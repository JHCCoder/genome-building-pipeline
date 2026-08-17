#!/bin/bash
#SBATCH -J finalize_annotation
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 01:00:00
#SBATCH --mem=32G
# --- Cluster-specific (Slurm on TSCC/UCSD) — adjust for your scheduler ---
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --no-requeue

# Load shared configuration (config.sh at the repo root)
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

source ~/.bashrc   # initialize conda (adjust for your setup)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
cd "$SCRIPT_DIR"

set -euo pipefail

# --- De-novo / degu-specific file paths (edit for your run) ---
OLD_INPUT="$PEAKS2UTR_DIR/hifiasm_041425_denovoEnhanced_sorted.gff"
OLD_OUTPUT="$PEAKS2UTR_DIR/hifiasm_041425_denovoEnhanced_peaks2utr_sorted_perfect.gff"
MERGED="output/hifiasm-041425-denovoEnhanced_merged.gff3"
SUBSET_OUT="output/peaks2utr_subset_out.gff3"
FINAL="output/hifiasm-041425-denovoEnhanced_peaks2utr_sorted.gff3"
AGAT_OUT="${FINAL%.gff3}.agat.gff3"

# 1. Assemble final annotation (gene models + peaks2utr UTRs).
#    The subset peaks2utr output must already exist (see peaks2utr_subset_rerun.sh).
test -s "$SUBSET_OUT" || { echo "ERROR: $SUBSET_OUT missing/empty"; exit 1; }
echo "=== step 1: assemble_final.py ==="
python3 assemble_final.py \
  --old-input "$OLD_INPUT" --old-output "$OLD_OUTPUT" \
  --merged "$MERGED" --subset-out "$SUBSET_OUT" --final "$FINAL"

# 2. AGAT normalize (rebuild gene->mRNA->CDS/exon linkage + add introns/start-stop)
echo "=== step 2: AGAT normalize ==="
conda activate "$ENV_AGAT"
agat_convert_sp_gxf2gxf.pl -g "$FINAL" -o "$AGAT_OUT"

# 3. verify the AGAT-normalized final
echo "=== step 3: verification ==="
export AGAT_OUT
python3 - <<'EOF'
import os
import re
from collections import Counter
path = os.environ["AGAT_OUT"]
gene_ids = set(); parentable = set(); utrs = []
mrna_orphan = mrna_none = 0
cdspar = Counter(); exonpar = Counter(); feat = Counter()
for line in open(path):
    if line.startswith("#"): continue
    p = line.rstrip().split("\t")
    if len(p) < 9: continue
    f, a = p[2], p[8]; feat[f] += 1
    idm = re.search(r"ID=([^;]+)", a); par = re.search(r"Parent=([^;]+)", a)
    if f == "gene":
        if idm: gene_ids.add(idm.group(1))
    elif f == "mRNA":
        if idm: parentable.add(idm.group(1))
        if par and par.group(1) not in gene_ids: mrna_orphan += 1
        elif not par: mrna_none += 1
    elif f in ("transcript","ncRNA","lnc_RNA","tRNA","rRNA","snRNA","snoRNA","C_gene_segment","V_gene_segment","D_gene_segment","J_gene_segment"):
        if idm: parentable.add(idm.group(1))
    elif f == "CDS":
        if par and par.group(1) in parentable: pass
        elif par: cdspar["BAD"] += 1
        else: cdspar["NONE"] += 1
    elif f == "exon":
        if par and par.group(1) in parentable: pass
        elif par: exonpar["BAD"] += 1
        else: exonpar["NONE"] += 1
    elif f in ("three_prime_UTR","five_prime_UTR"):
        if par: utrs.append(par.group(1))
dangling = sum(1 for t in utrs if t not in parentable)
print(f"genes={feat['gene']} mRNA={feat['mRNA']} CDS={feat['CDS']} exon={feat['exon']} intron={feat['intron']} start_codon={feat['start_codon']} stop_codon={feat['stop_codon']}")
print(f"5UTR={feat['five_prime_UTR']} 3UTR={feat['three_prime_UTR']}")
print(f"orphan_mRNA={mrna_orphan} noParent_mRNA={mrna_none} orphan_CDS={cdspar['NONE']+cdspar['BAD']} orphan_exon={exonpar['NONE']+exonpar['BAD']} dangling_UTR={dangling}")
EOF

echo "=== finalize DONE ==="
