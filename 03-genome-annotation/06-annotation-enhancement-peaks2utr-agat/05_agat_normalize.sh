#!/bin/bash
#SBATCH -J agat_normalize
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
#SBATCH --mail-type END
#SBATCH --mail-user=you@example.com   # TODO: your email

# Load shared configuration (config.sh at the repo root)
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

set -euo pipefail

# 04_assemble_final.py must have produced the UTR-annotated GFF.
test -s "$FINAL_UTR_GFF" || { echo "ERROR: $FINAL_UTR_GFF missing/empty (run 04_assemble_final.py first)"; exit 1; }

source ~/.bashrc   # initialize conda (adjust for your setup)
conda activate "$ENV_AGAT"

# Normalize the assembled annotation: rebuild gene(mRNA(CDS/exon/UTR/intron))
# parentage and add introns + start/stop codons.
agat_convert_sp_gxf2gxf.pl \
  -g "$FINAL_UTR_GFF" \
  -o "$FINAL_AGAT_GFF"

# Verify the AGAT-normalized final (no orphan features / dangling UTR Parents).
python3 - <<'EOF'
import re
from collections import Counter
import os
path = os.environ["FINAL_AGAT_GFF"]
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

echo "AGAT normalization DONE: $FINAL_AGAT_GFF"
