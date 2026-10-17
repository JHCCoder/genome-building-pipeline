#!/bin/bash
#SBATCH -J 101326_blast_linkedPeaks_octDeg1
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 6:00:00
#SBATCH --mem=32G
# --- Cluster-specific (Slurm on TSCC/UCSD) — adjust for your scheduler ---
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user=you@example.com   # TODO: your email
#SBATCH --no-requeue

# Step 2: check each candidate "new" linked peak (step 1: < 50% covered by the genome
# alignment) directly against the OctDeg1.0 sequence with BLASTN, because whole-genome
# minimap2 alignment misses collapsed and high-copy repeats.
#
# Full peak sequences (median ~870 bp) are searched with megablast, no DUST filtering (so
# low-complexity peaks are not silently dropped), keeping HSPs with >= 90% identity.
# Step 3 classifies each peak by the number of OctDeg1.0 loci covering >= 80% of it.
#
# Outputs: $FIG3D_CCRE_DIR/02_octdeg1_blastdb/  and  02_candidates_vs_octdeg1.blastn.tsv
set -euo pipefail

if [[ -n "${SLURM_JOB_ID:-}" && -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    _script_dir="$SLURM_SUBMIT_DIR"
else
    _script_dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
fi
_repo_root="$_script_dir"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

OUT="$FIG3D_CCRE_DIR"
QUERY="$OUT/01_union_unaligned_candidates.fa"
[[ -s "$QUERY" ]] || { echo "ERROR: run 01_classify_linked_peaks.sh first" >&2; exit 1; }

DB="$OUT/02_octdeg1_blastdb/octdeg1"
if [[ ! -s "$DB.nsq" && ! -s "$DB.00.nsq" ]]; then
  mkdir -p "$(dirname "$DB")"
  "$MAKEBLASTDB_BIN" -in "$OCTDEG1_FASTA" -dbtype nucl -parse_seqids -out "$DB" -title "OctDeg1.0 RefSeq genome"
fi

TSV="$OUT/02_candidates_vs_octdeg1.blastn.tsv"
[[ -e "$TSV" ]] && { echo "ERROR: $TSV exists; move it aside before re-running" >&2; exit 1; }
"$BLASTN_BIN" -task megablast -query "$QUERY" -db "$DB" \
  -dust no -perc_identity 90 -evalue 1e-20 \
  -max_target_seqs 100 -num_threads 8 \
  -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen" \
  > "$TSV"

echo "queries: $(grep -c '>' "$QUERY")  HSPs: $(wc -l < "$TSV")  queries with any HSP: $(cut -f1 "$TSV" | sort -u | wc -l)"
echo "blast finished: $(date)"
