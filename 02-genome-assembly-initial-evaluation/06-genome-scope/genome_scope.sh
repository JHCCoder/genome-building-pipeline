#!/bin/bash
#SBATCH -J genome_scope
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 6:00:00
#SBATCH --mem=175G
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

# ====== RUN PARAMETERS (edit these for each run) ====== #
KMER=21               # k-mer size
MIN_COVERAGE=1        # KMC -ci lower coverage cutoff
MAX_COVERAGE=1200000  # KMC -cs upper coverage cutoff
MEM=170               # KMC max memory (GB)
THREADS=8
SAMPLE="male403"      # label used in the output directory name

# Count k-mers from Illumina WGS reads (KMC) and emit a frequency histogram.
# The histogram is NOT fitted here — upload the resulting .histo to the
# GenomeScope web interface to get the genome-size estimate (see README).

out_dir="${GENOMESCOPE_OUT_DIR}/${SAMPLE}-${KMER}kmer"
mkdir -p "$out_dir"
cd "$out_dir" || exit 1

prefix="reads${MAX_COVERAGE}"

# Build the KMC input file list from the read glob
realpath $GENOMESCOPE_READS > FILES

"$KMC_BIN" -k"$KMER" -t"$THREADS" -m"$MEM" -ci"$MIN_COVERAGE" -cs"$MAX_COVERAGE" @FILES "$prefix" "$GENOMESCOPE_TEMP_DIR"

"$KMC_TOOLS_BIN" transform "$prefix" histogram "${prefix}.histo" -cx"$MAX_COVERAGE"

echo "Done. Upload ${out_dir}/${prefix}.histo to https://genomescope.org/ for the genome-size estimate."
