#!/bin/bash
#SBATCH -J evo_raxml
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 1-00:00:00
#SBATCH --mem=120G
# --- Cluster-specific (Slurm on TSCC/UCSD) — adjust for your scheduler ---
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p platinum
#SBATCH -q hcp-csd788
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user=you@example.com   # TODO: your email

# =============================================================================
# Step 5 of the evolutionary-placement pipeline.
#
# Infers a maximum-likelihood species tree from the concatenated amino-acid
# supermatrix with RAxML, using the JTT substitution model with gamma-distributed
# rate heterogeneity (PROTGAMMAJTT), 100 rapid bootstrap replicates followed by a
# search for the best-scoring ML tree (-f a -N 100 -T 16).
#
# RUN PARAMETERS
# --------------
#   MODEL           substitution model (default: PROTGAMMAJTT)
#   NBOOT           number of rapid bootstrap replicates (default: 100)
#   RAXML_THREADS   threads for the PTHREADS build (default: 16)
#   PARSIMONY_SEED  -p random seed (default: 12345; not specified in the methods)
#   BOOT_SEED       -x bootstrap seed (default: 12345; not specified in the methods)
#   RUN_NAME        -n run name used to prefix RAxML output files
# =============================================================================

# Load shared configuration (config.sh at the repo root)
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

MODEL="PROTGAMMAJTT"
NBOOT=100
RAXML_THREADS=16
PARSIMONY_SEED=12345
BOOT_SEED=12345
RUN_NAME="degu_tree"

RAXML_DIR="$EVO_OUT_DIR/05_raxml"
mkdir -p "$RAXML_DIR"

source ~/.bashrc                      # initialize conda (adjust for your setup)
conda activate "$ENV_EVOLUTIONARY_TREE"   # RAxML 8.2.12

cd "$RAXML_DIR"

# -f a : rapid bootstrap + search for best-scoring ML tree in one run
# -N   : number of bootstrap replicates
# -T   : threads (PTHREADS build)
# The final tree with bootstrap support is RAxML_bipartitions.<RUN_NAME>.
raxmlHPC-PTHREADS-AVX2 \
    -s "$EVO_OUT_DIR/04_supermatrix/supermatrix.phy" \
    -m "$MODEL" \
    -f a \
    -N "$NBOOT" \
    -T "$RAXML_THREADS" \
    -p "$PARSIMONY_SEED" \
    -x "$BOOT_SEED" \
    -n "$RUN_NAME" < /dev/null

echo "Done. Best tree with bootstrap support: $RAXML_DIR/RAxML_bipartitions.$RUN_NAME"
