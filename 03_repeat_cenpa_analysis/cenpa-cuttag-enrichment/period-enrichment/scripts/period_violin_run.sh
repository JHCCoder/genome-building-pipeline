#!/bin/bash
#SBATCH --job-name=trf_period_violin
#SBATCH --output=/tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH --error=/tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH --account=csd788
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=you@example.com

# ============================================================================
# period_violin_run.sh
# Run the Δ CENP-A signal (fg − permutation null) violin plot, all 10 bins.
# Reads raw signal + null files directly (independent of 3_run_analysis.sh).
# Output: plots/period_violin_cenpa.{pdf,png,svg}
# ============================================================================

source /tscc/nfs/home/jhc103/.bashrc
conda activate r-visualizations

set -euo pipefail

cd /tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/period-enrichment

echo "[$(date)] Starting period_violin.R"
Rscript scripts/period_violin.R
echo "[$(date)] Done"