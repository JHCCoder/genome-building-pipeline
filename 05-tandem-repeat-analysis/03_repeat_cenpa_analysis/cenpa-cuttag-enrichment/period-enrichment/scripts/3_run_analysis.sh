#!/bin/bash
#SBATCH --job-name=trf_period_R
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
# 3_run_analysis.sh
# Run the period-size stratified CENP-A enrichment R analysis.
# Input:  trf_signal_*.tsv + trf_bg_signal_*.tsv (from steps 2 and 2b)
# Output: plots/ and results/ in period-enrichment/
# ============================================================================

source /tscc/nfs/home/jhc103/.bashrc
conda activate r-visualizations

set -euo pipefail

cd /tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/period-enrichment

echo "[$(date)] Starting 3_analyze_period_enrichment.R"
Rscript scripts/3_analyze_period_enrichment.R
echo "[$(date)] Done"
