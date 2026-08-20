#!/bin/bash
#SBATCH -J 081926_interferon_singlecopy_featureCounts
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 12:00:00
#SBATCH --mem=32G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user jhc103@ucsd.edu
#SBATCH --no-requeue

# ===========================================================================
# featureCounts for the 3 single-copy interferon genes (Gbp1, Casp1, Syncrip)
#
# These three are in the "Cellular Response to Type II Interferon" GO term but
# have no paralogs, so they were absent from paralog_families.gff and therefore
# from paralog_families_counts.tsv. This job quantifies their expression against
# the SAME 29 tissue RNA-seq BAMs used for the paralog counts, with the SAME
# featureCounts parameters (-g ID -t gene -p --primary), for a directly
# comparable result.
# ===========================================================================

source /tscc/nfs/home/jhc103/.bashrc
conda activate transcriptome-mapping

set -e

WORKDIR="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/code/github-code-to-share/07-segmental-duplication-analysis/02_analysis"
BAM_DIR="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/paralog-alignment-visualization"
cd "$WORKDIR"

ANNOTATION="interferon_singlecopy_genes.gff"
OUTPUT="interferon_singlecopy_counts.tsv"

echo "Started at: $(date)"
echo "Annotation: $ANNOTATION"
echo "Gene features: $(awk -F'\t' '$3=="gene"' "$ANNOTATION" | wc -l)"

# 29 tissue RNA-seq BAMs (same set used for paralog_families_counts.tsv)
BAM_FILES=( "$BAM_DIR"/*_Aligned.sortedByCoord.out.bam )
echo "BAM files found: ${#BAM_FILES[@]}"

featureCounts \
    -a "$ANNOTATION" \
    -o "$OUTPUT" \
    -g ID \
    -t gene \
    -p \
    -T 16 \
    --primary \
    "${BAM_FILES[@]}"

rc=$?
echo "featureCounts exit code: $rc"
echo "Finished at: $(date)"
ls -lh "$OUTPUT" "${OUTPUT}.summary" 2>/dev/null
echo "DONE"
