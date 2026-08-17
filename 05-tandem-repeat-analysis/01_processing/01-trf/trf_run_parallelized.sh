#!/bin/bash
#SBATCH -J 061025_run_trf_jobNumber
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 12:00:00
#SBATCH --mem=32G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user you@example.com
#SBATCH --no-requeue

# Load required modules
## module purge
module load shared
module load gpu/0.17.3
module load parallel/20210922-hzmc3me
module load pyfasta/0.5.2

# Set variables
INFILE="assembly_final.sorted.headerRenamed.chrAssigned"  # OMIT .fa extension
OUTPUT_DIR="/tscc/lustre/ddn/scratch/jhc103/trf-output"  # All outputs go here
INPUT_DIR="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/code/command-line-script/genome-annotation/trf-tandem-repeat"
num_chunk="32"

# Step 1: Split input FASTA into 32 chunks
#pyfasta split -n $num_chunk "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked-chrNameAssigned/${INFILE}.fasta"
#mv "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked-chrNameAssigned/${INFILE}.*.fasta" .

# Step 2: Submit array of TRF jobs
echo "Submitting TRF jobs..."
JOBIDS=()
for chunk in "${INFILE}".??.fasta; do
  JOBIDS+=($(sbatch \
    --job-name="trf_${chunk##*.}" \
    --output="/tscc/nfs/home/jhc103/cluster-logs/trf_${chunk##*.}.%j.out" \
    --error="/tscc/nfs/home/jhc103/cluster-logs/trf_${chunk##*.}.%j.err" \
    --time=12:00:00 \
    --mem=64G \
    --partition=condo \
    --qos=condo \
    --account=csd788 \
    --no-requeue \
    --wrap="cd $OUTPUT_DIR && /tscc/projects/ps-renlab2/jhc103/toolshed/TRF-4.09.1/build/src/trf $INPUT_DIR/${chunk} 2 7 7 80 10 50 500 -d -l 20" \
    --parsable))
done

# Step 3: Monitor job completion
echo "Submitted ${#JOBIDS[@]} TRF jobs. Job IDs: ${JOBIDS[*]}"
echo "Waiting for jobs to complete..."
while squeue -u $USER -j ${JOBIDS// /,} | grep -q "TRF"; do
  sleep 60
done

# Step 4: Merge outputs
echo "Merging results..."
cat "${INFILE}".??.fa.*.mask > "${INFILE}.fa.mask"
cat "${INFILE}".??.fa.*.dat > "${INFILE}.fa.dat"

# Cleanup
echo "Cleaning up chunk files..."
rm "${INFILE}".??.fa "${INFILE}".??.fa.*.{dat,mask}

echo "TRF processing complete. Final outputs:"
echo "- Masked sequences: ${INFILE}.fa.mask"
echo "- Repeat annotations: ${INFILE}.fa.dat"


