#!/bin/bash
# samtools view -H 403_ear_deep_file_hifi_041425_1.bam | grep '^@SQ' | cut -f 2 | sed 's/SN://g' > chromosome_list.txt then delete while only leaving chr chromosomes
# This script submits an array job, one task per chromosome.

# Read the chromosome list into a Bash array
mapfile -t CHROMOSOMES < chromosome_list.txt

# Get the number of chromosomes
NUM_CHROM=${#CHROMOSOMES[@]}

# Create a temporary file to hold all the chromosome arguments, one per line.
# This is the key file for the job array.
printf "%s\n" "${CHROMOSOMES[@]}" > chrom_args_for_array.txt


# Submit the job array.
# The --array option will create a task for each line (1-N).
# The command for each task will be to run the SBATCH script and pass it the line from the file.
sbatch --array=1-${NUM_CHROM}%20 \
       --dependency=singleton \
       split_chromosome.sh \
       $(sed -n "${SLURM_ARRAY_TASK_ID}p" chrom_args_for_array.txt)

echo "finish submission"
