#!/bin/bash

file_list="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/file_list_111824_deepSeq_403Male.txt"

echo "Starting job submission..."
count=0

while IFS=$'\t' read -r read_file1 read_file2 sampleid; do
    count=$((count + 1))
    echo "Processing job $count:"
    echo "  R1: $read_file1"
    echo "  R2: $read_file2"
    echo "  Sample: $sampleid"
    
    job_id=$(sbatch 01.proc_bulk_hic_sbatch_for_multiSubmission_degu_coolFiles.sh "$read_file1" "$read_file2" "$sampleid")
    echo "  Submitted job: $job_id"
    echo "---"
done < "$file_list"

echo "Total jobs submitted: $count"
