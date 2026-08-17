#!/bin/bash

# Usage: ./count_transcripts.sh <chromosome_name> <bed_file> <anchor_file>

chrom=$1
bed_file=$2
anchor_file=$3

# Step 1: Extract all transcripts for the given chromosome from the BED file
transcripts=$(awk -v chrom="$chrom" '$1 == chrom {print $4}' "$bed_file" | sort -u)

# Step 2: Count how many of these transcripts appear in the anchor file
count=$(grep -Fwf <(echo "$transcripts") "$anchor_file" | wc -l)

# Step 3: Print results
echo "Chromosome $chrom has $(echo "$transcripts" | wc -l) unique transcripts in the BED file."
echo "$count of these transcripts appear in the anchor file."
