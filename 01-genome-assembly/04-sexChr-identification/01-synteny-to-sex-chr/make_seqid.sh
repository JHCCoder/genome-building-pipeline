#!/bin/bash

# Usage: ./generate_seqid.sh <species1.bed> <species2.bed> <anchor_file> <output.seqid>

bed1="$1"   # BED file for species 1 (e.g., mouse)
bed2="$2"   # BED file for species 2 (e.g., degus)
anchors="$3" # Anchor file
output="$4"  # Output seqid file

# Step 1: Extract all unique transcript IDs from the anchor file
transcripts=$(cut -f1-4 "$anchors" | tr '\t' '\n' | sort -u)

# Step 2: Write transcripts to a temporary file for grep
tmp_transcripts=$(mktemp)
echo "$transcripts" > "$tmp_transcripts"

# Step 3: Find chromosomes from species1 BED file with these transcripts
chroms1=$(grep -Fwf "$tmp_transcripts" "$bed1" | awk '{print $1}' | sort -u | tr '\n' ',' | sed 's/,$//')

# Step 4: Find chromosomes from species2 BED file with these transcripts
chroms2=$(grep -Fwf "$tmp_transcripts" "$bed2" | awk '{print $1}' | sort -u | tr '\n' ',' | sed 's/,$//')

# Step 5: Write to seqid file
echo "$chroms1" > "$output"
echo "$chroms2" >> "$output"

# Cleanup
rm "$tmp_transcripts"

echo "Generated seqid file: $output"
echo "Species 1 chromosomes: $chroms1"
echo "Species 2 chromosomes: $chroms2"
