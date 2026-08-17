#!/bin/bash
# Usage: ./highlight_chrom.sh <simple-file> <species1-bed> <species2-bed> <output-file> [chromosomes]

# Required arguments
simple_file="$1"
species1_bed="$2"
species2_bed="$3"
output_file="$4"

# Optional: Chromosomes to highlight (default: ChrX)
highlight_chroms="${5:-ChrX}"  # Can specify multiple: "ChrX ChrY group1"

# Temporary files
tmp1=$(mktemp)
tmp2=$(mktemp)

# Step 1: Create gene-to-chromosome maps
awk '{print $4,$1}' "$species1_bed" > "$tmp1"  # gene chrom
awk '{print $4,$1}' "$species2_bed" > "$tmp2"

# Step 2: Process anchors file
while read -r line; do
    # Extract first 4 genes
    genes=($(echo "$line" | cut -f1-4))
    
    # Check if any gene is on highlighted chromosomes
    highlight=0
    for gene in "${genes[@]}"; do
        chrom=$(awk -v g="$gene" '$1==g {print $2}' "$tmp1" "$tmp2" | head -1)
        for hc in $highlight_chroms; do
            if [[ "$chrom" == "$hc" ]]; then
                highlight=1
                break 2
            fi
        done
    done
    
    # Add g* prefix if needed
    [[ $highlight -eq 1 ]] && echo -n "g*"
    echo "$line"
done < "$simple_file" > "$output_file"

# Cleanup
rm "$tmp1" "$tmp2"

echo "Created highlighted anchors file: $output_file"
echo "Highlighted chromosomes: $highlight_chroms"
