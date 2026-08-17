#!/bin/bash
# Usage: ./highlight_chrom.sh <simple-file> <species1-bed> <species2-bed> <output-file> <chrom1:prefix> <chrom2:prefix> ...

# Required arguments
simple_file="$1"
species1_bed="$2"
species2_bed="$3"
output_file="$4"
shift 4

# Create associative array for chromosome prefixes
declare -A chrom_prefixes
for arg in "$@"; do
    IFS=':' read -r chrom prefix <<< "$arg"
    chrom_prefixes["$chrom"]="$prefix"
done

# Temporary files
tmp1=$(mktemp)
tmp2=$(mktemp)

# Create gene-to-chromosome maps
awk '{print $4,$1}' "$species1_bed" > "$tmp1"  # gene chrom
awk '{print $4,$1}' "$species2_bed" > "$tmp2"

# Process anchors file
while read -r line; do
    # Extract first 4 genes
    genes=($(echo "$line" | cut -f1-4))
    
    # Check chromosome for each gene
    prefixes=()
    for gene in "${genes[@]}"; do
        chrom=$(awk -v g="$gene" '$1==g {print $2}' "$tmp1" "$tmp2" | head -1)
        [[ -n "${chrom_prefixes[$chrom]}" ]] && prefixes+=("${chrom_prefixes[$chrom]}")
    done
    
    # Remove duplicate prefixes and sort for consistency
    unique_prefixes=$(printf "%s\n" "${prefixes[@]}" | sort -u | tr '\n' ' ' | sed 's/ $//')
    
    # Add prefixes if found
    if [[ -n "$unique_prefixes" ]]; then
        echo "${unique_prefixes// /}${line}"
    else
        echo "$line"
    fi
done < "$simple_file" > "$output_file"

# Cleanup
rm "$tmp1" "$tmp2"

echo "Created highlighted anchors file: $output_file"
echo "Chromosome prefixes:"
for chrom in "${!chrom_prefixes[@]}"; do
    echo "  $chrom → ${chrom_prefixes[$chrom]}"
done
