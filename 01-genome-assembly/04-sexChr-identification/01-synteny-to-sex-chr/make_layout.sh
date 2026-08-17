#!/bin/bash

# Usage: ./make_layout.sh <species1-bed> <species2-bed> <anchors-simple-file> <output.layout>

species1_bed="$1"  # e.g., mouse.bed
species2_bed="$2"  # e.g., degus.bed
anchors_simple="$3"
output="$4"

# Function to capitalize first letter and remove extension
get_label() {
    local filename=$(basename "$1")
    local name=${filename%.*}  # Remove extension
    echo "${name^}"           # Capitalize first letter
}

label1=$(get_label "$species1_bed")
label2=$(get_label "$species2_bed")

# Create the layout file
cat > "$output" <<EOF
# y, xstart, xend, rotation, color, label, va, bed
 0.7,    0.1,    0.8,       0,      , $label1, top, $species1_bed
 0.3,    0.1,    0.8,       0,      , $label2, bottom, $species2_bed
# edges
e, 0, 1, $anchors_simple
EOF

echo "Generated layout file: $output"
echo "Labels used: $label1 and $label2"
