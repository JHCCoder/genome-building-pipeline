#!/usr/bin/env python
import sys
from collections import defaultdict

def read_gff(gff_file, gene_interest):
    """Parse GFF to get gene-to-transcript mappings"""
    gene_to_transcripts = defaultdict(list)
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            attrs = fields[8]
            if 'gene_name='+gene_interest+';' in attrs or 'gene='+gene_interest+';' in attrs:
                transcript_id = None
                for attr in attrs.split(';'):
                    if attr.startswith('ID=') or attr.startswith('Parent='):
                        transcript_id = attr.split('=')[1]
                        break
                if transcript_id:
                    gene_to_transcripts[gene_interest].append(transcript_id)
    return gene_to_transcripts

def find_blocks(blocks_file, transcripts, output_file, up_num=0, down_num=0):
    """Find all blocks containing the target transcripts plus upstream/downstream blocks"""
    found_blocks = set()
    all_blocks = []
    target_indices = set()
    
    # First pass: read all blocks and identify target positions
    with open(blocks_file) as fin:
        for i, line in enumerate(fin):
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            all_blocks.append(line)
            # Check if either gene in the block matches our transcripts
            if any(t in parts[0] for t in transcripts) or any(t in parts[1] for t in transcripts):
                target_indices.add(i)
    
    # Second pass: write blocks including upstream/downstream
    with open(output_file, 'w') as fout:
        for i, line in enumerate(all_blocks):
            # Check if this block is within range of any target
            should_write = False
            for target_idx in target_indices:
                if i >= target_idx - up_num and i <= target_idx + down_num:
                    should_write = True
                    break
            
            if should_write:
                fout.write(line)
                parts = line.strip().split()
                if len(parts) >= 2:
                    found_blocks.add((parts[0], parts[1]))
    
    return found_blocks

def main():
    # Check if we have at least the minimum required arguments
    if len(sys.argv) < 5:
        sys.exit("Usage: {} <gff_file> <blocks_file> <gene_name> <output_file> [upstream_num] [downstream_num]\n"
                 "       upstream_num and downstream_num default to 0 if not provided".format(sys.argv[0]))
    
    gff_file = sys.argv[1]
    blocks_file = sys.argv[2]
    gene_name = sys.argv[3]  # In your case "App"
    output_file = sys.argv[4]
    
    # Set defaults for upstream and downstream numbers
    up_num = 0
    down_num = 0
    
    # Parse optional arguments if provided
    if len(sys.argv) >= 6:
        try:
            up_num = int(sys.argv[5])
        except ValueError:
            sys.exit("Error: upstream_num must be an integer")
    
    if len(sys.argv) >= 7:
        try:
            down_num = int(sys.argv[6])
        except ValueError:
            sys.exit("Error: downstream_num must be an integer")
    
    print(f"Using upstream blocks: {up_num}, downstream blocks: {down_num}")

    # 1. Find all transcripts for the gene
    gene_transcripts = read_gff(gff_file, gene_name)
    if not gene_transcripts:
        sys.exit("Error: No transcripts found for gene {}".format(gene_name))
    
    print("Found {} transcripts for gene {}: {}".format(
        len(gene_transcripts[gene_name]), gene_name, ", ".join(gene_transcripts[gene_name])))

    # 2. Find all blocks containing these transcripts (plus upstream/downstream)
    found = find_blocks(blocks_file, gene_transcripts[gene_name], output_file, up_num, down_num)
    print("Found {} blocks containing {} gene transcripts (including {} upstream and {} downstream blocks)".format(
        len(found), gene_name, up_num, down_num))
    print("Output written to {}".format(output_file))

if __name__ == '__main__':
    main()
