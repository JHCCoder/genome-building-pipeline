#!/usr/bin/env python
import sys
import re
from collections import defaultdict

def read_gff(gff_file, target_gene):
    """Parse GFF to get gene-to-transcript mappings with precise attribute parsing"""
    gene_to_transcripts = defaultdict(list)
    transcript_to_gene = {}
    
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
                
            attrs = fields[8]
            attr_dict = dict(re.findall(r'([^=;\s]+)=([^=;\s]+)', attrs))
            
            # Check both gene_name and gene attributes
            if attr_dict.get('gene_name') == target_gene or attr_dict.get('gene') == target_gene:
                transcript_id = attr_dict.get('ID') or attr_dict.get('Parent')
                if transcript_id:
                    gene_to_transcripts[target_gene].append(transcript_id)
                    transcript_to_gene[transcript_id] = target_gene
    
    return gene_to_transcripts, transcript_to_gene

def find_blocks(blocks_file, transcripts, output_file):
    """Find blocks containing target transcripts and their syntenic partners"""
    found_blocks = set()
    linked_transcripts = set(transcripts)  # Start with our target transcripts
    
    # First pass: find direct matches
    with open(blocks_file) as fin:
        for line in fin:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
                
            t1, t2 = parts[0], parts[1]
            if t1 in transcripts or t2 in transcripts:
                found_blocks.add((t1, t2))
                linked_transcripts.update([t1, t2])  # Add both transcripts
    
    # Second pass: find any blocks containing the expanded transcript set
    with open(blocks_file) as fin, open(output_file, 'w') as fout:
        for line in fin:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
                
            t1, t2 = parts[0], parts[1]
            if t1 in linked_transcripts or t2 in linked_transcripts:
                fout.write(line)
                
    return linked_transcripts

def main():
    if len(sys.argv) != 5:
        sys.exit("Usage: {} <gff_file> <blocks_file> <gene_name> <output_file>".format(sys.argv[0]))
    
    gff_file = sys.argv[1]
    blocks_file = sys.argv[2]
    gene_name = sys.argv[3]
    output_file = sys.argv[4]

    # 1. Find all transcripts for the target gene
    gene_transcripts, transcript_map = read_gff(gff_file, gene_name)
    if not gene_transcripts:
        sys.exit(f"Error: No transcripts found for gene {gene_name}")
    
    target_transcripts = gene_transcripts[gene_name]
    print(f"Found {len(target_transcripts)} transcripts for {gene_name}:")
    print(", ".join(target_transcripts))

    # 2. Find all blocks containing these transcripts and their partners
    linked_transcripts = find_blocks(blocks_file, target_transcripts, output_file)
    print(f"\nFound {len(linked_transcripts)} linked transcripts total")
    print(f"Output written to {output_file}")

    # 3. Generate companion BED files for visualization
    print("\nNext steps:")
    print(f"1. Filter BED files using:")
    print(f"   grep -f <(cut -f1 {output_file} | sort -u) mouse.bed > mouse_app.bed")
    print(f"   grep -f <(cut -f2 {output_file} | sort -u) degus.bed > degus_app.bed")

if __name__ == '__main__':
    main()
