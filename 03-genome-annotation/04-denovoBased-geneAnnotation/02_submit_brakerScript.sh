#!/bin/bash
#SBATCH -J braker3
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 48
#SBATCH -t 3-00:00:00
#SBATCH --mem=256G
# --- Cluster-specific (Slurm on TSCC/UCSD) — adjust for your scheduler ---
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user=you@example.com   # TODO: your email

# Load shared configuration (config.sh at the repo root)
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

module load singularitypro/3.11

# Braker3 runs inside its Singularity image; its working directory holds the
# masked assembly and protein evidence (both referenced by name below).
cd "$BRAKER3_DIR"

input_dir="$HISAT2_OUT_DIR/$HISAT2_INDEX_NAME/"
SRR_ids="SRR17216301,SRR17216302,SRR17216303,SRR17216304,SRR17216305,SRR17216306,SRR17216307,SRR17216308,SRR17216309,SRR17216310,SRR17216311,SRR17216312,SRR17216313,SRR17216314,SRR17216315,SRR17216293,SRR17216294,SRR17216295,SRR17216296,SRR17216297,SRR17216320,SRR17216319,SRR17216318,SRR17216317,SRR17216316,SRR17216299,SRR17216298,SRR17216321,SRR17216300"
input_SRR_bam=$(echo "$SRR_ids" | sed "s|,|.bam,$input_dir/|g; s|^|$input_dir/|; s|$|.bam|")
wd="$BRAKER3_DIR/braker-hifiasm-041425-vertebrate-plus-relatives-curated-genome/"
mkdir -p "$wd"

# De-novo gene prediction with Braker3: RNA-seq BAMs + vertebrate protein evidence
singularity exec -B "$SCRATCH_DIR/temp-dir-fast,$PROJ_ROOT" "$BRAKER3_SIF" braker.pl \
  --genome=assembly_final.sorted.headerRenamed.chrAssigned.masked.fasta \
  --prot_seq=Vertebrata_plus_relatives.fa \
  --bam="$input_SRR_bam" --gff3 --threads=48 --workingdir="$wd"
