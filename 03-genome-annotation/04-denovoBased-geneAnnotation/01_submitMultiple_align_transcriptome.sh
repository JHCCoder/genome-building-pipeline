#!/bin/bash
#SBATCH -J 042425_submit_hisat2_jobs
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH -t 02:00:00
#SBATCH --mem=32G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user jhc103@ucsd.edu

source /tscc/nfs/home/jhc103/.bashrc
conda activate genome-annotation

#export input_fa_dir="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/output/outputs-from-haphic-alignment/references_hifiasm_male403_hifiHiCMode_022425/04.build/"
#export assembly="hifiasm_022425_haphic"
export base_name="hifiasm_041425_haphic_masked_curated"
assembly_file='/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked-chrNameAssigned/assembly_final.sorted.headerRenamed.chrAssigned.masked.fasta' #'/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked/assembly_final.sorted.headerRenamed.fasta.masked' #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/output/outputs-from-haphic-alignment/references_hifiasm_male403_hifiHiCMode_022425/04.build/scaffolds.fa.masked" #Previous assembly: "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/output/outputs-from-repeatmasker/hifiasm-120624-haphic/scaffolds.fa.masked"


# Check if HISAT2 index exists
if [[ ! -f ${assembly_file}.1.ht2 ]]; then
    echo "HISAT2 index not found. Building index..."
    hisat2-build -p 12 --verbose ${assembly_file} ${base_name}
    if [[ $? -ne 0 ]]; then
        echo "Error: HISAT2 index build failed."
        exit 1
    fi
else
    echo "HISAT2 index already exists. Skipping index build."
fi

# Submit array job for alignment
SRR_list=(
    SRR17216301 SRR17216302 SRR17216303 SRR17216304 SRR17216305
    SRR17216306 SRR17216307 SRR17216308 SRR17216309 SRR17216310
    SRR17216311 SRR17216312 SRR17216313 SRR17216314 SRR17216315
    SRR17216293 SRR17216294 SRR17216295 SRR17216296 SRR17216297
    SRR17216320 SRR17216319 SRR17216318 SRR17216317 SRR17216316
    SRR17216299 SRR17216298 SRR17216321 SRR17216300
)

sbatch --array=0-$((${#SRR_list[@]} - 1)) hisat2_alignment.sh "${SRR_list[@]}"
