# Or if using STAR (better for splice junctions)
STAR --runMode genomeGenerate \
     --genomeDir star_index \
     --genomeFastaFiles ~/ps-renlab2-link/degu-genome-assembly-proj/data/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked-chrNameAssigned/assembly_final.sorted.headerRenamed.chrAssigned.mito.fasta \
     --sjdbGTFfile cct7_complete_annotations.gff \
     --genomeSAindexNbases 12
