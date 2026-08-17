module load parallel/20210922-hzmc3me
export DATADIR=`pwd`
INFILE="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked-chrNameAssigned/assembly_final.sorted.headerRenamed.chrAssigned.fasta" # OMIT EXTENSION!
export SCRIPTDIR=~
pyfasta split -n 32 $INFILE
parallel -j 32 qsub -q nice.q -v SEQFILE={} -sync y $SCRIPTDIR/trf.sh ::: $INFILE.??.fa
cat $INFILE.??.fa.*.mask > $INFILE.fa.mask
cat $INFILE.??.fa.*.dat > $INFILE.fa.dat
