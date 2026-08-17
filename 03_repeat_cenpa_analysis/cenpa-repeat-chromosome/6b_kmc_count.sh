#!/bin/bash
#SBATCH -J 0807_kmc_count
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 8:00:00
#SBATCH --mem=64G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user you@example.com
#SBATCH --no-requeue

## KMC 31-mer counting for the k-mer leg.
## Builds KMC DBs for:
##   genome (chr1-28,X,Y)              -- for probe specificity (genome counts)
##   each repeat family array seq      -- bin6(349) bin4(195) bin8(389)
##   each CUT&Tag library reads        -- XG_150/151/152/153 (trimmed fastq)
## Uses /tscc/lustre/ddn/scratch/jhc103 for KMC temp.

KMC=/tscc/projects/ps-renlab2/jhc103/toolshed/kmc/bin
SCRATCH=/tscc/lustre/ddn/scratch/jhc103/kmc_tmp
mkdir -p "$SCRATCH"

source /tscc/nfs/home/jhc103/.bashrc
conda activate bulk-HiC-processing

cd /tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-repeat-chromosome
source scripts/config.sh
init_dirs

K=31; THREADS=$SLURM_CPUS_PER_TASK

# ---- genome k-mer DB (chr only; big) ----
GEN_DB="${PROBE_DIR}/genome_k31"
if [[ ! -s "${GEN_DB}.kmc_suf" ]]; then
    log "Counting genome 31-mers (this is the big one)..."
    # subset fasta to chr1-28,X,Y
    FA_CHR="${PROBE_DIR}/assembly_chr.fa"
    if [[ ! -s "$FA_CHR" ]]; then
        # subset the assembly to chr1-28,X,Y with samtools faidx
        export PATH="/tscc/projects/ps-renlab2/jhc103/miniconda3-storage/envs/bulk-HiC-processing/bin:$PATH"
        samtools faidx "$ASSEMBLY_FASTA" $(awk '{print $1}' "$CHROM_SIZES") > "$FA_CHR" 2>/dev/null
        log "genome fasta (chr only): $(du -sh "$FA_CHR" | cut -f1)"
    fi
    "$KMC/kmc" -k"$K" -m48 -sm -t"$THREADS" -ci1 -fm "$FA_CHR" "$GEN_DB" "$SCRATCH" \
        > "${PROBE_DIR}/kmc_genome.log" 2>&1
    log "genome DB: $(grep 'unique counted' "${PROBE_DIR}/kmc_genome.log" || echo done)"
fi

# ---- family array-sequence DBs ----
for spec in "bin6:$PROBE_DIR/bin6_seq.fa" "bin4:$PROBE_DIR/bin4_seq.fa" "bin8:$PROBE_DIR/bin8_seq.fa"; do
    tag="${spec%%:*}"; fa="${spec##*:}"
    db="${PROBE_DIR}/${tag}_k31"
    if [[ ! -s "${db}.kmc_suf" ]]; then
        check_file "$fa" "$tag fasta (run 6a first)"
        log "Counting $tag 31-mers..."
        "$KMC/kmc" -k"$K" -m8 -sm -t"$THREADS" -ci1 -fm "$fa" "$db" "$SCRATCH" \
            > "${PROBE_DIR}/kmc_${tag}.log" 2>&1
        log "$tag DB done: $(grep 'unique counted' "${PROBE_DIR}/kmc_${tag}.log" | awk '{print $NF}')"
    fi
done

# ---- library read DBs ----
declare -A READ_ONE=( [XG_150]="$READ_150_1" [XG_151]="$READ_151_1" [XG_152]="$READ_152_1" )
for s in XG_150 XG_151 XG_152; do
    db="${PROBE_DIR}/${s}_reads_k31"
    if [[ ! -s "${db}.kmc_suf" ]]; then
        log "Counting ${s} read 31-mers..."
        "$KMC/kmc" -k"$K" -m32 -sm -t"$THREADS" -ci1 -fq "${READ_ONE[$s]}" "$db" "$SCRATCH" \
            > "${PROBE_DIR}/kmc_${s}.log" 2>&1
        log "${s} read DB done"
    fi
done

# XG_153 reads (R1 only, consistent with others)
DB="${PROBE_DIR}/XG_153_reads_k31"
if [[ ! -s "${DB}.kmc_suf" ]]; then
    check_file "$READ_153_1" "XG_153 reads"
    "$KMC/kmc" -k"$K" -m32 -sm -t"$THREADS" -ci1 -fq "$READ_153_1" "$DB" "$SCRATCH" \
        > "${PROBE_DIR}/kmc_XG_153.log" 2>&1
    log "XG_153 read DB done"
fi

log "=== 6b done ==="
