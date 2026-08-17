#!/bin/bash
#SBATCH -J 0807_kmer_probes
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 4:00:00
#SBATCH --mem=16G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user you@example.com
#SBATCH --no-requeue

## Prepare family k-mer dumps + genome counts for probe definition (6c).
##   kmc_dump  bin6/bin4/bin8 family DBs            -> per-family k-mer counts
##   intersect family DB x genome DB (-oc right)    -> genome count per family k-mer
##   kmc_dump  fam189/fam18 consensus 31-mers       -> L1 subfamily probes

KMC=/tscc/projects/ps-renlab2/jhc103/toolshed/kmc/bin
source /tscc/nfs/home/jhc103/.bashrc
conda activate bulk-HiC-processing

cd /tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-repeat-chromosome
source scripts/config.sh
init_dirs

THREADS=$SLURM_CPUS_PER_TASK
GEN_DB="${PROBE_DIR}/genome_k31"
check_file "${GEN_DB}.kmc_suf" "genome DB (run 6b first)"

# ---- per-family dumps + genome counts ----
for tag in bin6 bin4 bin8; do
    db="${PROBE_DIR}/${tag}_k31"
    check_file "${db}.kmc_suf" "$tag DB (run 6b first)"
    # family counts
    if [[ ! -s "${PROBE_DIR}/${tag}_counts.tsv" ]]; then
        "$KMC/kmc_dump" "$db" "${PROBE_DIR}/${tag}_counts.tsv" 2>/dev/null
        log "$tag counts: $(wc -l < "${PROBE_DIR}/${tag}_counts.tsv") k-mers"
    fi
    # genome count for each family k-mer (intersect with genome, keep right count)
    if [[ ! -s "${PROBE_DIR}/${tag}_genomecounts.tsv" ]]; then
        "$KMC/kmc_tools" simple "$db" "$GEN_DB" intersect "${PROBE_DIR}/${tag}_x_genome" -ocright 2>/dev/null
        "$KMC/kmc_dump" "${PROBE_DIR}/${tag}_x_genome" "${PROBE_DIR}/${tag}_genomecounts.tsv" 2>/dev/null
        log "$tag x genome: $(wc -l < "${PROBE_DIR}/${tag}_genomecounts.tsv") k-mers"
    fi
done

# ---- L1 subfamily probes from consensi ----
# Count 31-mers in each L1 consensus; these become fam-189 / fam-18 subfamily probes
for fam in "fam189:$L1_189_CONS" "fam18:$L1_18_CONS"; do
    tag="${fam%%:*}"; fa="${fam##*:}"
    if [[ ! -s "${PROBE_DIR}/${tag}_counts.tsv" ]]; then
        check_file "$fa" "$tag consensus"
        # count 31-mers in the consensus with KMC
        "$KMC/kmc" -k31 -m2 -sm -ci1 -fm "$fa" "${PROBE_DIR}/${tag}_k31" \
            /tscc/lustre/ddn/scratch/jhc103/kmc_tmp 2>/dev/null
        "$KMC/kmc_dump" "${PROBE_DIR}/${tag}_k31" "${PROBE_DIR}/${tag}_counts.tsv" 2>/dev/null
        log "$tag counts: $(wc -l < "${PROBE_DIR}/${tag}_counts.tsv") k-mers"
    fi
done

log "=== 6c probe setup done ==="
