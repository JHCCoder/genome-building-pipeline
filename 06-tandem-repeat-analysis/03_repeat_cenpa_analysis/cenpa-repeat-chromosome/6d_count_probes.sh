#!/bin/bash
#SBATCH -J 0807_count_probes
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

## Count each family probe set in each library's 31-mer read DB (KMC intersect,
## -oc right => report the read count per probe). Emits per-probe counts TSV.

KMC=/tscc/projects/ps-renlab2/jhc103/toolshed/kmc/bin
SCRATCH=/tscc/lustre/ddn/scratch/jhc103/kmc_tmp
source /tscc/nfs/home/jhc103/.bashrc
conda activate bulk-HiC-processing

cd /tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-repeat-chromosome
source scripts/config.sh
init_dirs

for s in XG_150 XG_151 XG_152 XG_153; do
    check_file "${PROBE_DIR}/${s}_reads_k31.kmc_suf" "$s read DB (run 6b first)"
done

for fam in 349 195 389; do
    probes_fa="${PROBE_DIR}/probes_${fam}.fa"
    check_file "$probes_fa" "probes_${fam} (run 6c first)"
    # build probe DB (each 31-mer once)
    "${KMC}/kmc" -k31 -m4 -sm -ci1 -fm "$probes_fa" \
        "${PROBE_DIR}/probes_${fam}_db" "$SCRATCH" >/dev/null 2>&1
    for s in XG_150 XG_151 XG_152 XG_153; do
        out="${PROBE_DIR}/probes_${fam}_x_${s}"
        if [[ ! -s "${out}.kmc_suf" ]]; then
            "${KMC}/kmc_tools" simple "${PROBE_DIR}/probes_${fam}_db" \
                "${PROBE_DIR}/${s}_reads_k31" intersect "$out" -ocright >/dev/null 2>&1
        fi
        "${KMC}/kmc_dump" "$out" "${PROBE_DIR}/probes_${fam}_${s}_counts.tsv" 2>/dev/null
        log "probes_${fam} x ${s}: $(wc -l < "${PROBE_DIR}/probes_${fam}_${s}_counts.tsv") probes with read hits"
    done
done

log "=== 6d done ==="
