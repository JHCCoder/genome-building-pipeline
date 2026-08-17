#!/bin/bash
#SBATCH -J 0807_kmer_seqs
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 2:00:00
#SBATCH --mem=8G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user you@example.com
#SBATCH --no-requeue

## Extract array sequences per repeat family from the assembly (for k-mer probes)
## Families: bin6=349, bin4=195, bin8=389

source /tscc/nfs/home/jhc103/.bashrc
conda activate bulk-HiC-processing

cd /tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-repeat-chromosome
source scripts/config.sh
init_dirs

for spec in "bin6:$BIN6_349" "bin4:$BIN4_195" "bin8:$BIN8_389"; do
    bin=$(echo "$spec" | cut -d: -f1)
    bed=$(echo "$spec" | cut -d: -f2)
    out="${PROBE_DIR}/${bin}_seq.fa"
    if [[ -s "$out" ]]; then
        log "SKIP $bin (exists)"
        continue
    fi
    check_file "$bed" "$bin arrays"
    bedtools getfasta -fi "$ASSEMBLY_FASTA" -bed "$bed" -name+ -s 2>/dev/null > "$out"
    log "$bin: $(grep -c '^>' "$out") sequences, $(du -sh "$out" | cut -f1)"
done

log "=== 6a done ==="
