#!/bin/bash
#SBATCH --job-name=349bp_bg
#SBATCH --output=/tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH --error=/tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH --account=csd788
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=you@example.com

# ============================================================================
# 2_generate_background.sh
# Generate chromosome- and length-matched shuffled intervals
# N=1000 iterations per foreground interval
# ============================================================================

source /tscc/nfs/home/jhc103/.bashrc
conda activate bulk-HiC-processing

set -euo

SCRIPT_DIR="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/349-bp"
source "${SCRIPT_DIR}/config_349bp.sh"
init_349_dirs

log "=== Step 2: Generate shuffled background ==="

check_file "${FOREGROUND_MERGED}" "foreground merged BED"
check_file "${ASSEMBLY_FAI}" "assembly FAI"

# ============================================================================
# Generate shuffled background using bedtools shuffle
# ============================================================================
if [[ -f "${BG_SHUFFLED}" ]]; then
    log "Background BED already exists: ${BG_SHUFFLED}"
    N_LINES=$(wc -l < "${BG_SHUFFLED}")
    log "  ${N_LINES} lines"
else
    log "Generating ${N_SHUFFLES} iterations of chromosome-matched shuffled intervals"

    # Create chromosome sizes file (chr1-28, chrX only)
    CHR_SIZES_349="${DATA_349_DIR}/chrom_sizes_349.txt"
    grep -w -f <(echo "${CHROMOSOMES_349}" | tr ' ' '\n') "${ASSEMBLY_FAI}" | \
        awk '{print $1, $2}' OFS='\t' > "${CHR_SIZES_349}"
    log "Chromosome sizes: $(wc -l < ${CHR_SIZES_349}) chromosomes"

    # Get gap BED for exclusion (reuse from main pipeline)
    GAP_BED="${EXCLUSION_DIR}/assembly_gaps.bed"
    GAP_BED_349="${DATA_349_DIR}/assembly_gaps_349.bed"
    if [[ -f "${GAP_BED}" ]]; then
        grep -w -f <(echo "${CHROMOSOMES_349}" | tr ' ' '\n') "${GAP_BED}" > "${GAP_BED_349}" || true
        log "Gap BED: $(wc -l < ${GAP_BED_349}) regions"
    else
        # Create empty exclusion file if no gap BED
        touch "${GAP_BED_349}"
        log "No gap BED found, shuffling without exclusion"
    fi

    # Extract foreground intervals (chr, start, end, id)
    FG_SIMPLE="${DATA_349_DIR}/foreground_simple.bed"
    cut -f1-4 "${FOREGROUND_MERGED}" > "${FG_SIMPLE}"
    N_FG=$(wc -l < "${FG_SIMPLE}")
    log "Foreground intervals: ${N_FG}"

    # Generate shuffled intervals using bedtools shuffle
    # We need to shuffle each interval N_SHUFFLES times
    # Strategy: use a loop with different seeds, concatenate
    BG_TMP_DIR="${DATA_349_DIR}/bg_tmp"
    mkdir -p "${BG_TMP_DIR}"

    # Use GNU parallel if available, otherwise sequential
    log "Shuffling ${N_FG} intervals x ${N_SHUFFLES} iterations..."

    # Create shuffled copies using R (more control over per-chromosome matching)
    conda activate r-visualizations
    Rscript --no-save - \
        "${FG_SIMPLE}" "${CHR_SIZES_349}" "${GAP_BED_349}" \
        "${BG_SHUFFLED}" "${N_SHUFFLES}" "${SHUFFLE_SEED}" << 'REOF'
args <- commandArgs(trailingOnly = TRUE)
fg_file <- args[1]
chr_sizes_file <- args[2]
gap_file <- args[3]
outfile <- args[4]
n_shuffles <- as.integer(args[5])
seed <- as.integer(args[6])

set.seed(seed)

message("Loading foreground intervals")
fg <- read.table(fg_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                 col.names = c("chrom", "start", "end", "id"))
message("  ", nrow(fg), " intervals")

message("Loading chromosome sizes")
chr_sizes <- read.table(chr_sizes_file, header = FALSE, sep = "\t",
                        col.names = c("chrom", "size"))
chr_max <- setNames(chr_sizes$size, chr_sizes$chrom)

message("Loading gaps")
gaps <- NULL
if (file.exists(gap_file) && file.info(gap_file)$size > 0) {
    gaps <- read.table(gap_file, header = FALSE, sep = "\t",
                       col.names = c("chrom", "start", "end"))
    message("  ", nrow(gaps), " gap regions")
}

# Build gap intervals per chromosome
gaps_by_chr <- split(gaps, gaps$chrom)

# Check if a shuffled interval overlaps any gap
overlaps_gap <- function(chrom, start, end) {
    if (is.null(gaps) || nrow(gaps) == 0) return(FALSE)
    chr_gaps <- gaps_by_chr[[chrom]]
    if (is.null(chr_gaps) || nrow(chr_gaps) == 0) return(FALSE)
    any(start < chr_gaps$end & end > chr_gaps$start)
}

# Generate shuffled positions
message("Generating ", n_shuffles, " shuffled backgrounds per interval...")
results <- vector("list", nrow(fg))
pb <- txtProgressBar(min = 0, max = nrow(fg), style = 3)

for (i in seq_len(nrow(fg))) {
    chrom <- fg$chrom[i]
    width <- fg$end[i] - fg$start[i]
    interval_id <- fg$id[i]
    chr_len <- chr_max[chrom]

    if (is.na(chr_len) || chr_len <= width) {
        setTxtProgressBar(pb, i)
        next  # skip if chr not in sizes or interval > chr
    }

    max_start <- chr_len - width
    iter_rows <- vector("list", n_shuffles)

    for (j in seq_len(n_shuffles)) {
        # Try up to 100 times to avoid gaps
        for (attempt in 1:100) {
            new_start <- sample.int(max_start, 1)
            new_end <- new_start + width
            if (!overlaps_gap(chrom, new_start, new_end)) break
        }
        iter_rows[[j]] <- data.frame(
            chrom = chrom,
            start = new_start,
            end = new_end,
            iter = sprintf("iter_%03d", j),
            interval_id = interval_id,
            stringsAsFactors = FALSE
        )
    }
    results[[i]] <- do.call(rbind, iter_rows)
    setTxtProgressBar(pb, i)
}
close(pb)

bg <- do.call(rbind, results)
message("\nTotal shuffled intervals: ", nrow(bg))
message("Writing: ", outfile)
write.table(bg, outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
message("Done")
REOF

    # Clean up temp
    rm -rf "${BG_TMP_DIR}"

    N_BG=$(wc -l < "${BG_SHUFFLED}")
    log "Background BED created: ${N_BG} lines"
    log "Expected: $((N_FG * N_SHUFFLES)) lines"
fi

# ============================================================================
# QC: Check shuffle uniformity
# ============================================================================
log "QC: Checking shuffle uniformity..."

conda activate r-visualizations
Rscript --no-save - "${BG_SHUFFLED}" "${FOREGROUND_MERGED}" << 'REOF'
args <- commandArgs(trailingOnly = TRUE)
bg_file <- args[1]
fg_file <- args[2]

message("Loading background...")
bg <- read.table(bg_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                 col.names = c("chrom", "start", "end", "iter", "interval_id"))
message("  Background intervals: ", nrow(bg))

fg <- read.table(fg_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                 col.names = c("chrom", "start", "end", "id", "size"))
message("  Foreground intervals: ", nrow(fg))

# Check: does each foreground interval have n_shuffles background copies?
n_iter <- length(unique(bg$iter))
message("  Unique iterations: ", n_iter)

interval_counts <- table(bg$interval_id)
expected <- n_iter
missing <- sum(interval_counts != expected)
if (missing > 0) {
    message("  WARNING: ", missing, " intervals have != ", expected, " background copies")
} else {
    message("  OK: All intervals have exactly ", expected, " background copies")
}

# Check chromosome distribution
for (chr in unique(fg$chrom)) {
    fg_count <- sum(fg$chrom == chr)
    bg_count <- sum(bg$chrom == chr) / n_iter
    if (fg_count > 0) {
        cat(sprintf("  %-6s: fg=%3d  bg/iter=%7.1f\n", chr, fg_count, bg_count))
    }
}
REOF

log "=== Step 2 complete ==="
