#!/bin/bash
# ============================================================
# Monitor HiCAT sub-chunk jobs: check completion, failures, efficiency.
# Reads subchunk_manifest.tsv and checks each job via sacct/seff.
#
# Usage: bash monitor_subchunks.sh [--summary]
#   --summary: only print summary counts, no per-job detail
# ============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
MANIFEST="${SCRIPT_DIR}/HiCAT_genome/subchunk_manifest.tsv"

if [[ ! -f "${MANIFEST}" ]]; then
    echo "ERROR: Manifest not found: ${MANIFEST}"
    exit 1
fi

SUMMARY_ONLY=0
[[ "${1:-}" == "--summary" ]] && SUMMARY_ONLY=1

COMPLETED=0
FAILED=0
RUNNING=0
PENDING=0
UNKNOWN=0

declare -a FAILED_JOBS

echo "[$(date)] HiCAT sub-chunk monitor"
echo "================================"

# Skip header line
tail -n +2 "${MANIFEST}" | while IFS=$'\t' read -r sub_chunk chr parent sub_idx job_id mem cpu fasta_path; do
    # Get job state from sacct
    state=$(sacct -j "${job_id}" --format=State --noheader -P 2>/dev/null | head -1 || echo "UNKNOWN")

    # Normalize state
    state="${state:-UNKNOWN}"

    case "${state}" in
        COMPLETED)
            # Get efficiency stats
            seff_out=$(seff "${job_id}" 2>/dev/null || echo "")
            mem_used=$(echo "${seff_out}" | grep "Memory Utilized" | grep -oP '[\d.]+ [GM]B' || echo "?")
            cpu_eff=$(echo "${seff_out}" | grep "CPU Efficiency" | grep -oP '[\d.]+%' || echo "?")
            walltime=$(echo "${seff_out}" | grep "Wall-clock" | grep -oP '[\d-]+:[\d:]+' || echo "?")

            # Check for output files
            out_dir="${SCRIPT_DIR}/HiCAT_genome/${sub_chunk}"
            has_out=0
            ls "${out_dir}"/out_* >/dev/null 2>&1 && has_out=1

            exit_code=$(sacct -j "${job_id}" --format=ExitCode --noheader -P 2>/dev/null | head -1 | cut -d: -f1 || echo "?")

            if [[ "${exit_code}" == "0" ]] && [[ ${has_out} -eq 1 ]]; then
                if [[ "${SUMMARY_ONLY}" != "1" ]]; then
                    echo "  OK  ${sub_chunk} (${job_id}): ${walltime}, mem=${mem_used}, cpu=${cpu_eff}"
                fi
                echo "COMPLETED" >> /tmp/hicat_monitor_$$
            else
                if [[ "${SUMMARY_ONLY}" != "1" ]]; then
                    echo "  FAIL ${sub_chunk} (${job_id}): exit=${exit_code}, out_files=${has_out}, mem=${mem_used}"
                fi
                echo "FAILED" >> /tmp/hicat_monitor_$$
            fi
            ;;
        RUNNING)
            echo "RUNNING" >> /tmp/hicat_monitor_$$
            ;;
        PENDING)
            echo "PENDING" >> /tmp/hicat_monitor_$$
            ;;
        FAILED|TIMEOUT|CANCELLED|OUT_OF_MEMORY|NODE_FAIL)
            if [[ "${SUMMARY_ONLY}" != "1" ]]; then
                echo "  FAIL ${sub_chunk} (${job_id}): state=${state}"
            fi
            echo "FAILED" >> /tmp/hicat_monitor_$$
            ;;
        *)
            echo "UNKNOWN" >> /tmp/hicat_monitor_$$
            ;;
    esac
done

# Count
COMPLETED=$(grep -c "COMPLETED" /tmp/hicat_monitor_$$ 2>/dev/null || echo 0)
FAILED=$(grep -c "FAILED" /tmp/hicat_monitor_$$ 2>/dev/null || echo 0)
RUNNING=$(grep -c "RUNNING" /tmp/hicat_monitor_$$ 2>/dev/null || echo 0)
PENDING=$(grep -c "PENDING" /tmp/hicat_monitor_$$ 2>/dev/null || echo 0)
rm -f /tmp/hicat_monitor_$$

echo ""
echo "Status: ${COMPLETED} done, ${RUNNING} running, ${PENDING} pending, ${FAILED} failed"
echo ""

# Flag failures for attention
if [[ "${FAILED}" -gt 0 ]]; then
    echo "⚠️  ${FAILED} jobs failed — check logs with:"
    echo "   grep FAIL /tmp/hicat_monitor_output  # (if saved)"
    echo "   ls /tscc/nfs/home/jhc103/cluster-logs/HiCAT_c*.err | xargs grep -l SIGKILL"
fi
