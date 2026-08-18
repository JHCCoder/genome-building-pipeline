# Contamination removal

Removal of contaminant sequences from the degu assembly. Contamination was
screened with **NCBI FCS-GX** (genome contaminants) and **FCS-adapter**
(adapter contamination), both run through the NCBI Galaxy portal, then each
flagged contig was traced through the scaffolded assembly to confirm none
landed in the final chromosomes.

## Workflow

1. Run FCS-GX and FCS-adapter on the assembly via the
   [NCBI Galaxy](https://galaxy.ncbi.nlm.nih.gov/) portal:
   - FCS-GX genome report → `contamination_action.txt` (`EXCLUDE`/`REVIEW` per contig).
   - FCS-adapter report → `adaptor_report.txt` (`ACTION_TRIM` per adapter hit).
2. `check_contam.ipynb` — verify none of the `EXCLUDE` contigs ended up in the
   30 chromosome scaffolds (traces contigs via the HapHiC AGP file).
3. `check_adapter.ipynb` — verify the adapter hit (`ptg000006l`) is a false
   positive (it lies inside exon 4 of the *Irf4* gene).
4. `contamination-timing-verification.txt` — full write-up of the findings and
   the early-vs-late removal comparison.

## Files

| File | Description |
|---|---|
| `check_contam.ipynb` | Verify no contaminant contig is in the 30 chromosomes |
| `check_adapter.ipynb` | Verify the adapter hit is a false positive |
| `contamination_action.txt` | FCS-GX genome report (EXCLUDE/REVIEW list) |
| `adaptor_report.txt` | FCS-adapter report (the single adapter hit) |
| `contamination-timing-verification.txt` | Findings + timing comparison |

## Key result

FCS-GX flagged 90 contigs (89 `EXCLUDE` + 1 `REVIEW`); **none** were
incorporated into the 30 chromosome scaffolds. The 88 prokaryotic contigs
removed in the final assembly are exactly the 88 from the early screen.

## Notes

- FCS-GX / FCS-adapter are web tools (Galaxy), not command-line scripts — the
  reports above are their inputs, and the notebooks verify the downstream
  assembly. See `contamination-timing-verification.txt` for the full reasoning.
- The notebooks read their input paths from a clearly-marked configuration cell
  at the top (mirrors `config.sh` at the repo root).
