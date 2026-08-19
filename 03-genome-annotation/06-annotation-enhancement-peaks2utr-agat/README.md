# 06-annotation-enhancement-peaks2utr-agat

Finalize the merged gene annotation: annotate **3′ UTRs** from RNA-seq coverage
peaks with [peaks2utr](https://github.com/haessar/peaks2utr), re-key the UTRs
onto the merged gene models, and **normalize** the result with AGAT.

This is the last step of `03-genome-annotation`. It takes the merged annotation
from `05-annotation-merging/` (which carries gene models but no UTR features)
and produces the final annotation file.

## RNA-seq evidence (scMultiome)

3′ UTRs are inferred from 10x Genomics Single-Cell Multiome (scMultiome ATAC +
Gene Expression) libraries generated from degu brain tissue. Four Gene
Expression (GEX) BAMs were merged into a single `merged.bam`:

| Sample | Region |
|---|---|
| `181_PFC`  | prefrontal cortex |
| `181_dHIP` | dorsal hippocampus |
| `6997_dHIP`| dorsal hippocampus |
| `7000_PFC` | prefrontal cortex |

The same GEX data also provides the expression track `merged.bw` (a `bamCoverage`
BigWig of `merged.bam`) used by the "add novel genes" filter in step 05.

## Steps (in order)

| # | Script | What it does |
|---|--------|--------------|
| 1 | `01_peaks2utr_annotate_utrs.sh` | Full-annotation 3′ UTR run (the "previous" UTR set reused in step 4) |
| 2 | `02_rebuild_changed_subset.py` | Find genes whose CDS span changed vs. the previous annotation |
| 3 | `03_peaks2utr_changed_subset.sh` | Re-run peaks2utr on only the changed-gene subset |
| 4 | `04_assemble_final.py` | Combine gene models + UTRs (unchanged ← step 1, changed ← step 3) |
| 5 | `05_agat_normalize.sh` | AGAT normalization + final GFF3 integrity check |

Run order and dependencies:

```
01_peaks2utr_annotate_utrs.sh          (full 3' UTR run on the previous annotation)
        │
02_rebuild_changed_subset.py            (identify changed genes)
        │
03_peaks2utr_changed_subset.sh          (subset re-run)
        │
04_assemble_final.py                    (re-key UTRs onto merged models)
        │
05_agat_normalize.sh                    (final AGAT normalization)
```

The Python steps (`02`, `04`) read their input/output paths from `config.sh`
(via the repo-root walk-up), so run them after `source config.sh`, e.g.:

```bash
source config.sh
python3 02_rebuild_changed_subset.py
python3 04_assemble_final.py
```

## Key behaviours (non-obvious, worth knowing)

* **peaks2utr generates ONLY 3′ UTRs.** 5′ UTRs in its output are carried
  through from the input GFF, never created by peaks2utr itself.
* **peaks2utr's gffutils output is lossy** for inputs with orphaned CDS/exon
  (`no Parent`): it keeps only the outermost transcript per gene. So step 4
  matches *changed* genes by **gene ID**, not CDS span.
* **UTR `Parent` re-keying** (step 4) tries, in order: direct transcript-ID
  match → CDS-span match → single-transcript fallback → nearest 3′ end. This
  keeps UTRs parented even for de-novo genes (Braker transcript IDs, CDS rows
  with no `Parent`) and ncRNA (no CDS).
* **AGAT is the last step** (after peaks2utr), rebuilding
  `gene → mRNA → CDS/exon/intron/UTR` parentage and adding introns and
  start/stop codons.

## Requirements

* Conda environments: `toolshed-peaks2utr`, `toolshed-agat`,
  `toolshed-deeptools` (all names are set in `config.sh`).
* A `merged.bam` of stranded, coordinate-sorted scMultiome GEX alignments
  (`samtools merge` of the four samples above).
* The `src/` GFF3 helper package from `../05-annotation-merging/` (used by
  steps 02 and 04).

## Notes

* Scripts are Slurm batch scripts written for TSCC at UCSD; adjust the `#SBATCH`
  scheduler lines for your own cluster.
* `peaks2utr` resolves its `.cache`/`.log` relative to the current working
  directory, so steps 01 and 03 `cd` into the peaks2utr working directory first
  (so cached stranded BAMs, MACS peaks and pileups are reused).
