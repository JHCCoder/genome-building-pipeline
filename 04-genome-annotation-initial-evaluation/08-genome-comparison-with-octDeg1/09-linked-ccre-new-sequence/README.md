# 09-linked-ccre-new-sequence

Which peak–gene linked cCREs from the OctDeg2.0 scMultiome data lie in sequence that
OctDeg1.0 does not contain, and which TF motifs are enriched in them. Planned as the
replacement for the "median genes per cell" half of Figure 3D (genes per cell barely
differs once both assemblies use matched gene models; see `../README.md`).

Inputs come from the matched-annotation Cell Ranger ARC reruns of the parent folder
(`../07_cellranger_arc_count_octDeg2.sh`, Liftoff-only OctDeg2.0 reference) and the
OctDeg1.0 ↔ OctDeg2.0 minimap2 alignment of Figure 3A (`$FIG3_ALIGNED_BED`). Large
outputs are written to `$FIG3D_CCRE_DIR`; small summary tables to `results/`.

## Steps

| # | Script | What it does | Runs |
|---|--------|--------------|------|
| 1 | `01_classify_linked_peaks.sh` | Collect the `peak-gene` linked peaks of the four samples (peak–peak links excluded), merge them into a union set, and flag peaks < 50% covered by the genome alignment as candidates | login node, < 1 min |
| 2 | `02_blast_candidates_octdeg1.sh` | BLASTN (megablast, no DUST, ≥ 90% identity) of each candidate's full sequence against OctDeg1.0, because minimap2 misses collapsed and high-copy repeats | Slurm, 8 cores |
| 3 | `03_classify_new_sequence.py` | Classify each union peak: `shared`, `resolved_repeat` (≥ 2 OctDeg1.0 loci), `partial`, or `new_sequence` (no OctDeg1.0 hit) | login node |
| 4 | `04_homer_known_motifs.sh` | HOMER known motifs (`-size 200 -nomotif`, no repeat masking) in all linked peaks not shared with OctDeg1.0 (new + partial + resolved repeat, 2,315) vs shared linked peaks (98,860) as GC-matched background | Slurm, 8 cores |

## Definitions

* **Linked peak:** a peak in a `peak-gene` row of Cell Ranger ARC `feature_linkage.bedpe`.
  These are about 10% of Cell Ranger's "feature linkages detected"; the rest are
  peak–peak links.
* **OctDeg1.0 locus (BLAST):** HSPs with ≥ 90% identity that together cover ≥ 80% of the
  peak, on one OctDeg1.0 sequence within 5 kb.
* Peaks are 220–1,960 bp wide (median ~870 bp); the BLAST check uses the full peak.

## Results so far

Step 1 (`results/01_alignment_overlap_summary.tsv`): 101,175 union linked peaks, of which
3,037 (3.0%) are < 50% covered by the genome alignment (1,488 with no aligned base);
per sample 3.1–3.8%.

Step 2–3 (`results/03_classification_summary.tsv`, BLAST job 12812647, 14 min): of the 3,037
candidates, 722 have exactly one full-length OctDeg1.0 locus (missed by the genome
alignment, reclassified as `shared`), 39 occur at ≥ 2 loci (`resolved_repeat`), 1,940 have
only partial hits (`partial`; repeat-rich, 5.7 M HSPs in total), and 336 have no hit
(`new_sequence`, 0.33% of 101,175 union linked peaks; 0.31–0.35% per sample).
Overlap with the genome alignment alone would have called all 3,037 "new".

Step 4 (`results/04_homer_knownResults_notShared_vs_shared.txt`, HOMER job 12812855, 6.5 min):
2,315 not-shared linked peaks vs 98,860 shared linked peaks (GC-matched, 200 bp centred).
**No known motif is enriched.** With HOMER's full library (1,049 motifs, including plant
motifs, because the genome is a custom FASTA) the best q is 0.051 (a plant motif); restricted
to the 473 vertebrate motifs (BH recomputed from the same p-values) the best are Srebp2 and
ISRE (p ≈ 1.2–1.4 × 10⁻³, q = 0.32; 1.7% vs 1.0% and 0.95% vs 0.46% of peaks), and no motif
reaches q ≤ 0.1. The linked peaks in OctDeg2.0-specific sequence have the same motif
composition as linked peaks shared with OctDeg1.0. A re-run should pass
`-mknown <homer>/data/knownTFs/vertebrates/known.motifs` to test vertebrate motifs only.
