# annotation-merging

Merge a **Liftoff**-transferred (reference-based) gene annotation with a
**de-novo** (e.g. Braker/AUGUSTUS) gene annotation into a single final
annotation, using three explicit enhancement rules.

This is the merge step used for the *Octodon degus* (degu) genome annotation
and is released here as a standalone, dependency-light tool.

---

## Why this merge?

A Liftoff annotation is a high-quality, reference-derived gene model transfer,
but it cannot discover genes absent from the reference. A de-novo annotation
(Braker + RNA-seq) can, but its gene structures are more error-prone and its
functional names are only as good as the BLAST/InterProScan evidence behind
them.

The merge keeps the reference-derived annotation by default and only
incorporates de-novo genes under explicit, conservative conditions.

## Two usage modes

**List mode (recommended, reproduces the original curated logic).** The
name-change / add / replace decisions are supplied as curated TSV lists —
exactly how the degu annotation was built. This keeps the analysis decisions
(which de-novo genes are expressed and well-supported) explicit and auditable.

**Automatic mode.** If no lists are provided, the three rules are derived
directly from the two annotations (a de-novo gene is considered when it has a
functional name, and the replace rule applies the concordance fix below). This
is a convenient fallback but does not reproduce the curated expression /
length filters of the original pipeline.

## The three rules

### 1. name-change

Rename Liftoff genes from a mapping (e.g. correcting an unannotated `LOC` name
based on a comparative analysis). The Liftoff gene structure is kept; only
`Name=`/`ID=` are updated.

*Input (list mode):* a two-column TSV (`old_name<TAB>new_name`), via
`--rename-map`.

### 2. replace — with the *concordance fix*

A Liftoff gene is replaced by a de-novo gene only when the Liftoff gene is an
unannotated `LOC`, or when the new name is **concordant** (case/suffix-insensitive)
with the Liftoff gene's own name:

| Liftoff gene | Decision | Reason |
|---|---|---|
| `LOC##########` | **replaced** | the de-novo gene supplies a real name |
| named, **same** name as the replacing de-novo gene | **replaced** | concordant — same locus |
| named, **different** name (e.g. a paralog) | **kept** | discordant — reference annotation preserved |

> **The concordance fix.** In the degu annotation, replacing a named Liftoff
> gene with a de-novo gene of a *different* name silently removed reference
> genes at paralog loci (CHMP3→RNF103, CLEC12A→CLEC4A, DNAH12→DNAH7,
> EIF4EBP3→ANKHD1, TUBG1→TUBG2). The fix keeps such discordant named genes —
> they carry the reference exon structure and should not be dropped for a
> differently-named de-novo prediction.

*Input (list mode):* a two-column TSV (`liftoff_gene<TAB>new_name`), via
`--replace-list`. The tool resolves the de-novo gene by name+overlap and applies
the fix (a discordant named Liftoff gene is kept, never removed).

### 3. add (novel genes)

A de-novo gene is added as a **novel** gene only when it is in the curated
`--add-list` (list mode), or (automatic mode) it has a functional name whose
CDS does not overlap any retained Liftoff gene's CDS and passes optional
expression / length filters.

### 4. gene naming (evidence-ranked, post-processing)

The final annotation is named so that every gene name is unique, following an
evidence-ranked convention.  Within a gene family (same suffix-stripped base
name), the gene that keeps the **plain parent name** is chosen by evidence
priority:

1. a **replace**-introduced de-novo gene whose replacement was *concordant*
   (it took over the reference locus that carried the same name);
2. a gene **renamed** by the name-change rule that holds the plain name;
3. a kept **Liftoff** copy with the plain name;
4. a **de-novo** (added) gene.

Every other copy is a paralog and is suffixed according to its evidence class —
`-RL` for replace copies, `-L` for renamed paralogs, `-DL` for de-novo copies,
and `-L` for kept Liftoff paralogs — **renumbered sequentially** within the
family by genomic order (a lone paralog is always `-L1`/`-DL1`/`-RL1`, never a
gapped number like `-l3`).

---

## Inputs & outputs

### Inputs
* **Liftoff GFF3** (`--liftOff`) — the reference-transferred annotation, e.g.
  from `liftoff`.
* **De-novo GFF3** (`--denovo`) — the de-novo annotation with functional gene
  names in the `Name=` attribute (e.g. Braker genes annotated by
  BLAST/InterProScan).

### Outputs (written to `--outdir` with `--prefix`)
| File | Contents |
|---|---|
| `<prefix>_merged.gff3` | final merged annotation (sorted) |
| `<prefix>_report.txt` | summary counts per rule |
| `<prefix>_decisions.tsv` | per-gene decision log (category, genes, detail) |
| `<prefix>_replaced_by_de_novo_LOC_rule.tsv` | Liftoff LOC genes replaced |
| `<prefix>_replaced_by_de_novo_concordant.tsv` | concordant named genes replaced |
| `<prefix>_discordant_kept.tsv` | named Liftoff genes kept (name differs) |
| `<prefix>_added_novel.tsv` | de-novo genes added as novel |
| `<prefix>_renamed.tsv` | genes renamed by the name-change rule |

## Usage

**List mode** (curated name-change / add / replace lists):

```bash
python merge_annotations.py \
    --liftOff       liftoff.gff \
    --denovo        braker.functional.gff \
    --rename-map    name_change.tsv     # old_name<TAB>new_name
    --add-list      add_genes.tsv       # de_novo_gene_id<TAB>new_name
    --replace-list  replace_genes.tsv   # liftoff_gene<TAB>new_name
    --outdir        output \
    --prefix        merged
```

**Automatic mode** (derive the three rules from the two annotations):

```bash
python merge_annotations.py \
    --liftOff liftoff.gff \
    --denovo  braker.functional.gff \
    --outdir  output \
    --prefix  merged
```

Run the unit tests:

```bash
python -m unittest discover -s tests -v
```

### Options

| Option | Default | Meaning |
|---|---|---|
| `--liftOff` | — | Liftoff GFF3 (required) |
| `--denovo` | — | de-novo GFF3 with `Name=` (required) |
| `--outdir` | `output` | output directory |
| `--prefix` | `merged` | output filename prefix |
| `--rename-map` | — | TSV `old_name<TAB>new_name` for rule 1 |
| `--require-name` / `--no-require-name` | on | require a functional name on de-novo genes |
| `--min-cds-bp` | 0 | min de-novo CDS length (bp) to consider |
| `--overlap-fraction` | 0.0 | min gene-span overlap fraction to trigger rule 2 |
| `--cds-overlap-fraction` | 0.0 | max CDS overlap fraction allowed for a novel add |
| `--no-skip-dup-base-name` | (dup check on) | allow adding a de-novo gene whose base name exists |
| `--bigwig` | — | optional BigWig expression track (requires `pyBigWig`) |
| `--min-expression` | 0.0 | min mean signal for the add rule |
| `--bigwig-chrom-style` | `ucsc` | `ucsc` (`chr1`) or `ensembl` (`1`) |

## Dependencies

Python ≥ 3.9, standard library only. `pyBigWig` is required **only** if an
expression track is supplied (`--bigwig`).

## Example

`example/` contains a tiny synthetic Liftoff and de-novo annotation that
exercise all three rules:

```bash
python merge_annotations.py \
    --liftOff example/liftoff.example.gff \
    --denovo  example/denovo.example.gff \
    --outdir  example/output \
    --prefix  example
```

Expected: `Chmp3` and `Discord` are kept (discordant names), `Rnf103` is
replaced by de-novo `RNF103` (concordant), `LOC111111111` is replaced by
`NAMED1`, `NOVEL1` is added, and de-novo `OTHER` is not placed (it overlaps only
a discordant named gene).

## Reproducing the degu (*Octodon degus*) annotation

The degu annotation was built with the list mode, using curated
name-change / add / replace lists derived from RNA-seq expression, CDS-length
and functional-annotation evidence:

```bash
python merge_annotations.py \
    --liftOff       hifiasm-041425-scaffolded-chrAssigned-mito.gff \
    --denovo        braker_peak2utr.gff3 \
    --rename-map    gene_LOC_nameChange_unique.tsv \
    --add-list      gene_gBraker_toAdd_unique.tsv \
    --replace-list  gene_LOC_replace_unique.tsv \
    --outdir        output \
    --prefix        hifiasm-041425-denovoEnhanced
```

Notes for the degu run:
* The de-novo input carries functional names assigned by BLAST (UniProt/NR) +
  InterProScan in the `Name=` attribute.
* The curated replace list still contains the five discordant entries
  (`Chmp3→RNF103`, `Clec12a→CLEC4A`, `Dnah12→DNAH7`, `Eif4ebp3→ANKHD1`,
  `Tubg1→TUBG2`); the tool's concordance fix keeps those Liftoff genes and
  records them under `replace_skipped_discordant_kept`.
* Post-merge formatting (sorting with AGAT, 3′ UTR extension with `peaks2utr`,
  paralog `-L`/`-DL`/`-RL` suffix naming) are separate, degu-specific steps and
  are not part of `merge_annotations.py` itself. The UTR step is, however,
  captured in this directory — see below.

## Finalization: 3′ UTR annotation (`peaks2utr`) + assembly

The merge output carries gene models but no UTR features. UTRs are restored in
two steps, both in this directory:

### 1. `peaks2utr_subset_rerun.sh` — 3′ UTR annotation with peaks2utr

[peaks2utr](https://github.com/haessar/peaks2utr) annotates 3′ UTRs from RNA-seq
coverage peaks. **It generates only 3′ UTRs** — 5′ UTRs are carried through from
its input GFF, never created. Run it inside the `peaks2utr` working directory
(it resolves `.cache`/`.log` relative to `cwd`) so the cached stranded BAMs,
MACS peaks and pileups are reused:

```bash
peaks2utr <input.gff3> merged.bam \
    --do-pseudo --keep-cache --extend-utr -p 8 --max-distance 1500 \
    -o <output.gff3>
```

### 2. `assemble_final.py` — combine gene models + UTRs

Takes the merged gene models and re-keys UTR features onto them from a previous
peaks2utr output (unchanged genes) plus a fresh subset peaks2utr output
(changed genes). Two non-obvious behaviours are baked in (see the module
docstring):

* **changed genes match by gene ID, not CDS span** — peaks2utr's gffutils output
  keeps only the outermost transcript per gene, so gene-level CDS spans disagree
  with the merged (multi-transcript) gene;
* **UTR `Parent` re-keying** tries, in order: direct transcript-ID match →
  CDS-span match → single-transcript fallback → nearest-3′-end. This keeps UTRs
  parented even for de-novo genes (Braker transcript IDs, CDS rows with no
  `Parent`) and ncRNA (no CDS), so the final GFF has no dangling UTR Parents.

## Share-repo notes

This tool lives at `03-genome-annotation/05-annotation-merging/` in the
`genome-building-pipeline` repository.

**Finalization workflow** (degus-specific, run in this order):

1. `list-building/` — build the three curated TSV lists (rename / replace / add).
2. `merge_annotations.py` (list mode) — `output/<prefix>_merged.gff3`.
3. `rebuild_subset.py` — the "changed-gene" subset for peaks2utr.
4. `peaks2utr_subset_rerun.sh` — 3′ UTRs for the changed genes.
5. `assemble_final.py` — combine gene models + UTRs.
6. `agat_normalize_final.sh` (or `finalize_annotation.sh`) — AGAT normalize + verify.

The degus-specific file paths were moved out of the code so the scripts are
reusable:

- `assemble_final.py` and `rebuild_subset.py` take `--*` CLI arguments (run with
  `--help` to see them);
- the three `.sh` orchestration scripts read conda env names and the peaks2utr
  working directory from `config.sh` (repo root), and keep their degus-specific
  input filenames as clearly-marked variables at the top of each script.

## License / citation

Released as part of the *Octodon degus* genome project. Please cite the
accompanying manuscript when using this tool.
