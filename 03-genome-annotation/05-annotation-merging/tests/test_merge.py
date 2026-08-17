"""Synthetic tests for the three annotation-merging rules.

Run from the repository root:

    python -m unittest discover -s tests -v
"""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.gffio import parse_gff, write_gff
from src.merge import MergeConfig, merge_annotations

LIFTOFF = """##gff-version 3
chrom1\tLiftoff\tgene\t100\t500\t.\t+\t.\tID=gene-Chmp3;Name=Chmp3;gene_biotype=protein_coding
chrom1\tLiftoff\tmRNA\t100\t500\t.\t+\t.\tID=rna-Chmp3-1;Parent=gene-Chmp3
chrom1\tLiftoff\tCDS\t150\t400\t.\t+\t0\tID=cds1;Parent=rna-Chmp3-1
chrom1\tLiftoff\tgene\t400\t700\t.\t+\t.\tID=gene-Rnf103;Name=Rnf103;gene_biotype=protein_coding
chrom1\tLiftoff\tmRNA\t400\t700\t.\t+\t.\tID=rna-Rnf103-1;Parent=gene-Rnf103
chrom1\tLiftoff\tCDS\t420\t680\t.\t+\t0\tID=cds2;Parent=rna-Rnf103-1
chrom1\tLiftoff\tgene\t1000\t1500\t.\t+\t.\tID=gene-LOC111111111;Name=LOC111111111;gene_biotype=protein_coding
chrom1\tLiftoff\tmRNA\t1000\t1500\t.\t+\t.\tID=rna-LOC1;Parent=gene-LOC111111111
chrom1\tLiftoff\tCDS\t1100\t1400\t.\t+\t0\tID=cds3;Parent=rna-LOC1
chrom1\tLiftoff\tgene\t2000\t2500\t.\t+\t.\tID=gene-Kept;Name=Kept;gene_biotype=protein_coding
chrom1\tLiftoff\tmRNA\t2000\t2500\t.\t+\t.\tID=rna-Kept-1;Parent=gene-Kept
chrom1\tLiftoff\tCDS\t2100\t2400\t.\t+\t0\tID=cds4;Parent=rna-Kept-1
chrom1\tLiftoff\tgene\t3000\t3300\t.\t+\t.\tID=gene-Discord;Name=Discord;gene_biotype=protein_coding
chrom1\tLiftoff\tmRNA\t3000\t3300\t.\t+\t.\tID=rna-Discord-1;Parent=gene-Discord
chrom1\tLiftoff\tCDS\t3050\t3250\t.\t+\t0\tID=cds5;Parent=rna-Discord-1
chrom1\tLiftoff\tgene\t6000\t6200\t.\t-\t.\tID=gene-LOC999999999;Name=LOC999999999;gene_biotype=protein_coding
chrom1\tLiftoff\tmRNA\t6000\t6200\t.\t-\t.\tID=rna-LOC2;Parent=gene-LOC999999999
chrom1\tLiftoff\tCDS\t6050\t6150\t.\t-\t0\tID=cds6;Parent=rna-LOC2
"""

DENOVO = """##gff-version 3
chrom1\tAUGUSTUS\tgene\t380\t720\t.\t+\t.\tID=g1;Name=RNF103
chrom1\tAUGUSTUS\tmRNA\t380\t720\t.\t+\t.\tID=g1.t1;Parent=g1
chrom1\tAUGUSTUS\tCDS\t400\t700\t.\t+\t0\tID=g1.t1.CDS1;Parent=g1.t1
chrom1\tAUGUSTUS\tgene\t950\t1550\t.\t+\t.\tID=g2;Name=NAMED1
chrom1\tAUGUSTUS\tmRNA\t950\t1550\t.\t+\t.\tID=g2.t1;Parent=g2
chrom1\tAUGUSTUS\tCDS\t1000\t1500\t.\t+\t0\tID=g2.t1.CDS1;Parent=g2.t1
chrom1\tAUGUSTUS\tgene\t3050\t3250\t.\t+\t.\tID=g3;Name=OTHER
chrom1\tAUGUSTUS\tmRNA\t3050\t3250\t.\t+\t.\tID=g3.t1;Parent=g3
chrom1\tAUGUSTUS\tCDS\t3070\t3230\t.\t+\t0\tID=g3.t1.CDS1;Parent=g3.t1
chrom1\tAUGUSTUS\tgene\t5000\t5300\t.\t+\t.\tID=g4;Name=NOVEL1
chrom1\tAUGUSTUS\tmRNA\t5000\t5300\t.\t+\t.\tID=g4.t1;Parent=g4
chrom1\tAUGUSTUS\tCDS\t5050\t5250\t.\t+\t0\tID=g4.t1.CDS1;Parent=g4.t1
"""


def _tmp(path: str, content: str) -> str:
    with open(path, "w") as fh:
        fh.write(content)
    return path


class TestMergeRules(unittest.TestCase):

    def setUp(self):
        self._dir = tempfile.mkdtemp()
        self.lo = _tmp(os.path.join(self._dir, "liftoff.gff"), LIFTOFF)
        self.dn = _tmp(os.path.join(self._dir, "denovo.gff"), DENOVO)

    def _names(self, genes):
        return {g.name.upper() for g in genes}

    def test_replace_and_add(self):
        result = merge_annotations(
            parse_gff(self.lo), parse_gff(self.dn), MergeConfig())
        names = self._names(result.final_genes)
        counts = result.counts
        sources = {g.name.upper(): g.source for g in result.final_genes}

        # Replace rule: discordant named Liftoff gene is KEPT ...
        self.assertIn("CHMP3", names)
        self.assertIn("DISCORD", names)
        # ... the concordant Liftoff gene (Rnf103) is gone and the de-novo
        #     RNF103 (AUGUSTUS source) took its place ...
        self.assertNotIn("LOC111111111", names)
        self.assertEqual(sources["RNF103"], "AUGUSTUS")
        # LOC gene replaced by a named de-novo gene
        self.assertIn("NAMED1", names)
        # de-novo gene overlapping only a discordant named gene is NOT placed
        self.assertNotIn("OTHER", names)
        # add rule: novel non-overlapping de-novo gene is added
        self.assertIn("NOVEL1", names)
        # no-overlap Liftoff genes (Kept, LOC999999999) stay
        self.assertIn("KEPT", names)
        self.assertIn("LOC999999999", names)

        self.assertEqual(counts["discordant_kept"], 2)      # CHMP3, DISCORD
        self.assertEqual(counts["replaced_concordant"], 1)  # Rnf103 <- RNF103
        self.assertEqual(counts["replaced_loc"], 1)         # LOC111111111 <- NAMED1
        self.assertEqual(counts["added"], 1)                # NOVEL1
        self.assertEqual(len(result.final_genes), 7)

    def test_name_change(self):
        result = merge_annotations(
            parse_gff(self.lo), parse_gff(self.dn),
            MergeConfig(rename_map={"LOC999999999": "REALGENE"}))
        names = self._names(result.final_genes)
        self.assertIn("REALGENE", names)
        self.assertNotIn("LOC999999999", names)
        self.assertEqual(result.counts.get("renamed", 0), 1)

    def test_roundtrip_write(self):
        result = merge_annotations(
            parse_gff(self.lo), parse_gff(self.dn), MergeConfig())
        out = os.path.join(self._dir, "merged.gff3")
        write_gff(result.final_genes, out)
        genes2 = parse_gff(out)
        self.assertEqual(len(genes2), len(result.final_genes))

    def test_list_mode_replace_fix(self):
        """List-driven merge keeps a discordant named Liftoff gene (the fix)."""
        from src.merge import merge_with_lists
        liftoff = parse_gff(self.lo)
        denovo = parse_gff(self.dn)
        # Chmp3 -> RNF103 is discordant (fix => keep); LOC111111111 -> NAMED1 is
        # a LOC replacement (allowed); Kept has no de-novo counterpart.
        result = merge_with_lists(
            liftoff, denovo, MergeConfig(),
            replace_list=[("gene-Chmp3", "RNF103"),
                          ("gene-LOC111111111", "NAMED1")],
            add_list={"g4": "NOVEL1"})
        names = {g.name.upper() for g in result.final_genes}
        self.assertIn("CHMP3", names)          # discordant named -> KEPT
        self.assertIn("RNF103", names)         # concordant member replaced in
        self.assertNotIn("LOC111111111", names)  # LOC -> replaced
        self.assertIn("NAMED1", names)
        self.assertIn("NOVEL1", names)         # add-list gene
        self.assertEqual(result.counts["replace_skipped_discordant_kept"], 1)

    def test_assign_names_evidence_priority(self):
        """assign_names follows the evidence-ranked naming convention."""
        from src.gffio import Gene
        from src.merge import assign_names

        def mk(source, name, start, evidence, conc=False):
            g = Gene(seqid="chr1", source=source, start=start, end=start + 1000,
                     strand="+", attributes=f"ID=gene-{name};Name={name}", name=name)
            g.evidence = evidence
            g.replace_concordant = conc
            return g

        # 1) replace-concordant de-novo gene is the plain parent; renamed paralog
        #    is renumbered sequentially (-l3 -> -L1, no gaps)
        genes = [mk("AUGUSTUS", "ALDH1A1", 1000, "replace", conc=True),
                 mk("Liftoff", "ALDH1A1-l3", 2000, "rename")]
        assign_names(genes)
        self.assertEqual({g.name.upper() for g in genes},
                         {"ALDH1A1", "ALDH1A1-L1"})

        # 2) no replace-concordant: renamed-plain is parent; de-novo -> -DL, replace -> -RL
        genes2 = [mk("Liftoff", "CLEC4A", 3000, "rename"),
                  mk("AUGUSTUS", "CLEC4A-dl1", 1000, "de_novo"),
                  mk("AUGUSTUS", "CLEC4A-rl1", 1500, "replace")]
        assign_names(genes2)
        self.assertEqual({g.name.upper() for g in genes2},
                         {"CLEC4A", "CLEC4A-DL1", "CLEC4A-RL1"})

        # 3) a single de-novo (novel) gene keeps its plain name
        genes3 = [mk("AUGUSTUS", "H2AB2", 5000, "de_novo")]
        assign_names(genes3)
        self.assertEqual(genes3[0].name.upper(), "H2AB2")

    def test_refseq_paralog_locus_preserved(self):
        """A RefSeq paralog locus tag (gene-Iqcn-2) keeps its original name."""
        from src.gffio import Gene
        from src.merge import assign_names
        g1 = Gene(seqid="chr3", source="Liftoff", start=105380351, end=105384824,
                  strand="+", attributes="ID=gene-Iqcn-2;Name=Iqcn", name="Iqcn",
                  gene_id="gene-Iqcn-2")
        g2 = Gene(seqid="chr3", source="Liftoff", start=105385591, end=105399646,
                  strand="+", attributes="ID=gene-Iqcn;Name=Iqcn", name="Iqcn",
                  gene_id="gene-Iqcn")
        assign_names([g1, g2])
        self.assertEqual(sorted(g.name for g in (g1, g2)), ["Iqcn", "Iqcn-2"])

    def test_parse_feature_before_gene(self):
        """Features preceding their gene row attach to that gene, not the previous."""
        from src.gffio import parse_gff
        gff = ("##gff-version 3\n"
               "chr10\tAUGUSTUS\tgene\t1000\t2000\t.\t+\t.\tID=g1;Name=A\n"
               "chr10\tAUGUSTUS\tmRNA\t1000\t2000\t.\t+\t.\tID=g1.t1;Parent=g1\n"
               "chr10\tAUGUSTUS\tstart_codon\t3000\t3002\t.\t+\t.\t.\n"
               "chr10\tAUGUSTUS\tmRNA\t3000\t4000\t.\t+\t.\t.\n"
               "chr10\tAUGUSTUS\tgene\t3000\t4000\t.\t+\t.\tID=g2;Name=B\n"
               "chr10\tAUGUSTUS\tmRNA\t3000\t4000\t.\t+\t.\tID=g2.t1;Parent=g2\n")
        import tempfile, os
        f = tempfile.NamedTemporaryFile("w", suffix=".gff", delete=False)
        f.write(gff); f.close()
        try:
            genes = parse_gff(f.name)
        finally:
            os.unlink(f.name)
        g1, g2 = genes
        self.assertEqual(g1.end, 2000)   # g1 must NOT absorb g2's pre-gene features
        self.assertEqual((g2.start, g2.end), (3000, 4000))

    def test_canonical_case(self):
        """Gene-name case matches the original annotation convention."""
        from src.merge import _canonical_name
        cases = {
            "RNF103": "Rnf103",
            "ALDH1A1": "Aldh1a1",
            "CLEC4A-RL1": "Clec4a-rl1",
            "DNAH7-RL2": "Dnah7-rl2",
            "MAGEA10-DL1": "Magea10-dl1",
            "ZNF709-L1": "Znf709-l1",
            "H2AB2": "H2ab2",
            "SLC22A13": "Slc22a13",
            "ND1": "ND1",            # mitochondrial: stays all-caps
            "COX1": "COX1",
            "C3": "C3",              # complement
            "F9": "F9",              # coagulation factor
            "LOC111814891": "LOC111814891",
            "Iqcn-2": "Iqcn-2",      # ref-paralog locus tag kept
        }
        for inp, exp in cases.items():
            self.assertEqual(_canonical_name(inp), exp, f"{inp} -> {_canonical_name(inp)}")

    def test_cli_end_to_end(self):
        """The CLI entry point runs on the example inputs and writes all outputs."""
        from merge_annotations import main
        exdir = os.path.join(self._dir, "cli_out")
        rc = main(["--liftOff", self.lo, "--denovo", self.dn,
                   "--outdir", exdir, "--prefix", "t"])
        self.assertEqual(rc, 0)
        for f in ("t_merged.gff3", "t_report.txt", "t_decisions.tsv",
                  "t_replaced_by_de_novo_concordant.tsv",
                  "t_discordant_kept.tsv", "t_added_novel.tsv"):
            self.assertTrue(os.path.exists(os.path.join(exdir, f)),
                            f"missing {f}")
        merged = parse_gff(os.path.join(exdir, "t_merged.gff3"))
        self.assertEqual(len(merged), 7)


if __name__ == "__main__":
    unittest.main()
