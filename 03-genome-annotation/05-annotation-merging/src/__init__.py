"""annotation-merging: merge a Liftoff-transferred annotation with a de-novo (e.g. Braker) annotation.

Three enhancement rules (see README.md):
  1. name-change : rename mis-assigned Liftoff genes (optional mapping file)
  2. add         : add de-novo genes with a functional name whose CDS does not
                   overlap any retained Liftoff gene
  3. replace     : replace a Liftoff gene with a de-novo gene only when they are
                   the *same* gene (concordant name) or the Liftoff gene is an
                   unannotated LOC; keep named Liftoff genes whose name differs
                   (avoid dropping reference-derived annotations on paralog loci)
"""
