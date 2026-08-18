# Inputs — edit for your environment (mirror config.sh at the repo root).
import pandas as pd, numpy as np, re
from intervaltree import IntervalTree
from collections import Counter, defaultdict

DUP = "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/output/outputs-from-purge-duplicate/hifiasm-041425-scaffolded-assembly/dups.bed"
GFF = "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/code/command-line-script/annotation-merging/output/hifiasm-041425-denovoEnhanced_peaks2utr_sorted.agat.gff3"

CATS = ['HAPLOTIG','HIGHCOV','JUNK','OVLP','REPEAT']

# dups.bed -> per-category interval trees
dup = pd.read_csv(DUP, sep='\t', header=None, names=['chrom','start','end','cat','extra'])
trees = {c: {} for c in CATS}
for _, r in dup.iterrows():
    if r['cat'] in trees:
        trees[r['cat']].setdefault(r['chrom'], IntervalTree())[int(r['start']):int(r['end'])] = True

# new GFF gene features + Name
gff = pd.read_csv(GFF, sep='\t', comment='#', header=None,
                  names=['chrom','src','type','gs','ge','sc','str','ph','attr'])
gff = gff[gff['type']=='gene'].copy()
gff['name'] = gff['attr'].str.extract(r'Name=([^;]+)', expand=False)

# classify each gene by Name suffix
SUF = re.compile(r'-(dl|l|rl)(\d+)$')
def base(name):
    m = SUF.search(name)
    return name[:m.start()] if m else name

names = gff['name'].dropna().tolist()
base_set = set(base(n) for n in names)
# base names that have at least one paralog copy
bases_with_copies = set(base(n) for n in names if SUF.search(n))
parent_set = set(bases_with_copies)  # any base that has copies is a "family parent"

def gtype(name):
    m = SUF.search(name)
    if m:
        return {'dl':'paralog_dl','l':'paralog_l','rl':'paralog_rl'}[m.group(1)]
    return 'parent' if name in bases_with_copies else 'singleton'

totals = Counter()
overlap = defaultdict(Counter)
for _, r in gff.iterrows():
    nm = r['name']
    if not isinstance(nm, str):
        continue
    gt = gtype(nm)
    totals[gt] += 1
    chrom = r['chrom']
    t_gs, t_ge = int(r['gs']), int(r['ge'])
    for cat in CATS:
        t = trees[cat].get(chrom)
        if t is not None and t.overlap(t_gs, t_ge):
            overlap[gt][cat] += 1

print("TOTAL_GENES =", len(gff))
print("\nTotals:")
for gt in ['singleton','parent','paralog_l','paralog_dl','paralog_rl']:
    print(f"  {gt}: {totals[gt]}")
print("\nOverlap counts (gene_type -> {cat: n}):")
for gt in ['singleton','parent','paralog_l','paralog_dl','paralog_rl']:
    print(f"  {gt}: " + repr({c: overlap[gt].get(c,0) for c in CATS}))
