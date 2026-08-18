#!/usr/bin/env python3
"""Execute pyCircularize_features_noMethylation.ipynb cells top-to-bottom with
the current interpreter (python-visualizations env). Handles the `%%time`
magic by stripping it. Chdirs to the notebook's directory so relative data
paths resolve. Stops on the first cell error and reports which cell failed."""
import os
import re
import sys
import traceback

import nbformat

NB = "pyCircularize_features_noMethylation.ipynb"
nb_dir = os.path.dirname(os.path.abspath(NB)) or "."
os.chdir(nb_dir)

nb = nbformat.read(NB, as_version=4)

MAGIC = re.compile(r"^\s*%%\w+")

ns = {"__name__": "__main__", "__file__": os.path.abspath(NB)}

for i, cell in enumerate(nb.cells):
    if cell.cell_type != "code":
        continue
    src = cell.source
    # strip Jupyter cell magics (e.g. %%time)
    src = "\n".join(
        line for line in src.splitlines() if not MAGIC.match(line)
    )
    print(f"=== executing code cell {i} ({len(src.splitlines())} lines) ===",
          flush=True)
    try:
        exec(compile(src, f"<cell {i}>", "exec"), ns)
    except Exception:
        print(f"!!! ERROR in cell {i} !!!", flush=True)
        traceback.print_exc()
        sys.exit(1)

print("ALL CELLS EXECUTED OK")
