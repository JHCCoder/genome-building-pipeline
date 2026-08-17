#!/usr/bin/env python3
"""Build .ipynb notebooks from the `# %%`-formatted .py scripts in this directory.

Each .py file uses the VSCode/Jupyter cell convention:

    # %% [markdown]     -> markdown cell (lines prefixed with ``# ``)
    # %%               -> code cell

Run from this directory:

    python build_notebooks.py

This regenerates the three .ipynb files from the three .py files (the .py files
are the source of truth, so they can be diffed and run headlessly).
"""
import json
import re
from pathlib import Path

HERE = Path(__file__).parent
SCRIPTS = ["01_add_novel_genes.py", "02_rename_loc_genes.py", "03_replace_loc_genes.py"]

MARKDOWN = re.compile(r"^#\s*%%\s*\[markdown\]\s*$")
CODE = re.compile(r"^#\s*%%\s*$")


def parse(src: str):
    cells = []
    lines = src.splitlines()
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        if MARKDOWN.match(line):
            i += 1
            buf = []
            while i < n and not CODE.match(lines[i]) and not MARKDOWN.match(lines[i]):
                raw = lines[i]
                if raw.startswith("# "):
                    buf.append(raw[2:])
                elif raw == "#":
                    buf.append("")
                else:
                    buf.append(raw)
                i += 1
            cells.append({"cell_type": "markdown", "source": buf})
        elif CODE.match(line):
            i += 1
            buf = []
            while i < n and not CODE.match(lines[i]) and not MARKDOWN.match(lines[i]):
                buf.append(lines[i])
                i += 1
            cells.append({"cell_type": "code", "source": buf})
        else:
            i += 1
    return cells


def to_notebook(cells):
    nb_cells = []
    for c in cells:
        if c["cell_type"] == "markdown":
            nb_cells.append({
                "cell_type": "markdown",
                "metadata": {},
                "source": [ln + "\n" for ln in c["source"]],
            })
        else:
            nb_cells.append({
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [ln + "\n" for ln in c["source"]],
            })
    return {
        "cells": nb_cells,
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python", "version": "3"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


def main():
    for name in SCRIPTS:
        py = HERE / name
        src = py.read_text()
        cells = parse(src)
        nb = to_notebook(cells)
        out = HERE / (py.stem + ".ipynb")
        out.write_text(json.dumps(nb, indent=1))
        n_md = sum(1 for c in nb["cells"] if c["cell_type"] == "markdown")
        n_code = sum(1 for c in nb["cells"] if c["cell_type"] == "code")
        print(f"{name} -> {out.name}  ({n_md} markdown, {n_code} code cells)")


if __name__ == "__main__":
    main()
