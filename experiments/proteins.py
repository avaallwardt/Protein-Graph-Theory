"""Protein contact-graph treewidth experiments. Needs local PDB files.

Download the structures first (RCSB), then run from the repo root:
    PYTHONPATH=. python experiments/proteins.py path/to/pdbs

Reproduces three things from the project:
  1. treewidth vs. size/degree across proteins of different folds,
  2. the Cα cutoff sweep on one structure (contact radius -> treewidth),
  3. the long-range-contact experiment (cap |i-j| -> treewidth collapses).

PDB IDs used in the writeup: 1GC1, 1UBQ, 1CRN, 1AKE, 1TIM.
"""

import sys

import pandas as pd

from src.contacts import contact_graph, protein_stats
from src.graphs import estimate_treewidth


def survey(pdb_dir: str):
    jobs = [("1UBQ", "A"), ("1CRN", "A"), ("1AKE", "A"), ("1TIM", "A"), ("1GC1", "G")]
    rows = []
    for pdb_id, chain in jobs:
        path = f"{pdb_dir}/{pdb_id}.pdb"
        rows.append(protein_stats(path, pdb_id, chain_ids=[chain], cutoff=15.0))
    print(pd.DataFrame(rows).to_string(index=False))


def cutoff_sweep(pdb_dir: str, pdb_id="1GC1", chain="G"):
    print(f"\nCutoff sweep on {pdb_id} chain {chain}:")
    rows = []
    for cutoff in range(5, 16):
        g = contact_graph(f"{pdb_dir}/{pdb_id}.pdb", chain_ids=[chain], cutoff=float(cutoff))
        deg = sum(d for _, d in g.degree()) / g.number_of_nodes()
        rows.append({"cutoff_A": cutoff, "edges": g.number_of_edges(),
                     "avg_degree": round(deg, 1), "treewidth": estimate_treewidth(g)})
    df = pd.DataFrame(rows)
    print(df.to_string(index=False))
    print("corr(avg_degree, treewidth) =",
          round(df["avg_degree"].corr(df["treewidth"]), 3))


def long_range(pdb_dir: str, pdb_id="1GC1", chain="G"):
    print(f"\nLong-range-contact experiment on {pdb_id} chain {chain} (15 A):")
    for sep in [2, 4, 8, 12, 16, 24, 32, 50, 10**9]:
        g = contact_graph(f"{pdb_dir}/{pdb_id}.pdb", chain_ids=[chain],
                          cutoff=15.0, min_seq_separation=0)
        # keep only edges within the sequence-separation cap
        capped = g.copy()
        capped.remove_edges_from([(i, j) for i, j in g.edges() if abs(i - j) > sep])
        label = "full" if sep > 10**8 else f"|i-j|<={sep}"
        print(f"  {label:12} edges={capped.number_of_edges():5d}  treewidth={estimate_treewidth(capped)}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: PYTHONPATH=. python experiments/proteins.py <dir with PDB files>")
    d = sys.argv[1]
    survey(d)
    cutoff_sweep(d)
    long_range(d)
