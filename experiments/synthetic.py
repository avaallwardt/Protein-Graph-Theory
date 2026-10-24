"""Synthetic treewidth experiments. No external data required.

Run from the repo root:  PYTHONPATH=. python experiments/synthetic.py

Produces:
  results/synthetic_reference.png   treewidth of path/cycle/grid/clique at n=20
  results/random_edges.png          treewidth vs. number of random chords added
and prints the underlying tables.
"""

import random

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from src.graphs import estimate_treewidth, path_with_random_edges, reference_graphs

N = 20


def reference():
    graphs = reference_graphs(N)
    names = list(graphs)
    tws = [estimate_treewidth(g) for g in graphs.values()]
    print("Reference graphs (n = 20):")
    for name, g, t in zip(names, graphs.values(), tws):
        print(f"  {name:6} edges={g.number_of_edges():3d}  treewidth={t}")

    fig, ax = plt.subplots(figsize=(5, 3.6))
    ax.bar(names, tws, color="#4C72B0")
    ax.set_ylabel("estimated treewidth")
    ax.set_title("Same size (n=20), different structure")
    for i, t in enumerate(tws):
        ax.text(i, t + 0.3, str(t), ha="center", fontsize=9)
    fig.tight_layout()
    fig.savefig("results/synthetic_reference.png", dpi=150)
    print("  wrote results/synthetic_reference.png")


def random_edges(trials=20):
    rng = random.Random(42)
    counts = list(range(0, 21, 2))
    means, sds = [], []
    for ne in counts:
        ws = [estimate_treewidth(path_with_random_edges(N, ne, rng)) for _ in range(trials)]
        means.append(np.mean(ws))
        sds.append(np.std(ws))
    print("\nRandom chords added to a path (n=20, %d trials each):" % trials)
    for ne, m, s in zip(counts, means, sds):
        print(f"  +{ne:2d} edges  mean treewidth={m:.2f}  sd={s:.2f}")

    fig, ax = plt.subplots(figsize=(5.4, 3.8))
    ax.errorbar(counts, means, yerr=sds, marker="o", capsize=4, color="#4C72B0")
    ax.set_xlabel("number of random edges added to a path")
    ax.set_ylabel("mean estimated treewidth ± 1 SD")
    ax.set_title("Connectivity drives treewidth at fixed size")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig("results/random_edges.png", dpi=150)
    print("  wrote results/random_edges.png")


if __name__ == "__main__":
    reference()
    random_edges()
