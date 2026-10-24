"""Synthetic graph builders and a treewidth helper.

Treewidth measures how close a graph is to a tree. It sets the cost of exact
inference on a graphical model: solving a Potts / MAP problem exactly is
exponential in treewidth, so a low-treewidth interaction graph is tractable and
a high-treewidth one is not. This module builds reference graphs of a fixed size
and estimates their treewidth, so structure can be compared independent of size.

Treewidth is NP-hard to compute exactly; `estimate_treewidth` uses NetworkX's
min-fill-in heuristic, which gives an upper bound.
"""

from __future__ import annotations

import networkx as nx
from networkx.algorithms.approximation import treewidth_min_fill_in


def estimate_treewidth(graph: nx.Graph) -> int:
    """Upper-bound estimate of treewidth via the min-fill-in heuristic."""
    width, _ = treewidth_min_fill_in(graph)
    return width


def reference_graphs(n: int = 20) -> dict[str, nx.Graph]:
    """Four graphs on ~n nodes spanning low to high treewidth.

    Path (treewidth 1) and cycle (2) are tree-like; a grid is intermediate; the
    complete graph (n-1) is the worst case. Same node budget, very different
    structure — and very different treewidth.
    """
    rows, cols = 4, max(1, n // 4)
    return {
        "Path": nx.path_graph(n),
        "Cycle": nx.cycle_graph(n),
        "Grid": nx.grid_2d_graph(rows, cols),
        "Clique": nx.complete_graph(n),
    }


def path_with_random_edges(n: int, extra_edges: int, rng) -> nx.Graph:
    """A path on n nodes plus `extra_edges` random chords.

    Models adding long-range contacts to an otherwise linear chain. `rng` is a
    `random.Random` so the draw is reproducible.
    """
    graph = nx.path_graph(n)
    candidates = [
        (i, j)
        for i in range(n)
        for j in range(i + 1, n)
        if not graph.has_edge(i, j)
    ]
    graph.add_edges_from(rng.sample(candidates, extra_edges))
    return graph


def graph_stats(graph: nx.Graph) -> dict[str, float]:
    """Vertices, edges, average degree, density, and estimated treewidth."""
    n = graph.number_of_nodes()
    avg_degree = sum(d for _, d in graph.degree()) / n if n else 0.0
    return {
        "vertices": n,
        "edges": graph.number_of_edges(),
        "avg_degree": avg_degree,
        "density": nx.density(graph),
        "treewidth": estimate_treewidth(graph),
    }
