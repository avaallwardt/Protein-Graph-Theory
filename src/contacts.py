"""Build residue-contact graphs from protein structures.

A protein's interaction graph has one node per residue and an edge between two
residues whose Cα atoms are within a distance cutoff. This is the graph whose
treewidth governs the cost of exact Potts-model sequence optimization, so its
structure — not just the number of residues — determines tractability.

Needs Biopython and a local PDB file; kept separate from `graphs.py` so the
synthetic experiments run with no external data.
"""

from __future__ import annotations

import networkx as nx
import numpy as np

from .graphs import estimate_treewidth


def _ca_residues(pdb_path: str, chain_ids: list[str] | None = None):
    from Bio.PDB import PDBParser, is_aa

    model = PDBParser(QUIET=True).get_structure("protein", pdb_path)[0]
    residues = []
    for chain in model:
        if chain_ids is not None and chain.id not in chain_ids:
            continue
        for residue in chain:
            if is_aa(residue) and "CA" in residue:
                residues.append(residue)
    return residues


def contact_graph(
    pdb_path: str,
    chain_ids: list[str] | None = None,
    cutoff: float = 15.0,
    min_seq_separation: int = 0,
) -> nx.Graph:
    """Cα contact graph for a structure.

    An edge joins residues i, j when their Cα atoms are within ``cutoff`` Å.
    ``min_seq_separation`` keeps only edges with ``|i - j| >`` that value, which
    is how the long-range-contact experiments isolate non-local structure.
    """
    residues = _ca_residues(pdb_path, chain_ids)
    coords = np.array([r["CA"].coord for r in residues])
    graph = nx.Graph()
    graph.add_nodes_from(range(len(residues)))
    for i in range(len(residues)):
        for j in range(i + 1, len(residues)):
            if j - i <= min_seq_separation:
                continue
            if np.linalg.norm(coords[i] - coords[j]) < cutoff:
                graph.add_edge(i, j)
    return graph


def protein_stats(pdb_path: str, pdb_id: str, **kwargs) -> dict:
    """Summary stats + estimated treewidth for a structure's contact graph."""
    graph = contact_graph(pdb_path, **kwargs)
    n = graph.number_of_nodes()
    return {
        "protein": pdb_id,
        "vertices": n,
        "edges": graph.number_of_edges(),
        "avg_degree": sum(d for _, d in graph.degree()) / n if n else 0.0,
        "density": nx.density(graph),
        "treewidth": estimate_treewidth(graph),
    }
