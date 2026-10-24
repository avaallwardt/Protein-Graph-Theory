"""Render results/protein_treewidth.png from the measured protein numbers.

These values come from running experiments/proteins.py on the downloaded PDB
structures (RCSB). They're transcribed here so the figure builds without the
large PDB files; re-run experiments/proteins.py to regenerate them.
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# 1GC1 chain G, 15 Å, cap on sequence separation |i-j|.  "full" = no cap.
SEP = [2, 4, 8, 12, 16, 24, 32, 50]
SEP_TW = [2, 4, 8, 12, 16, 16, 21, 48]
FULL_TW = 114

# 1GC1 chain G: Cα distance cutoff -> estimated treewidth.
CUTOFF = list(range(5, 16))
CUTOFF_TW = [7, 14, 34, 48, 51, 60, 76, 92, 94, 109, 114]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))

ax1.axhline(FULL_TW, ls="--", color="#C44E52", label=f"full graph (tw = {FULL_TW})")
ax1.plot(SEP, SEP_TW, marker="o", color="#4C72B0", label="capped graph")
ax1.set_xlabel("max sequence separation |i − j| kept")
ax1.set_ylabel("estimated treewidth")
ax1.set_title("Long-range contacts drive treewidth (1GC1)")
ax1.legend(frameon=False, fontsize=9)
ax1.grid(alpha=0.3)

ax2.plot(CUTOFF, CUTOFF_TW, marker="o", color="#4C72B0")
ax2.set_xlabel("Cα distance cutoff (Å)")
ax2.set_ylabel("estimated treewidth")
ax2.set_title("Contact radius vs. treewidth (1GC1)")
ax2.grid(alpha=0.3)

fig.tight_layout()
fig.savefig("results/protein_treewidth.png", dpi=150)
print("wrote results/protein_treewidth.png")
