import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys
import os

group         = sys.argv[1]
admixture_out = sys.argv[2]
plot_out      = sys.argv[3]
k_min         = int(sys.argv[4])
k_max         = int(sys.argv[5])

os.makedirs(plot_out, exist_ok=True)

# --- Load sample order from .fam file ---
fam = pd.read_csv(f"{admixture_out}/{group}_pruned_remapped.fam",
                  sep=" ", header=None)
samples = fam[1].values

# --- Parse CV errors ---
cv = {}
with open(f"{admixture_out}/cv_errors.txt") as f:
    for line in f:
        parts = line.strip().split()
        k_val  = int(parts[2].replace("(K=", "").replace("):", ""))
        cv_val = float(parts[-1])
        cv[k_val] = cv_val

best_k = min(cv, key=cv.get)
print(f"{group}: Best K = {best_k} (CV error = {cv[best_k]:.4f})")

# --- Build figure: CV error plot + one bar per K ---
n_plots = (k_max - k_min + 1) + 1  # one per K + CV error
fig = plt.figure(figsize=(14, n_plots * 1.5))
gs  = gridspec.GridSpec(n_plots, 1, hspace=0.6)

# CV error plot
ax_cv = fig.add_subplot(gs[0])
ax_cv.plot(list(cv.keys()), list(cv.values()), marker="o", color="black")
ax_cv.axvline(best_k, color="red", linestyle="--", label=f"Best K={best_k}")
ax_cv.set_xlabel("K")
ax_cv.set_ylabel("CV error")
ax_cv.set_title(f"{group} — CV error")
ax_cv.set_xticks(range(k_min, k_max + 1))
ax_cv.legend()

# One stacked bar chart per K
colors = plt.cm.tab10.colors

for idx, k in enumerate(range(k_min, k_max + 1)):
    q_file = f"{admixture_out}/{group}_pruned_remapped.{k}.Q"
    if not os.path.exists(q_file):
        print(f"  Warning: {q_file} not found, skipping K={k}")
        continue

    q = pd.read_csv(q_file, sep=" ", header=None)
    q.columns = [f"pop{i}" for i in range(k)]
    q["sample"] = samples

    # --- Sort by dominant population then by its proportion ---
    q["dominant_pop"] = q[[f"pop{i}" for i in range(k)]].idxmax(axis=1)
    q["dominant_prop"] = q[[f"pop{i}" for i in range(k)]].max(axis=1)
    q = q.sort_values(["dominant_pop", "dominant_prop"], ascending=[True, False]).reset_index(drop=True)

    ax = fig.add_subplot(gs[idx + 1])
    bottom = [0] * len(q)

    for pop_idx, pop in enumerate([f"pop{i}" for i in range(k)]):
        ax.bar(range(len(q)), q[pop], bottom=bottom,
               color=colors[pop_idx % len(colors)], width=1.0, edgecolor="none")
        bottom = [b + v for b, v in zip(bottom, q[pop])]

    # Add vertical lines between population groups
    group_sizes = q.groupby("dominant_pop", sort=True).size().cumsum().values
    for boundary in group_sizes[:-1]:
        ax.axvline(boundary, color="black", linewidth=0.8, linestyle="--")

    title = f"K={k}" if k != best_k else f"K={k} ← best"
    ax.set_title(title, fontsize=9)
    ax.set_xlim(0, len(q))
    ax.set_ylim(0, 1)
    ax.set_ylabel("Ancestry")
    ax.set_xticks([])

plt.savefig(f"{plot_out}/{group}_admixture.png", dpi=150, bbox_inches="tight")
print(f"Saved to {plot_out}/{group}_admixture.png")