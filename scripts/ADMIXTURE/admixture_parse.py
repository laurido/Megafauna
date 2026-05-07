import pandas as pd
import numpy as np
import sys
import os

# -----------------------------
# Manual K overrides
# None or 1 = treat as ONE population
# -----------------------------
k_overrides = {
    "Boselaphus_tragocamelus": 1,
    "Rhinoceros_unicornis": 1,
    "Panthera_tigris": 3,
    "Panthera_uncia": 1,
    "Panthera_leo": 2
}

# -----------------------------
# Inputs
# -----------------------------
group         = sys.argv[1]
admixture_out = sys.argv[2]
k_min         = int(sys.argv[3])
k_max         = int(sys.argv[4])

# -----------------------------
# Load sample order
# -----------------------------
fam = pd.read_csv(f"{admixture_out}/{group}_pruned_remapped.fam",
                  sep=" ", header=None)
samples = fam[1].values

# -----------------------------
# Parse CV errors (ignore K=1)
# -----------------------------
cv = {}
with open(f"{admixture_out}/cv_errors.txt") as f:
    for line in f:
        parts = line.strip().split()
        k_val  = int(parts[2].replace("(K=", "").replace("):", ""))
        cv_val = float(parts[-1])
        if k_val >= 2:          # <-- IMPORTANT
            cv[k_val] = cv_val

best_k_admixture = min(cv, key=cv.get)

# -----------------------------
# Determine final K
# -----------------------------
override = k_overrides.get(group, None)

if override is not None:
    if override >= 2:
        final_k = override
        print(f"{group}: Using OVERRIDDEN K = {final_k}")
    else:
        # override == 1 → treat as one population
        final_k = None
        print(f"{group}: Override requests ONE population (no ADMIXTURE split)")
else:
    final_k = best_k_admixture
    print(f"{group}: Using ADMIXTURE best K = {final_k} (CV = {cv[final_k]:.4f})")

# -----------------------------
# Case 1: final_k is None → ONE POPULATION
# -----------------------------
if final_k is None:
    out_dir = f"{admixture_out}/pop_files/K1"
    os.makedirs(out_dir, exist_ok=True)

    df = pd.DataFrame({"sample_id": samples, "population": "pop0"})
    df.to_csv(f"{out_dir}/{group}_population_assignments_K1.tsv",
              sep="\t", index=False)
    df["sample_id"].to_csv(f"{out_dir}/pop0.txt", index=False, header=False)

    print(f"{group}: Saved one-population assignment to {out_dir}")
    sys.exit(0)

# -----------------------------
# Case 2: final_k >= 2 → load Q file
# -----------------------------
q_file = f"{admixture_out}/{group}_pruned_remapped.{final_k}.Q"

if not os.path.exists(q_file):
    # fallback to one population
    print(f"WARNING: Q file for K={final_k} not found → treating as one population")

    out_dir = f"{admixture_out}/pop_files/K1"
    os.makedirs(out_dir, exist_ok=True)

    df = pd.DataFrame({"sample_id": samples, "population": "pop0"})
    df.to_csv(f"{out_dir}/{group}_population_assignments_K1.tsv",
              sep="\t", index=False)
    df["sample_id"].to_csv(f"{out_dir}/pop0.txt", index=False, header=False)
    sys.exit(0)

# Load Q file
q = pd.read_csv(q_file, sep=" ", header=None)
q.columns = [f"pop{i}" for i in range(final_k)]
q["sample_id"] = samples

# Assign each sample to its max‑probability cluster
q["population"] = q[[f"pop{i}" for i in range(final_k)]].idxmax(axis=1)

# -----------------------------
# Output directory
# -----------------------------
out_dir = f"{admixture_out}/pop_files/K{final_k}"
os.makedirs(out_dir, exist_ok=True)

# -----------------------------
# Write population files
# -----------------------------
for pop, group_df in q.groupby("population"):
    out_path = f"{out_dir}/{pop}.txt"
    group_df["sample_id"].to_csv(out_path, index=False, header=False)
    print(f"  {pop}: {len(group_df)} samples → {out_path}")

# Save full assignment table
q.to_csv(f"{out_dir}/{group}_population_assignments_K{final_k}.tsv",
         sep="\t", index=False)

print(f"\nDone. Population files saved to {out_dir}")
