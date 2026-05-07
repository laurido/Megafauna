import pandas as pd
import os
import sys

# ─── Configuration ───────────────────────────────────────────────
admixture_threshold = 0.85   # minimum ancestry to be considered non-admixed
fst_threshold       = 0.05   # minimum Fst to consider a split meaningful

# Manual overrides:
#  - value >= 2 → force that K
#  - value == 1 → treat as ONE population (no ADMIXTURE split)
#  - missing → use CV-best K
k_overrides = {
    "Boselaphus_tragocamelus": 1,
    "Rhinoceros_unicornis": 1,
    "Panthera_tigris": 3,
    "Panthera_uncia": 1,
    "Panthera_leo": 2
}
# ─────────────────────────────────────────────────────────────────
genus_list = ["Ceratotherium"]
genus_list = ["Loxodonta", "Elephas", "Boselaphus", "Panthera", "Rhinoceros", "Ceratotherium", "Diceros"]

# Load metadata
data = pd.concat([
    pd.read_table(f"/faststorage/project/megaFauna/sa_megafauna/metadata/samples_{genus}.txt")
    for genus in genus_list
], ignore_index=True)

data = data.reset_index(drop=True)
ref_folders = sorted(set(data.REFERENCE_FOLDER))

references = pd.read_table(f"/faststorage/project/megaFauna/sa_megafauna/metadata/references.txt")
references = references.loc[references.REFERENCE_FOLDER.isin(ref_folders)].reset_index(drop=True)

species_and_refs = pd.DataFrame({
    "FOLDER": data.FOLDER,
    "REFERENCE_FOLDER": data.REFERENCE_FOLDER,
    "GVCF_FOLDER": [data.GENUS.iloc[j] + "_" + data.SPECIES.iloc[j] for j in range(data.shape[0])]
}).drop_duplicates().reset_index(drop=True)

species_and_refs = species_and_refs.merge(references, how="left")

# ─────────────────────────────────────────────────────────────────

for i in range(species_and_refs.shape[0]):

    group      = species_and_refs.FOLDER[i]
    ref_folder = species_and_refs.REFERENCE_FOLDER[i]

    print(f"\n{'='*50}")
    print(f"{group}")
    print(f"{'='*50}")

    sample_df = pd.read_csv(
        f"/faststorage/project/megaFauna/sa_megafauna/data/{ref_folder}/ref/samples_coverage_stats_{group}_filtered.txt",
        sep="\t"
    )
    out_dir = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}"
    admixture_base = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/ADMIXTURE"
    cv_file        = f"{admixture_base}/cv_errors.txt"

    if not os.path.exists(cv_file):
        print("  No cv_errors.txt found, skipping")
        continue

    # --- Load CV errors (ignore K=1) ---
    cv = {}
    with open(cv_file) as f:
        for line in f:
            parts  = line.strip().split()
            k_val  = int(parts[2].replace("(K=", "").replace("):", ""))
            cv_val = float(parts[-1])
            if k_val >= 2:
                cv[k_val] = cv_val

    cv_best_k = min(cv, key=cv.get)

    # --- Apply override logic ---
    override = k_overrides.get(group, None)

    if override is not None:
        if override >= 2:
            best_k = override
            print(f"  Using OVERRIDDEN K={best_k} (CV best was {cv_best_k})")
        else:
            # override == 1 → one population
            best_k = None
            print(f"  Override requests ONE population (no ADMIXTURE split)")
    else:
        best_k = cv_best_k
        print(f"  Using CV best K={best_k}")

    # Print CV table
    print("\n  CV errors:")
    for k, v in sorted(cv.items()):
        marker = " ←" if k == cv_best_k else ""
        print(f"    K={k}: {v:.5f}{marker}")

    # ───────────────────────────────────────────────────────────────
    # CASE 1: ONE POPULATION (override=1 or missing Q file)
    # ───────────────────────────────────────────────────────────────
    if best_k is None:
        print("\n  Treating all samples as one population")

        q_clean = pd.DataFrame({
            "sample_id": sample_df["IND_ID"].drop_duplicates().tolist(),
            "population": "pop0"
        })

        os.makedirs(out_dir, exist_ok=True)

        q_clean["sample_id"].to_csv(f"{out_dir}/pop0.txt", index=False, header=False)
        q_clean.to_csv(f"{out_dir}/{group}_population_assignments.tsv",
                       sep="\t", index=False)
        # Write population list file
        with open(f"{out_dir}/population_list.txt", "w") as f:
            f.write("pop0\n")


        print(f"  Saved one-population files to {out_dir}")
        continue

    # ───────────────────────────────────────────────────────────────
    # CASE 2: MULTIPLE POPULATIONS (best_k >= 2)
    # ───────────────────────────────────────────────────────────────
    assignment_file = f"{admixture_base}/pop_files/K{best_k}/{group}_population_assignments_K{best_k}.tsv"

    if not os.path.exists(assignment_file):
        print(f"\n  Assignment file missing → treating as one population")

        q_clean = pd.DataFrame({
            "sample_id": sample_df["IND_ID"].drop_duplicates().tolist(),
            "population": "pop0"
        })

        os.makedirs(out_dir, exist_ok=True)

        q_clean["sample_id"].to_csv(f"{out_dir}/pop0.txt", index=False, header=False)
        q_clean.to_csv(f"{out_dir}/{group}_population_assignments.tsv",
                       sep="\t", index=False)
        continue

    # Load assignment file
    q = pd.read_csv(assignment_file, sep="\t")
    pop_cols = [f"pop{i}" for i in range(best_k)]

    # --- Parse Fst ---
    print("\n  Fst between populations:")
    fst_ok = True
    log_file = f"{admixture_base}/admixture_K{best_k}.log"
    
    if os.path.exists(log_file):
        reading_fst = False
        header_pops = None
    
        with open(log_file) as f:
            for line in f:
                stripped = line.strip()
    
                # Find the "Fst divergences" line and then expect the header next
                if "Fst divergences" in stripped:
                    print(f"    {stripped}")
                    reading_fst = True
                    continue
                
                if not reading_fst or not stripped:
                    continue
                
                parts = stripped.split()
    
                # HEADER ROW: something like "Pop0 Pop1 ..."
                if header_pops is None and parts[0].startswith("Pop") and len(parts) >= 1:
                    header_pops = parts  # e.g. ["Pop0", "Pop1"]
                    print(f"    Header populations: {', '.join(header_pops)}")
                    continue
                
                # DATA ROWS: "PopX <values...>"
                if parts[0].startswith("Pop"):
                    row_pop = parts[0]
                    values = parts[1:]
    
                    for j, val in enumerate(values):
                        try:
                            fst_val = float(val)
                            col_pop = header_pops[j]  # safe now, header_pops has all names
                            print(f"    {row_pop} vs {col_pop}: Fst = {fst_val:.4f}")
    
                            if fst_val < fst_threshold:
                                print(f"      WARNING: Fst={fst_val:.4f} < {fst_threshold}")
                                fst_ok = False
    
                        except (ValueError, IndexError):
                            # Value not numeric or header shorter than expected
                            pass
    else:
        print("  Log file not found")


    # --- Flag admixed individuals ---
    q["dominant_prop"] = q[pop_cols].max(axis=1)
    q["is_admixed"]    = q["dominant_prop"] < admixture_threshold

    print(f"\n  Admixture summary:")
    print(f"    Total samples:   {len(q)}")
    print(f"    Clean samples:   {len(q[~q['is_admixed']])}")
    print(f"    Admixed samples: {len(q[q['is_admixed']])}")

    # --- Filter out admixed ---
    q_clean = q[~q["is_admixed"]].copy()

    if len(q_clean) == 0:
        print("  ERROR: No samples left after filtering")
        continue

    # --- Find biggest population ---
    pop_sizes = q_clean.groupby("population").size().sort_values(ascending=False)
    biggest_pop = pop_sizes.index[0]

    # Write population list sorted by size (largest first)
    pop_list_path = f"{out_dir}/population_list.txt"
    with open(pop_list_path, "w") as f:
        for pop in pop_sizes.index:   # already sorted descending
            f.write(f"{pop}\n")

    print(f"  Population list saved to {pop_list_path}")


    print("\n  Population sizes:")
    for pop, size in pop_sizes.items():
        mark = " ← biggest" if pop == biggest_pop else ""
        print(f"    {pop}: {size}{mark}")

    # --- Save population files ---
    os.makedirs(out_dir, exist_ok=True)

    print("\n  Saving population files:")

    for pop, pop_df in q_clean.groupby("population"):
        out_path = f"{out_dir}/{pop}.txt"
        pop_df["sample_id"].to_csv(out_path, index=False, header=False)
        print(f"    {pop}: {len(pop_df)} → {out_path}")

    # Save clean assignment table
    q_clean.to_csv(f"{out_dir}/{group}_population_assignments.tsv",
                   sep="\t", index=False)
    print(f"  Saved clean assignments")
