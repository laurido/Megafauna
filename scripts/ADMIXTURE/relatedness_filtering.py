import pandas as pd
import os
import sys

# ─── Configuration ───────────────────────────────────────────────
relatedness_threshold = 0.125  # KING kinship threshold (1st degree relatives only)

# ─────────────────────────────────────────────────────────────────
genus_list = ["Ceratotherium"]
#genus_list = ["Boselaphus"]

genus_list = ["Loxodonta", "Elephas", "Boselaphus", "Panthera", "Rhinoceros", "Ceratotherium", "Diceros"]

data = pd.concat([pd.read_table(f"/faststorage/project/megaFauna/sa_megafauna/metadata/samples_{genus}.txt") 
                  for genus in genus_list], ignore_index=True)

groups = data[["FOLDER"]].drop_duplicates()["FOLDER"].tolist()

# ─────────────────────────────────────────────────────────────────

def filter_related_samples(relatedness_file, sample_list, threshold, group_name):
    """
    Smart relatedness filtering:
    Removes the minimum number of samples needed so that no pair has kinship > threshold.
    Uses a greedy vertex-cover approach.
    """
    if not os.path.exists(relatedness_file):
        print(f"\n  WARNING: Relatedness file not found at {relatedness_file}")
        return set()

    rel_df = pd.read_csv(relatedness_file, sep=r'\s+')
    rel_df = rel_df[rel_df['Kinship'] > threshold]

    # Keep only pairs where both samples are in this population
    rel_df = rel_df[
        rel_df['ID1'].isin(sample_list) &
        rel_df['ID2'].isin(sample_list)
    ]

    if rel_df.empty:
        print(f"  ✓ No related pairs above threshold {threshold}")
        return set()

    print(f"\n  Related pairs found: {len(rel_df)}")

    # Build conflict graph
    conflicts = []
    for _, row in rel_df.iterrows():
        conflicts.append((row['ID1'], row['ID2']))

    to_remove = set()

    # Greedy vertex cover
    while conflicts:
        # Count degrees
        counts = {}
        for a, b in conflicts:
            counts[a] = counts.get(a, 0) + 1
            counts[b] = counts.get(b, 0) + 1

        # Remove the sample with the most conflicts
        worst = max(counts, key=counts.get)
        to_remove.add(worst)

        # Remove all edges involving this sample
        conflicts = [(a, b) for (a, b) in conflicts if a != worst and b != worst]

    print(f"  Removing {len(to_remove)} samples to break all related pairs")
    for s in to_remove:
        print(f"    Removing {s}")

    return to_remove


for group in groups:
    print(f"\n{'='*50}")
    print(f"{group}")
    print(f"{'='*50}")

    with open(f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/population_list.txt") as f:
        pops = [line.strip() for line in f]

    for pop in pops:

        pop_file = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/{pop}.txt"
        
        # Load samples from pop.txt
        with open(pop_file) as f:
            sample_list = [line.strip() for line in f if line.strip()]

        print(f"  Loaded {len(sample_list)} samples from {pop}.txt")

        # Create dataframe
        q = pd.DataFrame({"sample_id": sample_list})

        # Filter relatedness
        relatedness_file = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/relatedness/{pop}/{group}_king.kin0"
        related_samples = filter_related_samples(relatedness_file, q['sample_id'].tolist(), relatedness_threshold, group)

        # Apply filter
        q_clean = q[~q["sample_id"].isin(related_samples)].copy()

        print(f"\n  Final sample counts:")
        print(f"    Total samples:             {len(q)}")
        print(f"    Removed for relatedness:   {len(related_samples)}")
        print(f"    After relatedness filter:  {len(q_clean)}")

        # Save filtered samples back to pop.txt
        out_dir = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}"

        pop_out_path = f"{out_dir}/{pop}_filtered.txt"
        q_clean["sample_id"].to_csv(pop_out_path, index=False, header=False)
        print(f"\n  Saved {len(q_clean)} samples to {pop_out_path}")

        # Save clean assignments
        q_clean.to_csv(f"{out_dir}/{group}_relatedness_filtered.tsv", sep="\t", index=False)
        print(f"  Saved filtered sample list to {out_dir}/{group}_relatedness_filtered.tsv")

        # ───────────────────────────────────────────────────────────────
    # Recompute population sizes AFTER relatedness filtering
    # ───────────────────────────────────────────────────────────────
    print("\n  Re-evaluating population order based on filtered sample counts")

    pop_sizes = {}

    for pop in pops:
        filtered_file = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/{pop}_filtered.txt"

        if os.path.exists(filtered_file):
            with open(filtered_file) as f:
                samples = [line.strip() for line in f if line.strip()]
            pop_sizes[pop] = len(samples)
        else:
            pop_sizes[pop] = 0

    # Sort populations by size (largest → smallest)
    sorted_pops = sorted(pop_sizes, key=lambda p: pop_sizes[p], reverse=True)

    # Write updated population_list.txt
    pop_list_path = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/population_list.txt"
    with open(pop_list_path, "w") as f:
        for pop in sorted_pops:
            f.write(f"{pop}\n")

    print(f"  Updated population_list.txt written with order:")
    for pop in sorted_pops:
        print(f"    {pop}: {pop_sizes[pop]} samples")
