import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import os
import math

genus_list      = ["Loxodonta", "Boselaphus", "Panthera", "Elephas", "Naja", "Rhinoceros", "Ceratotherium", "Diceros"]
data            = pd.concat([pd.read_table(f) for f in [f"/faststorage/project/megaFauna/sa_megafauna/metadata/samples_{genus}.txt" for genus in genus_list]], ignore_index=True)
data            = data.reset_index(drop=True)
ref_folders     = sorted(set(data.REFERENCE_FOLDER))

references      = pd.read_table(f"/faststorage/project/megaFauna/sa_megafauna/metadata/references.txt")
references      = references.loc[references.REFERENCE_FOLDER.isin(ref_folders)]
references      = references.reset_index(drop=True)

species_and_refs = pd.DataFrame({
    "FOLDER": data.FOLDER,
    "REFERENCE_FOLDER": data.REFERENCE_FOLDER,
    "GVCF_FOLDER": [data.GENUS.iloc[jj] + "_" + data.SPECIES.iloc[jj] for jj in range(data.shape[0])]
}).drop_duplicates().reset_index(drop=True)

species_and_refs = species_and_refs.merge(references, how="left")

def remove_iqr_outliers(df, columns, sample_col="sample_id", multiplier=5):
    mask = pd.Series(True, index=df.index)
    
    for col in columns:
        Q1 = df[col].quantile(0.25)
        Q3 = df[col].quantile(0.75)
        IQR = Q3 - Q1
        
        lower = Q1 - multiplier * IQR
        upper = Q3 + multiplier * IQR
        
        mask &= df[col].between(lower, upper)
        
        print(f"  {col}: Q1={Q1:.3f}, Q3={Q3:.3f}, IQR={IQR:.3f} → keeping [{lower:.3f}, {upper:.3f}]")
    
    df_clean   = df[mask].copy()   # <-- added .copy()
    df_removed = df[~mask].copy()  # <-- added .copy()

    print(f"  Samples before filtering: {len(df)}")
    print(f"  Samples removed:          {len(df_removed)}")
    print(f"  Samples after filtering:  {len(df_clean)}")

    if len(df_removed) > 0 and sample_col in df.columns:
        print(f"  Removed samples:")
        print(df_removed[[sample_col, "PC1", "PC2"]].to_string(index=False))

    return df_clean, df_removed

##########################################
#        --- PCA GRID PLOT ---
##########################################

pca_results   = []   # store (group, df, explained_var)
all_kept      = []   # track kept samples across all groups
all_removed   = []   # track removed samples across all groups

for i in range(species_and_refs.shape[0]):
    group      = species_and_refs.FOLDER[i]
    ref_folder = species_and_refs.REFERENCE_FOLDER[i]

    df = pd.read_csv(f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/PCA/pca_dataset_{group}.txt", sep='\t')

    with open(f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/PCA/pca_model_{group}.pkl", "rb") as f:
        loaded_model = pickle.load(f)

    explained_var = loaded_model.explained_variance_ratio_

    print(f"\n{group}")
    df_clean, df_removed = remove_iqr_outliers(df=df, columns=["PC1", "PC2"])

    # Tag with group so you know which species each sample belongs to
    df_clean["group"]   = group
    df_removed["group"] = group

    all_kept.append(df_clean)
    all_removed.append(df_removed)

    pca_results.append((group, df_clean, explained_var))

# -------------------------------
# Save kept/removed sample lists
# -------------------------------

kept_samples    = pd.concat(all_kept,    ignore_index=True)
removed_samples = pd.concat(all_removed, ignore_index=True)

out_dir = f"/faststorage/project/megaFauna/sa_megafauna/results"

kept_samples.to_csv(f"{out_dir}/kept_samples.txt",    sep="\t", index=False)
removed_samples.to_csv(f"{out_dir}/removed_samples.txt", sep="\t", index=False)

print(f"\n--- Summary ---")
print(f"Total kept:    {len(kept_samples)}")
print(f"Total removed: {len(removed_samples)}")
print(f"Saved to {out_dir}/kept_samples.txt and removed_samples.txt")

# -------------------------------
# Build grid plot
# -------------------------------

n = len(pca_results)
cols = 4
rows = math.ceil(n / cols)

fig, axes = plt.subplots(rows, cols, figsize=(18, 10))
axes = axes.flatten()

for idx, (group, df, explained_var) in enumerate(pca_results):
    ax = axes[idx]
    ax.scatter(df["PC1"], df["PC2"], edgecolor='black', s=40)
    ax.set_xlabel(f"PC1 ({explained_var[0]*100:.2f}%)")
    ax.set_ylabel(f"PC2 ({explained_var[1]*100:.2f}%)")
    ax.set_title(group)
    ax.grid(True)

# Hide any unused subplots
for j in range(idx + 1, len(axes)):
    axes[j].axis("off")

plt.tight_layout()
plt.savefig(f"{out_dir}/PCA_all_species_grid.png", dpi=300)
plt.show()