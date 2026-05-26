import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns

if os.getcwd().startswith("/home/lakrids"):
    path_prefix = "/home/lakrids/GenomeDK"
else:
    path_prefix = "/faststorage/project/"

genus_list      = ["Loxodonta", "Elephas", "Boselaphus", "Panthera", "Rhinoceros", "Ceratotherium", "Diceros"]
data            = pd.concat([pd.read_table(f) for f in [f"{path_prefix}/megaFauna/sa_megafauna/metadata/samples_{genus}.txt" for genus in genus_list]], ignore_index=True) # add all SRR accession from genus_list in a single data frame
data            = data.reset_index(drop=True)
ref_folders     = sorted(set(data.REFERENCE_FOLDER)) # list of references needed to map the SRR accessions

references      = pd.read_table(f"{path_prefix}/megaFauna/sa_megafauna/metadata/references.txt")
references      = references.loc[references.REFERENCE_FOLDER.isin(ref_folders)] # filter references to keep only the ones needed to map the SRR accessions
references      = references.reset_index(drop=True)

## make data frame that contains names of species-specific folders and the reference folders used to map the species
species_and_refs = pd.DataFrame({"FOLDER": data.FOLDER, "REFERENCE_FOLDER": data.REFERENCE_FOLDER, "GVCF_FOLDER": [data.GENUS.iloc[jj] + "_" + data.SPECIES.iloc[jj] for jj in range(data.shape[0])]}).drop_duplicates()
species_and_refs = species_and_refs.reset_index(drop=True)

## merge dataframes to have species-specific folder, reference folder and fastq ftps in same dataframe
species_and_refs = species_and_refs.merge(references, how = "left")

##########################################
#    	     --- Preparations ---
##########################################
colors = plt.cm.tab10.colors

for i in range(species_and_refs.shape[0]):
    group      = species_and_refs.FOLDER[i]
    ref_folder = species_and_refs.REFERENCE_FOLDER[i]
    inds_to_include = (data[data["FOLDER"] == group]["IND_ID"].drop_duplicates().tolist())

    merged_counts_df = pd.read_csv(f"{path_prefix}/megaFauna/sa_megafauna/data/{group}/VCF/snp_counts/merged_counts.txt", sep="\t")
    coverage_df = pd.read_csv(f"{path_prefix}/megaFauna/sa_megafauna/data/{ref_folder}/ref/samples_coverage_stats_{group}_filtered.txt", sep="\t")
    coverage_df = coverage_df[coverage_df.IND_ID.isin(inds_to_include)]
    merged_df = merged_counts_df.merge(coverage_df[['IND_ID', 'len_covered_raw_A']],
                                  left_on='IND_ID', right_on='IND_ID', how='left')
    merged_df['het_autosomal'] = merged_df['autosomal'] / merged_df['len_covered_raw_A']

    # --- Load population assignments (same logic as PCA script) ---
    assignment_file = f"{path_prefix}/megaFauna/sa_megafauna/results/{group}/{group}_population_assignments.tsv"
    pop_col = None

    if os.path.exists(assignment_file):
        assignments = pd.read_csv(assignment_file, sep="\t")[["sample_id", "population"]]
        merged_df = merged_df.merge(assignments, left_on="IND_ID", right_on="sample_id", how="left")
        pop_col = "population"

    # --- Sort by population then by heterozygosity within each population ---
    if pop_col and pop_col in merged_df.columns:
        merged_df = merged_df.sort_values([pop_col]).reset_index(drop=True)
    
    fig, ax = plt.subplots(figsize=(10, 5))

    if pop_col and pop_col in merged_df.columns:
        populations = sorted(merged_df[pop_col].dropna().unique())
        for pop_idx, pop in enumerate(populations):
            subset = merged_df[merged_df[pop_col] == pop]
            color = colors[pop_idx % len(colors)]
            ax.scatter(subset.index, subset["het_autosomal"],
                       label=f"{pop} (n={len(subset)})",
                       color=color, edgecolor='black', s=40)
            
            # Per-cluster mean
            pop_mean = subset["het_autosomal"].mean()
            ax.hlines(pop_mean, subset.index.min(), subset.index.max(),
                      color=color, linestyle='--', linewidth=1.5)
            ax.text(subset.index.max(), pop_mean, f"  {pop_mean:.5f}",
                    va='bottom', ha='left', fontsize=8, color=color)
        ax.legend(fontsize=8)
    else:
        ax.scatter(merged_df.index, merged_df["het_autosomal"], edgecolor='black', s=40)
        # Fallback: single global mean
        mean_val = merged_df['het_autosomal'].mean()
        ax.axhline(mean_val, color='black', linestyle='--')
        ax.text(len(merged_df) * 0.01, mean_val, f"mean = {mean_val:.5f}",
                va='bottom', ha='left', fontsize=9, color='black')

    ax.set_xticks([])
    ax.set_xlabel("Samples", fontsize=12)
    ax.set_ylabel("Autosome Heterozygosity", fontsize=12)
    ax.set_title(group)
    plt.tight_layout()

    os.makedirs(f"{path_prefix}/megaFauna/sa_megafauna/results/shared/heterozygosity/", exist_ok=True)
    os.makedirs(f"{path_prefix}/megaFauna/sa_megafauna/results/{group}/heterozygosity/", exist_ok=True)
    plt.savefig(f"{path_prefix}/megaFauna/sa_megafauna/results/shared/heterozygosity/heterozygosity_{group}.png")
    plt.savefig(f"{path_prefix}/megaFauna/sa_megafauna/results/{group}/heterozygosity/heterozygosity_{group}.png")
    plt.close()