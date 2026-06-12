import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns

if os.getcwd().startswith("/home/lakrids"):
    path_prefix = "/home/lakrids/GenomeDK"
else:
    path_prefix = "/faststorage/project/"

genus_list      = ["Boselaphus", "Ceratotherium", "Diceros", "Elephas", "Loxodonta", "Panthera", "Rhinoceros"]
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
    parts = group.split('_')
    formatted = f"$\it{{{parts[0][0]}. {' '.join(parts[1:])}}}$"
    ax.set_xticks([])
    ax.set_xlabel("Samples", fontsize=12)
    ax.set_ylabel("Autosome Heterozygosity", fontsize=12)
    ax.set_title(formatted)
    plt.tight_layout()

    os.makedirs(f"{path_prefix}/megaFauna/sa_megafauna/results/shared/heterozygosity/", exist_ok=True)
    os.makedirs(f"{path_prefix}/megaFauna/sa_megafauna/results/{group}/heterozygosity/", exist_ok=True)
    plt.savefig(f"{path_prefix}/megaFauna/sa_megafauna/results/shared/heterozygosity/heterozygosity_{group}.png")
    plt.savefig(f"{path_prefix}/megaFauna/sa_megafauna/results/{group}/heterozygosity/heterozygosity_{group}.png")
    plt.close()

# --- Combined boxplot across all species and populations ---
all_data = []

for i in range(species_and_refs.shape[0]):
    group      = species_and_refs.FOLDER[i]
    ref_folder = species_and_refs.REFERENCE_FOLDER[i]
    inds_to_include = (data[data["FOLDER"] == group]["IND_ID"].drop_duplicates().tolist())

    try:
        merged_counts_df = pd.read_csv(f"{path_prefix}/megaFauna/sa_megafauna/data/{group}/VCF/snp_counts/merged_counts.txt", sep="\t")
        coverage_df = pd.read_csv(f"{path_prefix}/megaFauna/sa_megafauna/data/{ref_folder}/ref/samples_coverage_stats_{group}_filtered.txt", sep="\t")
        coverage_df = coverage_df[coverage_df.IND_ID.isin(inds_to_include)]
        merged_df = merged_counts_df.merge(coverage_df[['IND_ID', 'len_covered_raw_A']],
                                      on='IND_ID', how='left')
        merged_df['het_autosomal'] = merged_df['autosomal'] / merged_df['len_covered_raw_A']

        assignment_file = f"{path_prefix}/megaFauna/sa_megafauna/results/{group}/{group}_population_assignments.tsv"
        if os.path.exists(assignment_file):
            assignments = pd.read_csv(assignment_file, sep="\t")[["sample_id", "population"]]
            merged_df = merged_df.merge(assignments, left_on="IND_ID", right_on="sample_id", how="left")
            merged_df['group'] = group
            all_data.append(merged_df[['het_autosomal', 'population', 'group']])
    except FileNotFoundError:
        continue

combined_df = pd.concat(all_data, ignore_index=True)

# Build labels with sample counts
combined_df['label'] = combined_df.apply(
    lambda r: f"{r['population']}\n(n={combined_df[(combined_df['population']==r['population']) & (combined_df['group']==r['group'])].shape[0]})",
    axis=1
)

# Get unique ordered groups and populations
groups = combined_df['group'].unique()
# Calculate width per subplot based on number of populations
n_pops_per_group = [len(sorted(combined_df[combined_df['group'] == group]['population'].dropna().unique())) 
                    for group in groups]
box_width = 0.4  # width per population in inches
subplot_widths = [max(n * box_width, box_width * 1.5) for n in n_pops_per_group]  # minimum 2 boxes wide

fig, axes = plt.subplots(1, len(groups), figsize=(sum(subplot_widths), 6),
                          sharey=True, gridspec_kw={'wspace': 0.0, 'width_ratios': subplot_widths})

if len(groups) == 1:
    axes = [axes]

for idx, (ax, group) in enumerate(zip(axes, groups)):
    subset = combined_df[combined_df['group'] == group]
    populations = sorted(subset['population'].dropna().unique())
    
    data_to_plot = [subset[subset['population'] == pop]['het_autosomal'].values 
                    for pop in populations]
    n_labels = [f"{pop}\n(n={len(subset[subset['population']==pop])})" 
                for pop in populations]
    
    # Tighter positions - spacing of 0.6 instead of default 1.0
    positions = [j * 0.4 for j in range(len(populations))]
    
    bp = ax.boxplot(data_to_plot, patch_artist=True, labels=n_labels,
                    widths=0.3, positions=positions)
    
    # Set x limits to avoid too much padding
    ax.set_xlim(positions[0] - 0.2, positions[-1] + 0.2)
    
    for patch, pop_idx in zip(bp['boxes'], range(len(populations))):
        patch.set_facecolor(colors[pop_idx % len(colors)])
    for median in bp['medians']:
        median.set_color('black')
    # Format species name as L. africana in italics
    parts = group.split('_')
    formatted = f"$\it{{{parts[0][0]}. {' '.join(parts[1:])}}}$"
    ax.set_title(formatted, fontsize=10, rotation=20)    
    ax.tick_params(axis='x', rotation=45, labelsize=7)
    ax.yaxis.grid(True, linestyle='--', alpha=0.7)
    ax.set_axisbelow(True)

    if idx < len(groups) - 1:
        ax.spines['right'].set_visible(True)
        ax.spines['right'].set_linewidth(1.5)
    
    if idx > 0:
        ax.spines['left'].set_visible(False)
        ax.tick_params(left=False)

axes[0].set_ylabel("Autosome Heterozygosity", fontsize=12)
#fig.suptitle("Heterozygosity by Population", fontsize=14, y=1.02)
plt.tight_layout()

os.makedirs(f"{path_prefix}/megaFauna/sa_megafauna/results/shared/heterozygosity/", exist_ok=True)
plt.savefig(f"{path_prefix}/megaFauna/sa_megafauna/results/shared/heterozygosity/heterozygosity_boxplot_combined.png",
            bbox_inches='tight')
plt.close()