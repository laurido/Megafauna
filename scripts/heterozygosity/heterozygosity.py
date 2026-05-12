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

for i in range(species_and_refs.shape[0]):
    # Initialising folders and variables for putting in the functions

    group      = species_and_refs.FOLDER[i]
    ref_folder = species_and_refs.REFERENCE_FOLDER[i]
    inds_to_include = (data[data["FOLDER"] == group]["IND_ID"].drop_duplicates().tolist())

    merged_counts_df = pd.read_csv(f"{path_prefix}/megaFauna/sa_megafauna/data/{group}/VCF/snp_counts/merged_counts.txt", sep="\t")
    coverage_df = pd.read_csv(f"{path_prefix}/megaFauna/sa_megafauna/data/{ref_folder}/ref/samples_coverage_stats_{group}_filtered.txt", sep="\t")
    coverage_df = coverage_df[coverage_df.IND_ID.isin(inds_to_include)]
    merged_df = merged_counts_df.merge(coverage_df[['IND_ID', 'len_covered_raw_A']],
                                  left_on='IND_ID', right_on='IND_ID', how='left')
    merged_df['het_autosomal'] = merged_df['autosomal'] / merged_df['len_covered_raw_A']

    sns.scatterplot(data=merged_df, x='IND_ID', y='het_autosomal', edgecolor='black')
    plt.xticks([])
    #plt.title(f"Individual autosomal heterozygosity ({species})")

    # get mean value
    mean_val = merged_df['het_autosomal'].mean()
    plt.axhline(mean_val, color='black', linestyle='--', label="mean")
    plt.xlabel("Samples", fontsize=12)
    plt.ylabel("Autosome Heterozygosity", fontsize=12)
    plt.legend()
    plt.tight_layout()
    os.makedirs(f"{path_prefix}/megaFauna/sa_megafauna/results/shared/heterozygosity/", exist_ok=True)
    os.makedirs(f"{path_prefix}/megaFauna/sa_megafauna/results/{group}/heterozygosity/", exist_ok=True)
    plt.savefig(f"{path_prefix}/megaFauna/sa_megafauna/results/shared/heterozygosity/heterozygosity_{group}.png")
    plt.savefig(f"{path_prefix}/megaFauna/sa_megafauna/results/{group}/heterozygosity/heterozygosity_{group}.png")
    plt.close()