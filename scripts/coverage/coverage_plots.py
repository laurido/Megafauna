import os
import pickle
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from plotnine import ggplot, aes, geom_point, scale_x_log10, scale_y_continuous, scale_color_manual, labs, theme_bw, theme, element_text
import matplotlib.gridspec as gridspec

# For running it locally
if os.getcwd().startswith("/home/lakrids"):
    path_prefix = "/home/lakrids/GenomeDK"
else:
    path_prefix = "/faststorage/project/"

genus_list      = ["Loxodonta", "Elephas", "Boselaphus", "Panthera", "Rhinoceros", "Ceratotherium", "Diceros"]
#genus_list = ["Elephas"]
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

    with open(f"{path_prefix}/megaFauna/sa_megafauna/results/{group}/parameters_{group}.pkl", "rb") as file:
        parameters = pickle.load(file)

    df_regions = pd.read_table(f"{path_prefix}/megaFauna/sa_megafauna/data/{ref_folder}/ref/regions_{ref_folder}_updated.txt")
    df_coverage = pd.read_table(f"{path_prefix}/megaFauna/sa_megafauna/data/{ref_folder}/ref/samples_coverage_stats.txt", delimiter="\t")
    df_coverage = df_coverage[df_coverage.IND_ID.isin(inds_to_include)]
    # edit dataframe so we have cov_A >= 10 and cov_len_A >= 0.95 --> these will get a true in the "cov_filter" column
    df_coverage['cov_filter'] = (df_coverage['cov_A'] >= 10) & (df_coverage['cov_len_A'] >= 0.90)

    # set plot style
    sns.set_theme(style='whitegrid')
    plt.rcParams.update({'font.size': 16})

    # figure & grid layout
    fig = plt.figure(figsize=(10, 6))
    gs = gridspec.GridSpec(2, 2, width_ratios=[4, 1], height_ratios=[1, 4], wspace=0.065, hspace=0.1) # adjust so nothing overlaps

    ### COVERAGE ###
    # Top histogram (cov_chrA)
    ax_top = plt.subplot(gs[0, 0])
    sns.histplot(data=df_coverage, x='cov_A', bins=50, ax=ax_top, legend=False, hue='cov_filter',
                    palette={True: 'orange', False: 'darkgrey'})
    ax_top.set_xlim(0, 85)
    #ax_top.set_xscale('log')
    ax_top.set_ylabel('Count')
    ax_top.set_xticks([])
    ax_top.set_xlabel('')
    ax_top.tick_params(bottom=False)

    ### FRACTION ###
    # Right histogram (frac_cov_chrA)
    ax_right = plt.subplot(gs[1, 1])
    sns.histplot(data=df_coverage, y='cov_len_A', bins=50, ax=ax_right, legend=False, hue='cov_filter',
                    palette={True: 'orange', False: 'darkgrey'})
    ax_right.set_ylim(-0.05, 1.05)
    ax_right.set_ylabel('')
    ax_right.set_yticks([])
    ax_right.set_xlabel('Count')
    ax_right.tick_params(left=False)

    ### SCATTER PLOT ###
    ax_main = plt.subplot(gs[1, 0])
    sns.scatterplot(data=df_coverage, x='cov_A', y='cov_len_A', ax=ax_main, hue='cov_filter',
                    palette={True: 'orange', False: 'darkgrey'})
    ax_main.set_xlim(0, 85)
    ax_main.set_ylim(-0.05, 1.05)
    ax_main.set_xlabel('Depth')
    ax_main.set_ylabel('Genome length covered')
    #fig.suptitle('Correlation between Coverage Depth and Genome Coverage (Panthera tigris)', fontsize=14)
    os.makedirs(f"{path_prefix}/megaFauna/sa_megafauna/results/shared/coverage/", exist_ok=True)
    os.makedirs(f"{path_prefix}/megaFauna/sa_megafauna/results/{group}/coverage/", exist_ok=True)
    plt.savefig(f"{path_prefix}/megaFauna/sa_megafauna/results/shared/coverage/Correlation_Coverage_Depth_and_Genome_Coverage_{group}.png")
    plt.savefig(f"{path_prefix}/megaFauna/sa_megafauna/results/{group}/coverage/Correlation_Coverage_Depth_and_Genome_Coverage_{group}.png")
    plt.close()