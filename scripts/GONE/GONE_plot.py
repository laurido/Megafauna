import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr
import pickle
import glob

# For running it locally
if os.getcwd().startswith("/home/lakrids"):
    path_prefix = "/home/lakrids/GenomeDK"
else:
    path_prefix = "/faststorage/project/"

genus_list      = ["Loxodonta", "Elephas", "Boselaphus", "Panthera", "Rhinoceros", "Ceratotherium", "Diceros"]
genus_list = ["Elephas"]
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


    with open(f"{path_prefix}/megaFauna/sa_megafauna/results/{group}/parameters_{group}.pkl", "rb") as f: 
        params = pickle.load(f)

    generation = float(params.get("generation"))
    colour = params.get("colour")

    ne_df = pd.read_csv(glob.glob(f"{path_prefix}/megaFauna/sa_megafauna/results/{group}/GONE/*_GONE2_Ne")[0], sep="\t")
    ne_df["time_years_ago"] = ne_df["Generation"] * generation
    ne_df[:10]

    # Plotting
    plt.figure(figsize=(8, 5))
    plt.plot(ne_df["time_years_ago"], ne_df["Ne_diploids"], linestyle='-', color=colour)
    plt.gca()  # Optional: makes the most recent generation appear on the left
    #plt.title(f"Past effective population sizes {species} (GONE2)")
    plt.xlabel("Years Ago", fontsize=14)
    plt.ylabel("Population Size", fontsize=14)
    plt.grid(True)
    plt.tight_layout()
    plt.gca() # so it matches the paper figure
    plt.savefig(f"{path_prefix}/megaFauna/sa_megafauna/results/shared/GONE/{group}_GONE.png")
    plt.show()