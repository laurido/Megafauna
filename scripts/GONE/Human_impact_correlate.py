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

    generation = round(float(params.get("generation")))
    biogeo = params.get("biogeography")
    if biogeo == "Palearctic":
        continue
    if biogeo == "Afrotropical":
        country = "Kenya"
    else:
        country = "India"

    ne_df = pd.read_csv(glob.glob(f"{path_prefix}/megaFauna/sa_megafauna/results/{group}/GONE/*_GONE2_Ne")[0], sep="\t")
    ne_df["time_years_ago"] = ne_df["Generation"] * generation

    human_df = pd.read_csv(f'{path_prefix}/megaFauna/sa_megafauna/data/human_impact/{biogeo}.csv')
    human_df.columns = ["Year", "Inhabitants", "Cropland", "Urban"] + list(human_df.columns[4:])

    os.makedirs(f"{path_prefix}/megaFauna/sa_megafauna/results/shared/human_impact/", exist_ok=True)

    # some data wrangling to get the years up to the same axis :)
    reference_year = 2023
    ne_df["Year"] = reference_year - ne_df["time_years_ago"]
    human_vars = ["Inhabitants", "Cropland", "Urban"]

    # idea: interpolate the human impact variables (they are smooth and slow changing, so it should work). You still need to be careful though.
    # Step 1: Create a full year range covering Ne_diploids years 
    full_years = pd.DataFrame({"Year": ne_df["Year"]})
    full_years

    # Step 2: Reindex human impact data to this full range
    impact_full = (
        human_df.set_index("Year")
        .reindex(full_years["Year"])  # this adds missing years as NaN
        .interpolate(method="index", limit_direction="both")  # interpolate + extrapolate
        .reset_index()
    )

    # Step 3: Merge with Ne_diploids
    aligned = ne_df[["Year", "Ne_diploids"]].merge(
        impact_full[["Year"] + human_vars],
        on="Year",
        how="left"
    )
    """
    Use scipy.stats.spearmanr to calculate the correlation between each variable in human_vars and Ne_diploids:
    """
    # Step 4: Drop any remaining NaNs
    aligned = aligned.dropna()

    # Step 5: Compute Spearman correlations
    correlations = {}

    for var in human_vars:
        rho, pval = spearmanr(aligned[var], aligned["Ne_diploids"])
        correlations[var] = (rho, pval)
