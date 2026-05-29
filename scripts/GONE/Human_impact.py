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
def significance_stars(p):
    if p < 0.01:
        return "**"   # very significant
    elif p < 0.05:
        return "*"    # significant
    else:
        return ""

ne_dict = {}
correlations_by_species = {}

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
    
    # Add to dictionary
    entry = {"species": group, "df": ne_df}
    if country not in ne_dict:
        ne_dict[country] = []
    ne_dict[country].append(entry)

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
    # Step 5: Store in dictionary
    correlations_by_species[group] = correlations

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 12), sharey=False)

    # --- Left plot: Population vs Ne_diploids ---
    ax1.plot(
        human_df["Year"],
        human_df["Inhabitants"],
        linestyle="--",
        color="#ff6bbc",
        linewidth=4,
        label="Population [inh]"
    )
    rho, p = correlations["Inhabitants"]
    ax1.text(
        0.02, 0.95,
        f"Spearman rho = {rho:.2f}\np = {p:.2e}",
        transform=ax1.transAxes,
        fontsize=16,
        verticalalignment="top",
        bbox=dict(facecolor="white", alpha=0.7, edgecolor="none")
    )
    ax1.set_xlabel("Year", fontsize = 20)
    ax1.set_ylabel(f"{country} - Population [inh]", fontsize = 20)

    # Secondary axis for Ne_diploids
    ax1b = ax1.twinx()
    ax1b.plot(
        ne_df["Year"],
        ne_df["Ne_diploids"],
        linewidth=4,
        label="Ne_diploids"
    )
    ax1b.set_ylabel("Ne_diploids", fontsize = 20)

    # Combine legends
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax1b.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper right", fontsize = 18)
    # --- Right plot: Urban Area & Grazing Land vs Ne_diploids ---
    colors = ["#017df0", "#7eab57"]  # orange, green
    for col, color in zip([f"Urban", f"Cropland"], colors):
        ax2.plot(
        human_df["Year"],
        human_df[col],
        linestyle="--",
        color=color,
        linewidth=4,
        label=col
    )
    # Urban
    rho_u, p_u = correlations["Urban"]
    ax2.text(
        0.02, 0.95,
        f"Urban: rho = {rho_u:.2f}, p = {p_u:.2e}",
        transform=ax2.transAxes,
        fontsize=16,
        verticalalignment="top",
        bbox=dict(facecolor="white", alpha=0.7, edgecolor="none")
    )

    # Cropland
    rho_c, p_c = correlations["Cropland"]
    ax2.text(
        0.02, 0.80,
        f"Cropland: rho = {rho_c:.2f}, p = {p_c:.2e}",
        transform=ax2.transAxes,
        fontsize=16,
        verticalalignment="top",
        bbox=dict(facecolor="white", alpha=0.7, edgecolor="none")
    )

    ax2.set_xlabel("Year", fontsize = 20)
    ax2.set_ylabel(f"{country} Land-use impact in km²", fontsize = 20)

    # Secondary axis for Ne_diploids
    ax2b = ax2.twinx()
    ax2b.plot(
        ne_df["Year"],
        ne_df["Ne_diploids"],
        linewidth=4,
        label="Ne_diploids"
    )
    ax2b.set_ylabel("Ne_diploids", fontsize = 20)

    # Combine legends
    lines1, labels1 = ax2.get_legend_handles_labels()
    lines2, labels2 = ax2b.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, loc="upper right", fontsize = 18)

    # Formatting
    for ax in [ax1, ax2]:
        ax1.set_xlim(reference_year, 800) # recent years on the left, older on the right
        ax2.set_xlim(reference_year, 800)
        xticks = np.arange(800, reference_year + 1, 100)
        ax.set_xticks(xticks)
        ax.tick_params(axis="x", rotation=45, labelsize = 18)
        ax.tick_params(axis="y", labelsize=18)

    # Apply to the twin y-axes as well
    ax1b.tick_params(axis="y", labelsize=18)
    ax2b.tick_params(axis="y", labelsize=18)

    plt.tight_layout()
    plt.gca() # so it matches the paper figure
    plt.savefig(f"{path_prefix}/megaFauna/sa_megafauna/results/shared/human_impact/{group}_human_impact")
    plt.close()

for country, entries in ne_dict.items():

    fig, ax1 = plt.subplots(figsize=(20, 10))

    # Plot inhabitants
    ax1.plot(
        human_df["Year"],
        human_df["Inhabitants"],
        color="#ff6bbc",
        linestyle="--",
        linewidth=4,
        label="Inhabitants"
    )

    ax1.set_xlabel("Year", fontsize=20)
    ax1.set_ylabel(f"{country} - Population [inh]", fontsize=20)

    # Twin axis for Ne curves
    ax2 = ax1.twinx()
    ax2.set_ylabel("Ne_diploids", fontsize=20)
    ax2.set_yscale("log")

    # Plot all Ne curves for this country
    for entry in entries:
        species = entry["species"]
        ne_df = entry["df"]

        ne_df["Year"] = reference_year - ne_df["time_years_ago"]

        ax2.plot(
            ne_df["Year"],
            ne_df["Ne_diploids"],
            linewidth=4,
            alpha=0.8,
            label=species
        )

    # Combine legends
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper right", fontsize=18)

    # Formatting
    ax1.set_xlim(reference_year, 825)
    xticks = np.arange(825, reference_year + 1, 100)
    ax1.set_xticks(xticks)
    ax1.tick_params(axis="x", rotation=45, labelsize=18)
    ax1.tick_params(axis="y", labelsize=18)
    ax2.tick_params(axis="y", labelsize=18)

    plt.tight_layout()
    plt.savefig(f"{path_prefix}/megaFauna/sa_megafauna/results/shared/human_impact/{country}_combined_Ne_population")
    plt.close()

rows = []

for species, corr_dict in correlations_by_species.items():
    row = {"Species": species}
    for var, (rho, p) in corr_dict.items():
        star = significance_stars(p)
        row[var] = f"{rho:.2f}{star}"
    rows.append(row)

rho_table = pd.DataFrame(rows)
rho_table = rho_table.rename(columns={
    "Inhabitants": "Inhabitants",
    "Cropland": "Cropland (km²)",
    "Urban": "Urban (km²)"
})
value_columns = ["Inhabitants", "Cropland (km²)", "Urban (km²)"]

import matplotlib.pyplot as plt
import numpy as np

rho_numeric = rho_table.copy()
for col in value_columns:
    rho_numeric[col] = rho_numeric[col].str.extract(r"([-+]?\d*\.\d+|\d+)").astype(float)

def cell_color(val):
    if np.isnan(val):
        return "white"
    if val > 0:
        return "#b6e3b6"   # light green
    if val < 0:
        return "#f4b6b6"   # light red
    return "white"


import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(12, len(rho_table) * 0.6 + 1))
ax.axis("off")

# Build cell colours
cell_colours = []
for i in range(len(rho_table)):
    row_colors = ["white"]  # Species column stays white
    for col in value_columns:
        row_colors.append(cell_color(rho_numeric.loc[i, col]))
    cell_colours.append(row_colors)

# Render table
tbl = ax.table(
    cellText=rho_table.values,
    colLabels=rho_table.columns,
    cellColours=cell_colours,
    loc="center",
    cellLoc="center"
)

tbl.auto_set_font_size(False)
tbl.set_fontsize(12)
tbl.scale(1, 1.5)

# Save PNG
png_path = f"{path_prefix}/megaFauna/sa_megafauna/results/shared/human_impact/species_correlation_table.png"
plt.savefig(png_path, dpi=300, bbox_inches="tight")
plt.close()
