import re
import pandas as pd
import os
import matplotlib.pyplot as plt

k_overrides = {
    "Boselaphus_tragocamelus": 1,
    "Rhinoceros_unicornis": 1,
    "Panthera_tigris": 3,
    "Loxodonta_africana": 3,
    "Loxodonta_cyclotis": 2,
    "Elephas_maximus": 2,
    "Panthera_leo": 3,
    "Panthera_pardus": 2,
    "Panthera_uncia": 2,
    "Ceratotherium_simum": 2,
    "Diceros_bicornis": 2}

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

fst_records = []

for i in range(species_and_refs.shape[0]):
    group = species_and_refs.FOLDER[i]
    k = k_overrides.get(group)
    if k == 1:
        continue

    log_file = f"{path_prefix}/megaFauna/sa_megafauna/results/{group}/ADMIXTURE/admixture_K{k}.log"
    
    with open(log_file) as f:
        content = f.read()
            
    # Extract Fst matrix
    lines = content.split('\n')
    fst_start = next((i for i, l in enumerate(lines) if 'Fst divergences' in l), None)
    if fst_start is None:
        continue
    
    # Parse population headers and Fst values
    pop_lines = []
    for line in lines[fst_start+1:]:
        if line.startswith('CV error') or line.strip() == '':
            break
        pop_lines.append(line)
    
    # Extract pairwise Fst values
    for line in pop_lines[1:]:  # skip header row
        parts = line.strip().split('\t')
        if len(parts) < 2:
            continue
        row_pop = parts[0]
        for col_idx, val in enumerate(parts[1:]):
            if val.strip():
                col_pop = f'Pop{col_idx}'
                fst_records.append({
                    'species': group,
                    'K': k,
                    'pair': f'{col_pop} x {row_pop}',
                    'Fst': float(val.strip())
                })

# Build table
fst_df = pd.DataFrame(fst_records)
print(fst_df.head())
print(fst_df.columns.tolist())
print(fst_df.dtypes)
fst_pivot = fst_df.pivot_table(index='species', columns='pair', values='Fst', aggfunc='first')
fst_pivot = fst_pivot.reset_index()

print(fst_pivot)

fig, ax = plt.subplots(figsize=(max(6, len(fst_pivot.columns) * 1.5), max(3, len(fst_pivot) * 0.5 + 1)))
ax.axis('off')

# Format species names as L. africana
def format_species(name):
    parts = name.split('_')
    return f"{parts[0][0]}. {' '.join(parts[1:])}"

fst_display = fst_pivot.copy()
fst_display['species'] = fst_display['species'].apply(format_species)

# Round Fst values
for col in fst_display.columns:
    if col != 'species':
        fst_display[col] = fst_display[col].apply(lambda x: f'{x:.3f}' if pd.notna(x) else '-')

table = ax.table(
    cellText=fst_display.values,
    colLabels=fst_display.columns,
    cellLoc='center',
    loc='center'
)

table.auto_set_font_size(False)
table.set_fontsize(10)
table.auto_set_column_width(col=list(range(len(fst_display.columns))))

# Style header row
for j in range(len(fst_display.columns)):
    table[0, j].set_facecolor('#2c7bb6')
    table[0, j].set_text_props(color='white', fontweight='bold')

# Alternate row colors
for row_idx in range(1, len(fst_display) + 1):
    table[row_idx, 0].set_text_props(fontstyle='italic')
    color = '#f0f0f0' if row_idx % 2 == 0 else 'white'
    for j in range(len(fst_display.columns)):
        table[row_idx, j].set_facecolor(color)

plt.tight_layout()
os.makedirs(f"{path_prefix}/megaFauna/sa_megafauna/results/shared/", exist_ok=True)
plt.savefig(f"{path_prefix}/megaFauna/sa_megafauna/results/shared/ADMIXTURE/admixture_fst_table.png",
            bbox_inches='tight', dpi=150)
plt.close()