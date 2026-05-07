from Bio import Entrez
import xml.etree.ElementTree as ET
import pandas as pd
import os
from geopy.geocoders import Nominatim
import folium
import seaborn as sns

Entrez.email = "laurids@bio.au.dk"

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

# Fetching metadata from the NCBI database
# Have a look at the added columns or go to the sample pages to decide how to make the geography column.
def fetch_selected_attributes(sample_id, target_attrs):
    results = {attr: None for attr in target_attrs}

    try:
        # Fetch BioSample XML
        handle = Entrez.efetch(db="biosample", id=sample_id, rettype="xml")
        xml_data = handle.read()
        root = ET.fromstring(xml_data)

        # Extract target attributes
        for attr_elem in root.findall(".//Attribute"):
            name = attr_elem.attrib.get("harmonized_name")
            if name in target_attrs:
                results[name] = attr_elem.text

    except Exception as e:
        print(f"Error fetching {sample_id}: {e}")

    return results

target_attrs = {"ecotype", "biomaterial_provider", "geo_loc_name"}


for i in range(species_and_refs.shape[0]):
    group      = species_and_refs.FOLDER[i]
    ref_folder = species_and_refs.REFERENCE_FOLDER[i]

    pca_df = pd.read_csv(f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/PCA/pca_dataset_{group}.txt", sep='\t')

    pca_df[list(target_attrs)] = pca_df["IND_ID"].apply(
    lambda sid: pd.Series(fetch_selected_attributes(sid, target_attrs)))

    pca_df.to_csv(f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/PCA/pca_dataset_{group}.txt", sep='\t', index=False)