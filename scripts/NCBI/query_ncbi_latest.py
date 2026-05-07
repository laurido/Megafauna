import argparse
import os
import time
from Bio import Entrez
import xml.etree.ElementTree as ET
import pandas as pd
import xml

def search_sra_and_download_summary(genus, output_dir="."):
    # Set your email for NCBI's reference
    Entrez.email = "laurids@math.au.dk"  # Replace with your email
    
    print(f"Searching NCBI SRA for genus: {genus}")
    
    # Construct search term with all filters
    
    search_term = f' ("{genus}"[Organism] AND ("biomol dna"[Properties] AND "strategy wgs"[Properties])'
    search_term = f' ("{genus}"[Organism])'

    try:
        # First, search to get the IDs
        print("Performing search with filters...")
        search_handle = Entrez.esearch(db="sra", term=search_term, retmax=99999) # retmax is the max queries you want to note 
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        # Get the IDs from the search
        id_list = search_results["IdList"]
        
        if not id_list:
            print(f"No results found for {genus} with the specified filters.")
            return
        
        print(f"Found {len(id_list)} results. Fetching summaries...")
        
        # Fetch the summaries for these IDs
        fetch_handle = Entrez.efetch(db="sra", id=",".join(id_list), rettype="summary", retmode="text")
        summary_data = fetch_handle.read()
        fetch_handle.close()
        # Create the output file path
        output_file = os.path.join(output_dir, f"{genus}_SRA_summary.txt")
        
        # Write the summary to file - fix: don't encode since summary_data is already bytes
        with open(output_file, 'wb') as f:
            f.write(summary_data)  # Removed the .encode('utf-8')
        
        print(f"Successfully downloaded summary to: {output_file}")
        print(f"Found {len(id_list)} SRA entries for {genus} with applied filters:")
        print("- SOURCE=DNA")
        print("- TYPE=GENOME") 
        print("- LIBRARY LAYOUT=PAIRED")
        print("- STRATEGY=GENOME")
        
    except Exception as e:
        print(f"Error: {str(e)}")
        import traceback
        traceback.print_exc()

# Function to print XML content
def print_xml_content(element, indent=0):
    print(" " * indent + f"<{element.tag}> {element.text.strip() if element.text else ''}")
    for child in element:
        print_xml_content(child, indent + 2)
    print(" " * indent + f"</{element.tag}>")

def parse_sra_xml(xml_file):
    # Parse XML file
    tree = ET.parse(xml_file)

    root = tree.getroot()

    #print_xml_content(root)
    
    # Lists to store data
    records = []
    
    # Process each experiment package
    for exp_package in root.findall('./EXPERIMENT_PACKAGE'):
        # Extract experiment info
        exp = exp_package.find('./EXPERIMENT')
        if exp is None:
            continue
        
        exp_accession = exp.get('accession')
        exp_title = exp.find('./TITLE').text if exp.find('./TITLE') is not None else None

        # Get library info
        lib_info = exp.find('.//LIBRARY_DESCRIPTOR')
        lib_strategy = lib_info.find('./LIBRARY_STRATEGY').text if lib_info is not None and lib_info.find('./LIBRARY_STRATEGY') is not None else None
        lib_source = lib_info.find('./LIBRARY_SOURCE').text if lib_info is not None and lib_info.find('./LIBRARY_SOURCE') is not None else None
        lib_selection = lib_info.find('./LIBRARY_SELECTION').text if lib_info is not None and lib_info.find('./LIBRARY_SELECTION') is not None else None
        library_layout = lib_info.find("LIBRARY_LAYOUT")
        
        if library_layout is not None:
            layout_type = list(library_layout)[0].tag  # Extracts SINGLE or PAIRED
        else:
            layout_type = "Unknown"

        lib_layout = lib_info.find('./LIBRARY_LAYOUT/PAIRED').text if lib_info is not None and lib_info.find('./LIBRARY_LAYOUT/PAIRED') is not None else None
        
        # Extract platform/instrument info
        instrument = exp.find('.//INSTRUMENT_MODEL')
        instrument_model = instrument.text if instrument is not None else None
        
        # Extract study information
        study = exp_package.find('./STUDY')
        study_title = None
        bioproject_id = None
        
        if study is not None:
            study_title_elem = study.find('./DESCRIPTOR/STUDY_TITLE')
            if study_title_elem is not None:
                study_title = study_title_elem.text
            # If not found in direct path, try another common path
            elif study.find('./STUDY_TITLE') is not None:
                study_title = study.find('./STUDY_TITLE').text
                
            # Extract BioProject ID
            # First look for external_id with "BioProject" namespace
            for ext_id in study.findall('.//EXTERNAL_ID'):
                if ext_id.get('namespace') == 'BioProject' and ext_id.text is not None:
                    bioproject_id = ext_id.text
                    break
            
            # If not found, try to find it in study links
            if bioproject_id is None:
                for db_link in study.findall('.//DB_XREF'):
                    if db_link.get('db') == 'BioProject' and db_link.text is not None:
                        bioproject_id = db_link.text
                        break
        
        # Handle case where study info might be in a different location
        if study_title is None and exp_package.find('.//STUDY_TITLE') is not None:
            study_title = exp_package.find('.//STUDY_TITLE').text
            
        # Look for BioProject ID in any location if not found yet
        if bioproject_id is None:
            for ext_id in exp_package.findall('.//EXTERNAL_ID'):
                if ext_id.get('namespace') == 'BioProject' and ext_id.text is not None:
                    bioproject_id = ext_id.text
                    break
            
            # Also check for a db_xref
            if bioproject_id is None:
                for db_link in exp_package.findall('.//DB_XREF'):
                    if db_link.get('db') == 'BioProject' and db_link.text is not None:
                        bioproject_id = db_link.text
                        break
        
        # Extract sample info
        sample = exp_package.find('./SAMPLE')
        sample_name = None
        sex = 'not recorded'
        
        if sample is not None:
            sample_accession = sample.get('accession')
            sample_title = sample.find('./TITLE').text if sample.find('./TITLE') is not None else None
            
            # Get BioSample ID
            biosample_id = None
            
            # Look for BioSample ID in external_id with "BioSample" namespace
            for ext_id in sample.findall('.//EXTERNAL_ID'):
                if ext_id.get('namespace') == 'BioSample' and ext_id.text is not None:
                    biosample_id = ext_id.text
                    break
            
            # Also check for a db_xref if not found
            if biosample_id is None:
                for db_link in sample.findall('.//DB_XREF'):
                    if db_link.get('db') == 'BioSample' and db_link.text is not None:
                        biosample_id = db_link.text
                        break
                        
            # If still not found, check if the sample accession itself is a BioSample ID
            if biosample_id is None and sample_accession and sample_accession.startswith('SAMN'):
                biosample_id = sample_accession
            
            # Get scientific name
            sci_name_elem = sample.find('.//SCIENTIFIC_NAME')
            scientific_name = sci_name_elem.text if sci_name_elem is not None else None
            
            # Get taxon ID
            taxon_elem = sample.find('.//TAXON_ID')
            taxon_id = taxon_elem.text if taxon_elem is not None else None
            
            # Extract sample_name - check multiple locations where it might be stored
            # First check in SAMPLE_NAME element
            sample_name_elem = sample.find('.//SAMPLE_NAME/ANONYMIZED_NAME')
            if sample_name_elem is not None and sample_name_elem.text:
                sample_name = sample_name_elem.text
            
            # If not found, check in SAMPLE_ATTRIBUTES
            if sample_name is None:
                for attr in sample.findall('.//SAMPLE_ATTRIBUTE'):
                    tag = attr.find('./TAG')
                    value = attr.find('./VALUE')
                    if tag is not None and tag.text and value is not None and value.text:
                        if tag.text.lower() == 'sample_name' or tag.text.lower() == 'sample name':
                            sample_name = value.text
                            break
            
            # Extract sex information from SAMPLE_ATTRIBUTES
            for attr in sample.findall('.//SAMPLE_ATTRIBUTE'):
                tag = attr.find('./TAG')
                value = attr.find('./VALUE')
                if tag is not None and tag.text and value is not None and value.text:
                    if tag.text.lower() == 'sex' or tag.text.lower() == 'gender':
                        sex = value.text
                        break
        else:
            sample_accession = sample_title = scientific_name = taxon_id = biosample_id = None
        
        # Process runs for this experiment
        runs = exp_package.findall('.//RUN')
        for run in runs:
            run_accession = run.get('accession')
            run_published = run.get('published')
            total_spots = run.get('total_spots')
            total_bases = run.get('total_bases')
            size_bytes = run.get('size')
            
            # Add to records
            record = {
                'experiment_accession': exp_accession,
                'experiment_title': exp_title,
                'study_title': study_title,
                'bioproject_id': bioproject_id,
                'sample_accession': sample_accession,
                'biosample_id': biosample_id,  # Added BioSample ID
                'sample_title': sample_title,
                'sample_name': sample_name,  # Added sample_name
                'sex': sex,  # Added sex information
                'scientific_name': scientific_name,
                'taxon_id': taxon_id,
                'run_accession': run_accession,
                'instrument_model': instrument_model,
                'library_strategy': lib_strategy,
                'library_source': lib_source,
                'library_selection': lib_selection,
                'library_layout': layout_type,
                'total_spots': total_spots,
                'total_bases': total_bases,
                'size_bytes': size_bytes,
                'published_date': run_published
            }
            records.append(record)
    
    # Create DataFrame
    df = pd.DataFrame(records)
    
    # Convert numeric columns
    numeric_cols = ['total_spots', 'total_bases', 'size_bytes', 'taxon_id']
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    
    # Convert date column
    if 'published_date' in df.columns:
        df['published_date'] = pd.to_datetime(df['published_date'], errors='coerce')
        
    # Add some derived columns
    if 'total_bases' in df.columns and 'total_spots' in df.columns:
        df['avg_read_length'] = df['total_bases'] / df['total_spots'] / 2  # Paired-end, so divide by 2
        df['gb_size'] = df['total_bases'] / 1e9  # Convert to gigabases
    
    # Add this after creating the DataFrame
    print(f"Found {len(df)} total runs")
    print(f"Number of unique experiment accessions: {df['experiment_accession'].nunique()}")
    print(f"Number of unique run accessions: {df['run_accession'].nunique()}")
    print(f"Number of unique BioProject IDs: {df['bioproject_id'].nunique()}")
    print(f"Number of unique BioSample IDs: {df['biosample_id'].nunique()}")
    print(f"Number of samples with sex information: {(df['sex'] != 'not recorded').sum()}")
    print(f"Number of samples with sample_name: {df['sample_name'].notna().sum()}")
    
    return df

Genus_list = ["Diceros", "Ceratotherium"]


rerun = True

#samples = pd.read_csv('Samples_currently_in.csv', sep=',')
for genus in Genus_list:
    if rerun:
        search_result = search_sra_and_download_summary(genus, '/faststorage/project/megaFauna/people/laurids/new_samples')
        try:
            df = parse_sra_xml(f'/faststorage/project/megaFauna/people/laurids/new_samples/{genus}_SRA_summary.txt')
            df.to_csv(f'/faststorage/project/megaFauna/people/laurids/new_samples/{genus}_all.csv')
        except:
            continue
        #filter
        #filtered_df = df[~df['sample_accession'].isin(samples['ENA-SRA_SAMPLE_ID'])]
        #filtered_df.to_csv(f'new_samples/{genus}_potential_new_samples.csv')
