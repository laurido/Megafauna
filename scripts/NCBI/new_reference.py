import os
if os.getcwd().startswith("/home/lakrids"):
    path_prefix = "/home/lakrids/GenomeDK"
else:
    path_prefix = "/faststorage/project"

reference_path = path_prefix + "/megaFauna/sa_megafauna/metadata/references.txt"

reference_folder = "Panthera_onca"
ref_genome_name = "mPanOnc1_haplotype_1"
genbank = "GCA_046562885.2"
ftp = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/046/562/885/GCF_046562885.1_mPanOnc1_haplotype_1/GCF_046562885.1_mPanOnc1_haplotype_1_genomic.fna.gz"

def add_reference(path, reference_folder, ref_genome_name, genbank, ftp):
    new_line = f"{reference_folder}\t{ref_genome_name}\t{genbank}\t{ftp}\n"

    # Check for duplicates
    with open(path, "r") as f:
        for line in f:
            if line.startswith(reference_folder + "\t") or genbank in line:
                return  # silently skip or raise an error if you prefer

    with open(path, "a") as f:
        f.write(new_line)

add_reference(path=reference_path, reference_folder=reference_folder, ref_genome_name=ref_genome_name, genbank=genbank, ftp=ftp)