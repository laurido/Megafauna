from gwf import Workflow
import subprocess
import glob
import pandas as pd
import pickle

# Define workflow and set default options for all targets

gwf = Workflow(defaults = {
    "cores": 1,
    "memory": "8g",
    "walltime": "01:00:00",
    "account": "megaFauna"})

##########################################
#    	     --- User inputs ---
##########################################

# Prerequisites: 
# 1 Finish running workflow.py
# 2 Run PCA analysis and identify populations in file: 
#   {results}/PCA/population_{population}_{species}.txt with a sample name in each line.
# 3 Create ref/parameters_{species}.txt mutation rate and generation time in the format
#   mu = "{mutation rate}"
#   generation = "{generation time}"

species = "Panthera_tigris"
population_name = "India_NE"

# Filtering for missingness. 10 means keeping only < than 10% Missingness SNPs
filter = "30"
#filter = ""

number_of_contigs_included = 2

##########################################
#    	     --- Preparations ---
##########################################

# Conda environments initiation
general_header = '''
eval "$(conda shell.bash hook)"
conda activate megafauna
'''

# Directory paths
data = f"/faststorage/project/megaFauna/sa_megafauna/data/{species}"
results = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}"
contigs = f"{data}/VCF/contigs/miss_{filter}"
gone_in = f"{data}/analysis_input/GONE"
gone_out = f"{results}/GONE"
done_path = f"{results}/done/"

# Initializing directories
subprocess.run(["mkdir", "-p", gone_in, gone_out, contigs, done_path])

# Load regions table
df = pd.read_csv(f"{data}/ref/regions_{species}_updated.txt", sep="\t")
# Extract only NC chromosomes with ploidy 2/2
chromosomes = df[
    df["region"].str.startswith("NC") &
    (df["female_ploidy"] == 2) &
    (df["male_ploidy"] == 2)
]["region"].tolist()
# Limit number of contigs if needed
chromosomes_count = len(chromosomes)
if chromosomes_count < number_of_contigs_included or number_of_contigs_included <1:
    number_of_contigs_included = chromosomes_count

chromosomes = chromosomes[:number_of_contigs_included]
contigs_included = f"_{len(chromosomes)}"

# Read sample IDs from file for a given population
samples_file = f"{results}/PCA/population_{population_name}_{species}.txt"
with open(samples_file) as f:
    sample_list = [line.strip() for line in f if line.strip()]
# Build the final string
population_samples = ",".join(sample_list)
sample_size = len(sample_list)

##########################################
#     --- GONE file preparation ---
##########################################

# A.1 Concatenate autosome VCFs
filtered_vcf_batches = f"{data}/VCF/batches_filtered"
chrA_concat = f"{data}/VCF/chrA_concat_{species}.vcf.gz"
job_id_concat_chrA = f"concatenate_autosomes_{species}"
done = done_path + job_id_concat_chrA

gwf.target(job_id_concat_chrA, 
           inputs=[filtered_vcf_batches], 
           outputs=[chrA_concat, done]
           ) << general_header + f"""
    bcftools concat {filtered_vcf_batches}/*fploidy_2_mploidy_2* -Oz -o {chrA_concat}
    touch {done}
"""

# A.2 Choose a population for analysis and do missingness filtering
chrA_pop = f"{data}/VCF/chrA_{population_name}_{filter}_{species}.vcf.gz"
job_id_pop_and_missingness = f"population_and_missingness_filtering_{population_name}_{filter}_{species}"
prev_done = done
done = done_path + job_id_pop_and_missingness

gwf.target(job_id_pop_and_missingness, 
           inputs=[chrA_concat, prev_done], 
           outputs=[chrA_pop, f"{chrA_pop}.csi", done]
           ) << general_header + f"""
    
    bcftools view -S {samples_file} {chrA_concat} \
    | bcftools +fill-tags -- -t AC,AN,F_MISSING \
    | bcftools view -i 'INFO/AC>0 && INFO/AN>0' \
    | bcftools view -i 'INFO/F_MISSING < {int(filter)/100}' \
    | sed 's/|/\\//g' \
    | bcftools view -Oz -o {chrA_pop}

    bcftools index {chrA_pop}

    touch {done}
"""

##########################################
#    	     --- GONE ---
##########################################

# B.2 - estimate population size
job_id_gone = f"GONE_{population_name}{contigs_included}_{filter}_{species}"
gone_estimate = f"{gone_out}/{job_id_gone}"
prev_done = done
done = done_path + job_id_gone

gwf.target(job_id_gone, 
           inputs=[chrA_pop, done], 
           outputs=[gone_estimate, done], 
            cores=8, memory="64g", 
            walltime="03:00:00", account="megaFauna"
            ) << f"""
    gone2 -r {chrA_pop} -o 
    touch {done}
"""