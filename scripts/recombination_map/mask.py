import subprocess
import sys
import os
import pandas as pd

sample = sys.argv[1]
chrom = sys.argv[2]
bed = sys.argv[3]
cov_min, min_het_ad, gq_min = sys.argv[4].split(",")
ref_folder = sys.argv[5]
cov_min    = int(cov_min)
min_het_ad = int(min_het_ad)
gq_min     = int(gq_min)

species = bed.split("/")[6]

df_regions = pd.read_csv(f"/faststorage/project/megaFauna/sa_megafauna/data/{ref_folder}/ref/regions_{ref_folder}_updated.txt", sep="\t")
df_cov = pd.read_table(f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/cov/{sample}.cov", header=None)

batch   = int(df_regions[df_regions.region == chrom]["batch"].iat[0])
cov     = float(df_cov[df_cov[0] == chrom][6].iat[0])
cov_max = cov * 2
cov_min = max(cov / 2, cov_min)
gvcf    = f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/gVCF/{species}_batch_{batch}_fploidy_2_mploidy_2_gt.gvcf.gz"

os.makedirs(os.path.dirname(bed), exist_ok=True)

print(f"Sample: {sample} | Chrom: {chrom} | Batch: {batch} | "
      f"min_dp: {cov_min} | max_dp: {cov_max:.1f} | GQ: {gq_min} | HET_AD: {min_het_ad}")

bash_cmd = f"""
    bcftools view -s {sample} -r {chrom} {gvcf} | \\
    bcftools filter \\
        -i "(GT='./.') | \\
            (GT='het' & FMT/AD[*:*] < {min_het_ad}) | \\
            FMT/DP < {cov_min} | \\
            FMT/DP > {cov_max} | \\
            (FMT/GQ != '.' & FMT/GQ < {gq_min})" | \\
    grep -v '^#' | \\
    awk 'BEGIN{{OFS="\\t"}} {{print $1, $2-1, $2}}' | \\
    sort -k1,1 -k2,2n | \\
    bedtools merge > {bed}
"""

result = subprocess.run(
    bash_cmd,
    shell=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    text=True
)

if result.stdout:
    print(result.stdout)
if result.stderr:
    print(f"STDERR: {result.stderr}", file=sys.stderr)
if result.returncode != 0:
    sys.exit(1)
# | \ (FORMAT/RGQ[0] != '.' & FORMAT/RGQ[0] < {gq_min})