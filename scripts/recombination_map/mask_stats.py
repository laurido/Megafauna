import subprocess
import argparse
import pandas as pd
from pathlib import Path
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--sample-beds", nargs="+")
parser.add_argument("--mappability-bed")
parser.add_argument("--cov-bed")
parser.add_argument("--final-bed")  # bgzipped
parser.add_argument("--out")
args = parser.parse_args()

def masked_bases_per_chrom(bed_path, bgzipped=False):
    """Return dict of {chrom: masked_bases} for a bed file."""
    if bgzipped:
        cmd = f"bgzip -d -c {bed_path}"
    else:
        cmd = f"cat {bed_path}"

    result = subprocess.run(
        f"{cmd} | awk '{{sum[$1] += $3-$2}} END {{for (c in sum) print c, sum[c]}}'",
        shell=True, capture_output=True, text=True
    )
    
    # Debug — print stderr if something went wrong
    if result.returncode != 0 or not result.stdout.strip():
        print(f"WARNING: no output for {bed_path}", file=sys.stderr)
        print(f"STDERR: {result.stderr}", file=sys.stderr)
        return {}
    
    chrom_counts = {}
    for line in result.stdout.strip().split("\n"):
        if line:
            chrom, count = line.split()
            chrom_counts[chrom] = int(count)
    return chrom_counts

# =========================================================
# Collect all data
# =========================================================
rows = {}

# Per sample rows
for bed in args.sample_beds:
    name = Path(bed).stem  # removes .bed extension
    # Extract sample ID - adjust the pattern to match your naming
    # e.g. "merge_mask_Panthera_tigris_SAMN20424160" -> "SAMN20424160"
    sample = name.split("_")[-1]
    rows[sample] = masked_bases_per_chrom(bed)


# Species level rows
rows["mappability"]  = masked_bases_per_chrom(args.mappability_bed)
rows["coverage"]     = masked_bases_per_chrom(args.cov_bed)
rows["final_merged"] = masked_bases_per_chrom(args.final_bed, bgzipped=True)

# =========================================================
# Build dataframe and write
# =========================================================
df = pd.DataFrame(rows).T  # rows = samples/masks, cols = chroms
df.index.name = "sample"

# Sort columns by chrom name
df = df.reindex(sorted(df.columns), axis=1)
df = df.fillna(0).astype(int)

df.to_csv(args.out, sep="\t")
print(f"Stats written to {args.out}")
print(df)
