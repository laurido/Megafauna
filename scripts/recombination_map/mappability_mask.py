import subprocess
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--fasta")
parser.add_argument("--mask-bed")
parser.add_argument("--index")
parser.add_argument("--genmap-out")
parser.add_argument("--nc-contigs")
args = parser.parse_args()

nc_contigs = set(args.nc_contigs.split(","))

def run(cmd, step):
    print(f"Running: {step}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.stdout:
        print(result.stdout)
    if result.stderr:
        print(f"STDERR: {result.stderr}", file=sys.stderr)
    if result.returncode != 0:
        print(f"ERROR in {step}", file=sys.stderr)
        sys.exit(1)

# Step 1 - Index
run(f"""
    genmap index \\
        -F {args.fasta} \\
        -I {args.index}
""", "genmap index")

# Step 2 - Compute mappability
run(f"""
    genmap map \\
        -K 100 \\
        -E 0 \\
        -I {args.index} \\
        -O {args.genmap_out}/ \\
        -bg
""", "genmap map")

# Step 3 - Extract non-uniquely mappable regions
run(f"""
    awk 'BEGIN{{OFS="\\t"}} $4 < 1.0 {{print $1, $2, $3}}' \\
        {args.genmap_out}/*.bedgraph | \\
        sort -k1,1 -k2,2n | \\
        bedtools merge -i stdin > {args.mask_bed}
""", "extract mappability mask")

# Step 4 - Filter to NC contigs only
print("Filtering to NC contigs...")
nc_only_bed = args.mask_bed.replace(".bed", "_NC_only.bed")

with open(args.mask_bed) as f_in, open(nc_only_bed, "w") as f_out:
    for line in f_in:
        if line.split("\t")[0] in nc_contigs:
            f_out.write(line)

total   = sum(1 for _ in open(args.mask_bed))
nc_only = sum(1 for _ in open(nc_only_bed))
print(f"Total intervals:   {total}")
print(f"NC-only intervals: {nc_only}")
print(f"Done: {args.mask_bed}")
print(f"Done: {nc_only_bed}")
