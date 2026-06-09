import subprocess
import argparse
import pandas as pd
from pathlib import Path
import sys
import matplotlib.pyplot as plt
import os

parser = argparse.ArgumentParser()
parser.add_argument("--sample-beds", nargs="+")
parser.add_argument("--roh-bed")
parser.add_argument("--mappability-bed")
parser.add_argument("--cov-bed")
parser.add_argument("--rep-bed")
parser.add_argument("--ass-bed")
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
rows["assembly"]  = masked_bases_per_chrom(args.ass_bed)
rows["ROH"]  = masked_bases_per_chrom(args.roh_bed)
rows["repeats"]  = masked_bases_per_chrom(args.rep_bed)
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

def compute_callable_lengths(bed_path, bgzipped=True):
    """Return list of callable region lengths (gaps between masked intervals)."""
    
    if bgzipped:
        cmd = f"bgzip -d -c {bed_path}"
    else:
        cmd = f"cat {bed_path}"

    awk_script = r"""
    BEGIN {prev_chr=""; prev_end=0}
    {
        if ($1 == prev_chr) {
            gap = $2 - prev_end;
            if (gap > 0) print gap;
        }
        prev_chr = $1;
        prev_end = $3;
    }
    """

    result = subprocess.run(
        f"{cmd} | sort -k1,1 -k2,2n | awk '{awk_script}'",
        shell=True, capture_output=True, text=True
    )

    if result.returncode != 0 or not result.stdout.strip():
        print("WARNING: no callable regions extracted", file=sys.stderr)
        print(result.stderr, file=sys.stderr)
        return []

    lengths = [int(x) for x in result.stdout.strip().split("\n") if x]
    return lengths

def build_histogram(lengths, bin_size=10000):
    """Return histogram as dict: {bin_start: count}"""
    hist = {}
    for l in lengths:
        b = (l // bin_size) * bin_size
        hist[b] = hist.get(b, 0) + 1
    return hist

import math

def build_log_histogram(lengths):
    """Return log-scale histogram"""
    hist = {}
    for l in lengths:
        if l > 0:
            b = int(math.log10(l))
            hist[b] = hist.get(b, 0) + 1
    return hist

# =========================================================
# Callable region histogram (from final mask)
# =========================================================

lengths = compute_callable_lengths(args.final_bed, bgzipped=True)

if lengths:
    print(f"\nComputed {len(lengths)} callable regions")

    # ---- standard histogram ----
    hist = build_histogram(lengths, bin_size=10000)

    hist_df = pd.DataFrame([
        {"bin_start": k, "count": v}
        for k, v in sorted(hist.items())
    ])

    hist_path = args.out + ".hist.tsv"
    hist_df.to_csv(hist_path, sep="\t", index=False)
    print(f"Histogram written to {hist_path}")

    # ---- log histogram ----
    log_hist = build_log_histogram(lengths)

    log_df = pd.DataFrame([
        {"log10_bin": f"1e{int(k)}", "count": v}
        for k, v in sorted(log_hist.items())
    ])

    log_path = args.out + ".loghist.tsv"
    log_df.to_csv(log_path, sep="\t", index=False)
    print(f"Log histogram written to {log_path}")

    # ---- quick summary ----
    print("\nCallable region summary:")
    print(f"≥10 kb:  {sum(1 for x in lengths if x >= 10000)}")
    print(f"≥50 kb:  {sum(1 for x in lengths if x >= 50000)}")
    print(f"≥100 kb: {sum(1 for x in lengths if x >= 100000)}")
    print(f"Max:     {max(lengths)} bp")

def get_group_name(path):
    # parent of 'smcpp' directory
    parts = Path(path).parts
    if "results" in parts:
        idx = parts.index("results")
        return parts[idx+1]
    return Path(path).stem

def plot_histogram(lengths, group, out_prefix):
    """Plot histogram of callable region sizes."""

    if not lengths:
        print("No lengths to plot", file=sys.stderr)
        return

    # Convert to numpy-like list
    lengths_sorted = sorted(lengths)

    # Create histogram
    plt.figure(figsize=(8,6))

    plt.hist(lengths_sorted, bins=100, log=True)

    plt.xlabel("Callable region length (bp)")
    plt.ylabel("Count (log scale)")
    plt.title(f"Callable region size distribution: {group}")

    plt.xscale("log")

    plt.tight_layout()

    out_file = f"{out_prefix}.{group}.hist.png"
    plt.savefig(out_file, dpi=300)
    plt.close()

    print(f"Histogram plot saved to {out_file}")

# =========================================================
# Callable region histogram + plot
# =========================================================

group = get_group_name(args.final_bed)

lengths = compute_callable_lengths(args.final_bed, bgzipped=True)

if lengths:
    print(f"\nComputed {len(lengths)} callable regions")

    # Save histogram table (optional)
    hist = build_histogram(lengths, bin_size=10000)
    hist_df = pd.DataFrame([
        {"bin_start": k, "count": v}
        for k, v in sorted(hist.items())
    ])
    hist_df.to_csv(args.out + f".{group}.hist.tsv", sep="\t", index=False)

    # ✅ Plot
    plot_histogram(lengths, group, args.out)