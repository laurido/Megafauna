import pandas as pd
import matplotlib.pyplot as plt

species = "Panthera_tigris"

# Load combined recombination map
df = pd.read_csv(f"~/GenomeDK/megaFauna/sa_megafauna/results/{species}/pyrho/chrA_recomb_map_2_8.10_10_{species}.txt", 
                 sep=r"\s+", names=["chrom", "start", "end", "rate"])

# Compute cM/Mb
df["cM_Mb"] = df["rate"] * 1e8

# Midpoint of each interval for plotting
df["mid"] = (df["start"] + df["end"]) / 2

# Build a genome-wide coordinate system
chrom_offsets = {}
offset = 0
for chrom in df["chrom"].unique():
    chrom_offsets[chrom] = offset
    offset += df[df["chrom"] == chrom]["end"].max()

df["genome_pos"] = df.apply(lambda row: row["mid"] + chrom_offsets[row["chrom"]], axis=1)

# Plot
plt.figure(figsize=(18, 4))
for i, chrom in enumerate(df["chrom"].unique()):
    sub = df[df["chrom"] == chrom]
    plt.plot(sub["genome_pos"], sub["cM_Mb"], 
             linestyle="none", marker=".", markersize=2,
             label=chrom if i < 15 else None,  # avoid clutter
             color="black" if i % 2 == 0 else "gray")

plt.ylabel("Recombination rate (cM/Mb)")
plt.xlabel("Genome position (bp)")
plt.title("Genome-wide recombination rate")
plt.tight_layout()
plt.show()
