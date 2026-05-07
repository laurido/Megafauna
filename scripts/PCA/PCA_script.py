import allel
import numpy as np
import pandas as pd
import os
import pickle
import sys

group        = sys.argv[1]
vcf_path     = sys.argv[2]
pca_out      = sys.argv[3]
pruned_snps  = sys.argv[4]  # path to .prune.in file from ADMIXTURE step

print("Running analysis for:", group)

# Load pruned SNP IDs
pruned_ids = set(pd.read_csv(pruned_snps, header=None)[0].values)
print(f"Loaded {len(pruned_ids)} pruned SNP IDs")

# Find all VCF files
vcfs = [x for x in os.listdir(vcf_path)
        if x.endswith("fploidy_2_mploidy_2_gt_snps.vcf.gz")]

gn_list      = []
sample_ids   = []

for vcf in vcfs:
    callset = allel.read_vcf(vcf_path + vcf,
                         fields=['calldata/GT', 'variants/CHROM', 'variants/POS', 
                                 'variants/REF', 'variants/ALT', 'samples'],
                         chunk_length=50000)

    # Skip empty or unreadable VCFs
    if callset is None:
        print(f"  Warning: {vcf} returned no data, skipping")
        continue

    # Build SNP IDs the same way PLINK did: CHROM:POS:REF:ALT
    chrom  = callset['variants/CHROM']
    pos    = callset['variants/POS']
    ref    = callset['variants/REF']
    alt    = callset['variants/ALT'][:, 0]  # first alt allele
    ids    = np.array([f"{c}:{p}:{r}:{a}" for c, p, r, a in zip(chrom, pos, ref, alt)])

    # Filter to pruned SNPs
    mask_pruned = np.isin(ids, list(pruned_ids))
    print(f"  {vcf}: {mask_pruned.sum()} / {len(ids)} SNPs kept after LD filter")

    gt = allel.GenotypeArray(callset['calldata/GT'])
    gt = gt[mask_pruned]
    gn = gt.to_n_alt()
    gn_list.append(gn)

    if not sample_ids:
        sample_ids.extend(callset['samples'])

gn_combined = np.vstack(gn_list)

# Remove invariant and non-finite sites
mask_finite = np.all(np.isfinite(gn_combined), axis=1)
gn_filtered = gn_combined[mask_finite]
stds        = np.std(gn_filtered, axis=1)
gn_filtered = gn_filtered[stds > 0]

print(f"SNPs after all filters: {gn_filtered.shape[0]}")

# Random subsample if still too many
#max_snps = 100000
#if gn_filtered.shape[0] > max_snps:
#    idx         = np.random.choice(gn_filtered.shape[0], max_snps, replace=False)
#    gn_filtered = gn_filtered[idx]

# PCA
coords, model = allel.pca(
    gn_filtered,
    n_components=5,
    copy=True,
    scaler='patterson',
    ploidy=2
)

# Build output dataframe
pca_df          = pd.DataFrame(coords, columns=[f'PC{i+1}' for i in range(coords.shape[1])])
pca_df['IND_ID'] = sample_ids[:coords.shape[0]]
cols            = ['IND_ID'] + [col for col in pca_df.columns if col != 'IND_ID']
pca_df          = pca_df[cols]

pca_df.to_csv(f"{pca_out}/pca_dataset_{group}.txt", sep='\t', index=False)

with open(f"{pca_out}/pca_model_{group}.pkl", "wb") as f:
    pickle.dump(model, f)

print(f"Saved PCA results to {pca_out}")