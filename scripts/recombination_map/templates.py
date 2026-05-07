from gwf import AnonymousTarget, Workflow
from gwf.executors import Conda

# Define workflow and set default options for all targets

gwf = Workflow()
default_options = {
    "cores": 1,
    "memory": "8g",
    "walltime": "06:00:00",
    "account": "megaFauna",
}

##########################################
#     --- smc++ file preparation ---
##########################################


# A.1 Concatenate autosome VCFs
def vcfconcat(vcfs, chrA_vcf, done):
    inputs = []
    outputs = [done]
    options = default_options.copy()
    executor = Conda("megafauna")
    spec = f"""
    bcftools concat {" ".join(vcfs)} -Oz -o {chrA_vcf}
    bcftools index {chrA_vcf}
    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=executor)


# A.2 Choose a population for analysis and do missingness filtering
def subset_and_filter(vcf_in, samples, filter, pop_vcf, done_prev, done):
    inputs = [done_prev]
    outputs = [done]
    options = default_options.copy()
    executor = Conda("megafauna")
    spec = """
    bcftools view -S {samples} {vcf_in} \
    | bcftools +fill-tags -- -t AC,AN,F_MISSING \
    | bcftools view -i 'INFO/AC>0 && INFO/AN>0 && INFO/F_MISSING < {f_missing}' \
    | sed 's/|/\\//g' \
    | bcftools view -Oz -o {pop_vcf}

    bcftools index {pop_vcf}

    touch {done}
    """.format(vcf_in=vcf_in, samples=samples, filter=int(filter) / 100, pop_vcf=pop_vcf, done=done)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=executor)

# A.2 Choose a population for analysis and do missingness filtering
def subset_and_filter(vcf_in, samples, filter, pop_vcf, done_prev, done):
    inputs = [done_prev]
    outputs = [done]
    options = default_options.copy()
    executor = Conda("megafauna")
    spec = """
    bcftools view {vcf_in} \
    | bcftools view -i 'INFO/AC>0 && INFO/AN>0' \
    | sed 's/|/\\//g' \
    | bcftools view -Oz -o {pop_vcf}

    bcftools index {pop_vcf}

    touch {done}
    """.format(vcf_in=vcf_in, filter=int(filter) / 100, pop_vcf=pop_vcf, done=done)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=executor)


# A.3 - Mask bed
def mask_beds(sample, chrom, bed, parameters, done_prev, done):
    inputs = [done_prev]
    outputs = [done]
    options = {"cores" : 1, 'memory': "8g", 'walltime': "06:00:00", 'account': "megaFauna"}
    executor = Conda("megafauna")
    spec = f"""
    python /faststorage/project/megaFauna/people/laurids/scripts/recombination_map/mask.py {sample} {chrom} {bed} {parameters}

    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=executor)

def sample_merge_mask(beds, merged_bed, done_prev, done):
    inputs = done_prev
    outputs = [done]
    options = default_options.copy()
    executor = Conda("megafauna")
    spec = f"""
    # Merge all per-contig mask files into one population mask
    sort -k1,1 -k2,2n {" ".join(beds)} | \
       bedtools merge -i stdin > {merged_bed}
    rm {" ".join(beds)}
       
    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=executor)

def merge_mask(beds, merged_bed, done_prev, done):
    inputs = done_prev
    outputs = [done]
    options = default_options.copy()
    executor = Conda("megafauna")
    spec = f"""
    # Merge all per-contig mask files into one population mask
    sort -k1,1 -k2,2n {" ".join(beds)} | \
       bedtools merge -i stdin > {merged_bed}
       
    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=executor)

def mappability_mask(fasta, mask_bed, index, genmap_out, done):
    inputs  = [fasta]
    outputs = [done]
    options = {"memory": "64g", "cores": 8, "walltime": "10:00:00", 'account': "megaFauna"}
    executor = Conda("megafauna")
    spec = f"""
    mkdir -p {genmap_out}
    # Step 1 - Index
    genmap index \\
        -F {fasta} \\
        -I {index}

    # Step 2 - Compute mappability
    genmap map \\
        -K 100 \\
        -E 1 \\
        -T 8 \\
        -I {index} \\
        -O {genmap_out}/ \\
        --bg 

    # Step 3 - Extract non-uniquely mappable regions
    awk 'BEGIN{{OFS="\\t"}} $4 < 1.0 {{print $1, $2, $3}}' \\
        {genmap_out}/*.bedgraph | \\
        sort -k1,1 -k2,2n | \\
        bedtools merge -i stdin > {mask_bed}

    touch {done}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=executor)

#def mappability_mask(fasta, mask_bed, index, genmap_out, nc_contigs, done):
#    inputs  = [fasta]
#    outputs = [done]
#    options = {"memory": "32g", "cores": 1, "walltime": "10:00:00", "account": "megaFauna"}
#    executor = Conda("megafauna")
#    nc_contigs_str = ",".join(nc_contigs)
#    spec = f"""
#    python mappability_mask.py \\
#        --fasta       {fasta} \\
#        --mask-bed    {mask_bed} \\
#        --index       {index} \\
#        --genmap-out  {genmap_out} \\
#        --nc-contigs  {nc_contigs_str}
#
#    touch {done}
#"""
#    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=executor)

def final_merge_mask(map_bed, cov_bed, final_bed, done_prev, done):
    inputs = done_prev
    outputs = [done]
    options = default_options.copy()
    executor = Conda("megafauna")
    spec = f"""
    # Merge mappability bed and coverage bed into one
    cat {map_bed} {cov_bed} | \
        sort -k1,1 -k2,2n | \
        bedtools merge -i stdin > {final_bed}

    # Compress and index
    bgzip -f {final_bed}
    tabix -p bed {final_bed}.gz

    touch {done}

    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=executor)

def mask_stats(sample_beds, mappability_bed, cov_bed_merged, final_bed, stats_file, done_prev, done):
    inputs  = [done_prev]
    outputs = [done]
    options = default_options.copy()
    executor = Conda("megafauna")
    beds_str = " ".join(sample_beds)
    spec = f"""
    python mask_stats.py \
        --sample-beds {beds_str} \
        --mappability-bed {mappability_bed} \
        --cov-bed {cov_bed_merged} \
        --final-bed {final_bed} \
        --out {stats_file}

    touch {done}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=executor)

##########################################
#    	     --- smc++ ---
##########################################

# B.1 - convert vcf file to smc files


def vcf2smc(vcf_in, mask, chrom, pop, smc_file, done_prev, done):
    inputs = done_prev
    outputs = [done]
    options = default_options.copy()
    executor = Conda("smcpp")
    spec = f"""
    smc++ vcf2smc {vcf_in} {smc_file} \
    {chrom} {pop} -m {mask}

    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=executor)
#-m {mask}

# B.2 - estimate population size
def smcpp_estimate(smc_files, mu, estimate_name, outdir, done_prev, done):
    inputs = [done_prev]
    outputs = [done]
    options = {"memory": "32g", "cores":  8, "walltime": "10:00:00", 'account': "megaFauna"}
    executor = Conda("smcpp")
    spec = f"""
    smc++ estimate --base {estimate_name} --em-iterations 30 --cores 8 \
    {mu} {" ".join(smc_files)} -o {outdir}

    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=executor)
#--knots 26
#--knots 26 --timepoints 40 90000
 

# B.3 - create plot of population trajectory
def smcpp_plot(estimate_json, plot_name, done_prev, done):
    inputs = [done_prev]
    outputs = [done]
    options = default_options.copy()
    executor = Conda("smcpp")
    spec = f"""
    smc++ plot --csv {plot_name} {estimate_json}

    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=executor)


##########################################
# 	  --- pyrho steps ---
##########################################
# A.3 Create a VCF for each contig
def make_contig_files(vcf_in, chrom, contig_vcf, done_prev, done):
    inputs = [done_prev]
    outputs = [done]
    options = default_options.copy()
    executor = Conda("megafauna")
    spec = f"""
    bcftools view -r {chrom} {vcf_in} -Oz -o {contig_vcf}
    bcftools index {contig_vcf}

    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=executor)


# C.1 - make pyrho lookup table
def pyrho_lookup(estimate_csv, sample_size, mu, pyrho_table, done_prev, done):
    inputs = [done_prev]
    outputs = [done]
    options = {"memory": "32g", "cores":  1, "walltime": "12:00:00", 'account': "megaFauna"}
    executor = Conda("pyrho")
    spec = f"""
    pyrho make_table --samplesize {sample_size} \
    --mu {mu} --outfile {pyrho_table} --smcpp_file {estimate_csv}

    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=executor)


# C.2 - adjust pyrho hyperparameters
def pyrho_hyperparam(estimate_csv, sample_size, pyrho_table, mu, hyperparam_file, done_prev, done):
    inputs = [done_prev]
    outputs = [done]
    options = {"memory": "32g", "cores":  1, "walltime": "12:00:00", 'account': "megaFauna"}
    executor = Conda("pyrho")
    spec = f"""
    pyrho hyperparam -n {sample_size} --tablefile {pyrho_table} \
    --mu {mu} --ploidy 2 \
    --smcpp_file {estimate_csv} --outfile {hyperparam_file}

    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=executor)

    # C.3 - pyrho recombination map estimation

def pyrho_optimize(contig_vcf, pyrho_table, rmap_out, done_prev, done):
    inputs = [done_prev]
    outputs = [done]
    options = default_options.copy()
    executor = Conda("pyrho")
    spec = f"""
    pyrho optimize --vcffile {contig_vcf} --windowsize 50 --blockpenalty 50 \
    --tablefile {pyrho_table} --ploidy 2 --outfile {rmap_out}

    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=executor)


# C.4 - compute r2 from the recombination map
def pyrho_compute(pyrho_table, r2_out, done_prev, done):
    inputs = [done_prev]
    outputs = [done]
    options = default_options.copy()
    executor = Conda("pyrho")
    spec = f"""
    pyrho compute_r2 --quantiles .25,.5,.75 --compute_mean --samplesize 14 \
    --tablefile {pyrho_table} --outfile {r2_out}

    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=executor)


# C.5 - combine recombination maps and convert to PLINK format
def combine_maps(rmaps, combined_map, done_prev, done):
    inputs = [done_prev]
    outputs = [done]
    options = default_options.copy()
    executor = Conda("megafauna")
    spec = f"""
        # Combine all rmap files
        for f in {" ".join(rmaps)}; do
            chr=$(basename "$f" .rmap)
            awk -v c="$chr" '{{print c, $0}}' "$f"
        done | sort -V -k1,1 -k2,2n > {combined_map}

        touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=executor)


def plink_map(combined_map, plink_map, done_prev, done):
    inputs = [done_prev]
    outputs = [done]
    options = default_options.copy()
    executor = Conda("megafauna")
    spec = f"""
        # Convert to PLINK map format (cM)
        awk '
        BEGIN {{cum=0}}
        {{
            chrom=$1; start=$2; end=$3; rate=$4
            inc=(end-start)*rate
            cum+=inc
            snp_id=chrom"_"end
            print chrom, snp_id, cum*100, end
        }}' {combined_map} > {plink_map}

        touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=executor)
