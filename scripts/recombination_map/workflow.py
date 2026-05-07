from gwf import Workflow
import subprocess
import pandas as pd
import pickle
from templates import *
import glob

# Define workflow and set default options for all targets

gwf = Workflow()

##########################################
#    	     --- User inputs ---
##########################################

# Prerequisites: 
# 1 Finish running workflow.py
# 2 Run PCA analysis and identify populations in file: 
#   {results}/PCA/population_{population}_{species}.txt with a sample name in each line.
# 3 Create ref/parameters_{species}.pkl mutation rate and generation time in the dictionary
#   {"mu":"{mutation rate}"
#   "generation":"{generation time}"}

species_list = ["Panthera_tigris"]
population_name = "selected"

# Filtering for missingness. 10 means keeping only < than 10% Missingness SNPs
filter = "10"
#filter = ""

number_of_contigs_included = 0

##########################################
#    	     --- Preparations ---
##########################################

for species in species_list:
    # Initializing directories
    subprocess.run(["mkdir", "-p", f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/analysis_input/smcpp", 
                    f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/smcpp", 
                    f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/analysis_input/pyrho", 
                    f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/pyrho", 
                    f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/VCF/contigs/miss_{filter}", 
                    f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/done/", 
                    f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/analysis_input/smcpp/masked_regions"])

    # Load regions table
    df_regions = pd.read_csv(f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/ref/regions_{species}_updated.txt", sep="\t")
    # Extract only NC chromosomes with ploidy 2/2
    chromosomes = df_regions[
        df_regions["region"].str.startswith("NC") &
        (df_regions["female_ploidy"] == 2) &
        (df_regions["male_ploidy"] == 2)
    ]["region"].tolist()
    
    autosome_batch_list = df_regions[
        df_regions["region"].str.startswith("NC") &
        (df_regions["female_ploidy"] == 2) &
        (df_regions["male_ploidy"] == 2)
    ]["batch"].tolist()

    # Limit number of contigs if needed
    if len(chromosomes) < number_of_contigs_included or number_of_contigs_included <1:
        number_of_contigs_included = len(chromosomes)

    chromosomes = chromosomes[:number_of_contigs_included]
    contigs_included = f"_{len(chromosomes)}"

    # Read sample IDs from file for a given population
    with open(f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/PCA/population_{population_name}_{species}.txt") as f:
        sample_list = [line.strip() for line in f if line.strip()]
    # Build the final string

    with open(f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/parameters_{species}.pkl", "rb") as f: 
        params = pickle.load(f)

    mu = params.get("mu")
    mu_name = mu[:mu.find("e")]

    ##########################################
    #     --- smc++ file preparation ---
    ##########################################

    # A.1 Concatenate autosome VCFs
    job_id_concat_chrA = f"concatenate_autosomes_{species}"
    gwf.target_from_template(job_id_concat_chrA, 
                             vcfconcat(vcfs     = [f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/VCF/batches/{species}_batch_{i}_fploidy_2_mploidy_2_gt_snps.vcf.gz" for i in autosome_batch_list],
                                       chrA_vcf = f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/VCF/chrA_concat_{species}.vcf.gz",
                                       done     = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/done/" + job_id_concat_chrA))

#    # A.2 Choose a population for analysis and do missingness filtering
#    job_id_pop_and_missingness = f"population_and_missingness_filtering_{population_name}_{filter}_{species}"
#    gwf.target_from_template(job_id_pop_and_missingness, 
#                             subset_and_filter(vcf_in       = f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/VCF/chrA_concat_{species}.vcf.gz",
#                                               samples      = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/PCA/population_{population_name}_{species}.txt",
#                                               filter       = filter,
#                                               pop_vcf      = f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/VCF/chrA_{population_name}_{filter}_{species}.vcf.gz",
#                                               done_prev    = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/done/" + job_id_concat_chrA,
#                                               done         = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/done/" + job_id_pop_and_missingness))
#
#    # A.3 - Mask bed
#    sample_beds_dones = []
#    sample_beds = []
#    min_cov     = 5
#    min_het_ad  = 3
#    gq_min      = 30
#    for sample in sample_list:
#        sample_chrom_beds_dones = []
#        sample_chrom_beds = []
#        for chrom in chromosomes:
#            job_id_sample_chrom_beds = f"coverage_mask_{species}_{sample}_{chrom}"
#            sample_chrom_beds.append(f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/analysis_input/smcpp/masked_regions/{job_id_sample_chrom_beds}.bed")
#            sample_chrom_beds_dones.append(f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/done/" + job_id_sample_chrom_beds)
#            gwf.target_from_template(job_id_sample_chrom_beds,
#                                     mask_beds(sample       = sample,
#                                               chrom        = chrom,
#                                               bed          = f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/analysis_input/smcpp/masked_regions/{job_id_sample_chrom_beds}.bed",
#                                               parameters   = f"{min_cov},{min_het_ad},{gq_min}",
#                                               done_prev    = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/done/" + f"PCA_{species}",
#                                               done         = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/done/" + job_id_sample_chrom_beds))
#        
#        jobid_merge_per_sample = f"merge_mask_{species}_{sample}"
#        sample_beds.append(f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/analysis_input/smcpp/masked_regions/{jobid_merge_per_sample}.bed")
#        sample_beds_dones.append(f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/done/" + jobid_merge_per_sample)
#        gwf.target_from_template(jobid_merge_per_sample,
#                                 sample_merge_mask(beds    = sample_chrom_beds,
#                                        merged_bed  = f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/analysis_input/smcpp/masked_regions/{jobid_merge_per_sample}.bed",
#                                        done_prev   = sample_chrom_beds_dones,
#                                        done        = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/done/" + jobid_merge_per_sample))
#    # A.4 - Merged coverage mask
#    job_id_mask_merge = f"merge_mask_{population_name}_{species}"
#    gwf.target_from_template(job_id_mask_merge,
#                             merge_mask(beds = sample_beds,
#                                        merged_bed  = f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/analysis_input/smcpp/masked_regions/merge_mask_{population_name}_{species}.bed",
#                                        done_prev   = sample_beds_dones,
#                                        done        = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/done/" + job_id_mask_merge))
#    
#    # A.4 - Mapability mask
#    job_id_mappability_mask = f"mapability_mask_{species}"
#    gwf.target_from_template(job_id_mappability_mask,
#                             mappability_mask(fasta = glob.glob(f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/ref/*genomic_LargerThan1000bp.fasta")[0],
#                                        mask_bed    = f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/masking/{species}_mappability_mask.bed",
#                                        index       = f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/masking/index/",
#                                        genmap_out  = f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/masking/",
#                                        done        = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/done/" + job_id_mappability_mask))
#    
#    # A.4 - Merged coverage mask
#    job_id_final_mask_merge = f"final_merge_mask_{population_name}_{species}"
#    gwf.target_from_template(job_id_final_mask_merge,
#                             final_merge_mask(map_bed   = f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/masking/{species}_mappability_mask.bed",
#                                            cov_bed     = f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/analysis_input/smcpp/masked_regions/merge_mask_{population_name}_{species}.bed",  
#                                            final_bed   = f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/analysis_input/smcpp/masked_regions/final_mask_{population_name}_{species}.bed",
#                                            done_prev   = [f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/done/" + job_id_mask_merge, f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/done/" + job_id_mappability_mask],
#                                            done        = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/done/" + job_id_final_mask_merge))
#    
#    # A.4 - Merged coverage mask
#    job_id_mask_stats = f"mask_stats_{population_name}_{species}"
#    gwf.target_from_template(job_id_mask_stats,
#                             mask_stats(sample_beds         = sample_beds,
#                                        cov_bed_merged  = f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/analysis_input/smcpp/masked_regions/merge_mask_{population_name}_{species}.bed",  
#                                        mappability_bed = f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/masking/{species}_mappability_mask.bed",
#                                        final_bed       = f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/analysis_input/smcpp/masked_regions/final_mask_{population_name}_{species}.bed.gz",
#                                        stats_file      = f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/analysis_input/smcpp/masked_regions/mask_stats_{species}.txt",
#                                        done_prev       = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/done/" + job_id_final_mask_merge,
#                                        done            = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/done/" + job_id_mask_stats))
#
#    ##########################################
#    #    	     --- smc++ ---
#    ##########################################
#
##    # B.1 - convert vcf file to smc files
#    job_id_pop_and_missingness = f"population_and_missingness_filtering_{population_name}_{filter}_{species}"
#    job_id_final_mask_merge = f"final_merge_mask_{population_name}_{species}"
#    smc_files = []
#    vcf2smc_dones = []
#    for chrom in chromosomes:
#        job_id_vcf2smc = f"vcf2smc_{population_name}_{chrom}_{filter}_{species}"
#        vcf2smc_dones.append(f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/done/" + job_id_vcf2smc)
#        smc_files.append(f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/analysis_input/smcpp/{population_name}_{chrom}_{filter}_{species}.smc.gz")
#        gwf.target_from_template(job_id_vcf2smc,
#                                 vcf2smc(vcf_in     = f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/VCF/chrA_{population_name}_{filter}_{species}.vcf.gz",
#                                         mask       = f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/analysis_input/smcpp/masked_regions/final_mask_{population_name}_{species}.bed.gz",
#                                         chrom      = chrom,
#                                         pop        = f"{population_name}:{",".join(sample_list)}",
#                                         smc_file   = f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/analysis_input/smcpp/{population_name}_{chrom}_{filter}_{species}.smc.gz",
#                                         done_prev  = [f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/done/" + job_id_final_mask_merge, 
#                                                       f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/done/" + job_id_pop_and_missingness],
#                                         done       = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/done/" + job_id_vcf2smc))
#
#    # B.2 - estimate population size
#    job_id_smcpp_estimate = f"smcpp_estimate_{population_name}{contigs_included}_{mu_name}_{filter}_{species}"
#    gwf.target_from_template(job_id_smcpp_estimate,
#                             smcpp_estimate(smc_files       = smc_files,
#                                            mu              = mu,
#                                            estimate_name   = job_id_smcpp_estimate,
#                                            outdir          = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/smcpp/",
#                                            done_prev       = vcf2smc_dones,
#                                            done            = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/done/" + job_id_smcpp_estimate))
#
#    # B.3 - create plot of population trajectory
#    job_id_smcpp_plot = f"smcpp_plot_{population_name}{contigs_included}_{mu_name}_{filter}_{species}"
#    gwf.target_from_template(job_id_smcpp_plot,
#                             smcpp_plot(estimate_json   = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/smcpp/{job_id_smcpp_estimate}.final.json",
#                                        plot_name       = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/smcpp/{job_id_smcpp_plot}.png",
#                                        done_prev       = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/done/" + job_id_smcpp_estimate,
#                                        done            = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/done/" + job_id_smcpp_plot))
#
#    ##########################################
#    #    	  --- pyrho steps ---
#    ##########################################
#
#
#        # A.3 Create a VCF for each contig
#    job_id_pop_and_missingness = f"population_and_missingness_filtering_{population_name}_{filter}_{species}"
#    contig_dones = []
#    contig_files = []
#    for chrom in chromosomes:
#        job_id_contig_files = f"contig_files_{population_name}_{chrom}_{filter}_{species}"
#        contig_dones.append(f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/done/" + job_id_contig_files)
#        contig_files.append(f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/VCF/contigs/miss_{filter}/{population_name}_{chrom}_{filter}_{species}.vcf.gz")
#        gwf.target_from_template(job_id_contig_files,
#                                 make_contig_files(vcf_in   = f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/VCF/chrA_{population_name}_{filter}_{species}.vcf.gz",
#                                              chrom         = chrom,
#                                              contig_vcf    = f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/VCF/contigs/miss_{filter}/{population_name}_{chrom}_{filter}_{species}.vcf.gz",
#                                              done_prev     = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/done/" + job_id_pop_and_missingness,
#                                              done          = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/done/" + job_id_contig_files))
#
#    # C.1 - make pyrho lookup table
#    job_id_smcpp_plot = f"smcpp_plot_{population_name}{contigs_included}_{mu_name}_{filter}_{species}"
#    job_id_pyrho_table = f"pyrho_lookup_table_{population_name}{contigs_included}_{mu_name}_{filter}_{species}"
#    gwf.target_from_template(job_id_pyrho_table,
#                             pyrho_lookup(estimate_csv  = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/smcpp/{job_id_smcpp_plot}.csv",
#                                          sample_size   = len(sample_list) * 2,
#                                          mu            = mu,
#                                          pyrho_table   = f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/analysis_input/pyrho/{job_id_pyrho_table}.hdf",
#                                          done_prev     = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/done/" + job_id_smcpp_plot,
#                                          done          = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/done/" + job_id_pyrho_table))
#
#        # C.2 - adjust pyrho hyperparameters
#    job_id_pyrho_hyperparam = f"pyrho_hyperparam_{population_name}{contigs_included}_{mu_name}_{filter}_{species}"
#    gwf.target_from_template(job_id_pyrho_hyperparam, 
#                             pyrho_hyperparam(estimate_csv  = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/smcpp/{job_id_smcpp_plot}.csv",
#                                          sample_size       = len(sample_list) * 2,
#                                          mu                = mu,
#                                          pyrho_table       = f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/analysis_input/pyrho/{job_id_pyrho_table}.hdf",
#                                          hyperparam_file  = f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/analysis_input/pyrho/{job_id_pyrho_hyperparam}.txt",
#                                          done_prev         = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/done/" + job_id_pyrho_table,
#                                          done              = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/done/" + job_id_pyrho_hyperparam))
#
#    # C.3 - pyrho recombination map estimation
#    optimize_dones = []
#    rmaps = []
#    for chrom in chromosomes:
#        job_id_pyrho_optimize = f"pyrho_optimize_{chrom}_{mu_name}_{filter}_{species}"
#        optimize_dones.append(done_path + job_id_pyrho_optimize)
#        rmaps.append(f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/pyrho/{chrom}_{mu_name}_{filter}_{species}.rmap")
#        gwf.target_from_template(job_id_pyrho_optimize,
#                                 pyrho_optimize(contig_vcf = f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/VCF/contigs/miss_{filter}/{population_name}_{chrom}_{filter}_{species}.vcf.gz",
#                                                pyrho_table = f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/analysis_input/pyrho/{job_id_pyrho_table}.hdf", 
#                                                rmap_out= f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/pyrho/{chrom}_{mu_name}_{filter}_{species}.rmap",
#                                                done_prev = done_path + job_id_pyrho_hyperparam, 
#                                                done = done_path + job_id_pyrho_optimize))
#
#    # C.4 - compute r2 from the recombination map
#    compute_dones = []
#    r2s = []
#    for chrom in chromosomes:
#        job_id_pyrho_r2 = f"pyrho_compute_r2_{chrom}_{mu_name}_{filter}_{species}"
#        compute_dones.append(done_path + job_id_pyrho_r2)
#        r2s.append(f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/pyrho/r2_{chrom}_{mu_name}_{filter}_{species}.txt")
#        gwf.target_from_template(job_id_pyrho_r2, 
#                                 pyrho_compute(pyrho_table = f"/faststorage/project/megaFauna/sa_megafauna/data/{species}/analysis_input/pyrho/{job_id_pyrho_table}.hdf", 
#                                               r2_out = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/pyrho/r2_{chrom}_{mu_name}_{filter}_{species}.txt", 
#                                               done_prev = done_path + job_id_pyrho_optimize, 
#                                               done = done_path + job_id_pyrho_r2))
#
#    # C.5 - combine recombination maps
#    job_id_combine_maps = f"combine_recomb_maps{contigs_included}_{mu_name}_{filter}_{species}"
#    gwf.target_from_template(job_id_combine_maps,
#                             combine_maps(rmaps = rmaps, 
#                                          combined_map = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/pyrho/chrA_recomb_map{contigs_included}_{mu_name}_{filter}_{species}.txt", 
#                                          done_prev = done_path + job_id_pyrho_r2, 
#                                          done = done_path + job_id_combine_maps))
#
#    # Make PLINK map
#    job_id_plink_map = f"plink_map{contigs_included}_{mu_name}_{filter}_{species}"
#    gwf.target_from_template(job_id_plink_map,
#                            plink_map(combined_map = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/pyrho/chrA_recomb_map{contigs_included}_{mu_name}_{filter}_{species}.txt", 
#                                      plink_map = f"/faststorage/project/megaFauna/sa_megafauna/results/{species}/pyrho/plink_map{contigs_included}_{mu_name}_{filter}_{species}.map",
#                                      done_prev = done_path + job_id_combine_maps, 
#                                      done = done_path + job_id_plink_map))