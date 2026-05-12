from gwf import Workflow
import subprocess
import pandas as pd
import pickle
from templates import *
import glob

gwf = Workflow()

##########################################
#    	     --- User inputs ---         #
##########################################

# Prerequisites: 
# 1 Finish running workflow.py in workflow_data/
# 2 Create ref/parameters_{group}.pkl mutation rate and generation time in the dictionary
#   {"mu":"{mutation rate}"
#   "generation":"{generation time}"}

genus_list      = ["Loxodonta", "Elephas", "Boselaphus", "Panthera", "Rhinoceros", "Ceratotherium", "Diceros"]
#genus_list =["Rhinoceros"] # -d 1300
#genus_list = ["Ceratotherium"] # -d 800
#genus_list = ["Diceros"] # - d 200
#genus_list = ["Boselaphus", "Panthera"] # -d 50 also uncia is here
#genus_list = ["Loxodonta"] # -d 75
#genus_list = ["Elephas","Panthera"] # -d 50 pardus, tigris and leo

# Filtering for missingness. 10 means keeping only < than 10% Missingness SNPs
filter = "10"

data            = pd.concat([pd.read_table(f) for f in [f"/faststorage/project/megaFauna/sa_megafauna/metadata/samples_{genus}.txt" for genus in genus_list]], ignore_index=True) # add all SRR accession from genus_list in a single data frame
data            = data.reset_index(drop=True)
ref_folders     = sorted(set(data.REFERENCE_FOLDER)) # list of references needed to map the SRR accessions

references      = pd.read_table(f"/faststorage/project/megaFauna/sa_megafauna/metadata/references.txt")
references      = references.loc[references.REFERENCE_FOLDER.isin(ref_folders)] # filter references to keep only the ones needed to map the SRR accessions
references      = references.reset_index(drop=True)

## make data frame that contains names of species-specific folders and the reference folders used to map the species
species_and_refs = pd.DataFrame({"FOLDER": data.FOLDER, "REFERENCE_FOLDER": data.REFERENCE_FOLDER, "GVCF_FOLDER": [data.GENUS.iloc[jj] + "_" + data.SPECIES.iloc[jj] for jj in range(data.shape[0])]}).drop_duplicates()
species_and_refs = species_and_refs.reset_index(drop=True)
## merge dataframes to have species-specific folder, reference folder and fastq ftps in same dataframe
species_and_refs = species_and_refs.merge(references, how = "left")

##########################################
#    	     --- Preparations ---        #
##########################################

for i in range(species_and_refs.shape[0]):
    # Initialising folders and variables for putting in the functions
    n_contigs_included = 0
    n_pops = 1

    group      = species_and_refs.FOLDER[i]
    ref_folder = species_and_refs.REFERENCE_FOLDER[i]

    # Load regions table
    df_regions = pd.read_csv(f"/faststorage/project/megaFauna/sa_megafauna/data/{ref_folder}/ref/regions_{ref_folder}_updated.txt", sep="\t")
    df_regions["end"] = pd.to_numeric(df_regions["end"], errors='coerce')     # Convert end column to integer (in case it was read as string)

    autosomal_regions = df_regions[(df_regions["female_ploidy"] == 2) & (df_regions["male_ploidy"] == 2) & (df_regions["end"] > 20000000)]
    chromosomes = autosomal_regions.sort_values("end", ascending=False)["region"].tolist()
    autosome_batch_list = autosomal_regions.sort_values("end", ascending=False)["batch"].drop_duplicates().tolist()

    # Limit number of contigs if needed
    if len(chromosomes) < n_contigs_included or n_contigs_included < 1:
        n_contigs_included = len(chromosomes)
    all_chromosomes = chromosomes
    subset_chromosomes = chromosomes[:n_contigs_included]

    #Load parameters such as mutation rate and generation time
    with open(f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/parameters_{group}.pkl", "rb") as f: 
        params = pickle.load(f)
    mu = params.get("mu")
    gen = params.get("generation")
    # Population list
    sample_dict = {}
    with open(f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/population_list.txt") as f:
        pops = [line.strip() for line in f]
    if len(pops) < n_pops or n_pops < 1:
        n_pops = len(pops)
    pops = pops[:n_pops]
    print(group, len(subset_chromosomes))

    ##########################################
    #   A  --- Population structure ---      #
    ##########################################
#    # A.1 Concatenate autosome VCFs
    job_id_concat_chrA = f"concatenate_autosomes_{group}"
    gwf.target_from_template(job_id_concat_chrA, 
                             vcfconcat(vcfs     = [f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/VCF/batches/{group}_batch_{i}_fploidy_2_mploidy_2_gt_snps.vcf.gz" for i in autosome_batch_list],
                                       contigs  = all_chromosomes,
                                       chrA_vcf = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/VCF/chrA_concat_{group}.vcf.gz",
                                       done     = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_concat_chrA))

#    # A.2 ADMIXTURE
#    job_id_ADMIXTURE = f"ADMIXTURE_{group}"
#    gwf.target_from_template(job_id_ADMIXTURE,
#                             ADMIXTURE(
#                                 group         = group,
#                                 vcf_path      = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/VCF/chrA_concat_{group}.vcf.gz",
#                                 admixture_out = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/ADMIXTURE/",
#                                 k_min         = 2,
#                                 k_max         = 6,
#                                 done_prev     = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_concat_chrA,
#                                 done          = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_ADMIXTURE))
#    
#    # A.3 PLOT ADMIXTURE
#    job_id_ADMIXTURE_PLOT = f"ADMIXTURE_PLOT_{group}"
#    gwf.target_from_template(job_id_ADMIXTURE_PLOT,
#                             ADMIXTURE_PLOT(
#                                 admixture_plot_script = f"/faststorage/project/megaFauna/people/laurids/scripts/ADMIXTURE/admixture_plot.py",
#                                 group                 = group,
#                                 admixture_out         = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/ADMIXTURE/",
#                                 plot_out              = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/ADMIXTURE/plots/",
#                                 k_min                 = 2,
#                                 k_max                 = 6,
#                                 done_prev             = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_ADMIXTURE,
#                                 done                  = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_ADMIXTURE_PLOT))
#    
#    # A.4 Parse "best" populations
#    # Inspect plots before continuing
#    job_id_ADMIXTURE_PARSE = f"ADMIXTURE_PARSE_{group}"
#    gwf.target_from_template(job_id_ADMIXTURE_PARSE,
#                             ADMIXTURE_PARSE(
#                                 admixture_parse_script = f"/faststorage/project/megaFauna/people/laurids/scripts/ADMIXTURE/admixture_parse.py",
#                                 group                  = group,
#                                 admixture_out          = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/ADMIXTURE",
#                                 k_min                  = 2,
#                                 k_max                  = 6,
#                                 done_prev              = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_ADMIXTURE_PLOT,
#                                 done                   = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_ADMIXTURE_PARSE))
#
#    # A.5 PCA
#    job_id_PCA = f"PCA_{group}"
#    gwf.target_from_template(job_id_PCA, 
#                             PCA(pca_script = f"/faststorage/project/megaFauna/people/laurids/scripts/PCA/PCA_script.py",
#                             group          = group,
#                             vcf_path       = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/VCF/batches/",
#                             pca_out        = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/PCA/",
#                             pruned_snps    = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/ADMIXTURE/{group}_pruned.prune.in",
#                             done_prev      = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_ADMIXTURE,
#                             done           = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_PCA))
#
#    # A.6 KING relatedness analysis for each population
#    for pop in pops:
#        with open(f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/{pop}.txt") as f:
#            sample_list = [line.strip() for line in f if line.strip()]
#        if len(sample_list) > 1:
#            job_id_relatedness = f"relatedness_{group}_{pop}"
#            gwf.target_from_template(
#                job_id_relatedness,
#                check_relatedness_king(
#                    group           = group,
#                    contigs         = all_chromosomes,
#                    samples         = sample_list,
#                    vcf_path        = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/VCF/chrA_concat_{group}.vcf.gz",
#                    relatedness_out = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/relatedness/{pop}",
#                    done_prev       = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_ADMIXTURE,
#                    done            = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_relatedness
#                )
#            )
#
    ##########################################
    #  B  --- mask and prep for smc++ ---    #
    ##########################################
    for pop in pops:
    # Read sample IDs from file for a given population SHOULD NEVER BE COMMENTED OUT
        with open(f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/{pop}_filtered.txt") as f:
            sample_list = [line.strip() for line in f if line.strip()]
        sample_dict[pop] = sample_list
        
        # B.1 Subset vcf file by population and filter missingness
        job_id_pop_and_missingness = f"population_and_missingness_filtering_{pop}_{group}"
        gwf.target_from_template(job_id_pop_and_missingness, 
                                 subset_and_filter(vcf_in       = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/VCF/chrA_concat_{group}.vcf.gz",
                                                   samples      = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/{pop}_filtered.txt",
                                                   filter       = filter,
                                                   pop_vcf      = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/VCF/chrA_{pop}_{group}.vcf.gz",
                                                   done_prev    = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_concat_chrA,
                                                   done         = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_pop_and_missingness))
        # B.2 - Mask beds per individual per chromosome
        sample_beds_dones = []
        sample_beds = []
        for sample in sample_dict[pop]:
            sample_chrom_beds_dones = []
            sample_chrom_beds = []
            for chrom in subset_chromosomes:
                job_id_sample_chrom_beds = f"coverage_mask_{group}_{sample}_{chrom}"
                sample_chrom_beds.append(f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/smcpp/masked_regions/{job_id_sample_chrom_beds}.bed")
                sample_chrom_beds_dones.append(f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_sample_chrom_beds)
                gwf.target_from_template(job_id_sample_chrom_beds,
                                         mask_beds(sample       = sample,
                                                   chrom        = chrom,
                                                   ref_folder   = ref_folder,
                                                   bed          = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/smcpp/masked_regions/{job_id_sample_chrom_beds}.bed",
                                                   done_prev    = [],
                                                   done         = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_sample_chrom_beds))
            # B.3 Merge to one bed per individual
            jobid_merge_per_sample = f"merge_mask_{group}_{sample}"
            sample_beds.append(f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/smcpp/masked_regions/{jobid_merge_per_sample}.bed")
            sample_beds_dones.append(f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + jobid_merge_per_sample)
            gwf.target_from_template(jobid_merge_per_sample,
                                     sample_merge_mask(beds = sample_chrom_beds,
                                            merged_bed      = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/smcpp/masked_regions/{jobid_merge_per_sample}.bed",
                                            done_prev       = sample_chrom_beds_dones,
                                            done            = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + jobid_merge_per_sample))
        # B.4 - Merge to one bed per population
        job_id_mask_merge = f"merge_mask_{pop}_{group}"
        gwf.target_from_template(job_id_mask_merge,
                                 merge_mask(beds = sample_beds,
                                            merged_bed  = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/smcpp/masked_regions/merge_mask_{pop}_{group}.bed",
                                            done_prev   = sample_beds_dones,
                                            done        = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_mask_merge))
    # B.5 - Mapability mask
    job_id_mappability_mask = f"mapability_mask_{group}"
    gwf.target_from_template(job_id_mappability_mask,
                             mappability_mask(fasta = glob.glob(f"/faststorage/project/megaFauna/sa_megafauna/data/{ref_folder}/ref/*LargerThan1000bp.fasta")[0],
                                        mask_bed    = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/masking/{ref_folder}_mappability_mask.bed",
                                        index       = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/masking/index/",
                                        genmap_out  = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/masking/",
                                        done        = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_mappability_mask))
    for pop in pops:
        # B.6 - Merged population bed with mappability bed for final bed 
        job_id_final_mask_merge = f"final_merge_mask_{pop}_{group}"
        gwf.target_from_template(job_id_final_mask_merge,
                                 final_merge_mask(map_bed   = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/masking/{ref_folder}_mappability_mask.bed",
                                                cov_bed     = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/smcpp/masked_regions/merge_mask_{pop}_{group}.bed",  
                                                final_bed   = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/smcpp/masked_regions/final_mask_{pop}_{group}.bed",
                                                done_prev   = [f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_mask_merge, f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_mappability_mask],
                                                done        = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_final_mask_merge))

        # B.7 - Compute mask statistics per individual and for merged beds
        job_id_mask_stats = f"mask_stats_{pop}_{group}"
        gwf.target_from_template(job_id_mask_stats,
                                 mask_stats(sample_beds         = sample_beds,
                                            cov_bed_merged  = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/smcpp/masked_regions/merge_mask_{pop}_{group}.bed",  
                                            mappability_bed = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/masking/{ref_folder}_mappability_mask.bed",
                                            final_bed       = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/smcpp/masked_regions/final_mask_{pop}_{group}.bed.gz",
                                            stats_file      = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/smcpp/masked_regions/mask_stats_{group}.txt",
                                            done_prev       = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_final_mask_merge,
                                            done            = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_mask_stats))

    ##########################################
    #      C     --- smc++ ---               #
    ##########################################
    for pop in pops:
        # C.1 - convert vcf file to smc files
        job_id_final_mask_merge = f"final_merge_mask_{pop}_{group}"
        smc_files = []
        vcf2smc_dones = []
        for chrom in chromosomes:
            job_id_vcf2smc = f"vcf2smc_{pop}_{chrom}_{group}"
            vcf2smc_dones.append(f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_vcf2smc)
            smc_files.append(f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/smcpp/smc_files/{pop}_{chrom}_{group}.smc.gz")
            gwf.target_from_template(job_id_vcf2smc,
                                     vcf2smc(vcf_in     = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/VCF/chrA_{pop}_{n_contigs_included}_{group}.vcf.gz",
                                             mask       = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/smcpp/masked_regions/final_mask_{pop}_{group}.bed.gz",
                                             chrom      = chrom,
                                             pop        = f"{pop}:{",".join(sample_dict[pop])}",
                                             smc_file   = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/smcpp/smc_files/{pop}_{chrom}_{group}.smc.gz",
                                             done_prev  = [f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_final_mask_merge, 
                                                           f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_pop_and_missingness],
                                             done       = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_vcf2smc))

        # C.2 - estimate population size
        job_id_smcpp_estimate = f"smcpp_estimate_{pop}_{n_contigs_included}_{group}"
        gwf.target_from_template(job_id_smcpp_estimate,
                                 smcpp_estimate(smc_files       = smc_files,
                                                mu              = mu,
                                                estimate_name   = job_id_smcpp_estimate,
                                                outdir          = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/smcpp/",
                                                done_prev       = vcf2smc_dones,
                                                done            = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_smcpp_estimate))

        # C.3 - create plot of population trajectory in years
        job_id_smcpp_plot = f"smcpp_plot_{pop}_{n_contigs_included}_{group}_{gen}"
        gwf.target_from_template(job_id_smcpp_plot,
                                 smcpp_plot(estimate_json   = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/smcpp/{job_id_smcpp_estimate}.final.json",
                                            generation      = gen,
                                            plot_name       = f"/faststorage/project/megaFauna/sa_megafauna/results/shared/smcpp/{job_id_smcpp_plot}.png",
                                            done_prev       = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_smcpp_estimate,
                                            done            = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_smcpp_plot))
                # C.4 - create plot of population trajectory in generations
        job_id_smcpp_plot = f"smcpp_plot_{pop}_{n_contigs_included}_{group}"
        gwf.target_from_template(job_id_smcpp_plot,
                                 smcpp_plot_generation(estimate_json   = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/smcpp/{job_id_smcpp_estimate}.final.json",
                                            plot_name       = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/smcpp/{job_id_smcpp_plot}.png",
                                            done_prev       = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_smcpp_estimate,
                                            done            = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_smcpp_plot))

    ##########################################
    #    	D      --- pyrho ---             #
    ##########################################
#    for pop in pops:
#        # D.1 - make pyrho lookup table
#        job_id_smcpp_plot = f"smcpp_plot_{pop}_{n_contigs_included}_{group}"
#        job_id_pyrho_table = f"pyrho_lookup_table_{pop}_{n_contigs_included}_{group}"
#        gwf.target_from_template(job_id_pyrho_table,
#                                 pyrho_lookup(estimate_csv  = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/smcpp/{job_id_smcpp_plot}.csv",
#                                              sample_size   = len(sample_dict[pop]) * 2,
#                                              mu            = mu,
#                                              pyrho_table   = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/pyrho/{job_id_pyrho_table}.hdf",
#                                              done_prev     = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_smcpp_plot,
#                                              done          = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_pyrho_table))
#
#        # D.2 - adjust pyrho hyperparameters
#        job_id_pyrho_hyperparam = f"pyrho_hyperparam_{pop}_{n_contigs_included}_{group}"
#        gwf.target_from_template(job_id_pyrho_hyperparam, 
#                                 pyrho_hyperparam(estimate_csv  = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/smcpp/{job_id_smcpp_plot}.csv",
#                                              sample_size       = len(sample_dict[pop]) * 2,
#                                              mu                = mu,
#                                              pyrho_table       = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/pyrho/{job_id_pyrho_table}.hdf",
#                                              hyperparam_file  = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/pyrho/{job_id_pyrho_hyperparam}.txt",
#                                              done_prev         = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_pyrho_table,
#                                              done              = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_pyrho_hyperparam))
#        
#        optimize_dones = []
#        rmaps = []
#        compute_dones = []
#        r2s = []
#        for chrom in subset_chromosomes:
#            # D.3 Create a VCF for each contig
#            job_id_contig_files = f"contig_files_{pop}_{chrom}_{group}"
#            gwf.target_from_template(job_id_contig_files,
#                                     make_contig_files(vcf_in   = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/VCF/chrA_{pop}_{group}.vcf.gz",
#                                                  chrom         = chrom,
#                                                  contig_vcf    = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/VCF/contigs/{pop}_{chrom}_{group}.vcf.gz",
#                                                  done_prev     = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_pop_and_missingness,
#                                                  done          = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_contig_files))
#            
#            # D.4 - pyrho recombination map estimation
#            job_id_pyrho_optimize = f"pyrho_optimize_{chrom}_{group}"
#            optimize_dones.append(f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_pyrho_optimize)
#            rmaps.append(f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/pyrho/{chrom}_{group}.rmap")
#            gwf.target_from_template(job_id_pyrho_optimize,
#                                     pyrho_optimize(contig_vcf = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/VCF/contigs/miss_{filter}/{pop}_{chrom}_{group}.vcf.gz",
#                                                    pyrho_table = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/pyrho/{job_id_pyrho_table}.hdf", 
#                                                    rmap_out= f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/pyrho/{chrom}_{group}.rmap",
#                                                    done_prev = [f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_pyrho_hyperparam, 
#                                                                 f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_contig_files], 
#                                                    done = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_pyrho_optimize))
#
#            # D.5 - compute r2
#            job_id_pyrho_r2 = f"pyrho_compute_r2_{chrom}_{group}"
#            compute_dones.append(f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_pyrho_r2)
#            r2s.append(f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/pyrho/r2_{chrom}_{group}.txt")
#            gwf.target_from_template(job_id_pyrho_r2, 
#                                     pyrho_compute(pyrho_table = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/pyrho/{job_id_pyrho_table}.hdf", 
#                                                   r2_out = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/pyrho/r2_{chrom}_{group}.txt", 
#                                                   done_prev = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_pyrho_optimize, 
#                                                   done = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_pyrho_r2))
#
#        # D.6 - combine recombination maps
#        job_id_combine_maps = f"combine_recomb_maps_{n_contigs_included}_{group}"
#        gwf.target_from_template(job_id_combine_maps,
#                                 combine_maps(rmaps = rmaps, 
#                                              combined_map = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/pyrho/chrA_recomb_map_{n_contigs_included}_{group}.txt", 
#                                              done_prev = optimize_dones, 
#                                              done = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_combine_maps))
#
#        # D.7 Make PLINK map
#        job_id_plink_map = f"plink_map_{n_contigs_included}_{group}"
#        gwf.target_from_template(job_id_plink_map,
#                                plink_map(combined_map = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/pyrho/chrA_recomb_map_{n_contigs_included}_{group}.txt", 
#                                          plink_map = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/pyrho/plink_map_{n_contigs_included}_{group}.map",
#                                          done_prev = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_combine_maps, 
#                                          done = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_plink_map))
#
#    ##########################################
#    #    	     --- GONE ---                #
#    ##########################################
#    for pop in pops:
#        # E.1 Make unzipped vcf for GONE 
#        job_id_unzip_vcf = f"unzip_vcf_{pop}_{n_contigs_included}_{group}"
#        gwf.target_from_template(job_id_pop_and_missingness, 
#                                 subset_and_filter(vcf_in       = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/VCF/chrA_{pop}_{group}.vcf.gz",
#                                                   chromosomes  = subset_chromosomes,
#                                                   subset_vcf   = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/VCF/chrA_{pop}_{n_contigs_included}_{group}.vcf",
#                                                   done_prev    = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_pop_and_missingness,
#                                                   done         = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_unzip_vcf))
#        
#        # E.2 - estimate population size    
#        job_id_gone = f"GONE_{pop}_{n_contigs_included}_{group}"
#        gwf.target_from_template(job_id_gone, 
#                                 GONE(chrA_pop      = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/VCF/chrA_{pop}{n_contigs_included}_{group}.vcf",
#                                      gone_estimate = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/GONE/{job_id_gone}",
#                                      done_prev     = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_unzip_vcf,
#                                      done          = f"/faststorage/project/megaFauna/sa_megafauna/results/{group}/done/" + job_id_gone))