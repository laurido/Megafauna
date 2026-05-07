"""
------------------------------------------------------------------------------------------------------------------------
This workflow handles raw files (reference genomes, paired-end fastqs), maps reads, calls and genotypes variants.
------------------------------------------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------------------------------------
Authors: Juraj Bergman, Vasili Pankratov, Bjarke M. Pedersen, Laurids Svenningsen, Antonia Mehlig
Date: 26/03/2025
------------------------------------------------------------------------------------------------------------------------
"""

import os
import numpy as np
import subprocess
from gwf import Workflow, AnonymousTarget
from templates import *
import pandas as pd
import re
gwf = Workflow()

genus_list      = ["Loxodonta", "Elephas", "Boselaphus", "Panthera", "Naja", "Rhinoceros", "Ceratotherium", "Diceros"]
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

"""
------------------------------------------------------------------------------------------------------------------------
A. REFERENCE-ASSOCIATED JOBS
------------------------------------------------------------------------------------------------------------------------
"""
# run jobs that download reference, make a reference that contains only contigs more or equal to 1000 bp in length and makes a referece-specific regions file (to be manually modified as necessary due to X, PAR, Y, mitochondrial ploidy, etc.
# submit jobs
#for i in range(references.shape[0]):
#    ## A.0. Initialize directories
#    subprocess.run(["mkdir", "-p", "/faststorage/project/megaFauna/sa_megafauna/data/{}/done".format(references.REFERENCE_FOLDER[i])])
#    subprocess.run(["mkdir", "-p", "/faststorage/project/megaFauna/sa_megafauna/data/{}/ref".format(references.REFERENCE_FOLDER[i])])
#    subprocess.run(["mkdir", "-p", "/faststorage/project/megaFauna/sa_megafauna/data/{}/fastq".format(references.REFERENCE_FOLDER[i])])
#    subprocess.run(["mkdir", "-p", "/faststorage/project/megaFauna/sa_megafauna/data/{}/bam".format(references.REFERENCE_FOLDER[i])])
#    subprocess.run(["mkdir", "-p", "/faststorage/project/megaFauna/sa_megafauna/data/{}/gVCF".format(references.REFERENCE_FOLDER[i])])
#    ## A.1. download reference
#    jobid_download_ref = "download_ref_" + references.ref_genome_name[i].replace("-","_")
#    gwf.target_from_template(jobid_download_ref,
#                             download_ref(ftp  = references.ftp[i],
#                                          out  = "/faststorage/project/megaFauna/sa_megafauna/data/{}/ref/{}".format(references.REFERENCE_FOLDER[i], str(references.ftp[i]).split("/")[-1][:-3]),
#                                          done = "/faststorage/project/megaFauna/sa_megafauna/data/{}/done/{}".format(references.REFERENCE_FOLDER[i],jobid_download_ref)))
#    # A.2. mask reference - uncomment if reference needs to be masked; "masked_regions.bed is expected in the reference folder
#    #jobid_mask_reference = "mask_reference_" + references.ref_genome_name[i].replace("-", "_")
#    #gwf.target_from_template(jobid_mask_reference,
#    #                         mask_reference(input="/faststorage/project/megaFauna/sa_megafauna/data/{}/ref/{}".format(references.REFERENCE_FOLDER[i], str(references.ftp[i]).split("/")[-1][:-3]),
#    #                                        bed="/faststorage/project/megaFauna/sa_megafauna/data/{}/ref/{}".format(references.REFERENCE_FOLDER[i], "masked_regions.bed"),
#    #                                        output="/faststorage/project/megaFauna/sa_megafauna/data/{}/ref/{}".format(references.REFERENCE_FOLDER[i], str(references.ftp[i]).split("/")[-1][:-7] + "_masked.fna"),
#    #                                        done_prev="/faststorage/project/megaFauna/sa_megafauna/data/{}/done/{}".format(references.REFERENCE_FOLDER[i], jobid_download_ref),
#    #                                        done="/faststorage/project/megaFauna/sa_megafauna/data/{}/done/{}".format(references.REFERENCE_FOLDER[i], jobid_mask_reference)))
#    # A.3. cut out contigs less than 1000 bp in length
#    jobid_cut_contigs = "cut_contigs_" + references.ref_genome_name[i].replace("-", "_")
#    gwf.target_from_template(jobid_cut_contigs,
#                             cut_contigs(input     = "/faststorage/project/megaFauna/sa_megafauna/data/{}/ref/{}".format(references.REFERENCE_FOLDER[i], str(references.ftp[i]).split("/")[-1][:-3]),
#                             # cut_contigs(input      = "/faststorage/project/megaFauna/sa_megafauna/data/{}/ref/{}".format(references.REFERENCE_FOLDER[i], str(references.ftp[i]).split("/")[-1][:-7]+"_masked.fna"), # if masked reference present comment out previous line and uncomment this line
#                                         output     = "/faststorage/project/megaFauna/sa_megafauna/data/{}/ref/{}".format(references.REFERENCE_FOLDER[i], str(references.ftp[i]).split("/")[-1][:-7]+"_LargerThan1000bp.fasta"),
#                                         min_length = 1000,
#                                         done_prev  = "/faststorage/project/megaFauna/sa_megafauna/data/{}/done/{}".format(references.REFERENCE_FOLDER[i], jobid_download_ref),
#                                         # done_prev="/faststorage/project/megaFauna/sa_megafauna/data/{}/done/{}".format(references.REFERENCE_FOLDER[i], jobid_mask_reference), # if masked reference present comment out previous line and uncomment this line
#                                         done       = "/faststorage/project/megaFauna/sa_megafauna/data/{}/done/{}".format(references.REFERENCE_FOLDER[i], jobid_cut_contigs)))
#    ## A.4. make reference-associated files
#    jobid_make_fasta = "make_fasta_" + references.ref_genome_name[i].replace("-","_")
#    gwf.target_from_template(jobid_make_fasta,
#                             make_fasta(input     = "/faststorage/project/megaFauna/sa_megafauna/data/{}/ref/{}".format(references.REFERENCE_FOLDER[i], str(references.ftp[i]).split("/")[-1][:-3]),
#                             # make_fasta(input="/faststorage/project/megaFauna/sa_megafauna/data/{}/ref/{}".format(references.REFERENCE_FOLDER[i], str(references.ftp[i]).split("/")[-1][:-7]+"_masked.fna"), ## if masked reference present comment out previous line and uncomment this line
#                                        done_prev = "/faststorage/project/megaFauna/sa_megafauna/data/{}/done/{}".format(references.REFERENCE_FOLDER[i], jobid_cut_contigs),
#                                        done      = "/faststorage/project/megaFauna/sa_megafauna/data/{}/done/{}".format(references.REFERENCE_FOLDER[i], jobid_make_fasta)))
#    ## A.5. make regions files
#    jobid_make_regions = "make_regions_" + references.ref_genome_name[i].replace("-","_")
#    gwf.target_from_template(jobid_make_regions,
#                            make_regions(refFolder = references.REFERENCE_FOLDER[i],
#                                         input     = str(references.ftp[i]).split("/")[-1][:-3],
#                                         # input=str(references.ftp[i]).split("/")[-1][:-7]+"_masked.fna", ## if masked reference present comment out previous line and uncomment this line
#                                         done_prev = "/faststorage/project/megaFauna/sa_megafauna/data/{}/done/{}".format(references.REFERENCE_FOLDER[i], jobid_make_fasta),
#                                         done      = "/faststorage/project/megaFauna/sa_megafauna/data/{}/done/{}".format(references.REFERENCE_FOLDER[i], jobid_make_regions)))
#"""
#------------------------------------------------------------------------------------------------------------------------
#B. SAMPLE DOWNLOAD AND PREPARATION FOR MAPPING
#------------------------------------------------------------------------------------------------------------------------
#"""
# submit jobs
#for i in range(species_and_refs.shape[0]):
#    group = species_and_refs.FOLDER[i]
#    inds  = list(data.loc[data.FOLDER == group].IND_ID.drop_duplicates())
#    srrs  = list(data.loc[data.FOLDER == group].fastq_ftp)
#  
#    ## B.0. Initialize directories
#    subprocess.run(["mkdir", "-p", f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done"])
#    subprocess.run(["mkdir", "-p", f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/fastq"])
#    subprocess.run(["mkdir", "-p", f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/bam"])
#    subprocess.run(["mkdir", "-p", f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/gVCF"])
#    for ind in inds:   
#        ## B.1. Download fastq files per individual
#        data_subset = data.loc[data.IND_ID == ind]
#        jobid_download_per_individual = f"download_per_individual_{group}_{ind}"
#        gwf.target_from_template(jobid_download_per_individual,
#                                download_per_individual(group          = group,
#                                                        ind            = ind,
#                                                        run_accessions = ",".join(list(data_subset.fastq_ftp)),
#                                                        md5s           = ",".join(list(data_subset.fastq_md5)),
#                                                        done           = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_download_per_individual}"))
#        ## B.2. Concatenate fastq files per individual if more than one accession are available, else rename the fastq files; prepare uBAMs for mapping
#        if data_subset.shape[0] > 2:
#            jobid_concat_or_rename_fastqs = f"concat_or_rename_fastqs_{group}_{ind}"
#            gwf.target_from_template(jobid_concat_or_rename_fastqs,
#                                    concatfastqs(group     = group,
#                                                ind       = ind,
#                                                fastqs_1  = ["/faststorage/project/megaFauna/sa_megafauna/data/" + group + "/fastq/" + data_subset.loc[data_subset.R1_or_R2 == "R1"].fastq_ftp.iloc[jj].split("/")[-1] for jj in range(data_subset.loc[data_subset.R1_or_R2 == "R1"].shape[0])],
#                                                fastqs_2  = ["/faststorage/project/megaFauna/sa_megafauna/data/" + group + "/fastq/" + data_subset.loc[data_subset.R1_or_R2 == "R2"].fastq_ftp.iloc[jj].split("/")[-1] for jj in range(data_subset.loc[data_subset.R1_or_R2 == "R2"].shape[0])],
#                                                prev_done = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_download_per_individual}",# +["/faststorage/project/megaFauna/sa_megafauna/data/" + group + "/done/download_pe2_" + group + "_" + srr.replace("/", "_") for srr in data_subset.fastq_ftp],
#                                                done      = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_concat_or_rename_fastqs}"))
#        else:
#            jobid_concat_or_rename_fastqs = f"concat_or_rename_fastqs_{group}_{ind}"
#            gwf.target_from_template(jobid_concat_or_rename_fastqs,
#                                    renamefastqs(group    = group,
#                                                ind       = ind,
#                                                fastq_1   = "/faststorage/project/megaFauna/sa_megafauna/data/" + group + "/fastq/" + data_subset.loc[data_subset.R1_or_R2 == "R1"].fastq_ftp.iloc[0].split("/")[-1],
#                                                fastq_2   = "/faststorage/project/megaFauna/sa_megafauna/data/" + group + "/fastq/" + data_subset.loc[data_subset.R1_or_R2 == "R2"].fastq_ftp.iloc[0].split("/")[-1],
#                                                prev_done = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_download_per_individual}",# +["/faststorage/project/megaFauna/sa_megafauna/data/" + group + "/done/download_pe2_" + group + "_" + srr.replace("/", "_") for srr in data_subset.fastq_ftp],
#                                                done      = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_concat_or_rename_fastqs}"))
#        ## B.3. Make uBAMs from concatenated or renamed fastqs
#        jobid_makeuBAM = f"makeuBAM_{group}_{ind}"
#        gwf.target_from_template(jobid_makeuBAM,
#                                makeuBAM(group=group,
#                                        ind       = ind,
#                                        fastq1    = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/fastq/{ind}_R1.fastq.gz",
#                                        fastq2    = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/fastq/{ind}_R2.fastq.gz",
#                                        prev_done = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_concat_or_rename_fastqs}",
#                                        done      = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_makeuBAM}"))
#        ## B.4. Split uBAMs
#        jobid_splituBAM = f"splituBAM_{group}_{ind}"
#        gwf.target_from_template(jobid_splituBAM,
#                                splituBAM(group     = group,
#                                        ind       = ind,
#                                        prev_done = [f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_makeuBAM}"],
#                                        done      = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_splituBAM}"))
#      
# """
# ------------------------------------------------------------------------------------------------------------------------
# C. MAPPING
# ------------------------------------------------------------------------------------------------------------------------
# """
 ## prequisite to start mapping is the existence of the f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_splituBAM}" done file

## submit jobs
#for i in range(species_and_refs.shape[0]):
#    group      = species_and_refs.FOLDER[i]
#    ref_folder = species_and_refs.REFERENCE_FOLDER[i]
#    ref_path   = f"/faststorage/project/megaFauna/sa_megafauna/data/{ref_folder}/ref/" + [file for file in os.listdir(f"/faststorage/project/megaFauna/sa_megafauna/data/{ref_folder}/ref/") if file.endswith("_LargerThan1000bp.fasta")][0]
#    inds       = list(data.loc[data.FOLDER == group].IND_ID.drop_duplicates())
#
#    for ind in inds:
#        if os.path.isfile(f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/bam/{ind}_nsplitubams.txt"):
#            jobid_splituBAM = f"splituBAM_{group}_{ind}" ## prequisite to start mapping is the existence of the f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_splituBAM}" done file
#            if os.path.isfile(f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_splituBAM}") and os.path.isfile(f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_splituBAM}"): ## checks if f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_splituBAM}" done file exists
#                with open(f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/bam/{ind}_nsplitubams.txt") as f:
#                    for l in f:
#                        nbams = int(l.strip())
#                with open(f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/bam/{ind}_bpreads.txt") as f:
#                    for l in f:
#                        nreads = float(l.strip())
#                        nbams = min(int(np.ceil(30/(48000000*nreads/4000000000))), nbams)
#                
#                for ishard in range(nbams):
#                    shard = shardstr(ishard)
#
#                
#                    # ## C.1. Mark the adapters in each shard uBAM file
#                    jobid_markadapt = f"markadapt_{group}_{ind}_{shard}"
#                    gwf.target_from_template(jobid_markadapt,
#                                                markadapt(group   = group,
#                                                        ind       = ind,
#                                                        shard     = shard,
#                                                        prev_done = [f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_splituBAM}"],
#                                                        done      = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_markadapt}"))
#
#                    ## C.2. Map each shard uBAM file
#                    jobid_mapBAM = f"mapBAM_{group}_{ind}_{shard}"
#                    gwf.target_from_template(jobid_mapBAM,
#                                                mapBAM(group  = group,
#                                                    ind       = ind,
#                                                    shard     = shard,
#                                                    ref       = ref_path,
#                                                    prev_done = [f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_markadapt}"],
#                                                    done      = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_mapBAM}"))
#                    
#                ## C.3 Merge mapped shard bams into single bam per individual
#                jobid_mergeBAMs = f"mergeBAMs_{group}_{ind}"
#                gwf.target_from_template(jobid_mergeBAMs,
#                                            mergeBAMs(group     = group,
#                                                    ind       = ind,
#                                                    bams      = [f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/bam/split_uBAM{ind}/shard_{shardstr(ishard)}_markadapt_mapped.bam" for ishard in range(nbams)],
#                                                    prev_done = [f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/mapBAM_{group}_{ind}_{shardstr(ishard)}" for ishard in range(nbams)],
#                                                    done      = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_mergeBAMs}"))
#                    
#                ## C.4. Mark and remove duplicates from bam
#                jobid_markduplicates = f"markduplicates_{group}_{ind}"
#                gwf.target_from_template(jobid_markduplicates,
#                                            markduplicates(group     = group,
#                                                        ind       = ind,
#                                                        prev_done = [f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_mergeBAMs}"],
#                                                        done      = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_markduplicates}"))
#
#                ## C.5. Sort bam by coordinates
#                jobid_coordsort = f"coordsort_{group}_{ind}"
#                gwf.target_from_template(jobid_coordsort,
#                                            coordsort(group     = group,
#                                                    ind       = ind,
#                                                    prev_done = [f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_markduplicates}"],
#                                                    done      = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_coordsort}"))
#
#                ## C.6. Get coverage statistics
#                regions      = pd.read_table(f"/faststorage/project/megaFauna/sa_megafauna/data/{ref_folder}/ref/regions_{ref_folder}.txt")
#                regions_list = list(regions.region)
#                chrom_list   = list(regions.chrom)
#                start_list   = [jj + 1 for jj in list(regions.start)]
#                end_list     = list(regions.end)
#
#                jobid_cov = f"cov_{group}_{ind}"
#                gwf.target_from_template(jobid_cov,
#                                            cov(group       = group,
#                                                ind         = ind,
#                                                regions     = regions_list,
#                                                chromosomes = chrom_list,
#                                                starts      = start_list,
#                                                ends        = end_list,
#                                                prev_done   = [f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_coordsort}"],
#                                                done        = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_cov}"))
#                    
#                ## C.7 Use this instead of C.6. if you need to get coverage statistics when there are too many contigs in reference
#                # regions      = pd.read_table(f"/faststorage/project/megaFauna/sa_megafauna/data/{ref_folder}/ref/regions_{ref_folder}.txt")
#                # regions_list = list(regions.region)
#                # chrom_list   = list(regions.chrom)
#                # start_list   = [jj + 1 for jj in list(regions.start)]
#                # end_list     = list(regions.end)
#    
#                # no_regions_per_batch = 10000
#                   
#                # cov_done_files = []
#                # for b in range(len(regions_list)//no_regions_per_batch):
#                #     jobid_cov_batched = f"cov_{group}_{ind}_batch_{b}"
#                #     cov_done_files.append(f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_cov_batched}")
#                #     gwf.target_from_template(jobid_cov_batched,
#                #                              cov_batched(group       = group,
#                #                                          ind         = ind,
#                #                                          batch       = b,
#                #                                          regions     = regions_list[(b*no_regions_per_batch):(b*no_regions_per_batch+no_regions_per_batch)],
#                #                                          chromosomes = chrom_list[(b*no_regions_per_batch):(b*no_regions_per_batch+no_regions_per_batch)],
#                #                                          starts      = start_list[(b*no_regions_per_batch):(b*no_regions_per_batch+no_regions_per_batch)],
#                #                                          ends        = end_list[(b*no_regions_per_batch):(b*no_regions_per_batch+no_regions_per_batch)],
#                #                                          prev_done   = [f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_coordsort}"],
#                #                                          done        = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_cov_batched}"))
#                # b = b + 1
#                # jobid_cov_batched = f"cov_{group}_{ind}_batch_{b}"
#                # cov_done_files.append(f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_cov_batched}")
#                # gwf.target_from_template(jobid_cov_batched,
#                #                          cov_batched(group       = group,
#                #                                      ind         = ind,
#                #                                      batch       = b,
#                #                                      regions     = regions_list[(b*no_regions_per_batch):],
#                #                                      chromosomes = chrom_list[(b*no_regions_per_batch):],
#                #                                      starts      = start_list[(b*no_regions_per_batch):],
#                #                                      ends        = end_list[(b*no_regions_per_batch):],
#                #                                      prev_done   = [f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_coordsort}"],
#                #                                      done        = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_cov_batched}"))
#                
#                # jobid_concatenate_cov_files= f"concatenate_cov_files_{group}_{ind}"
#                # gwf.target_from_template(jobid_concatenate_cov_files,
#                #                          concatenate_cov_files(files = [f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/cov/{ind}_batch_{b}.cov" for b in range(len(regions_list)//no_regions_per_batch + 1)],
#                #                                                out_file = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/cov/{ind}.cov", 
#                #                                                prev_done = cov_done_files, 
#                #                                                done = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/cov_{group}_{ind}"))
                        
# """
# ------------------------------------------------------------------------------------------------------------------------
# D. CALLING
# ------------------------------------------------------------------------------------------------------------------------
# """

###------------------------------------------------ pre-calling jobs ------------------------------------------------###
# D.1. find X contigs and update regions file
#for ref_folder in ref_folders:
#    sub_data = data[data.REFERENCE_FOLDER == ref_folder]
#    sub_data.to_csv(f"/faststorage/project/megaFauna/sa_megafauna/metadata/samples_{ref_folder}.txt", sep="\t", index=False)
#    jobid_find_chrX = f'find_chrX_{ref_folder}'
#    gwf.target_from_template(jobid_find_chrX,
#                                find_chrX(subset_file   = f"/faststorage/project/megaFauna/sa_megafauna/metadata/samples_{ref_folder}.txt",
#                                            ref_folder  = ref_folder,
#                                            contigs     = "/faststorage/project/megaFauna/sa_megafauna/main_workflow/additional_scripts/special_contigs.txt",
#                                            minlen      = 1e6,
#                                            prev_done   = [],
#                                            done        = f"/faststorage/project/megaFauna/sa_megafauna/data/{ref_folder}/done/{jobid_find_chrX}"))
#
#    ## D.2. make simplified batch file from regions file
#    jobid_make_simplified_batch_file = f'make_simplified_batch_file_{ref_folder}'
#    gwf.target_from_template(jobid_make_simplified_batch_file,
#                                make_simplified_batch_file(regions_file = f"/faststorage/project/megaFauna/sa_megafauna/data/{ref_folder}/ref/regions_{ref_folder}_updated.txt",
#                                                            ref_folder   = ref_folder,
#                                                            prev_done    = [f"/faststorage/project/megaFauna/sa_megafauna/data/{ref_folder}/done/{jobid_find_chrX}"],
#                                                            done         = f"/faststorage/project/megaFauna/sa_megafauna/data/{ref_folder}/done/{jobid_make_simplified_batch_file}"))

###-------  calling at individual level with sp_ssp designation (maintained so that done files are recognized) ------###
## prequisite to start calling and genotyping is the existence of the f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_cov}" and f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_make_loc_metadata}" done files
    
#for i in range(species_and_refs.shape[0]):
#    group           = species_and_refs.FOLDER[i]
#    ref_folder      = species_and_refs.REFERENCE_FOLDER[i]
#    inds = []
#    ref_path        = f"/faststorage/project/megaFauna/sa_megafauna/data/{ref_folder}/ref/" + [file for file in os.listdir(f"/faststorage/project/megaFauna/sa_megafauna/data/{ref_folder}/ref/") if file.endswith("_LargerThan1000bp.fasta")][0]
#    inds_to_include = list(data.loc[data.FOLDER == group]["IND_ID"].drop_duplicates())
#    inds_and_sexes  = pd.read_table(f"/faststorage/project/megaFauna/sa_megafauna/data/{ref_folder}/ref/samples_coverage_stats.txt")
#    inds_and_sexes  = inds_and_sexes[inds_and_sexes.IND_ID.isin(inds_to_include)]
#    inds            = list(inds_and_sexes.IND_ID)
#    sexes           = list(inds_and_sexes.gSEX)
#
#    regions = pd.read_table(f"/faststorage/project/megaFauna/sa_megafauna/data/{ref_folder}/ref/batches_{ref_folder}.txt")
#    batches = list(regions.batch)
#
#    # ## D.3. make metadata files
#    length_of_chunk = 2000000
#    jobid_make_simplified_batch_file = f'make_simplified_batch_file_{ref_folder}'
#    jobid_make_loc_metadata = f'make_loc_metadata_{group}'
#    gwf.target_from_template(jobid_make_loc_metadata,
#                                make_batch_metadata(loc          = length_of_chunk,
#                                                    group        = group,
#                                                    regions_file = f"/faststorage/project/megaFauna/sa_megafauna/data/{ref_folder}/ref/regions_{ref_folder}_updated.txt",
#                                                    inds         = ",".join(inds),
#                                                    sexes        = ",".join(sexes),
#                                                    ref_folder   = ref_folder,
#                                                    prev_done    = [f"/faststorage/project/megaFauna/sa_megafauna/data/{ref_folder}/done/{jobid_make_simplified_batch_file}"],
#                                                    done         = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_make_loc_metadata}"))
#
#    ## D.4. Call variants per individual considering region ploidy (as determined by the sex of the individual)
#    for j in range(len(inds)):
#        jobid_cov = f"cov_{group}_{inds[j]}" ## prequisite to start calling and genotyping is the existence of the f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_cov}" done file
#        if os.path.isfile(f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_cov}"): ## checks if  f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_cov}" done file exists
#            ## get ploidy of regions based on sex of individual
#            if sexes[j] == "M":
#                ploidies_list = list(regions.male_ploidy)
#            else:
#                ploidies_list = list(regions.female_ploidy)
#
#            for batch in sorted(list(set(regions.batch))):
#                iis = [index for index in range(len(batches)) if batches[index] == batch]
#                set_of_ploidies = sorted(list(set([ploidies_list[ii] for ii in iis if ploidies_list[ii] != 0])))
#                for pl in set_of_ploidies:
#                    jobid_call = f"call_with_bed_{group}_{inds[j]}_batch_{batch}_ploidy_{pl}"
#                    gwf.target_from_template(jobid_call,
#                                                call_batch_with_bed(group = group,
#                                                    ind                = inds[j],
#                                                    batch              = batch,
#                                                    bed                = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/gVCF/beds_and_intervals/{inds[j]}_batch_{batch}_ploidy_{pl}.bed",
#                                                    intervals          = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/gVCF/beds_and_intervals/{inds[j]}_batch_{batch}_ploidy_{pl}.intervals",
#                                                    ref                = ref_path,
#                                                    ploidy             = pl,
#                                                    prev_done          = [f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_cov}", f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_make_loc_metadata}"],
#                                                    done               = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_call}"))
#
# """
# ------------------------------------------------------------------------------------------------------------------------
# E. GENOTYPING
# ------------------------------------------------------------------------------------------------------------------------
# """

###------------------------------------------------ pre-genotyping jobs ------------------------------------------------###

#for gvcf_folder in species_and_refs.GVCF_FOLDER.drop_duplicates():
#    
#    subprocess.run(["mkdir", "-p", "/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done"])
#    subprocess.run(["mkdir", "-p", "/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/gVCF"])
#
#    length_of_chunk = 2000000
#    groups      = list(species_and_refs.loc[species_and_refs.GVCF_FOLDER == gvcf_folder].FOLDER)
#    rr_folders  = list(species_and_refs.loc[species_and_refs.GVCF_FOLDER == gvcf_folder].REFERENCE_FOLDER)
#
#    inds       = []
#    sexes      = []
#    folders    = []
#    
#    ### gather inidividuals, sexes and done files
#    for i in range(len(groups)):
#        group           = groups[i]
#        ref_folder      = rr_folders[i]
#        ref_path        = f"/faststorage/project/megaFauna/sa_megafauna/data/{ref_folder}/ref/" + [file for file in os.listdir(f"/faststorage/project/megaFauna/sa_megafauna/data/{ref_folder}/ref/") if file.endswith("_LargerThan1000bp.fasta")][0]
#        inds_to_include = list(data.loc[data.FOLDER == group]["IND_ID"].drop_duplicates())
#        inds_and_sexes  = pd.read_table(f"/faststorage/project/megaFauna/sa_megafauna/data/{ref_folder}/ref/samples_coverage_stats.txt")
#        inds_and_sexes  = inds_and_sexes[inds_and_sexes.IND_ID.isin(inds_to_include)]
#        group_inds      = list(inds_and_sexes.IND_ID)
#        group_sexes     = list(inds_and_sexes.gSEX)
#
#        regions = pd.read_table(f"/faststorage/project/megaFauna/sa_megafauna/data/{ref_folder}/ref/batches_{ref_folder}.txt")
#        batches = list(regions.batch)
#
#        for j in range(len(group_inds)):
#            inds.append(group_inds[j])
#            sexes.append(group_sexes[j])
#            folders.append(group)
#
#    ## E.1. Make metadata files for genotyping
#    length_of_chunk = 2000000
#    jobid_make_simplified_batch_file = f'make_simplified_batch_file_{ref_folder}'
#    jobid_make_geno_metadata = f'make_geno_metadata_{gvcf_folder}'
#    gwf.target_from_template(jobid_make_geno_metadata,
#                                make_geno_metadata(loc           = length_of_chunk,
#                                                    group       = gvcf_folder, # ",".join(folders),
#                                                    regions_file = f"/faststorage/project/megaFauna/sa_megafauna/data/{ref_folder}/ref/regions_{ref_folder}_updated.txt",
#                                                    inds         = ",".join(inds),
#                                                    sexes        = ",".join(sexes),
#                                                    ref_folder   = ref_folder,
#                                                    prev_done    = [f"/faststorage/project/megaFauna/sa_megafauna/data/{ref_folder}/done/{jobid_make_simplified_batch_file}"],
#                                                    done         = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/{jobid_make_geno_metadata}"))
#    
#    ## E.2. Make a GenomicsDB folder and cohort.sample_map file for each batch and ploidy category within batch
#    for batch in sorted(list(set(regions.batch))):
#        iis                 = [index for index in range(len(batches)) if batches[index] == batch]
#        set_of_ploidy_pairs = pd.DataFrame({"female_ploidies":[list(regions.female_ploidy)[ii] for ii in iis], "male_ploidies":[list(regions.male_ploidy)[ii] for ii in iis]}).drop_duplicates()
#        if "M" not in sexes:
#            set_of_ploidy_pairs.drop(set_of_ploidy_pairs[set_of_ploidy_pairs['female_ploidies'] == 0].index, inplace = True)
#        female_ploidies     = list(set_of_ploidy_pairs.female_ploidies)
#        male_ploidies       = list(set_of_ploidy_pairs.male_ploidies)
#
#        for j in range(set_of_ploidy_pairs.shape[0]):
#            if female_ploidies[j] == male_ploidies[j]:
#                jobid_make_genDB_folder_and_map = f'make_genDB_folder_and_map_{gvcf_folder}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{female_ploidies[j]}'
#                gwf.target_from_template(jobid_make_genDB_folder_and_map,
#                                            make_genDB_folder_and_map(group     = ",".join(folders),
#                                                                    batch       = batch,
#                                                                    inds        = ",".join(inds),
#                                                                    ploidies    = ",".join(len(inds)*[str(female_ploidies[j])]),
#                                                                    folder_name = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/GenomicsDB/batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}/",
#                                                                    prev_done   = [f"/faststorage/project/megaFauna/sa_megafauna/data/{folders[jj]}/done/call_with_bed_{folders[jj]}_{inds[jj]}_batch_{batch}_ploidy_{female_ploidies[j]}" for jj in range(len(inds))] + [f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/{jobid_make_geno_metadata}"],
#                                                                    done        = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/{jobid_make_genDB_folder_and_map}"))
#
#            elif female_ploidies[j] == 2 and male_ploidies[j] == 1:
#                jobid_make_genDB_folder_and_map = f'make_genDB_folder_and_map_{gvcf_folder}_batch_{batch}_fploidy_2_mploidy_1'
#                gwf.target_from_template(jobid_make_genDB_folder_and_map,
#                                            make_genDB_folder_and_map(group     = ",".join(folders),
#                                                                    batch       = batch,
#                                                                    inds        = ",".join(inds),
#                                                                    ploidies    = ",".join(["2" if sexes[jj] == "F" else "1" for jj in range(len(sexes))]),
#                                                                    folder_name = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/GenomicsDB/batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}/",
#                                                                    prev_done   = [f"/faststorage/project/megaFauna/sa_megafauna/data/{folders[jj]}/done/call_with_bed_{folders[jj]}_{inds[jj]}_batch_{batch}_ploidy_2" if sexes[jj] == "F" else f"/faststorage/project/megaFauna/sa_megafauna/data/{folders[jj]}/done/call_with_bed_{folders[jj]}_{inds[jj]}_batch_{batch}_ploidy_1" for jj in range(len(inds))] + [f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/{jobid_make_geno_metadata}"],
#                                                                    done        = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/{jobid_make_genDB_folder_and_map}"))
#            elif female_ploidies[j] == 1 and male_ploidies[j] == 2:
#                jobid_make_genDB_folder_and_map = f'make_genDB_folder_and_map_{gvcf_folder}_batch_{batch}_fploidy_1_mploidy_2'
#                gwf.target_from_template(jobid_make_genDB_folder_and_map,
#                                            make_genDB_folder_and_map(group     = ",".join(folders),
#                                                                    batch       = batch,
#                                                                    inds        = ",".join(inds),
#                                                                    ploidies    = ",".join(["1" if sexes[jj] == "F" else "2" for jj in range(len(sexes))]),
#                                                                    folder_name = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/GenomicsDB/batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}/",
#                                                                    prev_done   = [f"/faststorage/project/megaFauna/sa_megafauna/data/{folders[jj]}/done/call_with_bed_{folders[jj]}_{inds[jj]}_batch_{batch}_ploidy_1" if sexes[jj] == "F" else f"/faststorage/project/megaFauna/sa_megafauna/data/{folders[jj]}/done/call_with_bed_{folders[jj]}_{inds[jj]}_batch_{batch}_ploidy_2" for jj in range(len(inds))] + [f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/{jobid_make_geno_metadata}"],
#                                                                    done        = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/{jobid_make_genDB_folder_and_map}"))
#
#            elif female_ploidies[j] == 0 and male_ploidies[j] == 1:
#                jobid_make_genDB_folder_and_map = f'make_genDB_folder_and_map_{gvcf_folder}_batch_{batch}_fploidy_0_mploidy_1'
#                gwf.target_from_template(jobid_make_genDB_folder_and_map,
#                                            make_genDB_folder_and_map(group     = ",".join([folders[jj] for jj in range(len(folders)) if sexes[jj] == "M"]),
#                                                                    batch       = batch,
#                                                                    inds        = ",".join([inds[jj] for jj in range(len(inds)) if sexes[jj] == "M"]),
#                                                                    ploidies    = ",".join(sexes.count("M")*["1"]),
#                                                                    folder_name = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/GenomicsDB/batch_{batch}_fploidy_0_mploidy_1/",
#                                                                    prev_done   = [f"/faststorage/project/megaFauna/sa_megafauna/data/{folders[jj]}/done/call_with_bed_{folders[jj]}_{inds[jj]}_batch_{batch}_ploidy_1" for jj in range(len(inds)) if sexes[jj] == "M"] + [f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/{jobid_make_geno_metadata}"],
#                                                                    done        = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/{jobid_make_genDB_folder_and_map}"))
#
 ##------------------------------------------ genotyping at species level -------------------------------------------###

#for gvcf_folder in species_and_refs.GVCF_FOLDER.drop_duplicates():
#    
#    subprocess.run(["mkdir", "-p", "/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done"])
#    subprocess.run(["mkdir", "-p", "/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/gVCF"])
#
#    length_of_chunk = 2000000
#    groups      = list(species_and_refs.loc[species_and_refs.GVCF_FOLDER == gvcf_folder].FOLDER)
#    rr_folders  = list(species_and_refs.loc[species_and_refs.GVCF_FOLDER == gvcf_folder].REFERENCE_FOLDER)
#
#    inds       = []
#    sexes      = []
#    folders    = []
#
#    ### gather inidividuals, sexes and done files
#    for i in range(len(groups)):
#        group           = groups[i]
#        ref_folder      = rr_folders[i]
#        ref_path        = f"/faststorage/project/megaFauna/sa_megafauna/data/{ref_folder}/ref/" + [file for file in os.listdir(f"/faststorage/project/megaFauna/sa_megafauna/data/{ref_folder}/ref/") if file.endswith("_LargerThan1000bp.fasta")][0]
#        inds_to_include = list(data.loc[data.FOLDER == group]["IND_ID"].drop_duplicates())
#        inds_and_sexes  = pd.read_table(f"/faststorage/project/megaFauna/sa_megafauna/data/{ref_folder}/ref/samples_coverage_stats.txt")
#        inds_and_sexes  = inds_and_sexes[inds_and_sexes.IND_ID.isin(inds_to_include)]
#        group_inds      = list(inds_and_sexes.IND_ID)
#        group_sexes     = list(inds_and_sexes.gSEX)
#
#        regions = pd.read_table(f"/faststorage/project/megaFauna/sa_megafauna/data/{ref_folder}/ref/batches_{ref_folder}.txt")
#        batches = list(regions.batch)
#
#        for j in range(len(group_inds)):
#            inds.append(group_inds[j])
#            sexes.append(group_sexes[j])
#            folders.append(group)
#    
#    #### Make GenomicDB workspaces, genotype and index
#    for batch in sorted(list(set(regions.batch))):
#        iis                 = [index for index in range(len(batches)) if batches[index] == batch]
#        set_of_ploidy_pairs = pd.DataFrame({"female_ploidies":[list(regions.female_ploidy)[ii] for ii in iis], "male_ploidies":[list(regions.male_ploidy)[ii] for ii in iis]}).drop_duplicates()
#        if "M" not in sexes:
#            set_of_ploidy_pairs.drop(set_of_ploidy_pairs[set_of_ploidy_pairs['female_ploidies'] == 0].index, inplace = True)
#        female_ploidies     = list(set_of_ploidy_pairs.female_ploidies)
#        male_ploidies       = list(set_of_ploidy_pairs.male_ploidies)
#
#        for j in range(set_of_ploidy_pairs.shape[0]):
#            if female_ploidies[j] == male_ploidies[j]:
#                jobid_make_genDB_folder_and_map = f'make_genDB_folder_and_map_{gvcf_folder}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{female_ploidies[j]}'
#            elif female_ploidies[j] == 2 and male_ploidies[j] == 1:
#                jobid_make_genDB_folder_and_map = f'make_genDB_folder_and_map_{gvcf_folder}_batch_{batch}_fploidy_2_mploidy_1'
#            elif female_ploidies[j] == 1 and male_ploidies[j] == 2:
#                jobid_make_genDB_folder_and_map = f'make_genDB_folder_and_map_{gvcf_folder}_batch_{batch}_fploidy_1_mploidy_2'
#            elif female_ploidies[j] == 0 and male_ploidies[j] == 1:
#                jobid_make_genDB_folder_and_map = f'make_genDB_folder_and_map_{gvcf_folder}_batch_{batch}_fploidy_0_mploidy_1'
#
#            prev_done_files = [f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/make_geno_metadata_{gvcf_folder}"]
#            files_to_concatenate = []
#            if os.path.isfile(f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/gVCF/geno_beds_and_intervals/{gvcf_folder}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_loc_{length_of_chunk}_subbatches.intervals"):
#                subbatch_file = pd.read_table(f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/gVCF/geno_beds_and_intervals/{gvcf_folder}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_loc_{length_of_chunk}_subbatches.intervals", header = None)
#                for sb in range(subbatch_file.shape[0]):
#                    ## E.3.1. Make a GenomicsDB workspace for each batch and ploidy category within batch; separate regions into chunks of length_of_chunk if the are longer than length_of_chunk
#                    jobid_make_genDB = f'make_genDB_{gvcf_folder}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_loc_{length_of_chunk}_subbatch_{sb}'
#                    gwf.target_from_template(jobid_make_genDB,
#                                                make_genDB_subbatch_with_bed(group     = gvcf_folder,
#                                                                            batch      = batch,
#                                                                            fploidy    = female_ploidies[j],
#                                                                            mploidy    = male_ploidies[j],
#                                                                            subbatch   = str(sb),
#                                                                            interval   = subbatch_file.iloc[sb,0],
#                                                                            prev_done  = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/{jobid_make_genDB_folder_and_map}",
#                                                                            done       = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/{jobid_make_genDB}"))
#
#                    ## E.4.1. Genotype for each batch and ploidy category within batch; separate regions into chunks of length_of_chunk if the are longer than length_of_chunk
#                    jobid_GenotypeGVCFs = f'GenotypeGVCFs_{gvcf_folder}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_loc_{length_of_chunk}_subbatch_{sb}'
#                    prev_done_files.append(f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/{jobid_GenotypeGVCFs}")
#                    files_to_concatenate.append(f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/gVCF/{gvcf_folder}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_subbatch_{sb}_gt.gvcf.gz")
#                    jobid_make_genDB = f'make_genDB_{gvcf_folder}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_loc_{length_of_chunk}_subbatch_{sb}'
#                    gwf.target_from_template(jobid_GenotypeGVCFs,
#                                            GenotypeGVCFs_subbatch_with_bed_new(group      = gvcf_folder,
#                                                                            batch      = batch,
#                                                                            fploidy    = female_ploidies[j],
#                                                                            mploidy    = male_ploidies[j],
#                                                                            subbatch   = str(sb),
#                                                                            interval   = subbatch_file.iloc[sb,0],
#                                                                            out        = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/gVCF/{gvcf_folder}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_subbatch_{sb}_gt.gvcf.gz",
#                                                                            ref        = ref_path,
#                                                                            cores      = 1,
#                                                                            memory     = "32g",
#                                                                            prev_done  = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/{jobid_make_genDB}",
#                                                                            done       = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/{jobid_GenotypeGVCFs}"))
#
#            if os.path.isfile(f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/gVCF/geno_beds_and_intervals/{gvcf_folder}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_loc_{length_of_chunk}_short_segments.intervals"): ## genotype smaller regions shorter than length_of_chunk (if they exist)
#                short_segments_limit = 1000
#                short_segments_count = sum(1 for _ in open(f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/gVCF/geno_beds_and_intervals/{gvcf_folder}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_loc_{length_of_chunk}_short_segments.intervals"))-1
#                if short_segments_count > short_segments_limit:
#                    to_s = [min(start + short_segments_limit - 1, short_segments_count) for start in range(1, short_segments_count + 1, short_segments_limit)]
#                    for sb in range(short_segments_count//short_segments_limit + 1):
#                        ## E.3.2. Make a GenomicsDB workspace for short segments subbatch
#                        jobid_make_genDB = f'make_genDB_{gvcf_folder}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_loc_{length_of_chunk}_short_segments_subbatch_{sb}'
#                        gwf.target_from_template(jobid_make_genDB,
#                                                    make_genDB_short_segments_subbatch(group      = gvcf_folder,
#                                                                                    batch       = batch,
#                                                                                    fploidy     = female_ploidies[j],
#                                                                                    mploidy     = male_ploidies[j],
#                                                                                    subbatch    = str(sb),
#                                                                                    intervals   = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/gVCF/geno_beds_and_intervals/{gvcf_folder}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_loc_{length_of_chunk}_short_segments.intervals",
#                                                                                    fr          = sb*short_segments_limit + 1,
#                                                                                    to          = to_s[sb],
#                                                                                    prev_done   = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/{jobid_make_genDB_folder_and_map}",
#                                                                                    done        = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/{jobid_make_genDB}"))
#
#                        ## E.4.2. Genotype short segments subbatch
#                        jobid_GenotypeGVCFs = f'GenotypeGVCFs_{gvcf_folder}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_short_segments_subbatch_{sb}'
#                        prev_done_files.append(f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/{jobid_GenotypeGVCFs}")
#                        files_to_concatenate.append(f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/gVCF/{gvcf_folder}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_short_segments_subbatch_{sb}_gt.gvcf.gz")
#                        gwf.target_from_template(jobid_GenotypeGVCFs,
#                                                GenotypeGVCFs_short_segments_subbatch(group       = gvcf_folder,
#                                                                                    batch       = batch,
#                                                                                    fploidy     = female_ploidies[j],
#                                                                                    mploidy     = male_ploidies[j],
#                                                                                    subbatch    = str(sb),
#                                                                                    intervals   = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/gVCF/geno_beds_and_intervals/{gvcf_folder}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_loc_{length_of_chunk}_short_segments.intervals",
#                                                                                    fr          = sb*short_segments_limit + 1,
#                                                                                    to          = to_s[sb],
#                                                                                    out         = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/gVCF/{gvcf_folder}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_short_segments_subbatch_{sb}_gt.gvcf.gz",
#                                                                                    ref         = ref_path,
#                                                                                    cores       = 1,
#                                                                                    memory      = "32g",
#                                                                                    prev_done   = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/{jobid_make_genDB}",
#                                                                                    done        = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/{jobid_GenotypeGVCFs}"))
#                else:
#                    ## E.3.2. Make a GenomicsDB workspace for short segments
#                    jobid_make_genDB = f'make_genDB_{gvcf_folder}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_loc_{length_of_chunk}_short_segments'
#                    gwf.target_from_template(jobid_make_genDB,
#                                                make_genDB_with_bed(group       = gvcf_folder,
#                                                                    batch       = batch,
#                                                                    fploidy     = female_ploidies[j],
#                                                                    mploidy     = male_ploidies[j],
#                                                                    intervals   = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/gVCF/geno_beds_and_intervals/{gvcf_folder}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_loc_{length_of_chunk}_short_segments.intervals",
#                                                                    prev_done   = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/{jobid_make_genDB_folder_and_map}",
#                                                                    done        = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/{jobid_make_genDB}"))
#
#                    ## E.4.2. Genotype short segments
#                    jobid_GenotypeGVCFs = f'GenotypeGVCFs_{gvcf_folder}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_short_segments'
#                    prev_done_files.append(f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/{jobid_GenotypeGVCFs}")
#                    files_to_concatenate.append(f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/gVCF/{gvcf_folder}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_short_segments_gt.gvcf.gz")
#                    gwf.target_from_template(jobid_GenotypeGVCFs,
#                                            GenotypeGVCFs_with_bed(group       = gvcf_folder,
#                                                                    batch       = batch,
#                                                                    fploidy     = female_ploidies[j],
#                                                                    mploidy     = male_ploidies[j],
#                                                                    intervals   = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/gVCF/geno_beds_and_intervals/{gvcf_folder}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_loc_{length_of_chunk}_short_segments.intervals",
#                                                                    out         = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/gVCF/{gvcf_folder}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_short_segments_gt.gvcf.gz",
#                                                                    ref         = ref_path,
#                                                                    cores       = 1,
#                                                                    memory      = "32g",
#                                                                    prev_done   = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/{jobid_make_genDB}",
#                                                                    done        = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/{jobid_GenotypeGVCFs}"))
#            
#            if len(files_to_concatenate) > 1: ## if chunking was done, concatenate all chunks within a batch and ploidy category, then index the concatenated gvcf and finally remove individual-based GVCFs
#                ## E.5.1. Concatenate
#                jobid_picardconcat = f'picardconcat_{gvcf_folder}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}'
#                gwf.target_from_template(jobid_picardconcat,
#                                            picardconcat(vcfs      = files_to_concatenate,
#                                                        vcf       = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/gVCF/{gvcf_folder}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_gt.gvcf.gz",
#                                                        prev_done = prev_done_files,
#                                                        done      = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/{jobid_picardconcat}"))
#
#                ## E.6. Index
#                jobid_index = f'IndexGVCFs_{gvcf_folder}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}'
#                gwf.target_from_template(jobid_index,
#                                        IndexGVCFs(group      = gvcf_folder,
#                                                    batch     = batch,
#                                                    fploidy   = female_ploidies[j],
#                                                    mploidy   = male_ploidies[j],
#                                                    prev_done = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/{jobid_picardconcat}",
#                                                    done      = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/{jobid_index}"))
#
#            else:  ## if chunking was NOT done, rename the gvcf to follow convention, then index the renamed gvcf and finally remove individual-based GVCFs
#                ## E.5.2. Rename
#                jobid_renameGVCF = f'renameGVCF_{gvcf_folder}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}'
#                gwf.target_from_template(jobid_renameGVCF,
#                                            renameGVCF(vcf_in     = files_to_concatenate[0],
#                                                        vcf_out   = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/gVCF/{gvcf_folder}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_gt.gvcf.gz",
#                                                        prev_done = prev_done_files,  
#                                                        done      = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/{jobid_renameGVCF}"))
#
#                ## E.6. Index
#                jobid_index = f'IndexGVCFs_{gvcf_folder}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}'
#                gwf.target_from_template(jobid_index,
#                                        IndexGVCFs(group     = gvcf_folder,
#                                                    batch     = batch,
#                                                    fploidy   = female_ploidies[j],
#                                                    mploidy   = male_ploidies[j],
#                                                    prev_done = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/{jobid_renameGVCF}",
#                                                    done      = f"/faststorage/project/megaFauna/sa_megafauna/data/{gvcf_folder}/done/{jobid_index}"))

# """
# ------------------------------------------------------------------------------------------------------------------------
# F. Filtering
# ------------------------------------------------------------------------------------------------------------------------
# """
#
for i in range(species_and_refs.shape[0]):
    # Initialising folders and variables for putting in the functions

    group      = species_and_refs.FOLDER[i]
    ref_folder = species_and_refs.REFERENCE_FOLDER[i]

    subprocess.run(["mkdir", "-p", f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/VCF/batches"])
    subprocess.run(["mkdir", "-p", f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/VCF/snp_counts"])

    ## F.1 filter on sample level 
    ''' Filters the samples_coverage_stats.txt based on the total coverage and the fraction of the genome covered'''

    inds_to_include = (data[data["FOLDER"] == group]["IND_ID"].drop_duplicates().tolist())
    coverage = pd.read_table(f"/faststorage/project/megaFauna/sa_megafauna/data/{ref_folder}/ref/samples_coverage_stats.txt")
    coverage = coverage[coverage.IND_ID.isin(inds_to_include)]
    coverage_filtered = coverage[(coverage['cov_A'] > 10) & (coverage['cov_len_A'] > 0.90)]
    if group == "Boselaphus_tragocamelus":
        coverage_filtered = coverage[(coverage['cov_A'] > 10) & (coverage['cov_len_A'] > 0.84)]
    with open(f"/faststorage/project/megaFauna/sa_megafauna/data/{ref_folder}/ref/samples_coverage_stats_{group}_filtered.txt", 'w') as f:
        coverage_filtered.to_csv(f, sep="\t", index=False)
 
    ind_list   = list(coverage_filtered["IND_ID"])

    gvcfs = [x for x in os.listdir(f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/gVCF/") 
        if x.endswith("gvcf.gz") and x.startswith(f"{group}")]
    batch_dones = []
    filtered_vcfs = []
    for gvcf in gvcfs:
        jobid_snp = f"{gvcf.split(".")[0]}_snps"
        batch_dones.append(f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_snp}")
        filtered_vcfs.append(f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/VCF/batches/{gvcf.split(".")[0]}_snps.vcf.gz")
        ## F.2 Extracting and filtering SNPs from gVCFs into VCFs
        gwf.target_from_template(
            jobid_snp,
            snps_filtering(
                inds        = ",".join(ind_list),
                gvcf_in     = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/gVCF/{gvcf}",
                vcf_out     = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/VCF/batches/{gvcf.split(".")[0]}_snps.vcf.gz",
                prev_done   = [],
                done        = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid_snp}"
                )
            )
                
    # F.4 SNP counting and more filtering

    def get_chrom_type_from_ploidy(filename):
        """"
        Get the chromosome type (autosomal, X, Y) from the file name.
        """
        match = re.search(r"_fploidy_(\d+)_mploidy_(\d+)_", filename)
        if not match:
            return "unknown"
        fploidy = int(match.group(1))
        mploidy = int(match.group(2))

        if fploidy == 2 and mploidy == 2:
            return "autosomal"
        elif fploidy == 2 and mploidy == 1:
            return "chrX"
        elif fploidy == 0 and mploidy == 1:
            return "chrY"
        elif fploidy == 1 and mploidy == 1:
            return "chrM"
        else:
            return "unknown"

    snp_count_files = []
    snp_count_dones = []

    for vcf in filtered_vcfs:
        chrom_type = get_chrom_type_from_ploidy(vcf)
        jobid = f"{vcf.split('/')[-1].split('.')[0]}_count"
        outfile = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/VCF/snp_counts/{jobid}.txt"
        done = f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{jobid}"

        gwf.target_from_template(
            jobid,
            count_snps_per_vcf(
            vcf=vcf,
            inds=",".join(ind_list),
            chrom_type=chrom_type,
            outfile=outfile,
            prev_done=batch_dones,
            done=done
        ))
        snp_count_files.append(outfile)
        snp_count_dones.append(done)

    gwf.target_from_template(
    f"{group}_merge_snp_counts",
    merge_snp_counts_by_type(
        files=snp_count_files,
        outfile=f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/VCF/snp_counts/merged_counts.txt",
        prev_done=snp_count_dones,
        done=f"/faststorage/project/megaFauna/sa_megafauna/data/{group}/done/{group}_merge_snp_counts"))