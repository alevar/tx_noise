#!/usr/bin/env python

import pandas as pd
import numpy as np
import subprocess
import argparse
import random
from scipy import stats
import glob
import math
import csv
import sys
import os

def wrapper(args):
    # declarations
    base_dir_data = args.basedirdata
    base_dir_out = args.basedirout
    out_dir = args.outdir
    hg38_fa = args.ref

    genRNAseq = args.genRNAseq

    readlen = args.readlen
    num_samples = args.numSamples # number of samples to simulate

    gff3cols=["seqid","source","type","start","end","score","strand","phase","attributes"]

    # load base annotations
    # print(">>>loading base annotations")
    # real_baseDF = pd.read_csv(base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.true.gtf",sep="\t",names=gff3cols)
    # real_baseDF = real_baseDF[real_baseDF["type"]=="transcript"].reset_index(drop=True)
    # nonint_baseDF = pd.read_csv(base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.no.contained.non.intronic.gtf",sep="\t",names=gff3cols)
    # nonint_baseDF = nonint_baseDF[nonint_baseDF["type"]=="transcript"].reset_index(drop=True)
    # int_baseDF = pd.read_csv(base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.no.contained.intronic.gtf",sep="\t",names=gff3cols)
    # int_baseDF = int_baseDF[int_baseDF["type"]=="transcript"].reset_index(drop=True)
    # pol_baseDF = pd.read_csv(base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.no.contained.RNApol.gtf",sep="\t",names=gff3cols)
    # pol_baseDF = pol_baseDF[pol_baseDF["type"]=="transcript"].reset_index(drop=True)

    # # get all loci and transcript IDs
    # print(">>>getting loci IDs")
    # real_baseDF["lid"] = real_baseDF["attributes"].str.split("gene_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    # nonint_baseDF["lid"] = nonint_baseDF["attributes"].str.split("gene_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    # int_baseDF["lid"] = int_baseDF["attributes"].str.split("gene_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    # pol_baseDF["lid"] = pol_baseDF["attributes"].str.split("gene_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    # real_locs = set(real_baseDF["lid"])
    # nonint_locs = set(nonint_baseDF["lid"])
    # int_locs = set(int_baseDF["lid"])
    # pol_locs = set(pol_baseDF["lid"])
    # print("number of real locs: "+str(len(real_locs)))
    # print("number of nonint locs: "+str(len(nonint_locs)))
    # print("number of int locs: "+str(len(int_locs)))
    # print("number of pol locs: "+str(len(pol_locs)))

    # total_n_locs = len(real_locs)+len(nonint_locs)+len(int_locs)+len(pol_locs)

    # # perform cleanup, by removing any loci in 
    # #   1. int that are also in real
    # #   2. nonint that are not in real
    # #   3. pol that are in real
    # int_locs = int_locs - real_locs.intersection(int_locs)
    # assert(len(real_locs.intersection(int_locs))==0),"something wrong intronic"
    # int_baseDF = int_baseDF[int_baseDF["lid"].isin(int_locs)].reset_index(drop=True)

    # nonint_locs = nonint_locs - nonint_locs.difference(real_locs)
    # assert(len(nonint_locs) == len(real_locs.intersection(nonint_locs))),"something wrong non-intronic"
    # nonint_baseDF = nonint_baseDF[nonint_baseDF["lid"].isin(nonint_locs)].reset_index(drop=True)

    # pol_locs = pol_locs - real_locs.intersection(pol_locs)
    # assert(len(real_locs.intersection(pol_locs))==0),"something wrong polymerase"
    # pol_baseDF = pol_baseDF[pol_baseDF["lid"].isin(pol_locs)].reset_index(drop=True)

    # print("number of real locs: "+str(len(real_locs)))
    # print("number of nonint locs: "+str(len(nonint_locs)))
    # print("number of int locs: "+str(len(int_locs)))
    # print("number of pol locs: "+str(len(pol_locs)))

    # real_baseDF["tid"] = real_baseDF["attributes"].str.split("transcript_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    # real_baseDF = real_baseDF[["lid","tid"]]
    # nonint_baseDF["tid"] = nonint_baseDF["attributes"].str.split("transcript_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    # nonint_baseDF = nonint_baseDF[["lid","tid"]]
    # int_baseDF["tid"] = int_baseDF["attributes"].str.split("transcript_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    # int_baseDF = int_baseDF[["lid","tid"]]
    # pol_baseDF["tid"] = pol_baseDF["attributes"].str.split("transcript_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    # pol_baseDF = pol_baseDF[["lid","tid"]]


    # # get the averages
    # avg_real_locs = 0
    # avg_noise_locs = 0
    # avg_undefined_locs = 0
    # with open(base_dir_out+"res_distrib.loc_stats","r") as inFP:
    #     for line in inFP.readlines():
    #         type_loc,num_loc = line.split(":")
    #         if type_loc=="real":
    #             avg_real_locs = int(num_loc)
    #         elif type_loc=="noise":
    #             avg_noise_locs = int(num_loc)
    #         elif type_loc=="undefined":
    #             if(int(num_loc))>0:
    #                 avg_undefined_locs = int(num_loc)
    #         else:
    #             print("error in loading averages")
    # print("average real locs: "+str(avg_real_locs))
    # print("average noise locs: "+str(avg_noise_locs))
    # print("average undefined locs: "+str(avg_undefined_locs))

    # # randomly select loci of each type based on the averages
    # # noise loci are a combination of: intronic and polymerase
    # # this needs to be done proportionally to the consistent parts of each group (real(real,nonint) and noise(int,pol))

    # real_only_locs = real_locs - nonint_locs

    # # in case of inconsistencies associated with merging transcripts of different types under the same locus
    # # we need to scale the averages with respect to the new numbers
    # avg_total_locs = avg_real_locs+avg_noise_locs
    # # now scale these averages such that
    # # new_avg_real_locs <= len(real_only_locs and nonint_locs)
    # # new_avg_noise_locs <= len(int_locs and pol_locs)
    # sum_locs_to_choose_from = len(real_only_locs)+len(nonint_locs)+len(int_locs)+len(pol_locs)
    # # new_avg_real_locs = int((len(nonint_locs)+len(real_only_locs))*(avg_real_locs/avg_total_locs))
    # # new_avg_noise_locs = int((len(pol_locs)+len(int_locs))*(avg_noise_locs/avg_total_locs))

    # print("old avg real locs: "+str(avg_real_locs))
    # # print("new avg real locs: "+str(new_avg_real_locs))
    # print("old avg noise locs: "+str(avg_noise_locs))
    # # print("new avg noise locs: "+str(new_avg_noise_locs))

    # print("real locs: "+str(len(real_only_locs)))
    # print("nonint locs: "+str(len(nonint_locs)))
    # print("int locs: "+str(len(int_locs)))
    # print("pol locs: "+str(len(pol_locs)))

    # perc_real = len(real_only_locs)/len(real_locs)
    # perc_nonint = len(nonint_locs)/len(real_locs)
    # assert(perc_real+perc_nonint)==1,"wrong percent real and nonint"

    # perc_int = len(int_locs)/(len(int_locs)+len(pol_locs))
    # perc_pol = len(pol_locs)/(len(int_locs)+len(pol_locs))
    # assert((perc_int+perc_pol)==1),"wrong percent pol and int"
    # print("selecting n real locs: "+str(int(avg_real_locs*perc_real)))
    # print("selecting n nonint locs: "+str(int(avg_real_locs*perc_nonint)))
    # print("selecting n int locs: "+str(int(avg_noise_locs*perc_int)))
    # print("selecting n pol locs: "+str(int(avg_noise_locs*perc_pol)))

    # real_only_locs_rand = np.random.choice(list(real_only_locs),int(avg_real_locs*perc_real), replace=False)
    # nonint_locs_rand = np.random.choice(list(nonint_locs),int(avg_real_locs*perc_nonint), replace=False)
    # int_locs_rand = np.random.choice(list(int_locs),int(avg_noise_locs*perc_int), replace=False)
    # pol_locs_rand = np.random.choice(list(pol_locs),int(avg_noise_locs*perc_pol), replace=False)

    # real_only_baseDF = real_baseDF[real_baseDF["lid"].isin(real_only_locs_rand)].reset_index(drop=True)
    # real_nonint_baseDF = real_baseDF[(real_baseDF["lid"].isin(nonint_locs_rand))&~(real_baseDF["lid"].isin(real_only_locs_rand))].reset_index(drop=True)
    # nonint_baseDF = nonint_baseDF[nonint_baseDF["lid"].isin(nonint_locs_rand)].reset_index(drop=True)
    # int_baseDF = int_baseDF[int_baseDF["lid"].isin(int_locs_rand)].reset_index(drop=True)
    # pol_baseDF = pol_baseDF[pol_baseDF["lid"].isin(pol_locs_rand)].reset_index(drop=True)

    # # load the transcript tissue distribution
    # ttx_real = list()
    # ttx_nonint = list()
    # ttx_int = list()
    # ttx_pol = list()

    # tmp_tx = pd.read_csv(base_dir_out+"res_distrib.num_tx_per_tissue_loc2")
    # tmp_tx_real_loc = tmp_tx[(tmp_tx["real"]>0)][["total","real","nonintronic"]].reset_index(drop=True)
    # tmp_tx_real_loc["perc_real"] = tmp_tx_real_loc["real"]/tmp_tx_real_loc["total"]
    # tmp_tx_real_loc["perc_nonint"] = tmp_tx_real_loc["nonintronic"]/tmp_tx_real_loc["total"]
    # d_real = list(tmp_tx_real_loc["perc_real"])
    # d_nonint = list(tmp_tx_real_loc["perc_nonint"])

    # tmp_tx_int_loc = tmp_tx[(tmp_tx["intronic"]>0)][["total","intronic"]].reset_index(drop=True)
    # tmp_tx_int_loc["perc_int"] = tmp_tx_int_loc["intronic"]/tmp_tx_int_loc["total"]
    # d_int = list(tmp_tx_int_loc["perc_int"])

    # tmp_tx_pol_loc = tmp_tx[(tmp_tx["polymerase"]>0)][["total","polymerase"]].reset_index(drop=True)
    # tmp_tx_pol_loc["perc_pol"] = tmp_tx_pol_loc["polymerase"]/tmp_tx_pol_loc["total"]
    # d_pol = list(tmp_tx_pol_loc["perc_pol"])

    # # now that we have the loci, we can subset the gtfs and move on to working with the transcripts
    # # the question is: how do we allocate transcripts per locus?
    # # do we do this based on the distribution of the number of transcripts per each type of locus?
    # # are there any caveates in the co-dependencies such as between real and non-intronic?

    # # we should probably only select transcripts of each group within each type of locus dataframe

    # # one important consideration is that the total number of transcripts should be close to what is observed
    # #  in an average tissue (needs to be computed in the gtex_stats)
    # # this approach, however, might be overly complicated, and simply randomly selecting n transcripts for each group
    # #  based on distribution computed by gtex_stats might be sufficient, and should approximate well to desired values

    # real_only_txs = set()
    # for name, group in real_only_baseDF.groupby("lid"):
    #     tids = group["tid"].tolist()
    #     nt = math.ceil(d_real[random.randint(0,len(d_real)-1)]*len(tids))
    #     random.shuffle(tids)
    #     real_only_txs.update(tids[:nt])

    # real_nonint_txs = set()
    # for name, group in real_nonint_baseDF.groupby("lid"):
    #     tids = group["tid"].tolist()
    #     nt = math.ceil(d_real[random.randint(0,len(d_real)-1)]*len(tids))
    #     random.shuffle(tids)
    #     real_nonint_txs.update(tids[:nt])
        
    # nonint_txs = set()
    # for name, group in nonint_baseDF.groupby("lid"):
    #     tids = group["tid"].tolist()
    #     nt = math.ceil(d_nonint[random.randint(0,len(d_nonint)-1)]*len(tids))
    #     random.shuffle(tids)
    #     nonint_txs.update(tids[:nt])
        
    # int_txs = set()
    # for name, group in int_baseDF.groupby("lid"):
    #     tids = group["tid"].tolist()
    #     nt = math.ceil(d_int[random.randint(0,len(d_int)-1)]*len(tids))
    #     random.shuffle(tids)
    #     int_txs.update(tids[:nt])
        
    # pol_txs = set()
    # for name, group in pol_baseDF.groupby("lid"):
    #     tids = group["tid"].tolist()
    #     nt = math.ceil(d_pol[random.randint(0,len(d_pol)-1)]*len(tids))
    #     random.shuffle(tids)
    #     pol_txs.update(tids[:nt])

    # print("number of real only txs: "+str(len(real_only_txs)))
    # print("number of real nonint txs: "+str(len(real_nonint_txs)))
    # print("number of nonint txs: "+str(len(nonint_txs)))
    # print("number of int txs: "+str(len(int_txs)))
    # print("number of pol txs: "+str(len(pol_txs)))

    # # now that we have lists of transcripts - we can subset the annotation further
    # print("number of real only txs before: "+str(len(real_only_baseDF)))
    # real_only_baseDF = real_only_baseDF[real_only_baseDF["tid"].isin(real_only_txs)].reset_index(drop=True)
    # print("number of real nonint txs before: "+str(len(real_nonint_baseDF)))
    # real_nonint_baseDF = real_nonint_baseDF[real_nonint_baseDF["tid"].isin(real_nonint_txs)].reset_index(drop=True)
    # print("number of nonint txs before: "+str(len(nonint_baseDF)))
    # nonint_baseDF = nonint_baseDF[nonint_baseDF["tid"].isin(nonint_txs)].reset_index(drop=True)
    # print("number of int txs before: "+str(len(int_baseDF)))
    # int_baseDF = int_baseDF[int_baseDF["tid"].isin(int_txs)].reset_index(drop=True)
    # print("number of pol txs before: "+str(len(pol_baseDF)))
    # pol_baseDF = pol_baseDF[pol_baseDF["tid"].isin(pol_txs)].reset_index(drop=True)

    # # lastly we need to generate samples based on the underlying annotation
    # # for this we will need additional information, similar to the tissue_per_loc2 in order to tell
    # # how likely a sample is to contribute a transcript to a tissue

    # # load the transcript sample distribution
    # ttx_real = list()
    # ttx_nonint = list()
    # ttx_int = list()
    # ttx_pol = list()

    # tmp_tx = pd.read_csv(base_dir_out+"res_distrib.num_tx_per_sample_loc2")
    # tmp_tx_real_loc = tmp_tx[(tmp_tx["real"]>0)][["total","real","nonintronic"]].reset_index(drop=True)
    # tmp_tx_real_loc["perc_real"] = tmp_tx_real_loc["real"]/tmp_tx_real_loc["total"]
    # tmp_tx_real_loc["perc_nonint"] = tmp_tx_real_loc["nonintronic"]/tmp_tx_real_loc["total"]
    # d_real = list(tmp_tx_real_loc["perc_real"])
    # d_nonint = list(tmp_tx_real_loc["perc_nonint"])

    # tmp_tx_int_loc = tmp_tx[(tmp_tx["intronic"]>0)][["total","intronic"]].reset_index(drop=True)
    # tmp_tx_int_loc["perc_int"] = tmp_tx_int_loc["intronic"]/tmp_tx_int_loc["total"]
    # d_int = list(tmp_tx_int_loc["perc_int"])

    # tmp_tx_pol_loc = tmp_tx[(tmp_tx["polymerase"]>0)][["total","polymerase"]].reset_index(drop=True)
    # tmp_tx_pol_loc["perc_pol"] = tmp_tx_pol_loc["polymerase"]/tmp_tx_pol_loc["total"]
    # d_pol = list(tmp_tx_pol_loc["perc_pol"])

    # real_only_txs = [set() for x in range(num_samples)]
    # for name, group in real_only_baseDF.groupby("lid"):
    #     tids = group["tid"].tolist()
    #     for i in range(num_samples):
    #         nt = math.ceil(d_real[random.randint(0,len(d_real)-1)]*len(tids))
    #         random.shuffle(tids)
    #         real_only_txs[i].update(tids[:nt])

    # real_nonint_txs = [set() for x in range(num_samples)]
    # for name, group in real_nonint_baseDF.groupby("lid"):
    #     tids = group["tid"].tolist()
    #     for i in range(num_samples):
    #         nt = math.ceil(d_real[random.randint(0,len(d_real)-1)]*len(tids))
    #         random.shuffle(tids)
    #         real_nonint_txs[i].update(tids[:nt])

    # nonint_txs = [set() for x in range(num_samples)]
    # for name, group in nonint_baseDF.groupby("lid"):
    #     tids = group["tid"].tolist()
    #     for i in range(num_samples):
    #         nt = math.ceil(d_nonint[random.randint(0,len(d_nonint)-1)]*len(tids))
    #         random.shuffle(tids)
    #         nonint_txs[i].update(tids[:nt])

    # int_txs = [set() for x in range(num_samples)]
    # for name, group in int_baseDF.groupby("lid"):
    #     tids = group["tid"].tolist()
    #     for i in range(num_samples):
    #         nt = math.ceil(d_int[random.randint(0,len(d_int)-1)]*len(tids))
    #         random.shuffle(tids)
    #         int_txs[i].update(tids[:nt])
            
    # pol_txs = [set() for x in range(num_samples)]
    # for name, group in pol_baseDF.groupby("lid"):
    #     tids = group["tid"].tolist()
    #     for i in range(num_samples):
    #         nt = math.ceil(d_pol[random.randint(0,len(d_pol)-1)]*len(tids))
    #         random.shuffle(tids)
    #         pol_txs[i].update(tids[:nt])

    # # now that we have sample specific stuff - we can generate sample-specific GTFs
    # real_baseDF = pd.read_csv(base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.true.gtf",sep="\t",names=gff3cols)
    # nonint_baseDF = pd.read_csv(base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.no.contained.non.intronic.gtf",sep="\t",names=gff3cols)
    # int_baseDF = pd.read_csv(base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.no.contained.intronic.gtf",sep="\t",names=gff3cols)
    # pol_baseDF = pd.read_csv(base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.no.contained.RNApol.gtf",sep="\t",names=gff3cols)
    # real_baseDF["tid"] = real_baseDF["attributes"].str.split("transcript_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    # nonint_baseDF["tid"] = nonint_baseDF["attributes"].str.split("transcript_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    # int_baseDF["tid"] = int_baseDF["attributes"].str.split("transcript_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    # pol_baseDF["tid"] = pol_baseDF["attributes"].str.split("transcript_id \"",expand=True)[1].str.split("\"",expand=True)[0]

    # for i in range(num_samples):
    #     real_baseDF[real_baseDF["tid"].isin(real_only_txs[i].union(real_nonint_txs[i]))][gff3cols].to_csv(out_dir+"/res_distrib.real.sample"+str(i)+".gtf",sep="\t",index=False,header=False,quoting=csv.QUOTE_NONE)
    #     nonint_baseDF[nonint_baseDF["tid"].isin(nonint_txs[i])][gff3cols].to_csv(out_dir+"/res_distrib.nonint.sample"+str(i)+".gtf",sep="\t",index=False,header=False,quoting=csv.QUOTE_NONE)
    #     int_baseDF[int_baseDF["tid"].isin(int_txs[i])][gff3cols].to_csv(out_dir+"/res_distrib.int.sample"+str(i)+".gtf",sep="\t",index=False,header=False,quoting=csv.QUOTE_NONE)
    #     pol_baseDF[pol_baseDF["tid"].isin(pol_txs[i])][gff3cols].to_csv(out_dir+"/res_distrib.pol.sample"+str(i)+".gtf",sep="\t",index=False,header=False,quoting=csv.QUOTE_NONE)

    # # load distribution of the sum of expressions
    # tmp_exp_real = pd.read_csv(base_dir_out+"res_distrib.cov_sample_real")
    # d_exp_real = list(tmp_exp_real[tmp_exp_real["cov"]>=1]["cov"])

    # tmp_exp_nonint = pd.read_csv(base_dir_out+"res_distrib.cov_sample_nonintronic")
    # d_exp_nonint = list(tmp_exp_nonint[tmp_exp_nonint["cov"]>=1]["cov"])

    # tmp_exp_int = pd.read_csv(base_dir_out+"res_distrib.cov_sample_intronic")
    # d_exp_int = list(tmp_exp_int[tmp_exp_int["cov"]>=1]["cov"])

    # tmp_exp_pol = pd.read_csv(base_dir_out+"res_distrib.cov_sample_polymerase")
    # d_exp_pol = list(tmp_exp_pol[tmp_exp_pol["cov"]>=1]["cov"])

    # # now generate corresponding expression matrices for the transcripts in each sample for each type
    # for i in range(num_samples):
    #     realDF = pd.read_csv(out_dir+"/res_distrib.real.sample"+str(i)+".gtf",sep="\t",names=gff3cols)
    #     nt_real = len(realDF[realDF["type"]=="transcript"])
    #     exps_real = np.random.choice(d_exp_real, nt_real, replace=False)
    #     np.savetxt(out_dir+"/res_distrib.real.sample"+str(i)+".exp", exps_real, delimiter=',',fmt='%0.2f')
        
    #     nonintDF = pd.read_csv(out_dir+"/res_distrib.nonint.sample"+str(i)+".gtf",sep="\t",names=gff3cols)
    #     nt_nonint = len(nonintDF[nonintDF["type"]=="transcript"])
    #     exps_nonint = np.random.choice(d_exp_nonint, nt_nonint, replace=False)
    #     np.savetxt(out_dir+"/res_distrib.nonint.sample"+str(i)+".exp", exps_nonint, delimiter=',',fmt='%0.2f')
        
    #     intDF = pd.read_csv(out_dir+"/res_distrib.int.sample"+str(i)+".gtf",sep="\t",names=gff3cols)
    #     nt_int = len(intDF[intDF["type"]=="transcript"])
    #     exps_int = np.random.choice(d_exp_int, nt_int, replace=False)
    #     np.savetxt(out_dir+"/res_distrib.int.sample"+str(i)+".exp", exps_int, delimiter=',',fmt='%0.2f')
        
    #     polDF = pd.read_csv(out_dir+"/res_distrib.pol.sample"+str(i)+".gtf",sep="\t",names=gff3cols)
    #     nt_pol = len(polDF[polDF["type"]=="transcript"])
    #     exps_pol = np.random.choice(d_exp_pol, nt_pol, replace=False)
    #     np.savetxt(out_dir+"/res_distrib.pol.sample"+str(i)+".exp", exps_pol, delimiter=',',fmt='%0.2f')

    # # generate transcriptomes for each file with gffread
    # for i in range(num_samples):
    #     cmd_real = ["gffread",
    #                 "-w",out_dir+"res_distrib.real.sample"+str(i)+".fasta",
    #                 "-g",hg38_fa,
    #                 out_dir+"res_distrib.real.sample"+str(i)+".gtf"]
    #     subprocess.call(cmd_real)
        
    #     cmd_nonint = ["gffread",
    #                 "-w",out_dir+"res_distrib.nonint.sample"+str(i)+".fasta",
    #                 "-g",hg38_fa,
    #                 out_dir+"res_distrib.nonint.sample"+str(i)+".gtf"]
    #     subprocess.call(cmd_nonint)
        
    #     cmd_int = ["gffread",
    #                 "-w",out_dir+"res_distrib.int.sample"+str(i)+".fasta",
    #                 "-g",hg38_fa,
    #                 out_dir+"res_distrib.int.sample"+str(i)+".gtf"]
    #     subprocess.call(cmd_int)
        
    #     cmd_pol = ["gffread",
    #                 "-w",out_dir+"res_distrib.pol.sample"+str(i)+".fasta",
    #                 "-g",hg38_fa,
    #                 out_dir+"res_distrib.pol.sample"+str(i)+".gtf"]
    #     subprocess.call(cmd_pol)

    # # also extract fasta from the set of all real transcripts for the salmon/kallisto indices
    # cmd_all = ["gffread",
    #            "-w",out_dir+"res_distrib.all.real.fasta",
    #            "-g",hg38_fa,
    #            base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.true.gtf"]
    # subprocess.call(cmd_all)

    # for i in range(num_samples):
    #     if not os.path.exists(out_dir+"res_distrib.real.sample"+str(i)):
    #         os.makedirs(out_dir+"res_distrib.real.sample"+str(i))
    #     cmd_real = [genRNAseq,
    #                 out_dir+"res_distrib.real.sample"+str(i)+".fasta",
    #                 out_dir+"res_distrib.real.sample"+str(i)+".exp",
    #                 str(readlen),
    #                 out_dir+"res_distrib.real.sample"+str(i)+"/"]
    #     subprocess.call(cmd_real)
        
    #     if not os.path.exists(out_dir+"res_distrib.nonint.sample"+str(i)):
    #         os.makedirs(out_dir+"res_distrib.nonint.sample"+str(i))
    #     cmd_nonint = [genRNAseq,
    #                 out_dir+"res_distrib.nonint.sample"+str(i)+".fasta",
    #                 out_dir+"res_distrib.nonint.sample"+str(i)+".exp",
    #                 str(readlen),
    #                 out_dir+"res_distrib.nonint.sample"+str(i)+"/"]
    #     subprocess.call(cmd_nonint)
        
    #     if not os.path.exists(out_dir+"res_distrib.int.sample"+str(i)):
    #         os.makedirs(out_dir+"res_distrib.int.sample"+str(i))
    #     cmd_int = [genRNAseq,
    #                 out_dir+"res_distrib.int.sample"+str(i)+".fasta",
    #                 out_dir+"res_distrib.int.sample"+str(i)+".exp",
    #                 str(readlen),
    #                 out_dir+"res_distrib.int.sample"+str(i)+"/"]
    #     subprocess.call(cmd_int)
        
    #     if not os.path.exists(out_dir+"res_distrib.pol.sample"+str(i)):
    #         os.makedirs(out_dir+"res_distrib.pol.sample"+str(i))
    #     cmd_pol = [genRNAseq,
    #                 out_dir+"res_distrib.pol.sample"+str(i)+".fasta",
    #                 out_dir+"res_distrib.pol.sample"+str(i)+".exp",
    #                 str(readlen),
    #                 out_dir+"res_distrib.pol.sample"+str(i)+"/"]
    #     subprocess.call(cmd_pol)

    # # build salmon and hisat2 indices
    # f=open(out_dir+"res_distrib.all.real.exons","wb")
    # subprocess.call([args.hisat+"_extract_exons.py",
    #                  base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.true.gtf"],
    #                  stdout=f)
    # f.close()
    # f=open(out_dir+"res_distrib.all.real.ss","wb")
    # subprocess.call([args.hisat+"_extract_splice_sites.py",
    #                  base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.true.gtf"],
    #                  stdout=f)
    # f.close()
    # subprocess.call([args.hisat+"-build",
    #                  args.ref,
    #                  out_dir+"res_distrib.all.real",
    #                  "-p",str(args.threads),
    #                  "--ss",out_dir+"res_distrib.all.real.ss",
    #                  "--exon",out_dir+"res_distrib.all.real.exons"])

    # subprocess.call(["salmon","index",
    #                  "-t",out_dir+"res_distrib.all.real.fasta",
    #                  "-p",str(args.threads),
    #                  "-i",out_dir+"res_distrib.all.real"])

    # # shuffle generated reads and create merged datasets for downstream analysis
    # print(">>>shuffling reads and creating combinations")
    # for i in range(num_samples):
    #     print(i)
    #     # shuffle real
    #     shuffle_cmd = [args.shuffleReads,
    #                    out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_1.fasta",
    #                    out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_2.fasta",
    #                    out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_1.shuffled.fasta",
    #                    out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_2.shuffled.fasta"]
    #     subprocess.call(shuffle_cmd)
    #     # create a combination of real and nonint
    #     if not os.path.exists(out_dir+"res_distrib.real_nonint.sample"+str(i)):
    #         os.makedirs(out_dir+"res_distrib.real_nonint.sample"+str(i))
    #     with open(out_dir+"res_distrib.real_nonint.sample"+str(i)+"/sample_01_1.fasta", 'w+') as outfile:
    #         for fname in [out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_1.fasta",
    #                       out_dir+"res_distrib.nonint.sample"+str(i)+"/sample_01_1.fasta"]:
    #             with open(fname) as infile:
    #                 for line in infile:
    #                     outfile.write(line)
    #     with open(out_dir+"res_distrib.real_nonint.sample"+str(i)+"/sample_01_2.fasta", 'w+') as outfile:
    #         for fname in [out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_2.fasta",
    #                       out_dir+"res_distrib.nonint.sample"+str(i)+"/sample_01_2.fasta"]:
    #             with open(fname) as infile:
    #                 for line in infile:
    #                     outfile.write(line)
    #     shuffle_cmd = [args.shuffleReads,
    #                    out_dir+"res_distrib.real_nonint.sample"+str(i)+"/sample_01_1.fasta",
    #                    out_dir+"res_distrib.real_nonint.sample"+str(i)+"/sample_01_2.fasta",
    #                    out_dir+"res_distrib.real_nonint.sample"+str(i)+"/sample_01_1.shuffled.fasta",
    #                    out_dir+"res_distrib.real_nonint.sample"+str(i)+"/sample_01_2.shuffled.fasta"]
    #     subprocess.call(shuffle_cmd)
    #     # create a combination of real and int
    #     if not os.path.exists(out_dir+"res_distrib.real_int.sample"+str(i)):
    #         os.makedirs(out_dir+"res_distrib.real_int.sample"+str(i))
    #     with open(out_dir+"res_distrib.real_int.sample"+str(i)+"/sample_01_1.fasta", 'w+') as outfile:
    #         for fname in [out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_1.fasta",
    #                       out_dir+"res_distrib.int.sample"+str(i)+"/sample_01_1.fasta"]:
    #             with open(fname) as infile:
    #                 for line in infile:
    #                     outfile.write(line)
    #     with open(out_dir+"res_distrib.real_int.sample"+str(i)+"/sample_01_2.fasta", 'w+') as outfile:
    #         for fname in [out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_2.fasta",
    #                       out_dir+"res_distrib.int.sample"+str(i)+"/sample_01_2.fasta"]:
    #             with open(fname) as infile:
    #                 for line in infile:
    #                     outfile.write(line)
    #     shuffle_cmd = [args.shuffleReads,
    #                    out_dir+"res_distrib.real_int.sample"+str(i)+"/sample_01_1.fasta",
    #                    out_dir+"res_distrib.real_int.sample"+str(i)+"/sample_01_2.fasta",
    #                    out_dir+"res_distrib.real_int.sample"+str(i)+"/sample_01_1.shuffled.fasta",
    #                    out_dir+"res_distrib.real_int.sample"+str(i)+"/sample_01_2.shuffled.fasta"]
    #     subprocess.call(shuffle_cmd)
    #     # create a combination of real and polymerase
    #     if not os.path.exists(out_dir+"res_distrib.real_pol.sample"+str(i)):
    #         os.makedirs(out_dir+"res_distrib.real_pol.sample"+str(i))
    #     with open(out_dir+"res_distrib.real_pol.sample"+str(i)+"/sample_01_1.fasta", 'w+') as outfile:
    #         for fname in [out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_1.fasta",
    #                       out_dir+"res_distrib.pol.sample"+str(i)+"/sample_01_1.fasta"]:
    #             with open(fname) as infile:
    #                 for line in infile:
    #                     outfile.write(line)
    #     with open(out_dir+"res_distrib.real_pol.sample"+str(i)+"/sample_01_2.fasta", 'w+') as outfile:
    #         for fname in [out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_2.fasta",
    #                       out_dir+"res_distrib.pol.sample"+str(i)+"/sample_01_2.fasta"]:
    #             with open(fname) as infile:
    #                 for line in infile:
    #                     outfile.write(line)
    #     shuffle_cmd = [args.shuffleReads,
    #                    out_dir+"res_distrib.real_pol.sample"+str(i)+"/sample_01_1.fasta",
    #                    out_dir+"res_distrib.real_pol.sample"+str(i)+"/sample_01_2.fasta",
    #                    out_dir+"res_distrib.real_pol.sample"+str(i)+"/sample_01_1.shuffled.fasta",
    #                    out_dir+"res_distrib.real_pol.sample"+str(i)+"/sample_01_2.shuffled.fasta"]
    #     subprocess.call(shuffle_cmd)
    #     # create a combination of all reads
    #     if not os.path.exists(out_dir+"res_distrib.all.sample"+str(i)):
    #         os.makedirs(out_dir+"res_distrib.all.sample"+str(i))
    #     with open(out_dir+"res_distrib.all.sample"+str(i)+"/sample_01_1.fasta", 'w+') as outfile:
    #         for fname in [out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_1.fasta",
    #                       out_dir+"res_distrib.nonint.sample"+str(i)+"/sample_01_1.fasta",
    #                       out_dir+"res_distrib.int.sample"+str(i)+"/sample_01_1.fasta",
    #                       out_dir+"res_distrib.pol.sample"+str(i)+"/sample_01_1.fasta"]:
    #             with open(fname) as infile:
    #                 for line in infile:
    #                     outfile.write(line)
    #     with open(out_dir+"res_distrib.all.sample"+str(i)+"/sample_01_2.fasta", 'w+') as outfile:
    #         for fname in [out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_2.fasta",
    #                       out_dir+"res_distrib.nonint.sample"+str(i)+"/sample_01_2.fasta",
    #                       out_dir+"res_distrib.int.sample"+str(i)+"/sample_01_2.fasta",
    #                       out_dir+"res_distrib.pol.sample"+str(i)+"/sample_01_2.fasta"]:
    #             with open(fname) as infile:
    #                 for line in infile:
    #                     outfile.write(line)
    #     shuffle_cmd = [args.shuffleReads,
    #                    out_dir+"res_distrib.all.sample"+str(i)+"/sample_01_1.fasta",
    #                    out_dir+"res_distrib.all.sample"+str(i)+"/sample_01_2.fasta",
    #                    out_dir+"res_distrib.all.sample"+str(i)+"/sample_01_1.shuffled.fasta",
    #                    out_dir+"res_distrib.all.sample"+str(i)+"/sample_01_2.shuffled.fasta"]
    #     subprocess.call(shuffle_cmd)

    # now align with hisat, assemble with stringtie and quantify with salmon
    print(">>>aligning,assemblying,quantifying")
    for i in range(num_samples):
        print(i)
        # firtst perform analysis using real reads only
        if not os.path.exists(out_dir+"strg.real.sample"+str(i)):
            os.makedirs(out_dir+"strg.real.sample"+str(i))
        # align
        hisat_cmd = [args.hisat,
                     "-x",out_dir+"res_distrib.all.real",
                     "-p",str(args.threads),
                     "--rna-sensitive","-f",
                     "-1",out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_1.shuffled.fasta",
                     "-2",out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_2.shuffled.fasta",
                     "-S",out_dir+"strg.real.sample"+str(i)+"/hisat.sam"]
        subprocess.call(hisat_cmd)
        # sort
        sort_cmd = ["samtools","sort",
                    "-@",str(args.threads),
                    "-o",out_dir+"strg.real.sample"+str(i)+"/hisat.sorted.bam",
                    out_dir+"strg.real.sample"+str(i)+"/hisat.sam"]
        subprocess.call(sort_cmd)
        # assemble
        strg_cmd = ["stringtie",
                    out_dir+"strg.real.sample"+str(i)+"/hisat.sorted.bam",
                    "-t",
                    "-p",str(args.threads),
                    "-G",base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.true.gtf",
                    "-o",out_dir+"strg.real.sample"+str(i)+"/strg.gtf"]
        subprocess.call(strg_cmd)
        # quantify with salmon
        salmon_cmd = ["salmon","quant","--validateMappings","-l","A",
                      "-i",out_dir+"res_distrib.all.real",
                      "-p",str(args.threads),
                      "-1",out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_1.shuffled.fasta",
                      "-2",out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_2.shuffled.fasta",
                      "-o",out_dir+"slmn.real.sample"+str(i)]
        subprocess.call(salmon_cmd)

        # next perform the same analysis with real and non-intronic reads
        if not os.path.exists(out_dir+"strg.real_nonint.sample"+str(i)):
            os.makedirs(out_dir+"strg.real_nonint.sample"+str(i))
        # align
        hisat_cmd = [args.hisat,
                     "-x",out_dir+"res_distrib.all.real",
                     "-p",str(args.threads),
                     "--rna-sensitive","-f",
                     "-1",out_dir+"res_distrib.real_nonint.sample"+str(i)+"/sample_01_1.shuffled.fasta",
                     "-2",out_dir+"res_distrib.real_nonint.sample"+str(i)+"/sample_01_2.shuffled.fasta",
                     "-S",out_dir+"strg.real_nonint.sample"+str(i)+"/hisat.sam"]
        subprocess.call(hisat_cmd)
        # sort
        sort_cmd = ["samtools","sort",
                    "-@",str(args.threads),
                    "-o",out_dir+"strg.real_nonint.sample"+str(i)+"/hisat.sorted.bam",
                    out_dir+"strg.real_nonint.sample"+str(i)+"/hisat.sam"]
        subprocess.call(sort_cmd)
        # assemble
        strg_cmd = ["stringtie",
                    out_dir+"strg.real_nonint.sample"+str(i)+"/hisat.sorted.bam",
                    "-t",
                    "-p",str(args.threads),
                    "-G",base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.true.gtf",
                    "-o",out_dir+"strg.real_nonint.sample"+str(i)+"/strg.gtf"]
        subprocess.call(strg_cmd)
        # quantify with salmon
        salmon_cmd = ["salmon","quant","--validateMappings","-l","A",
                      "-i",out_dir+"res_distrib.all.real",
                      "-p",str(args.threads),
                      "-1",out_dir+"res_distrib.real_nonint.sample"+str(i)+"/sample_01_1.shuffled.fasta",
                      "-2",out_dir+"res_distrib.real_nonint.sample"+str(i)+"/sample_01_2.shuffled.fasta",
                      "-o",out_dir+"slmn.real_nonint.sample"+str(i)]
        subprocess.call(salmon_cmd)

        # next perform the same analysis with real and intronic reads
        if not os.path.exists(out_dir+"strg.real_int.sample"+str(i)):
            os.makedirs(out_dir+"strg.real_int.sample"+str(i))
        # align
        hisat_cmd = [args.hisat,
                     "-x",out_dir+"res_distrib.all.real",
                     "-p",str(args.threads),
                     "--rna-sensitive","-f",
                     "-1",out_dir+"res_distrib.real_int.sample"+str(i)+"/sample_01_1.shuffled.fasta",
                     "-2",out_dir+"res_distrib.real_int.sample"+str(i)+"/sample_01_2.shuffled.fasta",
                     "-S",out_dir+"strg.real_int.sample"+str(i)+"/hisat.sam"]
        subprocess.call(hisat_cmd)
        # sort
        sort_cmd = ["samtools","sort",
                    "-@",str(args.threads),
                    "-o",out_dir+"strg.real_int.sample"+str(i)+"/hisat.sorted.bam",
                    out_dir+"strg.real_int.sample"+str(i)+"/hisat.sam"]
        subprocess.call(sort_cmd)
        # assemble
        strg_cmd = ["stringtie",
                    out_dir+"strg.real_int.sample"+str(i)+"/hisat.sorted.bam",
                    "-t",
                    "-p",str(args.threads),
                    "-G",base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.true.gtf",
                    "-o",out_dir+"strg.real_int.sample"+str(i)+"/strg.gtf"]
        subprocess.call(strg_cmd)
        # quantify with salmon
        salmon_cmd = ["salmon","quant","--validateMappings","-l","A",
                      "-i",out_dir+"res_distrib.all.real",
                      "-p",str(args.threads),
                      "-1",out_dir+"res_distrib.real_int.sample"+str(i)+"/sample_01_1.shuffled.fasta",
                      "-2",out_dir+"res_distrib.real_int.sample"+str(i)+"/sample_01_2.shuffled.fasta",
                      "-o",out_dir+"slmn.real_int.sample"+str(i)]
        subprocess.call(salmon_cmd)

        # next perform the same analysis with real and polymerase reads
        if not os.path.exists(out_dir+"strg.real_pol.sample"+str(i)):
            os.makedirs(out_dir+"strg.real_pol.sample"+str(i))
        # align
        hisat_cmd = [args.hisat,
                     "-x",out_dir+"res_distrib.all.real",
                     "-p",str(args.threads),
                     "--rna-sensitive","-f",
                     "-1",out_dir+"res_distrib.real_pol.sample"+str(i)+"/sample_01_1.shuffled.fasta",
                     "-2",out_dir+"res_distrib.real_pol.sample"+str(i)+"/sample_01_2.shuffled.fasta",
                     "-S",out_dir+"strg.real_pol.sample"+str(i)+"/hisat.sam"]
        subprocess.call(hisat_cmd)
        # sort
        sort_cmd = ["samtools","sort",
                    "-@",str(args.threads),
                    "-o",out_dir+"strg.real_pol.sample"+str(i)+"/hisat.sorted.bam",
                    out_dir+"strg.real_pol.sample"+str(i)+"/hisat.sam"]
        subprocess.call(sort_cmd)
        # assemble
        strg_cmd = ["stringtie",
                    out_dir+"strg.real_pol.sample"+str(i)+"/hisat.sorted.bam",
                    "-t",
                    "-p",str(args.threads),
                    "-G",base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.true.gtf",
                    "-o",out_dir+"strg.real_pol.sample"+str(i)+"/strg.gtf"]
        subprocess.call(strg_cmd)
        # quantify with salmon
        salmon_cmd = ["salmon","quant","--validateMappings","-l","A",
                      "-i",out_dir+"res_distrib.all.real",
                      "-p",str(args.threads),
                      "-1",out_dir+"res_distrib.real_pol.sample"+str(i)+"/sample_01_1.shuffled.fasta",
                      "-2",out_dir+"res_distrib.real_pol.sample"+str(i)+"/sample_01_2.shuffled.fasta",
                      "-o",out_dir+"slmn.real_pol.sample"+str(i)]
        subprocess.call(salmon_cmd)

        # next perform the same analysis with all reads
        if not os.path.exists(out_dir+"strg.all.sample"+str(i)):
            os.makedirs(out_dir+"strg.all.sample"+str(i))
        # align
        hisat_cmd = [args.hisat,
                     "-x",out_dir+"res_distrib.all.real",
                     "-p",str(args.threads),
                     "--rna-sensitive","-f",
                     "-1",out_dir+"res_distrib.all.sample"+str(i)+"/sample_01_1.shuffled.fasta",
                     "-2",out_dir+"res_distrib.all.sample"+str(i)+"/sample_01_2.shuffled.fasta",
                     "-S",out_dir+"strg.all.sample"+str(i)+"/hisat.sam"]
        subprocess.call(hisat_cmd)
        # sort
        sort_cmd = ["samtools","sort",
                    "-@",str(args.threads),
                    "-o",out_dir+"strg.all.sample"+str(i)+"/hisat.sorted.bam",
                    out_dir+"strg.all.sample"+str(i)+"/hisat.sam"]
        subprocess.call(sort_cmd)
        # assemble
        strg_cmd = ["stringtie",
                    out_dir+"strg.all.sample"+str(i)+"/hisat.sorted.bam",
                    "-t",
                    "-p",str(args.threads),
                    "-G",base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.true.gtf",
                    "-o",out_dir+"strg.all.sample"+str(i)+"/strg.gtf"]
        subprocess.call(strg_cmd)
        # quantify with salmon
        salmon_cmd = ["salmon","quant","--validateMappings","-l","A",
                      "-i",out_dir+"res_distrib.all.real",
                      "-p",str(args.threads),
                      "-1",out_dir+"res_distrib.all.sample"+str(i)+"/sample_01_1.shuffled.fasta",
                      "-2",out_dir+"res_distrib.all.sample"+str(i)+"/sample_01_2.shuffled.fasta",
                      "-o",out_dir+"slmn.all.sample"+str(i)]
        subprocess.call(salmon_cmd)

    # todo: the problem is with the intronic loci... Those should only be inside the real loci which have been selected...
    # how do we do that???

    return

def main(args):
    parser = argparse.ArgumentParser(description='''Help Page''')

    parser.add_argument('--basedirdata',
                        required=True,
                        type=str,
                        help="base data directory")
    parser.add_argument('--basedirout',
                        required=True,
                        type=str,
                        help="base output directory")
    parser.add_argument('--outdir',
                        required=True,
                        type=str,
                        help="output directory")
    parser.add_argument('--ref',
                        required=True,
                        type=str,
                        help="reference fasta")
    parser.add_argument('--genRNAseq',
                        required=True,
                        type=str,
                        help="genRNAseq path")
    parser.add_argument('--hisat',
                        required=True,
                        type=str,
                        help="hisat2 path")
    parser.add_argument('--threads',
                        required=True,
                        type=int,
                        help="number of threads to use")
    parser.add_argument("--shuffleReads",
                        required=True,
                        type=str,
                        help="path to the shuffleReads script")
    parser.add_argument("--numSamples",
                        required=True,
                        type=int,
                        help="number of samples to simulate")
    parser.add_argument("--readlen",
                        required=True,
                        type=int,
                        help="read length")

    parser.set_defaults(func=wrapper)
    args=parser.parse_args()
    args.func(args)


if __name__=="__main__":
    main(sys.argv[1:])