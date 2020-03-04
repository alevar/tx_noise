#!/usr/bin/env python

# the difference from the regular wrapper is that here we sample TPMs for real and nonintronic conditionally, preserving relationships

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

def writeRSEM_stats(args):
    rsem_model_fp = open(outdir+"/")

def wrapper(args):
    # declarations
    base_dir_data = args.basedirdata
    base_dir_out = args.basedirout
    out_dir = os.path.abspath(args.outdir)+"/"
    hg38_fa = args.ref

    genRNAseq = args.genRNAseq

    readlen = args.readlen
    num_samples = args.numSamples # number of samples to simulate

    gff3cols=["seqid","source","type","start","end","score","strand","phase","attributes"]
    rsem_tpm_cols = ["transcript_id","gene_id","length","effective_length","expected_count","TPM","FPKM","IsoPct"]
    extension = ".fasta"
    if args.type == "rsem":
        extension = ".fa"
    # Stats holders
    stats = {
        "starting number of real and splicing noise loci":0,
        "starting number of splicing noise loci":0,
        "starting number of intronic noise loci":0,
        "starting number of polymerase noise loci": 0,
        "adjusted starting number of real only loci":0,
        "adjusted starting number of splicing noise loci":0,
        "adjusted starting number of intronic noise loci":0,
        "adjusted starting number of polymerase noise loci": 0,
        "average number of real loci":0,
        "average number of noise loci":0,
        "average number of undefined loci":0,
        "selecting n real locs":0,
        "selecting n splicing noise locs":0,
        "selecting n intronic locs":0,
        "selecting n polymerase locs":0,
    }

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # load base annotations
    print(">>>loading base annotations")
    real_baseDF = pd.read_csv(base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.true.gtf",sep="\t",names=gff3cols)
    real_baseDF = real_baseDF[(real_baseDF["type"]=="transcript") & (real_baseDF["strand"].isin(["+","-"]))].reset_index(drop=True)
    nonint_baseDF = pd.read_csv(base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.no.contained.non.intronic.gtf",sep="\t",names=gff3cols)
    nonint_baseDF = nonint_baseDF[(nonint_baseDF["type"]=="transcript") & (nonint_baseDF["strand"].isin(["+","-"]))].reset_index(drop=True)
    int_baseDF = pd.read_csv(base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.no.contained.intronic.gtf",sep="\t",names=gff3cols)
    int_baseDF = int_baseDF[(int_baseDF["type"]=="transcript") & (int_baseDF["strand"].isin(["+","-"]))].reset_index(drop=True)
    pol_baseDF = pd.read_csv(base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.no.contained.RNApol.gtf",sep="\t",names=gff3cols)
    pol_baseDF = pol_baseDF[(pol_baseDF["type"]=="transcript") & (pol_baseDF["strand"].isin(["+","-"]))].reset_index(drop=True)

    # get all loci and transcript IDs
    print(">>>getting loci IDs")
    real_baseDF["lid"] = real_baseDF["attributes"].str.split("gene_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    nonint_baseDF["lid"] = nonint_baseDF["attributes"].str.split("gene_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    int_baseDF["lid"] = int_baseDF["attributes"].str.split("gene_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    pol_baseDF["lid"] = pol_baseDF["attributes"].str.split("gene_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    real_locs = set(real_baseDF["lid"])
    nonint_locs = set(nonint_baseDF["lid"])
    int_locs = set(int_baseDF["lid"])
    pol_locs = set(pol_baseDF["lid"])
    stats["starting number of real and splicing noise loci"] = len(real_locs)
    stats["starting number of splicing noise loci"] = len(nonint_locs)
    stats["starting number of intronic loci"] = len(int_locs)
    stats["starting number of polymerase loci"] = len(pol_locs)

    total_n_locs = len(real_locs)+len(nonint_locs)+len(int_locs)+len(pol_locs)

    # perform cleanup, by removing any loci in 
    #   1. int that are also in real
    #   2. nonint that are not in real
    #   3. pol that are in real
    int_locs = int_locs - real_locs.intersection(int_locs)
    assert(len(real_locs.intersection(int_locs))==0),"something wrong intronic"
    int_baseDF = int_baseDF[int_baseDF["lid"].isin(int_locs)].reset_index(drop=True)

    nonint_locs = nonint_locs - nonint_locs.difference(real_locs)
    assert(len(nonint_locs) == len(real_locs.intersection(nonint_locs))),"something wrong non-intronic"
    nonint_baseDF = nonint_baseDF[nonint_baseDF["lid"].isin(nonint_locs)].reset_index(drop=True)

    pol_locs = pol_locs - real_locs.intersection(pol_locs)
    assert(len(real_locs.intersection(pol_locs))==0),"something wrong polymerase"
    pol_baseDF = pol_baseDF[pol_baseDF["lid"].isin(pol_locs)].reset_index(drop=True)

    real_baseDF["tid"] = real_baseDF["attributes"].str.split("transcript_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    real_baseDF = real_baseDF[["lid","tid"]]
    nonint_baseDF["tid"] = nonint_baseDF["attributes"].str.split("transcript_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    nonint_baseDF = nonint_baseDF[["lid","tid"]]
    int_baseDF["tid"] = int_baseDF["attributes"].str.split("transcript_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    int_baseDF = int_baseDF[["lid","tid"]]
    pol_baseDF["tid"] = pol_baseDF["attributes"].str.split("transcript_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    pol_baseDF = pol_baseDF[["lid","tid"]]

    # get the averages
    avg_real_locs = 0
    avg_noise_locs = 0
    avg_undefined_locs = 0
    with open(base_dir_out+"res_distrib.loc_stats","r") as inFP:
        for line in inFP.readlines():
            type_loc,num_loc = line.split(":")
            if type_loc=="real":
                avg_real_locs = int(num_loc)
            elif type_loc=="noise":
                avg_noise_locs = int(num_loc)
            elif type_loc=="undefined":
                if(int(num_loc))>0:
                    avg_undefined_locs = int(num_loc)
            else:
                print("error in loading averages")
    stats["average number of real loci"] = avg_real_locs
    stats["average number of noise loci"] = avg_noise_locs
    stats["average number of undefined loci"] = avg_undefined_locs

    # randomly select loci of each type based on the averages
    # noise loci are a combination of: intronic and polymerase
    # this needs to be done proportionally to the consistent parts of each group (real(real,nonint) and noise(int,pol))

    real_only_locs = real_locs - nonint_locs

    # in case of inconsistencies associated with merging transcripts of different types under the same locus
    # we need to scale the averages with respect to the new numbers
    avg_total_locs = avg_real_locs+avg_noise_locs
    sum_locs_to_choose_from = len(real_only_locs)+len(nonint_locs)+len(int_locs)+len(pol_locs)

    stats["adjusted number of real only loci"] = len(real_only_locs)
    stats["adjusted number of splicing noise loci"] = len(nonint_locs)
    stats["adjusted number of intronic loci"] = len(int_locs)
    stats["adjusted number of polymerase loci"] = len(pol_locs)

    perc_real = len(real_only_locs)/len(real_locs)
    perc_nonint = len(nonint_locs)/len(real_locs)
    assert(perc_real+perc_nonint)==1,"wrong percent real and nonint"

    perc_int = len(int_locs)/(len(int_locs)+len(pol_locs))
    perc_pol = len(pol_locs)/(len(int_locs)+len(pol_locs))
    assert((perc_int+perc_pol)==1),"wrong percent pol and int"
    stats["selecting n real locs"] = int(avg_real_locs*perc_real)
    stats["selecting n splicing noise locs"] = int(avg_real_locs*perc_nonint)
    stats["selecting n intronic locs"] = int(avg_noise_locs*perc_int)
    stats["selecting n polymerase locs"] = int(avg_noise_locs*perc_pol)

    real_only_locs_rand = np.random.choice(list(real_only_locs),int(avg_real_locs*perc_real), replace=False)
    nonint_locs_rand = np.random.choice(list(nonint_locs),int(avg_real_locs*perc_nonint), replace=False)
    int_locs_rand = np.random.choice(list(int_locs),int(avg_noise_locs*perc_int), replace=False)
    pol_locs_rand = np.random.choice(list(pol_locs),int(avg_noise_locs*perc_pol), replace=False)

    real_only_baseDF = real_baseDF[real_baseDF["lid"].isin(real_only_locs_rand)].reset_index(drop=True)
    real_nonint_baseDF = real_baseDF[(real_baseDF["lid"].isin(nonint_locs_rand))&~(real_baseDF["lid"].isin(real_only_locs_rand))].reset_index(drop=True)
    nonint_baseDF = nonint_baseDF[nonint_baseDF["lid"].isin(nonint_locs_rand)].reset_index(drop=True)
    int_baseDF = int_baseDF[int_baseDF["lid"].isin(int_locs_rand)].reset_index(drop=True)
    pol_baseDF = pol_baseDF[pol_baseDF["lid"].isin(pol_locs_rand)].reset_index(drop=True)

    # load the transcript tissue distribution
    ttx_real = list()
    ttx_nonint = list()
    ttx_int = list()
    ttx_pol = list()

    tmp_tx = pd.read_csv(base_dir_out+"res_distrib.num_tx_per_tissue_loc2")
    tmp_tx_real_loc = tmp_tx[(tmp_tx["real"]>0)][["total","real","nonintronic"]].reset_index(drop=True)
    tmp_tx_real_loc["perc_real"] = tmp_tx_real_loc["real"]/tmp_tx_real_loc["total"]
    tmp_tx_real_loc["perc_nonint"] = tmp_tx_real_loc["nonintronic"]/tmp_tx_real_loc["total"]
    d_real = list(tmp_tx_real_loc["perc_real"])
    d_nonint = list(tmp_tx_real_loc["perc_nonint"])

    tmp_tx_int_loc = tmp_tx[(tmp_tx["intronic"]>0)][["total","intronic"]].reset_index(drop=True)
    tmp_tx_int_loc["perc_int"] = tmp_tx_int_loc["intronic"]/tmp_tx_int_loc["total"]
    d_int = list(tmp_tx_int_loc["perc_int"])

    tmp_tx_pol_loc = tmp_tx[(tmp_tx["polymerase"]>0)][["total","polymerase"]].reset_index(drop=True)
    tmp_tx_pol_loc["perc_pol"] = tmp_tx_pol_loc["polymerase"]/tmp_tx_pol_loc["total"]
    d_pol = list(tmp_tx_pol_loc["perc_pol"])

    # now that we have the loci, we can subset the gtfs and move on to working with the transcripts
    # the question is: how do we allocate transcripts per locus?
    # do we do this based on the distribution of the number of transcripts per each type of locus?
    # are there any caveates in the co-dependencies such as between real and non-intronic?

    # we should probably only select transcripts of each group within each type of locus dataframe

    # one important consideration is that the total number of transcripts should be close to what is observed
    #  in an average tissue (needs to be computed in the gtex_stats)
    # this approach, however, might be overly complicated, and simply randomly selecting n transcripts for each group
    #  based on distribution computed by gtex_stats might be sufficient, and should approximate well to desired values

    real_only_txs = set()
    for name, group in real_only_baseDF.groupby("lid"):
        tids = group["tid"].tolist()
        nt = math.ceil(d_real[random.randint(0,len(d_real)-1)]*len(tids))
        random.shuffle(tids)
        real_only_txs.update(tids[:nt])

    real_nonint_txs = set()
    for name, group in real_nonint_baseDF.groupby("lid"):
        tids = group["tid"].tolist()
        nt = math.ceil(d_real[random.randint(0,len(d_real)-1)]*len(tids))
        random.shuffle(tids)
        real_nonint_txs.update(tids[:nt])
        
    nonint_txs = set()
    for name, group in nonint_baseDF.groupby("lid"):
        tids = group["tid"].tolist()
        nt = math.ceil(d_nonint[random.randint(0,len(d_nonint)-1)]*len(tids))
        random.shuffle(tids)
        nonint_txs.update(tids[:nt])
        
    int_txs = set()
    for name, group in int_baseDF.groupby("lid"):
        tids = group["tid"].tolist()
        nt = math.ceil(d_int[random.randint(0,len(d_int)-1)]*len(tids))
        random.shuffle(tids)
        int_txs.update(tids[:nt])
        
    pol_txs = set()
    for name, group in pol_baseDF.groupby("lid"):
        tids = group["tid"].tolist()
        nt = math.ceil(d_pol[random.randint(0,len(d_pol)-1)]*len(tids))
        random.shuffle(tids)
        pol_txs.update(tids[:nt])

    stats["number of real only txs"] = len(real_only_txs)
    stats["number of txs in real loci with splicing noise"] = len(real_nonint_txs)
    stats["number of nonint txs"] = len(nonint_txs)
    stats["number of int txs"] = len(int_txs)
    stats["number of pol txs"] = len(pol_txs)

    # now that we have lists of transcripts - we can subset the annotation further
    real_only_baseDF = real_only_baseDF[real_only_baseDF["tid"].isin(real_only_txs)].reset_index(drop=True)
    real_nonint_baseDF = real_nonint_baseDF[real_nonint_baseDF["tid"].isin(real_nonint_txs)].reset_index(drop=True)
    nonint_baseDF = nonint_baseDF[nonint_baseDF["tid"].isin(nonint_txs)].reset_index(drop=True)
    int_baseDF = int_baseDF[int_baseDF["tid"].isin(int_txs)].reset_index(drop=True)
    pol_baseDF = pol_baseDF[pol_baseDF["tid"].isin(pol_txs)].reset_index(drop=True)

    # lastly we need to generate samples based on the underlying annotation
    # for this we will need additional information, similar to the tissue_per_loc2 in order to tell
    # how likely a sample is to contribute a transcript to a tissue

    # load the transcript sample distribution
    ttx_real = list()
    ttx_nonint = list()
    ttx_int = list()
    ttx_pol = list()

    tmp_tx = pd.read_csv(base_dir_out+"res_distrib.num_tx_per_sample_loc2")
    tmp_tx_real_loc = tmp_tx[(tmp_tx["real"]>0)][["total","real","nonintronic"]].reset_index(drop=True)
    tmp_tx_real_loc["perc_real"] = tmp_tx_real_loc["real"]/tmp_tx_real_loc["total"]
    tmp_tx_real_loc["perc_nonint"] = tmp_tx_real_loc["nonintronic"]/tmp_tx_real_loc["total"]
    d_real = list(tmp_tx_real_loc["perc_real"])
    d_nonint = list(tmp_tx_real_loc["perc_nonint"])

    tmp_tx_int_loc = tmp_tx[(tmp_tx["intronic"]>0)][["total","intronic"]].reset_index(drop=True)
    tmp_tx_int_loc["perc_int"] = tmp_tx_int_loc["intronic"]/tmp_tx_int_loc["total"]
    d_int = list(tmp_tx_int_loc["perc_int"])

    tmp_tx_pol_loc = tmp_tx[(tmp_tx["polymerase"]>0)][["total","polymerase"]].reset_index(drop=True)
    tmp_tx_pol_loc["perc_pol"] = tmp_tx_pol_loc["polymerase"]/tmp_tx_pol_loc["total"]
    d_pol = list(tmp_tx_pol_loc["perc_pol"])

    real_only_txs = [set() for x in range(num_samples)]
    for name, group in real_only_baseDF.groupby("lid"):
        tids = group["tid"].tolist()
        for i in range(num_samples):
            nt = math.ceil(d_real[random.randint(0,len(d_real)-1)]*len(tids))
            random.shuffle(tids)
            real_only_txs[i].update(tids[:nt])

    real_nonint_txs = [set() for x in range(num_samples)]
    for name, group in real_nonint_baseDF.groupby("lid"):
        tids = group["tid"].tolist()
        for i in range(num_samples):
            nt = math.ceil(d_real[random.randint(0,len(d_real)-1)]*len(tids))
            random.shuffle(tids)
            real_nonint_txs[i].update(tids[:nt])

    nonint_txs = [set() for x in range(num_samples)]
    for name, group in nonint_baseDF.groupby("lid"):
        tids = group["tid"].tolist()
        for i in range(num_samples):
            nt = math.ceil(d_nonint[random.randint(0,len(d_nonint)-1)]*len(tids))
            random.shuffle(tids)
            nonint_txs[i].update(tids[:nt])

    int_txs = [set() for x in range(num_samples)]
    for name, group in int_baseDF.groupby("lid"):
        tids = group["tid"].tolist()
        for i in range(num_samples):
            nt = math.ceil(d_int[random.randint(0,len(d_int)-1)]*len(tids))
            random.shuffle(tids)
            int_txs[i].update(tids[:nt])
            
    pol_txs = [set() for x in range(num_samples)]
    for name, group in pol_baseDF.groupby("lid"):
        tids = group["tid"].tolist()
        for i in range(num_samples):
            nt = math.ceil(d_pol[random.randint(0,len(d_pol)-1)]*len(tids))
            random.shuffle(tids)
            pol_txs[i].update(tids[:nt])

    print(">>> Generating subsets for Intronic and Polymerase")
    # now that we have sample specific stuff - we can generate sample-specific GTFs
    int_baseDF = pd.read_csv(base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.no.contained.intronic.gtf",sep="\t",names=gff3cols)
    pol_baseDF = pd.read_csv(base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.no.contained.RNApol.gtf",sep="\t",names=gff3cols)
    int_baseDF["tid"] = int_baseDF["attributes"].str.split("transcript_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    pol_baseDF["tid"] = pol_baseDF["attributes"].str.split("transcript_id \"",expand=True)[1].str.split("\"",expand=True)[0]

        # load joint expressions
    print(">>> Getting subsets for Real and Nonint")
    real_baseDF = pd.read_csv(base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.true.gtf",sep="\t",names=gff3cols)
    nonint_baseDF = pd.read_csv(base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.no.contained.non.intronic.gtf",sep="\t",names=gff3cols)
    real_baseDF["tid"] = real_baseDF["attributes"].str.split("transcript_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    nonint_baseDF["tid"] = nonint_baseDF["attributes"].str.split("transcript_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    # load joint expressions
    joined_df = pd.read_csv(base_dir_out+"res_distrib.covs_joined")
    joined_df["code"] = joined_df["nt_real"].astype(int).astype(str)+"-"+joined_df["nt_nonint"].astype(int).astype(str)
    joined_df.drop(["nt_real","nt_nonint"],axis=1,inplace=True)

    for i in range(num_samples):
        int_baseDF[int_baseDF["tid"].isin(int_txs[i])][gff3cols].to_csv(out_dir+"/res_distrib.int.sample"+str(i)+".gtf",sep="\t",index=False,header=False,quoting=csv.QUOTE_NONE)
        pol_baseDF[pol_baseDF["tid"].isin(pol_txs[i])][gff3cols].to_csv(out_dir+"/res_distrib.pol.sample"+str(i)+".gtf",sep="\t",index=False,header=False,quoting=csv.QUOTE_NONE)

    # load distribution of the sum of expressions
    tmp_exp_int = pd.read_csv(base_dir_out+"res_distrib.cov_sample_intronic")
    d_exp_int = list(tmp_exp_int["cov"]+1)

    tmp_exp_pol = pd.read_csv(base_dir_out+"res_distrib.cov_sample_polymerase")
    d_exp_pol = list(tmp_exp_pol["cov"]+1)

    num_reads_real = []
    num_reads_nonint = []
    num_reads_int = []
    num_reads_pol = []

    tpm_contrib_real = []
    tpm_contrib_nonint = []
    tpm_contrib_int = []
    tpm_contrib_pol = []

    # now generate corresponding expression matrices for the transcripts in each sample for each type
    for i in range(num_samples):        
        intDF = pd.read_csv(out_dir+"/res_distrib.int.sample"+str(i)+".gtf",sep="\t",names=gff3cols)
        nt_int = len(intDF[intDF["type"]=="transcript"])
        exps_int = np.random.choice(d_exp_int, nt_int, replace=False)
        np.savetxt(out_dir+"/res_distrib.int.sample"+str(i)+".exp", exps_int, delimiter=',',fmt='%0.2f')

        polDF = pd.read_csv(out_dir+"/res_distrib.pol.sample"+str(i)+".gtf",sep="\t",names=gff3cols)
        nt_pol = len(polDF[polDF["type"]=="transcript"])
        exps_pol = np.random.choice(d_exp_pol, nt_pol, replace=False)
        np.savetxt(out_dir+"/res_distrib.pol.sample"+str(i)+".exp", exps_pol, delimiter=',',fmt='%0.2f')

        # need to generate data for the RSEM simulation
        tdf = int_baseDF[int_baseDF["tid"].isin(int_txs[i])]
        edf = tdf[tdf["type"]=="exon"].reset_index(drop=True)
        tdf = tdf[tdf["type"]=="transcript"][["tid","start","end"]].reset_index(drop=True) # intended for order
        tdf["len"] = tdf["end"]-tdf["start"]
        tdf.drop(["start","end"],axis=1,inplace=True)
        edf["elen"] = edf["end"]-edf["start"]
        edf = edf[["tid","elen"]]
        edf = edf.groupby("tid").agg({"elen":"sum"}).reset_index()
        # now organize by
        tdf = tdf.merge(edf,on="tid",how="left",indicator=True)
        assert len(tdf) == len(tdf[tdf["_merge"]=="both"]),"missing/extra exons"

        tdf["cov"] = exps_int
        tdf["cov"] = tdf["cov"].fillna(0)
        assert len(tdf[tdf["cov"]==0])==0,"missing TPMs"
        tdf["nr"] = (tdf["elen"]*tdf["cov"])/readlen
        num_reads_int.append(int(tdf["nr"].sum()))
        tdf["kelen"] = tdf["elen"]/1000
        tdf["rpk"] = tdf["nr"]/tdf["kelen"]
        tdf_int = tdf.copy(deep=True)
        # pm_sf = tdf["rpk"].sum()/1000000
        # tdf["tpm"] = tdf["rpk"]/pm_sf
        # tdf.drop(["_merge","nr","kelen","rpk","cov"],axis=1,inplace=True)
        # tdf.columns = ["transcript_id","length","effective_length","TPM"]
        # tdf["gene_id"] = tdf["transcript_id"]
        # tdf["expected_count"] = 0.00
        # tdf["FPKM"] = 0.00
        # tdf["IsoPct"] = 0.00
        # tdf[rsem_tpm_cols].to_csv(out_dir+"/res_distrib.int.sample"+str(i)+".rsem_tpms.tmp",sep="\t",index=False)
        # tpm_contrib_int.append(tdf["TPM"].sum())

        # need to generate data for the RSEM simulation
        tdf = pol_baseDF[pol_baseDF["tid"].isin(pol_txs[i])]
        edf = tdf[tdf["type"]=="exon"].reset_index(drop=True)
        tdf = tdf[tdf["type"]=="transcript"][["tid","start","end"]].reset_index(drop=True) # intended for order
        tdf["len"] = tdf["end"]-tdf["start"]
        tdf.drop(["start","end"],axis=1,inplace=True)
        edf["elen"] = edf["end"]-edf["start"]
        edf = edf[["tid","elen"]]
        edf = edf.groupby("tid").agg({"elen":"sum"}).reset_index()
        # now organize by
        tdf = tdf.merge(edf,on="tid",how="left",indicator=True)
        assert len(tdf) == len(tdf[tdf["_merge"]=="both"]),"missing/extra exons"

        tdf["cov"] = exps_pol
        tdf["cov"] = tdf["cov"].fillna(0)
        assert len(tdf[tdf["cov"]==0])==0,"missing coverages"
        tdf["nr"] = (tdf["elen"]*tdf["cov"])/readlen
        num_reads_pol.append(int(tdf["nr"].sum()))
        tdf["kelen"] = tdf["elen"]/1000
        tdf["rpk"] = tdf["nr"]/tdf["kelen"]
        tdf_pol = tdf.copy(deep=True)
        # pm_sf = tdf["rpk"].sum()/1000000
        # tdf["tpm"] = tdf["rpk"]/pm_sf
        # tdf.drop(["_merge","nr","kelen","rpk","cov"],axis=1,inplace=True)
        # tdf.columns = ["transcript_id","length","effective_length","TPM"]
        # tdf["gene_id"] = tdf["transcript_id"]
        # tdf["expected_count"] = 0.00
        # tdf["FPKM"] = 0.00
        # tdf["IsoPct"] = 0.00
        # tdf[rsem_tpm_cols].to_csv(out_dir+"/res_distrib.pol.sample"+str(i)+".rsem_tpms.tmp",sep="\t",index=False)
        # tpm_contrib_pol.append(tdf["TPM"].sum())

        # now need to somehow get the tpms for the samples
        # while intronic and polymerase are easy to sample
        # real and nonintronic are correlated
        # so we need to use the special distribution to sample

        # first need to be able to process real and nonintronic together
        # and output expression accordingly
        sub_real = real_baseDF[real_baseDF["tid"].isin(real_only_txs[i].union(real_nonint_txs[i]))].reset_index(drop=True)
        sub_real["gid"] = sub_real["attributes"].str.split("gene_id \"",expand=True)[1].str.split("\"",expand=True)[0]
        sub_nonint = nonint_baseDF[nonint_baseDF["tid"].isin(nonint_txs[i])].reset_index(drop=True)
        sub_nonint["gid"] = sub_nonint["attributes"].str.split("gene_id \"",expand=True)[1].str.split("\"",expand=True)[0]

        gsr = sub_real[sub_real["type"]=="transcript"][["gid","tid"]].groupby("gid").count().reset_index()
        gsr.columns = ["gid","ntr"]
        gsn = sub_nonint[sub_nonint["type"]=="transcript"][["gid","tid"]].groupby("gid").count().reset_index()
        gsn.columns = ["gid","ntn"]
        gs = gsr.merge(gsn,on="gid",how="outer",indicator=True)
        print("removing "+str(len(gs[gs["_merge"]=="right_only"]))+" number of loci due to presence of nonint only transcripts")
        gs = gs[~(gs["_merge"]=="right_only")].reset_index(drop=True)
        gs["ntn"] = gs["ntn"].fillna(0)
        gs["ntr"] = gs["ntr"].fillna(0)
        gs["code"] = gs["ntr"].astype(int).astype(str)+"-"+gs["ntn"].astype(int).astype(str)
        gs.sort_values(by="code",inplace=True)

        real_exp = open(out_dir+"/res_distrib.real.sample"+str(i)+".exp","w+")
        nonint_exp = open(out_dir+"/res_distrib.nonint.sample"+str(i)+".exp","w+")

        jd = joined_df.copy(deep=True)
        # compute the number of each code in gs
        tmp = gs.groupby(by = 'code')['gid'].count().reset_index()
        tmp.columns = ["code","nc"]
        jd = jd.merge(tmp,on="code",how="inner",indicator=True)
        # need to have a way of knowing which loci have no code matches...
        # and to maintain the order between the expressions and corresponding gtfs
        jd = jd[jd["_merge"]=="both"]
        gs = gs[gs["code"].isin(set(jd["code"]))]
        # now we can select th
        jd = jd.groupby('code', as_index=False).apply(lambda obj: obj.loc[np.random.choice(obj.index, obj.nc.iloc[0],replace=True),:]).reset_index(drop=True)

        gs.sort_values(by="code",inplace=True,ascending=False)
        gs.reset_index(drop=True,inplace=True)
        jd.sort_values(by="code",inplace=True,ascending=False)
        jd.reset_index(drop=True,inplace=True)
        gs[["real","nonint","code_jd"]] = jd[["real","nonint","code"]]
        assert len(gs)==len(gs[gs["code"]==gs["code_jd"]]),"codes do not match"

        # now we can split the expressions
        rgs = gs[gs["ntr"]>0][["gid","real"]].reset_index(drop=True)
        rgs["reals"] = rgs["real"].str.split(";")
        rgs = rgs.explode("reals")
        rgs.drop("real",axis=1,inplace=True)
        # we can save the expressions, however, now we need to make sure the gtf
        # is saved in the same order
        sub_real = sub_real[sub_real["gid"].isin(rgs["gid"])].reset_index(drop=True)
        assert len(set(rgs["gid"]).intersection(set(sub_real["gid"])))==len(set(sub_real["gid"])),"incompatible gids in sub_real and gs"
        # to finally ensure the same order between the two dataframes we can
        # get the order of genes in one dataframe
        # and then merge into it independently both expression and gtf data
        #   to guarantee identical ordering
        gene_order = rgs[["gid"]].reset_index(drop=True)
        gene_order.drop_duplicates(keep="first",inplace=True)
        gene_order.reset_index(drop=True,inplace=True)
        real_gtf = gene_order.merge(sub_real,on="gid",how="inner")
        real_exp = gene_order.merge(rgs,on="gid",how="inner")
        real_exp.reset_index(drop=True,inplace=True)
        real_exp["reals"] = real_exp["reals"].astype(float).round(2)+1.0
        real_gtf[gff3cols].to_csv(out_dir+"/res_distrib.real.sample"+str(i)+".gtf",sep="\t",index=False,header=False,quoting=csv.QUOTE_NONE)
        real_exp[["reals"]].to_csv(out_dir+"/res_distrib.real.sample"+str(i)+".exp",index=False,header=False)
        
        # now we can split the expressions for the nonint
        ngs = gs[gs["ntn"]>0][["gid","nonint"]].reset_index(drop=True)
        ngs["nonints"] = ngs["nonint"].str.split(";")
        ngs = ngs.explode("nonints")
        ngs.drop("nonint",axis=1,inplace=True)
        # we can save the expressions, however, now we need to make sure the gtf
        # is saved in the same order
        sub_nonint = sub_nonint[sub_nonint["gid"].isin(ngs["gid"])].reset_index(drop=True)
        assert len(set(ngs["gid"]).intersection(set(sub_nonint["gid"])))==len(set(sub_nonint["gid"])),"incompatible gids in sub_nonint and gs"
        # to finally ensure the same order between the two dataframes we can
        # get the order of genes in one dataframe
        # and then merge into it independently both expression and gtf data
        #   to guarantee identical ordering
        gene_order = ngs[["gid"]].reset_index(drop=True)
        gene_order.drop_duplicates(keep="first",inplace=True)
        gene_order.reset_index(drop=True,inplace=True)
        nonint_gtf = gene_order.merge(sub_nonint,on="gid",how="inner")
        nonint_exp = gene_order.merge(ngs,on="gid",how="inner")
        nonint_exp.reset_index(drop=True,inplace=True)
        nonint_exp["nonints"] = nonint_exp["nonints"].astype(float).round(2)+1.0
        nonint_gtf[gff3cols].to_csv(out_dir+"/res_distrib.nonint.sample"+str(i)+".gtf",sep="\t",index=False,header=False,quoting=csv.QUOTE_NONE)
        nonint_exp[["nonints"]].to_csv(out_dir+"/res_distrib.nonint.sample"+str(i)+".exp",index=False,header=False)
        
        # process for RSEM
        edf = real_gtf[real_gtf["type"]=="exon"].reset_index(drop=True)
        tdf = real_gtf[real_gtf["type"]=="transcript"][["tid","start","end"]].reset_index(drop=True) # intended for order
        tdf["len"] = tdf["end"]-tdf["start"]
        tdf.drop(["start","end"],axis=1,inplace=True)
        edf["elen"] = edf["end"]-edf["start"]
        edf = edf[["tid","elen"]]
        edf = edf.groupby("tid").agg({"elen":"sum"}).reset_index()
        # now organize by
        tdf = tdf.merge(edf,on="tid",how="left",indicator=True)
        assert len(tdf) == len(tdf[tdf["_merge"]=="both"]),"missing/extra exons"

        tdf["cov"] = real_exp["reals"]
        tdf["cov"] = tdf["cov"].fillna(0)
        assert len(tdf[tdf["cov"]==0])==0,"missing TPMs"
        tdf["nr"] = (tdf["elen"]*tdf["cov"])/readlen
        num_reads_real.append(int(tdf["nr"].sum()))
        tdf["kelen"] = tdf["elen"]/1000
        tdf["rpk"] = tdf["nr"]/tdf["kelen"]
        tdf_real = tdf.copy(deep=True)
        # pm_sf = tdf["rpk"].sum()/1000000
        # tdf["tpm"] = tdf["rpk"]/pm_sf
        # tdf.drop(["_merge","nr","kelen","rpk","cov"],axis=1,inplace=True)
        # tdf.columns = ["transcript_id","length","effective_length","TPM"]
        # tdf["gene_id"] = tdf["transcript_id"]
        # tdf["expected_count"] = 0.00
        # tdf["FPKM"] = 0.00
        # tdf["IsoPct"] = 0.00
        # tdf[rsem_tpm_cols].to_csv(out_dir+"/res_distrib.real.sample"+str(i)+".rsem_tpms.tmp",sep="\t",index=False)
        # tpm_contrib_real.append(tdf["TPM"].sum())

        edf = nonint_gtf[nonint_gtf["type"]=="exon"].reset_index(drop=True)
        tdf = nonint_gtf[nonint_gtf["type"]=="transcript"][["tid","start","end"]].reset_index(drop=True) # intended for order
        tdf["len"] = tdf["end"]-tdf["start"]
        tdf.drop(["start","end"],axis=1,inplace=True)
        edf["elen"] = edf["end"]-edf["start"]
        edf = edf[["tid","elen"]]
        edf = edf.groupby("tid").agg({"elen":"sum"}).reset_index()
        # now organize by
        tdf = tdf.merge(edf,on="tid",how="left",indicator=True)
        assert len(tdf) == len(tdf[tdf["_merge"]=="both"]),"missing/extra exons"

        tdf["cov"] = nonint_exp["nonints"]
        tdf["cov"] = tdf["cov"].fillna(0)
        assert len(tdf[tdf["cov"]==0])==0,"missing TPMs"
        tdf["nr"] = (tdf["elen"]*tdf["cov"])/readlen
        num_reads_nonint.append(int(tdf["nr"].sum()))
        tdf["kelen"] = tdf["elen"]/1000
        tdf["rpk"] = tdf["nr"]/tdf["kelen"]
        tdf_nonint = tdf.copy(deep=True)
        # pm_sf = tdf["rpk"].sum()/1000000
        # tdf["tpm"] = tdf["rpk"]/pm_sf
        # tdf.drop(["_merge","nr","kelen","rpk","cov"],axis=1,inplace=True)
        # tdf.columns = ["transcript_id","length","effective_length","TPM"]
        # tdf["gene_id"] = tdf["transcript_id"]
        # tdf["expected_count"] = 0.00
        # tdf["FPKM"] = 0.00
        # tdf["IsoPct"] = 0.00
        # tdf[rsem_tpm_cols].to_csv(out_dir+"/res_distrib.nonint.sample"+str(i)+".rsem_tpms.tmp",sep="\t",index=False)
        # tpm_contrib_nonint.append(tdf["TPM"].sum())

        # lastly need to process TPMs jointly here to get the relative abundances of each
        tdf = pd.concat([tdf_real,tdf_nonint,tdf_int,tdf_pol])
        tdf.reset_index(drop=True)
        pm_sf = tdf["rpk"].sum()/1000000
        
        tdf_real["tpm"] = tdf_real["rpk"]/pm_sf
        tdf_real.drop(["_merge","nr","kelen","rpk","cov"],axis=1,inplace=True)
        tdf_real.columns = ["transcript_id","length","effective_length","TPM"]
        tdf_real["gene_id"] = tdf_real["transcript_id"]
        tdf_real["expected_count"] = 0.00
        tdf_real["FPKM"] = 0.00
        tdf_real["IsoPct"] = 0.00
        tdf_real[rsem_tpm_cols].to_csv(out_dir+"/res_distrib.real.sample"+str(i)+".rsem_tpms.tmp",sep="\t",index=False)
        tpm_contrib_real.append(tdf_real["TPM"].sum())

        tdf_nonint["tpm"] = tdf_nonint["rpk"]/pm_sf
        tdf_nonint.drop(["_merge","nr","kelen","rpk","cov"],axis=1,inplace=True)
        tdf_nonint.columns = ["transcript_id","length","effective_length","TPM"]
        tdf_nonint["gene_id"] = tdf_nonint["transcript_id"]
        tdf_nonint["expected_count"] = 0.00
        tdf_nonint["FPKM"] = 0.00
        tdf_nonint["IsoPct"] = 0.00
        tdf_nonint[rsem_tpm_cols].to_csv(out_dir+"/res_distrib.nonint.sample"+str(i)+".rsem_tpms.tmp",sep="\t",index=False)
        tpm_contrib_nonint.append(tdf_nonint["TPM"].sum())

        tdf_int["tpm"] = tdf_int["rpk"]/pm_sf
        tdf_int.drop(["_merge","nr","kelen","rpk","cov"],axis=1,inplace=True)
        tdf_int.columns = ["transcript_id","length","effective_length","TPM"]
        tdf_int["gene_id"] = tdf_int["transcript_id"]
        tdf_int["expected_count"] = 0.00
        tdf_int["FPKM"] = 0.00
        tdf_int["IsoPct"] = 0.00
        tdf_int[rsem_tpm_cols].to_csv(out_dir+"/res_distrib.int.sample"+str(i)+".rsem_tpms.tmp",sep="\t",index=False)
        tpm_contrib_int.append(tdf_int["TPM"].sum())

        tdf_pol["tpm"] = tdf_pol["rpk"]/pm_sf
        tdf_pol.drop(["_merge","nr","kelen","rpk","cov"],axis=1,inplace=True)
        tdf_pol.columns = ["transcript_id","length","effective_length","TPM"]
        tdf_pol["gene_id"] = tdf_pol["transcript_id"]
        tdf_pol["expected_count"] = 0.00
        tdf_pol["FPKM"] = 0.00
        tdf_pol["IsoPct"] = 0.00
        tdf_pol[rsem_tpm_cols].to_csv(out_dir+"/res_distrib.pol.sample"+str(i)+".rsem_tpms.tmp",sep="\t",index=False)
        tpm_contrib_pol.append(tdf_pol["TPM"].sum())

    # Right now we need to also have information about the proportion that each type of transcript contributes to total expression in this dataset
    # these numbers can then be used to allocate the number of reads accordingly in or data


    # save simulation statistics
    stats_fp = open(out_dir+"res_distrib.sim_stats","w+")
    for k,v in stats.items():
        stats_fp.write(k+","+str(v)+"\n")

    stats_fp.close()

    print(">>> Extracting transcripts with gffread")
    for i in range(num_samples):
        cmd_real = ["gffread",
                    "-w",out_dir+"/res_distrib.real.sample"+str(i)+".fasta",
                    "-g",hg38_fa,
                    out_dir+"/res_distrib.real.sample"+str(i)+".gtf"]
        subprocess.call(cmd_real)
        
        cmd_nonint = ["gffread",
                    "-w",out_dir+"/res_distrib.nonint.sample"+str(i)+".fasta",
                    "-g",hg38_fa,
                    out_dir+"/res_distrib.nonint.sample"+str(i)+".gtf"]
        subprocess.call(cmd_nonint)
        
        cmd_int = ["gffread",
                    "-w",out_dir+"/res_distrib.int.sample"+str(i)+".fasta",
                    "-g",hg38_fa,
                    out_dir+"/res_distrib.int.sample"+str(i)+".gtf"]
        subprocess.call(cmd_int)
        
        cmd_pol = ["gffread",
                    "-w",out_dir+"/res_distrib.pol.sample"+str(i)+".fasta",
                    "-g",hg38_fa,
                    out_dir+"/res_distrib.pol.sample"+str(i)+".gtf"]
        subprocess.call(cmd_pol)

    # also extract fasta from the set of all real transcripts for the salmon/kallisto indices
    cmd_all = ["gffread",
               "-w",out_dir+"/res_distrib.all.real.fasta",
               "-g",hg38_fa,
               base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.true.gtf"]
    subprocess.call(cmd_all)

    print(">>> Building reference with RSEM")
    for i in range(num_samples):
        if not os.path.exists(out_dir+"rsem.real.sample"+str(i)):
            os.makedirs(out_dir+"rsem.real.sample"+str(i))
        cmd_real = ["rsem-prepare-reference",
                    out_dir+"res_distrib.real.sample"+str(i)+".fasta",
                    out_dir+"rsem.real.sample"+str(i)+"/rsem_ref"]
        subprocess.call(cmd_real)

        if not os.path.exists(out_dir+"rsem.nonint.sample"+str(i)):
            os.makedirs(out_dir+"rsem.nonint.sample"+str(i))
        cmd_nonint = ["rsem-prepare-reference",
                      out_dir+"res_distrib.nonint.sample"+str(i)+".fasta",
                      out_dir+"rsem.nonint.sample"+str(i)+"/rsem_ref"]
        subprocess.call(cmd_nonint)

        if not os.path.exists(out_dir+"rsem.int.sample"+str(i)):
            os.makedirs(out_dir+"rsem.int.sample"+str(i))
        cmd_int = ["rsem-prepare-reference",
                    out_dir+"res_distrib.int.sample"+str(i)+".fasta",
                    out_dir+"rsem.int.sample"+str(i)+"/rsem_ref"]
        subprocess.call(cmd_int)

        if not os.path.exists(out_dir+"rsem.pol.sample"+str(i)):
            os.makedirs(out_dir+"rsem.pol.sample"+str(i))
        cmd_pol = ["rsem-prepare-reference",
                    out_dir+"res_distrib.pol.sample"+str(i)+".fasta",
                    out_dir+"rsem.pol.sample"+str(i)+"/rsem_ref"]
        subprocess.call(cmd_pol)

    # now order expressions for rsem such that they appear in the same order as in the ref_name index
    for i in range(num_samples):
        order_real = []
        with open(out_dir+"rsem.real.sample0/rsem_ref.transcripts.fa","r") as inFP:
            for line in inFP.readlines():
                if line[0] == ">":
                    order_real.append(line[1:-1])
        odf_real = pd.DataFrame(order_real)
        odf_real.columns = ["transcript_id"]
        exp_real = pd.read_csv(out_dir+"/res_distrib.real.sample"+str(i)+".rsem_tpms.tmp",sep="\t")
        exp_real = odf_real.merge(exp_real,on="transcript_id",how="left",indicator=True)
        assert len(exp_real[exp_real["_merge"]=="both"])==len(exp_real),"incorrect tids"
        exp_real[rsem_tpm_cols].to_csv(out_dir+"/res_distrib.real.sample"+str(i)+".rsem_tpms",sep="\t",index=False)

        order_nonint = []
        with open(out_dir+"rsem.nonint.sample0/rsem_ref.transcripts.fa","r") as inFP:
            for line in inFP.readlines():
                if line[0] == ">":
                    order_nonint.append(line[1:-1])
        odf_nonint = pd.DataFrame(order_nonint)
        odf_nonint.columns = ["transcript_id"]
        exp_nonint = pd.read_csv(out_dir+"/res_distrib.nonint.sample"+str(i)+".rsem_tpms.tmp",sep="\t")
        exp_nonint = odf_nonint.merge(exp_nonint,on="transcript_id",how="left",indicator=True)
        assert len(exp_nonint[exp_nonint["_merge"]=="both"])==len(exp_nonint),"incorrect tids"
        exp_nonint[rsem_tpm_cols].to_csv(out_dir+"/res_distrib.nonint.sample"+str(i)+".rsem_tpms",sep="\t",index=False)

        order_int = []
        with open(out_dir+"rsem.int.sample0/rsem_ref.transcripts.fa","r") as inFP:
            for line in inFP.readlines():
                if line[0] == ">":
                    order_int.append(line[1:-1])
        odf_int = pd.DataFrame(order_int)
        odf_int.columns = ["transcript_id"]
        exp_int = pd.read_csv(out_dir+"/res_distrib.int.sample"+str(i)+".rsem_tpms.tmp",sep="\t")
        exp_int = odf_int.merge(exp_int,on="transcript_id",how="left",indicator=True)
        assert len(exp_int[exp_int["_merge"]=="both"])==len(exp_int),"incorrect tids"
        exp_int[rsem_tpm_cols].to_csv(out_dir+"/res_distrib.int.sample"+str(i)+".rsem_tpms",sep="\t",index=False)

        order_pol = []
        with open(out_dir+"rsem.pol.sample0/rsem_ref.transcripts.fa","r") as inFP:
            for line in inFP.readlines():
                if line[0] == ">":
                    order_pol.append(line[1:-1])
        odf_pol = pd.DataFrame(order_pol)
        odf_pol.columns = ["transcript_id"]
        exp_pol = pd.read_csv(out_dir+"/res_distrib.pol.sample"+str(i)+".rsem_tpms.tmp",sep="\t")
        exp_pol = odf_pol.merge(exp_pol,on="transcript_id",how="left",indicator=True)
        assert len(exp_pol[exp_pol["_merge"]=="both"])==len(exp_pol),"incorrect tids"
        exp_pol[rsem_tpm_cols].to_csv(out_dir+"/res_distrib.pol.sample"+str(i)+".rsem_tpms",sep="\t",index=False)

    if args.type == "rsem":
        # another solution is that we concatenate all transcripts and expressions together
        # then simulate with RSEM
        # and then break them apart into individual groups to proceed forward
        # for i in range(num_samples):
        #     with open(out_dir+"all_combined.sample"+str(i)+".gtf", 'w+') as outfile:
        #         for fname in [out_dir+"res_distrib.real.sample"+str(i)+".gtf",
        #                       out_dir+"res_distrib.nonint.sample"+str(i)+".gtf",
        #                       out_dir+"res_distrib.int.sample"+str(i)+".gtf",
        #                       out_dir+"res_distrib.pol.sample"+str(i)+".gtf"]:
        #             with open(fname) as infile:
        #                 for line in infile:
        #                     outfile.write(line)

        #     with open(out_dir+"all_combined.sample"+str(i)+".rsem_tpms", 'w+') as outfile:
        #         for fname in [out_dir+"res_distrib.real.sample"+str(i)+".rsem_tpms",
        #                       out_dir+"res_distrib.nonint.sample"+str(i)+".rsem_tpms",
        #                       out_dir+"res_distrib.int.sample"+str(i)+".rsem_tpms",
        #                       out_dir+"res_distrib.pol.sample"+str(i)+".rsem_tpms"]:
        #             with open(fname) as infile:
        #                 for line in infile:
        #                     outfile.write(line)

        # now we can simulate RSEM reads
        # for this we need to sample from the distribution of the number of reads per sample
        print(">>> Simulating reads with RSEM")
        # load distribution of the number of reads per sample
        nrdist_df = pd.read_csv(args.nrdist)
        nrdist = list(nrdist_df["readlen"].tolist())
        for i in range(num_samples):
            total_nr = random.choice(nrdist)
            print(i,tpm_contrib_real[i],tpm_contrib_nonint[i],tpm_contrib_int[i],tpm_contrib_pol[i])
            print(i,num_reads_real[i],num_reads_nonint[i],num_reads_int[i],num_reads_pol[i])
            total_tpm = tpm_contrib_real[i]+tpm_contrib_nonint[i]+tpm_contrib_int[i]+tpm_contrib_pol[i]

            if not os.path.exists(out_dir+"res_distrib.real.sample"+str(i)):
                os.makedirs(out_dir+"res_distrib.real.sample"+str(i))
            frac_real = tpm_contrib_real[i]/total_tpm
            real_nr = int(frac_real*total_nr)
            print("number of real reads: "+str(real_nr)+"\t"+str(num_reads_real[i]))
            cmd_real = ["rsem-simulate-reads",
                        out_dir+"rsem.real.sample"+str(i)+"/rsem_ref",
                        args.rsem_model,
                        out_dir+"/res_distrib.real.sample"+str(i)+".rsem_tpms",
                        str(0),
                        str(num_reads_real[i]),
                        out_dir+"res_distrib.real.sample"+str(i)+"/sample_01"]
            subprocess.call(cmd_real)

            if not os.path.exists(out_dir+"res_distrib.nonint.sample"+str(i)):
                os.makedirs(out_dir+"res_distrib.nonint.sample"+str(i))
            frac_nonint = tpm_contrib_nonint[i]/total_tpm
            nonint_nr = int(frac_nonint*total_nr)
            print("number of nonint reads: "+str(nonint_nr)+"\t"+str(num_reads_nonint[i]))
            cmd_nonint = ["rsem-simulate-reads",
                        out_dir+"rsem.nonint.sample"+str(i)+"/rsem_ref",
                        args.rsem_model,
                        out_dir+"/res_distrib.nonint.sample"+str(i)+".rsem_tpms",
                        str(0),
                        str(num_reads_nonint[i]),
                        out_dir+"res_distrib.nonint.sample"+str(i)+"/sample_01"]
            subprocess.call(cmd_nonint)

            if not os.path.exists(out_dir+"res_distrib.int.sample"+str(i)):
                os.makedirs(out_dir+"res_distrib.int.sample"+str(i))
            frac_int = tpm_contrib_int[i]/total_tpm
            int_nr = int(frac_int*total_nr)
            print("number of int reads: "+str(int_nr)+"\t"+str(num_reads_int[i]))
            cmd_int = ["rsem-simulate-reads",
                        out_dir+"rsem.int.sample"+str(i)+"/rsem_ref",
                        args.rsem_model,
                        out_dir+"/res_distrib.int.sample"+str(i)+".rsem_tpms",
                        str(0),
                        str(num_reads_int[i]),
                        out_dir+"res_distrib.int.sample"+str(i)+"/sample_01"]
            subprocess.call(cmd_int)

            if not os.path.exists(out_dir+"res_distrib.pol.sample"+str(i)):
                os.makedirs(out_dir+"res_distrib.pol.sample"+str(i))
            frac_pol = tpm_contrib_pol[i]/total_tpm
            pol_nr = int(frac_pol*total_nr)
            print("number of pol reads: "+str(pol_nr)+"\t"+str(num_reads_pol[i]))
            cmd_pol = ["rsem-simulate-reads",
                        out_dir+"rsem.pol.sample"+str(i)+"/rsem_ref",
                        args.rsem_model,
                        out_dir+"/res_distrib.pol.sample"+str(i)+".rsem_tpms",
                        str(0),
                        str(num_reads_pol[i]),
                        out_dir+"res_distrib.pol.sample"+str(i)+"/sample_01"]
            subprocess.call(cmd_pol)

    else:
        print(">>> Simulating reads with Polyester")
        for i in range(num_samples):
            if not os.path.exists(out_dir+"res_distrib.real.sample"+str(i)):
                os.makedirs(out_dir+"res_distrib.real.sample"+str(i))
            cmd_real = [genRNAseq,
                        out_dir+"res_distrib.real.sample"+str(i)+".fasta",
                        out_dir+"res_distrib.real.sample"+str(i)+".exp",
                        str(readlen),
                        out_dir+"res_distrib.real.sample"+str(i)+"/",
                        str(int(args.paired))]
            subprocess.call(cmd_real)
            
            if not os.path.exists(out_dir+"res_distrib.nonint.sample"+str(i)):
                os.makedirs(out_dir+"res_distrib.nonint.sample"+str(i))
            cmd_nonint = [genRNAseq,
                        out_dir+"res_distrib.nonint.sample"+str(i)+".fasta",
                        out_dir+"res_distrib.nonint.sample"+str(i)+".exp",
                        str(readlen),
                        out_dir+"res_distrib.nonint.sample"+str(i)+"/",
                        str(int(args.paired))]
            subprocess.call(cmd_nonint)
            
            if not os.path.exists(out_dir+"res_distrib.int.sample"+str(i)):
                os.makedirs(out_dir+"res_distrib.int.sample"+str(i))
            cmd_int = [genRNAseq,
                        out_dir+"res_distrib.int.sample"+str(i)+".fasta",
                        out_dir+"res_distrib.int.sample"+str(i)+".exp",
                        str(readlen),
                        out_dir+"res_distrib.int.sample"+str(i)+"/",
                        str(int(args.paired))]
            subprocess.call(cmd_int)
            
            if not os.path.exists(out_dir+"res_distrib.pol.sample"+str(i)):
                os.makedirs(out_dir+"res_distrib.pol.sample"+str(i))
            cmd_pol = [genRNAseq,
                        out_dir+"res_distrib.pol.sample"+str(i)+".fasta",
                        out_dir+"res_distrib.pol.sample"+str(i)+".exp",
                        str(readlen),
                        out_dir+"res_distrib.pol.sample"+str(i)+"/",
                        str(int(args.paired))]
            subprocess.call(cmd_pol)

    # build salmon and hisat2 indices
    f=open(out_dir+"res_distrib.all.real.exons","wb")
    subprocess.call([args.hisat+"_extract_exons.py",
                     base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.true.gtf"],
                     stdout=f)
    f.close()
    f=open(out_dir+"res_distrib.all.real.ss","wb")
    subprocess.call([args.hisat+"_extract_splice_sites.py",
                     base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.true.gtf"],
                     stdout=f)
    f.close()
    subprocess.call([args.hisat+"-build",
                     args.ref,
                     out_dir+"res_distrib.all.real",
                     "-p",str(args.threads),
                     "--ss",out_dir+"res_distrib.all.real.ss",
                     "--exon",out_dir+"res_distrib.all.real.exons"])

    subprocess.call(["salmon","index",
                     "-t",out_dir+"res_distrib.all.real.fasta",
                     "-p",str(args.threads),
                     "-i",out_dir+"res_distrib.all.real"])

    subprocess.call(["kallisto","index",
                     "-i",out_dir+"res_distrib.all.real.kallisto",
                     out_dir+"res_distrib.all.real.fasta"])

    # shuffle generated reads and create merged datasets for downstream analysis
    print(">>>shuffling reads and creating combinations")
    for i in range(num_samples):
        print(i)
        # shuffle real
        shuffle_cmd = []
        if(args.paired):
            shuffle_cmd = [args.shuffleReads,
                           out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_1"+extension,
                           out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_2"+extension,
                           out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_1.shuffled.fasta",
                           out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_2.shuffled.fasta"]
        else:
            shuffle_cmd = [args.shuffleReads,
                           out_dir+"res_distrib.real.sample"+str(i)+"/sample_01"+extension,
                           out_dir+"res_distrib.real.sample"+str(i)+"/sample_01.shuffled.fasta"]
        subprocess.call(shuffle_cmd)
        # create a combination of real and nonint
        if not os.path.exists(out_dir+"res_distrib.real_nonint.sample"+str(i)):
            os.makedirs(out_dir+"res_distrib.real_nonint.sample"+str(i))
        if args.paired:
            with open(out_dir+"res_distrib.real_nonint.sample"+str(i)+"/sample_01_1"+extension, 'w+') as outfile:
                for fname in [out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_1"+extension,
                              out_dir+"res_distrib.nonint.sample"+str(i)+"/sample_01_1"+extension]:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
            with open(out_dir+"res_distrib.real_nonint.sample"+str(i)+"/sample_01_2"+extension, 'w+') as outfile:
                for fname in [out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_2"+extension,
                              out_dir+"res_distrib.nonint.sample"+str(i)+"/sample_01_2"+extension]:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
        else:
            with open(out_dir+"res_distrib.real_nonint.sample"+str(i)+"/sample_01"+extension, 'w+') as outfile:
                for fname in [out_dir+"res_distrib.real.sample"+str(i)+"/sample_01"+extension,
                              out_dir+"res_distrib.nonint.sample"+str(i)+"/sample_01"+extension]:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
        shuffle_cmd = []
        if(args.paired):
            shuffle_cmd = [args.shuffleReads,
                           out_dir+"res_distrib.real_nonint.sample"+str(i)+"/sample_01_1"+extension,
                           out_dir+"res_distrib.real_nonint.sample"+str(i)+"/sample_01_2"+extension,
                           out_dir+"res_distrib.real_nonint.sample"+str(i)+"/sample_01_1.shuffled.fasta",
                           out_dir+"res_distrib.real_nonint.sample"+str(i)+"/sample_01_2.shuffled.fasta"]
        else:
            shuffle_cmd = [args.shuffleReads,
                           out_dir+"res_distrib.real_nonint.sample"+str(i)+"/sample_01"+extension,
                           out_dir+"res_distrib.real_nonint.sample"+str(i)+"/sample_01.shuffled.fasta"]
        subprocess.call(shuffle_cmd)
        # create a combination of real and int
        if not os.path.exists(out_dir+"res_distrib.real_int.sample"+str(i)):
            os.makedirs(out_dir+"res_distrib.real_int.sample"+str(i))
        if args.paired:
            with open(out_dir+"res_distrib.real_int.sample"+str(i)+"/sample_01_1"+extension, 'w+') as outfile:
                for fname in [out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_1"+extension,
                              out_dir+"res_distrib.int.sample"+str(i)+"/sample_01_1"+extension]:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
            with open(out_dir+"res_distrib.real_int.sample"+str(i)+"/sample_01_2"+extension, 'w+') as outfile:
                for fname in [out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_2"+extension,
                              out_dir+"res_distrib.int.sample"+str(i)+"/sample_01_2"+extension]:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
        else:
            with open(out_dir+"res_distrib.real_int.sample"+str(i)+"/sample_01"+extension, 'w+') as outfile:
                for fname in [out_dir+"res_distrib.real.sample"+str(i)+"/sample_01"+extension,
                              out_dir+"res_distrib.int.sample"+str(i)+"/sample_01"+extension]:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
        shuffle_cmd = []
        if(args.paired):
            shuffle_cmd = [args.shuffleReads,
                           out_dir+"res_distrib.real_int.sample"+str(i)+"/sample_01_1"+extension,
                           out_dir+"res_distrib.real_int.sample"+str(i)+"/sample_01_2"+extension,
                           out_dir+"res_distrib.real_int.sample"+str(i)+"/sample_01_1.shuffled.fasta",
                           out_dir+"res_distrib.real_int.sample"+str(i)+"/sample_01_2.shuffled.fasta"]
        else:
            shuffle_cmd = [args.shuffleReads,
                           out_dir+"res_distrib.real_int.sample"+str(i)+"/sample_01"+extension,
                           out_dir+"res_distrib.real_int.sample"+str(i)+"/sample_01.shuffled.fasta"]
        subprocess.call(shuffle_cmd)
        # create a combination of real and polymerase
        if not os.path.exists(out_dir+"res_distrib.real_pol.sample"+str(i)):
            os.makedirs(out_dir+"res_distrib.real_pol.sample"+str(i))
        if args.paired:
            with open(out_dir+"res_distrib.real_pol.sample"+str(i)+"/sample_01_1"+extension, 'w+') as outfile:
                for fname in [out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_1"+extension,
                              out_dir+"res_distrib.pol.sample"+str(i)+"/sample_01_1"+extension]:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
            with open(out_dir+"res_distrib.real_pol.sample"+str(i)+"/sample_01_2"+extension, 'w+') as outfile:
                for fname in [out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_2"+extension,
                              out_dir+"res_distrib.pol.sample"+str(i)+"/sample_01_2"+extension]:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
        else:
            with open(out_dir+"res_distrib.real_pol.sample"+str(i)+"/sample_01"+extension, 'w+') as outfile:
                for fname in [out_dir+"res_distrib.real.sample"+str(i)+"/sample_01"+extension,
                              out_dir+"res_distrib.pol.sample"+str(i)+"/sample_01"+extension]:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
        shuffle_cmd = []
        if(args.paired):
            shuffle_cmd = [args.shuffleReads,
                           out_dir+"res_distrib.real_pol.sample"+str(i)+"/sample_01_1"+extension,
                           out_dir+"res_distrib.real_pol.sample"+str(i)+"/sample_01_2"+extension,
                           out_dir+"res_distrib.real_pol.sample"+str(i)+"/sample_01_1.shuffled.fasta",
                           out_dir+"res_distrib.real_pol.sample"+str(i)+"/sample_01_2.shuffled.fasta"]
        else:
            shuffle_cmd = [args.shuffleReads,
                           out_dir+"res_distrib.real_pol.sample"+str(i)+"/sample_01"+extension,
                           out_dir+"res_distrib.real_pol.sample"+str(i)+"/sample_01.shuffled.fasta"]
        subprocess.call(shuffle_cmd)
        # create a combination of all reads
        if not os.path.exists(out_dir+"res_distrib.all.sample"+str(i)):
            os.makedirs(out_dir+"res_distrib.all.sample"+str(i))
        if args.paired:
            with open(out_dir+"res_distrib.all.sample"+str(i)+"/sample_01_1"+extension, 'w+') as outfile:
                for fname in [out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_1"+extension,
                              out_dir+"res_distrib.nonint.sample"+str(i)+"/sample_01_1"+extension,
                              out_dir+"res_distrib.int.sample"+str(i)+"/sample_01_1"+extension,
                              out_dir+"res_distrib.pol.sample"+str(i)+"/sample_01_1"+extension]:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
            with open(out_dir+"res_distrib.all.sample"+str(i)+"/sample_01_2"+extension, 'w+') as outfile:
                for fname in [out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_2"+extension,
                              out_dir+"res_distrib.nonint.sample"+str(i)+"/sample_01_2"+extension,
                              out_dir+"res_distrib.int.sample"+str(i)+"/sample_01_2"+extension,
                              out_dir+"res_distrib.pol.sample"+str(i)+"/sample_01_2"+extension]:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
        else:
            with open(out_dir+"res_distrib.all.sample"+str(i)+"/sample_01"+extension, 'w+') as outfile:
                for fname in [out_dir+"res_distrib.real.sample"+str(i)+"/sample_01"+extension,
                              out_dir+"res_distrib.nonint.sample"+str(i)+"/sample_01"+extension,
                              out_dir+"res_distrib.int.sample"+str(i)+"/sample_01"+extension,
                              out_dir+"res_distrib.pol.sample"+str(i)+"/sample_01"+extension]:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
        shuffle_cmd = []
        if(args.paired):
            shuffle_cmd = [args.shuffleReads,
                           out_dir+"res_distrib.all.sample"+str(i)+"/sample_01_1"+extension,
                           out_dir+"res_distrib.all.sample"+str(i)+"/sample_01_2"+extension,
                           out_dir+"res_distrib.all.sample"+str(i)+"/sample_01_1.shuffled.fasta",
                           out_dir+"res_distrib.all.sample"+str(i)+"/sample_01_2.shuffled.fasta"]
        else:
            shuffle_cmd = [args.shuffleReads,
                           out_dir+"res_distrib.all.sample"+str(i)+"/sample_01"+extension,
                           out_dir+"res_distrib.all.sample"+str(i)+"/sample_01.shuffled.fasta"]
        subprocess.call(shuffle_cmd)

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
                     "-S",out_dir+"strg.real.sample"+str(i)+"/hisat.sam"]
        if args.paired:
            hisat_cmd.extend(["-1",out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_1.shuffled.fasta",
                              "-2",out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_2.shuffled.fasta"])
        else:
            hisat_cmd.extend(["-U",out_dir+"res_distrib.real.sample"+str(i)+"/sample_01.shuffled.fasta"])
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
                    "-p",str(args.threads),
                    "-G",base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.true.gtf",
                    "-o",out_dir+"strg.real.sample"+str(i)+"/strg.gtf"]
        subprocess.call(strg_cmd)
        # quantify with salmon
        salmon_cmd = ["salmon","quant","--validateMappings","-l","A",
                      "-i",out_dir+"res_distrib.all.real",
                      "-p",str(args.threads),
                      "-o",out_dir+"slmn.real.sample"+str(i)]
        if args.paired:
            salmon_cmd.extend(["-1",out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_1.shuffled.fasta",
                               "-2",out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_2.shuffled.fasta"])
        else:
            salmon_cmd.extend(["-r",out_dir+"res_distrib.real.sample"+str(i)+"/sample_01.shuffled.fasta"])
        subprocess.call(salmon_cmd)

        # next perform the same analysis with real and non-intronic reads
        if not os.path.exists(out_dir+"strg.real_nonint.sample"+str(i)):
            os.makedirs(out_dir+"strg.real_nonint.sample"+str(i))
        # align
        hisat_cmd = [args.hisat,
                     "-x",out_dir+"res_distrib.all.real",
                     "-p",str(args.threads),
                     "--rna-sensitive","-f",
                     "-S",out_dir+"strg.real_nonint.sample"+str(i)+"/hisat.sam"]
        if args.paired:
            hisat_cmd.extend(["-1",out_dir+"res_distrib.real_nonint.sample"+str(i)+"/sample_01_1.shuffled.fasta",
                              "-2",out_dir+"res_distrib.real_nonint.sample"+str(i)+"/sample_01_2.shuffled.fasta"])
        else:
            hisat_cmd.extend(["-U",out_dir+"res_distrib.real_nonint.sample"+str(i)+"/sample_01.shuffled.fasta"])
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
                    "-p",str(args.threads),
                    "-G",base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.true.gtf",
                    "-o",out_dir+"strg.real_nonint.sample"+str(i)+"/strg.gtf"]
        subprocess.call(strg_cmd)
        # quantify with salmon
        salmon_cmd = ["salmon","quant","--validateMappings","-l","A",
                      "-i",out_dir+"res_distrib.all.real",
                      "-p",str(args.threads),
                      "-o",out_dir+"slmn.real_nonint.sample"+str(i)]
        if args.paired:
            salmon_cmd.extend(["-1",out_dir+"res_distrib.real_nonint.sample"+str(i)+"/sample_01_1.shuffled.fasta",
                               "-2",out_dir+"res_distrib.real_nonint.sample"+str(i)+"/sample_01_2.shuffled.fasta"])
        else:
            salmon_cmd.extend(["-r",out_dir+"res_distrib.real_nonint.sample"+str(i)+"/sample_01.shuffled.fasta"])
        subprocess.call(salmon_cmd)

        # next perform the same analysis with real and intronic reads
        if not os.path.exists(out_dir+"strg.real_int.sample"+str(i)):
            os.makedirs(out_dir+"strg.real_int.sample"+str(i))
        # align
        hisat_cmd = [args.hisat,
                     "-x",out_dir+"res_distrib.all.real",
                     "-p",str(args.threads),
                     "--rna-sensitive","-f",
                     "-S",out_dir+"strg.real_int.sample"+str(i)+"/hisat.sam"]
        if args.paired:
            hisat_cmd.extend(["-1",out_dir+"res_distrib.real_int.sample"+str(i)+"/sample_01_1.shuffled.fasta",
                              "-2",out_dir+"res_distrib.real_int.sample"+str(i)+"/sample_01_2.shuffled.fasta"])
        else:
            hisat_cmd.extend(["-U",out_dir+"res_distrib.real_int.sample"+str(i)+"/sample_01.shuffled.fasta"])
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
                    "-p",str(args.threads),
                    "-G",base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.true.gtf",
                    "-o",out_dir+"strg.real_int.sample"+str(i)+"/strg.gtf"]
        subprocess.call(strg_cmd)
        # quantify with salmon
        salmon_cmd = ["salmon","quant","--validateMappings","-l","A",
                      "-i",out_dir+"res_distrib.all.real",
                      "-p",str(args.threads),
                      "-o",out_dir+"slmn.real_int.sample"+str(i)]
        if args.paired:
            salmon_cmd.extend(["-1",out_dir+"res_distrib.real_int.sample"+str(i)+"/sample_01_1.shuffled.fasta",
                               "-2",out_dir+"res_distrib.real_int.sample"+str(i)+"/sample_01_2.shuffled.fasta"])
        else:
            salmon_cmd.extend(["-r",out_dir+"res_distrib.real_int.sample"+str(i)+"/sample_01.shuffled.fasta"])
        subprocess.call(salmon_cmd)

        # next perform the same analysis with real and polymerase reads
        if not os.path.exists(out_dir+"strg.real_pol.sample"+str(i)):
            os.makedirs(out_dir+"strg.real_pol.sample"+str(i))
        # align
        hisat_cmd = [args.hisat,
                     "-x",out_dir+"res_distrib.all.real",
                     "-p",str(args.threads),
                     "--rna-sensitive","-f",
                     "-S",out_dir+"strg.real_pol.sample"+str(i)+"/hisat.sam"]
        if args.paired:
            hisat_cmd.extend(["-1",out_dir+"res_distrib.real_pol.sample"+str(i)+"/sample_01_1.shuffled.fasta",
                              "-2",out_dir+"res_distrib.real_pol.sample"+str(i)+"/sample_01_2.shuffled.fasta"])
        else:
            hisat_cmd.extend(["-U",out_dir+"res_distrib.real_pol.sample"+str(i)+"/sample_01.shuffled.fasta"])
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
                    "-p",str(args.threads),
                    "-G",base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.true.gtf",
                    "-o",out_dir+"strg.real_pol.sample"+str(i)+"/strg.gtf"]
        subprocess.call(strg_cmd)
        # quantify with salmon
        salmon_cmd = ["salmon","quant","--validateMappings","-l","A",
                      "-i",out_dir+"res_distrib.all.real",
                      "-p",str(args.threads),
                      "-o",out_dir+"slmn.real_pol.sample"+str(i)]
        if args.paired:
            salmon_cmd.extend(["-1",out_dir+"res_distrib.real_pol.sample"+str(i)+"/sample_01_1.shuffled.fasta",
                               "-2",out_dir+"res_distrib.real_pol.sample"+str(i)+"/sample_01_2.shuffled.fasta"])
        else:
            salmon_cmd.extend(["-r",out_dir+"res_distrib.real_pol.sample"+str(i)+"/sample_01.shuffled.fasta"])
        subprocess.call(salmon_cmd)

        # next perform the same analysis with all reads
        if not os.path.exists(out_dir+"strg.all.sample"+str(i)):
            os.makedirs(out_dir+"strg.all.sample"+str(i))
        # align
        hisat_cmd = [args.hisat,
                     "-x",out_dir+"res_distrib.all.real",
                     "-p",str(args.threads),
                     "--rna-sensitive","-f",
                     "-S",out_dir+"strg.all.sample"+str(i)+"/hisat.sam"]
        if args.paired:
            hisat_cmd.extend(["-1",out_dir+"res_distrib.all.sample"+str(i)+"/sample_01_1.shuffled.fasta",
                              "-2",out_dir+"res_distrib.all.sample"+str(i)+"/sample_01_2.shuffled.fasta"])
        else:
            hisat_cmd.extend(["-U",out_dir+"res_distrib.all.sample"+str(i)+"/sample_01.shuffled.fasta"])
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
                    "-p",str(args.threads),
                    "-G",base_dir_data+"ALL.combined.IDs.olny.in.ALL.combined.true.gtf",
                    "-o",out_dir+"strg.all.sample"+str(i)+"/strg.gtf"]
        subprocess.call(strg_cmd)
        # quantify with salmon
        salmon_cmd = ["salmon","quant","--validateMappings","-l","A",
                      "-i",out_dir+"res_distrib.all.real",
                      "-p",str(args.threads),
                      "-o",out_dir+"slmn.all.sample"+str(i)]
        if args.paired:
            salmon_cmd.extend(["-1",out_dir+"res_distrib.all.sample"+str(i)+"/sample_01_1.shuffled.fasta",
                               "-2",out_dir+"res_distrib.all.sample"+str(i)+"/sample_01_2.shuffled.fasta"])
        else:
            salmon_cmd.extend(["-r",out_dir+"res_distrib.all.sample"+str(i)+"/sample_01.shuffled.fasta"])
        subprocess.call(salmon_cmd)

    # now running kallisto
    for i in range(num_samples):
        klst_cmd = ["kallisto","quant",
                    "-i",out_dir+"res_distrib.all.real.kallisto",
                    "-o",out_dir+"klst.real.sample"+str(i),
                    "-t",str(args.threads),
                    "-l",str(250),"-s",str(25)]
        if args.paired:
            klst_cmd.extend([out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_1.shuffled.fasta",out_dir+"res_distrib.real.sample"+str(i)+"/sample_01_2.shuffled.fasta"])
        else:
            klst_cmd.extend(["--single",out_dir+"res_distrib.real.sample"+str(i)+"/sample_01.shuffled.fasta"])
        subprocess.call(klst_cmd)
        klst_cmd = ["kallisto","quant",
                    "-i",out_dir+"res_distrib.all.real.kallisto",
                    "-o",out_dir+"klst.real_nonint.sample"+str(i),
                    "-t",str(args.threads),
                    "-l",str(250),"-s",str(25)]
        if args.paired:
            klst_cmd.extend([out_dir+"res_distrib.real_nonint.sample"+str(i)+"/sample_01_1.shuffled.fasta",out_dir+"res_distrib.real_nonint.sample"+str(i)+"/sample_01_2.shuffled.fasta"])
        else:
            klst_cmd.extend(["--single",out_dir+"res_distrib.real_nonint.sample"+str(i)+"/sample_01.shuffled.fasta"])
        subprocess.call(klst_cmd)
        klst_cmd = ["kallisto","quant",
                    "-i",out_dir+"res_distrib.all.real.kallisto",
                    "-o",out_dir+"klst.real_int.sample"+str(i),
                    "-t",str(args.threads),
                    "-l",str(250),"-s",str(25)]
        if args.paired:
            klst_cmd.extend([out_dir+"res_distrib.real_int.sample"+str(i)+"/sample_01_1.shuffled.fasta",out_dir+"res_distrib.real_int.sample"+str(i)+"/sample_01_2.shuffled.fasta"])
        else:
            klst_cmd.extend(["--single",out_dir+"res_distrib.real_int.sample"+str(i)+"/sample_01.shuffled.fasta"])
        subprocess.call(klst_cmd)
        klst_cmd = ["kallisto","quant",
                    "-i",out_dir+"res_distrib.all.real.kallisto",
                    "-o",out_dir+"klst.real_pol.sample"+str(i),
                    "-t",str(args.threads),
                    "-l",str(250),"-s",str(25)]
        if args.paired:
            klst_cmd.extend([out_dir+"res_distrib.real_pol.sample"+str(i)+"/sample_01_1.shuffled.fasta",out_dir+"res_distrib.real_pol.sample"+str(i)+"/sample_01_2.shuffled.fasta"])
        else:
            klst_cmd.extend(["--single",out_dir+"res_distrib.real_pol.sample"+str(i)+"/sample_01.shuffled.fasta"])
        subprocess.call(klst_cmd)
        klst_cmd = ["kallisto","quant",
                    "-i",out_dir+"res_distrib.all.real.kallisto",
                    "-o",out_dir+"klst.all.sample"+str(i),
                    "-t",str(args.threads),
                    "-l",str(250),"-s",str(25)]
        if args.paired:
            klst_cmd.extend([out_dir+"res_distrib.all.sample"+str(i)+"/sample_01_1.shuffled.fasta",out_dir+"res_distrib.all.sample"+str(i)+"/sample_01_2.shuffled.fasta"])
        else:
            klst_cmd.extend(["--single",out_dir+"res_distrib.all.sample"+str(i)+"/sample_01.shuffled.fasta"])
        subprocess.call(klst_cmd)

    # lastly run same samples using sim2sam
    if not args.paired:
        for i in range(num_samples):
            plst_cmd = [args.sim2sam,
                        "-i",args.ref+".fai",
                        "-o",out_dir+"strg.real.sample"+str(i)+".sam",
                        "-s",out_dir+"res_distrib.real.sample"+str(i)+"/sample_01"+extension,
                        "-g",out_dir+"res_distrib.real.sample"+str(i)+".gtf"]
            if args.type == "rsem":
                plst_cmd.extend(["-t","rsem","-r",out_dir+"res_distrib.real.sample"+str(i)+"/sample_01.sim.isoforms.results"])
            else:
                plst_cmd.extend(["-t","poly"])
            subprocess.call(plst_cmd)
            sort_cmd = ["samtools","sort",
                        "-@",str(args.threads),
                        "-o",out_dir+"strg.real.sample"+str(i)+".sorted.bam",
                        out_dir+"strg.real.sample"+str(i)+".sam"]
            subprocess.call(sort_cmd)

            plst_cmd = [args.sim2sam,
                        "-i",args.ref+".fai",
                        "-o",out_dir+"strg.nonint.sample"+str(i)+".sam",
                        "-s",out_dir+"res_distrib.nonint.sample"+str(i)+"/sample_01"+extension,
                        "-g",out_dir+"res_distrib.nonint.sample"+str(i)+".gtf"]
            if args.type == "rsem":
                plst_cmd.extend(["-t","rsem","-r",out_dir+"res_distrib.nonint.sample"+str(i)+"/sample_01.sim.isoforms.results"])
            else:
                plst_cmd.extend(["-t","poly"])
            subprocess.call(plst_cmd)
            sort_cmd = ["samtools","sort",
                        "-@",str(args.threads),
                        "-o",out_dir+"strg.nonint.sample"+str(i)+".sorted.bam",
                        out_dir+"strg.nonint.sample"+str(i)+".sam"]
            subprocess.call(sort_cmd)

            plst_cmd = [args.sim2sam,
                        "-i",args.ref+".fai",
                        "-o",out_dir+"strg.int.sample"+str(i)+".sam",
                        "-s",out_dir+"res_distrib.int.sample"+str(i)+"/sample_01"+extension,
                        "-g",out_dir+"res_distrib.int.sample"+str(i)+".gtf"]
            if args.type == "rsem":
                plst_cmd.extend(["-t","rsem","-r",out_dir+"res_distrib.int.sample"+str(i)+"/sample_01.sim.isoforms.results"])
            else:
                plst_cmd.extend(["-t","poly"])
            subprocess.call(plst_cmd)
            sort_cmd = ["samtools","sort",
                        "-@",str(args.threads),
                        "-o",out_dir+"strg.int.sample"+str(i)+".sorted.bam",
                        out_dir+"strg.int.sample"+str(i)+".sam"]
            subprocess.call(sort_cmd)

            plst_cmd = [args.sim2sam,
                        "-i",args.ref+".fai",
                        "-o",out_dir+"strg.pol.sample"+str(i)+".sam",
                        "-s",out_dir+"res_distrib.pol.sample"+str(i)+"/sample_01"+extension,
                        "-g",out_dir+"res_distrib.pol.sample"+str(i)+".gtf"]
            if args.type == "rsem":
                plst_cmd.extend(["-t","rsem","-r",out_dir+"res_distrib.pol.sample"+str(i)+"/sample_01.sim.isoforms.results"])
            else:
                plst_cmd.extend(["-t","poly"])
            subprocess.call(plst_cmd)
            sort_cmd = ["samtools","sort",
                        "-@",str(args.threads),
                        "-o",out_dir+"strg.pol.sample"+str(i)+".sorted.bam",
                        out_dir+"strg.pol.sample"+str(i)+".sam"]
            subprocess.call(sort_cmd)

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
    parser.add_argument("--paired",
                        required=False,
                        action='store_true',
                        help="Whether to simulate paired reads")
    parser.add_argument("--sim2sam",
                        required=True,
                        type=str,
                        help="path to the sim2sam executable")
    parser.add_argument("--rsem_model",
                        required=True,
                        type=str,
                        help="path to the RSEM model parameters file")
    parser.add_argument("--type",
                        required=True,
                        type=str,
                        choices=['polyester','rsem'],
                        help="which simulation tool to use")
    parser.add_argument("--nrdist",
                        required=True,
                        type=str,
                        help="path to the file that contains the distribution of the number of reads per sample")

    parser.set_defaults(func=wrapper)
    args=parser.parse_args()
    args.func(args)


if __name__=="__main__":
    main(sys.argv[1:])


# Need to make sure the expressions start out from TPMs
# this way we can later enforce the number of reads for each group
# based on averages across samples, which needs to be computed by gtex_stats