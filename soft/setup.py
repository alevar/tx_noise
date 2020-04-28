#!/usr/bin/env python

# Performs typing and filtering and stats aggregation of the input dataset

from scipy import stats
import pandas as pd
import numpy as np
import subprocess
import argparse
import random
import glob
import math
import csv
import sys
import os

# ./setup.py --base_dir_data /ccb/salz8-1/avaraby/tx_noise/data/ --base_dir_out /ccb/salz8-1/avaraby/tx_noise/data/gtex_aggs/ --gtex_stats /ccb/salz8-1/avaraby/tx_noise/soft/gtex_stats/gtex_stats --filters /ccb/salz8-1/avaraby/tx_noise/soft/pipeline/filters.sh --threads 30 --seed 101 > /ccb/salz8-1/avaraby/tx_noise/data/setup.log

def setup(args):
    base_dir_data = os.path.abspath(args.base_dir_data)+"/"

    base_dir_out = os.path.abspath(args.base_dir_out)+"/"
    if not os.path.exists(base_dir_out):
        os.makedirs(base_dir_out)

    gff3cols=["seqid","source","type","start","end","score","strand","phase","attributes"]
    tmapcols=["ref_gene_id","ref_id","class_code","qry_gene_id","qry_id","num_exons","FPKM","TPM","cov","len","major_iso_id","ref_match_len"]

    # first run the filters.sh script to perform initial typing of transcripts and filtering
    # subprocess.call([args.filters,base_dir_data])

    # # extract a subset of CHESS comprised of these loci only - this will serve as a truth set
    # # first read in all the chess geneIDs
    # real_tm = pd.read_csv(base_dir_data+"real.chess.locs",names=tmapcols,sep="\t",usecols=["ref_id","qry_id"])
    # # only include those with CHS type IDs in the matching reference transcript
    # real_tm = real_tm[real_tm["ref_id"].str.contains("CHS\.[0-9]*\.[0-9]*")].reset_index(drop=True)
    # real_tm["gid"] = real_tm["ref_id"].str.split("\.",expand=True)[1].astype(int)
    # # convert chess GFF with matching transcript to GTF as real.gtf
    # assert len(real_tm[real_tm.duplicated("qry_id",keep="first")])==0,"duplicated ALL IDs in real"
    # real_tm.drop_duplicates("ref_id",keep="first",inplace=True)
    # assert len(real_tm[real_tm.duplicated("ref_id",keep="first")])==0,"duplicated CHESS IDs in real"
    # print("number of transcripts in real: "+str(len(set(real_tm["ref_id"]))))
    # real_tm.to_csv(base_dir_data+"real_chess2all.tgids",index=False)
    
    # chessDF = pd.read_csv(base_dir_data+"chess2.2_assembly.gff",sep="\t",comment="#",names=gff3cols)
    # chess_types = [x for x in set(chessDF["type"]) if x not in ["gene","CDS"]]
    # chessDF = chessDF[chessDF["type"].isin(chess_types)].reset_index(drop=True)
    # chessDF["gid"] = chessDF["attributes"].str.split("CHS.",expand=True,n=1)[1].str.split("\.|;",expand=True,n=1)[0].astype(int)
    # chessDF["tid"] = np.where(chessDF["type"]!="exon",
    #             chessDF["attributes"].str.split("ID=",expand=True,n=1)[1].str.split(";",expand=True,n=1)[0],
    #             chessDF["attributes"].str.split("Parent=",expand=True,n=1)[1].str.split(";",expand=True,n=1)[0])
    # chessDF = chessDF[chessDF["gid"].isin(real_tm["gid"])].reset_index(drop=True)
    # chessDF = chessDF.merge(real_tm,left_on="tid",right_on="ref_id",how="right",indicator=True)
    # assert len(chessDF[chessDF["_merge"]=="both"])==len(chessDF),"not all are merged for real"
    # assert len(real_tm["ref_id"])==len(set(chessDF["tid"])),"length of the dataframe has changed in real"
    # chessDF["attributes"] = "transcript_id \""+chessDF["tid"]+\
    #                         "\"; gene_id \"CHS."+chessDF["gid_x"].astype(str)+\
    #                         "\"; old_transcript_id_old \""+chessDF["qry_id"]+"\";"
    # chessDF[gff3cols].to_csv(base_dir_data+"chess_real_matches.gtf",header=False,index=False,sep="\t",quoting=csv.QUOTE_NONE)
    
    # # next we need to perform the same operation for the splicing noise
    # # first need to also retrieve chess IDs
    # chess_spliceloc_cmd = """awk -v pat='c|k|m|n|j|e' -F '\t' '$3~pat {print}' """+base_dir_data+"""nonreal.tmap > """+base_dir_data+"splice.chess.locs"
    # subprocess.call(chess_spliceloc_cmd,shell=True)

    # # now need to obtain a mapping from chess to ALL for the splicing noise transcripts
    # sp_tm = pd.read_csv(base_dir_data+"splice.chess.locs",names=tmapcols,sep="\t",usecols=["ref_id","qry_id"])
    # # only include those that belong to valid transcripts in CHESS
    # sp_tm = sp_tm[sp_tm["ref_id"].isin(real_tm["ref_id"])].reset_index(drop=True)
    # # only include those with CHS type IDs in the matching reference transcript
    # sp_tm = sp_tm[sp_tm["ref_id"].str.contains("CHS\.[0-9]*\.[0-9]*")].reset_index(drop=True)
    # sp_tm["gid"] = sp_tm["ref_id"].str.split("\.",expand=True)[1].astype(int)
    # # make sure the IDs exist in the set of real genes
    # sp_tm = sp_tm[sp_tm["gid"].isin(real_tm["gid"])]
    # sp_tm.to_csv(base_dir_data+"splicing_chess2all.tgids",index=False)

    # # next we need to perform the same operation for the intronic noise
    # # first need to also retrieve chess IDs
    # chess_intronloc_cmd = """awk -v pat='i' -F '\t' '$3~pat {print}' """+base_dir_data+"""nonreal.tmap > """+base_dir_data+"intron.chess.locs"
    # subprocess.call(chess_intronloc_cmd,shell=True)

    # # now need to obtain a mapping from chess to ALL for the intronic noise transcripts
    # int_tm = pd.read_csv(base_dir_data+"intron.chess.locs",names=tmapcols,sep="\t",usecols=["ref_id","qry_id"])
    # # only include those that belong to valid transcripts in CHESS
    # int_tm = int_tm[int_tm["ref_id"].isin(real_tm["ref_id"])].reset_index(drop=True)
    # # only include those with CHS type IDs in the matching reference transcript
    # int_tm = int_tm[int_tm["ref_id"].str.contains("CHS\.[0-9]*\.[0-9]*")].reset_index(drop=True)
    # int_tm["gid"] = int_tm["ref_id"].str.split("\.",expand=True)[1].astype(int)
    # int_tm.to_csv(base_dir_data+"intronic_chess2all.tgids",index=False)

    # all_qry_ids = set(real_tm["qry_id"]).union(set(sp_tm["qry_id"]).union(int_tm["qry_id"]))

    # # need to build a dictionary from qry_id to ref_id
    # all_tm = pd.concat([real_tm,sp_tm,int_tm],axis=0)
    # id_dict = pd.Series(all_tm.ref_id.values,index=all_tm.qry_id).to_dict()
    # real_dict = pd.Series(real_tm.ref_id.values,index=real_tm.qry_id).to_dict()
    # sp_dict = pd.Series(sp_tm.ref_id.values,index=sp_tm.qry_id).to_dict()
    # int_dict = pd.Series(int_tm.ref_id.values,index=int_tm.qry_id).to_dict()

    # assert set(id_dict)==set(all_tm["qry_id"]),"qry sets not identical"
    # assert set(id_dict.values())==set(all_tm["ref_id"]),"ref sets not identical"
    # assert len(list(id_dict))==len(list(all_tm["qry_id"].tolist())),"qry lists of different lengths"
    # assert len(list(id_dict.values()))==len(list(all_tm["ref_id"].tolist())),"ref lists of different lengths"

    # count=0
    # real_gtf_fp = open(base_dir_data+"real_1.gtf","w+")
    # sp_gtf_fp = open(base_dir_data+"splicing_1.gtf","w+")
    # int_gtf_fp = open(base_dir_data+"intronic_1.gtf","w+")

    # with open(base_dir_data+"ALL.merged.gtf","r") as inFP:
    #     for line in inFP.readlines():
    #         line = line.strip()
    #         lineCols = line.split("\t")
    #         count+=1
    #         if count%1000000==0:
    #             print(count/1000000)
                
    #         # get transcript and gene IDs
    #         attrs = lineCols[8].split("\"")
    #         tid = attrs[1]
    #         gid = attrs[3]
    #         ntid = ""
    #         ngid = ""
    #         # now find the corresponding chess IDs
    #         if not tid in id_dict: # nothing found
    #             continue
                
    #         if tid in real_dict:
    #             ntid = real_dict[tid]
    #             ngid = "CHS."+ntid.split(".")[1]
    #             new_line = "\t".join(lineCols[:8])\
    #                              +"\ttranscript_id \""+ntid\
    #                              +"\"; gene_id \""+ngid\
    #                              +"\"; old_transcript_id_old \""+tid\
    #                              +"\"; old_gene_id \""+gid+"\";\n"
    #             real_gtf_fp.write(new_line)
                
            
    #         elif tid in sp_dict:
    #             ntid = sp_dict[tid]
    #             ngid = "CHS."+ntid.split(".")[1]
    #             new_line = "\t".join(lineCols[:8])\
    #                              +"\ttranscript_id \""+tid\
    #                              +"\"; gene_id \""+ngid\
    #                              +"\"; old_transcript_id_old \""+tid\
    #                              +"\"; old_gene_id \""+gid+"\";\n"
    #             sp_gtf_fp.write(new_line)
            
    #         elif tid in int_dict:
    #             ntid = int_dict[tid]
    #             ngid = "CHS."+ntid.split(".")[1]
    #             new_line = "\t".join(lineCols[:8])\
    #                              +"\ttranscript_id \""+tid\
    #                              +"\"; gene_id \""+ngid\
    #                              +"\"; old_transcript_id_old \""+tid\
    #                              +"\"; old_gene_id \""+gid+"\";\n"
    #             int_gtf_fp.write(new_line)
                
    #         else: # nothing found
    #             print("ERROR")
    #             break;
            
    # real_gtf_fp.close()
    # sp_gtf_fp.close()
    # int_gtf_fp.close()

    # # lastly need to concatenate all 4 to create a new set of transcripts
    # # interesting to know how many remain in total
    # cat_cmd = "cat "+base_dir_data+"real_1.gtf "+base_dir_data+"splicing_1.gtf "+base_dir_data+"intronic_1.gtf "+base_dir_data+"RNApol_tmp.gtf > "+base_dir_data+"ALL.merged.filtered.gtf"
    # subprocess.call(cat_cmd,shell=True)

    # # now we need to modify intergenic set to replcae all undefined strands with real strand information
    # polDF = pd.read_csv(base_dir_data+"RNApol.gtf",sep="\t",comment="#",names=gff3cols)
    # polDF["tid"] = polDF["attributes"].str.split("transcript_id \"",expand=True,n=1)[1].str.split("\"",expand=True,n=1)[0]
    # undef_strand_tids = polDF[(polDF["type"]=="transcript")&(polDF["strand"]==".")].tid.tolist()
    # def_strand_tids = polDF[(polDF["type"]=="transcript")&~(polDF["strand"]==".")].tid.tolist()
    # assert len(undef_strand_tids)==len(set(undef_strand_tids)),"duplicates in undefined"
    # assert len(def_strand_tids)==len(set(def_strand_tids)),"duplicates in defined"
    # assert len(set(undef_strand_tids).intersection(set(def_strand_tids)))==0,"same transcripts in both groups"
    # random.seed(args.seed)
    # random_strands = random.choices(["+","-"],k=len(undef_strand_tids))
    # undef_strand_df = pd.concat([pd.DataFrame(undef_strand_tids,columns=["tid"]),\
    #                              pd.DataFrame(random_strands,columns=["strand_new"])],axis=1)
    # polDF = polDF.merge(undef_strand_df,on="tid",how="left")
    # polDF["strand"] = np.where(polDF["strand"]==".",polDF["strand_new"],polDF["strand"])
    # polDF[gff3cols].to_csv(base_dir_data+"RNApol.gtf",header=False,index=False,sep="\t",quoting=csv.QUOTE_NONE)

    # # first we need to load the RNApol set in order to know which transcripts to output for that set
    # polDF = pd.read_csv(base_dir_data+"RNApol.gtf",sep="\t",comment="#",names=gff3cols)
    # polDF = polDF[polDF["type"]=="transcript"].reset_index(drop=True)
    # polDF["tid"] = polDF["attributes"].str.split("transcript_id \"",expand=True,n=1)[1].str.split("\"",expand=True,n=1)[0]
    # pol_set = set(polDF["tid"])

    # # now we need to modify any other files with the updated information
    # # namely, we need to update the tracking file with the updated information
    # count=0

    # real_tr_fp = open(base_dir_data+"real.tracking","w+")
    # sp_tr_fp = open(base_dir_data+"splicing.tracking","w+")
    # int_tr_fp = open(base_dir_data+"intronic.tracking","w+")
    # pol_tr_fp = open(base_dir_data+"RNApol.tracking","w+")

    # with open(base_dir_data+"ALL.tracking","r") as inFP:
    #     for line in inFP.readlines():
    #         line = line.strip()
    #         lineCols = line.split("\t")
    #         count+=1
    #         if count%1000000==0:
    #             print(count/1000000)
                
    #         # get transcript and gene IDs
    #         tid = lineCols[0]
    #         # now find the corresponding chess IDs
    #         if not tid in id_dict and not tid in pol_set: # nothing found
    #             continue
                
    #         if tid in real_dict:
    #             ntid = real_dict[tid]
    #             ngid = "CHS."+ntid.split(".")[1]
    #             new_line = ntid+"\t"+ngid+"\t"+"\t".join(lineCols[2:])+"\n"
    #             real_tr_fp.write(new_line)
                
            
    #         elif tid in sp_dict:
    #             ntid = sp_dict[tid]
    #             ngid = "CHS."+ntid.split(".")[1]
    #             new_line = tid+"\t"+ngid+"\t"+"\t".join(lineCols[2:])+"\n"
    #             sp_tr_fp.write(new_line)
            
    #         elif tid in int_dict:
    #             ntid = int_dict[tid]
    #             ngid = "CHS."+ntid.split(".")[1]
    #             new_line = tid+"\t"+ngid+"\t"+"\t".join(lineCols[2:])+"\n"
    #             int_tr_fp.write(new_line)
        
    #         elif tid in pol_set:
    #             pol_tr_fp.write(line+"\n")
                
    #         else: # nothing found
    #             print("ERROR")
    #             break;
            
    # real_tr_fp.close()
    # sp_tr_fp.close()
    # int_tr_fp.close()
    # pol_tr_fp.close()

    # # lastly need to concatenate all 4 to create a new tracking file
    # cat_tr_cmd = "cat "+base_dir_data+"real.tracking "+\
    #           base_dir_data+"splicing.tracking "+\
    #           base_dir_data+"intronic.tracking "+\
    #           base_dir_data+"RNApol.tracking > "+base_dir_data+"ALL.merged.filtered.tracking"
    # subprocess.call(cat_tr_cmd,shell=True)

    # # deal with missing strand for real loci now
    # # load base annotations
    # print(">>>loading base annotations")
    # real_baseDF = pd.read_csv(base_dir_data+"real_1.gtf",sep="\t",names=gff3cols)
    # splice_baseDF = pd.read_csv(base_dir_data+"splicing_1.gtf",sep="\t",names=gff3cols)
    # int_baseDF = pd.read_csv(base_dir_data+"intronic_1.gtf",sep="\t",names=gff3cols)

    # # get all loci and transcript IDs
    # print(">>>getting loci IDs")
    # real_baseDF["lid"] = real_baseDF["attributes"].str.split("gene_id \"",expand=True,n=1)[1].str.split("\"",expand=True,n=1)[0]
    # splice_baseDF["lid"] = splice_baseDF["attributes"].str.split("gene_id \"",expand=True,n=1)[1].str.split("\"",expand=True,n=1)[0]
    # int_baseDF["lid"] = int_baseDF["attributes"].str.split("gene_id \"",expand=True,n=1)[1].str.split("\"",expand=True,n=1)[0]

    # chessDF = pd.read_csv(base_dir_data+"chess2.2_assembly.gff",sep="\t",comment="#",names=gff3cols)
    # chessDF = chessDF[chessDF["type"]=="gene"].reset_index(drop=True)
    # chessDF["lid"] = "CHS."+chessDF["attributes"].str.split("CHS.",expand=True,n=1)[1].str.split("\.|;",expand=True,n=1)[0]

    # print("getting real strands")
    # st_df = chessDF[["lid","strand"]].drop_duplicates(keep="first")
    # assert set(st_df["strand"])==set(["+","-"]),"wrong strands"
    # st_df.columns = ["lid","strand_new"]

    # print("setting real")
    # real_baseDF2=real_baseDF.merge(st_df,on="lid",how="left")
    # real_baseDF2["strand"] = real_baseDF2["strand_new"]
    # assert set(real_baseDF2["strand"]==set(["+","-"])),"wrong strands real"

    # print("setting intronic")
    # int_baseDF2=int_baseDF.merge(st_df,on="lid",how="left")
    # int_baseDF2["strand"] = int_baseDF2["strand_new"]
    # assert set(int_baseDF2["strand"])==set(["+","-"]),"wrong strands int"

    # print("setting splicing")
    # splice_baseDF2=splice_baseDF.merge(st_df,on="lid",how="left")
    # splice_baseDF2["strand"] = splice_baseDF2["strand_new"]
    # assert set(splice_baseDF2["strand"])==set(["+","-"]),"wrong strands splice"

    # real_baseDF2[gff3cols].to_csv(base_dir_data+"real.gtf",header=False,index=False,sep="\t",quoting=csv.QUOTE_NONE)
    # splice_baseDF2[gff3cols].to_csv(base_dir_data+"splicing.gtf",header=False,index=False,sep="\t",quoting=csv.QUOTE_NONE)
    # int_baseDF2[gff3cols].to_csv(base_dir_data+"intronic.gtf",header=False,index=False,sep="\t",quoting=csv.QUOTE_NONE)

    # mv_cmd = """mv """+base_dir_data+"""ALL.merged.filtered.gtf """+base_dir_data+"""ALL.merged.filtered_old.gtf"""
    # subprocess.call(mv_cmd,shell=True)

    # # lastly need to concatenate all 4 to create a new set of transcripts again
    # cat_cmd = "cat "+base_dir_data+"real.gtf "+base_dir_data+"splicing.gtf "+base_dir_data+"intronic.gtf "+base_dir_data+"RNApol.gtf > "+base_dir_data+"ALL.merged.filtered.gtf"
    # subprocess.call(cat_cmd,shell=True)

    print("running gtex_stats")
    agg_cmd = args.gtex_stats+" \
               -t "+base_dir_data+"tissues/tissues.lst \
               -a "+base_dir_data+"ALL.merged.filtered.tracking \
               -r "+base_dir_data+"real.gtf \
               -g "+base_dir_data+"ALL.merged.filtered.gtf \
               -s "+base_dir_data+"splicing.gtf \
               -i "+base_dir_data+"intronic.gtf \
               -p "+base_dir_data+"RNApol.gtf \
               -o "+base_dir_out+"res"
    subprocess.call(agg_cmd,shell=True)

def main(args):
    parser = argparse.ArgumentParser(description='''Help Page''')

#===========================================
#===================BUILD===================
#===========================================
    parser.add_argument('--base_dir_data',
                        required=True,
                        type=str,
                        help="path to the organized directory with base data. The results of filtering will be stored there as well")
    parser.add_argument('--base_dir_out',
                        required=True,
                        type=str,
                        help="path where the results of gtex_stats are to be stored")
    parser.add_argument('--gtex_stats',
                        required=True,
                        type=str,
                        help="path to gtex_stats executable")
    parser.add_argument('--filters',
                        required=True,
                        type=str,
                        help="path to the filters.sh script")
    parser.add_argument("--threads",
                        required=True,
                        type=int,
                        default=1,
                        help="number of threads permitted to use")
    parser.add_argument("--seed",
                        required=True,
                        type=int,
                        default=101,
                        help="seed for random number generator")

    parser.set_defaults(func=setup)
    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    main(sys.argv[1:])
