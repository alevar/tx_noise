#!/usr/bin/env bash

base_dir_data=${1}

if [ -z "$base_dir_data" ]
then
      echo "need to know base_dir_data"
fi

cd ${base_dir_data}

# cluster into loci
gffread ALL.gtf --cluster-only -T -o ALL.merged.gtf

# run gff compare between chess2.2_assembly and ALL
gffcompare -o chess2ALL -r chess2.2_assembly.gff ALL.merged.gtf

remove unnecessary data
rm chess2ALL.ALL.merged.gtf.refmap chess2ALL.annotated.gtf chess2ALL.tracking chess2ALL.loci chess2ALL.stats

# get all the real loci codes
awk -F '\t|\"' 'FNR==NR{ids[$1]++;next} {if($10 in ids) print $0}' <(awk -v pat='=' -F '\t' '$3~pat {print $5}' chess2ALL.ALL.merged.gtf.tmap) ALL.merged.gtf > real_tmp.gtf

# next get a list of all loci that contain real transcripts
awk -F '\t|locus "' '$3=="transcript" {print $10}' real_tmp.gtf | cut -d'"' -f1 | sort | uniq > real.locs

# extract all transcript IDs that are in the set of real loci
grep -w -F -f real.locs ALL.merged.gtf | awk -F '\t|\"' '$3=="transcript" {print $10}' | sort | uniq > real_loc.all_tids

# extract all transcript IDs
awk -F '\t|\"' '$3=="transcript" {print $10}' ALL.merged.gtf | sort | uniq > all.tids

# get tids of potentially RNApol
diff --new-line-format="" --unchanged-line-format=""  all.tids real_loc.all_tids > RNApol.tids

# get related RNApol tmap
grep -w -F -f RNApol.tids chess2ALL.ALL.merged.gtf.tmap > RNApol.tmap

# next, everything that is not in this set is RNApol noise
awk -F '\t|\"' 'FNR==NR{ids[$1]++;next} {if($10 in ids) print $0}' <(awk -v pat='u' -F '\t' '$3~pat {print $5}' RNApol.tmap) ALL.merged.gtf > RNApol_tmp.gtf

# next we need to separate out intronic and splicing noise
# first get a set of real tids
awk -F '\t|\"' '$3=="transcript" {print $10}' real_tmp.gtf | sort | uniq > real.tids

# next we can subtract this from the set of all tids in real locs
diff --new-line-format="" --unchanged-line-format=""  real_loc.all_tids real.tids > real_loc.nonreal_tids

# now we can generate the corresponding tmap file
grep -w -F -f real_loc.nonreal_tids chess2ALL.ALL.merged.gtf.tmap > nonreal.tmap

# now we can get intronic transcripts by following the code "i"
awk -F '\t|\"' 'FNR==NR{ids[$1]++;next} {if($10 in ids) print $0}' <(awk -v pat='i' -F '\t' '$3~pat {print $5}' nonreal.tmap) ALL.merged.gtf > intronic_tmp.gtf

# lastly we can extract splicing noise by following code "c,k,m,n,j,e"
awk -F '\t|\"' 'FNR==NR{ids[$1]++;next} {if($10 in ids) print $0}' <(awk -v pat='c|k|m|n|j|e' -F '\t' '$3~pat {print $5}' nonreal.tmap) ALL.merged.gtf > splicing_tmp.gtf

# now we can focus on figuring out how to extract real transcripts/loci
# and adding noise transcripts to them

# begin by extracting a list of chess IDs of the transcripts that match from ALL to CHESS
# and convert chess transcript IDs to gene IDs
awk -v pat='=' -F '\t' '$3~pat {print}' chess2ALL.ALL.merged.gtf.tmap > real.chess.locs