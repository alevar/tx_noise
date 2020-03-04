#!/usr/bin/env bash

if [ "$#" -eq 2 ]; then
	echo "shuffling single-end mode"
	if1=${1}
	of1=${2}
	paste ${if1} | paste - - | shuf | awk -v a=${of1} -F'\t' '{OFS="\n"; print $1,$2 > a}'
elif [ "$#" -eq 4 ]; then
	echo "shuffling paired-end mode"
	if1=${1}
	if2=${2}
	of1=${3}
	of2=${4}
	paste <(cat ${if1}) <(cat ${if2}) | paste - - | shuf | awk -v a=${of1} -v b=${of2} -F'\t' '{OFS="\n"; print $1,$3 > a; print $2,$4 > b}'
else
	echo "incorrect number of parameters"
fi