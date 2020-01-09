#!/usr/bin/env bash

if1=${1}
if2=${2}
of1=${3}
of2=${4}

paste <(cat ${if1}) <(cat ${if2}) | paste - - | shuf | awk -v a=${of1} -v b=${of2} -F'\t' '{OFS="\n"; print $1,$3 > a; print $2,$4 > b}'
