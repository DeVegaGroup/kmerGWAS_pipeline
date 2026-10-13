#!/bin/bash
#SBATCH --mem 50G -c 1 -p ei-medium
file=$1
cat ${file} | sed 's/$/£/g' | tr '1' '2' | sed 's/ /@/1' | tr -d ' ' | tr '@' '\t' | awk '{print $1 ", X, X, " $2}' | awk '{gsub(/.{1}/,"&, ",$4)}1' | sed 's/, £,//g' > ${file}.bimbam
