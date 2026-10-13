#!/bin/bash
#SBATCH --mem 1G -c 1

myfile=$1

zcat ${myfile} | awk '{if ($14<0.01) print $0}' | gzip > ${myfile}.LRT0-01.txt.gz
