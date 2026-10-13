#!/bin/bash
#SBATCH --mem 50G -c 4 -J kmap -p ei-medium

CPUS=4
ref=Urochloa_basilisk_9B9D_20251016
output=pheno_apo.k31.LRT0-01.Urochloa_basilisk_9B9D_20251016

source bowtie2-2.4.1_CBG
source samtools-1.11_CBG

zcat pheno_apo.txt_k31_bimbam_all.assoc.txt.assoc.txt.gz.LRT0-01.txt.gz | cut -f2 | bowtie2 -a --threads ${CPUS} --end-to-end --no-unal -x ${ref} -r -U - | samtools view -@ ${CPUS} -Sb > ${output}.bowtie.bam





