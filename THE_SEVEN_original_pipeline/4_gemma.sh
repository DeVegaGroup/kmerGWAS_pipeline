#!/bin/bash
#SBATCH --mem 200G -c 1 -p ei-medium

#for i in {1..91}; do
#~/bin/gemma -gk 1 -g k51_bimbam_${SLURM_ARRAY_TASK_ID} -p phenotypes.txt -n ${i} -outdir assoc_${SLURM_ARRAY_TASK_ID} -o kinship_${SLURM_ARRAY_TASK_ID}_pheno_${i}.txt
#~/bin/gemma -lmm 4 -g k51_bimbam_${SLURM_ARRAY_TASK_ID} -p phenotypes.txt -n ${i} -k assoc_${SLURM_ARRAY_TASK_ID}/kinship_${SLURM_ARRAY_TASK_ID}_pheno_${i}.txt.cXX.txt -outdir assoc_${SLURM_ARRAY_TASK_ID} -o assoc_${SLURM_ARRAY_TASK_ID}_pheno_${i}.assoc.txt -maf 0.1
#gzip assoc_${SLURM_ARRAY_TASK_ID}/kinship_${SLURM_ARRAY_TASK_ID}_pheno_${i}.txt.cXX.txt
#cat assoc_${SLURM_ARRAY_TASK_ID}/assoc_${SLURM_ARRAY_TASK_ID}_pheno_${i}.assoc.txt.assoc.txt | cut -f2,15 | gzip > assoc_${SLURM_ARRAY_TASK_ID}/assoc_${SLURM_ARRAY_TASK_ID}_pheno_${i}.assoc.txt.assoc.txt.gz && rm assoc_${SLURM_ARRAY_TASK_ID}/assoc_${SLURM_ARRAY_TASK_ID}_pheno_${i}.assoc.txt.assoc.txt
#done && touch k51_bimbam_${SLURM_ARRAY_TASK_ID}.done


genotype=$1
phenotype=pheno_apo.txt

~/bin/gemma -gk 1 -g ${genotype} -p ${phenotype} -outdir . -o ${genotype}.kinship.txt
~/bin/gemma -lmm 4 -maf 0.01 -g ${genotype} -p ${phenotype} -k ${genotype}.kinship.txt.cXX.txt -outdir . -o ${phenotype}_${genotype}.assoc.txt && gzip ${phenotype}_${genotype}.*.assoc.txt



 

