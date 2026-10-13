#!/bin/bash
#SBATCH --mem 200G -c 1

#singularity exec /ei/software/testing/kmtricks/1.2.1/kmtricks_1_2_1.img kmtricks pipeline --hard-min 2 --recurrence-min 5 --kmer-size 31 --file input_libs.txt --threads 4 --run-dir try1 --mode kmer:pa:bin

#singularity exec /ei/software/testing/kmtricks/1.2.1/kmtricks_1_2_1.img kmtricks aggregate --pa-matrix kmer --format text --run-dir try1 > k31.aggregate.txt



singularity exec ~/thecastle/singularities/kmtricks_151.img kmtricks pipeline --hard-min 4 --recurrence-min 4 --kmer-size 31 --file input.txt --threads 1 --run-dir try1 --mode kmer:pa:bin && singularity exec ~/thecastle/singularities/kmtricks_151.img kmtricks aggregate --pa-matrix kmer --format text --run-dir try1 > k31.aggregate.txt && touch k31.aggregate.txt.done;


#singularity exec /ei/software/testing/kmtricks/1.2.1/kmtricks_1_2_1.img kmtricks pipeline --hard-min 5 --recurrence-min 5 --kmer-size 31 --file input.txt --threads 8 --run-dir try1 --mode kmer:pa:bin

#singularity exec /ei/software/testing/kmtricks/1.2.1/kmtricks_1_2_1.img kmtricks aggregate --pa-matrix kmer --format text --run-dir try1 > k31.aggregate.txt && touch k31.aggregate.txt.done;




