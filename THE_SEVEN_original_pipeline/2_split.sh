#!/bin/bash
#SBATCH --mem 200G -c1 -p ei-medium

split -d -n l/10 k31.aggregate.txt.bimbam k31_bimbam_

