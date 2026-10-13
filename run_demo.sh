#!/bin/bash
# Runs the three scripts on the demo data. Demo-scale overrides (tiny genomes, 6x reads):
#   HARD_MIN/REC_MIN=2 (6x is too low for 4), CHUNK_ROWS=400000 (-> 3 chunks),
#   MIN_SEQ_LEN=10000 (drops the 8 kb contig), WIN=10000 (10 kb windows on 200-250 kb chr).
# Needs kmtricks, gemma, bowtie2, samtools, python3 (+pandas, matplotlib) in PATH,
# or set KMTRICKS= / GEMMA= as for the real scripts.
set -euo pipefail
cd "$(dirname "$0")"
[ -s input.txt ] || python3 make_demo_data.py
export THREADS=${THREADS:-4}
HARD_MIN=2 REC_MIN=2 CHUNK_ROWS=400000 bash kmergwas_1-4.sh input.txt pheno_apo.txt
MIN_SEQ_LEN=10000 WIN=10000 OUT=post bash kmergwas_5-7.sh kmergwas_k31/pheno_apo_k31_apo.assoc.txt.gz genA.fa genB_v2.fa.gz
MIN_LOGP=2.5 MIN_COUNT=3 PER_CHR=genA_chr2,genB-v2_chr1 bash kmergwas_plot.sh post/pheno_apo_k31_apo.LRT0-01.genA+genB-v2.windows10kb.tsv.gz
