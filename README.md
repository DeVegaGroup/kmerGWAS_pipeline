# kmerGWAS pipeline
kmer GWAS pipeline (from reads to Manhattan).

JDV's original Kmer GWAS pipeline has been reimplemented as three scripts. The scripts have been commented on by AI, and those comments have been checked for accuracy. The prototype, which grew as projects were reused into "the-seven" scripts used in several projects, it is included in its own folder. Another folder includes plotting scripts and examples of plot outputs. 

Note: Jamie has a reimplementation of the-seven as a Snakemake pipeline (email me).


## Inputs and outputs. 

The first script takes as input a table of paths to read files for each accession/sample (one or multiple, as fq.gz) using the kmtricks format for the input txt table, as well as one phenotype file (with a header and rownames matching those in the first column of the paths file). It can also take a single phenotype file in GEMMA format (no header, no rownames, with samples in some order). It outputs a file of kmers associated with the phenotype, with adjusted p-values in GEMMA format.

The second script takes the GEMMA output (one per phenotype in the input) and one or more genomes in FASTA format, and produces a table of genomic regions (100 Kbp by default) and the significant kmers in each (grouped by adjusted p-value, in bins of one decimal point), one per phenotype. 

The third script does the plotting.

Tested with kmtricks 1.6, GEMMA 0.98.5, bowtie2 and samtools



## How to run

```
# steps 1-4 (single job; resubmit the same command to resume)
KMTRICKS="singularity exec $HOME/thecastle/singularities/kmtricks_151.img kmtricks" GEMMA=$HOME/bin/gemma \
  sbatch kmergwas_1-4.sh input.txt pheno_apo.txt

#   or in parallel: chunks first, then one array task per chunk (the exact range is printed in the log)
STOP_AFTER=2 sbatch kmergwas_1-4.sh input.txt pheno_apo.txt
sbatch --array=0-N kmergwas_1-4.sh input.txt pheno_apo.txt

# steps 5-7, once per phenotype (refs are merged by default; NAME=file.fa sets the chr prefix)
sbatch kmergwas_5-7.sh kmergwas_k31/pheno_apo_k31_apo.assoc.txt.gz AAmalSCBG=malSCBG.fa BB=BB.fa SS=SS.fa

# plots
sbatch kmergwas_plot.sh pheno_apo_k31_apo.LRT0-01.AAmalSCBG+BB+SS.windows100kb.tsv.gz
```

Every setting at the top of each script (KMER=51, CHUNK_ROWS, MAF, PMAX, WIN, MERGE_REFS, MIN_SEQ_LEN, MIN_COUNT, PER_CHR, and so on) can be set on the command line in front of sbatch, as in the examples, so you don’t need to edit the files.

kmergwas_5-7.sh has a commented source bowtie2…; source samtools… line near the top. Uncomment it if those tools aren’t already in your PATH; the script stops with a clear message if they’re missing.


## What's going on in each script

### Script 1 (1–4).

* Phenotype check first: it matches IDs (the text before “:” in input.txt) before kmtricks starts, so ID mismatches fail. It reorders rows to the input order, fills missing samples with NA, and skips phenotypes with fewer than 10 non-NA individuals.
* Recovery: every step leaves a .done file.
* Chunking: default 9M k-mers per chunk, enforced below 10M.
* GEMMA: kinship plus -lmm 4 runs per chunk and per phenotype. Chunks are merged into one file per phenotype, which is the input to script 2.
* Parallel option: run once with STOP_AFTER=2, then submit a sbatch --array; each task converts and runs one chunk. A lock stops two tasks merging at the same time.
* Safety guard: it refuses to continue if the phenotype file changed since the last run.

### Script 2 (5–7).

* Filter: keeps k-mers with p_lrt (column 14) below 0.01.
* References: merged by default as in the banana paper. Sequences of 1 Mb or more are kept and renamed <ref>_<chr>; MERGE_REFS=0 gives one BAM per reference instead.
* a short Python helper does the join, windowing and counting. Output is chr, start, stop, -log10p, count, plus a per-hit table and chromosome lengths.
* Shared index: the bowtie2 index is locked, so phenotypes submitted in parallel wait for it rather than both building it.

### Plot script. 
* cumulative Manhattan in alternating blues, and per-chromosome panels with the orange-to-red gradient, sized by k-mer count. 
* It splits merged genomes back apart at the first _, and can draw single chromosomes with PER_CHR=

### Things to know when you run it:

* Singularity kmtricks: set the path with $HOME, not ~
* Lots of Disk: the merged GEMMA output is roughly 5–8 GB per phenotype. Per-chunk outputs are deleted by default after a verified merge unless you set KEEP_CHUNK_ASSOC=1.
* Point sizes: with correct counts, window counts reach the thousands, so the point-size scale is capped at the 99th percentile to keep it readable.