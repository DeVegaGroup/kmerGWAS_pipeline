# kmerGWAS demo: inputs and outputs of each script

Everything in `example_files/` was produced by running `run_demo.sh` with the three scripts in this folder. Files ending in `.head` are the first lines of files that are too big to ship; everything else is complete.

The demo data are simulated by `make_demo_data.py`, which is deterministic and takes about a minute.

- **References:** two genomes, `genA` (chr1–chr3 plus an 8 kb contig) and `genB_v2` (chr1–chr2).
- **Samples:** 48, carrying 60 presence/absence segments of 800 bp.
- **Reads:** 6x single-end 100 bp.
- **Causal segments:** phenotype `apo` is caused by genA chr2:119,000–119,800 and genB chr1:110,000–110,800. Those are the two peaks you should see.

To regenerate and rerun everything (about 8 minutes on 4 cores, 120 MB of reads):

```bash
bash run_demo.sh        # needs kmtricks, gemma, bowtie2, samtools, python3 + pandas/matplotlib in PATH
```

`run_demo.sh` uses demo-scale settings. They are all environment overrides, and the scripts' defaults are for real data:

| Setting | Demo value | Why |
|---|---|---|
| `HARD_MIN`, `REC_MIN` | 2 | 6x coverage is too low for the default 4 |
| `CHUNK_ROWS` | 400000 | makes 3 chunks |
| `MIN_SEQ_LEN` | 10000 | drops the 8 kb contig |
| `WIN` | 10000 | 10 kb windows suit 200–250 kb chromosomes |

---

## Script 1: `kmergwas_1-4.sh` (kmtricks → chunks → BIMBAM → GEMMA)

```bash
HARD_MIN=2 REC_MIN=2 CHUNK_ROWS=400000 bash kmergwas_1-4.sh input.txt pheno_apo.txt
# real data:  sbatch kmergwas_1-4.sh input.txt pheno_apo.txt
```

### Inputs

**`example_files/inputs/input.txt`** is the kmtricks file-of-files: `ID: file; file`, one sample per line. The ID is what the phenotype table must use.

```
16000: /path/to/kmergwas_demo/reads/16000_0.fq.gz
16001: /path/to/kmergwas_demo/reads/16001_0.fq.gz; /path/to/kmergwas_demo/reads/16001_1.fq.gz
16002: /path/to/kmergwas_demo/reads/16002_0.fq.gz
```

**`example_files/inputs/pheno_apo.txt`** is the *extended* phenotype format. It has a header, IDs in the first column and one column per phenotype, and is tab, comma or space separated.

```
accession	apo	plant height	rare
16031	NA	6.03	1.0
16035	1.997	7.81	1.0
16016	2.174	9.71	1.0
16000	2.124	9.50	1.0
88888	1.000	10.00	1.0
```

This file was built to exercise the phenotype module:
- Rows are shuffled.
- Three fof samples are absent (`16003`, `16010`, `3_SX18`) and become NA.
- Two IDs are not in the fof (`88888`, `99999`) and are ignored.
- `apo` has some NA values.
- The column name `plant height` contains a space and becomes `plant-height` in file names.
- `rare` has only 5 non-NA values, so it is skipped (`MIN_IND=10`).

**`example_files/inputs/pheno_gemmastyle.txt`** holds the same `apo` values in plain GEMMA format: no header, no IDs, one line per fof sample in fof order. The script detects which format it was given.

```
2.124
0.730
0.544
NA
```

### Outputs (in `kmergwas_k31/`; full tree in `outputs_1-4/directory_tree.txt`)

| File | What it is |
|---|---|
| `pheno/pheno_apo.gemma.txt` | phenotype table rewritten for GEMMA: fof order, no header, NA for missing |
| `pheno/pheno_apo.columns.tsv` | column number → phenotype name → non-NA individuals |
| `k31.aggregate.txt.gz` | kmtricks presence/absence matrix (step 1 + gzip step) |
| `chunks/chunks.list`, `chunks/k31_0000.bimbam.gz` … | chunks of < `CHUNK_ROWS` k-mers in BIMBAM (step 2 + 3) |
| `assoc/k31_0000/kin_pheno_apo.apo.cXX.txt.gz` | kinship per chunk and per phenotype |
| `assoc/k31_0000/pheno_apo.apo.log.txt` | GEMMA log (also `.stdout`) |
| **`pheno_apo_k31_apo.assoc.txt.gz`** | **all chunks merged, one file per phenotype → input of script 2** |
| `*.done` | resume markers |

`pheno/pheno_apo.columns.tsv`: `rare` is listed but not run.

```
1	apo	41
2	plant-height	45
3	rare	5
```

`pheno/pheno_apo.gemma.txt` has its first rows in fof order (16000, 16001, 16002, 16003, …). Row 4 is `16003`, which is absent from the table, so all NA:

```
2.124	9.50	1.0
0.730	9.02	NA
0.544	9.87	NA
NA	NA	NA
4.070	8.47	NA
```

`k31.aggregate.txt` (`outputs_1-4/k31.aggregate.txt.head`) has one k-mer per line and one 0/1 per sample in fof order:

```
AAAAAAAAATCCTGACCGGTACGCTGAATAA 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 0 1 1 1 0 ...
AAAAAAAACAATTATGATGGTACAAGAGTTC 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ...
```

`chunks/k31_0000.bimbam.gz` (`outputs_1-4/chunks/k31_0000.bimbam.head`) is GEMMA mean-genotype format, with presence coded as dosage 2:

```
AAAAAAAAATCCTGACCGGTACGCTGAATAA, X, X, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, ...
AAAAAAAACAATTATGATGGTACAAGAGTTC, X, X, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, ...
```

`assoc/k31_0000/kin_pheno_apo.apo.cXX.txt.gz` is a 48 × 48 matrix covering all samples, as GEMMA always outputs (top-left 5 × 5 shown):

```
0.3274340421	-0.000527061833	-0.007083983162	-0.008997102012	-0.01852837926
-0.000527061833	0.3205978365	-0.006824339437	-0.01579251571	-0.001206063627
```

`assoc/k31_0000/pheno_apo.apo.log.txt` contains the exact GEMMA command. 41 individuals were analysed for `apo`: 45 matched minus 4 NA.

```
## Command Line Input = gemma -lmm 4 -maf 0.01 -g kmergwas_k31/chunks/k31_0000.bimbam.gz -p kmergwas_k31/pheno/pheno_apo.gemma.txt -n 1 -k kmergwas_k31/assoc/k31_0000/kin_pheno_apo.apo.cXX.txt.gz -outdir kmergwas_k31/assoc/k31_0000 -o pheno_apo.apo
## number of analyzed individuals = 41
```

**`pheno_apo_k31_apo.assoc.txt.gz`** is the GEMMA `-lmm 4` output with 15 columns and one header. Column 2 is the k-mer and column 14 is `p_lrt`. The shipped copy in `outputs_1-4/` keeps only the rows with `p_lrt < 0.05` to stay small; script 2 gives identical results from it, so you can try script 2 without running kmtricks or GEMMA.

```
chr	rs	ps	n_miss	allele1	allele0	af	beta	se	logl_H1	l_remle	l_mle	p_wald	p_lrt	p_score
-9	AAAAAAAAATCCTGACCGGTACGCTGAATAA	-9	0	X	X	0.902	-4.185550e-01	3.679927e-01	-6.555006e+01	1.000000e+05	1.000000e+05	2.623109e-01	2.473955e-01	2.582368e-01
```

---

## Script 2: `kmergwas_5-7.sh` (filter → bowtie2 → window counts)

```bash
MIN_SEQ_LEN=10000 WIN=10000 OUT=post bash kmergwas_5-7.sh \
    kmergwas_k31/pheno_apo_k31_apo.assoc.txt.gz  genA.fa  genB_v2.fa.gz
# real data:  sbatch kmergwas_5-7.sh kmergwas_k31/pheno_apo_k31_apo.assoc.txt.gz AAmalSCBG=malSCBG.fa BB=BB.fa SS=SS.fa
```

### Inputs

The inputs are the GEMMA file from script 1 plus one or more references (plain or `.gz`). `example_files/inputs/genA.fa.head` shows that only the first word of a header is kept:

```
>chr1 some description
GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCG
```

With several references, the default `MERGE_REFS=1` merges them into one index. Only sequences of at least `MIN_SEQ_LEN` are kept, and they are renamed `<prefix>_<chr>`. The prefix is the file name, with `_` turned into `-` (`genB_v2.fa.gz` → `genB-v2`), or the name you give with `NAME=file.fa`.

### Outputs (in `post/`; all in `example_files/outputs_5-7/`)

| File | What it is |
|---|---|
| `pheno_apo_k31_apo.LRT0-01.txt.gz` | step 5: GEMMA rows with p_lrt < 0.01 (8,971 k-mers) |
| `bt2_index/genA+genB-v2/genA+genB-v2.fa.chroms.tsv` | merged reference: new name, genome, original name, length (shared index dir) |
| `pheno_apo_k31_apo.LRT0-01.genA+genB-v2.bam` + `.csi` | step 6: sorted, indexed BAM; read name = k-mer |
| `….bowtie2.log` | bowtie2 alignment summary |
| **`….windows10kb.tsv.gz`** | **step 7: chr, start, stop, −log10(p) rounded to 0.1, number of k-mer hits → input of the plot script** |
| `….hits.tsv.gz` | step 7: every hit with its p-value (useful to pull the k-mers under a peak) |
| `….genome.tsv` | chromosome lengths, used by the plot script |

`genA+genB-v2.fa.chroms.tsv`: `ctgSmall` (8 kb) was dropped.

```
genA_chr1	genA	chr1	250000
genA_chr2	genA	chr2	250000
genA_chr3	genA	chr3	200000
genB-v2_chr1	genB-v2	chr1	200000
genB-v2_chr2	genB-v2	chr2	150000
```

BAM records (`samtools view`). The read name (column 1) is the k-mer. In the second record, flag 16 means a reverse-strand hit, and SEQ (column 10) is reverse-complemented. This is why the old join on column 10 lost these hits, while joining on the read name keeps them:

```
CGGTAGACCAGAGGTCATCAGCTTCGGTAGC	0	genA_chr2	119201	255	31M	*	0	0	CGGTAGACCAGAGGTCATCAGCTTCGGTAGC	...
GGCTACCGAAGCTGATGACCTCTGGTCTACC	16	genA_chr2	119202	255	31M	*	0	0	GGTAGACCAGAGGTCATCAGCTTCGGTAGCC	...
```

**`….windows10kb.tsv.gz`** has no header, the same layout as your old `100kbwindows.txt`, and 1,003 rows. The first rows, then the peak window on genA chr2:

```
genA_chr1	0	10000	2.0	8
genA_chr1	0	10000	2.1	8
genA_chr1	0	10000	2.2	22
...
genA_chr2	110000	120000	7.0	3
genA_chr2	110000	120000	7.2	37
genA_chr2	110000	120000	7.3	3
genA_chr2	110000	120000	7.4	7
genA_chr2	110000	120000	7.5	9
genA_chr2	110000	120000	8.3	2
```

`….hits.tsv.gz` has a header and 8,865 hits. The best hits are the causal segment on genA chr2:

```
chr	start	end	kmer	strand	p	log10p
genA_chr2	119237	119268	AGTCTACGGTAGAATAGTCAGCAGGGCCACA	+	5.32435e-09	8.3
genA_chr2	119238	119269	TTGTGGCCCTGCTGACTATTCTACCGTAGAC	-	5.32435e-09	8.3
genA_chr2	119230	119261	CCTGGATAGTCTACGGTAGAATAGTCAGCAG	+	3.02258e-08	7.5
```

---

## Script 3: `kmergwas_plot.sh`

```bash
MIN_LOGP=2.5 MIN_COUNT=3 PER_CHR=genA_chr2,genB-v2_chr1 \
    bash kmergwas_plot.sh post/pheno_apo_k31_apo.LRT0-01.genA+genB-v2.windows10kb.tsv.gz
# real data defaults: MIN_LOGP=2.5 MIN_COUNT=10 (your "min10" plots), no single-chromosome plots
```

### Input

The input is the windows table from script 2. Old bedtools-era tables also work; their repeated window/p rows are summed. If a matching `….genome.tsv` sits next to the table, it supplies chromosome lengths.

### Outputs (in `example_files/outputs_plot/`; PDF plus a PNG preview of each)

Chromosome names contain `_`, so the plots are split back per genome:

| File | Plot |
|---|---|
| `….genA.manhattan.pdf`, `….genB-v2.manhattan.pdf` | whole genome, cumulative position, alternating blues, size = k-mer count |
| `….genA.chromosomes.pdf`, `….genB-v2.chromosomes.pdf` | one panel per chromosome, x in Mb, colour = −log10(p) |
| `….genA.chr2.pdf`, `….genB-v2.chr1.pdf` | single chromosomes requested with `PER_CHR` |

The log also lists the top windows per genome, for example:

```
  genA: plotted 168 rows; top windows:
      chr  start   stop  logp  count genome
genA_chr2 110000 120000   7.5      9   genA
genA_chr2 110000 120000   7.4      7   genA
  genB-v2: plotted 105 rows; top windows:
genB-v2_chr1 110000 120000   6.3      7 genB-v2
genB-v2_chr1 110000 120000   6.1     31 genB-v2
```
