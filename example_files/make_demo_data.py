#!/usr/bin/env python3
"""Simulated kmerGWAS demo data (deterministic, ~1 min, ~150 MB of reads).
Two reference genomes (genA: 3 chr + one 8 kb contig, genB: 2 chr), 48 samples with 60
presence/absence segments of 800 bp, 6x single-end 100 bp reads (no errors).
Phenotype 'apo' is caused by two segments: genA chr2:119000-119800 and genB chr1:110000-110800."""
import gzip, os, random
random.seed(7)
rc = str.maketrans("ACGT", "TGCA")
rseq = lambda n: ''.join(random.choice("ACGT") for _ in range(n))
refs = {"genA": {"chr1": rseq(250000), "chr2": rseq(250000), "chr3": rseq(200000), "ctgSmall": rseq(8000)},
        "genB": {"chr1": rseq(200000), "chr2": rseq(150000)}}
for g, fn in (("genA", "genA.fa"), ("genB", "genB_v2.fa.gz")):
    with (gzip.open if fn.endswith("gz") else open)(fn, "wt") as f:
        for c, s in refs[g].items():
            f.write(f">{c} some description\n" + "\n".join(s[i:i+60] for i in range(0, len(s), 60)) + "\n")
segs = []
for g in refs:
    for c in [c for c in refs[g] if c.startswith("chr")]:
        L = len(refs[g][c]); segs += [(g, c, p, p + 800) for p in sorted(random.sample(range(5000, L - 5000, 3000), 12))]
causal = ([i for i, s in enumerate(segs) if s[:2] == ("genA", "chr2")][3],
          [i for i, s in enumerate(segs) if s[:2] == ("genB", "chr1")][6])
print("causal segments:", [segs[i] for i in causal])
os.makedirs("reads", exist_ok=True)
ids = [str(16000 + i) for i in range(44)] + [f"{i}_SX18" for i in range(2, 6)]
fof, ph = [], {}
for sid in ids:
    pres = [random.random() < 0.5 for _ in segs]
    genome = []
    for g in refs:
        for c, s in refs[g].items():
            last, parts = 0, []
            for a, b in [(a, b) for (gg, cc, a, b), p in zip(segs, pres) if gg == g and cc == c and not p]:
                parts.append(s[last:a]); last = b
            genome.append(''.join(parts) + s[last:])
    nf = 2 if sid.endswith("1") else 1                      # some samples have 2 read files
    paths = [os.path.abspath(f"reads/{sid}_{j}.fq.gz") for j in range(nf)]
    fs = [gzip.open(p, "wt", compresslevel=1) for p in paths]
    k = 0
    for s in genome:
        for _ in range(len(s) * 6 // 100):
            p = random.randrange(0, len(s) - 100); r = s[p:p+100]
            if random.random() < .5: r = r.translate(rc)[::-1]
            fs[k % nf].write(f"@r{k}\n{r}\n+\n{'I'*100}\n"); k += 1
    for f in fs: f.close()
    fof.append(f"{sid}: " + "; ".join(paths))
    ph[sid] = (2.0 * pres[causal[0]] + 1.5 * pres[causal[1]] + random.gauss(0, .5), random.gauss(10, 2))
open("input.txt", "w").write("\n".join(fof) + "\n")
# extended phenotype table: header, shuffled rows, 3 fof samples absent, 2 extra IDs, NAs, a mostly-NA column
rows = [s for s in ids if s not in ("16003", "16010", "3_SX18")] + ["99999", "88888"]
random.shuffle(rows)
with open("pheno_apo.txt", "w") as f:
    f.write("accession\tapo\tplant height\trare\n")
    for i, s in enumerate(rows):
        a, h = ph.get(s, (1.0, 10.0))
        f.write(f"{s}\t{'NA' if i % 13 == 0 else f'{a:.3f}'}\t{h:.2f}\t{'NA' if i > 5 else '1.0'}\n")
# the same 'apo' values in plain GEMMA format (no header/IDs, fof order)
d = {l.split('\t')[0]: l.split('\t')[1] for l in open("pheno_apo.txt").read().splitlines()[1:]}
open("pheno_gemmastyle.txt", "w").write(''.join(d.get(s, 'NA') + '\n' for s in ids))
