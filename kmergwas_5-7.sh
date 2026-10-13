#!/bin/bash
#SBATCH -J kgwas5-7 --mem 50G -c 4 -p ei-medium
#
# kmerGWAS steps 5-7, ONE phenotype per run (one GEMMA output from kmergwas_1-4.sh)
#   5  keep k-mers with p_lrt (col 14) < PMAX                      -> <base>.LRT0-01.txt.gz
#   6  bowtie2 every k-mer (read name = k-mer, so no SEQ/strand join problems)
#        MERGE_REFS=1: refs (seqs >= MIN_SEQ_LEN) merged in ONE index, chr renamed <ref>_<chr> -> one BAM
#        MERGE_REFS=0: one index and one BAM per reference
#   7  python counter (replaces sort/join/bedtools): per WIN window, number of k-mer hits at each
#      -log10(p) rounded to 0.1                                   -> <...>.windows100kb.tsv.gz
#      (chr, start, stop, -log10p, count  = the table the plots use) + per-hit table + chr lengths
#
# usage  sbatch kmergwas_5-7.sh <pheno>.assoc.txt.gz ref1.fa[.gz] [ref2.fa ...]
#        a reference can be given as NAME=path.fa to choose its chr prefix (default: file name;
#        in merged mode '_' in the prefix becomes '-' so plots can split genomes at the first '_')
# resume resubmit the same command. Indexes go to IDX_DIR and are shared (and locked) between phenotypes.
# k-mer length (31/51) is taken from the k-mers themselves.

# source bowtie2-2.4.1_CBG; source samtools-1.11_CBG     # uncomment / adapt to your module system
set -euo pipefail

PMAX=${PMAX:-0.01}
PCOL=${PCOL:-14}                        # 13 p_wald, 14 p_lrt, 15 p_score (GEMMA -lmm 4)
WIN=${WIN:-100000}
MERGE_REFS=${MERGE_REFS:-1}
MIN_SEQ_LEN=${MIN_SEQ_LEN:-1000000}     # 0 keeps every sequence
BT2_OPTS=${BT2_OPTS:-"-a --end-to-end --no-unal"}
THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-4}}
IDX_DIR=${IDX_DIR:-bt2_index}
OUT=${OUT:-.}

[ $# -ge 2 ] || { sed -n '3,20p' "$0"; exit 1; }
IN=$1; shift
log(){ echo "[$(date '+%F %T')] $*" >&2; }
die(){ log "ERROR: $*"; exit 1; }
need(){ command -v "$1" >/dev/null 2>&1 || die "'$1' not found in PATH"; }
step(){ local m=$1; shift; [ -e "$m" ] && return 0; log "run: $*"; "$@"; touch "$m"; }
for b in bowtie2 bowtie2-build samtools python3; do need $b; done
[ -s "$IN" ] || die "no such file: $IN"

case $PCOL in 13) PT=WALD;; 14) PT=LRT;; 15) PT=SCORE;; *) PT=col$PCOL;; esac
WLAB=$WIN; (( WIN % 1000 )) || WLAB=$((WIN/1000))kb
BASE=$(basename "$IN"); BASE=${BASE%.gz}; BASE=${BASE%.txt}; BASE=${BASE%.assoc}
mkdir -p "$OUT" "$IDX_DIR"
FILT=$OUT/$BASE.$PT${PMAX//./-}.txt.gz

TAGS=(); FAS=()
for r in "$@"; do
  if [[ $r == *=* ]]; then t=${r%%=*}; f=${r#*=}; else f=$r; t=$(basename "$f"); t=${t%.gz}; t=${t%.*}; fi
  [ -s "$f" ] || die "reference not found: $f"
  if [ "$MERGE_REFS" = 1 ] && [ $# -gt 1 ]; then t=${t//_/-}; fi
  TAGS+=("$t"); FAS+=("$(readlink -f "$f")")
done

# ---- python helpers ----
PY_PREP=$(cat <<'PY'
import sys, gzip
out, minlen, rename = sys.argv[1], int(sys.argv[2]), sys.argv[3] == '1'
seen, kept = set(), 0
with open(out, 'w') as o, open(out + '.chroms.tsv', 'w') as m:
    for pair in sys.argv[4:]:
        tag, path = pair.split('=', 1)
        recs, name, buf = [], None, []
        def flush():
            global kept
            L = sum(len(x) for x in buf)
            if name is None or L < minlen: return
            new = f"{tag}_{name}" if rename else name
            if new in seen: sys.exit(f"duplicated sequence name {new}")
            seen.add(new); kept += 1
            o.write(f">{new}\n"); o.writelines(x + '\n' for x in buf)
            m.write(f"{new}\t{tag}\t{name}\t{L}\n")
        with (gzip.open if path.endswith('.gz') else open)(path, 'rt') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    flush(); name, buf = line[1:].split()[0], []
                elif line: buf.append(line)
            flush()
if not kept: sys.exit(f"no sequence >= {minlen} bp")
PY
)
PY_COUNT=$(cat <<'PY'
import sys, gzip, math
filt, win, pcol, out_w, out_h, out_g = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]) - 1, *sys.argv[4:7]
pv = {}
with gzip.open(filt, 'rt') as f:                       # k-mer -> p  (the old 'join -1 10 -2 2')
    for l in f:
        t = l.split(); pv[t[1]] = float(t[pcol])
chrs, ln, cnt, nhit, nmiss = [], {}, {}, 0, 0
with gzip.open(out_h, 'wt', compresslevel=3) as h:
    h.write("chr\tstart\tend\tkmer\tstrand\tp\tlog10p\n")
    for l in sys.stdin:                                # samtools view -h
        if l[0] == '@':
            if l.startswith('@SQ'):
                d = dict(x.split(':', 1) for x in l.rstrip('\n').split('\t')[1:])
                chrs.append(d['SN']); ln[d['SN']] = int(d['LN'])
            continue
        q, flag, c, pos = l.split('\t', 4)[:4]
        flag = int(flag)
        if flag & 4: continue
        p = pv.get(q)
        if p is None: nmiss += 1; continue
        s = int(pos) - 1
        lp = f"{-math.log10(max(p, 1e-300)):.1f}"      # the old awk printf "%.1f"
        k = (c, s // win, lp); cnt[k] = cnt.get(k, 0) + 1; nhit += 1
        h.write(f"{c}\t{s}\t{s+len(q)}\t{q}\t{'-' if flag & 16 else '+'}\t{p:.6g}\t{lp}\n")
order = {c: i for i, c in enumerate(chrs)}
with gzip.open(out_w, 'wt') as o:                      # chr, start, stop, -log10p, count
    for (c, w, lp), n in sorted(cnt.items(), key=lambda x: (order[x[0][0]], x[0][1], float(x[0][2]))):
        o.write(f"{c}\t{w*win}\t{min((w+1)*win, ln[c])}\t{lp}\t{n}\n")
with open(out_g, 'w') as g:
    g.writelines(f"{c}\t{ln[c]}\n" for c in chrs)
print(f"count: {len(pv)} k-mers, {nhit} hits, {len(cnt)} window/p rows, {nmiss} hits without p", file=sys.stderr)
PY
)

# ---- steps ----
s5_filter(){
  zcat "$IN" | awk -v c="$PCOL" -v p="$PMAX" '$1!="chr" && $c ~ /^[0-9]/ && $c+0 < p' | gzip > "$FILT.tmp"
  mv -f "$FILT.tmp" "$FILT"
  log "$(zcat "$FILT" | wc -l) k-mers with $PT p < $PMAX"
}
s6_index(){   # $1 index name, $2 rename(0/1), rest NAME=fasta
  local x=$IDX_DIR/$1/$1 rn=$2; shift 2
  python3 -c "$PY_PREP" "$x.fa" "$MIN_SEQ_LEN" "$rn" "$@"
  bowtie2-build --threads "$THREADS" "$x.fa" "$x" > "$x.bowtie2-build.log" 2>&1
}
index(){      # build once, lock so parallel phenotype jobs wait instead of racing
  local d=$IDX_DIR/$1
  mkdir -p "$d"
  until mkdir "$d/.lock" 2>/dev/null; do
    [ -e "$d/$1.done" ] && return 0
    log "waiting for $d/.lock (another job is building it; delete the lock if that job died)"; sleep 120
  done
  trap "rmdir '$d/.lock' 2>/dev/null" EXIT; trap 'exit 143' TERM INT
  step "$d/$1.done" s6_index "$@"
  rmdir "$d/.lock"; trap - EXIT TERM INT
}
s6_map(){     # $1 index name -> $2.bam (sorted + indexed)
  local x=$IDX_DIR/$1/$1
  zcat "$FILT" | awk '{print ">"$2"\n"$2}' \
    | bowtie2 $BT2_OPTS -p "$THREADS" -f -x "$x" -U - 2> "$2.bowtie2.log" \
    | samtools sort -@ "$THREADS" -T "$2.tmp" -o "$2.tmp.bam" -
  mv -f "$2.tmp.bam" "$2.bam"; samtools index -c "$2.bam"
}
s7_count(){   # $1 prefix
  samtools view -h "$1.bam" | python3 -c "$PY_COUNT" "$FILT" "$WIN" "$PCOL" \
    "$1.windows$WLAB.tsv.gz" "$1.hits.tsv.gz" "$1.genome.tsv"
}

# ---- main ----
step "$FILT.done" s5_filter
if [ "$MERGE_REFS" = 1 ]; then
  SETS=("${REFSET:-$(IFS=+; echo "${TAGS[*]}")}")
  PAIRS=(); for i in "${!TAGS[@]}"; do PAIRS+=("${TAGS[$i]}=${FAS[$i]}"); done
  index "${SETS[0]}" "$([ ${#TAGS[@]} -gt 1 ] && echo 1 || echo 0)" "${PAIRS[@]}"
else
  SETS=("${TAGS[@]}")
  for i in "${!TAGS[@]}"; do index "${TAGS[$i]}" 0 "${TAGS[$i]}=${FAS[$i]}"; done
fi
for s in "${SETS[@]}"; do
  P=$OUT/$BASE.$PT${PMAX//./-}.$s
  step "$P.bam.done"   s6_map "$s" "$P"
  step "$P.count.done" s7_count "$P"
  log "done: $P.windows$WLAB.tsv.gz   (plot: sbatch kmergwas_plot.sh $P.windows$WLAB.tsv.gz)"
done
