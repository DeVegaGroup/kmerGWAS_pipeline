#!/bin/bash
#SBATCH -J kgwas1-4 --mem 200G -c 1 -p ei-medium
#
# kmerGWAS steps 1-4
#   1   kmtricks presence/absence matrix (text)   -> k31.aggregate.txt
#   1b  gzip it                                   -> k31.aggregate.txt.gz
#   2   split into chunks of < 10M k-mers         -> chunks/k31_000.txt.gz ...
#   3   chunk -> BIMBAM (kmer, X, X, 0/2 ...)     -> chunks/k31_000.bimbam.gz ...
#   4   GEMMA, per chunk and per phenotype: kinship -gk 1 (per phenotype, so NA individuals are
#       dropped from the SNP filters) + -lmm 4 (wald, lrt, score). When every chunk of a phenotype
#       is done they are merged into  <phenofile>_k31_<pheno>.assoc.txt.gz  (input of kmergwas_5-7.sh)
#
# usage   sbatch kmergwas_1-4.sh input.txt pheno.txt
#   input.txt  kmtricks fof, one sample per line: "ID: r1.fq.gz; r2.fq.gz"
#   pheno.txt  a) GEMMA style: no header/rownames, one line per input.txt sample (same order), NA = missing
#              b) extended: header + first column = ID (text before ':' in input.txt), one column per
#                 phenotype; tab, comma or space separated; samples may be missing/extra/in any order
# resume  resubmit the same command; finished steps leave a .done file and are skipped
# faster  STOP_AFTER=2 sbatch kmergwas_1-4.sh input.txt pheno.txt            # kmtricks + chunks only
#         sbatch --array=0-<N-1> kmergwas_1-4.sh input.txt pheno.txt         # 1 task = 1 chunk (steps 3+4)
#         (the exact --array range is printed after step 2)
# All settings below can be overridden from the environment, e.g.  KMER=51 sbatch kmergwas_1-4.sh ...

set -euo pipefail

KMER=${KMER:-31}                        # 31 or 51
HARD_MIN=${HARD_MIN:-4}
REC_MIN=${REC_MIN:-4}
CHUNK_ROWS=${CHUNK_ROWS:-9000000}       # k-mers per chunk (must be < 10M)
MAF=${MAF:-0.01}
MIN_IND=${MIN_IND:-10}                  # skip phenotypes with fewer non-NA individuals
THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-1}}
KMTRICKS=${KMTRICKS:-kmtricks}          # e.g. KMTRICKS="singularity exec $HOME/thecastle/singularities/kmtricks_151.img kmtricks"
GEMMA=${GEMMA:-gemma}                   # e.g. GEMMA=$HOME/bin/gemma
OUT=${OUT:-kmergwas_k${KMER}}
STOP_AFTER=${STOP_AFTER:-4}             # 1, 2, 3 or 4
KEEP_CHUNK_ASSOC=${KEEP_CHUNK_ASSOC:-0}  # 1 = keep per-chunk GEMMA outputs after the per-phenotype merge

[ $# -eq 2 ] || { sed -n '3,25p' "$0"; exit 1; }
FOF=$(readlink -f "$1"); PHENO=$(readlink -f "$2")
PB=$(basename "$PHENO"); PB=${PB%.*}    # phenotype file name, used in output names
RUN=$OUT/kmtricks_run; AGG=$OUT/k$KMER.aggregate.txt; CH=$OUT/chunks
PG=$OUT/pheno/$PB.gemma.txt; PC=$OUT/pheno/$PB.columns.tsv
mkdir -p "$CH" "$OUT/pheno" "$OUT/assoc"

log(){ echo "[$(date '+%F %T')] $*" >&2; }
die(){ log "ERROR: $*"; exit 1; }
need(){ command -v "$1" >/dev/null 2>&1 || die "'$1' not found in PATH"; }
step(){ local m=$1; shift; [ -e "$m" ] && return 0; log "run: $*"; "$@"; touch "$m"; }
GZ=gzip; command -v pigz >/dev/null && GZ="pigz -p $THREADS"
(( CHUNK_ROWS < 10000000 )) || die "CHUNK_ROWS must be < 10M"
[ -s "$FOF" ] && [ -s "$PHENO" ] || die "missing input.txt or pheno file"

# ---- phenotype module: any pheno table -> GEMMA pheno file (fof order, no header) + column index ----
pheno_prep(){
  python3 - "$FOF" "$PHENO" "$PG.tmp.$$" "$PC.tmp.$$" <<'PY'
import sys, re
fof, pf, og, oc = sys.argv[1:5]
ids = [l.split(':', 1)[0].strip() for l in open(fof) if l.strip()]
if len(set(ids)) != len(ids): sys.exit("pheno: duplicated IDs in the fof")
NA = {'', 'na', 'nan', '.', '-9', 'null', 'none'}
def num(x):
    if x.lower() in NA: return True
    try: float(x); return True
    except ValueError: return False
lines = [l.rstrip('\r\n') for l in open(pf) if l.strip()]
sep = '\t' if '\t' in lines[0] else (',' if ',' in lines[0] else None)
rows = [[t.strip().strip('"') for t in l.split(sep)] for l in lines]
fix = lambda x: 'NA' if x.lower() in NA else x
if all(num(t) for r in rows for t in r):                      # GEMMA style
    if sum(r[0] in ids for r in rows) > len(ids) / 2:
        sys.exit("pheno: first column looks like sample IDs but there is no header line; add one")
    if len(rows) != len(ids):
        sys.exit(f"pheno: GEMMA-style file has {len(rows)} rows, fof has {len(ids)} samples")
    names = [f"p{j+1}" for j in range(len(rows[0]))]
    mat = [[fix(t) for t in r] for r in rows]
else:                                                         # header + rownames
    head, body = rows[0], rows[1:]
    nc = len(body[0]) - 1
    names, d = head[-nc:], {}
    for r in body:
        if len(r) != nc + 1: sys.exit(f"pheno: row '{r[0]}' has {len(r)} fields, expected {nc+1}")
        if r[0] in d: sys.exit(f"pheno: duplicated ID {r[0]}")
        bad = [t for t in r[1:] if not num(t)]
        if bad: sys.exit(f"pheno: non-numeric value(s) for {r[0]}: {bad[:3]}")
        d[r[0]] = [fix(t) for t in r[1:]]
    miss = [i for i in ids if i not in d]; extra = [i for i in d if i not in set(ids)]
    if len(miss) == len(ids): sys.exit("pheno: no fof ID found in the first column of the phenotype file")
    print(f"pheno: {len(ids)-len(miss)}/{len(ids)} fof samples matched; {len(miss)} without phenotype -> NA "
          f"{miss[:5]}; {len(extra)} phenotype rows not in fof ignored {extra[:5]}", file=sys.stderr)
    mat = [d.get(i, ['NA'] * nc) for i in ids]
names = [re.sub(r'[^A-Za-z0-9.-]+', '-', n).strip('-') or f"p{j+1}" for j, n in enumerate(names)]
names = [n if names.count(n) == 1 else f"{n}-{j+1}" for j, n in enumerate(names)]
open(og, 'w').write(''.join('\t'.join(r) + '\n' for r in mat))
open(oc, 'w').write(''.join(f"{j+1}\t{n}\t{sum(r[j] != 'NA' for r in mat)}\n" for j, n in enumerate(names)))
PY
  for f in "$PG" "$PC"; do   # first run: install; later runs: refuse a silently changed phenotype file
    if [ -e "$f" ] && ! cmp -s "$f" "$f.tmp.$$"; then
      rm -f "$f.tmp.$$"; die "$PHENO changed since the last run: remove $OUT/pheno/$PB.* and $OUT/assoc/*/*$PB.* first"
    fi
    mv -f "$f.tmp.$$" "$f"
  done
}

# ---- step functions ----
s1_kmtricks(){
  if [ -d "$RUN" ]; then log "removing incomplete kmtricks run $RUN"; rm -rf "$RUN"; fi
  $KMTRICKS pipeline --hard-min "$HARD_MIN" --recurrence-min "$REC_MIN" --kmer-size "$KMER" --file "$FOF" \
            --threads "$THREADS" --run-dir "$RUN" --mode kmer:pa:bin
}
s1_aggregate(){ $KMTRICKS aggregate --pa-matrix kmer --format text --threads "$THREADS" --run-dir "$RUN" > "$AGG"; }
s1_gzip(){ if [ -f "$AGG" ]; then $GZ -f "$AGG"; fi; gzip -t "$AGG.gz"; }
s2_split(){
  [ -z "$(ls -A "$OUT/assoc")" ] || die "re-splitting would invalidate $OUT/assoc: remove it first"
  rm -f "$CH"/k"${KMER}"_*
  zcat "$AGG.gz" | split -d -a 4 -l "$CHUNK_ROWS" --additional-suffix=.txt --filter="$GZ > \$FILE.gz" - "$CH/k${KMER}_"
  ls "$CH" | sed -n 's/\.txt\.gz$//p' | sort > "$CH/chunks.list"
}
s3_bimbam(){  # "KMER 0 1 1 0" -> "KMER, X, X, 0, 2, 2, 0"  (presence = dosage 2; columns = fof order)
  zcat "$CH/$1.txt.gz" | sed -e 's/1/2/g' -e 's/ /, /g' -e 's/, /, X, X, /' | $GZ > "$CH/$1.bimbam.gz"
}
s4_gemma(){   # $1 chunk, $2 pheno column, $3 pheno name
  local g=$CH/$1.bimbam.gz od=$OUT/assoc/$1 o=$PB.$3
  mkdir -p "$od"
  if [ ! -e "$od/kin_$o.done" ]; then
    "$GEMMA" -gk 1 -g "$g" -p "$PG" -n "$2" -outdir "$od" -o "kin_$o" > "$od/kin_$o.stdout" 2>&1 \
      || die "GEMMA kinship failed: $od/kin_$o.stdout"
    [ -s "$od/kin_$o.cXX.txt" ] || die "GEMMA kinship gave no output: $od/kin_$o.stdout"
    $GZ -f "$od/kin_$o.cXX.txt"; touch "$od/kin_$o.done"
  fi
  rm -f "$od/$o.assoc.txt"
  "$GEMMA" -lmm 4 -maf "$MAF" -g "$g" -p "$PG" -n "$2" -k "$od/kin_$o.cXX.txt.gz" -outdir "$od" -o "$o" \
    > "$od/$o.stdout" 2>&1 || die "GEMMA lmm failed: $od/$o.stdout"
  [ -s "$od/$o.assoc.txt" ] || die "GEMMA lmm gave no output: $od/$o.stdout"
  $GZ -f "$od/$o.assoc.txt"
}
merge_pheno(){  # all chunks of one phenotype -> one file with one header
  local f=$OUT/${PB}_k${KMER}_$1.assoc.txt.gz c
  [ -e "$f.done" ] && return 0
  for c in "${CHUNKS[@]}"; do [ -e "$OUT/assoc/$c/$PB.$1.done" ] || return 0; done
  if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then rmdir "$f.lock" 2>/dev/null || true; fi   # stale lock
  mkdir "$f.lock" 2>/dev/null || { log "another task is merging $f"; return 0; }
  if [ ! -e "$f.done" ]; then
    log "merging chunks -> $f"
    for c in "${CHUNKS[@]}"; do zcat "$OUT/assoc/$c/$PB.$1.assoc.txt.gz"; done \
      | awk 'NR==1 || $1!="chr"' | $GZ > "$f.tmp"
    gzip -t "$f.tmp"; mv -f "$f.tmp" "$f"; touch "$f.done"
    if [ "$KEEP_CHUNK_ASSOC" != 1 ]; then rm -f "$OUT"/assoc/*/"$PB.$1.assoc.txt.gz"; fi   # .done files stay
  fi
  rmdir "$f.lock"
}

# ---- main ----
need python3
pheno_prep
mapfile -t PHENOS < <(awk -v m="$MIN_IND" '$3>=m{print $1":"$2}' "$PC")
awk -v m="$MIN_IND" '$3<m{print "[skip] phenotype "$2": only "$3" non-NA individuals" > "/dev/stderr"}' "$PC"
log "$PB: ${#PHENOS[@]} phenotype(s) to run: ${PHENOS[*]}"

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
  if [ ! -e "$RUN.done" ] || [ ! -e "$AGG.done" ]; then
    need "${KMTRICKS%% *}"; $KMTRICKS infos >/dev/null 2>&1 || die "kmtricks does not run: $KMTRICKS"
  fi
  step "$RUN.done"      s1_kmtricks
  step "$AGG.done"      s1_aggregate
  step "$AGG.gz.done"   s1_gzip
  [ "$STOP_AFTER" -ge 2 ] || exit 0
  step "$CH/chunks.list.done" s2_split
else
  [ -e "$CH/chunks.list.done" ] || die "array mode needs steps 1-2: run once without --array (STOP_AFTER=2)"
fi
mapfile -t CHUNKS < "$CH/chunks.list"
log "${#CHUNKS[@]} chunk(s); parallel option: sbatch --array=0-$(( ${#CHUNKS[@]} - 1 )) kmergwas_1-4.sh $1 $2"
[ "$STOP_AFTER" -ge 3 ] || exit 0

if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
  [ "$SLURM_ARRAY_TASK_ID" -lt "${#CHUNKS[@]}" ] || die "array index $SLURM_ARRAY_TASK_ID >= ${#CHUNKS[@]} chunks"
  MINE=("${CHUNKS[$SLURM_ARRAY_TASK_ID]}")
else
  MINE=("${CHUNKS[@]}")
fi
[ "$STOP_AFTER" -lt 4 ] || need "$GEMMA"

for c in "${MINE[@]}"; do
  step "$CH/$c.bimbam.gz.done" s3_bimbam "$c"
  rm -f "$CH/$c.txt.gz"                                   # raw chunk no longer needed (all in $AGG.gz)
  [ "$STOP_AFTER" -ge 4 ] || continue
  ns=$(zcat "$CH/$c.bimbam.gz" 2>/dev/null | awk -F', ' 'NR==1{print NF-3; exit}' || true)
  [ "$ns" = "$(wc -l < "$PG")" ] || die "$c has $ns samples but the phenotype file has $(wc -l < "$PG") rows"
  for p in "${PHENOS[@]}"; do
    step "$OUT/assoc/$c/$PB.${p#*:}.done" s4_gemma "$c" "${p%%:*}" "${p#*:}"
  done
done
[ "$STOP_AFTER" -ge 4 ] || exit 0
for p in "${PHENOS[@]}"; do merge_pheno "${p#*:}"; done

n_done=$(ls "$OUT"/"${PB}"_k"${KMER}"_*.assoc.txt.gz.done 2>/dev/null | wc -l || true)
log "so far $n_done/${#PHENOS[@]} phenotypes merged into $OUT/${PB}_k${KMER}_<pheno>.assoc.txt.gz"
log "next: sbatch kmergwas_5-7.sh $OUT/${PB}_k${KMER}_<pheno>.assoc.txt.gz ref1.fa [ref2.fa ...]"
