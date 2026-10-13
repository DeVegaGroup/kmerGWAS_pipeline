#!/bin/bash
#SBATCH -J kgwas-plot --mem 16G -c 1 -p ei-medium
#
# kmerGWAS plots from the step-7 window table (chr, start, stop, -log10p, count; plain or .gz,
# also the old bedtools tables: repeated window/p rows are summed).
#   <prefix>[.<genome>].manhattan.pdf    whole genome, cumulative position, size = k-mer count
#   <prefix>[.<genome>].chromosomes.pdf  one panel per chromosome (x in Mb, colour = -log10p)
#   <prefix>[.<genome>].<chr>.pdf        single chromosomes (PER_CHR=all or PER_CHR=chr02,chr07)
# Merged references (chr = <genome>_<chr>, e.g. AAmalSCBG_chr1) are split back into one set of
# plots per genome (SPLIT=1, at the first '_'; only if every chr name has one).
# Chromosome lengths come from <prefix>.genome.tsv when present (written by kmergwas_5-7.sh).
#
# usage  sbatch kmergwas_plot.sh table.windows100kb.tsv.gz [more tables ...]
#        MIN_LOGP=4 MIN_COUNT=3 PER_CHR=chr02 sbatch kmergwas_plot.sh table.tsv.gz

set -euo pipefail
MIN_LOGP=${MIN_LOGP:-2.5}     # keep rows with -log10p >  MIN_LOGP
MIN_COUNT=${MIN_COUNT:-10}    # keep rows with count    >= MIN_COUNT
YMAX=${YMAX:-}                # empty = automatic
SPLIT=${SPLIT:-1}
PER_CHR=${PER_CHR:-}          # '' | all | comma list of chr names
NCOL=${NCOL:-3}
FMT=${FMT:-pdf}               # pdf | png | svg
export MIN_LOGP MIN_COUNT YMAX SPLIT PER_CHR NCOL FMT

[ $# -ge 1 ] || { sed -n '3,16p' "$0"; exit 1; }
command -v python3 >/dev/null || { echo "python3 not in PATH" >&2; exit 1; }
python3 -c 'import pandas, matplotlib' 2>/dev/null || { echo "python3 needs pandas + matplotlib" >&2; exit 1; }

python3 - "$@" <<'PY'
import os, re, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.ticker import MaxNLocator

E = os.environ
MIN_LOGP, MIN_COUNT = float(E["MIN_LOGP"]), float(E["MIN_COUNT"])
YMAX, SPLIT, PER_CHR, NCOL, FMT = E["YMAX"], E["SPLIT"] == "1", E["PER_CHR"], int(E["NCOL"]), E["FMT"]
CMAP = LinearSegmentedColormap.from_list("kgwas", ["#fdd49e", "#990000"])   # R: low/high of scale_colour_gradient
BLUES = ["#276FBF", "#183059"]
natkey = lambda s: [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", str(s))]

def sizer(counts, smin, smax):                     # area proportional to count (like ggplot scale_size),
    cmin, cmax = counts.min(), max(np.quantile(counts, .99), counts.min() + 1)   # capped at the 99th pct
    f = lambda c: smin + (np.clip(np.asarray(c, float), cmin, cmax) - cmin) / (cmax - cmin) * (smax - smin)
    f.lim = (int(cmin), int(cmax)); return f

def size_legend(ax, f, **kw):
    cmin, cmax = f.lim
    vals = [v for v in MaxNLocator(4, integer=True).tick_values(cmin, cmax) if cmin <= v < cmax] + [cmax]
    h = [ax.scatter([], [], s=float(f(v)), c="#333333", alpha=.8) for v in vals]
    labs = [str(int(v)) for v in vals[:-1]] + [f">={cmax}"]
    return ax.legend(h, labs, title="k-mers", frameon=False, **kw)

def ylims(d):
    top = float(YMAX) if YMAX else max(d.logp.max() + .5, MIN_LOGP + 1)
    return MIN_LOGP, top

def manhattan(d, lens, out, title):
    chrs = sorted(lens, key=natkey)
    off = dict(zip(chrs, np.concatenate([[0], np.cumsum([lens[c] for c in chrs])[:-1]])))
    d = d.sort_values("logp")
    x = (d.start + d.stop) / 2 + d.chr.map(off)
    col = d.chr.map({c: BLUES[i % 2] for i, c in enumerate(chrs)})
    f = sizer(d["count"], 3, 160)
    fig, ax = plt.subplots(figsize=(14, 6))
    ax.scatter(x, d.logp, s=f(d["count"]), c=col, alpha=.7, linewidths=0, rasterized=len(d) > 50000)
    ax.set_xticks([off[c] + lens[c] / 2 for c in chrs])
    ax.set_xticklabels([lab(c) for c in chrs], rotation=90 if len(chrs) > 12 else 0, fontsize=8)
    ax.set_xlim(0, sum(lens.values())); ax.set_ylim(*ylims(d))
    ax.set_ylabel("-log10(p)"); ax.set_title(title, fontsize=10, loc="left")
    ax.grid(alpha=.3); [ax.spines[s].set_visible(False) for s in ("top", "right")]
    size_legend(ax, f, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    fig.tight_layout(); fig.savefig(out, dpi=200); plt.close(fig)

def panel(ax, d, f, norm, raster):
    d = d.sort_values("logp")
    return ax.scatter((d.start + d.stop) / 2e6, d.logp, s=f(d["count"]), c=d.logp, cmap=CMAP, norm=norm,
                      alpha=.8, linewidths=0, rasterized=raster)

def chromosomes(d, lens, out, title):
    chrs = sorted(lens, key=natkey); n = len(chrs); nc = min(NCOL, n); nr = int(np.ceil(n / nc))
    norm = Normalize(MIN_LOGP, d.logp.max()); f = sizer(d["count"], 2, 110)
    fig, axs = plt.subplots(nr, nc, figsize=(4.3 * nc + 1.8, 2.3 * nr + .8), sharex=True, sharey=True, squeeze=False)
    for ax, c in zip(axs.flat, chrs):
        sc = panel(ax, d[d.chr == c], f, norm, len(d) > 20000)
        ax.set_title(lab(c), fontsize=9); ax.grid(alpha=.3)
        [ax.spines[s].set_visible(False) for s in ("top", "right")]
    for ax in axs.flat[n:]: ax.set_visible(False)
    for i in range(max(0, n - nc), n): axs.flat[i].tick_params(labelbottom=True)   # x labels on the last row of each column
    axs.flat[0].set_xlim(0, max(lens.values()) / 1e6); axs.flat[0].set_ylim(*ylims(d))
    fig.supxlabel("position (Mb)"); fig.supylabel("-log10(p)"); fig.suptitle(title, fontsize=10, x=.01, ha="left")
    fig.tight_layout(rect=(0, 0, .88, 1))
    fig.colorbar(sc, cax=fig.add_axes([.9, .1, .012, .3]), label="-log10(p)")
    size_legend(fig.add_axes([.89, .5, .1, .4], frameon=False, xticks=[], yticks=[]), f, loc="center left")
    fig.savefig(out, dpi=200); plt.close(fig)

def one_chr(d, c, L, out, title):
    norm = Normalize(MIN_LOGP, d.logp.max()); f = sizer(d["count"], 6, 300)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    sc = panel(ax, d, f, norm, len(d) > 20000)
    ax.set_xlim(0, L / 1e6); ax.set_ylim(*ylims(d)); ax.grid(alpha=.3)
    ax.set_xlabel(f"{lab(c)} position (Mb)"); ax.set_ylabel("-log10(p)"); ax.set_title(title, fontsize=10, loc="left")
    [ax.spines[s].set_visible(False) for s in ("top", "right")]
    fig.colorbar(sc, ax=ax, label="-log10(p)", shrink=.5, anchor=(0, 0))
    size_legend(ax, f, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    fig.tight_layout(); fig.savefig(out, dpi=200); plt.close(fig)

for tab in sys.argv[1:]:
    pre = re.sub(r"\.(tsv|txt)(\.gz)?$", "", tab)
    df = pd.read_csv(tab, sep="\t", header=None, names=["chr", "start", "stop", "logp", "count"],
                     dtype={"chr": str, "start": np.int64, "stop": np.int64, "logp": float, "count": np.int64})
    df = df.groupby(["chr", "start", "stop", "logp"], as_index=False)["count"].sum()
    gfile = re.sub(r"\.windows[^.]*$", "", pre) + ".genome.tsv"
    if os.path.exists(gfile):
        lens = dict(pd.read_csv(gfile, sep="\t", header=None, dtype={0: str}).values)
    else:
        lens = df.groupby("chr").stop.max().to_dict()
    split = SPLIT and all("_" in c for c in lens)
    gen = (lambda c: c.split("_", 1)[0]) if split else (lambda c: "")
    lab = (lambda c: c.split("_", 1)[1]) if split else (lambda c: c)
    df["genome"] = df.chr.map(gen)
    keep = df[(df.logp > MIN_LOGP) & (df["count"] >= MIN_COUNT)]
    print(f"{tab}: {len(df)} rows, {len(keep)} with -log10p > {MIN_LOGP:g} and count >= {MIN_COUNT:g}", file=sys.stderr)
    for g in sorted({gen(c) for c in lens}, key=natkey):
        d = keep[keep.genome == g]; L = {c: l for c, l in lens.items() if gen(c) == g}
        tag = f"{pre}.{g}" if g else pre
        title = f"{os.path.basename(pre)}  {g}   (-log10p > {MIN_LOGP:g}, k-mers per window/p >= {MIN_COUNT:g})"
        if d.empty: print(f"  {g or 'all'}: nothing above the thresholds", file=sys.stderr); continue
        manhattan(d, L, f"{tag}.manhattan.{FMT}", title)
        chromosomes(d, L, f"{tag}.chromosomes.{FMT}", title)
        want = set(L) if PER_CHR == "all" else {c for c in L if c in PER_CHR.split(",") or lab(c) in PER_CHR.split(",")}
        for c in sorted(want, key=natkey):
            if (d.chr == c).any(): one_chr(d[d.chr == c], c, L[c], f"{tag}.{lab(c)}.{FMT}", f"{title}  {lab(c)}")
        top = d.sort_values(["logp", "count"], ascending=False).head(5)
        print(f"  {g or 'all'}: plotted {len(d)} rows; top windows:\n" + top.to_string(index=False), file=sys.stderr)
PY
