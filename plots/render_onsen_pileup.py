#!/usr/bin/env python3
"""IGV-style read-pileup figure for full-length retrotransposon loci.

Designed to be run *after* `nanoaleseq_pipeline.sh` against the pipeline's
output directory. By default the eight Arabidopsis ATCOPIA78 (ONSEN)
full-length loci are visualised, but the locus list is configurable.

Usage:
    render_onsen_pileup.py \
        --out-dir   /path/to/pipeline/out_dir/ \
        --te-gff    /path/to/TAIR10_TEs.gff \
        --output    figure2_pileup.png

Each panel shows, for one locus:
    top    : per-base coverage profile across the displayed window
    middle : stacked primary alignments (one bar per read; coloured by strand)
    bottom : TAIR10 transposable_element annotation track, focal locus orange
"""
from __future__ import annotations
import argparse, os, re
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches
import numpy as np
import pysam

# ---------------------------------------------------------------------------
DEFAULT_LOCI = [
    "AT1TE59755", "AT1TE12295", "AT3TE92525", "AT5TE15240",
    "AT1TE71045", "AT3TE89830", "AT3TE54550", "AT1TE24850",
]
DEFAULT_FLANK_BP = 1500
DEFAULT_MAX_DRAW_READS = 400

def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out-dir", required=True,
                   help="Pipeline output directory containing aligned/all.bam")
    p.add_argument("--te-gff",  required=True,
                   help="TE annotation GFF (chrom names must match the BAM)")
    p.add_argument("--output",  default="figure_pileup.png",
                   help="Output PNG path (default: figure_pileup.png)")
    p.add_argument("--loci",    default=",".join(DEFAULT_LOCI),
                   help="Comma-separated list of TE feature IDs to render "
                        "(default: 8 Arabidopsis ATCOPIA78 FL loci)")
    p.add_argument("--flank-bp", type=int, default=DEFAULT_FLANK_BP,
                   help=f"Flank around each TE body to display (default: {DEFAULT_FLANK_BP})")
    p.add_argument("--max-draw-reads", type=int, default=DEFAULT_MAX_DRAW_READS,
                   help="Cap on stacked-read drawing per panel (default: 400)")
    p.add_argument("--bam",     default=None,
                   help="Override BAM path (default: <out-dir>/aligned/all.bam)")
    p.add_argument("--feature-type", default="transposable_element",
                   help="GFF feature type to match (column 3). TAIR10 uses "
                        "'transposable_element'; other annotations may use "
                        "'LTR_retrotransposon', 'transposon_fragment', etc.")
    p.add_argument("--title",   default="IGV-style nanoALE-seq read pileup",
                   help="Figure title")
    return p.parse_args()

def parse_gff(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        p = line.rstrip().split("\t")
        if len(p) < 9: continue
        attr = p[8]
        m_id    = re.search(r"ID=([^;]+)", attr)
        m_alias = re.search(r"Alias=([^;]+)", attr)
        rows.append({
            "chrom": p[0], "feature_type": p[2],
            "start": int(p[3]), "end": int(p[4]), "strand": p[6],
            "id":    m_id.group(1)    if m_id    else "",
            "alias": m_alias.group(1) if m_alias else "",
        })
    return rows

def find_locus(rows, locus_id, feature_type):
    for r in rows:
        if r["id"] == locus_id and r["feature_type"] == feature_type:
            return r
    return None

def pack_reads(reads):
    rows_end = []
    out = []
    for r in reads:
        placed = False
        for i, end in enumerate(rows_end):
            if r["start"] >= end + 50:
                rows_end[i] = r["end"]
                out.append((i, r))
                placed = True
                break
        if not placed:
            out.append((len(rows_end), r))
            rows_end.append(r["end"])
    return out, len(rows_end)

# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    out_dir = Path(args.out_dir).expanduser().resolve()
    bam_path = Path(args.bam) if args.bam else out_dir / "aligned" / "all.bam"
    if not bam_path.exists():
        raise SystemExit(f"BAM not found: {bam_path}")

    loci_ids = [s.strip() for s in args.loci.split(",") if s.strip()]
    gff_rows = parse_gff(args.te_gff)
    bam = pysam.AlignmentFile(str(bam_path), "rb")

    # 4-column layout, rows = ceil(n / 4)
    n = len(loci_ids)
    ncols = min(4, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.5*ncols, 5*nrows), squeeze=False)
    plt.subplots_adjust(hspace=0.55, wspace=0.35)

    for idx, locus_id in enumerate(loci_ids):
        ax = axes[idx // ncols][idx % ncols]
        locus = find_locus(gff_rows, locus_id, args.feature_type)
        if locus is None:
            ax.text(0.5, 0.5, f"{locus_id}\n(not in GFF)", ha="center", va="center")
            ax.axis("off"); continue
        chrom = locus["chrom"]
        te_start, te_end = locus["start"], locus["end"]
        rs, re_ = max(1, te_start - args.flank_bp), te_end + args.flank_bp
        region_len = re_ - rs

        reads = []
        for read in bam.fetch(chrom, rs, re_):
            if read.is_unmapped or read.is_secondary or read.is_supplementary: continue
            if read.reference_end is None: continue
            reads.append({"start": read.reference_start, "end": read.reference_end,
                          "strand": "-" if read.is_reverse else "+"})
        n_reads_window = len(reads)
        n_reads_on_te  = sum(1 for r in reads if r["end"] > te_start and r["start"] < te_end)
        reads.sort(key=lambda r: (r["start"], r["end"]))
        if len(reads) > args.max_draw_reads:
            sub = np.linspace(0, len(reads)-1, args.max_draw_reads).astype(int)
            draw_reads = [reads[i] for i in sub]
        else:
            draw_reads = reads
        placed, nrows_used = pack_reads(draw_reads)

        cov = np.zeros(region_len, dtype=int)
        for r in reads:
            s = max(0, r["start"] - rs); e = min(region_len, r["end"] - rs)
            if e > s: cov[s:e] += 1
        x_cov = np.arange(region_len) + rs

        cov_top    = nrows_used + 8
        cov_height = max(6, nrows_used * 0.6)
        cov_norm   = cov / max(cov.max(), 1) * cov_height
        ax.fill_between(x_cov, cov_top, cov_top + cov_norm, color="#446994", alpha=0.85, linewidth=0)
        ax.plot([rs, re_], [cov_top, cov_top], color="#888", linewidth=0.5)
        ax.text(rs, cov_top + cov_height + 0.3, f"max cov = {cov.max()}",
                fontsize=7, color="#446994")

        for row_idx, r in placed:
            color = "#1f3a93" if r["strand"] == "+" else "#c0392b"
            ax.add_patch(Rectangle((r["start"], row_idx + 0.15),
                                   r["end"] - r["start"], 0.7,
                                   facecolor=color, edgecolor="none"))

        annot_y, annot_h = -2.5, 1.2
        for r in gff_rows:
            if r["chrom"] != chrom: continue
            if r["feature_type"] != args.feature_type: continue
            if r["end"] < rs or r["start"] > re_: continue
            is_target = (r["id"] == locus_id)
            face = "#d97a3b" if is_target else "#cccccc"
            edge = "#a25212" if is_target else "#888"
            ax.add_patch(Rectangle((r["start"], annot_y),
                                   r["end"] - r["start"], annot_h,
                                   facecolor=face, edgecolor=edge, linewidth=0.6))
            if is_target:
                ax.text((r["start"]+r["end"])/2, annot_y - 0.6,
                        f"{r['alias']}  {(r['end']-r['start']+1)/1000:.1f} kb",
                        ha="center", va="top", fontsize=7.5, color="#a25212")

        ax.set_title(f"{locus_id}  ({chrom}:{te_start:,}–{te_end:,})\n"
                     f"{n_reads_on_te} primary reads on TE  ·  "
                     f"{n_reads_window} in window (±{args.flank_bp} bp)",
                     fontsize=9.5, loc="left", color="#222")
        ax.set_xlim(rs, re_)
        ax.set_ylim(annot_y - 2, cov_top + cov_height + 1.5)
        ax.set_yticks([])
        xt = np.linspace(rs, re_, 5)
        ax.set_xticks(xt)
        ax.set_xticklabels([f"{x/1000:.1f} kb" for x in xt], fontsize=8)
        ax.spines[["top","right","left"]].set_visible(False)
        ax.tick_params(axis="x", labelsize=8, length=2)

    # Hide unused panels
    for idx in range(n, nrows*ncols):
        axes[idx // ncols][idx % ncols].axis("off")

    fig.legend(handles=[
        mpatches.Patch(facecolor="#1f3a93", label="Read on + strand"),
        mpatches.Patch(facecolor="#c0392b", label="Read on − strand"),
        mpatches.Patch(facecolor="#446994", label="Coverage profile"),
        mpatches.Patch(facecolor="#d97a3b", label="Focal TE locus"),
        mpatches.Patch(facecolor="#cccccc", label="Other TE in window"),
    ], loc="lower center", ncol=5, frameon=False, fontsize=10, bbox_to_anchor=(0.5, -0.01))

    fig.suptitle(args.title, fontsize=13, x=0.04, ha="left", y=1.0)
    plt.savefig(args.output, dpi=300, bbox_inches="tight")
    print(f"Wrote {args.output}")

if __name__ == "__main__":
    main()
