#!/usr/bin/env python3
"""TE superfamily breakdown + per-chromosome density + per-source length distribution.

Designed to be run *after* `nanoaleseq_pipeline.sh`. Reads the pipeline's
strict-cutadapt outputs and produces a 3-panel summary figure.

Defaults are tuned for *Arabidopsis* TAIR10 (NCBI accessions Chr1–5 + Mt + Pt
+ ATCOPIA*/ATGP*/ATHILA* family aliases). For other species, override the
`--nuclear-regex`, `--mt-contig`, `--pt-contig` and `--chrom-labels` CLI flags
accordingly. Panel (a) family-superfamily mapping is a TAIR10 conservative
whitelist; non-matching families fall back to "Other / unclassified".

Usage:
    render_breakdown.py \
        --out-dir   /path/to/pipeline/out_dir/ \
        --te-gff    /path/to/TAIR10_TEs.gff \
        --output    figure3_te_breakdown.png

Three panels:
    (a) Top 10 TE superfamilies by read–TE overlap, colour-coded by class.
        ATCOPIA78 / ONSEN highlighted as a positive-control row.
    (b) Per-chromosome mapped-read density on log scale.
    (c) Per-source read length distribution with mean lines.
"""
from __future__ import annotations
import argparse, json, os, re, subprocess
from pathlib import Path
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# TAIR10 defaults
DEFAULT_NUCLEAR_REGEX = r"^NC_00(3070\.9|3071\.7|3074\.8|3075\.7|3076\.8)$"
DEFAULT_MT_CONTIG     = "NC_037304.1"
DEFAULT_PT_CONTIG     = "NC_000932.1"
DEFAULT_CHROM_LABELS  = {
    "NC_003070.9": "Chr1", "NC_003071.7": "Chr2", "NC_003074.8": "Chr3",
    "NC_003075.7": "Chr4", "NC_003076.8": "Chr5",
    "NC_037304.1": "Mt",   "NC_000932.1": "Pt",
}

# ---------------------------------------------------------------------------
COLOR_LTR_COPIA   = "#1f3a93"
COLOR_LTR_GYPSY   = "#2e7d8c"
COLOR_NON_LTR     = "#d97a3b"
COLOR_DNA         = "#7a5c9e"
COLOR_HELITRON    = "#3aa05a"
COLOR_ONSEN_HL    = "#c0392b"
COLOR_OTHER       = "#888888"

CLASS_OF = {
    "LTR/Copia":    ("LTR retro / Copia",          COLOR_LTR_COPIA),
    "LTR/Gypsy":    ("LTR retro / Gypsy",          COLOR_LTR_GYPSY),
    "non-LTR/LINE": ("non-LTR retro / LINE",       COLOR_NON_LTR),
    "non-LTR/SINE": ("non-LTR retro / SINE",       COLOR_NON_LTR),
    "DNA/MULE":     ("DNA transposon / MULE",      COLOR_DNA),
    "DNA/CACTA":    ("DNA transposon / CACTA",     COLOR_DNA),
    "Helitron":     ("Helitron / rolling-circle",  COLOR_HELITRON),
}
NON_CANON = {"Unannotated", "transposable_element_gene"}

def superfamily(alias):
    """Conservative whitelist mapping. ATHILA collapsed into Gypsy."""
    a = alias.upper()
    if "ATCOPIA"  in a: return "LTR/Copia"
    if "ATGP"     in a: return "LTR/Gypsy"
    if "ATHILA"   in a: return "LTR/Gypsy"      # Athila is a Gypsy clade
    if "ATLINE"   in a: return "non-LTR/LINE"
    if "ATSINE"   in a: return "non-LTR/SINE"
    if "VANDAL"   in a: return "DNA/MULE"
    if "ATENSPM"  in a: return "DNA/CACTA"
    if "HELITRON" in a: return "Helitron"
    if "TRANSPOSABLE_ELEMENT_GENE" in a: return "transposable_element_gene"
    if a in ("UNKNOWN", "UNASSIGNED", ""): return "Unannotated"
    return f"Other / unclassified — {alias}"

def is_athila(alias): return "ATHILA" in alias.upper()

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out-dir", required=True,
                   help="Pipeline output directory")
    p.add_argument("--te-gff",  required=True,
                   help="TE annotation GFF (chrom names must match BAM)")
    p.add_argument("--output",  default="figure_breakdown.png",
                   help="Output PNG (default: figure_breakdown.png)")
    p.add_argument("--bam",     default=None,
                   help="Override BAM path (default: <out-dir>/aligned/all.bam)")
    p.add_argument("--nuclear-regex", default=DEFAULT_NUCLEAR_REGEX,
                   help=f"Regex matching nuclear-chromosome contig SNs in the BAM "
                        f"header (default = TAIR10: {DEFAULT_NUCLEAR_REGEX})")
    p.add_argument("--mt-contig", default=DEFAULT_MT_CONTIG,
                   help=f"Mitochondrion contig SN; pass empty string '' to skip "
                        f"(default = TAIR10 {DEFAULT_MT_CONTIG})")
    p.add_argument("--pt-contig", default=DEFAULT_PT_CONTIG,
                   help=f"Plastid contig SN; pass empty string '' to skip "
                        f"(default = TAIR10 {DEFAULT_PT_CONTIG})")
    p.add_argument("--chrom-labels", default=None,
                   help="Optional JSON mapping of contig SN → friendly label "
                        "for panel (b). Default = TAIR10 NC_*→Chr1-5/Mt/Pt; "
                        "unmatched contigs use the raw SN.")
    p.add_argument("--title",   default="nanoALE-seq read distribution and TE family coverage",
                   help="Figure title")
    return p.parse_args()

def main():
    args = parse_args()
    out_dir = Path(args.out_dir).expanduser().resolve()
    bam = Path(args.bam) if args.bam else out_dir / "aligned" / "all.bam"
    cov_gff = out_dir / "target" / "te_overlap_counts.gff"      # produced by pipeline (bedtools -c)
    if not cov_gff.exists():
        # Fall back: compute on the fly
        cov_gff = out_dir / "target" / "_te_coverage_tmp.gff"
        cmd = (f"bedtools intersect -a {args.te_gff} -b {bam} -c "
               f"> {cov_gff}")
        print(f"Computing TE-coverage GFF on the fly: {cmd}")
        subprocess.check_call(cmd, shell=True)

    # ---- load coverage + family aggregation ----
    fam = defaultdict(lambda: dict(n=0, ov=0))
    super_fam = defaultdict(lambda: dict(n=0, ov=0))
    athila_only = dict(n=0, ov=0)
    atc78 = dict(n=0, ov=0)
    with open(cov_gff) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            p = line.rstrip().split("\t")
            if len(p) < 10: continue
            try: ov = int(p[-1])
            except ValueError: continue
            attr = p[8]
            m = (re.search(r"Alias=([^;]+)", attr) or re.search(r"Note=([^;]+)", attr)
                 or re.search(r"Name=([^;]+)", attr))
            a = m.group(1) if m else "UNKNOWN"
            sf = superfamily(a)
            fam[a]["n"]  += 1; fam[a]["ov"]  += ov
            super_fam[sf]["n"]  += 1; super_fam[sf]["ov"] += ov
            if is_athila(a):
                athila_only["n"]  += 1; athila_only["ov"] += ov
            if a == "ATCOPIA78":
                atc78["n"]  += 1; atc78["ov"] += ov

    # ---- figure ----
    fig = plt.figure(figsize=(13, 9.5))
    gs  = fig.add_gridspec(2, 2, height_ratios=[1.05, 1], hspace=0.5, wspace=0.3)

    # Panel (a)
    ax_a = fig.add_subplot(gs[0, :])
    top = [(sf, d) for sf, d in super_fam.items() if d["ov"] > 0 and sf not in NON_CANON]
    top.sort(key=lambda kv: -kv[1]["ov"])
    top = top[:10]
    rows_data = top + ([("ATCOPIA78 (ONSEN)", atc78)] if atc78["ov"] > 0 else [])
    rows_data = list(reversed(rows_data))
    labels = [r[0] for r in rows_data]
    counts = [r[1]["ov"] for r in rows_data]
    colors = []
    for sf, _ in rows_data:
        if sf == "ATCOPIA78 (ONSEN)":   colors.append(COLOR_ONSEN_HL)
        elif sf in CLASS_OF:            colors.append(CLASS_OF[sf][1])
        else:                           colors.append(COLOR_OTHER)
    y = np.arange(len(rows_data))
    ax_a.barh(y, counts, color=colors, edgecolor="white", linewidth=0.8)
    ax_a.set_yticks(y); ax_a.set_yticklabels(labels, fontsize=10)
    ax_a.set_xlabel("Read–TE overlaps", fontsize=11)
    ax_a.set_title("(a) TE superfamily distribution — colour-coded by class", fontsize=12, loc="left")
    ax_a.spines[["top","right"]].set_visible(False)
    for i, c in enumerate(counts):
        txt = f"  {c:,}"
        if labels[i] == "ATCOPIA78 (ONSEN)" and atc78["n"] > 0:
            txt += f"  (n={atc78['n']} loci)"
        elif labels[i] == "LTR/Gypsy" and athila_only["ov"] > 0:
            txt += f"  (incl. {athila_only['ov']:,} Athila clade)"
        ax_a.text(c, i, txt, va="center", fontsize=9)
    ax_a.legend(handles=[
        mpatches.Patch(facecolor=COLOR_LTR_COPIA, label="LTR retro / Copia"),
        mpatches.Patch(facecolor=COLOR_LTR_GYPSY, label="LTR retro / Gypsy (incl. Athila)"),
        mpatches.Patch(facecolor=COLOR_NON_LTR,   label="non-LTR retro (LINE / SINE)"),
        mpatches.Patch(facecolor=COLOR_DNA,       label="DNA transposon (MULE / CACTA)"),
        mpatches.Patch(facecolor=COLOR_HELITRON,  label="Helitron / rolling-circle"),
        mpatches.Patch(facecolor=COLOR_ONSEN_HL,  label="ATCOPIA78 / ONSEN positive control"),
    ], loc="lower right", frameon=False, fontsize=9, ncol=2)

    # Panel (b)
    ax_b = fig.add_subplot(gs[1, 0])
    chrom_labels = dict(DEFAULT_CHROM_LABELS)
    if args.chrom_labels:
        chrom_labels.update(json.loads(args.chrom_labels))
    nuclear_re = re.compile(args.nuclear_regex)
    out = subprocess.check_output(["samtools","idxstats", str(bam)]).decode().splitlines()
    chrom_table = []
    for line in out:
        cols = line.split("\t")
        if len(cols) < 3 or cols[0] == "*": continue
        size, nmap = int(cols[1]), int(cols[2])
        if size == 0 or nmap == 0: continue
        sn = cols[0]
        label = chrom_labels.get(sn, sn)
        if nuclear_re.match(sn):
            order = (0, label)
        elif sn == args.mt_contig or label == "Mt":
            order = (1, label)
        elif sn == args.pt_contig or label == "Pt":
            order = (2, label)
        else:
            order = (3, label)
        chrom_table.append((label, size, nmap, order))
    chrom_table.sort(key=lambda x: x[3])
    cls   = [t[0] for t in chrom_table]
    rpkb  = [t[2]/(t[1]/1000) for t in chrom_table]
    def _bar_color(l):
        if l == "Mt": return COLOR_NON_LTR
        if l == "Pt": return COLOR_HELITRON
        # Nuclear: any label that came from the nuclear regex bucket
        for entry in chrom_table:
            if entry[0] == l and entry[3][0] == 0:
                return "#3a4d7c"
        return "#888"
    bcols = [_bar_color(l) for l in cls]
    ax_b.bar(cls, rpkb, color=bcols, edgecolor="white", linewidth=0.5)
    ax_b.set_yscale("log")
    ax_b.set_ylabel("Mapped reads / kb (log)", fontsize=10)
    ax_b.set_title("(b) Per-chromosome read density", fontsize=12, loc="left")
    ax_b.spines[["top","right"]].set_visible(False)
    for i, d in enumerate(rpkb):
        ax_b.text(i, d*1.15, f"{d:.1f}", ha="center", va="bottom", fontsize=8)

    # Panel (c) — length distribution per source.
    # Sources are derived from BAM header using --nuclear-regex / --mt-contig / --pt-contig.
    ax_c = fig.add_subplot(gs[1, 1])
    bam_header = subprocess.check_output(["samtools","view","-H", str(bam)]).decode().splitlines()
    bam_contigs = [ln.split("\t")[1][3:] for ln in bam_header if ln.startswith("@SQ\t")]
    nuclear_contigs = [c for c in bam_contigs if nuclear_re.match(c)]
    sources = []
    if nuclear_contigs:
        sources.append(("Nuclear mapped", " ".join(nuclear_contigs), "-F", COLOR_LTR_COPIA))
    if args.mt_contig and args.mt_contig in bam_contigs:
        sources.append(("Mt mapped",      args.mt_contig, "-F", COLOR_NON_LTR))
    if args.pt_contig and args.pt_contig in bam_contigs:
        sources.append(("Pt mapped",      args.pt_contig, "-F", COLOR_HELITRON))
    sources.append(("Unmapped",           "",             "-f", "#888888"))
    arrays = []
    for label, refs, flag, color in sources:
        if flag == "-F":
            cmd = f"samtools view -F 0x904 {bam} {refs} | awk '{{print length($10)}}'"
        else:
            cmd = f"samtools view -f 4 {bam} | awk '{{print length($10)}}'"
        out = subprocess.check_output(cmd, shell=True).decode().split()
        arr = np.array([int(x) for x in out], dtype=int)
        if len(arr): arrays.append((arr, label, color))
    if arrays:
        max_len = max(a.max() for a,_,_ in arrays)
        bins = np.logspace(np.log10(50), np.log10(max_len), 60)
        for arr, label, color in arrays:
            m = arr.mean()
            ax_c.hist(arr, bins=bins, histtype="step", linewidth=1.6, color=color,
                      label=f"{label}  (n={len(arr):,}, μ={m:.0f} bp)", density=True)
            ax_c.axvline(m, color=color, linewidth=0.7, linestyle=":", alpha=0.55)
        ax_c.set_xscale("log")
        ax_c.set_xlabel("Read length (bp, log)", fontsize=10)
        ax_c.set_ylabel("Density", fontsize=10)
        ax_c.set_title("(c) Read length distribution by source", fontsize=12, loc="left")
        ax_c.legend(frameon=False, fontsize=8, loc="upper right")
        ax_c.spines[["top","right"]].set_visible(False)

    fig.suptitle(args.title, fontsize=13, x=0.04, ha="left", y=1.0)
    plt.savefig(args.output, dpi=300, bbox_inches="tight")
    print(f"Wrote {args.output}")

if __name__ == "__main__":
    main()
