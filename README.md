# nanoALE-seq

**Long-read profiling of active LTR retrotransposons via linear extrachromosomal DNA.**

`nanoALE-seq` is the Oxford Nanopore long-read adaptation of the **ALE-seq** method (Cho et al. 2019, *Nature Plants*) — it amplifies the LTR–primer-binding-site (PBS) region of linear extrachromosomal DNA (ecDNA) through T7 in vitro transcription and PBS-anchored reverse transcription, and sequences the products on Oxford Nanopore. The pipeline in this repository takes a raw Nanopore FASTQ as input.

**Core pipeline outputs** (generic, species-independent):

- a strict-cutadapt-trimmed FASTQ of canonical RT-derived reads,
- a minimap2-aligned BAM,
- a five-stage **funnel report** narrowing reads from raw → on-target retrotransposon (TSV).

**Benchmark / troubleshooting utilities** (provided alongside the pipeline; defaults are *Arabidopsis* TAIR10 — see CLI flags for retargeting):

- `plots/render_breakdown.py` — TE family / per-chromosome / read-length breakdown figure. Defaults assume TAIR10 (NC_* accessions for Chr1–5 + Mt + Pt). Override via `--nuclear-regex`, `--mt-contig`, `--pt-contig`, `--chrom-labels` for other species. Family-superfamily mapping in panel (a) is a TAIR10 conservative whitelist; non-matching families fall through to "Other / unclassified".
- `plots/render_onsen_pileup.py` — IGV-style positive-control pileup at the eight canonical full-length ATCOPIA78 (ONSEN) loci. Designed for the *Arabidopsis* ONSEN benchmark; pass `--loci` and `--feature-type` to retarget to a different positive-control set.

The companion *Methods in Molecular Biology* chapter (Ha, Kang & Cho, 2026) describes the wet-lab protocol and benchmarks the pipeline on heat-treated *Arabidopsis thaliana* Col-0 leaf, where it recovers all eight full-length ATCOPIA78 (ONSEN) reference loci as a positive control.

## Repository layout

```
nanoaleseq/
├── README.md                          # this file
├── LICENSE                            # MIT
├── environment.yml                    # conda environment
├── nanoaleseq_pipeline.sh             # main end-to-end pipeline
└── plots/
    ├── render_onsen_pileup.py         # IGV-style coverage figure
    └── render_breakdown.py            # TE family / per-chrom / length figure
```

## Install

```bash
git clone https://github.com/k821209/nanoaleseq.git
cd nanoaleseq
conda env create -f environment.yml
conda activate nanoaleseq
chmod +x nanoaleseq_pipeline.sh
```

That installs the four binaries the pipeline calls (`cutadapt`, `minimap2`, `samtools`, `bedtools`, `seqkit`) plus the Python plotting dependencies (`pysam`, `matplotlib`, `numpy`).

## Quick start (Arabidopsis ONSEN benchmark)

```bash
# Reference + TE annotation
wget -O TAIR10.fa.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz
gunzip TAIR10.fa.gz
wget -O TAIR10_TEs.gff 'https://www.arabidopsis.org/api/download-files/download?filePath=Genes%2FTAIR10_genome_release%2FTAIR10_gff3%2FTAIR10_GFF3_genes_transposons.gff'

# End-to-end pipeline (core)
./nanoaleseq_pipeline.sh \
    --fastq    raw.fastq.gz \
    --ref      TAIR10.fa \
    --te-gff   TAIR10_TEs.gff \
    --out      out_dir/ \
    --threads  16

# Optional: Arabidopsis ONSEN benchmark plots
python plots/render_breakdown.py \
    --out-dir out_dir/ --te-gff TAIR10_TEs.gff \
    --output  figure3_breakdown.png

# Optional: ONSEN positive-control pileup (Arabidopsis-specific)
python plots/render_onsen_pileup.py \
    --out-dir out_dir/ --te-gff TAIR10_TEs.gff \
    --output  figure2_pileup.png
```

## Inputs

| Flag | Description | Default |
|------|-------------|---------|
| `--fastq` | Raw Nanopore FASTQ (gzipped OK) | required |
| `--ref` | Reference genome FASTA | required |
| `--te-gff` | TE annotation GFF (chrom names matching `--ref`) | required |
| `--out` | Output directory (will be created) | required |
| `--threads` | Parallel threads | 8 |
| `--pcr-f` | 5' technical primer | `ACACGACGCTCTTCCGATCT` |
| `--pcr-r` | 3' technical primer (reverse-complemented automatically) | `ACGCTCGACTAACTTGTACC` |
| `--fl-threshold` | LTR retro length (bp) considered "full-length" | 4500 |
| `--ltr-regex` | Regex matching LTR-retro family aliases in the GFF | `^(ATCOPIA\|ATGP\|ATHILA\|ATLANTYS)` |
| `--nuclear-regex` | Regex matching nuclear-chromosome contig names | `^NC_00(3070\.9\|3071\.7\|3074\.8\|3075\.7\|3076\.8)$` |

The defaults target the *Arabidopsis* ALE-seq design from Cho et al. 2019. Override `--pcr-f` / `--pcr-r` for a different primer design, `--ltr-regex` and `--nuclear-regex` for a different organism / annotation.

## Outputs

```
out_dir/
├── clean/
│   ├── trimmed.fastq.gz             # strict-linked-adapter-trimmed reads
│   └── cutadapt.log
├── aligned/
│   ├── all.bam                      # full minimap2 alignment
│   └── all.bam.bai
├── target/
│   ├── nuclear.bam                  # primary, nuclear chromosomes only
│   ├── te_overlap.bam               # primary reads overlapping any TE
│   ├── ltr_retro.bam                # primary reads overlapping LTR retro family
│   ├── ltr_retro.gff                # filtered LTR retro features
│   ├── fl_ltr_retro.bam             # primary reads on full-length LTR retros
│   ├── fl_ltr_retro.gff             # FL LTR retro features
│   ├── final_target.fastq.gz        # **on-target FASTQ**
│   └── funnel_report.tsv            # **funnel table**
└── logs/
    └── pipeline.log
```

The **funnel report** is the headline product: every row is a strict subset of the row above, counted as unique primary reads (`bedtools intersect -u`). On the heat-treated Col-0 benchmark we obtain:

| Stage | Reads | % of raw | % of previous |
|---|---|---|---|
| 1. Raw Nanopore reads | 9,707,795 | 100.00 | — |
| 2. Strict-trimmed | 235,175 | 2.42 | 2.42 |
| 3. Primary aligned to TAIR10 | 81,614 | 0.84 | 34.70 |
| 4. Nuclear (Chr1–5) | 73,405 | 0.76 | 89.94 |
| 5. Overlapping any TE | 23,672 | 0.24 | 32.25 |
| 6. Overlapping LTR retro family | 16,616 | 0.17 | 70.19 |
| 7. Overlapping full-length LTR retro | 11,380 | 0.12 | 68.49 |
| 8. Overlapping FL ATCOPIA78 / ONSEN | 135 | 0.0014 | 1.19 |

## Why the strict cutadapt mode matters

`nanoaleseq_pipeline.sh` runs cutadapt in **linked-adapter mode** (`-g 'PCR-F…PCR-R-RC' --discard-untrimmed`), which requires every retained read to carry **both** the 5' PCR-F and the 3' PCR-R (reverse-complemented) adapter. This discards two artefact classes:

1. **Adapter-only concatemers** (~92% of the loose-trim pool) — small fractions of T7 adapter ligate head-to-tail and get amplified into long, repetitive, adapter-only reads.
2. **Residual ds-DNA reads** (~4% of the loose-trim pool) — incompletely DNase-digested T7 template carries the palindromic adapter at both ends; loose `-g … -a …` cutadapt lets these through.

Switching from loose to strict drops the trimmed pool from 8.06 × 10⁶ to 2.35 × 10⁵ reads (97% drop) **but raises the primary mapping rate from 7.67% to 34.70%** — i.e. the discarded reads were almost entirely artefacts. See the chapter's Note 9 for the full diagnostic.

## Citation

If you use this pipeline, please cite both:

1. **Ha, S., Kang, Y. J. & Cho, J.** (2026). *nanoALE-seq: long-read profiling of active LTR retrotransposons via linear extrachromosomal DNA.* In: *Methods in Molecular Biology*. Springer. (this chapter)
2. **Cho, J., Benoit, M., Catoni, M. et al.** (2019). *Sensitive detection of pre-integration intermediates of long terminal repeat retrotransposons in crop plants.* Nature Plants 5: 26–33. doi:[10.1038/s41477-018-0320-9](https://doi.org/10.1038/s41477-018-0320-9)

The original ALE-seq Illumina protocol chapter is also a useful reference:

3. **Wang, L., Kim, E. Y. & Cho, J.** (2021). *High-Throughput Profiling of Extrachromosomal Linear DNAs of Long Terminal Repeat Retrotransposons by ALE-seq.* In: *Plant Transposable Elements*, Methods in Molecular Biology. doi:[10.1007/978-1-0716-1134-0_9](https://doi.org/10.1007/978-1-0716-1134-0_9)

## License

MIT — see [LICENSE](LICENSE).

## Contact

- **Yang Jae Kang** (co-corresponding) — kangyangjae@gnu.ac.kr
- **Jungnam Cho** (co-corresponding) — jungnam.cho@durham.ac.uk

Issues and pull requests welcome via GitHub.
