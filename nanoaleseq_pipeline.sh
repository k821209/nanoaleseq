#!/usr/bin/env bash
# nanoaleseq_pipeline.sh
#
# End-to-end nanoALE-seq analysis pipeline.
# Takes a Nanopore FASTQ and produces:
#   1) cleaned reads (cutadapt-trimmed FASTQ)
#   2) genome-aligned BAM
#   3) successively filtered subset BAMs at every funnel stage
#   4) the final target FASTQ (reads overlapping full-length LTR retrotransposons)
#   5) a TSV funnel report (read counts per stage)
#
# Usage:
#   nanoaleseq_pipeline.sh \
#       --fastq    raw.fastq.gz \
#       --ref      TAIR10.fa \
#       --te-gff   TAIR_TEs.gff   (with TE annotation; Alias= attribute used for family match) \
#       --out      out_dir/ \
#       --threads  16
#
# Designed for: Arabidopsis nanoALE-seq with PCR-F=ACACGACGCTCTTCCGATCT, PCR-R=ACGCTCGACTAACTTGTACC.
# Override the primer sequences below if running on a different species/design.
#
# Optional flags (with defaults):
#   --pcr-f SEQ        5' technical primer                     (ACACGACGCTCTTCCGATCT)
#   --pcr-r SEQ        3' technical primer; auto-RC computed   (ACGCTCGACTAACTTGTACC)
#   --min-len N        cutadapt --minimum-length               (50)
#   --fl-threshold N   bp threshold for "full-length" LTR retro (4500)
#   --ltr-regex RX     family-name regex for LTR retro         (^(ATCOPIA|ATGP|ATHILA|ATLANTYS))
#   --nuclear-regex RX BAM-header SN regex for nuclear contigs (TAIR10 NCBI Chr1-5)
#   --threads N                                                (8)
#
# Dependencies (in $PATH): cutadapt >=4, minimap2, samtools, bedtools, awk, python3.

set -euo pipefail

# ----------------------------------------------------------------------
# Defaults / CLI parsing
PCR_F='ACACGACGCTCTTCCGATCT'
PCR_R='ACGCTCGACTAACTTGTACC'
PCR_R_RC='GGTACAAGTTAGTCGAGCGT'           # reverse-complement of PCR_R, hits the 3' end of reads
MIN_LEN=50                                 # cutadapt --minimum-length
FL_THRESHOLD=4500                          # bp; used to flag full-length LTR retros
LTR_FAM_REGEX='^(ATCOPIA|ATGP|ATHILA|ATLANTYS)'
NUCLEAR_REGEX='^NC_00(3070\.9|3071\.7|3074\.8|3075\.7|3076\.8)$'   # TAIR10 NCBI Chr1-5

THREADS=8
FASTQ=
REF=
TE_GFF=
OUT=

while [[ $# -gt 0 ]]; do
    case "$1" in
        --fastq)         FASTQ="$2"; shift 2 ;;
        --ref)           REF="$2"; shift 2 ;;
        --te-gff)        TE_GFF="$2"; shift 2 ;;
        --out)           OUT="$2"; shift 2 ;;
        --threads)       THREADS="$2"; shift 2 ;;
        --pcr-f)         PCR_F="$2"; shift 2 ;;
        --pcr-r)         PCR_R="$2"; PCR_R_RC=$(python3 -c "import sys; s='$2'; tr=str.maketrans('ACGTN','TGCAN'); print(s.translate(tr)[::-1])"); shift 2 ;;
        --min-len)       MIN_LEN="$2"; shift 2 ;;
        --fl-threshold)  FL_THRESHOLD="$2"; shift 2 ;;
        --ltr-regex)     LTR_FAM_REGEX="$2"; shift 2 ;;
        --nuclear-regex) NUCLEAR_REGEX="$2"; shift 2 ;;
        -h|--help) sed -n '1,/^# Dependencies/p' "$0"; exit 0 ;;
        *) echo "Unknown arg: $1" >&2; exit 2 ;;
    esac
done

for v in FASTQ REF TE_GFF OUT; do
    if [[ -z "${!v}" ]]; then
        flag=$(echo "$v" | tr '[:upper:]_' '[:lower:]-')
        echo "Missing --$flag" >&2; exit 2
    fi
done

mkdir -p "$OUT"/{clean,aligned,target,logs}
LOG="$OUT/logs/pipeline.log"
exec > >(tee -a "$LOG") 2>&1
log() { echo "[$(date '+%F %T')] $*"; }

log "=== nanoALE-seq pipeline start ==="
log "fastq=$FASTQ  ref=$REF  te_gff=$TE_GFF  out=$OUT  threads=$THREADS"

# ----------------------------------------------------------------------
log "Step 1/6  cutadapt — STRICT linked-adapter trim (require both PCR-F at 5' AND PCR-R-RC at 3')"
# Linked-adapter syntax: -g 'ADAPTER1...ADAPTER2' enforces that BOTH adapters be
# found in the same read; reads with only one adapter (e.g. residual ds-DNA
# carrying PCR-F at both ends due to incomplete TURBO DNase digestion) are
# discarded. This was a bug fix from the original loose `-g … -a …` form.
cutadapt -j "$THREADS" \
         -g "${PCR_F}...${PCR_R_RC}" \
         --minimum-length "$MIN_LEN" --discard-untrimmed \
         -o "$OUT/clean/trimmed.fastq.gz" \
         "$FASTQ" \
         > "$OUT/logs/cutadapt.log" 2>&1
log "  trimmed reads: $(seqkit stats -T "$OUT/clean/trimmed.fastq.gz" | tail -1 | awk '{print $4}')"

# ----------------------------------------------------------------------
log "Step 2/6  minimap2 — map-ont alignment to reference"
minimap2 -ax map-ont -t "$THREADS" "$REF" "$OUT/clean/trimmed.fastq.gz" 2> "$OUT/logs/minimap2.log" \
    | samtools sort -@ "$THREADS" -o "$OUT/aligned/all.bam" -
samtools index "$OUT/aligned/all.bam"

# ----------------------------------------------------------------------
log "Step 3/6  Filter to nuclear chromosomes (regex: $NUCLEAR_REGEX)"
NUCLEAR_REFS=$(samtools view -H "$OUT/aligned/all.bam" \
    | awk -v rx="$NUCLEAR_REGEX" '$1=="@SQ"{for(i=2;i<=NF;i++) if(match($i,/^SN:/)) {sn=substr($i,4); if(sn ~ rx) print sn}}')
if [[ -z "$NUCLEAR_REFS" ]]; then
    log "  ERROR: no contigs matched NUCLEAR_REGEX. Header SQ:"
    samtools view -H "$OUT/aligned/all.bam" | grep '^@SQ' | head
    exit 3
fi
log "  nuclear contigs: $NUCLEAR_REFS"
# shellcheck disable=SC2086
samtools view -F 0x904 -b "$OUT/aligned/all.bam" $NUCLEAR_REFS \
    > "$OUT/target/nuclear.bam"
samtools index "$OUT/target/nuclear.bam"

# ----------------------------------------------------------------------
log "Step 4/6  Filter to TE-overlapping reads (any TE feature)"
bedtools intersect -a "$OUT/target/nuclear.bam" -b "$TE_GFF" -u \
    > "$OUT/target/te_overlap.bam" 2> /dev/null
samtools index "$OUT/target/te_overlap.bam"

# ----------------------------------------------------------------------
log "Step 5/6  Filter to LTR retrotransposon family overlap (regex: $LTR_FAM_REGEX)"
awk -F'\t' -v rx="$LTR_FAM_REGEX" '
  match($9, /Alias=[^;]+/) {
    a = substr($9, RSTART+6, RLENGTH-6);
    if (a ~ rx) print
  }' "$TE_GFF" > "$OUT/target/ltr_retro.gff"
log "  LTR retro features: $(wc -l < "$OUT/target/ltr_retro.gff")"

bedtools intersect -a "$OUT/target/nuclear.bam" -b "$OUT/target/ltr_retro.gff" -u \
    > "$OUT/target/ltr_retro.bam" 2> /dev/null
samtools index "$OUT/target/ltr_retro.bam"

# Sub-filter: full-length LTR retros only
awk -F'\t' -v fl="$FL_THRESHOLD" '($5-$4) >= fl' "$OUT/target/ltr_retro.gff" \
    > "$OUT/target/fl_ltr_retro.gff"
log "  full-length LTR retro features (>= $FL_THRESHOLD bp): $(wc -l < "$OUT/target/fl_ltr_retro.gff")"

bedtools intersect -a "$OUT/target/nuclear.bam" -b "$OUT/target/fl_ltr_retro.gff" -u \
    > "$OUT/target/fl_ltr_retro.bam" 2> /dev/null
samtools index "$OUT/target/fl_ltr_retro.bam"

# Final target FASTQ — reads overlapping any FL LTR retrotransposon
samtools fastq "$OUT/target/fl_ltr_retro.bam" 2> /dev/null \
    | gzip > "$OUT/target/final_target.fastq.gz"

# ----------------------------------------------------------------------
log "Step 6/6  Funnel report"
count_unique_reads_in_bam () {  # primary records only
    samtools view -c -F 0x904 "$1"
}
RAW=$(seqkit stats -T "$FASTQ" | tail -1 | awk '{print $4}')
TRIMMED=$(seqkit stats -T "$OUT/clean/trimmed.fastq.gz" | tail -1 | awk '{print $4}')
S3=$(count_unique_reads_in_bam "$OUT/aligned/all.bam")
S4=$(count_unique_reads_in_bam "$OUT/target/nuclear.bam")
S5=$(count_unique_reads_in_bam "$OUT/target/te_overlap.bam")
S6=$(count_unique_reads_in_bam "$OUT/target/ltr_retro.bam")
S7=$(count_unique_reads_in_bam "$OUT/target/fl_ltr_retro.bam")

python3 - <<PYEND > "$OUT/target/funnel_report.tsv"
stages = [
  ('1', 'Raw Nanopore reads',                                                                                  $RAW),
  ('2', 'After cutadapt PCR-F/R linked-adapter trim (>= ${MIN_LEN} bp, --discard-untrimmed)',                  $TRIMMED),
  ('3', 'Primary aligned to reference (any contig, -F 0x904)',                                                  $S3),
  ('4', 'Primary aligned to nuclear chromosomes only (regex match: ${NUCLEAR_REGEX})',                          $S4),
  ('5', 'Overlapping any TE feature in --te-gff',                                                               $S5),
  ('6', 'Overlapping an LTR retrotransposon family (regex: ${LTR_FAM_REGEX})',                                  $S6),
  ('7', 'Overlapping a full-length LTR retrotransposon (annotated length >= ${FL_THRESHOLD} bp) — FINAL TARGET',$S7),
]
print('stage\tdescription\tread_count\tpct_of_raw\tpct_of_previous')
prev = stages[0][2]; raw = stages[0][2]
for tag, desc, n in stages:
    pct_raw  = 100*n/raw  if raw  else 0
    pct_prev = 100*n/prev if prev else 0
    print(f'{tag}\t{desc}\t{n}\t{pct_raw:.4f}%\t{pct_prev:.2f}%')
    prev = n
PYEND

cat "$OUT/target/funnel_report.tsv"

log "=== DONE ==="
log "Final target reads (FL LTR retro):  $OUT/target/final_target.fastq.gz"
log "All stage BAMs:                     $OUT/target/*.bam"
log "Funnel report:                       $OUT/target/funnel_report.tsv"
