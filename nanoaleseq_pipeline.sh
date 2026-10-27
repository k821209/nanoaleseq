#!/usr/bin/env bash
# nanoaleseq_pipeline.sh
#
# End-to-end nanoALE-seq analysis pipeline.
# Takes a Nanopore FASTQ and produces:
#   1) cleaned reads (cutadapt-trimmed, adapter-concatemer-filtered FASTQ)
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
#   --te-feature TYPE  GFF column-3 type counted as a TE feature (transposable_element)
#   --no-concatemer-filter      keep adapter arrays intact (default: split them out)
#   --min-insert N              shortest insert kept after splitting an array (100)
#   --min-qual N                mean-base-quality floor on kept reads (10; 0 disables)
#   --concatemer-mismatch N     mismatches allowed when finding an internal adapter (0)
#   --concatemer-probe SEQ      extra internal-adapter probe; repeatable
#                               (default probes: PCR-F and PCR-R-RC)
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
TE_FEATURE='transposable_element'          # GFF column-3 type counted as a "TE feature" at stage 5
LTR_FAM_REGEX='^(ATCOPIA|ATGP|ATHILA|ATLANTYS)'
NUCLEAR_REGEX='^NC_00(3070\.9|3071\.7|3074\.8|3075\.7|3076\.8)$'   # TAIR10 NCBI Chr1-5

CONCATEMER_FILTER=1                        # drop adapter self-ligation concatemers
CONCATEMER_MM=0                            # mismatches allowed when locating internal adapter
MIN_INSERT=100                             # shortest insert worth keeping after splitting an array
MIN_QUAL=10                                # mean-base-quality floor; adapter remnants sit near Q5
CONCATEMER_PROBES=()                       # extra probes; PCR_F / PCR_R_RC are always included

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
        --te-feature)    TE_FEATURE="$2"; shift 2 ;;
        --ltr-regex)     LTR_FAM_REGEX="$2"; shift 2 ;;
        --nuclear-regex) NUCLEAR_REGEX="$2"; shift 2 ;;
        --no-concatemer-filter)  CONCATEMER_FILTER=0; shift 1 ;;
        --concatemer-mismatch)   CONCATEMER_MM="$2"; shift 2 ;;
        --min-insert)            MIN_INSERT="$2"; shift 2 ;;
        --min-qual)              MIN_QUAL="$2"; shift 2 ;;
        --concatemer-probe)      CONCATEMER_PROBES+=("$2"); shift 2 ;;
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
log "Step 1/7  cutadapt — STRICT linked-adapter trim (require both PCR-F at 5' AND PCR-R-RC at 3')"
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
log "Step 2/7  Split ADAPTER-CONCATEMER arrays, keep the embedded insert"
# A large fraction of every library measured so far is not sequence at all: the
# ligation adapter (PCR-F + spacer + T7 promoter) self-ligates into head-to-tail
# tandem arrays which are then amplified by the PCR step. Measured 2026-10-26 on
# four libraries, as a share of reads SURVIVING the linked-adapter trim above:
#
#     Cuscuta japonica Muju SD (22 PCR cycles)   88.7%
#     Cuscuta japonica Muju ST (26 PCR cycles)   98.3%
#     Arabidopsis Col-0 heat   (benchmark)       94.2%
#     Arabidopsis ddm1         (benchmark)       97.1%
#
# They pass stage 1 because the array begins with a genuine PCR-F at the 5' end
# and ends with PCR-R-RC, which is exactly what the linked adapter asks for.
# They survive the 0.5x bead clean-up because a 10-unit array is ~460 bp, not 46.
# More PCR cycles make it worse: the short template out-competes real inserts,
# which is why the 26-cycle arm is the dirtier of the matched pair.
#
# Signature, and why "internal adapter" is a safe test: consecutive PCR-F hits in
# these reads sit at an exact 46 bp period (225,241 of ~240,000 observed gaps),
# the unit being ACACGACGCTCTTCCGATCT AGAG AGTAATACGACTCACTATAGGG. A genuine
# amplicon carries PCR-F ONLY at the 5' end, which stage 1 has already removed,
# so any remaining internal occurrence means a second adapter copy in the molecule.
#
# The splitter also cuts at the T7 promoter and applies a mean-quality floor.
# Both were added after the first version let adapter remnants through: an array
# whose PCR-F copies had drifted survived as a >=100 bp "insert" made entirely of
# adapter, sitting near Q5 against the library's Q25, and such reads aligned 59 bp
# at a time to whatever matched. In a parasite experiment they produced an
# apparent 6-fold enrichment of host-derived sequence in the attached arm that
# was nothing but this artefact.
if [[ "$CONCATEMER_FILTER" -eq 1 ]]; then
    PROBE_ARGS=(--probe "$PCR_F" --probe "$PCR_R_RC")
    for pr in ${CONCATEMER_PROBES+"${CONCATEMER_PROBES[@]}"}; do PROBE_ARGS+=(--probe "$pr"); done

    # Diagnostic: the modal spacing between consecutive adapter hits IS the adapter
    # unit length when the reads are concatemers. Logged so a changed adapter design
    # shows up instead of silently passing the filter.
    seqkit locate -j "$THREADS" -m "$CONCATEMER_MM" -p "$PCR_F" -p "$PCR_R_RC" \
        "$OUT/clean/trimmed.fastq.gz" 2> "$OUT/logs/concatemer_locate.log" \
      | tail -n +2 | sort -k1,1 -k5,5n \
      | awk -F'\t' '$1==prev{d=$5-ps; if(d>0 && d<1000) g[d]++} {prev=$1; ps=$5}
            END{best=0; bn=0; tot=0; for(k in g){tot+=g[k]; if(g[k]>bn){bn=g[k]; best=k}}
                if(tot>0) printf "  modal adapter period: %d bp (%.1f%% of %d gaps)\n", best, 100*bn/tot, tot;
                else      print "  modal adapter period: n/a (no read carried 2+ hits)"}'

    "$(dirname "$0")/strip_adapter_concatemer.py" \
        --in "$OUT/clean/trimmed.fastq.gz" \
        --out "$OUT/clean/noconcat.fastq.gz" \
        --min-insert "$MIN_INSERT" \
        --min-qual "$MIN_QUAL" \
        --stats "$OUT/clean/concatemer_stats.tsv" \
        "${PROBE_ARGS[@]}" 2>&1 | sed 's/^/  /'
    READS_FOR_ALIGNMENT="$OUT/clean/noconcat.fastq.gz"
else
    log "  SKIPPED (--no-concatemer-filter). Expect the alignment stages to be"
    log "  dominated by insert-free adapter arrays; the funnel below will not be"
    log "  comparable to runs that applied the filter."
    cp "$OUT/clean/trimmed.fastq.gz" "$OUT/clean/noconcat.fastq.gz"
    READS_FOR_ALIGNMENT="$OUT/clean/noconcat.fastq.gz"
fi
log "  reads carried to alignment: $(seqkit stats -T "$READS_FOR_ALIGNMENT" | tail -1 | awk '{print $4}')"

# ----------------------------------------------------------------------
log "Step 3/7  minimap2 — map-ont alignment to reference"
minimap2 -ax map-ont -t "$THREADS" "$REF" "$READS_FOR_ALIGNMENT" 2> "$OUT/logs/minimap2.log" \
    | samtools sort -@ "$THREADS" -o "$OUT/aligned/all.bam" -
samtools index "$OUT/aligned/all.bam"

# ----------------------------------------------------------------------
log "Step 4/7  Filter to nuclear chromosomes (regex: $NUCLEAR_REGEX)"
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
log "Step 5/7  Filter to TE-overlapping reads (column-3 type == '$TE_FEATURE')"
# Intersect against ONLY the TE features of the GFF, not every record. A combined
# genes+transposons GFF (e.g. TAIR10_GFF3_genes_transposons.gff) also contains
# genome-spanning 'chromosome', 'gene' and 'exon' records; intersecting against
# those would make essentially every nuclear read "overlap a feature" and inflate
# this stage to ~100%. Restricting to the TE feature type keeps stage 5 meaning
# "overlaps an annotated TE", matching the manuscript funnel.
awk -F'\t' -v ft="$TE_FEATURE" '$3==ft' "$TE_GFF" > "$OUT/target/te_features.gff"
log "  TE features ($TE_FEATURE): $(wc -l < "$OUT/target/te_features.gff")"
bedtools intersect -a "$OUT/target/nuclear.bam" -b "$OUT/target/te_features.gff" -u \
    > "$OUT/target/te_overlap.bam" 2> /dev/null
samtools index "$OUT/target/te_overlap.bam"

# ----------------------------------------------------------------------
log "Step 6/7  Filter to LTR retrotransposon family overlap (regex: $LTR_FAM_REGEX)"
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
log "Step 7/7  Funnel report"
count_unique_reads_in_bam () {  # primary records only
    samtools view -c -F 0x904 "$1"
}
RAW=$(seqkit stats -T "$FASTQ" | tail -1 | awk '{print $4}')
TRIMMED=$(seqkit stats -T "$OUT/clean/trimmed.fastq.gz" | tail -1 | awk '{print $4}')
NOCONCAT=$(seqkit stats -T "$READS_FOR_ALIGNMENT" | tail -1 | awk '{print $4}')
if [[ "$CONCATEMER_FILTER" -eq 1 ]]; then
    CONCAT_DESC="After splitting adapter-concatemer arrays (insert >= ${MIN_INSERT} bp, mean Q >= ${MIN_QUAL})"
else
    CONCAT_DESC="Adapter-concatemer filter DISABLED (--no-concatemer-filter) - count unchanged"
fi
S3=$(count_unique_reads_in_bam "$OUT/aligned/all.bam")
S4=$(count_unique_reads_in_bam "$OUT/target/nuclear.bam")
S5=$(count_unique_reads_in_bam "$OUT/target/te_overlap.bam")
S6=$(count_unique_reads_in_bam "$OUT/target/ltr_retro.bam")
S7=$(count_unique_reads_in_bam "$OUT/target/fl_ltr_retro.bam")

python3 - <<PYEND > "$OUT/target/funnel_report.tsv"
stages = [
  ('1', 'Raw Nanopore reads',                                                                                  $RAW),
  ('2', 'After cutadapt PCR-F/R linked-adapter trim (>= ${MIN_LEN} bp, --discard-untrimmed)',                  $TRIMMED),
  ('3', '${CONCAT_DESC}',                                                                                       $NOCONCAT),
  ('4', 'Primary aligned to reference (any contig, -F 0x904)',                                                  $S3),
  ('5', 'Primary aligned to nuclear chromosomes only (regex match: ${NUCLEAR_REGEX})',                          $S4),
  ('6', 'Overlapping any TE feature in --te-gff',                                                               $S5),
  ('7', 'Overlapping an LTR retrotransposon family (regex: ${LTR_FAM_REGEX})',                                  $S6),
  ('8', 'Overlapping a full-length LTR retrotransposon (annotated length >= ${FL_THRESHOLD} bp) — FINAL TARGET',$S7),
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
log "Concatemer-free reads:              $OUT/clean/noconcat.fastq.gz"
log "Concatemer split stats:             $OUT/clean/concatemer_stats.tsv"
log "All stage BAMs:                     $OUT/target/*.bam"
log "Funnel report:                       $OUT/target/funnel_report.tsv"
