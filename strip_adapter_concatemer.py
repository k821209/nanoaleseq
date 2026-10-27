#!/usr/bin/env python3
"""Split adapter self-ligation concatemers out of nanoALE-seq reads.

The ligation adapter (PCR-F + spacer + T7 promoter) self-ligates into head-to-tail
tandem arrays that are then amplified by PCR. They pass the linked-adapter trim
because the array starts with a genuine PCR-F and ends with PCR-R-RC.

Discarding such a read whole is wrong: 42.6% of them still carry a genuine genomic
fragment embedded in the array (measured on Arabidopsis Col-0, 2026-10-26), and
those fragments were the bulk of the old pipeline's final target. So instead of
dropping the read, cut it at every adapter occurrence and keep the longest
remaining segment, if it reaches --min-insert.

Reads with no internal adapter pass through untouched.

Two things the first version of this script got wrong, both found by following a
false positive to its source (2026-10-27):

  * It cut only at the PCR primers, not at the T7 promoter that the adapter also
    carries. An array whose PCR-F copies had drifted therefore survived as a
    "recovered insert" made entirely of adapter.
  * It applied no quality floor. Those survivors had a median read quality of
    about Q5, against Q25 for the library as a whole, and they aligned 59 bp at
    a time to whatever they happened to match. In one comparison they produced
    an apparent 6-fold enrichment of host-derived sequence in the attached arm
    of a parasite experiment, which was entirely this artefact.

So the T7 promoter is a default probe, and --min-qual drops reads whose mean
base quality is below the floor. Pass --min-qual 0 to disable.
"""
import argparse, gzip, sys

def opener(p, mode="rt"):
    return gzip.open(p, mode) if p.endswith(".gz") else open(p, mode)

def rc(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]

def cut_positions(seq, probes):
    """All [start,end) intervals matching any probe, on either strand."""
    hits = []
    for p in probes:
        for pat in (p, rc(p)):
            i = seq.find(pat)
            while i != -1:
                hits.append((i, i + len(pat)))
                i = seq.find(pat, i + 1)
    return sorted(hits)

def longest_segment(seq, qual, hits):
    """Longest stretch of sequence not covered by an adapter hit."""
    best = (0, 0)
    cur = 0
    for s, e in hits:
        if s - cur > best[1] - best[0]:
            best = (cur, s)
        cur = max(cur, e)
    if len(seq) - cur > best[1] - best[0]:
        best = (cur, len(seq))
    return seq[best[0]:best[1]], qual[best[0]:best[1]]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--probe", action="append", required=True,
                    help="adapter sequence to cut at; repeatable")
    ap.add_argument("--min-insert", type=int, default=100,
                    help="drop recovered segments shorter than this (default 100)")
    ap.add_argument("--min-qual", type=float, default=10.0,
                    help="drop reads whose mean base quality is below this (default 10; 0 disables)")
    ap.add_argument("--no-default-probes", action="store_true",
                    help="use only the probes given with --probe, omitting the T7 promoter")
    ap.add_argument("--stats", default=None)
    a = ap.parse_args()

    T7 = "TAATACGACTCACTATAGGG"          # the adapter carries this beside the PCR primer
    probes = [p.upper() for p in a.probe]
    if not a.no_default_probes and T7 not in probes:
        probes.append(T7)
    n = clean = split = recovered = dropped = lowq = 0
    rec_bp = []

    def meanq(q):
        return sum(ord(c) - 33 for c in q) / len(q) if q else 0.0

    with opener(a.inp) as fh, opener(a.out, "wt") as out:
        while True:
            h = fh.readline()
            if not h:
                break
            seq = fh.readline().rstrip("\n")
            fh.readline()
            qual = fh.readline().rstrip("\n")
            n += 1
            if a.min_qual > 0 and meanq(qual) < a.min_qual:
                lowq += 1
                continue
            hits = cut_positions(seq.upper(), probes)
            if not hits:
                clean += 1
                out.write(f"{h}{seq}\n+\n{qual}\n")
                continue
            split += 1
            s2, q2 = longest_segment(seq, qual, hits)
            if len(s2) >= a.min_insert and (a.min_qual <= 0 or meanq(q2) >= a.min_qual):
                recovered += 1
                rec_bp.append(len(s2))
                name = h[1:].split()[0]
                out.write(f"@{name} concatemer_split=1 segments={len(hits)+1}\n{s2}\n+\n{q2}\n")
            else:
                dropped += 1

    rec_bp.sort()
    med = rec_bp[len(rec_bp)//2] if rec_bp else 0
    lines = [
        f"reads_in\t{n}",
        f"dropped_low_quality\t{lowq}",
        f"min_qual\t{a.min_qual}",
        f"no_internal_adapter\t{clean}",
        f"carried_internal_adapter\t{split}",
        f"insert_recovered\t{recovered}",
        f"dropped_below_min_insert\t{dropped}",
        f"min_insert\t{a.min_insert}",
        f"recovered_median_bp\t{med}",
        f"reads_out\t{clean + recovered}",
    ]
    sys.stderr.write("\n".join(lines) + "\n")
    if a.stats:
        open(a.stats, "w").write("\n".join(lines) + "\n")

if __name__ == "__main__":
    main()
