#!/usr/bin/env python3
"""Check nanoALE-seq inserts for the primer-binding site, and optionally keep only those that carry it.

The method reverse-transcribes from a primer annealed at the PBS, so a genuine
product ends at the PBS. Measured on four libraries (2026-10-28), the share of
inserts whose 3' terminal 12 bp is TGGTATCAGAGC within two mismatches:

    Cuscuta japonica Muju SD   89.3%
    Cuscuta japonica Muju ST   77.3%
    Arabidopsis Col-0 heat     84.3%
    Arabidopsis ddm1           82.8%

TGGTATCAGAGC is the plant LTR-retrotransposon PBS, complementary to the 3' end
of initiator tRNA-Met. In those libraries it is 90-93% conserved per base while
the twelve bases before it run 31-63%, so the boundary is the motif's own.

What this does and does not establish: the PBS at a defined end shows the
reverse transcription was primed where the method says it primes, so the
templates carried a PBS. It does not by itself show the template was
extrachromosomal rather than a genomic fragment; that needs the positional
evidence (reads falling on the LTR repeat rather than the element interior).
"""
import argparse, gzip, sys

PBS_DEFAULT = "TGGTATCAGAGC"

def opener(p, m="rt"):
    return gzip.open(p, m) if p.endswith(".gz") else open(p, m)

def rc(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]

def mismatches(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)

def carries(seq, pbs, mm, where):
    L = len(pbs)
    if len(seq) < L:
        return False
    if where in ("3", "both") and mismatches(seq[-L:].upper(), pbs) <= mm:
        return True
    if where in ("5", "both") and mismatches(seq[:L].upper(), rc(pbs)) <= mm:
        return True
    return False

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", help="write kept reads here; omit to only report")
    ap.add_argument("--pbs", default=PBS_DEFAULT,
                    help=f"primer-binding site sequence (default {PBS_DEFAULT})")
    ap.add_argument("--pbs-mismatch", type=int, default=2)
    ap.add_argument("--pbs-end", choices=["3", "5", "both"], default="3",
                    help="which insert end must carry it (default 3)")
    ap.add_argument("--require-pbs", action="store_true",
                    help="drop inserts without the PBS; otherwise all are kept and only counted")
    ap.add_argument("--stats")
    a = ap.parse_args()

    pbs = a.pbs.upper()
    n = hit = 0
    out = opener(a.out, "wt") if a.out else None
    with opener(a.inp) as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            seq = fh.readline().rstrip("\n")
            fh.readline()
            q = fh.readline().rstrip("\n")
            n += 1
            ok = carries(seq, pbs, a.pbs_mismatch, a.pbs_end)
            if ok:
                hit += 1
            if out and (ok or not a.require_pbs):
                out.write(f"{h}{seq}\n+\n{q}\n")
    if out:
        out.close()

    kept = hit if a.require_pbs else n
    lines = [f"pbs_sequence\t{pbs}", f"pbs_end\t{a.pbs_end}",
             f"pbs_mismatch_allowed\t{a.pbs_mismatch}",
             f"inserts_in\t{n}", f"inserts_with_pbs\t{hit}",
             f"pbs_rate\t{hit/n:.4f}" if n else "pbs_rate\tNA",
             f"require_pbs\t{int(a.require_pbs)}", f"inserts_out\t{kept}"]
    sys.stderr.write("\n".join(lines) + "\n")
    if a.stats:
        open(a.stats, "w").write("\n".join(lines) + "\n")

if __name__ == "__main__":
    main()
