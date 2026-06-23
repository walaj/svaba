#!/usr/bin/env python3
"""
benchmark.py -- score an SV caller's output against the simulator truth set.

Matches caller breakend-pairs to truth.bedpe within a position tolerance
(orientation-agnostic) and reports recall by SV type and size bin, plus overall
precision. Accepts svaba `bps.txt[.gz]` or a VCF (symbolic <DEL>/<DUP>/<INV>/<INS>
and paired BND records).

Usage:
  benchmark.py --truth truth.bedpe --calls sample.bps.txt.gz [--format auto]
               [--tol 100] [--pass-only] [--somatic-only]

Notes on precision: if the matched normal is a REAL blood BAM, the caller will
also (correctly) emit germline calls that are NOT in the somatic truth. Use
--somatic-only to score only the caller's somatic calls against the somatic
truth, which is the fair comparison for that setup.
"""
import argparse, bisect, gzip, sys
from collections import defaultdict

SIZE_BINS = [(0, 50, "<50bp"), (50, 1000, "50bp-1kb"), (1000, 10000, "1-10kb"),
             (10000, 100000, "10-100kb"), (100000, 1000000, "100kb-1Mb"),
             (1000000, 10**12, ">1Mb")]


def open_any(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def size_bin(svtype, svlen):
    if svtype in ("TRA", "BND"):
        return "interchrom"
    a = abs(svlen)
    for lo, hi, name in SIZE_BINS:
        if lo <= a < hi:
            return name
    return ">1Mb"


# ---------------------------------------------------------------------------
def load_truth(path):
    out = []
    with open_any(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            c1, s1, _, c2, s2, _, name, svtype, svlen, som = f[:10]
            out.append(dict(b1=(c1, int(s1)), b2=(c2, int(s2)), name=name,
                            svtype=svtype, svlen=int(svlen),
                            somatic=(som == "somatic"), matched=False))
    return out


# ---------------------------------------------------------------------------
def load_bps(path, pass_only, somatic_only):
    calls = []
    with open_any(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name.lstrip("#"): i for i, name in enumerate(header)}
        def col(*names):
            for n in names:
                if n in idx:
                    return idx[n]
            return None
        ic1, ip1, ic2, ip2 = col("chr1"), col("pos1"), col("chr2"), col("pos2")
        iconf, isom = col("conf", "confidence"), col("somatic")
        isl, iml, ity = col("somlod"), col("maxlod"), col("type", "evidence")
        icc1, icc2 = col("contig_conf1"), col("contig_conf2")
        get = lambda f, i: f[i] if (i is not None and i < len(f)) else ""
        def ccmin(f):  # min contig_conf over the two breakends ("" if absent)
            try:
                return str(min(float(get(f, icc1)), float(get(f, icc2))))
            except ValueError:
                return ""
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(ip2, iconf or 0):
                continue
            conf = f[iconf] if iconf is not None else "PASS"
            som = f[isom] if isom is not None else ""
            if pass_only and conf != "PASS":
                continue
            if somatic_only and som not in ("1", "SOMATIC", "true", "TRUE"):
                continue
            calls.append(dict(b1=(f[ic1], int(f[ip1])), b2=(f[ic2], int(f[ip2])),
                              conf=conf, somatic=som, evtype=get(f, ity),
                              somlod=get(f, isl), maxlod=get(f, iml), cc=ccmin(f)))
    return calls


def load_vcf(path, pass_only, somatic_only):
    calls = []
    bnd = {}  # id -> record
    with open_any(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, vid, ref, alt, qual, filt, info = f[:8]
            pos = int(pos)
            ikv = {}
            for kv in info.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1); ikv[k] = v
                else:
                    ikv[kv] = True
            is_som = ("SOMATIC" in ikv)
            if pass_only and filt not in ("PASS", "."):
                continue
            svtype = ikv.get("SVTYPE", "")
            if svtype == "BND" or "[" in alt or "]" in alt:
                # parse mate locus from ALT  e.g.  N[chr5:123[  or  ]chr5:123]N
                import re
                m = re.search(r'[\[\]]([^\[\]:]+):(\d+)[\[\]]', alt)
                if not m:
                    continue
                mate = (m.group(1), int(m.group(2)))
                rec = dict(b1=(chrom, pos), b2=mate, conf=filt, somatic=is_som,
                           somlod=ikv.get("SOMLOD", ""), maxlod=ikv.get("MAXLOD", ""),
                           evtype=ikv.get("EVDNC", ""))
                calls.append(rec)  # each BND row already carries both ends
            else:
                end = int(ikv.get("END", pos))
                calls.append(dict(b1=(chrom, pos), b2=(chrom, end),
                                  conf=filt, somatic=is_som,
                                  somlod=ikv.get("SOMLOD", ""), maxlod=ikv.get("MAXLOD", ""),
                                  evtype=ikv.get("EVDNC", "")))
    if somatic_only:
        calls = [c for c in calls if c["somatic"]]
    return calls


# ---------------------------------------------------------------------------
def match(truth, calls, tol):
    # index call breakends by chrom -> sorted [(pos, call_idx, which_end)]
    idx = defaultdict(list)
    for ci, c in enumerate(calls):
        idx[c["b1"][0]].append((c["b1"][1], ci))
        idx[c["b2"][0]].append((c["b2"][1], ci))
    for c in idx:
        idx[c].sort()

    def near(chrom, pos):
        lst = idx.get(chrom)
        if not lst:
            return []
        ps = [p for p, _ in lst]
        lo = bisect.bisect_left(ps, pos - tol)
        hi = bisect.bisect_right(ps, pos + tol)
        return [lst[k][1] for k in range(lo, hi)]

    def bp_match(a, b):
        return a[0] == b[0] and abs(a[1] - b[1]) <= tol

    matched_calls = set()
    call_to_truth = {}              # call idx -> truth name it matched (first wins)
    for t in truth:
        cand = set(near(*t["b1"])) | set(near(*t["b2"]))
        for ci in cand:
            c = calls[ci]
            if (bp_match(t["b1"], c["b1"]) and bp_match(t["b2"], c["b2"])) or \
               (bp_match(t["b1"], c["b2"]) and bp_match(t["b2"], c["b1"])):
                t["matched"] = True
                matched_calls.add(ci)       # credit EVERY matching call (no break)
                call_to_truth.setdefault(ci, t["name"])
    # NB: no break -> a call matching ANY truth is a TP, not an FP. svaba emits
    # BOTH junctions of an inversion (often 1 bp apart) for a single truth row;
    # breaking after the first match wrongly counted the 2nd correct junction FP.
    return matched_calls, call_to_truth


def write_events(path, truth, calls, matched_calls, call_to_truth):
    """Per-CALL event TSV (+ FN truths) so a viewer can re-score at any
    somlod/maxlod cutoff exactly:
      - TP : one row per call that matched a truth -> carries that truth's stable
             id/coords/type + the CALL's somlod/maxlod.
      - FP : one row per call that matched nothing -> the call's somlod/maxlod.
      - FN : one row per truth no call matched (no score).
    At cutoff (sc,mc): recall = |distinct truth ids over passing TP rows| / |truth|;
    precision = |passing TP rows| / (|passing TP| + |passing FP|)."""
    truth_by_name = {t["name"]: t for t in truth}
    with open(path, "w") as fh:
        fh.write("status\ttruth_id\tchrom1\tpos1\tchrom2\tpos2\t"
                 "svtype\tsvlen\tsize_bin\tsomlod\tmaxlod\tevtype\tcontig_conf\n")
        for ci, c in enumerate(calls):
            sl = c.get("somlod", ""); ml = c.get("maxlod", ""); ev = c.get("evtype", ""); cc = c.get("cc", "")
            if ci in matched_calls:
                t = truth_by_name[call_to_truth[ci]]
                fh.write("\t".join([
                    "TP", t["name"], t["b1"][0], str(t["b1"][1]), t["b2"][0], str(t["b2"][1]),
                    t["svtype"], str(t["svlen"]), size_bin(t["svtype"], t["svlen"]),
                    str(sl), str(ml), str(ev), str(cc)]) + "\n")
            else:
                fid = "FP_%s_%d_%s_%d" % (c["b1"][0], c["b1"][1], c["b2"][0], c["b2"][1])
                fh.write("\t".join([
                    "FP", fid, c["b1"][0], str(c["b1"][1]), c["b2"][0], str(c["b2"][1]),
                    "", "", "", str(sl), str(ml), str(ev), str(cc)]) + "\n")
        for t in truth:
            if not t["matched"]:
                fh.write("\t".join([
                    "FN", t["name"], t["b1"][0], str(t["b1"][1]), t["b2"][0], str(t["b2"][1]),
                    t["svtype"], str(t["svlen"]), size_bin(t["svtype"], t["svlen"]),
                    "", "", "", ""]) + "\n")


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True)
    ap.add_argument("--calls", required=True)
    ap.add_argument("--format", choices=["auto", "bps", "vcf"], default="auto")
    ap.add_argument("--tol", type=int, default=100, help="breakend tolerance bp [100]")
    ap.add_argument("--pass-only", action="store_true")
    ap.add_argument("--somatic-only", action="store_true",
                    help="score only somatic truth vs somatic calls")
    ap.add_argument("--events-out", default="",
                    help="write a per-event TSV (TP/FN per truth, FP per call) "
                         "for the cross-commit event matrix viewer")
    args = ap.parse_args()

    fmt = args.format
    if fmt == "auto":
        fmt = "vcf" if (".vcf" in args.calls) else "bps"

    truth = load_truth(args.truth)
    if args.somatic_only:
        truth = [t for t in truth if t["somatic"]]
    calls = (load_bps if fmt == "bps" else load_vcf)(
        args.calls, args.pass_only, args.somatic_only)

    matched_calls, call_to_truth = match(truth, calls, args.tol)

    if args.events_out:
        write_events(args.events_out, truth, calls, matched_calls, call_to_truth)

    # ---- recall by type ----
    by_type = defaultdict(lambda: [0, 0])
    by_bin = defaultdict(lambda: [0, 0])
    for t in truth:
        by_type[t["svtype"]][0] += 1
        by_type[t["svtype"]][1] += int(t["matched"])
        b = size_bin(t["svtype"], t["svlen"])
        by_bin[b][0] += 1
        by_bin[b][1] += int(t["matched"])

    print("== benchmark: %s (fmt=%s, tol=%dbp%s%s) ==" % (
        args.calls.split("/")[-1], fmt, args.tol,
        ", PASS-only" if args.pass_only else "",
        ", somatic-only" if args.somatic_only else ""))
    print("\nRecall by SV type:")
    print("  %-10s %8s %8s %8s" % ("type", "truth", "found", "recall"))
    for tp in sorted(by_type):
        n, k = by_type[tp]
        print("  %-10s %8d %8d %7.1f%%" % (tp, n, k, 100.0 * k / n if n else 0))

    order = [b[2] for b in SIZE_BINS] + ["interchrom"]
    print("\nRecall by size bin:")
    print("  %-12s %8s %8s %8s" % ("bin", "truth", "found", "recall"))
    for b in order:
        if b in by_bin:
            n, k = by_bin[b]
            print("  %-12s %8d %8d %7.1f%%" % (b, n, k, 100.0 * k / n if n else 0))

    ntru = len(truth); nfound = sum(t["matched"] for t in truth)
    ncalls = len(calls); ntp = len(matched_calls)
    print("\nOverall:")
    print("  truth=%d  detected=%d  recall=%.1f%%" % (
        ntru, nfound, 100.0 * nfound / ntru if ntru else 0))
    print("  calls=%d  matching-truth=%d  precision=%.1f%%%s" % (
        ncalls, ntp, 100.0 * ntp / ncalls if ncalls else 0,
        "  (note: real-normal germline calls count as FP unless --somatic-only)"
        if not args.somatic_only else ""))


if __name__ == "__main__":
    main()
