#!/usr/bin/env python3
"""
simulate_sv_genome.py  --  svaba SV simulator (genome construction + truth set)

Generates thousands of NON-OVERLAPPING structural variants spanning the full
size range (small indels -> inter-chromosomal rearrangements), applies them to a
reference genome to build a diploid DONOR genome (hap1 + hap2 FASTA), and emits a
ground-truth set (VCF + BEDPE + TSV) for benchmarking.

Design goals
------------
* Self-contained: Python 3 standard library + `samtools` only (no pysam/numpy).
* Reproducible: single --seed drives all randomness (random.Random).
* Non-overlapping: every SV reserves a padded footprint; placements that would
  overlap an existing footprint, an N-gap, or an --avoid-bed region are rejected.
* Diploid + somatic-aware: each SV is het (one haplotype) or hom (both), and is
  flagged germline (normal+tumor) or somatic (tumor only). Default = all somatic
  (the tumor/normal somatic-caller use case). --germline-frac>0 also builds a
  matched normal donor.
* Correct coordinates: truth breakends are emitted in REFERENCE coordinates
  (reads align back to the reference), directly comparable to a caller's output.

SV models (all edits are LOCAL splices into the reference -> composable, uncapped)
---------------------------------------------------------------------------------
* DEL  affected ref interval [a,b): bases a..b-1 deleted; junction joins a-1 to b.
       Truth breakends: (chrom,a)-(chrom,b).
* INS  at p: L novel random bases inserted between p-1 and p. Breakends: (p)-(p).
* DUP  tandem of [a,b): a copy is placed immediately after b. Breakends: (a)-(b).
* INV  [a,b) reverse-complemented in place. Breakends: (a)-(b).
* TRA  inter-chromosomal templated insertion (unbalanced translocation): a segment
       chrB[q:q+L] (forward or reverse-complemented) is spliced into chrA at p.
       This yields TWO inter-chromosomal junctions, emitted as two BND truth rows:
         proximal  (chrA,p)-(chrB,q)      distal  (chrA,p)-(chrB,q+L)
       Reads fully inside the inserted segment also add coverage at chrB[q:q+L]
       (the donor locus), exactly like a real templated insertion. Uncapped:
       a chromosome may host/donate many translocations.

Use sim/benchmark.py to score a caller against truth with a position tolerance
(default 100 bp, orientation-agnostic), so breakend-convention differences between
callers don't matter.
"""

import argparse
import bisect
import math
import os
import random
import subprocess
import sys
from collections import defaultdict

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def log(msg):
    sys.stderr.write(msg + "\n"); sys.stderr.flush()


def read_fai(fai_path):
    out = []
    with open(fai_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out.append((f[0], int(f[1])))
    return out


def fetch_chrom(ref, chrom):
    p = subprocess.run(["samtools", "faidx", ref, chrom],
                       capture_output=True, text=True, check=True)
    return "".join(p.stdout.split("\n")[1:]).upper()


def fetch_region(ref, chrom, start, end):
    """0-based half-open [start,end) -> uppercase sequence (samtools is 1-based incl)."""
    reg = "%s:%d-%d" % (chrom, start + 1, end)
    p = subprocess.run(["samtools", "faidx", ref, reg],
                       capture_output=True, text=True, check=True)
    return "".join(p.stdout.split("\n")[1:]).upper()


def callable_intervals(seq, min_len):
    out = []; n = len(seq); i = 0
    while i < n:
        if seq[i] == "N":
            i += 1; continue
        j = i
        while j < n and seq[j] != "N":
            j += 1
        if j - i >= min_len:
            out.append((i, j))
        i = j
    return out


def load_avoid_bed(path):
    av = defaultdict(list)
    if not path:
        return av
    opener = open
    if path.endswith(".gz"):
        import gzip; opener = gzip.open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                av[f[0]].append((int(f[1]), int(f[2])))
            except ValueError:
                continue
    for c in av:
        av[c].sort()
    return av


def overlaps_sorted(intervals, lo, hi):
    if not intervals:
        return False
    starts = [s for s, _ in intervals]
    idx = bisect.bisect_right(starts, hi)
    for k in range(max(0, idx - 1), min(len(intervals), idx + 1)):
        s, e = intervals[k]
        if s < hi and lo < e:
            return True
    return False


def logu(rng, lo, hi):
    if lo >= hi:
        return lo
    return int(math.exp(rng.uniform(math.log(lo), math.log(hi))))


def rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


class SV:
    __slots__ = ("idx", "svtype", "chrom", "a", "b", "chrom2", "q", "L", "rev",
                 "svlen", "novseq", "haps", "somatic", "name")

    def __init__(self, idx, svtype, chrom, a, b, svlen,
                 chrom2=None, q=None, L=0, rev=False,
                 novseq="", haps=(1,), somatic=True):
        self.idx = idx
        self.svtype = svtype          # DEL DUP INV INS TRA
        self.chrom = chrom
        self.a = a                    # primary ref coord / interval start
        self.b = b                    # interval end (== a for INS/TRA insert point)
        self.chrom2 = chrom2          # TRA donor chrom
        self.q = q                    # TRA donor segment start (0-based)
        self.L = L                    # TRA donor segment length
        self.rev = rev                # TRA segment reverse-complemented
        self.svlen = svlen
        self.novseq = novseq          # INS novel bases (small INS only)
        self.haps = tuple(haps)
        self.somatic = somatic
        self.name = "SV%07d_%s" % (idx, svtype)


def plan_svs(args, chrom_lengths):
    rng = random.Random(args.seed)
    ref = args.ref; target = args.chroms; pad = args.min_gap

    log("[plan] computing callable (non-N) intervals for %d chroms ..." % len(target))
    callable_by_chrom = {}
    for c in target:
        seq = fetch_chrom(ref, c)
        callable_by_chrom[c] = callable_intervals(seq, min_len=args.min_callable)
        del seq
    avoid = load_avoid_bed(args.avoid_bed)

    chrom_weight = {c: sum(e - s for s, e in callable_by_chrom[c]) for c in target}
    tot = sum(chrom_weight.values())
    if tot == 0:
        log("ERROR: no callable sequence in target chroms"); sys.exit(1)
    weighted = [c for c in target
                for _ in range(max(1, int(200 * chrom_weight[c] / tot)))]

    reserved = defaultdict(list)
    svs = []; idx = 0

    def reserve(chrom, lo, hi):
        bisect.insort(reserved[chrom], (lo, hi))

    def free_region(chrom, span):
        ivs = callable_by_chrom[chrom]
        if not ivs:
            return None
        for _ in range(args.max_attempts):
            s, e = rng.choice(ivs)
            need = span + 2 * pad
            if e - s <= need:
                continue
            a = rng.randint(s + pad, e - span - pad)
            lo, hi = a - pad, a + span + pad
            if overlaps_sorted(reserved[chrom], lo, hi):
                continue
            if overlaps_sorted(avoid.get(chrom, []), lo, hi):
                continue
            return a
        return None

    def assign_gt():
        somatic = (rng.random() >= args.germline_frac)
        haps = (1, 2) if rng.random() < args.hom_frac else (rng.choice([1, 2]),)
        return haps, somatic

    plan = [
        ("DEL", args.n_del_small, args.del_small_min, args.del_small_max),
        ("INS", args.n_ins_small, args.ins_small_min, args.ins_small_max),
        ("DEL", args.n_del,       args.del_min,       args.del_max),
        ("DUP", args.n_dup,       args.dup_min,       args.dup_max),
        ("INV", args.n_inv,       args.inv_min,       args.inv_max),
    ]
    for svtype, n, smin, smax in plan:
        placed = 0
        for _ in range(n):
            size = logu(rng, smin, smax)
            chrom = rng.choice(weighted)
            span = 1 if svtype == "INS" else size
            a = free_region(chrom, span)
            if a is None:
                continue
            haps, somatic = assign_gt()
            if svtype == "DEL":
                sv = SV(idx, "DEL", chrom, a, a + size, -size, haps=haps, somatic=somatic)
            elif svtype == "DUP":
                sv = SV(idx, "DUP", chrom, a, a + size, size, haps=haps, somatic=somatic)
            elif svtype == "INV":
                sv = SV(idx, "INV", chrom, a, a + size, size, haps=haps, somatic=somatic)
            else:
                sv = SV(idx, "INS", chrom, a, a, size,
                        novseq=rand_seq(rng, size), haps=haps, somatic=somatic)
            reserve(chrom, a - pad, a + max(span, 1) + pad)
            svs.append(sv); idx += 1; placed += 1
        log("[plan] %s size[%d,%d]: placed %d/%d" % (svtype, smin, smax, placed, n))

    # inter-chromosomal templated insertions (uncapped)
    placed = 0
    for _ in range(args.n_tra):
        if len(target) < 2:
            break
        haps, somatic = assign_gt()
        L = logu(rng, args.tra_min, args.tra_max)
        cA = rng.choice(weighted)
        cB = rng.choice(weighted)
        if cB == cA:
            continue
        p = free_region(cA, 1)            # insertion point on chrA
        q = free_region(cB, L)            # donor segment [q,q+L) on chrB
        if p is None or q is None:
            continue
        sv = SV(idx, "TRA", cA, p, p, L, chrom2=cB, q=q, L=L,
                rev=(rng.random() < 0.5), haps=haps, somatic=somatic)
        reserve(cA, p - pad, p + pad)
        reserve(cB, q - pad, q + L + pad)
        svs.append(sv); idx += 1; placed += 1
    log("[plan] TRA (inter-chrom templated insertion): placed %d/%d" % (placed, args.n_tra))

    log("[plan] total SVs placed: %d  (somatic=%d germline=%d)" % (
        len(svs), sum(s.somatic for s in svs), sum(not s.somatic for s in svs)))
    return svs


def build_chrom_seq(seq, edits):
    """Apply sorted non-overlapping edits to `seq`. edits: (kind,a,b,nov)."""
    edits = sorted(edits, key=lambda e: e[1])
    out = []; pos = 0
    for kind, a, b, nov in edits:
        out.append(seq[pos:a])
        if kind == "DEL":
            pos = b
        elif kind == "INS":      # also used for TRA templated insertion
            out.append(nov); pos = a
        elif kind == "DUP":
            out.append(seq[a:b]); out.append(seq[a:b]); pos = b
        elif kind == "INV":
            out.append(revcomp(seq[a:b])); pos = b
    out.append(seq[pos:])
    return "".join(out)


def write_fasta_record(fh, name, seq, width=60):
    fh.write(">%s\n" % name)
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width]); fh.write("\n")


def build_haplotype(ref, target, svs, hap, present, out_fa):
    """Write one haplotype FASTA. `present(sv)` gates SV membership (e.g.
    germline-only for the normal donor)."""
    intra = defaultdict(list)
    for sv in svs:
        if hap not in sv.haps or not present(sv):
            continue
        if sv.svtype == "TRA":
            seg = fetch_region(ref, sv.chrom2, sv.q, sv.q + sv.L)
            if sv.rev:
                seg = revcomp(seg)
            intra[sv.chrom].append(("INS", sv.a, sv.a, seg))
        else:
            intra[sv.chrom].append((sv.svtype, sv.a, sv.b, sv.novseq))

    with open(out_fa, "w") as fh:
        for c in target:
            if c in intra:
                seq = fetch_chrom(ref, c)
                es = build_chrom_seq(seq, intra[c])
                del seq
            else:
                es = fetch_chrom(ref, c)
            write_fasta_record(fh, "%s_h%d" % (c, hap), es)
            del es
    log("[build] hap%d -> %s  (%d edited chroms)" % (
        hap, os.path.basename(out_fa), len(intra)))


def emit_truth(svs, chrom_lengths, args):
    od = args.out_dir

    with open(os.path.join(od, "truth.bedpe"), "w") as fh:
        fh.write("#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tsvtype\t"
                 "svlen\tsomatic\tgt\thaps\n")
        for s in svs:
            gt = "hom" if s.haps == (1, 2) else "het"
            som = "somatic" if s.somatic else "germline"
            haps = ",".join(map(str, s.haps))
            if s.svtype == "TRA":
                # two inter-chromosomal junctions (proximal + distal)
                for tag, qpos in (("prox", s.q), ("dist", s.q + s.L)):
                    fh.write("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%s\t%s\t%s\n" % (
                        s.chrom, s.a, s.a + 1, s.chrom2, qpos, qpos + 1,
                        s.name + "_" + tag, "TRA", s.L, som, gt, haps))
            else:
                fh.write("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%s\t%s\t%s\n" % (
                    s.chrom, s.a, s.a + 1, s.chrom, s.b, s.b + 1,
                    s.name, s.svtype, s.svlen, som, gt, haps))

    with open(os.path.join(od, "truth.sv_table.tsv"), "w") as fh:
        fh.write("name\tsvtype\tchrom\ta\tb\tchrom2\tq\tL\trev\tsvlen\tsomatic\tgt\thaps\n")
        for s in svs:
            gt = "hom" if s.haps == (1, 2) else "het"
            fh.write("%s\t%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\n" % (
                s.name, s.svtype, s.chrom, s.a, s.b, s.chrom2 or ".",
                str(s.q) if s.q is not None else ".", s.L, int(s.rev),
                s.svlen, int(s.somatic), gt, ",".join(map(str, s.haps))))

    with open(os.path.join(od, "truth.vcf"), "w") as fh:
        fh.write("##fileformat=VCFv4.2\n##source=svaba_sim/simulate_sv_genome.py\n")
        for c, L in chrom_lengths:
            fh.write("##contig=<ID=%s,length=%d>\n" % (c, L))
        for a, d in (("DEL", "Deletion"), ("DUP", "Tandem duplication"),
                     ("INV", "Inversion"), ("INS", "Insertion")):
            fh.write('##ALT=<ID=%s,Description="%s">\n' % (a, d))
        for line in (
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
            '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">',
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">',
            '##INFO=<ID=CHR2,Number=1,Type=String,Description="Partner chrom">',
            '##INFO=<ID=MATEID,Number=1,Type=String,Description="Mate BND id">',
            '##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic event">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'):
            fh.write(line + "\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDONOR\n")
        for s in sorted(svs, key=lambda s: (s.chrom, s.a)):
            gt = "1|1" if s.haps == (1, 2) else ("1|0" if s.haps == (1,) else "0|1")
            som = ";SOMATIC" if s.somatic else ""
            if s.svtype == "TRA":
                for tag, qpos in (("prox", s.q), ("dist", s.q + s.L)):
                    a_id = s.name + "_" + tag + "_A"
                    b_id = s.name + "_" + tag + "_B"
                    fh.write("%s\t%d\t%s\tN\tN[%s:%d[\t.\tPASS\tSVTYPE=BND;CHR2=%s;MATEID=%s%s\tGT\t%s\n" % (
                        s.chrom, s.a + 1, a_id, s.chrom2, qpos + 1, s.chrom2, b_id, som, gt))
                    fh.write("%s\t%d\t%s\tN\t]%s:%d]N\t.\tPASS\tSVTYPE=BND;CHR2=%s;MATEID=%s%s\tGT\t%s\n" % (
                        s.chrom2, qpos + 1, b_id, s.chrom, s.a + 1, s.chrom, a_id, som, gt))
            else:
                end = s.b if s.svtype != "INS" else s.a
                fh.write("%s\t%d\t%s\tN\t<%s>\t.\tPASS\tSVTYPE=%s;SVLEN=%d;END=%d%s\tGT\t%s\n" % (
                    s.chrom, s.a + 1, s.name, s.svtype, s.svtype, s.svlen, end + 1, som, gt))
    log("[truth] wrote truth.vcf / truth.bedpe / truth.sv_table.tsv to %s" % od)


def main():
    ap = argparse.ArgumentParser(
        description="Simulate non-overlapping SVs and build a diploid donor genome + truth set.")
    ap.add_argument("--ref", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--chroms", default="chr1-22,chrX,chrY")
    ap.add_argument("--avoid-bed", default="")
    ap.add_argument("--min-gap", type=int, default=50000)
    ap.add_argument("--min-callable", type=int, default=200000)
    ap.add_argument("--max-attempts", type=int, default=200)
    ap.add_argument("--hom-frac", type=float, default=0.3)
    ap.add_argument("--germline-frac", type=float, default=0.0)
    ap.add_argument("--n-del-small", type=int, default=800)
    ap.add_argument("--n-ins-small", type=int, default=800)
    ap.add_argument("--n-del", type=int, default=600)
    ap.add_argument("--n-dup", type=int, default=400)
    ap.add_argument("--n-inv", type=int, default=400)
    ap.add_argument("--n-tra", type=int, default=200)
    ap.add_argument("--del-small-min", type=int, default=1)
    ap.add_argument("--del-small-max", type=int, default=50)
    ap.add_argument("--ins-small-min", type=int, default=1)
    ap.add_argument("--ins-small-max", type=int, default=50)
    ap.add_argument("--del-min", type=int, default=50)
    ap.add_argument("--del-max", type=int, default=3000000)
    ap.add_argument("--dup-min", type=int, default=200)
    ap.add_argument("--dup-max", type=int, default=1000000)
    ap.add_argument("--inv-min", type=int, default=200)
    ap.add_argument("--inv-max", type=int, default=1000000)
    ap.add_argument("--tra-min", type=int, default=1000, help="inter-chrom segment min bp")
    ap.add_argument("--tra-max", type=int, default=500000, help="inter-chrom segment max bp")
    ap.add_argument("--no-normal", action="store_true")
    args = ap.parse_args()

    chroms = []
    for tok in args.chroms.split(","):
        tok = tok.strip()
        if not tok:
            continue
        if "-" in tok and tok.lower().startswith("chr") and tok[3:].split("-")[0].isdigit():
            lo, hi = tok[3:].split("-")
            chroms += ["chr%d" % i for i in range(int(lo), int(hi) + 1)]
        else:
            chroms.append(tok)
    args.chroms = chroms

    if not os.path.exists(args.ref + ".fai"):
        log("ERROR: %s.fai not found (samtools faidx %s)" % (args.ref, args.ref)); sys.exit(1)
    os.makedirs(args.out_dir, exist_ok=True)

    fai = read_fai(args.ref + ".fai")
    present = set(c for c, _ in fai)
    missing = [c for c in args.chroms if c not in present]
    if missing:
        log("ERROR: chroms not in reference: %s" % ",".join(missing)); sys.exit(1)
    fai_d = dict(fai)
    chrom_lengths = [(c, fai_d[c]) for c in args.chroms]

    log("[main] ref=%s chroms=%d seed=%d min_gap=%d" % (
        args.ref, len(args.chroms), args.seed, args.min_gap))
    svs = plan_svs(args, chrom_lengths)

    log("[build] constructing TUMOR donor (germline + somatic) ...")
    build_haplotype(args.ref, args.chroms, svs, 1, lambda s: True,
                    os.path.join(args.out_dir, "tumor.hap1.fa"))
    build_haplotype(args.ref, args.chroms, svs, 2, lambda s: True,
                    os.path.join(args.out_dir, "tumor.hap2.fa"))

    n_germ = sum(1 for s in svs if not s.somatic)
    if n_germ > 0 and not args.no_normal:
        log("[build] constructing NORMAL donor (germline only, %d SVs) ..." % n_germ)
        build_haplotype(args.ref, args.chroms, svs, 1, lambda s: not s.somatic,
                        os.path.join(args.out_dir, "normal.hap1.fa"))
        build_haplotype(args.ref, args.chroms, svs, 2, lambda s: not s.somatic,
                        os.path.join(args.out_dir, "normal.hap2.fa"))
    else:
        log("[build] no germline SVs -> normal donor not written "
            "(use a real/clean normal BAM for the T/N pair)")

    emit_truth(svs, chrom_lengths, args)
    log("[done] donor genome(s) + truth written to %s" % args.out_dir)


if __name__ == "__main__":
    main()
