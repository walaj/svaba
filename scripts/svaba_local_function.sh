# svaba_local_function.sh — sourceable bash helpers for common svaba tasks.
#
# Designed to be dot-sourced into an interactive shell (bash/zsh), NOT
# executed as a script. Usage:
#
#   # in your ~/.bashrc or ~/.zshrc (or run once per shell session):
#   source ~/git/svaba/scripts/svaba_local_function.sh
#
#   # then use any of the svaba_* helpers:
#   svaba_help
#   svaba_grep_contig fermi4.r2c.pass.txt.gz c_fermi_chr19_48755001_48780001_13C
#
# Every function is prefixed `svaba_` to stay out of the way of anything
# else on your PATH. Every function responds to `-h` / `--help` with its
# own usage block, so you don't need to remember signatures — just type
# the name and `-h`.
#
# Set `SVABA_FUNCTIONS_QUIET=1` before sourcing to suppress the banner.

# ---------------------------------------------------------------------
# svaba_help — list all svaba_* functions defined here with one-liners.
# ---------------------------------------------------------------------
svaba_help() {
  cat <<'EOF'
svaba_* helper functions (from scripts/svaba_local_function.sh).
Call any function with -h / --help for detailed usage.

  svaba_help
      this message

  svaba_grep_contig FILE CNAME [OUTFILE]
      pull one contig's rows out of an r2c TSV, preserving the header
      so the result loads in viewer/r2c_explorer.html

  svaba_pass_cnames FILE
      print cnames of all PASS variants in a bps.txt[.gz]

  svaba_pass_somatic_cnames FILE
      same, restricted to PASS with somlod >= 1

  svaba_bps_cols
      print the bps.txt.gz column-index reference

  svaba_igv LOCUS [PORT]
      tell a running IGV (default port 60151) to jump to LOCUS

  svaba_swap_bps IN.bps.txt.gz OUT.bps.txt.gz
      sanity-check helper: randomly (50%/row) swap bp1<->bp2 fields
      across all paired columns, so downstream compare code can be
      tested for orientation-invariance

  svaba_bam_grep IN.bam SEQUENCE OUT.bam [REGION]
      pull reads whose SEQ contains SEQUENCE (or its reverse
      complement) into a sub-BAM. Optional REGION restricts the scan
      to a samtools-style locus first (requires indexed IN.bam)
EOF
}

# ---------------------------------------------------------------------
# svaba_grep_contig — the most frequent ask: pull one contig's rows out
# of an r2c TSV, preserving the header row so the sliced file loads in
# viewer/r2c_explorer.html without falling back to the canonical schema.
# ---------------------------------------------------------------------
svaba_grep_contig() {
  if [[ "${1:-}" == "-h" || "${1:-}" == "--help" || $# -lt 2 ]]; then
    cat <<'EOF'
svaba_grep_contig FILE CNAME [OUTFILE]

  Extract the "contig" row plus all "read" rows for a given cname from
  an r2c TSV produced by svaba run --dump-reads (or the PASS / PASS-
  somatic subsets from svaba_postprocess.sh). The header line is always
  preserved so the output loads directly in viewer/r2c_explorer.html.

  FILE      input r2c TSV — .txt.gz (auto-decompressed) or plain .txt
  CNAME     contig name (awk regex; a prefix or unique substring works)
  OUTFILE   output path   (default: CNAME.txt — note: plain text, not .gz)

Examples:
  svaba_grep_contig fermi4_full.r2c.pass.txt.gz c_fermi_chr19_48755001_48780001_13C
  svaba_grep_contig fermi4.r2c.pass.txt.gz c_fermi_chr7  /tmp/chr7_contigs.txt
EOF
    return 1
  fi
  local file="$1" cname="$2" out="${3:-${2}.txt}"
  if [[ ! -r "$file" ]]; then
    echo "svaba_grep_contig: cannot read $file" >&2
    return 2
  fi
  # Decompress if .gz, cat otherwise. `gzip -dc` works on macOS where
  # plain `zcat` expects a .Z suffix rather than .gz.
  local reader
  case "$file" in
    *.gz) reader="gzip -dc" ;;
    *)    reader="cat" ;;
  esac
  $reader "$file" | awk -v cname="$cname" 'NR==1 || $0 ~ cname' > "$out"
  local n
  n=$(wc -l < "$out" | awk '{print $1}')
  echo "svaba_grep_contig: wrote $out  ($n lines incl. header)" >&2
}

# ---------------------------------------------------------------------
# svaba_pass_cnames — one-liner for pulling the cnames of PASS variants
# out of bps.txt[.gz]. Useful for feeding into svaba_grep_contig in a
# loop, or for sanity-checking counts.
# ---------------------------------------------------------------------
svaba_pass_cnames() {
  if [[ "${1:-}" == "-h" || "${1:-}" == "--help" || $# -lt 1 ]]; then
    cat <<'EOF'
svaba_pass_cnames FILE

  Print cnames (bps.txt col 30 = contig_and_region) for every PASS row
  (col 32 == "PASS") in FILE. Output is sorted + unique.

  FILE   ${ID}.bps.txt.gz (or plain .txt)

Example:
  svaba_pass_cnames tumor_normal.bps.txt.gz | wc -l
EOF
    return 1
  fi
  local file="$1" reader
  case "$file" in
    *.gz) reader="gzip -dc" ;;
    *)    reader="cat" ;;
  esac
  $reader "$file" | awk -F'\t' 'NR>1 && $32=="PASS" {print $30}' | LC_ALL=C sort -u
}

# ---------------------------------------------------------------------
# svaba_pass_somatic_cnames — same as above but restricted to the
# PASS + somlod >= 1 subset (i.e. what ends up in r2c.pass.somatic).
# ---------------------------------------------------------------------
svaba_pass_somatic_cnames() {
  if [[ "${1:-}" == "-h" || "${1:-}" == "--help" || $# -lt 1 ]]; then
    cat <<'EOF'
svaba_pass_somatic_cnames FILE

  Print cnames for every PASS + somatic variant (col 32 == "PASS" AND
  col 37 somlod >= 1) in FILE. Non-numeric somlod values ("NA",
  missing) coerce to 0 and are excluded. Output is sorted + unique.

  FILE   ${ID}.bps.txt.gz (or plain .txt)
EOF
    return 1
  fi
  local file="$1" reader
  case "$file" in
    *.gz) reader="gzip -dc" ;;
    *)    reader="cat" ;;
  esac
  $reader "$file" | awk -F'\t' 'NR>1 && $32=="PASS" && $37+0 >= 1 {print $30}' | LC_ALL=C sort -u
}

# ---------------------------------------------------------------------
# svaba_bps_cols — quick printable reference for the bps.txt column
# layout. Useful when you're writing one-off awks and need to remember
# "which column is somlod again?"
# ---------------------------------------------------------------------
svaba_bps_cols() {
  cat <<'EOF'
svaba bps.txt column layout (from BreakPoint::toFileString).
Column numbers are 1-indexed (awk $1 == chr1).

  col  name                   notes
  ---  ---------------------  -----------------------------------------
   1   chr1
   2   pos1
   3   strand1
   4   chr2
   5   pos2
   6   strand2
   7   ref
   8   alt
   9   span
  10   split
  11   alt_count
  12   cov
  13   cigar
  14   cigar_near
  15   dmq1                   discordant-MAPQ, side 1
  16   dmq2                   discordant-MAPQ, side 2
  17   dcn                    discordant count, normal
  18   dct                    discordant count, tumor
  19   mapq1                  contig-alignment MAPQ, side 1
  20   mapq2
  21   nm1
  22   nm2
  23   as1
  24   as2
  25   sub1
  26   sub2
  27   homol
  28   insert
  29   repeat
  30   contig_and_region      <- cname, used as r2c join key
  31   naln
  32   conf                   <- "PASS" / "LOWLOD" / etc.
  33   type                   (older files: named "evidence")
  34   qual
  35   2ndary                 (older files: named "secondary")
  36   somatic
  37   somlod                 <- "somatic LOD"
  38   maxlod                 <- max over per-sample LO
  39   dbsnp
  40   contig_conf1
  41   contig_conf2
  42   cpos1                  v2 only (post-SvABA2.0 refilter round-trip)
  43   cpos2                  v2 only
  44   lmatch                 v2 only
  45   rmatch                 v2 only
  46   scov1                  v2 only
  47   scov2                  v2 only
  48   local1                 v2 only
  49   local2                 v2 only
  50   ctglen                 v2 only
  51   flipped                v2 only
  52   bp_id                  v3 only — unique per-BP identifier,
                              format "bpTTTNNNNNNNN" (thread TTT,
                              counter NNNNNNNN). Cross-references to
                              r2c.txt.gz's split_bps / disc_bps
                              columns so you can tell which BP on a
                              contig a given read actually supports.
                              ("." on very old v2 files rewritten
                              through a v3 emitter that left id unset.)
  53+  per-sample blocks      FORMAT: GT:AD:DP:SR:DR:GQ:PL:LO:LO_n
                              (one block per BAM, order from the header row)
EOF
}

# ---------------------------------------------------------------------
# svaba_igv — nudge a running IGV to a locus via its HTTP command port.
# Requires IGV's "Enable port" preference to be on (default port 60151).
# Useful in loops while eyeballing variants.
# ---------------------------------------------------------------------
svaba_igv() {
  if [[ "${1:-}" == "-h" || "${1:-}" == "--help" || $# -lt 1 ]]; then
    cat <<'EOF'
svaba_igv LOCUS [PORT]

  POST a goto request to a running IGV instance so it navigates to
  LOCUS. IGV must have the "Enable port" preference turned on.

  LOCUS   anything IGV understands: chr:pos, chr:start-end, gene name, ...
  PORT    IGV command port (default 60151)

Examples:
  svaba_igv chr7:55249071
  svaba_igv chr17:7675000-7687000
EOF
    return 1
  fi
  local locus="$1" port="${2:-60151}"
  # Use curl (widely available on macOS + linux). encodeURIComponent-style
  # minimal escaping: swap spaces for %20 so gene/locus queries work.
  local enc
  enc=$(printf '%s' "$locus" | sed 's/ /%20/g')
  if ! curl -sS --max-time 3 "http://localhost:${port}/goto?locus=${enc}" >/dev/null; then
    echo "svaba_igv: IGV not reachable at localhost:${port} (is IGV running with 'Enable port' on?)" >&2
    return 3
  fi
  echo "svaba_igv: IGV -> $locus"
}

# ---------------------------------------------------------------------
# svaba_swap_bps — sanity-check helper. Randomly (50% per row) swap all
# bp1-side <-> bp2-side paired columns in a bps.txt[.gz], producing an
# orientation-randomized file. Downstream compare code that's supposed
# to be orientation-invariant should treat in.bps.txt.gz and its
# swapped version as equivalent.
# ---------------------------------------------------------------------
svaba_swap_bps() {
  if [[ "${1:-}" == "-h" || "${1:-}" == "--help" || $# -lt 2 ]]; then
    cat <<'EOF'
svaba_swap_bps IN.bps.txt[.gz] OUT.bps.txt.gz

  Read IN, preserve the header, then for each data row independently
  flip a fair coin: with prob 1/2, swap every bp1-side column with its
  bp2-side partner. Output is gzip-compressed regardless of the input
  form.

  Paired columns swapped together (so rows stay internally consistent):
     1/4   chr1   <-> chr2
     2/5   pos1   <-> pos2
     3/6   strand1 <-> strand2
    15/16  dmq1   <-> dmq2
    19/20  mapq1  <-> mapq2
    21/22  nm1    <-> nm2
    23/24  as1    <-> as2
    25/26  sub1   <-> sub2
    40/41  contig_conf1 <-> contig_conf2
    42/43  cpos1  <-> cpos2          (v2 only — auto-skipped on v1 files)
    44/45  lmatch <-> rmatch          (v2 only)
    46/47  scov1  <-> scov2           (v2 only)
    48/49  local1 <-> local2          (v2 only)

  Non-paired side-agnostic columns (ref, alt, span, somlod, ...) are
  left alone.

Example:
  svaba_swap_bps fermi4_full.bps.txt.gz fermi4_full.bp_swapped.bps.txt.gz
EOF
    return 1
  fi
  local in="$1" out="$2" reader
  case "$in" in
    *.gz) reader="gzip -dc" ;;
    *)    reader="cat" ;;
  esac
  $reader "$in" | awk 'BEGIN{
      FS=OFS="\t"; srand()
      split("1,4 2,5 3,6 15,16 19,20 21,22 23,24 25,26 40,41 42,43 44,45 46,47 48,49", P, " ")
    }
    /^#/ { print; next }
    {
      if (rand() < 0.5) {
        for (k in P) {
          split(P[k], q, ",")
          a = q[1] + 0; b = q[2] + 0
          if (a <= NF && b <= NF) {
            t = $a; $a = $b; $b = t
          }
        }
      }
      print
    }' | gzip -c > "$out"
  echo "svaba_swap_bps: wrote $out" >&2
}

# ---------------------------------------------------------------------
# svaba_bam_grep — pull reads whose SEQ contains a target sequence (or
# its reverse complement) into a sub-BAM. Useful for fishing for reads
# carrying a specific motif, breakpoint junction, inserted sequence,
# etc. Decodes BAM -> SAM, filters col 10 (SEQ) via awk's index(), then
# re-encodes to BAM.
# ---------------------------------------------------------------------
svaba_bam_grep() {
  if [[ "${1:-}" == "-h" || "${1:-}" == "--help" || $# -lt 3 ]]; then
    cat <<'EOF'
svaba_bam_grep IN.bam SEQUENCE OUT.bam [REGION]

  Filter IN.bam to reads whose SEQ field contains SEQUENCE or its
  reverse complement. Writes OUT.bam (with the original header). Match
  is case-insensitive (SEQUENCE is uppercased, reverse-complemented
  using the standard A<->T / C<->G / N<->N mapping); both strands are
  checked in one pass so you don't need to submit both orientations.

  IN.bam      input BAM (reads will be decompressed + re-encoded;
              coordinate-sorted input yields coordinate-sorted output)
  SEQUENCE    DNA motif to grep for. ACGTN only.
  OUT.bam     output path.
  REGION      optional samtools-style region (e.g. "chr7:55000000-56000000")
              to restrict the scan. Requires IN.bam to be indexed (.bai).

  Depends on `samtools` being on PATH. Scanning is linear in the number
  of reads in the region (full-BAM scan of a 30x WGS takes ~15-30 min
  on a single core; use REGION to localize when you can). Short
  sequences (< ~8 bp) will match nearly every read — a warning is
  printed if you pass one.

Examples:
  # Find reads carrying a specific junction/insert sequence anywhere
  svaba_bam_grep tumor.bam TTAGGGTTAGGGTTAGGG telomere_reads.bam

  # Same, but only scan chr17 to speed things up
  svaba_bam_grep tumor.bam CATCGATCGATC chr17_hits.bam chr17

  # Short sequence — will match lots; use a longer motif if you can
  svaba_bam_grep tumor.bam ACGT lots.bam
EOF
    return 1
  fi
  local in="$1" seq="$2" out="$3" region="${4:-}"
  if [[ ! -r "$in" ]]; then
    echo "svaba_bam_grep: cannot read $in" >&2
    return 2
  fi
  if ! command -v samtools >/dev/null 2>&1; then
    echo "svaba_bam_grep: samtools not on PATH" >&2
    return 3
  fi
  if [[ -z "$seq" ]]; then
    echo "svaba_bam_grep: empty SEQUENCE" >&2
    return 4
  fi
  # Validate: ACGTN only, any case.
  if ! printf '%s' "$seq" | grep -qE '^[ACGTNacgtn]+$'; then
    echo "svaba_bam_grep: SEQUENCE must contain only ACGTN (got: $seq)" >&2
    return 5
  fi
  # Uppercase via tr (works on bash 3.2, unlike ${seq^^}).
  local upper_seq rc_seq
  upper_seq=$(printf '%s' "$seq" | tr '[:lower:]' '[:upper:]')
  # Reverse-complement: complement via tr, then reverse.
  rc_seq=$(printf '%s' "$upper_seq" | tr 'ACGTN' 'TGCAN' | rev)

  if [[ ${#upper_seq} -lt 8 ]]; then
    echo "svaba_bam_grep: warning — sequence length ${#upper_seq} is short; expect many incidental matches" >&2
  fi

  echo "svaba_bam_grep: matching '$upper_seq' (or rc '$rc_seq') in $in${region:+  region=$region}" >&2

  # samtools view -h emits header (@) + records; awk keeps header as-is
  # and data rows only if SEQ (col 10) contains either orientation.
  # index() is a substring search — fast, no regex compile per row.
  # Re-encode via `samtools view -b -o OUT -` so output is a valid BAM.
  # If the user passed a REGION, include it after the input path; that's
  # samtools' native locus-restriction syntax and it only works if the
  # input has a .bai (or .csi) alongside.
  if [[ -n "$region" ]]; then
    samtools view -h "$in" "$region"
  else
    samtools view -h "$in"
  fi | awk -v fwd="$upper_seq" -v rc="$rc_seq" '
      /^@/ { print; next }
      index($10, fwd) || index($10, rc) { print }
    ' | samtools view -b -o "$out" -

  # Report count. Try to build a .bai if the output ended up coord-sorted;
  # silently skip if it didn't (e.g. name-sorted input BAM).
  local n
  n=$(samtools view -c "$out" 2>/dev/null)
  if samtools index "$out" 2>/dev/null; then
    echo "svaba_bam_grep: wrote $out  ($n matching reads, indexed)" >&2
  else
    echo "svaba_bam_grep: wrote $out  ($n matching reads; not coord-sorted, no .bai built)" >&2
  fi
}

# ---------------------------------------------------------------------
# Banner — suppressible with SVABA_FUNCTIONS_QUIET=1. Sent to stderr so
# a non-interactive shell piping output isn't polluted.
# ---------------------------------------------------------------------
if [[ "${SVABA_FUNCTIONS_QUIET:-0}" != "1" ]]; then
  echo "svaba_local_function.sh sourced. Try 'svaba_help' for the list." >&2
fi
