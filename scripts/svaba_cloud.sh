#!/usr/bin/env bash
#
# svaba_cloud.sh — scatter svaba across GCP VMs, one partition per VM.
#
# Partitions the genome into roughly equal chunks (by base-pair count),
# spins up one c3d VM per partition with a shared read-only data disk,
# runs svaba on each partition in parallel, collects the per-partition
# bps.txt.gz outputs to a GCS bucket, and (optionally) merges +
# postprocesses them into a single call set on the local machine.
#
# The BAMs, reference, and svaba binary + blacklist should already live
# on a GCP persistent disk (see --data-disk). That disk is attached
# read-only to every worker, so all VMs share the same data without
# copying.
#
# Prerequisites:
#   - gcloud CLI authenticated with compute + storage permissions
#   - A GCP persistent disk with: BAM(s), reference (.fa + .fai + .bwt*),
#     svaba binary, blacklist BED. Create one manually or use --prep-disk.
#   - A GCS bucket for collecting outputs.
#   - jemalloc installed on the worker image (apt install libjemalloc2).
#     The default Debian/Ubuntu GCP images have it available; the startup
#     script installs it automatically.
#
# Usage:
#   svaba_cloud.sh [options]
#
# Required:
#   --tumor PATH          path to tumor BAM *on the data disk*
#   --normal PATH         path to normal BAM *on the data disk*
#   --ref PATH            path to reference FASTA *on the data disk*
#   --svaba PATH          path to svaba binary *on the data disk*
#   --blacklist PATH      path to blacklist BED *on the data disk*
#
#   Paths are relative to the data disk's filesystem root. The script
#   strips any leading / for safety.  If the data disk was formatted
#   with files at its root (e.g. /tumor.bam), pass --tumor tumor.bam.
#   If the files are in subdirectories (e.g. /data/bams/tumor.bam),
#   pass --tumor data/bams/tumor.bam.
#   --data-disk NAME      name of the GCP persistent disk
#   --bucket gs://...     GCS bucket (or bucket/prefix) for outputs
#   --id STR              analysis ID prefix
#
# Optional:
#   --machine TYPE        GCP machine type           (default: c3d-highcpu-30)
#   --zone ZONE           GCP zone                   (default: us-central1-a)
#   --threads N           svaba -p per VM            (default: 14)
#   --partitions N        number of VMs / partitions (default: 6)
#   --mount-point DIR     where data disk mounts     (default: /mnt/data)
#   --extra-args STR      extra args passed to svaba run (quoted)
#   --jemalloc            LD_PRELOAD jemalloc on workers (default: on)
#   --no-jemalloc         disable jemalloc preload
#   --merge               after all workers finish, merge + postprocess locally
#   --keep-vms            don't delete worker VMs on completion
#   --boot-disk-size STR  boot disk size             (default: 50GB)
#   --dry-run             print commands without executing
#   -h, --help            this message
#
# Genome partitions (hg38, roughly equal Mb per partition):
#   Default 6-way split balances by total base pairs so each VM does
#   similar wall-clock work. Override with --regions-file FILE (one
#   comma-separated region list per line, one line per partition).
#
# Examples:
#   # Full WGS, 6 VMs, merge at the end
#   svaba_cloud.sh \
#     --tumor tumor.bam --normal normal.bam \
#     --ref Homo_sapiens_assembly38.fasta \
#     --svaba svaba --blacklist hg38.combined_blacklist.bed \
#     --data-disk svaba-data --bucket gs://my-bucket/svaba-run \
#     --id my_run --merge
#
#   # Dry run to see what would happen
#   svaba_cloud.sh \
#     --tumor tumor.bam --normal normal.bam \
#     --ref ref.fa --svaba svaba \
#     --blacklist blacklist.bed --data-disk svaba-data \
#     --bucket gs://my-bucket/run1 --id run1 --dry-run
#
#   # Custom regions file (3 partitions)
#   svaba_cloud.sh ... --partitions 3 --regions-file my_regions.txt
#

set -euo pipefail

# ----------------------------------------------------------- defaults ---
MACHINE="c3d-highcpu-30"
ZONE="us-central1-a"
THREADS=14
PARTITIONS=6
MOUNT="/mnt/data"
BOOT_DISK_SIZE="10GB"
IMAGE="svaba-worker-image"
USE_JEMALLOC=1
DO_MERGE=0
KEEP_VMS=0
DRY_RUN=0
EXTRA_ARGS=""
REGIONS_FILE=""

# Required (must be set by flags)
TUMOR=""
NORMAL=""
REF=""
SVABA_BIN=""
BLACKLIST=""
DATA_DISK=""
BUCKET=""
ID=""

# ---------------------------------------------------------------- help ---
print_usage() {
  sed -n '2,70p' "$0" | sed 's/^# \{0,1\}//'
}

# -------------------------------------------------------------- parse ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    --tumor)          TUMOR="$2";          shift 2 ;;
    --normal)         NORMAL="$2";         shift 2 ;;
    --ref)            REF="$2";            shift 2 ;;
    --svaba)          SVABA_BIN="$2";      shift 2 ;;
    --blacklist)      BLACKLIST="$2";      shift 2 ;;
    --data-disk)      DATA_DISK="$2";      shift 2 ;;
    --bucket)         BUCKET="$2";         shift 2 ;;
    --id)             ID="$2";            shift 2 ;;
    --machine)        MACHINE="$2";        shift 2 ;;
    --zone)           ZONE="$2";           shift 2 ;;
    --threads)        THREADS="$2";        shift 2 ;;
    --partitions)     PARTITIONS="$2";     shift 2 ;;
    --mount-point)    MOUNT="$2";          shift 2 ;;
    --extra-args)     EXTRA_ARGS="$2";     shift 2 ;;
    --regions-file)   REGIONS_FILE="$2";   shift 2 ;;
    --jemalloc)       USE_JEMALLOC=1;      shift ;;
    --no-jemalloc)    USE_JEMALLOC=0;      shift ;;
    --merge)          DO_MERGE=1;          shift ;;
    --keep-vms)       KEEP_VMS=1;          shift ;;
    --boot-disk-size) BOOT_DISK_SIZE="$2"; shift 2 ;;
    --image)          IMAGE="$2";          shift 2 ;;
    --dry-run)        DRY_RUN=1;           shift ;;
    -h|--help)        print_usage;         exit 0 ;;
    --)               shift; break ;;
    -*)
      echo "svaba_cloud.sh: unknown option '$1'" >&2
      echo "run 'svaba_cloud.sh --help' for usage" >&2
      exit 2
      ;;
    *)  break ;;
  esac
done

# ------------------------------------------------------- validate args ---
missing=()
[[ -z "$TUMOR" ]]     && missing+=("--tumor")
[[ -z "$REF" ]]       && missing+=("--ref")
[[ -z "$SVABA_BIN" ]] && missing+=("--svaba")
[[ -z "$BLACKLIST" ]] && missing+=("--blacklist")
[[ -z "$DATA_DISK" ]] && missing+=("--data-disk")
[[ -z "$BUCKET" ]]    && missing+=("--bucket")
[[ -z "$ID" ]]        && missing+=("--id")

if [[ ${#missing[@]} -gt 0 ]]; then
  echo "svaba_cloud.sh: missing required flags: ${missing[*]}" >&2
  exit 2
fi

BUCKET=${BUCKET%/}  # strip trailing slash

# Sanitize ID for GCP resource names: lowercase, replace _ with -, strip
# anything that isn't [a-z0-9-].
GCP_ID=$(echo "$ID" | tr '[:upper:]' '[:lower:]' | tr '_' '-' | sed 's/[^a-z0-9-]//g')
if [[ "$GCP_ID" != "$ID" ]]; then
  echo "note: sanitized ID for GCP names: '$ID' -> '$GCP_ID'"
fi

# Strip leading slashes from disk-relative paths. Users will naturally
# pass absolute paths from their local machine (e.g. /mnt/ssd/tumor.bam)
# but on the worker these are relative to the data-disk mount point.
# SVABA_BIN and BLACKLIST are NOT stripped — they may live on the boot
# disk (baked into the worker image) rather than the data disk.
strip_leading_slash() { echo "${1#/}"; }
TUMOR=$(strip_leading_slash "$TUMOR")
NORMAL=$(strip_leading_slash "$NORMAL")
REF=$(strip_leading_slash "$REF")

# ---------------------------------------- default hg38 6-way partition ---
# Balanced by approximate Mb so each VM does similar wall-clock work.
# Total autosomal + X ≈ 3,088 Mb → ~515 Mb per partition.
default_regions() {
  local n=$1
  case $n in
    1)
      echo "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX"
      ;;
    2)
      echo "chr1,chr2,chr3,chr4,chr5,chr6,chr7"
      echo "chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX"
      ;;
    3)
      echo "chr1,chr2,chr3"
      echo "chr4,chr5,chr6,chr7,chr8"
      echo "chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX"
      ;;
    4)
      echo "chr1,chr2"
      echo "chr3,chr4,chr5"
      echo "chr6,chr7,chr8,chr9,chr10"
      echo "chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX"
      ;;
    5)
      echo "chr1,chr2"
      echo "chr3,chr4,chr5"
      echo "chr6,chr7,chr8,chr9"
      echo "chr10,chr11,chr12,chr13,chr14"
      echo "chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX"
      ;;
    6)
      echo "chr1,chr2"
      echo "chr3,chr4,chr5"
      echo "chr6,chr7,chr8,chr9"
      echo "chr10,chr11,chr12,chr13"
      echo "chr14,chr15,chr16,chr17,chr18"
      echo "chr19,chr20,chr21,chr22,chrX"
      ;;
    *)
      # For >6 partitions, fall back to one chromosome per partition
      # (user should supply --regions-file for fine-grained control)
      for c in chr{1..22} chrX; do
        echo "$c"
      done | head -n "$n"
      ;;
  esac
}

# Load region partitions into an array.
declare -a REGIONS
if [[ -n "$REGIONS_FILE" ]]; then
  mapfile -t REGIONS < "$REGIONS_FILE"
  PARTITIONS=${#REGIONS[@]}
else
  mapfile -t REGIONS < <(default_regions "$PARTITIONS")
fi

if [[ ${#REGIONS[@]} -ne $PARTITIONS ]]; then
  echo "svaba_cloud.sh: got ${#REGIONS[@]} region lines but expected $PARTITIONS" >&2
  exit 2
fi

# ------------------------------------------------------- jemalloc preamble ---
PRELOAD=""
if [[ $USE_JEMALLOC -eq 1 ]]; then
  PRELOAD="LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so.2"
fi

# ------------------------------------------------------- normal flag ---
NORMAL_FLAG=""
if [[ -n "$NORMAL" ]]; then
  NORMAL_FLAG="-n ${MOUNT}/${NORMAL}"
fi

# ----------------------------------------------------------- run/exec ---
run_cmd() {
  if [[ $DRY_RUN -eq 1 ]]; then
    echo "[dry-run] $*"
  else
    "$@"
  fi
}

# ================================================================
# STEP 1: Create worker VMs in parallel
# ================================================================
echo "============================================================"
echo "svaba_cloud.sh: launching $PARTITIONS workers"
echo "  machine=$MACHINE zone=$ZONE threads=$THREADS"
echo "  data_disk=$DATA_DISK bucket=$BUCKET id=$ID"
echo "  jemalloc=$USE_JEMALLOC merge=$DO_MERGE"
echo ""
echo "  resolved paths (on worker VMs):"
echo "    tumor:     ${MOUNT}/${TUMOR}  (data disk)"
[[ -n "$NORMAL" ]] && \
echo "    normal:    ${MOUNT}/${NORMAL}  (data disk)"
echo "    ref:       ${MOUNT}/${REF}  (data disk)"
echo "    svaba:     ${SVABA_BIN}  (boot image)"
echo "    blacklist: ${BLACKLIST}  (boot image)"
echo "============================================================"

# Temp dir for startup scripts (--metadata-from-file needs real files;
# --metadata=startup-script= chokes on commas in the script body because
# gcloud interprets them as key=value separators).
SCRIPT_TMPDIR=$(mktemp -d)
trap 'rm -rf "$SCRIPT_TMPDIR"' EXIT

VM_NAMES=()
for i in $(seq 1 "$PARTITIONS"); do
  vm="svaba-${GCP_ID}-worker-${i}"
  VM_NAMES+=("$vm")

  # Write startup script to a temp file so we can use --metadata-from-file
  # (avoids gcloud's comma-parsing issues with --metadata=).
  SCRIPT_FILE="${SCRIPT_TMPDIR}/startup_${i}.sh"
  cat > "$SCRIPT_FILE" <<STARTUP_EOF
#!/bin/bash
set -euo pipefail

# --- mount data disk ---
mkdir -p ${MOUNT} /mnt/output
# The data disk is the first non-boot disk (sdb on most GCP images).
# Wait briefly for the device to appear.
for attempt in \$(seq 1 30); do
  [[ -b /dev/sdb ]] && break
  sleep 2
done
if [[ ! -b /dev/sdb ]]; then
  echo "FATAL: /dev/sdb not found after 60s" >&2
  exit 1
fi
# Try partition first (sdb1), fall back to whole-disk (sdb).
# Most GCP persistent disks have a partition table; raw-formatted
# disks don't. Use noload to skip ext4 journal recovery (required
# when the source disk was detached without a clean unmount, which
# is the normal case for a read-only shared disk).
MOUNT_DEV=/dev/sdb
[[ -b /dev/sdb1 ]] && MOUNT_DEV=/dev/sdb1
mount -o ro,noload \${MOUNT_DEV} ${MOUNT} || mount -o ro \${MOUNT_DEV} ${MOUNT}

# --- install runtime dependencies (only if no custom image) ---
if ! ldconfig -p | grep -q libhts; then
  apt-get update -qq
  PKGS=(libhts3 libhtscodecs2)
  [[ ${USE_JEMALLOC} -eq 1 ]] && PKGS+=(libjemalloc2)
  apt-get install -y -qq "\${PKGS[@]}" >/dev/null 2>&1 || true
fi

# --- run svaba ---
PART_ID="${ID}_part${i}"
cd /mnt/output

${PRELOAD} \\
  ${SVABA_BIN} run \\
    -t ${MOUNT}/${TUMOR} \\
    ${NORMAL_FLAG} \\
    -G ${MOUNT}/${REF} \\
    -p ${THREADS} \\
    -k ${REGIONS[$((i-1))]} \\
    -a \${PART_ID} \\
    --blacklist ${BLACKLIST} \\
    ${EXTRA_ARGS} \\
  2>&1 | tee \${PART_ID}.log

# --- upload results ---
gsutil -m cp \\
  \${PART_ID}.bps.txt.gz \\
  \${PART_ID}.log \\
  \${PART_ID}.contigs.bam \\
  \${PART_ID}.runtime.txt \\
  ${BUCKET}/

# Upload optional outputs if they exist (--dump-reads)
for f in \${PART_ID}.discordant.bam \${PART_ID}.corrected.bam \${PART_ID}.r2c.txt.gz; do
  [[ -f "\$f" ]] && gsutil cp "\$f" ${BUCKET}/ || true
done

# Signal completion
echo "DONE" | gsutil cp - ${BUCKET}/.done_part${i}
STARTUP_EOF

  echo "[${i}/${PARTITIONS}] creating $vm  regions=${REGIONS[$((i-1))]}"
  run_cmd gcloud compute instances create "$vm" \
    --zone="$ZONE" \
    --machine-type="$MACHINE" \
    --image="$IMAGE" \
    --disk="name=${DATA_DISK},mode=ro,device-name=svaba-data" \
    --boot-disk-size="$BOOT_DISK_SIZE" \
    --scopes=storage-rw \
    --metadata-from-file=startup-script="$SCRIPT_FILE" \
    --no-address \
    &
done
wait
echo "all $PARTITIONS VMs created"

# ================================================================
# STEP 2: Wait for all workers to finish
# ================================================================
echo "waiting for workers to complete..."
echo "(checking ${BUCKET}/.done_part* every 60s)"

poll_done() {
  local expected=$1
  while true; do
    local count
    count=$(gsutil ls "${BUCKET}/.done_part*" 2>/dev/null | grep -c '\.done_part' || true)
    if [[ $count -ge $expected ]]; then
      echo "all $expected partitions complete"
      return 0
    fi
    echo "  $(date +%H:%M:%S)  $count / $expected done"
    sleep 60
  done
}

if [[ $DRY_RUN -eq 0 ]]; then
  poll_done "$PARTITIONS"
fi

# ================================================================
# STEP 3: Tear down workers (unless --keep-vms)
# ================================================================
if [[ $KEEP_VMS -eq 0 ]]; then
  echo "deleting worker VMs..."
  for vm in "${VM_NAMES[@]}"; do
    run_cmd gcloud compute instances delete "$vm" \
      --zone="$ZONE" --quiet &
  done
  wait
  echo "all workers deleted"
else
  echo "keeping VMs (--keep-vms): ${VM_NAMES[*]}"
fi

# Clean up done markers
if [[ $DRY_RUN -eq 0 ]]; then
  gsutil -m rm "${BUCKET}/.done_part*" 2>/dev/null || true
fi

# ================================================================
# STEP 4: Merge + postprocess (optional, --merge)
# ================================================================
if [[ $DO_MERGE -eq 1 ]]; then
  if [[ $DRY_RUN -eq 1 ]]; then
    echo "[dry-run] would merge ${PARTITIONS} partitions from ${BUCKET}/ and run svaba_postprocess.sh"
    echo "svaba_cloud.sh: finished"
    exit 0
  fi

  echo "============================================================"
  echo "MERGE: downloading and postprocessing"
  echo "============================================================"

  MERGE_DIR="${ID}_cloud_merge"
  mkdir -p "$MERGE_DIR"
  cd "$MERGE_DIR"

  # Download all bps.txt.gz files
  for i in $(seq 1 "$PARTITIONS"); do
    gsutil cp "${BUCKET}/${ID}_part${i}.bps.txt.gz" .
  done

  # Concatenate (gzip is concat-safe per RFC 1952)
  echo "concatenating ${PARTITIONS} bps.txt.gz files..."
  cat ${ID}_part*.bps.txt.gz > "${ID}.bps.txt.gz"

  # Merge contigs BAMs
  echo "downloading and merging contigs BAMs..."
  CONTIG_BAMS=()
  for i in $(seq 1 "$PARTITIONS"); do
    gsutil cp "${BUCKET}/${ID}_part${i}.contigs.bam" .
    CONTIG_BAMS+=("${ID}_part${i}.contigs.bam")
  done
  samtools merge -f "${ID}.contigs.bam" "${CONTIG_BAMS[@]}"
  rm -f "${CONTIG_BAMS[@]}"

  # Merge runtime.txt (just cat, first file has the header)
  echo "merging runtime.txt..."
  gsutil cp "${BUCKET}/${ID}_part1.runtime.txt" "${ID}.runtime.txt"
  for i in $(seq 2 "$PARTITIONS"); do
    gsutil cp "${BUCKET}/${ID}_part${i}.runtime.txt" tmp_rt.txt
    tail -n +2 tmp_rt.txt >> "${ID}.runtime.txt"
    rm -f tmp_rt.txt
  done

  # Merge optional r2c.txt.gz if present
  R2C_FILES=()
  for i in $(seq 1 "$PARTITIONS"); do
    if gsutil -q stat "${BUCKET}/${ID}_part${i}.r2c.txt.gz" 2>/dev/null; then
      gsutil cp "${BUCKET}/${ID}_part${i}.r2c.txt.gz" .
      R2C_FILES+=("${ID}_part${i}.r2c.txt.gz")
    fi
  done
  if [[ ${#R2C_FILES[@]} -gt 0 ]]; then
    echo "merging ${#R2C_FILES[@]} r2c.txt.gz files..."
    cat "${R2C_FILES[@]}" > "${ID}.r2c.txt.gz"
    rm -f "${R2C_FILES[@]}"
  fi

  # Postprocess: sort + dedup bps, filter r2c, sort + dedup + index
  # contigs BAM. Skip the per-read BAM dedup (discordant/corrected/weird
  # are too large to download for most WGS runs).
  SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
  if [[ -x "${SCRIPT_DIR}/svaba_postprocess.sh" ]]; then
    echo "running svaba_postprocess.sh..."
    "${SCRIPT_DIR}/svaba_postprocess.sh" -t "$THREADS" -m 4G "$ID"
  else
    echo "svaba_postprocess.sh not found at ${SCRIPT_DIR}; skipping postprocess"
    echo "run manually: scripts/svaba_postprocess.sh -t 8 -m 4G ${ID}"
  fi

  echo "============================================================"
  echo "DONE. Merged outputs in: $(pwd)"
  echo "  ${ID}.bps.sorted.dedup.txt.gz"
  echo "  ${ID}.contigs.bam"
  echo "  ${ID}.runtime.txt"
  echo "============================================================"
fi

echo "svaba_cloud.sh: finished"
