#!/usr/bin/env bash
#
# svaba_cloud.sh — scatter svaba across GCP VMs, one partition per VM.
#
# Takes the svaba run command you'd normally run locally and distributes
# it across N worker VMs, each handling a chromosome partition. The
# script manages -k (regions) and -a (analysis ID) per worker — you
# provide everything else exactly as you would on the command line.
#
# Workers boot from a pre-built image (--image) that has svaba, htslib,
# jemalloc, and the reference genome baked in. A shared read-only data
# disk provides the BAM files.
#
# Usage:
#   svaba_cloud.sh [cloud options] -- <svaba run args without -k and -a>
#
# Cloud options (before the --):
#   --id STR              analysis ID prefix (required)
#   --data-disk NAME      GCP persistent disk name for BAMs (required)
#   --bucket gs://...     GCS bucket for outputs (required)
#   --partitions N        number of VMs / partitions     (default: 6)
#   --machine TYPE        GCP machine type               (default: n2-highcpu-16)
#   --zone ZONE           GCP zone                       (default: us-central1-a)
#   --image NAME          GCP image for workers           (default: svaba-worker-image)
#   --mount-point DIR     where data disk mounts          (default: /mnt/data)
#   --boot-disk-size STR  boot disk size                  (default: 50GB)
#   --regions-file FILE   custom partition file (one comma-sep line per partition)
#   --merge               merge + postprocess after all workers finish
#   --keep-vms            don't delete worker VMs on completion
#   --dry-run             print commands without executing
#   -h, --help            this message
#
# Everything after -- is passed verbatim to `svaba run` on each worker,
# with -k and -a injected by the script. Do NOT include -k or -a in
# your svaba args.
#
# Examples:
#   # 5-way scatter, n2-highcpu-16, merge at the end
#   svaba_cloud.sh \
#     --id my_run --data-disk svaba-data \
#     --bucket gs://my-bucket/svaba-run \
#     --partitions 5 --machine n2-highcpu-16 --merge \
#     -- \
#     svaba run \
#       -t /mnt/data/tumor.bam -n /mnt/data/normal.bam \
#       -G /home/user/ref/hg38.fa -p 14 \
#       --blacklist /home/user/tracks/blacklist.bed
#
#   # Dry run
#   svaba_cloud.sh --id test --data-disk d1 \
#     --bucket gs://b/r --dry-run \
#     -- svaba run -t /mnt/data/t.bam -G /ref/hg38.fa -p 8
#

set -euo pipefail

# ----------------------------------------------------------- defaults ---
MACHINE="n2-highcpu-16"
ZONE="us-central1-a"
PARTITIONS=6
MOUNT="/mnt/data"
BOOT_DISK_SIZE="50GB"
IMAGE="svaba-worker-image"
DO_MERGE=0
KEEP_VMS=0
DRY_RUN=0
REGIONS_FILE=""

# Required
DATA_DISK=""
BUCKET=""
ID=""

# ---------------------------------------------------------------- help ---
print_usage() {
  sed -n '2,50p' "$0" | sed 's/^# \{0,1\}//'
}

# -------------------------------------------------------------- parse ---
# Everything before -- is cloud options; everything after is svaba args.
SVABA_ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --id)             ID="$2";             shift 2 ;;
    --data-disk)      DATA_DISK="$2";      shift 2 ;;
    --bucket)         BUCKET="$2";         shift 2 ;;
    --partitions)     PARTITIONS="$2";     shift 2 ;;
    --machine)        MACHINE="$2";        shift 2 ;;
    --zone)           ZONE="$2";           shift 2 ;;
    --image)          IMAGE="$2";          shift 2 ;;
    --mount-point)    MOUNT="$2";          shift 2 ;;
    --boot-disk-size) BOOT_DISK_SIZE="$2"; shift 2 ;;
    --regions-file)   REGIONS_FILE="$2";   shift 2 ;;
    --merge)          DO_MERGE=1;          shift ;;
    --keep-vms)       KEEP_VMS=1;          shift ;;
    --dry-run)        DRY_RUN=1;           shift ;;
    -h|--help)        print_usage;         exit 0 ;;
    --)               shift; SVABA_ARGS=("$@"); break ;;
    -*)
      echo "svaba_cloud.sh: unknown option '$1'" >&2
      echo "run 'svaba_cloud.sh --help' for usage" >&2
      exit 2
      ;;
    *)
      # No -- separator; assume everything from here on is svaba args
      SVABA_ARGS=("$@"); break ;;
  esac
done

# ------------------------------------------------------- validate args ---
missing=()
[[ -z "$ID" ]]        && missing+=("--id")
[[ -z "$DATA_DISK" ]] && missing+=("--data-disk")
[[ -z "$BUCKET" ]]    && missing+=("--bucket")
[[ ${#SVABA_ARGS[@]} -eq 0 ]] && missing+=("svaba run args (after --)")

if [[ ${#missing[@]} -gt 0 ]]; then
  echo "svaba_cloud.sh: missing required: ${missing[*]}" >&2
  exit 2
fi

BUCKET=${BUCKET%/}

# Sanitize ID for GCP resource names
GCP_ID=$(echo "$ID" | tr '[:upper:]' '[:lower:]' | tr '_' '-' | sed 's/[^a-z0-9-]//g')
if [[ "$GCP_ID" != "$ID" ]]; then
  echo "note: sanitized ID for GCP names: '$ID' -> '$GCP_ID'"
fi

# Strip "svaba run" prefix if present (we add it back ourselves)
if [[ "${SVABA_ARGS[0]:-}" == "svaba" ]]; then
  SVABA_ARGS=("${SVABA_ARGS[@]:1}")
fi
if [[ "${SVABA_ARGS[0]:-}" == "run" ]]; then
  SVABA_ARGS=("${SVABA_ARGS[@]:1}")
fi

# Build the svaba args string for embedding in the startup script
SVABA_ARGS_STR=""
for arg in "${SVABA_ARGS[@]}"; do
  # Quote args that contain spaces
  if [[ "$arg" == *" "* ]]; then
    SVABA_ARGS_STR+="\"${arg}\" "
  else
    SVABA_ARGS_STR+="${arg} "
  fi
done

# Extract -p value from svaba args for the merge step (default 8)
THREADS=8
for i in "${!SVABA_ARGS[@]}"; do
  if [[ "${SVABA_ARGS[$i]}" == "-p" ]]; then
    THREADS="${SVABA_ARGS[$((i+1))]}"
    break
  fi
done

# ---------------------------------------- default hg38 6-way partition ---
default_regions() {
  local n=$1
  case $n in
    1) echo "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX" ;;
    2)
      echo "chr1,chr2,chr3,chr4,chr5,chr6,chr7"
      echo "chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX" ;;
    3)
      echo "chr1,chr2,chr3"
      echo "chr4,chr5,chr6,chr7,chr8"
      echo "chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX" ;;
    4)
      echo "chr1,chr2"
      echo "chr3,chr4,chr5"
      echo "chr6,chr7,chr8,chr9,chr10"
      echo "chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX" ;;
    5)
      echo "chr1,chr2"
      echo "chr3,chr4,chr5"
      echo "chr6,chr7,chr8,chr9"
      echo "chr10,chr11,chr12,chr13,chr14"
      echo "chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX" ;;
    6)
      echo "chr1,chr2"
      echo "chr3,chr4,chr5"
      echo "chr6,chr7,chr8,chr9"
      echo "chr10,chr11,chr12,chr13"
      echo "chr14,chr15,chr16,chr17,chr18"
      echo "chr19,chr20,chr21,chr22,chrX" ;;
    *)
      for c in chr{1..22} chrX; do echo "$c"; done | head -n "$n" ;;
  esac
}

declare -a REGIONS
if [[ -n "$REGIONS_FILE" ]]; then
  while IFS= read -r line; do REGIONS+=("$line"); done < "$REGIONS_FILE"
  PARTITIONS=${#REGIONS[@]}
else
  while IFS= read -r line; do REGIONS+=("$line"); done < <(default_regions "$PARTITIONS")
fi

if [[ ${#REGIONS[@]} -ne $PARTITIONS ]]; then
  echo "svaba_cloud.sh: got ${#REGIONS[@]} region lines but expected $PARTITIONS" >&2
  exit 2
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
# STEP 1: Create worker VMs
# ================================================================
echo "============================================================"
echo "svaba_cloud.sh: launching $PARTITIONS workers"
echo "  machine=$MACHINE zone=$ZONE image=$IMAGE"
echo "  data_disk=$DATA_DISK bucket=$BUCKET id=$ID"
echo "  svaba args: ${SVABA_ARGS_STR}"
echo "============================================================"

SCRIPT_TMPDIR=$(mktemp -d)
trap 'rm -rf "$SCRIPT_TMPDIR"' EXIT

VM_NAMES=()
for i in $(seq 1 "$PARTITIONS"); do
  vm="svaba-${GCP_ID}-worker-${i}"
  VM_NAMES+=("$vm")

  SCRIPT_FILE="${SCRIPT_TMPDIR}/startup_${i}.sh"
  cat > "$SCRIPT_FILE" <<STARTUP_EOF
#!/bin/bash
exec > /var/log/svaba_startup.log 2>&1
set -euxo pipefail

echo "=== svaba worker ${i} starting at \$(date) ==="

# --- mount data disk ---
mkdir -p ${MOUNT} /mnt/output

echo "waiting for data disk..."
for attempt in \$(seq 1 30); do
  [[ -b /dev/sdb ]] && break
  echo "  attempt \${attempt}: /dev/sdb not yet available"
  sleep 2
done

if [[ ! -b /dev/sdb ]]; then
  echo "FATAL: /dev/sdb not found after 60s"
  echo "available block devices:"
  lsblk
  exit 1
fi

MOUNT_DEV=/dev/sdb
[[ -b /dev/sdb1 ]] && MOUNT_DEV=/dev/sdb1
echo "mounting \${MOUNT_DEV} -> ${MOUNT}"
mount -o ro,noload \${MOUNT_DEV} ${MOUNT} || mount -o ro \${MOUNT_DEV} ${MOUNT}

echo "data disk mounted, contents:"
ls ${MOUNT}/

# --- run svaba ---
PART_ID="${ID}_part${i}"
cd /mnt/output

echo "=== starting svaba at \$(date) ==="
svaba run \\
  ${SVABA_ARGS_STR} \\
  -k ${REGIONS[$((i-1))]} \\
  -a \${PART_ID} \\
  2>&1 | tee \${PART_ID}.startup.log

echo "=== svaba finished at \$(date) ==="
echo "output files:"
ls -lh /mnt/output/

# --- upload results ---
echo "uploading to ${BUCKET}/"
gsutil -m cp \\
  \${PART_ID}.bps.txt.gz \\
  \${PART_ID}.log \\
  \${PART_ID}.startup.log \\
  \${PART_ID}.contigs.bam \\
  \${PART_ID}.runtime.txt \\
  ${BUCKET}/

# Upload optional outputs if they exist (--dump-reads)
for f in \${PART_ID}.discordant.bam \${PART_ID}.corrected.bam \${PART_ID}.r2c.txt.gz; do
  [[ -f "\$f" ]] && gsutil cp "\$f" ${BUCKET}/ || true
done

# Signal completion
echo "=== uploading done marker at \$(date) ==="
echo "DONE" | gsutil cp - ${BUCKET}/.done_part${i}
echo "=== worker ${i} complete ==="
STARTUP_EOF

  echo "[${i}/${PARTITIONS}] creating $vm  regions=${REGIONS[$((i-1))]}"
  run_cmd gcloud compute instances create "$vm" \
    --zone="$ZONE" \
    --machine-type="$MACHINE" \
    --image="$IMAGE" \
    --disk="name=${DATA_DISK},mode=ro,device-name=svaba-data" \
    --boot-disk-size="$BOOT_DISK_SIZE" \
    --scopes=storage-rw \
    --metadata-from-file=startup-script="$SCRIPT_FILE"

  # Stagger VM creation to avoid disk-attach races
  [[ $i -lt $PARTITIONS ]] && sleep 5
done
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
    echo "[dry-run] would merge ${PARTITIONS} partitions from ${BUCKET}/"
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

  # Postprocess
  SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
  if [[ -x "${SCRIPT_DIR}/svaba_postprocess.sh" ]]; then
    echo "running svaba_postprocess.sh..."
    "${SCRIPT_DIR}/svaba_postprocess.sh" -t "$THREADS" -m 4G "$ID"
  else
    echo "svaba_postprocess.sh not found at ${SCRIPT_DIR}; skipping"
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
