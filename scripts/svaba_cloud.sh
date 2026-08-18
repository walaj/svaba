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
#   --chunk-mb M          fixed-size M-Mb shards across hg38 (1 shard/window;
#                         recreates a "1Mb-per-shard" scatter). Overrides --partitions.
#   --max-concurrent N    cap simultaneously-running VMs; process in waves of N
#                         (default 0 = launch all at once). Use for big shard counts
#                         to respect vCPU quota and avoid overloading the BAM disk.
#   --spot                use Spot/preemptible VMs (~70% cheaper; can be killed —
#                         re-run to fill any shard left without a .done_part marker)
#   --bam-params-gcs gs://...  each worker fetches this insert-size cache and runs
#                         --bam-params (precompute once with: svaba run --learn-only
#                         --bam-params p.tsv ... ; gsutil cp p.tsv gs://...)
#   --merge               merge + postprocess after all workers finish
#   --keep-vms            keep worker VMs running after their shard finishes
#                         (default: each worker self-deletes the moment its
#                         done-marker is written, to stop billing immediately)
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
SPOT=0              # --spot: use preemptible/Spot VMs (~70% cheaper, can be killed)
MAX_CONCURRENT=0    # --max-concurrent N: cap simultaneously-running workers (0=all at once)
CHUNK_MB=0          # --chunk-mb M: fixed-size M-Mb shards across hg38 (instead of chrom groups)
BAM_PARAMS_GCS=""   # --bam-params-gcs gs://...: each worker fetches this and runs --bam-params

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
    --spot)           SPOT=1;              shift ;;
    --max-concurrent) MAX_CONCURRENT="$2"; shift 2 ;;
    --chunk-mb)       CHUNK_MB="$2";       shift 2 ;;
    --bam-params-gcs) BAM_PARAMS_GCS="$2"; shift 2 ;;
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

# --chunk-mb M : fixed-size M-Mb shards across the hg38 primary assembly (one
# shard per window). Recreates the "1Mb-per-shard" style scatter. Uses embedded
# GRCh38 chrom lengths so it works without a local .fai. For non-hg38, use
# --regions-file instead.
hg38_windows() {
  local w=$(( $1 * 1000000 ))
  local names=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
               chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)
  local lens=(248956422 242193529 198295559 190214555 181538259 170805979 159345973 \
              145138636 138394717 133797422 135086622 133275309 114364328 107043718 \
              101991189 90338345 83257441 80373285 58617616 64444167 46709983 50818468 \
              156040895 57227415)
  local i s e
  for i in "${!names[@]}"; do
    s=0
    while [[ $s -lt ${lens[$i]} ]]; do
      e=$(( s + w )); [[ $e -gt ${lens[$i]} ]] && e=${lens[$i]}
      echo "${names[$i]}:$(( s + 1 ))-${e}"
      s=$e
    done
  done
}

declare -a REGIONS
if [[ -n "$REGIONS_FILE" ]]; then
  while IFS= read -r line; do REGIONS+=("$line"); done < "$REGIONS_FILE"
  PARTITIONS=${#REGIONS[@]}
elif [[ "$CHUNK_MB" -gt 0 ]]; then
  while IFS= read -r line; do REGIONS+=("$line"); done < <(hg38_windows "$CHUNK_MB")
  PARTITIONS=${#REGIONS[@]}
  echo "svaba_cloud.sh: --chunk-mb ${CHUNK_MB} -> ${PARTITIONS} fixed-size shards"
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

# --- spot + bam-params plumbing ---
SPOT_ARGS=()
[[ $SPOT -eq 1 ]] && SPOT_ARGS=(--provisioning-model=SPOT --instance-termination-action=DELETE --no-restart-on-failure)
BP_FETCH=""; BP_ARG=""
if [[ -n "$BAM_PARAMS_GCS" ]]; then
  BP_FETCH="echo 'fetching bam-params'; gsutil cp ${BAM_PARAMS_GCS} /tmp/bam_params.tsv"
  BP_ARG="--bam-params /tmp/bam_params.tsv"
fi

# create one worker VM for partition $1
create_worker() {
  local i=$1
  local vm="svaba-${GCP_ID}-worker-${i}"
  local SCRIPT_FILE="${SCRIPT_TMPDIR}/startup_${i}.sh"
  # by default each worker deletes itself once its shard's done-marker is written,
  # so a finished shard stops billing immediately (no waiting for slow siblings,
  # and it works even if this orchestrator isn't running). --keep-vms opts out.
  local SELF_DELETE=""
  if [[ $KEEP_VMS -eq 0 ]]; then
    SELF_DELETE="echo 'shard done -- self-deleting to stop billing'; gcloud compute instances delete \$(hostname) --zone=${ZONE} --quiet || sudo shutdown -h now"
  fi
  cat > "$SCRIPT_FILE" <<STARTUP_EOF
#!/bin/bash
exec > /var/log/svaba_startup.log 2>&1
set -euxo pipefail
echo "=== svaba worker ${i} starting at \$(date) ==="
mkdir -p ${MOUNT} /mnt/output
# Mount the data disk by its STABLE device-name (svaba-data, set in --disk), not
# /dev/sdb: the kernel name varies (sdb/sdc/nvme...) and the filesystem may be on
# a partition. GCP exposes /dev/disk/by-id/google-<device-name>[-part1].
echo "waiting for + mounting data disk (device-name svaba-data)..."
DEVBASE=/dev/disk/by-id/google-svaba-data
mok=0
for a in \$(seq 1 60); do
  for d in "\${DEVBASE}-part1" "\$DEVBASE"; do
    [ -b "\$d" ] || continue
    if mount -o ro "\$d" ${MOUNT} 2>/dev/null || mount -o ro,noload "\$d" ${MOUNT} 2>/dev/null; then mok=1; break; fi
  done
  [ \$mok -eq 1 ] && break
  sleep 3
done
[ \$mok -eq 1 ] || { echo "FATAL: data disk (device-name svaba-data) not mountable"; lsblk; ls -l /dev/disk/by-id/ || true; exit 1; }
echo "mounted: \$(mount | grep ${MOUNT})"
ls ${MOUNT}/
${BP_FETCH}
PART_ID="${ID}_part${i}"
cd /mnt/output
echo "=== starting svaba at \$(date) ==="
stdbuf -oL -eL svaba run \\
  ${SVABA_ARGS_STR} ${BP_ARG} \\
  -k ${REGIONS[$((i-1))]} \\
  -a \${PART_ID} \\
  2>&1 | tee \${PART_ID}.startup.log
echo "=== svaba finished at \$(date) ==="
ls -lh /mnt/output/
gsutil -m cp \\
  \${PART_ID}.bps.txt.gz \${PART_ID}.log \${PART_ID}.startup.log \\
  \${PART_ID}.contigs.bam \${PART_ID}.runtime.txt ${BUCKET}/
for f in \${PART_ID}.discordant.bam \${PART_ID}.corrected.bam \${PART_ID}.r2c.db; do
  [[ -f "\$f" ]] && gsutil cp "\$f" ${BUCKET}/ || true
done
echo "DONE" | gsutil cp - ${BUCKET}/.done_part${i}
echo "=== worker ${i} complete ==="
${SELF_DELETE}
STARTUP_EOF
  echo "  creating $vm  regions=${REGIONS[$((i-1))]}"
  run_cmd gcloud compute instances create "$vm" \
    --zone="$ZONE" --machine-type="$MACHINE" --image="$IMAGE" \
    --disk="name=${DATA_DISK},mode=ro,device-name=svaba-data" \
    --boot-disk-size="$BOOT_DISK_SIZE" --scopes=storage-rw,compute-rw \
    ${SPOT_ARGS[@]+"${SPOT_ARGS[@]}"} \
    --metadata-from-file=startup-script="$SCRIPT_FILE"
}

delete_worker() {
  # tolerant: workers self-delete on completion, so this is a fallback that may
  # find the VM already gone (--keep-vms, a failed shard, or a slow self-delete).
  run_cmd gcloud compute instances delete "svaba-${GCP_ID}-worker-${1}" --zone="$ZONE" --quiet || true
}

# wait until ${BUCKET}/.done_part<i> exists for every i in "$@"
poll_parts() {
  while true; do
    local have; have=$(gsutil ls "${BUCKET}/.done_part*" 2>/dev/null || true)
    local missing=0 i
    for i in "$@"; do echo "$have" | grep -q "/\.done_part${i}\$" || missing=$((missing+1)); done
    [[ $missing -eq 0 ]] && { echo "  all $# done"; return 0; }
    echo "  $(date +%H:%M:%S)  $missing of $# not yet done"
    sleep 60
  done
}

if [[ "$MAX_CONCURRENT" -le 0 ]]; then
  # ---- all-at-once (original behavior) ----
  for i in $(seq 1 "$PARTITIONS"); do create_worker "$i"; [[ $i -lt $PARTITIONS ]] && sleep 5; done
  echo "all $PARTITIONS VMs created"
  if [[ $DRY_RUN -eq 0 ]]; then
    echo "waiting for workers (checking .done_part* every 60s)..."
    poll_parts $(seq 1 "$PARTITIONS")
    if [[ $KEEP_VMS -eq 0 ]]; then
      echo "deleting workers..."; for i in $(seq 1 "$PARTITIONS"); do delete_worker "$i" & done; wait
    fi
  fi
else
  # ---- bounded-concurrency: waves of <= MAX_CONCURRENT VMs ----
  echo "wave mode: $PARTITIONS shards, <= $MAX_CONCURRENT VMs running at once"
  s=1
  while [[ $s -le $PARTITIONS ]]; do
    e=$(( s + MAX_CONCURRENT - 1 )); [[ $e -gt $PARTITIONS ]] && e=$PARTITIONS
    wave=($(seq "$s" "$e"))
    echo "--- wave: workers ${s}..${e} (of $PARTITIONS) ---"
    for i in "${wave[@]}"; do create_worker "$i"; sleep 3; done
    if [[ $DRY_RUN -eq 0 ]]; then
      poll_parts "${wave[@]}"
      if [[ $KEEP_VMS -eq 0 ]]; then for i in "${wave[@]}"; do delete_worker "$i" & done; wait; fi
    fi
    s=$(( e + 1 ))
  done
fi

# Clean up done markers
[[ $DRY_RUN -eq 0 ]] && gsutil -m rm "${BUCKET}/.done_part*" 2>/dev/null || true

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
