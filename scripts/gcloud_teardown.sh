#!/usr/bin/env bash
#
# gcloud_teardown.sh — kill svaba cloud workers and clean up the bucket.
#
# Usage:
#   gcloud_teardown.sh [options]
#
#   --id STR          analysis ID (default: my_run)
#   --partitions N    number of workers (default: 5)
#   --zone ZONE       GCP zone (default: us-central1-a)
#   --bucket gs://... GCS bucket to clean (default: gs://osteosarc-private/svaba-run)
#   --keep-bucket     don't delete bucket contents
#   -h, --help        this message
#

set -euo pipefail

ID="my_run"
PARTITIONS=5
ZONE="us-central1-a"
BUCKET="gs://osteosarc-private/svaba-run"
CLEAN_BUCKET=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --id)          ID="$2";         shift 2 ;;
    --partitions)  PARTITIONS="$2"; shift 2 ;;
    --zone)        ZONE="$2";       shift 2 ;;
    --bucket)      BUCKET="$2";     shift 2 ;;
    --keep-bucket) CLEAN_BUCKET=0;  shift ;;
    -h|--help)     sed -n '2,14p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *)             echo "unknown option: $1" >&2; exit 2 ;;
  esac
done

GCP_ID=$(echo "$ID" | tr '[:upper:]' '[:lower:]' | tr '_' '-' | sed 's/[^a-z0-9-]//g')

# 1. Delete worker VMs
echo "deleting $PARTITIONS workers (${GCP_ID})..."
for i in $(seq 1 "$PARTITIONS"); do
  vm="svaba-${GCP_ID}-worker-${i}"
  if gcloud compute instances describe "$vm" --zone="$ZONE" &>/dev/null; then
    echo "  deleting $vm"
    gcloud compute instances delete "$vm" --zone="$ZONE" --quiet &
  else
    echo "  $vm not found, skipping"
  fi
done
wait

# 2. Clean up the bucket
if [[ $CLEAN_BUCKET -eq 1 ]]; then
  echo "cleaning bucket ${BUCKET}/..."
  gsutil -m rm "${BUCKET}/**" 2>/dev/null || true
  echo "bucket cleaned"
else
  echo "keeping bucket contents (--keep-bucket)"
fi

echo "teardown complete"
