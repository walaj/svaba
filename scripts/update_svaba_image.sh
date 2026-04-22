#!/usr/bin/env bash
#
# update_svaba_image.sh — rebuild the svaba worker image from the builder VM.
#
# Stops the builder, deletes the old image, creates a new one, restarts
# the builder. Takes ~2-3 minutes.
#
# Usage:
#   update_svaba_image.sh [options]
#
#   --builder NAME    builder VM name     (default: svaba-image-builder)
#   --image NAME      image name          (default: svaba-worker-image)
#   --zone ZONE       GCP zone            (default: us-central1-a)
#   -h, --help        this message
#

set -euo pipefail

BUILDER="svaba-image-builder"
IMAGE="svaba-worker-image"
ZONE="us-central1-a"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --builder) BUILDER="$2"; shift 2 ;;
    --image)   IMAGE="$2";   shift 2 ;;
    --zone)    ZONE="$2";    shift 2 ;;
    -h|--help) sed -n '2,14p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *)         echo "unknown option: $1" >&2; exit 2 ;;
  esac
done

echo "=== updating image '${IMAGE}' from '${BUILDER}' ==="

echo "1. stopping ${BUILDER}..."
gcloud compute instances stop "$BUILDER" --zone="$ZONE"

echo "2. deleting old image '${IMAGE}'..."
gcloud compute images delete "$IMAGE" --quiet 2>/dev/null || echo "   (no existing image to delete)"

echo "3. creating new image..."
gcloud compute images create "$IMAGE" \
  --source-disk="$BUILDER" \
  --source-disk-zone="$ZONE"

echo "4. restarting ${BUILDER}..."
gcloud compute instances start "$BUILDER" --zone="$ZONE"

echo "=== done. image '${IMAGE}' is ready ==="
