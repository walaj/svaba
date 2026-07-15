#!/usr/bin/env bash
#
# build_svaba_image.sh — build a svaba worker image from a specific git commit.
#
# Boots the builder VM, checks out the requested commit/branch/tag, rebuilds
# (submodules + svaba at -O3 -march=native, the cloud-perf build), VERIFIES it
# compiled, installs svaba onto PATH, then snapshots the builder's boot disk to a
# new image. The data disk is untouched (workers attach it read-only at runtime).
#
# GCP image names can't contain underscores, so the default name uses hyphens:
#   svaba-<short-hash>   (override with --image).
#
# Usage:
#   build_svaba_image.sh --git <hash|branch|tag> [options]
#
#   --git REF          commit / branch / tag to build               (required)
#   --image NAME       image name to create        (default: svaba-<short-ref>)
#   --builder NAME     builder VM                   (default: svaba-image-builder)
#   --zone ZONE        GCP zone                     (default: us-central1-a)
#   --project PROJECT  GCP project                  (default: gcloud config)
#   --repo PATH        svaba repo on the builder    (default: /home/jeremiahwala/git/svaba)
#   --family NAME      image family                 (default: svaba)
#   --keep-running     leave the builder VM running afterwards
#   --dry-run          print the commands without running them
#   -h, --help         this message
#
set -euo pipefail

GIT=""; IMAGE=""; BUILDER="svaba-image-builder"; ZONE="us-central1-a"
PROJECT=""; REPO="/home/jeremiahwala/git/svaba"; FAMILY="svaba"
KEEP_RUNNING=0; DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --git)          GIT="$2";       shift 2 ;;
    --image)        IMAGE="$2";     shift 2 ;;
    --builder)      BUILDER="$2";   shift 2 ;;
    --zone)         ZONE="$2";      shift 2 ;;
    --project)      PROJECT="$2";   shift 2 ;;
    --repo)         REPO="$2";      shift 2 ;;
    --family)       FAMILY="$2";    shift 2 ;;
    --keep-running) KEEP_RUNNING=1; shift ;;
    --dry-run)      DRY_RUN=1;      shift ;;
    -h|--help)      sed -n '2,30p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done

[[ -z "$GIT" ]] && { echo "ERROR: --git <hash|branch|tag> is required" >&2; exit 2; }

# default image name from the ref (sanitized to a valid GCP name)
SHORT=$(printf '%s' "$GIT" | cut -c1-12 | tr '[:upper:]_./' '[:lower:]---' | tr -cd 'a-z0-9-')
[[ -z "$IMAGE" ]] && IMAGE="svaba-${SHORT}"
PROJ=(); [[ -n "$PROJECT" ]] && PROJ=(--project "$PROJECT")

run() { if [[ $DRY_RUN -eq 1 ]]; then echo "[dry-run] $*"; else "$@"; fi; }

echo "============================================================"
echo "build_svaba_image.sh"
echo "  git ref   : $GIT"
echo "  -> image  : $IMAGE   (family $FAMILY)"
echo "  builder   : $BUILDER  zone $ZONE  repo $REPO"
echo "============================================================"

# remote build script: checkout the ref, rebuild submodules + svaba at -O3
# -march=native, verify, install onto PATH. \$(nproc) stays literal (runs on VM);
# $REPO / $GIT are filled in here.
REMOTE=$(cat <<EOF
set -euo pipefail
cd "$REPO"
echo "== fetch + checkout $GIT =="
git fetch --all --tags --prune
git checkout --force "$GIT"
git submodule update --init --recursive || true
echo "== rebuild bwa + fermi-lite at -O3 -march=native =="
make -C SeqLib/bwa        clean || true
make -C SeqLib/fermi-lite clean || true
make -C SeqLib/bwa        -j \$(nproc) CFLAGS="-g -O3 -march=native -fno-omit-frame-pointer -Wall -Wno-unused-function"
make -C SeqLib/fermi-lite -j \$(nproc) CFLAGS="-g -O3 -march=native -fno-omit-frame-pointer -Wall -Wno-unused-function"
echo "== cmake build svaba (-O3 -march=native) =="
cmake -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo \
      -DCMAKE_CXX_FLAGS_RELWITHDEBINFO="-O3 -march=native -g -DNDEBUG -fno-omit-frame-pointer" >/dev/null
cmake --build build -j \$(nproc)
echo "== verify =="
test -x build/svaba
./build/svaba --version
sudo install -m 0755 build/svaba /usr/local/bin/svaba    # workers call bare 'svaba'
echo "BUILD_OK commit=\$(git rev-parse --short HEAD)"
EOF
)

echo "1/4  starting builder ${BUILDER} ..."
run gcloud compute instances start "$BUILDER" --zone "$ZONE" "${PROJ[@]+"${PROJ[@]}"}"

echo "2/4  building ${GIT} on ${BUILDER} (verifies it compiles) ..."
run gcloud compute ssh "$BUILDER" --zone "$ZONE" "${PROJ[@]+"${PROJ[@]}"}" --command "$REMOTE"
# if the build failed, set -e aborts here and no image is created.

echo "3/4  stopping builder + creating image ${IMAGE} ..."
run gcloud compute instances stop "$BUILDER" --zone "$ZONE" "${PROJ[@]+"${PROJ[@]}"}"
run gcloud compute images create "$IMAGE" \
  --source-disk "$BUILDER" --source-disk-zone "$ZONE" \
  --family "$FAMILY" "${PROJ[@]+"${PROJ[@]}"}"

if [[ $KEEP_RUNNING -eq 1 ]]; then
  echo "4/4  restarting builder ..."
  run gcloud compute instances start "$BUILDER" --zone "$ZONE" "${PROJ[@]+"${PROJ[@]}"}"
else
  echo "4/4  leaving builder stopped (pass --keep-running to restart)"
fi

echo "============================================================"
echo "DONE: image '${IMAGE}' built from ${GIT}."
echo "Use it in a run with:  svaba_cloud.sh --image ${IMAGE} ..."
echo "============================================================"
