#!/usr/bin/env bash
set -euo pipefail

# Download and unpack Eigen into external/eigen.
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXTERNAL_DIR="$ROOT_DIR/external"
EIGEN_DIR="$EXTERNAL_DIR/eigen"
EIGEN_VERSION="3.4.0"
ARCHIVE_URL="https://gitlab.com/libeigen/eigen/-/archive/${EIGEN_VERSION}/eigen-${EIGEN_VERSION}.tar.gz"

mkdir -p "$EXTERNAL_DIR"

if [ -d "$EIGEN_DIR/Eigen" ]; then
  echo "Eigen headers already present in $EIGEN_DIR"
  exit 0
fi

tmp_dir="$(mktemp -d)"
cleanup() {
  rm -rf "$tmp_dir"
}
trap cleanup EXIT

archive_path="$tmp_dir/eigen.tar.gz"

echo "Downloading Eigen ${EIGEN_VERSION}..."
curl -L "$ARCHIVE_URL" -o "$archive_path"

echo "Extracting Eigen..."
tar -xzf "$archive_path" -C "$tmp_dir"

extracted_dir="$(find "$tmp_dir" -maxdepth 1 -type d -name "eigen-*" | head -n 1)"
if [ -z "$extracted_dir" ]; then
  echo "Failed to locate extracted Eigen directory." >&2
  exit 1
fi

mv "$extracted_dir" "$EIGEN_DIR"
echo "Eigen installed to $EIGEN_DIR"
