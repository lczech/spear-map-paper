#!/usr/bin/env bash

# Convenience script to build all independent CMake subprojects in this repo.
# Each subproject keeps its own build directory and can also be built standalone.

set -euo pipefail

REPO_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
JOBS="${JOBS:-$(nproc)}"

PROJECTS=(
    "eval-align-libs"
    "eval-map/tools/generate_fragments"
)

for proj in "${PROJECTS[@]}"; do
    echo "==> Building ${proj}"
    cmake -S "${REPO_DIR}/${proj}" -B "${REPO_DIR}/${proj}/build" -DCMAKE_BUILD_TYPE=RELEASE
    cmake --build "${REPO_DIR}/${proj}/build" -j "${JOBS}"
done
