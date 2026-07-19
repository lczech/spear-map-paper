#!/usr/bin/env bash

# Convenience script to build all independent CMake subprojects in this repo.
# Each subproject keeps its own build directory and can also be built standalone.

set -euo pipefail

REPO_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
JOBS="${JOBS:-$(nproc)}"

# Activating a conda/micromamba clang env does not set CC/CXX or shadow cc/gcc, so CMake's
# default detection would otherwise silently pick up system GCC regardless of the active env.
# Default to clang; override by exporting CC/CXX yourself before calling this script.
CC="${CC:-clang}"
CXX="${CXX:-clang++}"

PROJECTS=(
    "eval-align-libs"
    "eval-map/tools/generate_fragments"
)

for proj in "${PROJECTS[@]}"; do
    echo "==> Building ${proj}"
    # Always reconfigure from scratch: CMake caches the compiler choice on first configure
    # and keeps it on reconfigure even if CC/CXX or the active environment change later, so a
    # stale build/ dir from a previous toolchain would otherwise keep using the old compiler.
    rm -rf "${REPO_DIR}/${proj}/build"
    CC="${CC}" CXX="${CXX}" cmake -S "${REPO_DIR}/${proj}" -B "${REPO_DIR}/${proj}/build" -DCMAKE_BUILD_TYPE=RELEASE
    cmake --build "${REPO_DIR}/${proj}/build" -j "${JOBS}"
done
