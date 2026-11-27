#!/usr/bin/env bash
# Reproducible independent verification suite for REVIEW.md (post-PR #2 main audit).
# Requires only g++ (C++20) and python3 (standard library).
# Compiled binaries are placed in bin/ (gitignored).
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

mkdir -p bin

echo "=== Compiling independent 80-bit long double auditors into bin/ ==="
CXX="${CXX:-g++}"
CXXFLAGS="-O3 -std=c++20 -pthread"

$CXX $CXXFLAGS review/audit_4b.cpp  -o bin/review_audit_4b
$CXX $CXXFLAGS review/audit_1pn.cpp -o bin/review_audit_1pn
echo "Compilation succeeded."
echo

echo "=== Running built-in repository self-tests ==="
if [ ! -x bin/search_1pn ] || [ ! -x bin/search_4body ]; then
    make -j4 bin/search_1pn bin/search_4body
fi
./bin/search_1pn --verify-physics
./bin/search_4body --verify-eig
python3 scripts/classify_topology.py --verify
echo

python3 review/audit_brake_and_choreo.py "$@"
echo
echo "=== Post-PR #2 independent audit complete ==="
