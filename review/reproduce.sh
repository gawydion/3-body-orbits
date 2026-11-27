#!/usr/bin/env bash
# Reproducible independent verification suite for REVIEW.md.
# Requires only g++ (C++20) and python3 (standard library).
# Compiled binaries and intermediate JSONL are placed in bin/ (gitignored).
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

mkdir -p bin

echo "=== Compiling independent 80-bit long double verification tools into bin/ ==="
CXX="${CXX:-g++}"
CXXFLAGS="-O3 -std=c++20 -pthread"

$CXX $CXXFLAGS review/audit_orbits.cpp           -o bin/review_audit_orbits
$CXX $CXXFLAGS review/verify_collinear_ultra.cpp -o bin/review_verify_collinear_ultra
$CXX $CXXFLAGS review/test_2000T.cpp             -o bin/review_test_2000T
$CXX $CXXFLAGS review/bhh_vertical_tune.cpp      -o bin/review_bhh_vertical_tune
$CXX $CXXFLAGS review/bhh_planar_floquet.cpp     -o bin/review_bhh_planar_floquet
$CXX $CXXFLAGS review/cont_perp_to_equal.cpp     -o bin/review_cont_perp_to_equal
echo "Compilation succeeded."
echo

python3 review/summarize_audit.py "$@"

echo
echo "=== [4/4] Origin proofs: Planar BHH 5:2 vertical resonance (3D-Delta14) & 2D-DEM-006 continuation ==="
echo "--- Planar BHH vertical tune sweep at m = [1, 1, 0.5] (locating nu_z = 2/5 = 0.40000) ---"
./bin/review_bhh_vertical_tune
echo
echo "--- Planar BHH 1-period and 5-period Floquet indices at vy2 = 0.8156 ---"
./bin/review_bhh_planar_floquet
echo
echo "--- Mass continuation of 2D-DEM-006 (AAAbbb shuttle) to equal masses [1, 1, 1] ---"
./bin/review_cont_perp_to_equal

echo
echo "=== Audit complete. Full 124-orbit JSONL saved to bin/review_audit_results.jsonl ==="
