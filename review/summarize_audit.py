#!/usr/bin/env python3
"""
Independent audit driver & report generator for 3-body-orbits.

Feeds catalog.json into the compiled 80-bit long double C++20 verification tools:
  - bin/review_audit_orbits           (all 124 orbits: exact variational monodromy,
                                       QR eigenvalues, Procrustes RPO/cover scan,
                                       adaptive S3 free-group words)
  - bin/review_verify_collinear_ultra (17-digit half-period symplectic involution
                                       refinement + full-period fallback on flagged
                                       collinear_L0 orbits)
  - bin/review_test_2000T             (2,000-period empirical perturbation test)

Usage:
  python3 review/summarize_audit.py [--workers N] [--out bin/review_audit_results.jsonl]
"""
import argparse
import concurrent.futures
import importlib.util
import json
import os
import subprocess
import sys
from collections import Counter


def to_3d_state(o):
    y = o["y0"]
    if o["dim"] == 3:
        return y
    return [
        y[0], y[1], 0.0, y[2], y[3], 0.0, y[4], y[5], 0.0,
        y[6], y[7], 0.0, y[8], y[9], 0.0, y[10], y[11], 0.0,
    ]


def format_audit_line(o):
    m = o["masses"]
    y3 = to_3d_state(o)
    precess = o.get("precess_rad", 0.0)
    vals = [
        o["id"],
        str(o["dim"]),
        o["subspace"],
        f"{m[0]:.17g}", f"{m[1]:.17g}", f"{m[2]:.17g}",
        f"{o['T']:.17g}", f"{precess:.17g}",
    ] + [f"{v:.17g}" for v in y3]
    return " ".join(vals) + "\n"


def run_audit_chunk(bin_path, lines):
    proc = subprocess.run(
        [bin_path],
        input="".join(lines),
        text=True,
        capture_output=True,
        check=True,
    )
    return [json.loads(ln) for ln in proc.stdout.splitlines() if ln.strip()]


def extract_nontrivial_k_3d(eigs_raw):
    """Extract the 3 nontrivial reciprocal Broucke indices k = lam + 1/lam from 18x18 inertial eigs."""
    eigs = [complex(r, i) for r, i in eigs_raw]
    rem = list(eigs)
    for _ in range(12):
        idx = min(range(len(rem)), key=lambda i: abs(rem[i] - 1.0))
        rem.pop(idx)
    ks = []
    used = [False] * len(rem)
    for i in range(len(rem)):
        if used[i]:
            continue
        best_j = min(
            (j for j in range(i + 1, len(rem)) if not used[j]),
            key=lambda j: abs(rem[i] * rem[j] - 1.0),
        )
        used[i] = used[best_j] = True
        ks.append((rem[i] + rem[best_j]).real)
    return sorted(ks, reverse=True)


def run_repo_classifier(orbits):
    """Load scripts/classify_topology.py and compute its reported canonical words for all 2D-DEM orbits."""
    path = os.path.join("scripts", "classify_topology.py")
    if not os.path.exists(path):
        return {}
    spec = importlib.util.spec_from_file_location("repo_classify_topology", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    out = {}
    for o in orbits:
        if o["dim"] == 2 and o.get("democratic"):
            _, canon = mod.compute_word(o)
            out[o["id"]] = canon
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", default="catalog.json")
    ap.add_argument("--workers", type=int, default=min(os.cpu_count() or 4, 8))
    ap.add_argument("--out", default="bin/review_audit_results.jsonl")
    args = ap.parse_args()

    with open(args.catalog) as f:
        data = json.load(f)
    orbits = data["orbits"] if isinstance(data, dict) else data
    cat_by_id = {o["id"]: o for o in orbits}

    print(f"=== [1/4] Auditing all {len(orbits)} orbits with 80-bit analytical variational QR ({args.workers} workers) ===")
    chunks = [[] for _ in range(args.workers)]
    for idx, o in enumerate(orbits):
        chunks[idx % args.workers].append(format_audit_line(o))
    chunks = [c for c in chunks if c]

    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=len(chunks)) as ex:
        futs = [ex.submit(run_audit_chunk, "bin/review_audit_orbits", ch) for ch in chunks]
        for fu in futs:
            results.extend(fu.result())

    order = {o["id"]: i for i, o in enumerate(orbits)}
    results.sort(key=lambda r: order[r["id"]])
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as f:
        for r in results:
            f.write(json.dumps(r) + "\n")
    print(f"Wrote {len(results)} audit records to {args.out}")

    # 1. 3D Audit Summary
    print("\n--- 3D Orbits: Exact Monodromy, Closure, and Swaps ---")
    for r in results:
        if r["dim"] != 3:
            continue
        o = cat_by_id[r["id"]]
        if r["subspace"] == "3d_rpo_rotating":
            ks_str = f"cat=({o['k1']:+.4f}, {o['k2']:+.4f}, {o['k3']:+.4f}) max|lam|={r['max_eig_3d']:.6f}"
        else:
            ks = extract_nontrivial_k_3d(r["eigs_3d"])
            margin = min(2.0 - abs(k) for k in ks)
            ks_str = (
                f"exact_k=({ks[0]:+.4f}, {ks[1]:+.4f}, {ks[2]:+.4f}) margin={margin:.4f} "
                f"[cat=({o['k1']:+.4f}, {o['k2']:+.4f}, {o['k3']:+.4f}) margin={o['stab_margin']:.4f}]"
            )
        print(
            f"  {r['id']} ({r['subspace']:18s}): ret={r['ret_err_ld']:.2e} symp={r['symp_err']:.2e} "
            f"swaps={r['exact_swaps']:3d} (cat={o['partner_swaps']:3d}) E={r['E']:.5f} | {ks_str}"
        )

    # 2. 2D perp_L_nonzero Audit Summary
    perp = [r for r in results if r["subspace"] == "perp_L_nonzero"]
    perp_hie = [r for r in perp if r["id"].startswith("2D-HIE")]
    perp_dem = [r for r in perp if r["id"].startswith("2D-DEM")]
    hie_rpo = sum(1 for r in perp_hie if r["q_rpo"] > 1)
    hie_covers = [r["id"] for r in perp_hie if r["q_inert"] > 1]
    dem_rpo = [r for r in perp_dem if r["q_rpo"] > 1]
    dem_prim = [r for r in perp_dem if r["q_rpo"] == 1]

    print(f"\n--- 2D perp_L_nonzero ({len(perp)} total = {len(perp_hie)} HIE + {len(perp_dem)} DEM) ---")
    print(f"  2D-HIE q_rpo > 1 (BHH closures): {hie_rpo}/{len(perp_hie)}")
    print(f"  2D-HIE exact 2x inertial covers (q_inert > 1): {len(hie_covers)} -> {hie_covers}")
    print(f"  2D-DEM q_rpo > 1 (BHH closures): {len(dem_rpo)}/{len(perp_dem)} "
          f"(S3 words: {dict(Counter(r['word_can_s3'] for r in dem_rpo))})")
    print(f"  2D-DEM q_rpo == 1 (primitive): {len(dem_prim)}/{len(perp_dem)}:")
    for r in dem_prim:
        print(f"    {r['id']} masses={r['masses']} word_can={r['word_can']} word_s3={r['word_can_s3']} "
              f"max|lam|={r['max_eig_3d']:.6f} rmin={min(r['rmin_pair']):.4f}")

    # 3. Free-group classifier comparison against scripts/classify_topology.py
    repo_words = run_repo_classifier(orbits)
    if repo_words:
        wrong_total = wrong_col = total_checked = col_checked = 0
        for r in results:
            if r["id"] in repo_words:
                total_checked += 1
                rw = repo_words[r["id"]] or "1"
                if r["subspace"] == "collinear_L0":
                    col_checked += 1
                if rw != r["word_can"]:
                    wrong_total += 1
                    if r["subspace"] == "collinear_L0":
                        wrong_col += 1
        print(f"\n--- Topology Classifier Audit vs scripts/classify_topology.py ---")
        print(f"  Corrupted words (same basis): {wrong_total}/{total_checked} "
              f"({100.0 * wrong_total / max(1, total_checked):.1f}%), "
              f"of which collinear_L0: {wrong_col}/{col_checked}")
        for showcase in ["2D-DEM-006", "2D-DEM-013", "2D-DEM-016", "2D-DEM-017", "2D-DEM-018", "2D-DEM-022", "2D-DEM-048"]:
            r = next(x for x in results if x["id"] == showcase)
            print(f"    {showcase}: repo='{repo_words.get(showcase, '?')}' -> exact='{r['word_can']}' (S3='{r['word_can_s3']}')")

    # 4. Collinear L=0 ultra-precision half-period symplectic verification & 2,000T test
    col = [r for r in results if r["subspace"] == "collinear_L0"]
    covers_col = [r["id"] for r in col if r["q_inert"] > 1]
    print(f"\n--- 2D collinear_L0 ({len(col)} total) ---")
    print(f"  Exact inertial multiple covers (q_inert > 1): {covers_col}")

    # Run ultra-precision half-period symplectic refinement on the 13 flagged/borderline
    # collinear_L0 orbits plus Figure-8 (2D-DEM-067) as an exact 16-digit benchmark
    target_ids = [
        "2D-DEM-021", "2D-DEM-030", "2D-DEM-033", "2D-DEM-034", "2D-DEM-035",
        "2D-DEM-036", "2D-DEM-037", "2D-DEM-039", "2D-DEM-040", "2D-DEM-042",
        "2D-DEM-049", "2D-DEM-057", "2D-DEM-067", "2D-DEM-069",
    ]
    flagged_orbits = [cat_by_id[oid] for oid in target_ids if oid in cat_by_id]

    print(f"\n=== [2/4] Half-period symplectic monodromy refinement (tol=1e-17L) on {len(flagged_orbits)} collinear_L0 orbits ===")
    ultra_in = []
    for o in flagged_orbits:
        y = o["y0"]
        ultra_in.append(f"{o['id']} {o['masses'][2]:.17g} {y[6]:.17g} {y[7]:.17g} {o['T']:.17g} {o['k1']:.17g} {o['k2']:.17g}\n")
    proc_u = subprocess.run(
        ["bin/review_verify_collinear_ultra"],
        input="".join(ultra_in),
        text=True,
        capture_output=True,
        check=True,
    )
    unstable_ids = []
    parabolic_ids = []
    sub_ids = []
    for ln in proc_u.stdout.splitlines():
        print(" ", ln)
        if "VERDICT=UNSTABLE" in ln:
            unstable_ids.append(ln.split()[0])
        elif "VERDICT=NEAR-PARABOLIC" in ln:
            parabolic_ids.append(ln.split()[0])
        elif "VERDICT=SUB<1e-4" in ln:
            sub_ids.append(ln.split()[0])

    print(f"\nConfirmed hyperbolically unstable (k1 > 2):  {len(unstable_ids)} -> {unstable_ids}")
    print(f"Ultra-near-parabolic (|k1-2| <= 1e-6):       {len(parabolic_ids)} -> {parabolic_ids}")
    print(f"Sub-threshold elliptic (margin < 1e-4):      {len(sub_ids)} -> {sub_ids}")

    print(f"\n=== [3/4] 2,000-period empirical perturbation test on top unstable collinear_L0 orbits ===")
    kick_targets = [oid for oid in ["2D-DEM-037", "2D-DEM-057", "2D-DEM-069", "2D-DEM-049"] if oid in cat_by_id]
    kick_in = []
    for oid in kick_targets:
        o = cat_by_id[oid]
        y = o["y0"]
        kick_in.append(f"{oid} {o['masses'][2]:.17g} {y[6]:.17g} {y[7]:.17g} {o['T']:.17g}\n")
    proc_k = subprocess.run(
        ["bin/review_test_2000T"],
        input="".join(kick_in),
        text=True,
        capture_output=True,
        check=True,
    )
    for ln in proc_k.stdout.splitlines():
        print(" ", ln)


if __name__ == "__main__":
    main()
