# Submission Cover Note — `threebodyorbits.com` (`info@threebodyorbits.com`)

**Proposed Family Name:** `AAAbbb` Shuttle-Exchange (Mass-Stabilized Non-Hierarchical $L \neq 0$ Exchange Family)  
**Primary Showcase Orbit:** `AAAbbb-DEM002` (`2D-DEM-002`, masses $m_1 : m_2 : m_3 = 1 : 0.9 : 0.8$, primitive $q_{\text{RPO}} = 1$, exact inertial period $T \approx 15.49514153$)  
**Additional Verified Benchmarks:** `AAAbbb-PEAK` ($[1, 0.94, 0.71]$, 3D margin `1.5963`), `AAAbbb-DEEP` ($[1, 0.94, 0.70]$, 3D margin `1.4762`), `AAAbbb-ISO89` ($[1, 0.89, 0.89]$, exact half-period exchange symmetry, 3D margin `1.0909`), `AAAbbb-WIDE` ($[1, 0.70, 0.70]$, $r_{\min} = 0.2429$), `AAAbbb-DEM019` ($[1, 0.8, 0.6]$), and `AAAbbb-EQ-RPO13` ($[1, 1, 1]$ full-3D-stable 13-fold inertial rosette closure `(AAAbbb)^13`).

---

## Why This Family Is Distinct and Notable

1. **Absent from Published Archives and the Three-Body Orbits Atlas:**
   - Verified absent from the live `threebodyorbits.com` Atlas (3,942 orbits in 41 families) and `/hunt/garden` (82,032 loops; shortest loop there has 16 syzygies/letters).
   - Verified absent from Li & Liao (2017, 695 equal-mass $L=0$), Li, Jing & Liao (2018, 1,349 unequal-mass $L=0$), Li & Liao (2020, 598,996 BHH $A^k$), Li, Li & Liao (2021, 135,445 non-hierarchical unequal-mass Moth I `ABAbaBABab` orbits), Janković & Dmitrašinović (2016), Šuvakov & Dmitrašinović (2013, 2014), Sheen (2016, all 17 orbits audited), Broucke (1975 A/R), Broucke & Boggs (1975 S-series, including `S4` and `S122`, both audited as `AAA`), Simó (2001), and Hristov et al. (2025).

2. **Physical Geometry (`AAAbbb` 3 + 3 Co-Rotating Shuttle Exchange):**
   - Planar Newtonian 3-body orbit ($G = 1$, nonzero angular momentum $L \neq 0$, scale-invariant $\Lambda = |L|\sqrt{|E|} \in [1.28, 2.14]$), with collinear perpendicular crossings at $t = 0$ and $t = T/2$.
   - Inertial pairwise windings are $(w_{12}, w_{13}, w_{23}) = (-5, -5, -2)$: while the lighter companions $(m_2, m_3)$ complete 2 revolutions around each other, the heaviest body $m_1$ shuttles back and forth between them, completing **3 prograde co-rotating loops around $m_2$ (`AAA`) and 3 opposite-chirality loops around $m_3$ (`bbb`)**.
   - At equal companion masses ($m_2 = m_3$), the orbit possesses exact half-period exchange + $180^\circ$ rotation symmetry and passes through an exact isosceles triangle ($r_{12} = r_{13}$) at $t = T/4$ and $3T/4$.

3. **Remarkable Mass-Stabilization Origin:**
   - Continuing $m = [\varepsilon, 1, 1]$ down to $\varepsilon \to 0$ connects smoothly to an exact 2:1 resonant periodic orbit of the equal-mass **Elliptic Restricted Three-Body Problem** ($e = 0.092555, T = 2 T_{\text{Kep}}$) and from there ($e \to 0$) to a violently unstable doubly-symmetric 3+3 prograde exchange orbit of the **Copenhagen CR3BP** ($C_J = 3.696245$, Hénon stability index $a = +168.90$, Floquet multiplier $|\lambda| \approx 337.8$).
   - Even at equal masses $[1, 1, 1]$ ($\alpha = 0$), the inertial orbit remains hyperbolically unstable ($k_1 = +3.5689$, $|\lambda| = 3.268$).
   - Full-3D linear stability ($|k_1|, |k_2|, |k_z| < 2$) emerges **exclusively** either (a) when the mean companion mass falls below the shuttle mass ($(m_2 + m_3)/2 < m_1$), opening a broad 2D mass-space stability island peaking at 3D margin **`1.5963`** ($[1, 0.94, 0.71]$; with narrow peninsulas reaching $[1, 1, 0.88]$ and $[1, 1.06, 0.88]$) before the isosceles branch folds back at $m_2 = m_3 \approx 0.45746$, or (b) on the equal-mass $[1, 1, 1]$ relative-periodic branch at $\alpha/\pi \in [0.061, 0.096]$, yielding full-3D-stable equal-mass inertial rosettes `(AAAbbb)^q` for $q \in \{11, 12, 13, 14, 15, 16\}$.

4. **34-Digit Independent Verification & 10,000-Period Full-3D Survival:**
   - All initial conditions in `aaabbb_submission_34digits.txt` and `aaabbb_quad32.json` were converged in 128-bit quad precision (`__float128`) using an independent 8-stage 16th-order exact symplectic Gauss–Legendre integrator (`GL16`, self-tested to $< 4 \times 10^{-34}$ order-16 / symplecticity error) with step-halving cross-check (`h` vs `h/2`), achieving return errors $< 10^{-30}$ and canonical variational symplecticity errors $< 10^{-28}$.
   - All stable representatives were empirically verified over **10,000 full periods** in full 3D (`GL8` symplectic) under zero-momentum 3D velocity perturbations $\delta = 10^{-5}$ and $10^{-4}$, remaining tightly bound with zero escape, while the unstable equal-mass control $[1, 1, 1]$ disintegrates rapidly under the same $\delta = 10^{-5}$ kick.

**Attached Machine-Readable Files:**
- `data/aaabbb/aaabbb_submission_34digits.txt` — human-readable 34-digit initial conditions table
- `data/aaabbb/aaabbb_quad32.json` — full 34-digit JSON with stability indices, windings, and verification bounds
- `data/aaabbb/aaabbb_longterm_3d.json` — 10,000-period full-3D perturbation survival log
- `src/confirm_aaabbb.cpp` — self-contained C++20 (`__float128`) verification source code
