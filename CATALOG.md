# 3-Body & 4-Body Steady (Linearly Stable) Periodic Orbit Discovery Catalog
Generated: `2026-11-27T16:33:53.174127`  
**Total unique 80-bit-audited primitive Newtonian 3-body steady orbits:** **133** (from 164 raw detections; 22 covers/unstable/near-parabolic/sub-threshold candidates strictly filtered)  
**Additional 80-bit-audited extension catalogs:** **21 Post-Newtonian (1PN ADM) relativistic steady orbits** · **3 Newtonian 4-body (full-3D-stable 24×24 variational) steady orbits** (8 weakly hyperbolic 3D 4-body candidates rejected)

- **3D Non-Planar Newtonian 3-Body Steady Orbits:** 16 (4 non-hierarchical 3D helical choreographies, 10 BHH 3D vertical subharmonic orbits [`3D-Δ14` & `3D-PD40/48`], 2 inclined 3:1 prograde hierarchical triple)
- **2D Non-Hierarchical (Exchange) Steady Orbits (`2D-DEM-*`):** 55 (33 collinear $L=0$ braids, 14 zero-velocity free-fall brake orbits, 6 2D rotating choreographies, 2 non-hierarchical `AAAbbb` shuttle-exchange orbits)
- **2D Hierarchical / BHH / Co-Orbital Steady Orbits (`2D-HIE-*`):** 62 (58 primitive $q$-fold BHH retrograde inertial closures, 2 BHH satellites, 2 1:1 eccentric planetary co-orbitals)

## Verification Criteria (80-Bit Analytical Variational Standard)
Every Newtonian 3-body orbit in this catalog passed **all five** independent audit criteria in `./bin/audit_catalog` (`long double`, 64-bit significand):
1. **Periodicity to high precision:** $|\mathbf{z}(T) - \mathbf{z}(0)| < 10^{-7}$ (most $10^{-11}\text{--}10^{-9}$) using 80-bit adaptive Dormand–Prince 5(4).
2. **Exact analytical variational monodromy & canonical symplecticity:** Variational equations $\dot{M} = J_f(\mathbf{y}) M$ integrated analytically (replacing finite differences), enforcing canonical symplecticity error $\|M_{\text{can}}^\top J M_{\text{can}} - J\|_\infty < 10^{-6}$ (and half-period symplectic involution $P\mathbf{y}(T/2) = \mathbf{y}(0)$ at $\text{tol}=10^{-17}$ for all `collinear_L0` orbits).
3. **Strict Floquet linear stability margin:** All non-trivial stability indices satisfy $|k_j| < 2 - 10^{-4}$ (confirmed by both Newton trace invariants and Hessenberg–Francis QR eigenvalues).
4. **Strict primitive inertial period (`q_inert = 1`):** Verified by Procrustes / divisor scan over $q = 2\ldots 36$ (all exact multiple covers removed).
5. **100-period empirical perturbation survival & true 3D thickness:** Survives 100/100 full periods after a momentum-preserving $10^{-5}$ kick; 3D orbits additionally satisfy $h_{\text{3D}} > 0.05$.

> **Note on planar $L = 0$ vertical stability:** For any planar $L = 0$ orbit (`collinear_L0` and `free_fall_brake`), infinitesimal spatial tilts preserve $\mathbf{L} = 0$, making all six eigenvalues of the vertical monodromy block $M_z$ identically $+1$ ($k_z = \operatorname{tr}(M_z) - 4 \equiv +2$). This is exact symplectic degeneracy (`L=0 degen`), **not** exponential vertical instability.

## 3D Non-Planar Newtonian 3-Body Steady Orbits
All stability indices $(k_1, k_2, k_3)$ and margins are exact 80-bit analytical variational values at the orbit's stated period $T$ (using the rotating-frame trivial trace $6 + 6\cos(m\theta)$ for RPOs and full $T_{\text{inert}}$ for inertial closures).

| ID | Subspace | Class | Masses | Period $T$ | $|\mathbf{L}|$ | $r_{\min}$ | $h_{\text{3D}}$ | Swaps | Stability $(k_1, k_2, k_3)$ | Margin | Symp Err |
| :--- | :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `3D-001` | `3d_choreo_helical` | `choreography_3d_helical` | `[1, 1, 1]` | 15.5432 | 0.564 | 2.0000 | 0.2101 | 3 | (+1.752, +1.003, -1.084) | **0.2477** | `1.6e-12` |
| `3D-002` | `3d_inertial_5fold` | `bhh_3d_subharmonic` | `[1, 1, 0.5]` | 91.4642 | 0.854 | 0.8438 | 0.2562 | 70 | (+1.800, +1.715, -0.237) | **0.2002** | `3.4e-11` |
| `3D-003` | `3d_vertical_pd_7fold` | `bhh_3d_subharmonic` | `[1, 1, 0.415]` | 335.2194 | 0.831 | 0.7463 | 0.1873 | 280 | (+1.545, -0.283, -1.897) | **0.1029** | `1.6e-09` |
| `3D-004` | `3d_choreo_helical_24fold` | `choreography_3d_helical` | `[1, 1, 1]` | 373.1700 | 0.530 | 2.0000 | 0.6979 | 72 | (+0.833, +0.424, -1.902) | **0.0983** | `5.3e-11` |
| `3D-005` | `3d_choreo_helical_11fold` | `choreography_3d_helical` | `[1, 1, 1]` | 170.9493 | 0.577 | 2.0000 | 0.6969 | 33 | (+1.905, -0.034, -1.661) | **0.0948** | `4.8e-11` |
| `3D-006` | `3d_inertial_12fold` | `bhh_3d_subharmonic` | `[1, 1, 0.5]` | 219.6972 | 0.833 | 0.8468 | 0.1833 | 168 | (+1.567, -1.187, -1.917) | **0.0827** | `4.3e-10` |
| `3D-007` | `3d_rpo_vertical_pd` | `bhh_3d_subharmonic` | `[1, 1, 0.415]` | 47.8586 | 0.836 | 0.7462 | 0.2111 | 40 | (+1.953, +1.953, +0.990) | **0.0466** | `1.4e-11` |
| `3D-008` | `3d_choreo_helical_15fold` | `choreography_3d_helical` | `[1, 1, 1]` | 233.4570 | 0.427 | 2.0000 | 0.6998 | 45 | (+1.313, +0.829, -1.980) | **0.0199** | `3.4e-11` |
| `3D-009` | `3d_vertical_pd_11fold` | `bhh_3d_subharmonic` | `[1, 1, 0.415]` | 526.5652 | 0.834 | 0.7460 | 0.1994 | 440 | (+1.777, -0.091, -1.981) | **0.0190** | `5.3e-09` |
| `3D-010` | `3d_vertical_pd_8fold` | `bhh_3d_subharmonic` | `[1, 1, 0.415]` | 385.4703 | 0.821 | 0.7534 | 0.1109 | 384 | (+1.983, +0.362, -1.754) | **0.0167** | `4.8e-09` |
| `3D-011` | `3d_inertial_7fold` | `bhh_3d_subharmonic` | `[1, 1, 0.5]` | 128.2675 | 0.822 | 0.8524 | 0.1292 | 70 | (+1.984, +0.816, -0.386) | **0.0161** | `7.3e-11` |
| `3D-012` | `3d_vertical_pd_3fold` | `bhh_3d_subharmonic` | `[1, 1, 0.415]` | 144.0555 | 0.825 | 0.7492 | 0.1459 | 120 | (+1.991, -1.524, -1.948) | **0.0091** | `2.8e-10` |
| `3D-013` | `3d_rpo_vertical_pd` | `bhh_3d_subharmonic` | `[1, 1, 0.369]` | 49.9526 | 0.842 | 0.7215 | 0.2007 | 40 | (+1.991, +0.539, +0.310) | **0.0090** | `1.8e-11` |
| `3D-014` | `3d_rpo_rotating` | `bhh_3d_subharmonic` | `[1, 1, 0.5]` | 18.2929 | 0.854 | 0.8438 | 0.1993 | 14 | (+1.991, -0.035, -1.490) | **0.0085** | `5.6e-13` |
| `3D-015` | `3d_distinct_15fold` | `hierarchical_3d_triple` | `[1, 0.7, 0.4]` | 163.6378 | 1.241 | 0.5539 | 0.0922 | 0 | (+1.993, +1.648, -0.993) | **0.0068** | `6.7e-10` |
| `3D-016` | `3d_rpo_rotating` | `hierarchical_3d_triple` | `[1, 0.7, 0.4]` | 11.0328 | 1.239 | 0.5574 | 0.0955 | 0 | (+1.998, +1.903, -0.705) | **0.0024** | `3.4e-12` |

## 2D Non-Hierarchical (Exchange) Steady Orbits
Strictly non-hierarchical 2D exchange orbits (`collinear_L0` braids, `free_fall_brake`, `2d_choreo_*`, and `AAAbbb` shuttle-exchange), ranked by exact 80-bit analytical variational margin $2 - \max(|k_1|, |k_2|)$ (compare Figure-8 margin $= 0.0028$):

| ID | Subspace | Masses | $S_3$ Word | Period $T$ | $L$ | $r_{\min}$ | Swaps | Stability $(k_1, k_2)$ | Margin | 3D Vertical ($k_z$) | Symp Err |
| :--- | :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `2D-DEM-001` | `2d_choreo_5petal_5fold` | `[1, 1, 1]` | `(ABab)^5` | 68.4160 | -0.214 | 0.51853 | 30 | (+0.4103, -0.0012) | **1.5897** | YES (+0.54) | `2.4e-10` |
| `2D-DEM-002` | `perp_L_nonzero` | `[1, 0.9, 0.8]` | `AAAbbb` | 15.4951 | -1.610 | 0.18932 | 2 | (+0.7940, +0.2189) | **1.2060** | YES (-0.26) | `5.8e-11` |
| `2D-DEM-003` | `collinear_L0` | `[1, 1, 1.75]` | `AAbaBAbaBBab` | 13.7455 | +0.000 | 0.00856 | 14 | (+0.4410, -1.1200) | **0.8800** | L=0 degen (+2) | `1.0e-09` |
| `2D-DEM-004` | `2d_choreo_8swap_3fold` | `[1, 1, 1]` | `1` | 31.3931 | -1.354 | 0.09491 | 27 | (+1.1786, -0.4641) | **0.8214** | YES (-0.69) | `1.3e-08` |
| `2D-DEM-005` | `2d_choreo_fig8_2fold` | `[1, 1, 1]` | `(ABab)^2` | 22.6351 | -0.265 | 0.33410 | 12 | (-0.5991, -1.2973) | **0.7027** | YES (+1.60) | `5.6e-11` |
| `2D-DEM-006` | `2d_choreo_fig8_4fold` | `[1, 1, 1]` | `(ABab)^4` | 73.2001 | -0.141 | 0.89915 | 24 | (-1.0300, -1.4056) | **0.5944** | YES (+1.58) | `1.1e-10` |
| `2D-DEM-007` | `2d_choreo_fig8_12fold` | `[1, 1, 1]` | `(ABab)^12` | 159.0740 | -0.222 | 0.48449 | 72 | (+1.4081, -0.2821) | **0.5919** | YES (-1.99) | `1.5e-09` |
| `2D-DEM-008` | `collinear_L0` | `[1, 1, 0.55]` | `AABaaBABaabAAbab` | 30.1457 | +0.000 | 0.07035 | 18 | (-0.9407, -1.5224) | **0.4776** | L=0 degen (+2) | `2.9e-13` |
| `2D-DEM-009` | `collinear_L0` | `[1, 1, 0.6]` | `ABABaBABabABaba…` | 24.5559 | +0.000 | 0.00429 | 32 | (+1.5482, +0.0859) | **0.4518** | L=0 degen (+2) | `7.6e-10` |
| `2D-DEM-010` | `collinear_L0` | `[1, 1, 0.75]` | `ABABabABAbabABa…` | 20.4871 | +0.000 | 0.00166 | 28 | (+0.6434, -1.5734) | **0.4266** | L=0 degen (+2) | `2.8e-10` |
| `2D-DEM-011` | `collinear_L0` | `[1, 1, 1.1]` | `ABAbaBAbaBABabABab` | 26.1606 | +0.000 | 0.03138 | 26 | (-0.6241, -1.6646) | **0.3354** | L=0 degen (+2) | `1.2e-12` |
| `2D-DEM-012` | `collinear_L0` | `[1, 1, 1]` | `ABAbaBAbaBABabABab` | 28.6693 | +0.000 | 0.04970 | 22 | (+1.6978, +1.5534) | **0.3022** | L=0 degen (+2) | `2.4e-13` |
| `2D-DEM-013` | `collinear_L0` | `[1, 1, 0.4]` | `AAbaBABaaBAbabA…` | 28.8540 | +0.000 | 0.00257 | 26 | (+1.7143, -0.0980) | **0.2857** | L=0 degen (+2) | `8.4e-10` |
| `2D-DEM-014` | `free_fall_brake` | `[1, 1, 0.167778568146]` | `1 (3-2-1-2-1-3)` | 3.0095 | +0.000 | 0.00076 | 4 | (-1.6697, -1.7684) | **0.2316** | L=0 degen (+2) | `1.6e-08` |
| `2D-DEM-015` | `free_fall_brake` | `[1, 1, 0.187778568146]` | `1 (3-2-1-2-1-3)` | 2.9483 | +0.000 | 0.00445 | 4 | (+1.2987, -1.7949) | **0.2051** | L=0 degen (+2) | `1.5e-09` |
| `2D-DEM-016` | `collinear_L0` | `[1, 1, 3]` | `ABABABababab` | 16.8941 | +0.000 | 0.02955 | 12 | (+0.2484, -1.8087) | **0.1913** | L=0 degen (+2) | `3.3e-12` |
| `2D-DEM-017` | `collinear_L0` | `[1, 1, 1.05]` | `ABAbaBAbaBABabABab` | 27.9295 | +0.000 | 0.04477 | 26 | (+1.8159, +0.5650) | **0.1841** | L=0 degen (+2) | `2.8e-13` |
| `2D-DEM-018` | `2d_choreo_trefoil_3fold` | `[1, 1, 1]` | `(ABab)^3` | 46.6426 | -0.181 | 0.66990 | 18 | (+1.8348, -1.3390) | **0.1652** | YES (+1.59) | `7.6e-11` |
| `2D-DEM-019` | `perp_L_nonzero` | `[1, 0.8, 0.6]` | `AAAbbb` | 19.0549 | -1.375 | 0.19703 | 2 | (+1.8402, -1.3431) | **0.1598** | YES (+0.12) | `3.6e-11` |
| `2D-DEM-020` | `free_fall_brake` | `[1, 1, 0.162778568146]` | `1 (3-2-1-2-1-3)` | 3.0177 | +0.000 | 0.00053 | 4 | (-1.4535, -1.8554) | **0.1446** | L=0 degen (+2) | `2.5e-08` |
| `2D-DEM-021` | `free_fall_brake` | `[1, 1, 0.173778568146]` | `1 (3-2-1-2-1-3)` | 2.9982 | +0.000 | 0.00115 | 4 | (-1.5932, -1.8747) | **0.1253** | L=0 degen (+2) | `8.7e-09` |
| `2D-DEM-022` | `collinear_L0` | `[1, 1, 1.5]` | `AAbaBBaaBAbb` | 13.4697 | +0.000 | 0.01164 | 13 | (-0.7745, -1.8800) | **0.1200** | L=0 degen (+2) | `4.7e-12` |
| `2D-DEM-023` | `collinear_L0` | `[1, 1, 1.1]` | `ABABabaBABAbabA…` | 29.5753 | +0.000 | 0.00550 | 36 | (+1.7818, -1.8975) | **0.1025** | L=0 degen (+2) | `2.7e-11` |
| `2D-DEM-024` | `collinear_L0` | `[1, 1, 0.65]` | `ABABaBABAbABABa…` | 23.5337 | +0.000 | 0.00139 | 34 | (+1.8991, +0.2971) | **0.1009** | L=0 degen (+2) | `3.8e-09` |
| `2D-DEM-025` | `free_fall_brake` | `[1, 1, 0.175]` | `1 (3-2-1-2-1-3)` | 2.9956 | +0.000 | 0.00125 | 4 | (-1.5418, -1.9068) | **0.0932** | L=0 degen (+2) | `7.6e-09` |
| `2D-DEM-026` | `free_fall_brake` | `[1, 0.998, 0.174278568146]` | `1 (3-2-1-2-1-3)` | 3.0006 | +0.000 | 0.00115 | 4 | (-1.5345, -1.9200) | **0.0800** | L=0 degen (+2) | `8.0e-09` |
| `2D-DEM-027` | `collinear_L0` | `[1, 1, 1.05]` | `ABABabaBABAbaba…` | 22.7751 | +0.000 | 0.00470 | 30 | (+1.5106, -1.9219) | **0.0781** | L=0 degen (+2) | `1.5e-11` |
| `2D-DEM-028` | `free_fall_brake` | `[1, 0.998, 0.175]` | `1 (3-2-1-2-1-3)` | 2.9990 | +0.000 | 0.00121 | 4 | (-1.5071, -1.9340) | **0.0660** | L=0 degen (+2) | `7.8e-09` |
| `2D-DEM-029` | `collinear_L0` | `[1, 1, 2]` | `ABABabab` | 12.6489 | +0.000 | 0.04589 | 8 | (-0.8754, -1.9346) | **0.0654** | L=0 degen (+2) | `7.8e-13` |
| `2D-DEM-030` | `free_fall_brake` | `[1, 1, 0.176278568146]` | `1 (3-2-1-2-1-3)` | 2.9928 | +0.000 | 0.00137 | 4 | (-1.4794, -1.9362) | **0.0638** | L=0 degen (+2) | `6.8e-09` |
| `2D-DEM-031` | `collinear_L0` | `[1, 1, 0.9]` | `ABAbaBABabaBAbab` | 18.1242 | +0.000 | 0.00519 | 16 | (+1.8301, -1.9384) | **0.0616** | L=0 degen (+2) | `3.5e-11` |
| `2D-DEM-032` | `free_fall_brake` | `[1, 1, 0.185278568146]` | `1 (3-2-1-2-1-3)` | 2.9659 | +0.000 | 0.00292 | 4 | (-0.3748, -1.9588) | **0.0412** | L=0 degen (+2) | `2.4e-09` |
| `2D-DEM-033` | `free_fall_brake` | `[1, 0.997, 0.175]` | `1 (3-2-1-2-1-3)` | 3.0007 | +0.000 | 0.00119 | 4 | (-1.4740, -1.9605) | **0.0395** | L=0 degen (+2) | `7.8e-09` |
| `2D-DEM-034` | `free_fall_brake` | `[1, 1, 0.177778568146]` | `1 (3-2-1-2-1-3)` | 2.9893 | +0.000 | 0.00153 | 4 | (-1.3919, -1.9646) | **0.0354** | L=0 degen (+2) | `5.8e-09` |
| `2D-DEM-035` | `collinear_L0` | `[1, 1, 0.45]` | `AAbaBAbabABab` | 11.2440 | +0.000 | 0.00194 | 14 | (+1.7420, -1.9655) | **0.0345** | L=0 degen (+2) | `1.2e-10` |
| `2D-DEM-036` | `collinear_L0` | `[1, 1, 1.75]` | `ABABAbabaBABABabab` | 23.9790 | +0.000 | 0.00424 | 18 | (+1.9671, -0.9860) | **0.0329** | L=0 degen (+2) | `1.6e-10` |
| `2D-DEM-037` | `collinear_L0` | `[1, 1, 2]` | `ABABABababab` | 11.7761 | +0.000 | 0.00799 | 12 | (+1.9806, -1.0825) | **0.0194** | L=0 degen (+2) | `1.4e-11` |
| `2D-DEM-038` | `collinear_L0` | `[1, 1, 0.95]` | `ABAbaBABabaBABa…` | 33.4089 | +0.000 | 0.00433 | 28 | (+1.9847, -1.7357) | **0.0153** | L=0 degen (+2) | `3.5e-10` |
| `2D-DEM-039` | `free_fall_brake` | `[1, 0.996, 0.174278568146]` | `1 (3-2-1-2-1-3)` | 3.0041 | +0.000 | 0.00110 | 4 | (-1.4547, -1.9847) | **0.0153** | L=0 degen (+2) | `8.6e-09` |
| `2D-DEM-040` | `collinear_L0` | `[1, 1, 1.1]` | `ABABabaBABababA…` | 20.6946 | +0.000 | 0.01089 | 20 | (+1.9875, +1.4623) | **0.0125** | L=0 degen (+2) | `3.5e-12` |
| `2D-DEM-041` | `collinear_L0` | `[1, 1, 0.5]` | `AAbaBAbabABab` | 16.8187 | +0.000 | 0.01089 | 14 | (+1.9908, -1.3504) | **0.0092** | L=0 degen (+2) | `4.2e-12` |
| `2D-DEM-042` | `free_fall_brake` | `[1, 1, 0.180278568146]` | `1 (3-2-1-2-1-3)` | 2.9829 | +0.000 | 0.00185 | 4 | (-1.1973, -1.9946) | **0.0054** | L=0 degen (+2) | `4.5e-09` |
| `2D-DEM-043` | `collinear_L0` | `[1, 1, 0.65]` | `ABAbaBAbaBABabABab` | 21.2959 | +0.000 | 0.00753 | 18 | (+1.9960, -1.9646) | **0.0040** | L=0 degen (+2) | `1.8e-11` |
| `2D-DEM-044` | `collinear_L0` | `[1, 1, 0.65]` | `ABAbaBABab` | 9.7913 | +0.000 | 0.00340 | 10 | (+1.9968, +1.7562) | **0.0032** | L=0 degen (+2) | `8.2e-12` |
| `2D-DEM-045` | `free_fall_brake` | `[1, 1, 0.182778568146]` | `1 (3-2-1-2-1-3)` | 2.9754 | +0.000 | 0.00228 | 4 | (-0.8997, -1.9972) | **0.0028** | L=0 degen (+2) | `3.4e-09` |
| `2D-DEM-046` | `collinear_L0` | `[1, 1, 1]` | `ABab` | 6.3259 | +0.000 | 0.69053 | 6 | (+1.9972, -0.5952) | **0.0028** | L=0 degen (+2) | `8.1e-16` |
| `2D-DEM-047` | `collinear_L0` | `[1, 1, 1]` | `ABAbaBABab` | 14.8943 | +0.000 | 0.07947 | 10 | (+1.0821, -1.9974) | **0.0026** | L=0 degen (+2) | `2.2e-14` |
| `2D-DEM-048` | `collinear_L0` | `[1, 1, 0.55]` | `AABaBBaabAbb` | 9.1522 | +0.000 | 0.00265 | 12 | (-1.3756, -1.9976) | **0.0024** | L=0 degen (+2) | `1.9e-11` |
| `2D-DEM-049` | `collinear_L0` | `[1, 1, 1.05]` | `ABABabaBABababA…` | 20.9067 | +0.000 | 0.00762 | 20 | (+1.9985, +1.9389) | **0.0015** | L=0 degen (+2) | `2.9e-11` |
| `2D-DEM-050` | `collinear_L0` | `[1, 1, 0.85]` | `ABAbaBABabaBAbab` | 18.4364 | +0.000 | 0.00726 | 16 | (+1.9872, -1.9986) | **0.0014** | L=0 degen (+2) | `1.1e-11` |
| `2D-DEM-051` | `collinear_L0` | `[1, 1, 0.6]` | `ABAbaBABab` | 9.7737 | +0.000 | 0.00432 | 10 | (+1.9989, +1.8298) | **0.0011** | L=0 degen (+2) | `1.3e-11` |
| `2D-DEM-052` | `collinear_L0` | `[1, 1, 0.65]` | `ABAbaBAbaBAbaBA…` | 32.1591 | +0.000 | 0.00579 | 26 | (+1.9989, +1.3315) | **0.0011** | L=0 degen (+2) | `6.2e-11` |
| `2D-DEM-053` | `collinear_L0` | `[1, 1, 3]` | `ABABABABabababa…` | 31.7109 | +0.000 | 0.00266 | 30 | (+1.9998, -1.6523) | **0.0002** | L=0 degen (+2) | `2.2e-09` |
| `2D-DEM-054` | `collinear_L0` | `[1, 1, 4]` | `ABABABABABababa…` | 18.5440 | +0.000 | 0.00510 | 20 | (+1.9999, -1.1270) | **0.0001** | L=0 degen (+2) | `2.9e-10` |
| `2D-DEM-055` | `collinear_L0` | `[1, 1, 1]` | `ABABabaBABAbab` | 13.8671 | +0.000 | 0.00834 | 14 | (+1.9999, -1.4470) | **0.0001** | L=0 degen (+2) | `5.8e-12` |

## 2D Hierarchical, BHH Closure & Co-Orbital Steady Orbits (Top 25 of 62)
Primitive ($q_{\text{inert}}=1$) Broucke–Hadjidemetriou–Hénon (1975–1976) retrograde $q$-fold inertial closures (`A^q`), BHH satellites (Janković & Dmitrašinović 2016), and 1:1 eccentric planetary co-orbitals:

| ID | Class | Masses | $S_3$ Word | $q_{\text{RPO}}$ | Period $T$ | $L$ | $r_{\min}$ | Stability $(k_1, k_2)$ | Margin | 3D Vertical ($k_z$) |
| :--- | :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `2D-HIE-001` | `bhh_qfold_closure` | `[1, 1, 0.5]` | `A^7` | 7 | 24.3591 | -0.843 | 0.9824 | (-0.1184, -0.1275) | **1.8725** | YES (-1.92) |
| `2D-HIE-002` | `bhh_qfold_closure` | `[1, 0.8, 0.6]` | `A^7` | 7 | 26.1898 | -0.805 | 1.0499 | (+0.0576, -0.2971) | **1.7029** | YES (-1.90) |
| `2D-HIE-003` | `bhh_qfold_closure` | `[1, 1, 0.6]` | `A^7` | 7 | 23.0378 | +0.826 | 0.9911 | (-0.0123, -0.3561) | **1.6439** | YES (-1.68) |
| `2D-HIE-004` | `bhh_qfold_closure` | `[1, 1, 2]` | `A^11` | 11 | 23.2637 | +0.818 | 0.8967 | (-0.0760, -0.3656) | **1.6344** | YES (+0.66) |
| `2D-HIE-005` | `bhh_qfold_closure` | `[1, 1, 0.75]` | `A^8` | 8 | 25.2759 | +0.514 | 1.0462 | (+0.4307, +0.3518) | **1.5693** | YES (+1.11) |
| `2D-HIE-006` | `bhh_qfold_closure` | `[1, 1, 0.75]` | `A^7` | 7 | 21.4896 | -0.805 | 1.0010 | (+0.0971, -0.5923) | **1.4077** | YES (-1.24) |
| `2D-HIE-007` | `bhh_qfold_closure` | `[1, 0.8, 0.6]` | `A^7` | 7 | 26.6191 | +0.626 | 1.0430 | (+0.0903, -0.6205) | **1.3795** | YES (-0.82) |
| `2D-HIE-008` | `bhh_qfold_closure` | `[1, 1, 0.3]` | `A^7` | 7 | 27.8692 | -0.889 | 0.9619 | (+0.6621, -0.6401) | **1.3379** | YES (-1.75) |
| `2D-HIE-009` | `bhh_qfold_closure` | `[1, 1, 3]` | `A^8` | 8 | 14.6059 | +0.801 | 0.8059 | (+0.6951, +0.6124) | **1.3049** | YES (+1.51) |
| `2D-HIE-010` | `bhh_qfold_closure` | `[1, 1, 1]` | `A^7` | 7 | 19.3442 | +0.775 | 0.9906 | (+0.2024, -0.7958) | **1.2042** | YES (-0.54) |
| `2D-HIE-011` | `coorbital_1to1` | `[1, 0.6, 0.3]` | `A` | 1 | 16.8823 | +0.859 | 0.6504 | (+0.8107, +0.5115) | **1.1893** | YES (+0.67) |
| `2D-HIE-012` | `bhh_qfold_closure` | `[1, 1, 0.5]` | `A^8` | 8 | 31.2320 | -0.138 | 1.2175 | (+0.8643, +0.0956) | **1.1357** | YES (+1.93) |
| `2D-HIE-013` | `bhh_qfold_closure` | `[1, 1, 1.25]` | `A^7` | 7 | 17.2466 | -0.747 | 0.9155 | (+0.2668, -0.8747) | **1.1253** | YES (-0.01) |
| `2D-HIE-014` | `bhh_qfold_closure` | `[1, 1, 1]` | `A^5` | 5 | 13.6109 | -0.802 | 0.9581 | (+1.0579, +1.0352) | **0.9421** | YES (+0.86) |
| `2D-HIE-015` | `bhh_qfold_closure` | `[1, 1, 1.25]` | `A^5` | 5 | 11.6253 | +0.765 | 0.8951 | (+1.0971, +0.9893) | **0.9029** | YES (+1.14) |
| `2D-HIE-016` | `bhh_qfold_closure` | `[1, 0.9, 0.8]` | `A^5` | 5 | 16.0591 | +0.575 | 1.0549 | (+1.1302, +0.9467) | **0.8698** | YES (+1.29) |
| `2D-HIE-017` | `bhh_qfold_closure` | `[1, 1, 0.75]` | `A^5` | 5 | 14.9522 | -0.826 | 0.9626 | (+1.1413, +1.0011) | **0.8587** | YES (+0.42) |
| `2D-HIE-018` | `bhh_qfold_closure` | `[1, 1, 2]` | `A^5` | 5 | 10.5962 | +0.760 | 0.8476 | (+1.1810, +1.0045) | **0.8190** | YES (+1.54) |
| `2D-HIE-019` | `bhh_qfold_closure` | `[1, 1, 0.6]` | `A^5` | 5 | 17.9147 | +0.275 | 1.1701 | (+1.3227, +0.7859) | **0.6773** | YES (+1.84) |
| `2D-HIE-020` | `coorbital_1to1` | `[1, 0.5, 0.25]` | `A` | 1 | 11.0047 | -0.639 | 0.5110 | (+1.3394, -0.2881) | **0.6606** | YES (+0.95) |
| `2D-HIE-021` | `bhh_qfold_closure` | `[1, 1, 2]` | `A^13` | 13 | 23.3760 | +0.742 | 0.7970 | (-1.3321, -1.3736) | **0.6264** | YES (-0.15) |
| `2D-HIE-022` | `bhh_qfold_closure` | `[1, 1, 0.3]` | `A^7` | 7 | 33.1795 | +0.265 | 1.2718 | (+1.3703, -1.3828) | **0.6172** | YES (+1.57) |
| `2D-HIE-023` | `bhh_qfold_closure` | `[1, 0.6, 0.3]` | `A^5` | 5 | 28.6794 | -0.554 | 1.1870 | (+1.4280, +0.7700) | **0.5720** | YES (+0.28) |
| `2D-HIE-024` | `bhh_qfold_closure` | `[1, 0.8, 0.4]` | `A^11` | 11 | 37.0411 | +0.780 | 0.8746 | (+1.7397, +0.8341) | **0.2603** | YES (+0.26) |
| `2D-HIE-025` | `bhh_qfold_closure` | `[1, 1, 3]` | `A^7` | 7 | 12.7807 | +1.021 | 0.8452 | (+1.8060, +1.6487) | **0.1940** | YES (+1.87) |

## Post-Newtonian (1PN ADM) Relativistic Steady Orbits (21 Verified)
Solved under the exact symplectic 1PN ADM / Einstein–Infeld–Hoffmann Hamiltonian $H = H_0 + \beta H_1$ ($\beta = 1/c^2$, $G=1$; cf. Imai, Chiba & Asada 2007; Lousto & Nakano 2008), audited in 80-bit `long double` at the true stated period $T$ ($T_{\text{inert}}$ with $\theta = 0$ for all snapped $q$-fold closures; $T_{\text{rpo}}$ for RPOs).

| ID | Family | Dim | Masses | $\beta = 1/c^2$ | $\langle v^2\rangle/c^2$ | Period $T$ | Precess $2\alpha$ | True Margin | Symp Err | Return Err | 100T Drift |
| :--- | :--- | :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `1PN-001` | `1pn_3d_delta14_12fold` | 3D | `[1, 1, 0.5]` | 0.00450 | 0.263% | 222.215 | -3.6652 | **0.10701** | `5.3e-07` | `2.7e-10` | 0.00323 |
| `1PN-002` | `1pn_3d_delta14_17fold` | 3D | `[1, 1, 0.5]` | 0.00400 | 0.234% | 314.304 | -3.6960 | **0.09198** | `4.2e-06` | `2.0e-10` | 0.00243 |
| `1PN-003` | `1pn_3d_delta14_5fold` | 3D | `[1, 1, 0.5]` | 0.00100 | 0.060% | 91.707 | -3.7699 | **0.06310** | `1.5e-07` | `5.9e-11` | 0.00614 |
| `1PN-004` | `1pn_3d_delta14_stabilized` | 3D | `[1, 1, 0.5]` | 0.00450 | 0.263% | 18.513 | -3.6833 | **0.04527** | `6.8e-09` | `2.3e-11` | 0.00244 |
| `1PN-005` | `1pn_3d_delta14` | 3D | `[1, 1, 0.5]` | 0.00450 | 0.263% | 18.516 | -3.6713 | **0.03316** | `9.2e-09` | `4.0e-11` | 0.00292 |
| `1PN-006` | `1pn_3d_delta14_stabilized` | 3D | `[1, 1, 0.5]` | 0.00300 | 0.177% | 18.440 | -3.7168 | **0.03205** | `8.1e-09` | `9.0e-12` | 0.00257 |
| `1PN-007` | `1pn_3d_delta14` | 3D | `[1, 1, 0.5]` | 0.00300 | 0.177% | 18.442 | -3.7040 | **0.03025** | `7.4e-09` | `3.8e-12` | 0.00288 |
| `1PN-008` | `1pn_3d_delta14` | 3D | `[1, 1, 0.5]` | 0.00150 | 0.089% | 18.368 | -3.7366 | **0.02082** | `1.2e-08` | `3.9e-12` | 0.00251 |
| `1PN-009` | `1pn_3d_delta14_stabilized` | 3D | `[1, 1, 0.5]` | 0.00150 | 0.089% | 18.366 | -3.7504 | **0.01743** | `1.2e-08` | `7.7e-12` | 0.00274 |
| `1PN-010` | `1pn_fig8_m9995_stabilized` | 2D | `[1, 1, 0.9995]` | 0.00500 | 0.408% | 6.506 | +0.0000 | **1.38778** | `3.7e-09` | `2.1e-11` | 0.11190 |
| `1PN-011` | `1pn_fig8` | 2D | `[1, 1, 1]` | 0.00500 | 0.409% | 6.501 | +0.0000 | **1.38718** | `7.3e-09` | `1.4e-11` | 0.11196 |
| `1PN-012` | `1pn_2d_dem001` | 2D | `[1, 1, 0.5]` | 0.00200 | 0.119% | 24.526 | -0.0742 | **1.29432** | `7.5e-09` | `3.6e-11` | 0.19672 |
| `1PN-013` | `1pn_fig8_m95_stabilized` | 2D | `[1, 1, 0.95]` | 0.00400 | 0.316% | 6.809 | +0.0000 | **0.97219** | `4.9e-09` | `1.8e-11` | 0.14724 |
| `1PN-014` | `1pn_fig8_m95_stabilized` | 2D | `[1, 1, 0.95]` | 0.00600 | 0.466% | 6.846 | +0.0000 | **0.80490** | `8.0e-09` | `2.0e-11` | 0.17630 |
| `1PN-015` | `1pn_2d_dem001` | 2D | `[1, 1, 0.5]` | 0.00400 | 0.236% | 24.696 | -0.1507 | **0.74870** | `1.4e-08` | `3.8e-12` | 0.20289 |
| `1PN-016` | `1pn_fig8` | 2D | `[1, 1, 1]` | 0.00750 | 0.598% | 6.591 | -0.0000 | **0.69665** | `6.9e-09` | `2.8e-11` | 0.16560 |
| `1PN-017` | `1pn_fig8` | 2D | `[1, 1, 1]` | 0.00250 | 0.209% | 6.413 | +0.0000 | **0.50135** | `3.3e-09` | `2.2e-11` | 0.08495 |
| `1PN-018` | `1pn_fig8_m9995_stabilized` | 2D | `[1, 1, 0.9995]` | 0.00250 | 0.209% | 6.421 | +0.0000 | **0.50040** | `4.6e-09` | `2.5e-11` | 0.08494 |
| `1PN-019` | `1pn_2d_dem001_25fold` | 2D | `[1, 1, 0.5]` | 0.00654 | 0.379% | 622.823 | -0.2513 | **0.26145** | `2.6e-05` | `7.4e-10` | 0.21153 |
| `1PN-020` | `1pn_fig8` | 2D | `[1, 1, 1]` | 0.01000 | 0.779% | 6.681 | +0.0000 | **0.03547** | `8.3e-09` | `1.3e-11` | 0.32854 |
| `1PN-021` | `1pn_2d_dem001_10fold` | 2D | `[1, 1, 0.5]` | 0.01516 | 0.835% | 256.760 | -0.6283 | **0.00172** | `6.9e-06` | `4.8e-10` | 0.25397 |

## 4-Body Steady Periodic Orbits (3 Verified Full-3D Linearly Stable)
Audited in 80-bit `long double` via exact $24\times 24$ analytical variational monodromy, canonical symplecticity $\|M_{\text{can}}^\top J M_{\text{can}} - J\|_\infty < 10^{-6}$, full Francis-QR spectrum, and trace-invariant 6th stability index $k_6$ (immune to Jordan-block splitting; 8 weakly hyperbolic 3D candidates removed).

| ID | Family | Dim | Masses | Period $T$ | $r_{\min}$ | Swaps | 24D Symp Err | Max $|\lambda|-1$ | 24D Margin (6/6 Ell.) | Full 3D Stable | 100T Drift |
| :--- | :--- | :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `4B-001` | `4b_2d_dem001_1to1` | 2D | `[1, 1, 0.5, 0.001]` | 24.3165 | 0.9811 | 14 | `1.5e-12` | `2.31e-14` | **0.007334** | YES | 0.45843 |
| `4B-002` | `4b_2d_dem001_1to1` | 2D | `[1, 1, 0.5, 0.003]` | 24.2309 | 0.9785 | 14 | `1.1e-12` | `2.30e-14` | **0.004210** | YES | 0.45011 |
| `4B-003` | `4b_2d_dem001_1to1` | 2D | `[1, 1, 0.5, 0.005]` | 24.1438 | 0.9759 | 14 | `1.3e-12` | `2.30e-14` | **0.001869** | YES | 0.44323 |

## Mass-Continuation Stability Windows ($m_1 = m_2 = 1$, varying $m_3$)
Segmented on branch jumps ($|\Delta T| > 1.0$ or $|\Delta \mathbf{v}| > 0.10$) and unstable gaps so only **contiguous same-branch linearly stable sub-intervals** on the primary seed branch are reported.

### Braid-1.10 (m3 around 1.10)
- **Primary-branch traced range:** $m_3 \in [1.0200,\; 1.1600]$ (29 pts) (primary branch; 4 jump segment(s) excluded)
- **Contiguous linearly stable sub-interval(s):** $[1.0450,\; 1.1350]$ (19 pts, peak 1.0999)
- **Peak stability:** at $m_3^* = 1.0750$ with margin **1.0999** ($k_1 = +0.8738,\; k_2 = -0.9001$, $T = 27.177367$)

### Heavy-Center-3.00 (m3 around 3.00)
- **Primary-branch traced range:** $m_3 \in [2.8200,\; 3.6000]$ (40 pts)
- **Contiguous linearly stable sub-interval(s):** $[2.8400,\; 3.0600]$ (12 pts, peak 0.8753), $[3.1000,\; 3.4200]$ (17 pts, peak 1.6079)
- **Peak stability:** at $m_3^* = 3.3000$ with margin **1.6079** ($k_1 = -0.3824,\; k_2 = -0.3921$, $T = 16.741939$)

### Heavy-Braid-1.75 (m3 around 1.75)
- **Primary-branch traced range:** $m_3 \in [1.7200,\; 1.7800]$ (5 pts)
- **Contiguous linearly stable sub-interval(s):** $[1.7200,\; 1.7500]$ (3 pts, peak 1.5913)
- **Peak stability:** at $m_3^* = 1.7350$ with margin **1.5913** ($k_1 = +0.4087,\; k_2 = -0.2474$, $T = 13.799897$)

### Light-Braid-0.55 (m3 around 0.55)
- **Primary-branch traced range:** $m_3 \in [0.5500,\; 0.5650]$ (4 pts) (primary branch; 1 jump segment(s) excluded)
- **Contiguous linearly stable sub-interval(s):** $[0.5500,\; 0.5650]$ (4 pts, peak 0.5090)
- **Peak stability:** at $m_3^* = 0.5600$ with margin **0.5090** ($k_1 = -1.1704,\; k_2 = -1.4910$, $T = 30.556849$)

### Figure-8 (m3 around 1.00)
- **Primary-branch traced range:** $m_3 \in [0.9850,\; 1.0000]$ (31 pts)
- **Contiguous linearly stable sub-interval(s):** $[1.0000,\; 1.0000]$ (1 pts, peak 0.0028)
- **Peak stability:** at $m_3^* = 1.0000$ with margin **0.0028** ($k_1 = +1.9972,\; k_2 = -0.5952$, $T = 6.325914$)

### Free-Fall Brake (m3 around 0.174)
- **Primary-branch traced range:** $m_3 \in [0.1558,\; 0.1878]$ (65 pts)
- **Contiguous linearly stable sub-interval(s):** $[0.1558,\; 0.1813]$ (52 pts, peak 0.2663), $[0.1823,\; 0.1878]$ (12 pts, peak 0.2051)
- **Peak stability:** at $m_3^* = 0.1693$ with margin **0.2663** ($k_1 = -1.7274,\; k_2 = -1.7337$, $T = 3.006813$)


## Exact Initial Conditions (Top Showcase Orbits)
### `3D-001` — `3d_choreo_helical` (`choreography_3d_helical`), masses `[1, 1, 1]`, $T = 15.543227271965$
- **Margin:** `0.247654`, **Return err:** `1.13e-10`, **Symp err:** `1.56e-12`
- **`y0`:** `[-1, 1.550285591865, 0, 1, 1.550285591865, 0, 0, -3.100571183731, 0, -0.1296499913256, 0.6420183647434, 0.06, -0.1296499913256, -0.6420183647434, 0.06, 0.2592999826512, 0, -0.12]`

### `3D-002` — `3d_inertial_5fold` (`bhh_3d_subharmonic`), masses `[1, 1, 0.5]`, $T = 91.464179118420$
- **Margin:** `0.200169`, **Return err:** `1.37e-10`, **Symp err:** `3.42e-11`
- **`y0`:** `[-1, 0, 0, 0.9479359541975, 0, 0, 0.104128091605, 0, 0, 0, 0.2298078731035, 0, 0, -0.4553004786205, 0.6734752656928, 0, 0.450985211034, -1.346950531386]`

### `3D-003` — `3d_vertical_pd_7fold` (`bhh_3d_subharmonic`), masses `[1, 1, 0.415]`, $T = 335.219354730900$
- **Margin:** `0.102887`, **Return err:** `3.76e-10`, **Symp err:** `1.60e-09`
- **`y0`:** `[-1, 0, 0, 0.9255864672514, 0, 0, 0.1793097174665, 0, 0, 0, 0.2244366458536, 0, 0, -0.6774303109278, 0.4191081742431, 0, 1.091551000179, -1.009899215044]`

### `3D-004` — `3d_choreo_helical_24fold` (`choreography_3d_helical`), masses `[1, 1, 1]`, $T = 373.170027060800$
- **Margin:** `0.098344`, **Return err:** `3.63e-11`, **Symp err:** `5.28e-11`
- **`y0`:** `[-1, 1.552047596447, 0, 1, 1.552047596447, 0, 0, -3.104095192894, 0, -0.1305269595906, 0.6423243546499, 0.05647873181652, -0.1305269595906, -0.6423243546499, 0.05647873181652, 0.2610539191811, 0, -0.112957463633]`

### `3D-005` — `3d_choreo_helical_11fold` (`choreography_3d_helical`), masses `[1, 1, 1]`, $T = 170.949280091550$
- **Margin:** `0.094767`, **Return err:** `3.86e-11`, **Symp err:** `4.80e-11`
- **`y0`:** `[-1, 1.549525685916, 0, 1, 1.549525685916, 0, 0, -3.099051371831, 0, -0.1292721854049, 0.641886394804, 0.06145193034628, -0.1292721854049, -0.641886394804, 0.06145193034628, 0.2585443708098, 0, -0.1229038606926]`

### `3D-006` — `3d_inertial_12fold` (`bhh_3d_subharmonic`), masses `[1, 1, 0.5]`, $T = 219.697222998700$
- **Margin:** `0.082745`, **Return err:** `2.48e-10`, **Symp err:** `4.34e-10`
- **`y0`:** `[-1, 0, 0, 0.9489288586447, 0, 0, 0.1021422827106, 0, 0, 0, 0.1770020266891, 0, 0, -0.5362198938625, 0.6164617513438, 0, 0.7184357343468, -1.232923502688]`

### `3D-007` — `3d_rpo_vertical_pd` (`bhh_3d_subharmonic`), masses `[1, 1, 0.415]`, $T = 47.858628766730$
- **Margin:** `0.046637`, **Return err:** `1.36e-10`, **Symp err:** `1.39e-11`
- **`y0`:** `[-1, 0, 0, 0.925551337217, 0, 0, 0.1793943681517, 0, 0, 0, 0.2347540779857, 0, 0, -0.654233000299, 0.452134728954, 0, 1.010792583888, -1.089481274588]`

### `3D-008` — `3d_choreo_helical_15fold` (`choreography_3d_helical`), masses `[1, 1, 1]`, $T = 233.456950433450$
- **Margin:** `0.019950`, **Return err:** `2.88e-11`, **Symp err:** `3.39e-11`
- **`y0`:** `[-1, 1.556853825916, 0, 1, 1.556853825916, 0, 0, -3.113707651832, 0, -0.1329257756037, 0.6431589252141, 0.0454091230715, -0.1329257756037, -0.6431589252141, 0.0454091230715, 0.2658515512074, 0, -0.090818246143]`

### `3D-009` — `3d_vertical_pd_11fold` (`bhh_3d_subharmonic`), masses `[1, 1, 0.415]`, $T = 526.565200439700$
- **Margin:** `0.018977`, **Return err:** `6.76e-10`, **Symp err:** `5.28e-09`
- **`y0`:** `[-1, 0, 0, 0.9255193441271, 0, 0, 0.1794714599347, 0, 0, 0, 0.2294907539663, 0, 0, -0.6659597842357, 0.4360602801258, 0, 1.051732603059, -1.050747662954]`

### `3D-010` — `3d_vertical_pd_8fold` (`bhh_3d_subharmonic`), masses `[1, 1, 0.415]`, $T = 385.470267624000$
- **Margin:** `0.016723`, **Return err:** `5.96e-10`, **Symp err:** `4.82e-09`
- **`y0`:** `[-1, 0, 0, 0.927665190062, 0, 0, 0.1743007468386, 0, 0, 0, 0.1967049195836, 0, 0, -0.7448463306793, 0.2841493146552, 0, 1.320822677339, -0.6846971437475]`

### `2D-DEM-001` — `2d_choreo_5petal_5fold` (`choreography_2d_rot`), masses `[1, 1, 1]`, $T = 68.416042643200$
- **Margin:** `1.589721`, **Return err:** `2.64e-10`, **Symp err:** `2.43e-10`
- **`y0`:** `[-1, 0.8818291939664, 1, 0.8818291939664, 0, -1.763658387933, -0.1816437405585, 0.587293465153, -0.1816437405585, -0.587293465153, 0.3632874811169, 0]`

### `2D-DEM-002` — `perp_L_nonzero` (`shuttle_exchange`), masses `[1, 0.9, 0.8]`, $T = 15.495141531780$
- **Margin:** `1.206048`, **Return err:** `3.96e-09`, **Symp err:** `5.79e-11`
- **`y0`:** `[-1, 0, -0.05951905193104, 0, 1.316958933422, 0, 0, 0.7812388652309, 0, -0.1613813787143, 0, -0.7949945304851]`

### `2D-DEM-003` — `collinear_L0` (`braid_L0_exchange`), masses `[1, 1, 1.75]`, $T = 13.745509466030$
- **Margin:** `0.880029`, **Return err:** `3.29e-09`, **Symp err:** `1.01e-09`
- **`y0`:** `[-1, 0, 1, 0, 0, 0, 0.1644176545216, 0.9211164935806, 0.1644176545216, 0.9211164935806, -0.1879058908819, -1.052704564092]`

### `2D-DEM-004` — `2d_choreo_8swap_3fold` (`choreography_2d_rot`), masses `[1, 1, 1]`, $T = 31.393053773730$
- **Margin:** `0.821382`, **Return err:** `3.63e-09`, **Symp err:** `1.28e-08`
- **`y0`:** `[-1, 0.2499223071152, 1, 0.2499223071152, 0, -0.4998446142305, 0.2098053508812, 0.5199031067188, 0.2098053508812, -0.5199031067188, -0.4196107017624, 0]`

### `2D-DEM-005` — `2d_choreo_fig8_2fold` (`choreography_2d_rot`), masses `[1, 1, 1]`, $T = 22.635128707730$
- **Margin:** `0.702694`, **Return err:** `9.19e-11`, **Symp err:** `5.61e-11`
- **`y0`:** `[-1, 0.7702753302814, 1, 0.7702753302814, 0, -1.540550660563, -0.1893779654758, 0.5698797625933, -0.1893779654758, -0.5698797625933, 0.3787559309515, 0]`

### `2D-DEM-006` — `2d_choreo_fig8_4fold` (`choreography_2d_rot`), masses `[1, 1, 1]`, $T = 73.200063705070$
- **Margin:** `0.594435`, **Return err:** `3.19e-10`, **Symp err:** `1.11e-10`
- **`y0`:** `[-1, 1.080173098475, 1, 1.080173098475, 0, -2.160346196949, -0.1670484345861, 0.6115848926766, -0.1670484345861, -0.6115848926766, 0.3340968691721, 0]`

### `2D-DEM-007` — `2d_choreo_fig8_12fold` (`choreography_2d_rot`), masses `[1, 1, 1]`, $T = 159.073997034800$
- **Margin:** `0.591879`, **Return err:** `1.72e-09`, **Symp err:** `1.47e-09`
- **`y0`:** `[-1, -0.8623273098799, 1, -0.8623273098799, 0, 1.72465461976, 0.1830627910545, 0.5844308183197, 0.1830627910545, -0.5844308183197, -0.3661255821089, 0]`

### `2D-DEM-008` — `collinear_L0` (`braid_L0_exchange`), masses `[1, 1, 0.55]`, $T = 30.145715748740$
- **Margin:** `0.477593`, **Return err:** `2.77e-09`, **Symp err:** `2.88e-13`
- **`y0`:** `[-1, 0, 1, 0, 0, 0, 0.3009665258277, 0.2527187743142, 0.3009665258277, 0.2527187743142, -1.094423730282, -0.9189773611427]`

### `2D-DEM-009` — `collinear_L0` (`braid_L0_exchange`), masses `[1, 1, 0.6]`, $T = 24.555889462460$
- **Margin:** `0.451750`, **Return err:** `1.15e-08`, **Symp err:** `7.65e-10`
- **`y0`:** `[-1, 0, 1, 0, 0, 0, 0.198933087672, 0.09648911003918, 0.198933087672, 0.09648911003918, -0.6631102922401, -0.3216303667973]`

### `2D-DEM-010` — `collinear_L0` (`braid_L0_exchange`), masses `[1, 1, 0.75]`, $T = 20.487096418130$
- **Margin:** `0.426619`, **Return err:** `1.39e-08`, **Symp err:** `2.83e-10`
- **`y0`:** `[-1, 0, 1, 0, 0, 0, 0.190246168878, 0.1147396346484, 0.190246168878, 0.1147396346484, -0.507323117008, -0.3059723590624]`

### `2D-DEM-011` — `collinear_L0` (`braid_L0_exchange`), masses `[1, 1, 1.1]`, $T = 26.160569266350$
- **Margin:** `0.335436`, **Return err:** `3.75e-09`, **Symp err:** `1.23e-12`
- **`y0`:** `[-1, 0, 1, 0, 0, 0, 0.3658879746975, 0.5486388556991, 0.3658879746975, 0.5486388556991, -0.6652508630864, -0.9975251921803]`

### `2D-DEM-012` — `collinear_L0` (`braid_L0_exchange`), masses `[1, 1, 1]`, $T = 28.669270915580$
- **Margin:** `0.302189`, **Return err:** `3.47e-09`, **Symp err:** `2.40e-13`
- **`y0`:** `[-1, 0, 1, 0, 0, 0, 0.4391659178758, 0.4529676431928, 0.4391659178758, 0.4529676431928, -0.8783318357516, -0.9059352863856]`

### `2D-HIE-001` — `perp_L_nonzero` (`bhh_qfold_closure`), masses `[1, 1, 0.5]`, $T = 24.359107221400$
- **Margin:** `1.872483`, **Return err:** `9.95e-11`, **Symp err:** `5.12e-12`
- **`y0`:** `[-1, 0, 0.9941458062522, 0, 0.01170838749561, 0, 0, 0.07446872307688, 0, -0.7818621047497, 0, 1.414786763346]`

### `2D-HIE-002` — `perp_L_nonzero` (`bhh_qfold_closure`), masses `[1, 0.8, 0.6]`, $T = 26.189849194790$
- **Margin:** `1.702879`, **Return err:** `2.75e-10`, **Symp err:** `5.73e-12`
- **`y0`:** `[-1, 0, 1.164247281754, 0, 0.1143369576612, 0, 0, 0.7328402560966, 0, 0.01399748905863, 0, -1.240063745572]`

### `2D-HIE-003` — `perp_L_nonzero` (`bhh_qfold_closure`), masses `[1, 1, 0.6]`, $T = 23.037812500710$
- **Margin:** `1.643869`, **Return err:** `1.87e-10`, **Symp err:** `4.54e-12`
- **`y0`:** `[-1, 0, 0.9966642760782, 0, 0.005559539869611, 0, 0, -0.00317809355291, 0, 0.8298856875384, 0, -1.377845989976]`

