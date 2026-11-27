# Post-PR #2 Scientific & Numerical Review of `main` (`01d4f7e`) & `threebodyorbits.com` Eligibility Audit

**Target Revision:** `main` at `01d4f7e` (*Merge pull request #2 from gawydion/feat/next-phase-roadmap*, merging `2ea8423` + `8ac4005`, 165 cataloged orbits: 133 Newtonian 3-body + 21 1PN ADM relativistic + 11 Newtonian 4-body).

**Questions Addressed:**
1. *Did PR #2 (`8ac4005`) fully resolve the prior audit, and do all five newly introduced roadmap tracks (Tracks 1–5) hold up under independent 80-bit `long double` verification and literature cross-checking?*
2. *Can any of the genuinely new families on `main` be posted to [threebodyorbits.com](https://www.threebodyorbits.com), and if so, which ones and via what process?*

All numbers in this review are reproduced by the zero-dependency verification suite in [`review/`](review/) (`bash review/reproduce.sh --long-kick`).

---

## Executive Verdict

### 1. State of `main` (`01d4f7e`)

Commit `8ac4005` cleanly fixed **every** critical Newtonian 3-body issue flagged in the earlier audit (PR #1): all **133** Newtonian 3-body orbits now pass 80-bit analytical variational auditing (`src/audit_catalog.cpp`), have `q_inert = 1` (zero multiple covers), `symp_err <= 2.55e-8`, `ret_err < 1e-7`, `stab_margin > 1e-4`, and `0/117` S3 word mismatches in `scripts/classify_topology.py --verify`. The 1PN ADM Hamiltonian in `src/search_1pn.cpp` also matches the exact Damour–Schäfer / Lousto–Nakano 1PN ADM Hamiltonian term-for-term, and all **21** 1PN orbits are genuinely linearly stable in independent 80-bit verification.

However, **four substantive issues remain across the newly added roadmap tracks (Tracks 1–5)**:

| Severity | Scope | Finding Summary |
| :--- | :--- | :--- |
| **Critical** | **8 orbits** (`4B-001`..`4B-008`) | **All 8 3D 4-body orbits are weakly hyperbolically unstable ($k_6 > 2$), not linearly stable.** In independent 80-bit analytical variational 24×24 monodromy (`symp_err = 4.6e-13..1.8e-12`), every 3D 4-body orbit has a real reciprocal hyperbolic pair ($\lambda \in [1.000199, 1.002409]$, $k_6 - 2 \in [+2.38\times 10^{-8}, +5.80\times 10^{-6}]$), confirmed empirically on `4B-008` where a $10^{-6}$ kick amplifies by **2.2 million×** over 4,000 periods. Conversely, the **3 planar 4-body orbits (`4B-009`..`4B-011`)** pass full 24×24 80-bit audit and are genuinely **full-3D linearly stable**. |
| **High** | **10 orbits** (`choreography_2d_rot` & `choreography_3d_helical`) | **Rediscoveries of classical Chenciner–Féjoz–Montgomery (2005) / Marchal (2000) families + closure mislabels.** `fig8_rot_2d.dat` and `fig8_helical_3d.dat` continue directly from the Figure-8 into the celebrated $\Gamma_z$ and $\Gamma_y$ ($P_{12}$) families (missing from `NEW_FAMILIES.md` references). `choreo_12swap_2d.dat` is the mirror reflection ($y \to -y, v_x \to -v_x$) of `fig8_rot_2d.dat`; `2D-DEM-007` (`12swap_6fold`) is a **12-fold** $\Gamma_z$ closure `(ABab)^12`; `3D-004` (`12fold`) is a **24-fold** closure (`precess = -11π/12`, `q_rpo = 24`). |
| **Medium** | **14 orbits** (`free_fall_brake`: `2D-DEM-014`..`045`) | **Coarse uniform sampling (`NS = 2400`) in `src/audit_catalog.cpp` misreports `rmin` (by up to 20.3×) and misses 8 of 12 syzygies.** Dense adaptive DOPRI5 tracing shows `2D-DEM-020` has true `rmin = 0.00053` (catalog: `0.01077`), `2D-DEM-014` has true `rmin = 0.00076` (catalog: `0.00298`), and all 14 orbits have **12 syzygies** with identical half-period Tanikawa–Mikkola sequence **`3-2-1-2-1-3`**. All 14 are samples of **one** continuous collisionless mass family (cleanest showcase: `2D-DEM-015` at $m_3 = 0.1878$, true `rmin = 0.00445`). |
| **Medium** | **5 orbits** (`1PN-002, 006, 009, 019, 020`) | **`src/search_1pn.cpp` writes 1-RPO monodromy & `swaps *= q_fold` into snapped $T_{\text{inert}}$ closure records.** All 21 1PN orbits are genuinely linearly stable in 80-bit verification, but the 5 snapped closures report 1-RPO stability margins instead of true $T_{\text{inert}}$ margins — notably `1PN-019` (10-fold), whose reported margin `0.37710` overstates the true $T_{\text{inert}}$ margin `0.00172` by **219×**. Also missing citations to Imai, Chiba & Asada (2007) and Lousto & Nakano (2008). |

### 2. Posting New Families to `threebodyorbits.com`

- **Submission channel:** `threebodyorbits.com` has **no manual upload form or public submission API**. Its interactive [`/hunt`](https://www.threebodyorbits.com/hunt) runs server-defined browser WASM searches (feeding `/garden`), whereas the curated atlas (~3,900+ orbits, CC BY 4.0) is maintained by its creator (`u/laddi_macchiato` on Reddit, **`info@threebodyorbits.com`**). Posting external discoveries requires sending a verified dataset directly to **`info@threebodyorbits.com`**.
- **Scope:** Strictly **planar (2D), Newtonian, 3-body, exact inertial-frame periodic orbits** (mostly $L = 0$ collinear and free-fall, plus choreographies and selected $L \neq 0$ planar sequences such as Broucke and Sheen).
- **Best candidates to pitch:**
  1. **`free_fall_brake`** (`2D-DEM-015`, `2D-DEM-014`, `2D-DEM-026` — half-period syzygy sequence **`3-2-1-2-1-3`**): fits the site's dedicated unequal-mass Free-fall section natively.
  2. **`AAAbbb` shuttle-exchange** (`2D-DEM-002` at `[1, 0.9, 0.8]`, `2D-DEM-019` at `[1, 0.8, 0.6]`): strongest non-hierarchical $L \neq 0$ planar candidate (`rmin ≈ 0.19`, full-3D stable).
  3. **6 unmatched intermediate-$m_3$ `collinear_L0` braids** (`2D-DEM-009, 022, 024, 048, 053, 054`): direct format match, though incremental and worth checking against `/garden` first.
- **Not eligible:** All 3D families (`3D-001`..`016`), all 1PN relativistic orbits (`1PN-001`..`021`), and all 4-body orbits (`4B-001`..`011`) are outside the site's planar Newtonian 3-body scope; `choreography_2d_rot` (`2D-DEM-001, 005, 006, 007, 018`) are closures of the already-published $\Gamma_z$ rotating Figure-8.

---

## 1. Critical Audit of Track 5: 4-Body Orbits (`4B-001`..`4B-011`)

### Root Cause in `src/search_4body.cpp` and `scripts/build_catalog.py`
1. `src/search_4body.cpp:325-411` computes the 24×24 monodromy matrix by `double` central finite differences (`eps = 1e-7`) without analytical variational equations or canonical symplecticity checks, and accepts orbits with `max_rad_err < 2e-3`, `max_triv < 5e-3`, and `max_k < 2.0 - 1e-6`.
2. `scripts/build_catalog.py:518-526` skips `audit_catalog` and **does not enforce `stab_margin > 1e-4`** on `four_body`. Consequently, 7 of 8 3D 4-body orbits were admitted with catalog margins of `2e-6`–`2.1e-5` — **50× to 100× smaller** than their own `trivial_eig_err` (`3.9e-4`–`2.6e-3`).
3. On `4B-001`, `max_eig_rad_err = 7.18e-4` and `trivial_eig_err = 1.51e-3`; the greedy trivial-eigenvalue stripper (`search_4body.cpp:362-373`) matched a split eigenvalue `> 1` as trivial and paired two `< 1` eigenvalues (`0.99928 + 0.99928 = 1.99856`), fabricating a bogus `stab_margin = 0.001439`.

### Independent 80-Bit Analytical Variational Audit (`review/audit_4b.cpp`)
Integrating the exact 24×24 analytical variational equations $\dot{M} = J_f(\mathbf{y}) M$ in 80-bit `long double` (`tol = 1e-16L`, canonical symplecticity $\|M_{\text{can}}^\top J M_{\text{can}} - J\|_\infty \in [4.6\times 10^{-13}, 1.8\times 10^{-12}]$) and computing the exact trace-invariant 6th stability index $k_6 = \operatorname{tr}(M) - (6 + 6\cos\theta) - \sum_{j=1}^5 k_j$ (immune to Jordan-block splitting of the trivial $+1$ block) yields:

| ID | Dim | Family (`subspace`) | Catalog Margin | 80-Bit Symp Err | 80-Bit $\max|\lambda|$ | Exact Trace $k_6 - 2$ (or 6/6 Ell. Margin) | Independent Verdict |
| :--- | :---: | :--- | :---: | :---: | :---: | :---: | :--- |
| `4B-001` | 3D | `4b_3d_dem001_r21_vz35` | `0.001439` | `4.6e-13` | `1.0003418` | **`+1.17e-07`** | **HYPERBOLIC UNSTABLE** |
| `4B-002` | 3D | `4b_3d_double_binary_N8_i35` | `0.000021` | `1.8e-12` | `1.0002986` | **`+9.11e-08`** | **HYPERBOLIC UNSTABLE** |
| `4B-003` | 3D | `4b_3d_dem001_r21_vz25` | `0.000009` | `4.6e-13` | `1.0001985` | **`+2.38e-08`** | **HYPERBOLIC UNSTABLE** |
| `4B-004` | 3D | `4b_3d_double_binary_N7_i35` | `0.000008` | `1.4e-12` | `1.0002398` | **`+5.80e-08`** | **HYPERBOLIC UNSTABLE** |
| `4B-005` | 3D | `4b_3d_double_binary_N7_i36` | `0.000008` | `1.4e-12` | `1.0002253` | **`+5.13e-08`** | **HYPERBOLIC UNSTABLE** |
| `4B-006` | 3D | `4b_3d_double_binary_N8_i33` | `0.000008` | `1.8e-12` | `1.0002688` | **`+7.30e-08`** | **HYPERBOLIC UNSTABLE** |
| `4B-007` | 3D | `4b_3d_dem001_r21_vz30` | `0.000004` | `1.3e-12` | `1.0002457` | **`+2.82e-08`** | **HYPERBOLIC UNSTABLE** |
| `4B-008` | 3D | `4b_3d_double_binary_N5_i40` | `0.000002` | `5.8e-13` | **`1.0024085`** | **`+5.80e-06`** | **HYPERBOLIC UNSTABLE** |
| `4B-009` | 2D | `4b_2d_dem001_1to1` ($m_4=0.001$) | `0.007334` | `1.5e-12` | `1.0000000` | `margin = 0.007334` (6/6 ell.) | **FULL-3D LINEARLY STABLE** |
| `4B-010` | 2D | `4b_2d_dem001_1to1` ($m_4=0.003$) | `0.004211` | `1.1e-12` | `1.0000000` | `margin = 0.004210` (6/6 ell.) | **FULL-3D LINEARLY STABLE** |
| `4B-011` | 2D | `4b_2d_dem001_1to1` ($m_4=0.005$) | `0.001868` | `1.3e-12` | `1.0000000` | `margin = 0.001869` (6/6 ell.) | **FULL-3D LINEARLY STABLE** |

Notice the exact algebraic consistency $(\lambda - 1)^2 \approx k_6 - 2$ (e.g., for `4B-008`, $(2.4085\times 10^{-3})^2 = 5.801\times 10^{-6}$ vs trace-invariant $k_6 - 2 = 5.804\times 10^{-6}$, five orders of magnitude above the $5.8\times 10^{-13}$ symplecticity floor). Empirical 4,000-period integration of `4B-008` from a $10^{-6}$ momentum-preserving kick confirms steady exponential blowup: `drift(500T) = 1.28e-2` → `drift(2000T) = 2.09e-1` → `drift(4000T) = 2.21` (a **2.2 million×** growth missed by the repo's 100-period test).

---

## 2. Audit of Track 3: Rotating 2D & 3D Choreographies (`src/search_choreo_rot.cpp`)

1. **Classical literature origin (missing citations):**
   - `data/choreo_rot/fig8_rot_2d.dat` begins directly at the planar Figure-8 ($\alpha = 0$, $T = 31.1815$, $k_1 = 1.9972$, $k_2 = -0.5952$) and continues in rotation angle $\alpha$. This is the classical **$\Gamma_z$ planar rotating Figure-8 family** proved by **Chenciner, Féjoz & Montgomery (2005, *Nonlinearity* 18, 1407)** (stability mapped by Broucke, Elipe & Riaguas 2006; Kapela & Simó 2007).
   - `data/choreo_rot/fig8_helical_3d.dat` likewise begins directly at the planar Figure-8 ($v_z = 0.002$, $T = 31.1814$, $k_1 = 1.99717$, $k_2 = -0.59544$) and tilts out of the plane with half-period RPO symmetry ($T^*_{\text{rpo}} = 4.645 \approx \frac{1}{2} T^*_{\text{Fig-8}}$). This is the celebrated **Marchal (2000) $P_{12}$ / Chenciner–Féjoz–Montgomery (2005) $\Gamma_y$** 3D rotating-eight family connecting the Figure-8 to the Lagrange relative equilibrium.
   - Neither Chenciner–Féjoz–Montgomery (2005), Marchal (2000), Broucke–Elipe–Riaguas (2006), nor Kapela–Simó (2007) is cited in `NEW_FAMILIES.md`, which lists all 10 orbits under `"Noteworthy / New"`.
2. **Duplicate/mirror branch and mislabeled fold counts:**
   - `data/choreo_rot/choreo_12swap_2d.dat` ($\alpha \in [0.7922, 1.1922]$, $y_0 \in [-0.9079, -0.6991]$) is the exact spatial mirror reflection ($y \to -y, v_x \to -v_x$) of the same $\Gamma_z$ branch in `fig8_rot_2d.dat`. Consequently, **`2D-DEM-007`** (labeled `2d_choreo_12swap_6fold`) is a **12-fold** closure `(ABab)^12` (`q_rpo = 12`, $T_{\text{rpo}} = 13.256$, $T^*_{\text{rpo}} = 9.1144$, $\Lambda = 0.1957$) sitting right beside `2D-DEM-001` ($T^*_{\text{rpo}} = 9.1244$, $\Lambda = 0.1865$) on the same $\Gamma_z$ curve — not a separate multi-loop choreography.
   - **`3D-004`** (labeled `3d_choreo_helical_12fold`) has `precess/π = -0.916667 = -11/12`. After 12 RPO periods, the accumulated rotation is $-11\pi \equiv \pi \pmod{2\pi}$, so inertial closure requires **24** RPO periods (`q_rpo = 24`, $T = 373.1700$). It is a **24-fold** closure, not 12-fold.
   - **`2D-DEM-004`** (`2d_choreo_8swap_3fold`, `q_rpo = 3`, $T^*_{\text{rpo}} = 23.1500$, $\Lambda = 1.7648$, `true_rmin = 0.09491`) is distinct from $\Gamma_z$; its 18 syzygies occur in immediately cancelling pairs `11-33-22-11-33-22-11-33-22`, giving contractible $S_3$ word `1` (warranting cross-check against Simó's rotating choreographies).

---

## 3. Audit of Track 1: Free-Fall Brake Family (`2D-DEM-014`..`045`)

All 14 `free_fall_brake` entries are genuinely linearly stable in 80-bit analytical variational arithmetic (`symp_err = 1.5e-9..2.5e-8`, `ret_err <= 1.5e-9`) and belong to **one** continuous 2-parameter mass family ($m = [1, 0.996\text{–}1.000, 0.163\text{–}0.188]$).

However, `src/audit_catalog.cpp:1189-1233` measures `rmin` and `syzygies` on a fixed uniform grid of `NS = 2400` points (`dt = T / 2400 ≈ 0.00125`), stepping right across the `dt ≈ 8e-5` binary slingshot passages of Body 3 around Bodies 1 and 2:

| ID | Masses | Catalog `rmin` | True `rmin` (Adaptive) | Overstatement | Catalog `syz` | True `syz` | Half-Period Syzygy Word |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| `2D-DEM-014` | `[1, 1, 0.16778]` | `0.00298` | **`0.00076`** | 3.9× | `4` | **`12`** | `321213` |
| `2D-DEM-015` | `[1, 1, 0.18778]` | `0.00769` | **`0.00445`** | 1.7× | `12` | **`12`** | `321213` |
| `2D-DEM-020` | `[1, 1, 0.16278]` | `0.01077` | **`0.00053`** | **20.3×** | `4` | **`12`** | `321213` |
| `2D-DEM-026` | `[1, 0.998, 0.17428]` | `0.00601` | **`0.00115`** | 5.2× | `4` | **`12`** | `321213` |
| `2D-DEM-030` | `[1, 1, 0.17628]` | `0.01075` | **`0.00137`** | 7.8× | `4` | **`12`** | `321213` |

**Classification note:** Every time-reversible brake orbit automatically has contractible reduced $F_2$ word $w \cdot w^{-1} = 1$. In free-fall literature (Li & Liao 2019; Hristov et al. 2024/2025), brake families are identified by their **half-period Tanikawa–Mikkola syzygy sequence** — here consistently **`3-2-1-2-1-3`** (six syzygies from rest at $t=0$ to rest at $t=T/2$, then retracing `3-1-2-1-2-3`). Because `2D-DEM-015` ($m_3 = 0.1878$) maintains `true_rmin = 0.00445` (nearly 6× wider than `2D-DEM-014`) with margin `0.2051`, it is the cleanest collisionless representative of this family.

---

## 4. Audit of Track 4: 1PN ADM Relativistic Orbits (`1PN-001`..`1PN-021`)

1. **Physics & linear stability confirmed:** Direct inspection of `src/search_1pn.cpp:25-148` confirms the 1PN ADM Hamiltonian matches Damour–Schäfer / Lousto–Nakano (2008) exactly. Our independent 80-bit `long double` canonical audit (`review/audit_1pn.cpp`) confirms that **all 21** 1PN orbits preserve $H_{\text{1PN}}$ to $1.8\times 10^{-15}\text{–}2.1\times 10^{-13}$, close to $3.2\times 10^{-12}\text{–}7.2\times 10^{-10}$, satisfy canonical symplecticity (`symp_err = 2.5e-9..1.1e-8` for all 16 RPOs; `1.6e-7..2.9e-5` for long $q$-fold closures up to $T = 622.8$), and are **genuinely linearly stable** — confirming the notable GR stabilization of Newtonian-unstable unequal-mass Figure-8 ($m_3 = 0.95$, `1PN-013/014`) and Newtonian-unstable `3D-Δ14` (`1PN-001/004/008`).
2. **Metadata bug on all 5 snapped inertial closures:** In `src/search_1pn.cpp:1096-1105`, `1158-1167`, and `1228-1238`, `--snap-2d` and `--snap-3d` call `analyze_stab_1pn` and `measure_geom` at `T_rpo` rather than `T_inert = q_fold * T_rpo` (and multiply `swaps *= q_fold`), writing 1-RPO indices into the $T_{\text{inert}}$ records:

| ID | Family | $q$ | $T_{\text{inert}}$ | Catalog Margin (1-RPO) | **True $T_{\text{inert}}$ Margin (80-Bit)** | True $T_{\text{inert}}$ $(k_1, k_2, k_3/k_z)$ | Symp Err |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| `1PN-002` | `1pn_3d_delta14_17fold` | 17 | `314.30` | `0.04095` | **`0.09197`** | `(+1.5633, -1.2311, -1.9080)` | `4.5e-6` |
| `1PN-006` | `1pn_3d_delta14_12fold` | 12 | `222.21` | `0.02375` | **`0.10701`** | `(-0.5533, -1.4304, -1.8930)` | `5.5e-7` |
| `1PN-009` | `1pn_3d_delta14_5fold` | 5 | `91.71` | `0.00254` | **`0.06310`** | `(+1.9369, +1.3037, -0.6449)` | `1.6e-7` |
| `1PN-019` | `1pn_2d_dem001_10fold` | 10 | `256.76` | **`0.37710`** *(219× high!)* | **`0.00172`** | `(+1.9983, +1.5561), kz=-0.0781` | `6.7e-6` |
| `1PN-020` | `1pn_2d_dem001_25fold` | 25 | `622.82` | `0.18304` | **`0.26145`** | `(+0.4297, -1.7385), kz=+1.8945` | `2.9e-5` |

3. **Missing literature citations:** Cite **Imai, Chiba & Asada (2007, *Phys. Rev. Lett.* 98, 201102)** (1PN choreographic Figure-8) and **Lousto & Nakano (2008, *Class. Quantum Grav.* 25, 145018)** (1PN/2PN ADM three-body equations of motion) in `NEW_FAMILIES.md`.

---

## 5. Posting Eligibility Assessment for `https://www.threebodyorbits.com`

### How `threebodyorbits.com` Accepts Orbits
- **No manual submission form or public API:** The site does not provide a web form or REST endpoint to upload external initial conditions directly.
- **Automated browser discovery (`/hunt` → `/garden`):** `https://www.threebodyorbits.com/hunt` runs a client-side Rust/WASM worker scanning server-defined planar Newtonian subspaces (e.g., $1:\mu:1$). Loops closed by the worker undergo server verification; unique braid words appear in `/garden` and grant naming rights.
- **Direct curator submission (`info@threebodyorbits.com` / `u/laddi_macchiato`):** The main atlas (~3,900+ orbits, CC BY 4.0) is compiled by its creator from published papers and external compilations (Šuvakov–Dmitrašinović, Li–Liao 2017/2018/2019, Hristov et al., Simó, Broucke, Sheen 2016, Reusser). Submitting external families requires contacting **`info@threebodyorbits.com`** directly with high-precision (ideally ≥16–32 digit) initial conditions, mass ratios, period, topology/syzygy classification, and stability/verification data.
- **Technical scope:** **Planar (2D), Newtonian, 3-body, exact inertial-frame periodic orbits** (predominantly $L = 0$ collinear and free-fall, plus choreographies and notable $L \neq 0$ planar exceptions like Broucke and Sheen).

### Family-by-Family Eligibility Matrix

| `main` Family | Fits Atlas Scope? | Worth Pitching to `info@threebodyorbits.com`? | Pre-Submission Requirements |
| :--- | :---: | :--- | :--- |
| **`free_fall_brake`** (`2D-DEM-015`, `014`, `026`) | **Yes** (2D, Newtonian 3B, $L=0$, free-fall) | **Best candidate** | Matches the site's ~311-orbit unequal-mass Free-fall section. Lead with `2D-DEM-015` ($m=[1,1,0.1878]$, clean `true_rmin = 0.00445`) and `2D-DEM-026` ($m=[1,0.998,0.1743]$). Include half-period syzygy word **`3-2-1-2-1-3`**, re-converge `y0` to ≥16–32 digits, and cross-check against Li & Liao (2019) / Hristov et al. (2024–2025). |
| **`shuttle_exchange` `AAAbbb`** (`2D-DEM-002`, `019`) | **Yes** (2D, Newtonian 3B, inertial periodic, $L \neq 0$) | **Strong candidate** | Non-hierarchical (`q_rpo = 1`), clean `rmin ≈ 0.19`, full-3D linearly stable, includes `[1, 0.8, 0.6]` mass ratio already present on the site. Cross-check older Broucke/Šuvakov $L \neq 0$ archives first. |
| **6 intermediate-$m_3$ `collinear_L0` braids** (`2D-DEM-009, 022, 024, 048, 053, 054`) | **Yes** (native $1:\mu:1$, $L=0$ format) | **Yes (incremental)** | Direct match to the bulk of the atlas; check whether their words already appear in `/garden` or Hristov et al. (2025). |
| **`choreography_2d_rot`** (`2D-DEM-001, 004, 005, 006, 007, 018`) | **Yes** (2D, Newtonian 3B, $L \neq 0$) | **Low priority** | `001/005/006/007/018` are closures of the known Chenciner–Féjoz–Montgomery (2005) $\Gamma_z$ family; cross-check `2D-DEM-004` against Simó's choreographies already in the atlas. |
| **3D Newtonian** (`3D-001`..`016`) | **No** (3D; several are RPOs) | **No** | Outside the site's planar 2D scope (`3D-001/004/005/008` are also Marchal $P_{12}$ / CFM $\Gamma_y$). |
| **1PN (`1PN-001`..`021`) & 4-Body (`4B-001`..`011`)** | **No** (non-Newtonian / 4-body) | **No** | Outside the site's Newtonian 3-body scope (and `4B-001`..`008` are hyperbolically unstable). |

---

## 6. Actionable Recommendations for `main`

1. **Audit 4-body in 80-bit analytical variational arithmetic & remove `4B-001`..`4B-008`:** Upgrade `src/search_4body.cpp` (or add a 4-body mode to `src/audit_catalog.cpp`) to use 80-bit analytical variational monodromy, canonical symplecticity `< 1e-6`, and `stab_margin > 1e-4` in `scripts/build_catalog.py`. Keep `4B-009`..`4B-011` (verified full-3D stable) and re-search 3D 4-body islands with analytical variational stability.
2. **Track `rmin` and syzygies inside adaptive steps in `src/audit_catalog.cpp`:** Replace the coarse `NS = 2400` uniform loop with per-adaptive-step tracking so close-encounter free-fall orbits report exact `rmin` and all 12 syzygies (plus half-period Tanikawa–Mikkola sequence for brake orbits).
3. **Evaluate true $T_{\text{inert}}$ monodromy for snapped 1PN closures in `src/search_1pn.cpp`:** Run `analyze_stab_1pn` and `measure_geom` at `T_inert` (with `two_alpha = 0`), add canonical symplecticity checking, and cite Imai–Chiba–Asada (2007) and Lousto–Nakano (2008).
4. **Update Track 3 attribution & labels in `NEW_FAMILIES.md` / `CATALOG.md`:** Cite Chenciner–Féjoz–Montgomery (2005), Marchal (2000), Broucke–Elipe–Riaguas (2006), and Kapela–Simó (2007) for $\Gamma_z$ and $\Gamma_y/P_{12}$; rename `2D-DEM-007` (`12fold` on $\Gamma_z$) and `3D-004` (`24fold` on $\Gamma_y/P_{12}$).
5. **Prepare a high-precision submission bundle for `info@threebodyorbits.com`:** Package `free_fall_brake` (`2D-DEM-015`, `014`, `026` with half-period word `3-2-1-2-1-3`), `AAAbbb` (`2D-DEM-002`, `019`), and the 6 intermediate-$m_3$ `collinear_L0` orbits re-converged to ≥16–32 digits.

---

## 7. Reproducing This Review

From the repository root (requires only `g++` with C++20 and `python3` stdlib):

```bash
bash review/reproduce.sh --long-kick
```
