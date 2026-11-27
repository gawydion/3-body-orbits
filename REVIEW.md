# Rigorous Scientific & Numerical Review of `3-body-orbits`

**Question:** *Are the claimed "new" 3-body periodic orbits in this repository genuinely new, linearly stable, and correctly classified?*

**Review Method:** Every one of the **124** cataloged orbits in [`catalog.json`](catalog.json) was independently audited using an **80-bit `long double` Dormand–Prince 5(4) integrator with exact analytical variational equations** (eliminating finite-difference truncation error), **direct Hessenberg–Francis QR eigendecomposition** of the full $12\times 12$, $6\times 6$, and $18\times 18$ monodromy matrices (with explicit canonical symplecticity verification $\|M_{\text{can}}^T J M_{\text{can}} - J\|_\infty$), **rotating-frame Procrustes relative-period and permutation-symmetry reduction**, **adaptive equator-shifted Montgomery shape-sphere word computation** across all $S_3$ puncture projections, **invariable-plane kinematic decomposition**, **2,000-period empirical perturbation tests**, and **direct cross-matching against published catalogs** — including Broucke (1975), Hadjidemetriou (1975), Hénon (1976, 1977), Šuvakov & Dmitrašinović (2013), Janković & Dmitrašinović (2016), Li & Liao (2017, 695 equal-mass orbits), Li, Jing & Liao (2018, 1,349 unequal-mass orbits), Li, Li & Liao (2021, 135,445 non-hierarchical orbits), and Li & Liao (2025, 10,059 spatial 3D orbits).

All verification code to reproduce every number in this report is included in [`review/`](review/).

---

## Executive Verdict

**Most of the orbits highlighted as "new families" are rediscoveries of classical or previously published families — primarily rational $q$-fold inertial closures of the 1975–1976 Broucke–Hadjidemetriou–Hénon (BHH) retrograde family, published Li & Liao (2017, 2018) $L = 0$ families, and exact multiple covers.** In addition, the verification pipeline contains severe numerical and topological bugs that admit **13 `collinear_L0` orbits violating the catalog's stability criteria** (**5 confirmed hyperbolically unstable**, **4 ultra-near-parabolic with $|k_1 - 2| \lesssim 10^{-6}$**, and **4 sub-threshold with margin $< 10^{-4}$**), corrupt **44 of 74 (59.5%) free-group words**, and mislabel hierarchical retrograde and planetary orbits as "democratic."

That said, the computational infrastructure is impressive and **a small, specific subset of orbits warrants genuine scientific interest** once stripped of overclaimed novelty and bugs:

1. **`3D-Δ14` (`3D-001`, `3D-002`, `3D-003`, `3D-005`) — Classical parent, potentially novel explicit high-order branch:**
   - **What it actually is:** Mathematically proven below to be the **5:2 vertical subharmonic bifurcation** ($\nu_z = 2/5$, coupled to horizontal mode $\nu_1 \approx 2/5$) of the **1975–1976 Broucke–Hadjidemetriou–Hénon (BHH) planar retrograde family**, NOT an independent "3D democratic" family. Dynamically, it is a **retrograde hierarchical triple** (permanent retrograde inner binary $(1,3)$ with $r_{13} \in [1.038, 1.309]$ plus nearly circular prograde outer body $2$ with $R_2 \in [0.890, 0.976]$); the "14 swaps" are shallow distance grazes, not partner exchanges.
   - **What IS real and valuable:** All 4 catalog points are genuine, high-precision, linearly stable 3D solutions (`ret_err ~ 1e-10`, $\max|\lambda| = 1.000000$). This specific 5:2 vertical branch at $m = [1, 1, 0.5]$ does **not** appear in Li & Liao's (2025) 10,059-orbit spatial database (which searches a different isosceles symmetry subspace), and high-order ($n=5$) vertical subharmonic bifurcations of unequal-mass BHH orbits are rarely tabulated explicitly in the literature (earlier studies such as Hénon 1977 and Ichtiaroglou & Michalodimitrakis 1980 focused on $n=1, 2$ vertical critical orbits). However, `3D-001/002/003/005` are **4 points along a single 1-parameter slice**, not 4 distinct families.
2. **`2D-DEM-006` & `2D-DEM-024` (`m = [1, 0.9, 0.8]` and `[1, 0.8, 0.6]`) — Strongest genuinely non-hierarchical $L \neq 0$ candidate:**
   - Misclassified as `AABaaB` by a bug in `classify_topology.py`; its true $S_3$-canonical free-group word is **`AAAbbb`** ($a^{-3} b^3$, inverse-equivalent to $a^3 b^{-3}$). Body 1 genuinely shuttles back and forth, making 3 close prograde revolutions around Body 2 ($r_{12,\min} = 0.23$) and 3 around Body 3 ($r_{13,\min} = 0.19$). Continuing it to equal masses $m = [1, 1, 1]$ yields an **unstable** $T/2$-symmetric orbit ($k_1 = 3.569$), meaning unequal masses **stabilize** it into a deep full-3D KAM island (`margin = 1.2061`, $k_z = -0.2568$). It does not appear in Li & Liao (2017, 2018, 2021), though priority still requires checking older $L \neq 0$ Broucke/Hadjidemetriou/Šuvakov–Dmitrašinović archives.
3. **6 genuinely stable $L = 0$ collinear words not verbatim in Li & Liao (2017, 2018):**
   - `2D-DEM-016` ($m_3=0.60$, len 32), `2D-DEM-025` ($m_3=0.65$, len 34), `2D-DEM-026` ($m_3=1.50$, len 12), `2D-DEM-045` ($m_3=3.00$, len 30, true margin $1.6\times 10^{-4}$), `2D-DEM-054` ($m_3=0.55$, len 12), and `2D-DEM-061` ($m_3=4.00$, $(AB)^5(ab)^5$). Five of these six sit at intermediate $m_3$ values not sampled in Li, Jing & Liao (2018), which only tabulated $m_3 \in \{0.5, 0.75, 2, 4, 5, 8, 10\}$.

---

## 1. Complete Reconciled Audit of All 124 Catalog Orbits

Every orbit in `catalog.json` falls into exactly one row below (totals sum to **124**):

| Catalog Bucket | Count | True Dynamical / Literature Classification | Genuinely New? |
| :--- | :---: | :--- | :--- |
| **3D — `3D-Δ14` slice** (`3D-001, 002, 003, 005`) | **4** | 4 points along **one** 1-parameter slice ($m=[1,1,0.5]$) of the **5:2 vertical subharmonic bifurcation** of the 1975–1976 planar **BHH retrograde family** (retrograde hierarchical, not democratic) | Classical parent & mechanism; explicit $n=5$ branch at $m_3=0.5$ not in Li & Liao (2025) |
| **3D — Distinct mass** (`3D-004`) | **1** | Standard inclined 3:1 prograde hierarchical binary (`swaps = 0`, $r_{23}\in[0.557, 0.641]$) | **No** (standard hierarchical triple) |
| **2D Hierarchical** (`2D-HIE-001..045`) | **45** | **45/45 (100%)** are $q$-fold rational inertial closures ($q\in\{3,4,5,6,7,8,11\}$, $S_3$-word $A^q$) of the 1975–1976 **BHH family**; **7** are exact 2× inertial multiple covers (`012, 013, 014, 015, 016, 018, 019`), 4 duplicating catalog entries | **No** (BHH family + 7 double-covers) |
| **2D `perp` "Democratic" — BHH closures** | **20** | $q$-fold rational inertial closures ($q\in\{2,5,7,8,9,13\}$, $S_3$-word $A^q$) of the 1975–1976 **BHH retrograde family** — includes ALL `A^7` (`2D-DEM-001..004, 008`) and `(Ab)^7` (`2D-DEM-013`) showcase orbits | **No** (1975–1976 BHH family) |
| **2D `perp` "Democratic" — BHH satellites** (`2D-DEM-029, 063`) | **2** | Period-multiplying BHH topological satellites ($q_{\text{rpo}}=1$, words $a^4, a^7$) of the class discovered by **Janković & Dmitrašinović (2016, *PRL* 116, 064301)** | **No** (known BHH satellite class) |
| **2D `perp` "Democratic" — 1:1 co-orbital** (`2D-DEM-007, 011`) | **2** | Classic 1:1 anti-aligned eccentric planetary/co-orbital orbits around $m_1=1$ ($r_{23}\in[2.7, 5.6]$, word $A$); mislabeled "democratic" due to alternating perihelion | **No** (classic 1:1 eccentric planetary resonance) |
| **2D `perp` "Democratic" — Shuttle exchange** (`2D-DEM-006, 024`) | **2** | Two mass points of **one** genuinely non-hierarchical, 3D-stable binary-exchange shuttle family; true $S_3$-word **`AAAbbb`** ($a^{-3} b^3$, misclassified as `AABaaB`) | **Candidate** (not in Li & Liao 2017–2021; check Broucke/Šuvakov $L\neq 0$ archives) |
| **2D `collinear_L0` — Multiple cover** (`2D-DEM-064`) | **1** | Exact **5× inertial multiple cover of the Figure-8** (`2D-DEM-067`, $T = 5 \times 6.3259 = 31.6296$) | **No** (duplicate of Moore 1993 / Chenciner–Montgomery 2000) |
| **2D `collinear_L0` — Published words** | **38** | Exact word match to **Šuvakov & Dmitrašinović (2013)**, **Li & Liao (2017)**, and/or **Li, Jing & Liao (2018)** (includes `2D-DEM-009, 015, 018, 020, 022, 023, 031, 067, 068`; **3** of these 38 are **confirmed unstable**, **4** are **ultra-near-parabolic** $|k_1-2|\lesssim 10^{-6}$, and **3** fail `margin > 1e-4`) | **No** (already published) |
| **2D `collinear_L0` — Unmatched, unstable / sub-threshold** (`2D-DEM-035, 037, 049`) | **3** | `037` is a hyperbolic 3× cover/period-tripling of `(AB)^3(ab)^3` ($k_1 = 2.000803$); `049` is hyperbolic $(AB)^4(ab)^4$ ($k_1 = 2.000002$); `035` fails `margin > 1e-4` ($4.8\times 10^{-5}$) | **No** (unstable / fails catalog criteria) |
| **2D `collinear_L0` — Unmatched, verified stable** (`2D-DEM-016, 025, 026, 045, 054, 061`) | **6** | Verified linearly stable ($\text{margin} > 10^{-4}$) $L=0$ words not found verbatim in Li & Liao (2017, 2018); 5 are at intermediate $m_3$ not scanned in LL2018, and `061` is $(AB)^5(ab)^5$ at $m_3=4$ | **Minor / incremental** (unequal-mass $L=0$ extensions) |
| **TOTAL** | **124** | — | — |

---

## 2. Deep Audit of `3D-Δ14` (`3D-001`, `3D-002`, `3D-003`, `3D-005`) and `3D-004`

### 2.1 Exact Mathematical Identification: 5:2 Vertical Subharmonic Bifurcation of the 1975–1976 Planar BHH Retrograde Family

Transforming `3D-001` ($m = [1, 1, 0.5]$, $T_{\text{rpo}} = 18.292871$, $2\alpha = -3.769217\text{ rad}$) into the **invariable plane** ($\hat{\mathbf{u}} = \mathbf{L}/|\mathbf{L}| = \hat{\mathbf{z}}'$, tilted $41.69^\circ$ from the solver's arbitrary $z$-axis) reveals its exact dynamical architecture:

- **Body 1 ($m_1=1$) and Body 3 ($m_3=0.5$)** form a **permanently bound, nearly circular retrograde inner binary**:
  $$r_{13}(t) \in [1.038,\; 1.309] \quad (r_{13,\max}/r_{13,\min} = 1.26,\; e \approx 0.11), \qquad \Delta\phi_{13} = -2.600\text{ turns.}$$
- **Body 2 ($m_2=1$)** moves on a **nearly circular prograde outer orbit** around the $(1,3)$ binary center of mass:
  $$R_{2,xy}(t) \in [0.890,\; 0.976] \quad (e \approx 0.04), \qquad \Delta\phi_{12} = \Delta\phi_{23} = +2.400\text{ turns.}$$
- **Synodic resonance count:** In the frame co-rotating with outer Body 2, the inner retrograde binary completes **exactly $5$ synodic revolutions** per relative period ($\Delta\phi_{12} - \Delta\phi_{13} = +2.400 - (-2.600) = 5.000$), while executing vertical oscillations (`z_cross = 4` for Body 2, `6` for Bodies 1 & 3).

This is the **exact same retrograde hierarchical architecture** as `2D-DEM-001` ($m = [1, 1, 0.5]$, $7$-fold planar BHH closure: $r_{13} \in [1.01, 1.13]$, $R_2 \in [0.98, 0.99]$, $\Delta\phi_{13} = -4$, $\Delta\phi_{12} = +3$, $7$ synodic turns). Continuing the **1975–1976 planar BHH retrograde family** at $m = [1, 1, 0.5]$ from `2D-DEM-001` ($\Lambda = 0.7036$) toward `2D-DEM-056` ($\Lambda = 0.6765$) proves that its single-period vertical Floquet tune $\nu_z = \arccos(k_z^{(1)}/2)/(2\pi)$ crosses **$\nu_z = 2/5 = 0.40000$** ($k_z^{(1)} = 2\cos(4\pi/5) = -\varphi \approx -1.618034$, so $k_z^{(5)} = +2.0000$) at `vy2 = 0.81560`. Comparing the 5th iterate of that planar BHH orbit with `3D-Δ14` as its invariable-plane out-of-plane velocity $v_{z2,\text{inv}} \to 0$ gives an exact quantitative match across **every invariant and all three Floquet modes**:

| Quantity | Planar BHH 5:2 Vertical Resonance ($5 \times T_1$, `vy2=0.8156`) | `3D-Δ14` (`vz2_inv = -0.1997`, `3D-001`) | `3D-Δ14` approaching planar limit (`vz2_inv = -0.0790`) |
| :--- | :---: | :---: | :---: |
| $x_2$ | `0.95395` | `0.94793` | `0.95167` $\to$ `0.9540` |
| $v_{y1,\text{inv}}$ | `-0.10020` | `-0.17134` | `-0.11227` $\to$ `-0.1002` |
| $v_{y2,\text{inv}}$ | `+0.81560` | `+0.78809` | `+0.81339` $\to$ `+0.8156` |
| Period $T$ | `18.3456` ($5 \times 3.6691$) | `18.2929` | `18.3330` $\to$ `18.3456` |
| Precession $2\alpha$ | `-3.4921 rad` ($5\cdot 2\alpha_1 - 6\pi$) | `-3.7692 rad` | `-3.5548 rad` $\to$ `-3.4921` |
| $\Lambda = |\mathbf{L}|\sqrt{-E}$ | `0.6799` | `0.7134` | `0.6847` $\to$ `0.6799` |
| $T^* = T|E|^{3/2}$ | `10.753` ($5 \times 2.1505$) | `10.656` | `10.731` $\to$ `10.753` |
| Floquet $k_2^{(5)}$ | **`-1.9775`** ($\nu_2 = 0.4952$) | `-1.4901` | **`-1.9686`** $\to$ `-1.9775` |
| Floquet $k_1^{(5)}, k_3^{(5)}$ | **`+1.8488`, `+2.0000`** ($\nu_1 = 0.3875, \nu_z = 0.4000$) | `+1.9915`, `-0.0346` | **`+2.003`, `+2.000`** (5:2 horiz–vert resonance coupling!) |

**Physical interpretation:** By Poincaré–Birkhoff / Birkhoff–Lewis bifurcation theory (Hénon 1973, 1977; Ichtiaroglou & Michalodimitrakis 1980), whenever a vertically elliptic planar family's vertical tune $\nu_z$ sweeps monotonically across rationals $m/n$, an $n$-fold 3D subharmonic resonant branch bifurcates at every rational $\nu_z = m/n$. Here $\nu_z = 2/5$ occurs right next to the horizontal resonance $\nu_1 \approx 2/5$, producing a narrow 5:2 resonant 3D island (`k1 ∈ [1.987, 1.999]`, stability margin $\le 0.0132$) coupled to the planar BHH parent.

### 2.2 Why `3D-Δ14` Is Not "Democratic"

The codebase defines `democratic = (partner_swaps >= 2)` (`src/campaign_2d.cpp:186`, `src/rpo_3d_fixedvz.cpp:433`). Because the inner retrograde binary $(1,3)$ has radius $r_{13} \in [1.04, 1.31]$ and outer Body 2 orbits at $R_2 \approx 0.94$ (distance $\sim 1.5$ from the $(1,3)$ COM), $r_{23}$ briefly dips to $0.844$ ($< r_{13}$) during each of the 5 synodic conjunctions — generating $2 \times 5 = 10$ transient closest-pair flips (`(1,3) → (2,3) → (1,3)`). For $v_{z2} > 0.583$, $r_{12}$ (`min = 1.188`) also grazes slightly below $r_{13}$ (`max = 1.309`) twice per period, adding 4 brief flips ($10 + 4 = 14$). **No pair ever unbinds or exchanges membership.** Below $v_{z2} = 0.583$, `r12` never becomes the closest pair at all (`closest_counts = [0, 1851, 549]` on `3D-005`).

### 2.3 Factual & Numerical Errors in `3D-Δ14` and 3D Catalog Reporting

1. **Inertial closures (`3D-002`, `3D-003`, `3D-005`) have copied/hand-estimated metadata instead of computed inertial Floquet indices:**
   In `data/discoveries_3d_inertial.jsonl`, the 1-RPO-period stability indices (or rough eyeball interpolations from `fine_vz_dn.dat`), approximate $E, L$, and copied `drift_100T` were written into the $q$-fold inertial records rather than computing the true $T_{\text{inert}}$-period monodromy matrix:

   | ID | Closure | Catalog $(k_1, k_2, k_3)$ & Margin | Exact $T_{\text{inert}}$ $(k_1, k_2, k_3)$ & Margin | Catalog vs Exact Swaps | Catalog vs Exact $E$ |
   | :--- | :--- | :---: | :---: | :---: | :---: |
   | `3D-001` | 1-RPO | `(+1.9915, -0.0346, -1.4901)` · `0.0085` | `(+1.9915, -0.0346, -1.4901)` · `0.0085` | 14 vs 14 | `-0.69752` vs `-0.69752` |
   | `3D-002` | 5-fold | `(+1.9915, -0.0346, -1.4901)` · `0.0085` *(copied from `3D-001`)* | **`(+1.7998, +1.7155, -0.2375)` · `0.2002`** | 70 vs 70 | `-0.69752` vs `-0.69751` |
   | `3D-003` | 12-fold | `(+1.9968, +1.5200, -1.9180)` · `0.0032` | **`(+1.5672, -1.1865, -1.9173)` · `0.0827`** | 168 vs 168 | `-0.69830` vs `-0.69872` |
   | `3D-005` | 7-fold | `(+1.9990, +1.8200, -1.9450)` · `0.0010` *(wrong `vz2` row!)* | **`(+1.9839, +0.8157, -0.3864)` · `0.0161`** | **98 vs 70** *(hardcoded $7\times 14$!)* | `-0.69870` vs `-0.69941` |

   *(Note: All three inertial closures ARE genuinely closed to $1.4\times 10^{-10}\text{–}2.5\times 10^{-10}$ and genuinely linearly stable with $\max|\lambda| = 1.000000$; their true inertial margins are actually larger than the 1-RPO values entered in the table.)*

2. **Upper stability boundary misidentified (`NEW_FAMILIES.md:41`):**
   `NEW_FAMILIES.md` states the upper boundary at $v_{z2} \approx 0.678$ is a *"Hamiltonian–Hopf (Krein collision) where $k_2, k_3$ collide on the unit circle."* However, the repo's own `data/family_3d/fine_vz_up.dat` line 2 (`vz2 = 0.6782`) shows **`k1 = +2.006053 > 2`, `k2 = -1.265276`, `k3 = -0.331393`** — stability is lost via a **tangent bifurcation ($k_1$ crossing $+2$)** BEFORE any Krein collision occurs between `0.6782` and `0.6832`.

3. **`3D-004` (`m = [1, 0.7, 0.4]`) is a standard 3:1 prograde hierarchical binary:**
   In the invariable plane, $(2,3)$ is a tight prograde binary ($r_{23} \in [0.557, 0.641]$, $+3.935$ turns) orbited by outer Body 1 ($R_1 \in [0.999, 1.046]$, $+0.935$ turns, exact 3:1 synodic resonance, `swaps = 0`, `democratic = false`).

4. **Scrambled IDs across `NEW_FAMILIES.md`, `scripts/generate_gifs.py`, and `catalog.json`:**
   - `NEW_FAMILIES.md` & `generate_gifs.py` label **5-fold** as `3D-001`, **7-fold** as `3D-004`, and **distinct-mass `[1,0.7,0.4]`** as `3D-005`.
   - In `catalog.json` & `CATALOG.md`, `3D-001` is the 1-RPO orbit ($T=18.29$), `3D-002` is 5-fold ($T=91.46$), `3D-004` is `[1,0.7,0.4]` ($T=11.03$), and `3D-005` is 7-fold ($T=128.27$).
   - Consequently, running `scripts/generate_gifs.py` renders the **wrong orbits** into all three 3D GIFs!

---

## 3. Audit of All 71 Planar $L \neq 0$ (`perp_L_nonzero`) Orbits

### 3.1 65 of 71 Orbits Are Rational Inertial Closures of the 1975–1976 BHH Family

On the Montgomery shape sphere (which quotients out planar rotations $\mathrm{SO}(2)$), any $q$-fold rational inertial closure ($T = q T_{\text{rel}}$, precession $\Delta\theta = 2\pi p / q$) of a period-1 relative periodic orbit with word $a$ automatically produces the word **$a^q$** (or $(Ab)^q$ if the binary collision is placed at $\infty$ in stereographic projection).

Our Procrustes rotating-frame period scan proves that **100% of `2D-HIE` (45/45)** and **76.9% of `2D-DEM` `perp_L_nonzero` (20/26)** satisfy $T = q T_{\text{rel}}$ with exact closure ($< 10^{-9}$) at $T_{\text{rel}} = T / q$, and all have $S_3$-canonical word $A^q$:

| Showcase / Group | Catalog IDs | $q_{\text{rpo}}$ ($2\pi p / q$) | $T^*_{\text{rel}}$ & $\Lambda$ | True Physical Identity |
| :--- | :--- | :---: | :---: | :--- |
| **`A^7` Deepest Island** | `2D-DEM-001, 002, 003, 004, 008` | `7` ($3/7$) | `1.46–4.23`, `0.50–0.89` | **1975–1976 BHH retrograde family** ($3/7$ inertial closure). On `2D-DEM-001`, $(1,3)$ is a permanent retrograde binary ($r_{13}\in[1.01,1.13]$, $-4$ turns) orbited by circular Body 2 ($R_2\in[0.98,0.99]$, $+3$ turns). Never exchanges partners! |
| **`(Ab)^7` Wide-Separation** | `2D-DEM-013` (`m=[1,1,0.3]`) | `7` ($-3/7$, perm `14`) | `1.62`, `0.185` | **1975–1976 BHH retrograde circumbinary family**! Permanent prograde inner binary $(1,2)$ ($r_{12}\in[1.48,1.52]$, $+4$ turns) + circular retrograde outer Body 3 ($R_3\in[1.71,1.75]$, $-3$ turns). Under $S_3$ relabeling, `(Ab)^7` $\equiv$ `A^7`. |
| **Other `A^q` closures** | `2D-DEM-005` ($q=8$), `010, 014` ($q=5$), `012` ($q=13$), `048, 052, 056, 058, 059, 060, 066, 071, 072` ($q=2$), `073` ($q=9$) | `2, 5, 8, 9, 13` | e.g. `059` ($m=[1,1,1]$): $T^*_{\text{rel}}=3.664, \Lambda=0.784$ | **1975–1976 BHH retrograde family** (`2D-DEM-059` is literally Hénon's 1976 equal-mass 1:1 retrograde orbit; `073` is the $4/9$ closure on the same curve!). |
| **All `2D-HIE`** | `2D-HIE-001..045` | `3, 4, 5, 6, 7, 8, 11` | BHH retrograde/prograde | **45/45** are BHH $q$-fold closures; **7** (`012, 013, 014, 015, 016, 018, 019`) are exact **2× inertial multiple covers**, with **4** duplicating catalog entries (`012→022`, `013→024`, `015→031`, `016→032`). |

*(Note: Why did `2D-DEM-001` have "668× deeper stability margin than Figure-8"? Because Hénon 1976 already showed that the planar BHH retrograde family possesses a massive linear stability island — and iterating $M_{\text{rel}}^7$ happened to rotate both elliptic Floquet angles near $\pm\pi/2$, yielding $k_1, k_2 \approx -0.12$.)*

### 3.2 The Remaining 6 `perp_L_nonzero` Orbits (`q_rpo = 1`)

- **`2D-DEM-029` ($a^4$, $m=[1,1,0.5]$) & `2D-DEM-063` ($a^7$, $m=[1,1,1]$):** Genuine period-multiplying ($q_{\text{rpo}} = 1$) **BHH topological satellites** — precisely the $a^k$ satellite family discovered and analyzed by **Janković & Dmitrašinović (2016, *Phys. Rev. Lett.* 116, 064301)** (see also Dmitrašinović et al. 2018).
- **`2D-DEM-007` ($m=[1,0.6,0.3]$) & `2D-DEM-011` ($m=[1,0.5,0.25]$):** Word `A`. Bodies 2 and 3 orbit the heaviest mass $m_1 = 1$ on anti-aligned **1:1 eccentric planetary/co-orbital orbits** ($r_{23} \in [3.91, 5.65]$ and $[2.73, 4.04]$; both complete $+1$ turn around $m_1$ and $0$ net turns around each other). Mislabeled "democratic" only because Bodies 2 and 3 alternate reaching perihelion around $m_1$.
- **`2D-DEM-006` ($m=[1,0.9,0.8]$) & `2D-DEM-024` ($m=[1,0.8,0.6]$):** Two mass points along **one genuinely non-hierarchical binary-exchange shuttle family**. Body 1 makes 3 close prograde revolutions around Body 2 ($r_{12,\min} = 0.232$), transfers to Body 3 for 3 close prograde revolutions ($r_{13,\min} = 0.190$), and returns, while $r_{23} \in [1.29, 1.95]$. Its true $S_3$-canonical word is **`AAAbbb`** ($a^{-3} b^3$), NOT `AABaaB`. Continuing it to $m = [1, 1, 1]$ gives a $T/2$-symmetric equal-mass orbit ($T^* = 27.925, \Lambda = 2.372$) that is **linearly unstable** ($k_1 = +3.5689, k_2 = +0.9233$), demonstrating genuine mass-asymmetry stabilization into a full-3D stable island (`margin = 1.2061`, $k_z = -0.2568$). **This is the most scientifically interesting 2D orbit in the repository.**

*(Stability verification note: Unlike the `collinear_L0` subset, all 71 `perp_L_nonzero` orbits are genuinely linearly stable in both 2D and 3D; our analytical variational indices match `catalog.json` to $\max |\Delta k| = 2.8\times 10^{-4}$.)*

---

## 4. Audit of All 48 Planar $L = 0$ (`collinear_L0`) Orbits

### 4.1 Literature Cross-Match Against Šuvakov–Dmitrašinović (2013) and Li & Liao (2017, 2018)

Using our corrected adaptive `long double` Montgomery free-group words (Section 5.2), we cross-matched all **48** `collinear_L0` orbits directly against the published datasets of **Li & Liao (2017)** (`three-body-free-group-word.md`, 695 equal-mass orbits) and **Li, Jing & Liao (2018)** (`three-body-unequal-mass-free-group-word.md`, 1,349 unequal-mass orbits):

- **38 of 48** match published free-group words in Li & Liao (2017, 2018), including **every single `collinear_L0` orbit featured in `NEW_FAMILIES.md` and `CATALOG.md` continuations**:
  - `2D-DEM-067` ($m_3=1$, `ABab`): **Figure-8** (Moore 1993; Chenciner & Montgomery 2000; `I.A^{i.c.}_1`).
  - `2D-DEM-064` ($m_3=1$, `(ABab)^5`): **Exact 5× multiple cover of the Figure-8** ($T = 31.6296 = 5 \times 6.3259$, $k_1 = 2\cos(5\theta_1) = 1.9304$, $k_2 = 2\cos(5\theta_2) = -1.9964$).
  - `2D-DEM-068`, `043`, `044` (`ABAbaBABab`): **Moth I** (Šuvakov & Dmitrašinović 2013; `I.B^{i.c.}_1`).
  - `2D-DEM-031` ($m_3=2$, `ABABabab`): **Butterfly I** (Šuvakov & Dmitrašinović 2013; LL2018 `I.A^{i.c.}_1 (2)`).
  - `2D-DEM-022` ($m_3=3$, featured as "`AABABABabaB` Heavy-Center"), `047` ($m_3=2$), `065` ($m_3=1.75$): True word is **`ABABABababab` = `(AB)^3(ab)^3`**, i.e., **Li, Jing & Liao (2018) `I.A^{i.c.}_2 (2)`**!
  - `2D-DEM-015` ($m_3=0.55$, featured as "20-Letter Unequal-Mass Braid"): **`ABaBAbABaBabAbaBabAb`** = exact mass continuation of **Li, Jing & Liao (2018) `I.A^{i.c.}_5 (0.5)`**!
  - `2D-DEM-018, 020, 023, 062` (`Braid-1.10`, `ABAbaBAbaBABabABab`): **Li & Liao (2017) `I.B^{i.c.}_5`** / **LL2018 `I.B^{i.c.}_6 (0.5)`**.
  - `2D-DEM-009` (`Heavy-Braid-1.75`, `AAbaBAbaBBab`): **Li, Jing & Liao (2018) `II.C^{i.c.}_5 (2)`**.
- **Only 10 of 48** have words not found verbatim in Li & Liao (2017, 2018):
  - 1 is the **5× Figure-8 duplicate** (`2D-DEM-064`).
  - 2 are **genuinely unstable** (`2D-DEM-037` = hyperbolic 3× `(AB)^3(ab)^3` with $k_1 = 2.000803$; `2D-DEM-049` = hyperbolic `(AB)^4(ab)^4` with $k_1 = 2.000002$).
  - 1 **fails the catalog's own `stab_margin > 1e-4` filter** (`2D-DEM-035`, true margin $4.8\times 10^{-5}$).
  - The remaining **6** genuinely stable orbits (`2D-DEM-016` $m_3=0.60$, `025` $m_3=0.65$, `026` $m_3=1.50$, `045` $m_3=3.00$, `054` $m_3=0.55$, `061` $m_3=4.00$ = $(AB)^5(ab)^5$) represent modest unequal-mass additions mostly at $m_3$ values not gridded in LL2018 (which scanned only $m_3 \in \{0.5, 0.75, 2, 4, 5, 8, 10\}$).

---

## 5. Critical Numerical, Topological, and Methodological Bugs

### 5.1 Finite-Difference Monodromy Corruption: 5 Confirmed Unstable, 4 Ultra-Near-Parabolic, and 4 Sub-Threshold Orbits

`src/campaign_2d.cpp:200-218` computes the $12\times 12$ monodromy matrix by central finite differences with fixed $\epsilon = 10^{-7}$ in IEEE-754 `double` and extracts $(k_1, k_2)$ solely from $\operatorname{tr}(M)$ and $\operatorname{tr}(M^2)$, without verifying symplecticity or eigenvalues. On `collinear_L0` orbits with close binary encounters ($r_{\min} = 0.008\text{–}0.10$, where gravitational Hessians reach $\sim 10^6$ and third derivatives $\sim 10^8$), $\epsilon = 10^{-7}$ enters the nonlinear regime and systematically biases $\operatorname{tr}(M)$ downward by $0.02\text{–}0.22$.

Re-refining the half-period spatial-involution symmetry $P(\theta)\,y(T/2) = y(0)$ (where $P^2 = I$ is symplectic, so $M(T) = (P M(T/2))^2$ and $k_j = h_j^2 - 2$, reducing the monodromy condition number from $\sim 10^4$ to $\sim 10\text{–}100$) in **80-bit `long double`** (`tol = 1e-17L`, exact analytical variational equations, canonical symplecticity error $3.8\times 10^{-11}\text{–}7.6\times 10^{-9}$; full-period refinement for asymmetric `2D-DEM-037`) and tracking iteration-to-iteration spread across converged Newton steps exposes **13 `collinear_L0` orbits violating the catalog's `linearly stable, stab_margin > 1e-4` filter**:

| ID | $m_3$ | $r_{\min}$ | Catalog $(k_1, k_2)$ & Margin | Exact Variational $(k_1, k_2)$ (`sympH`) | True Verdict & 2,000-Period Kick Test |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **`2D-DEM-037`** | 2.00 | 0.0705 | `(+1.9557, +1.8106)` · `0.0443` | **`(+2.0008033, +1.7727046)`** (`1.0e-9`) | **CONFIRMED UNSTABLE** ($\lambda = 1.02875$); drift grows $223\times$ by $500T$ (`0.0006 → 0.137 → 0.519`) |
| **`2D-DEM-057`** | 0.85 | 0.1018 | `(+1.9921, +1.5959)` · `0.0079` | **`(+2.0003437, +1.5884058)`** (`1.4e-10`) | **CONFIRMED UNSTABLE** ($\lambda = 1.01872$); **destroyed at period $546T$** ($574T$ for $v_y$ kick) |
| **`2D-DEM-069`** | 0.75 | 0.0207 | `(+1.9974, -1.2865)` · `0.0026` | **`(+2.0001553, -1.2867890)`** (`3.8e-11`) | **CONFIRMED UNSTABLE** ($\lambda = 1.01254$); **destroyed at period $811T$** ($783T$ for $v_y$ kick) |
| **`2D-DEM-049`** | 2.50 | 0.0867 | `(+1.9801, +1.3681)` · `0.0199` | **`(+2.0000020, +1.3595712)`** (`1.1e-10`) | **CONFIRMED UNSTABLE** ($\lambda = 1.00141$; $k_1-2 = +1.984\times 10^{-6} \pm 1.3\times 10^{-8}$) |
| **`2D-DEM-040`** | 0.90 | 0.0373 | `(+1.9666, -1.5454)` · `0.0334` | **`(+2.0000002, -1.5491107)`** (`5.6e-10`) | **CONFIRMED UNSTABLE** ($k_1-2 = +2.40\times 10^{-7} \pm 2.5\times 10^{-8}$, consistently $> 2$ across all 5 iterations) |
| **`2D-DEM-030`** | 1.05 | 0.0486 | `(+1.9287, +1.6222)` · `0.0713` | `(+2.0000000 ± 2.8e-7, +1.6132872)` (`1.3e-9`) | **ULTRA-NEAR-PARABOLIC** ($|k_1-2| \lesssim 3\times 10^{-7}$; fails `> 1e-4` by $>300\times$, overstated $>10^5\times$) |
| **`2D-DEM-034`** | 0.55 | 0.0606 | `(+1.9471, -0.4677)` · `0.0529` | `(+2.0000001 ± 4.6e-7, -0.4793621)` (`3.6e-10`) | **ULTRA-NEAR-PARABOLIC** ($|k_1-2| \lesssim 5\times 10^{-7}$; fails `> 1e-4` by $>200\times$, overstated $>10^5\times$) |
| **`2D-DEM-036`** | 0.85 | 0.0351 | `(+1.9510, +1.4192)` · `0.0490` | `(+1.9999998 ± 9.8e-7, +1.2208875)` (`7.6e-9`) | **ULTRA-NEAR-PARABOLIC** ($|k_1-2| \lesssim 10^{-6}$; fails `> 1e-4` by $>100\times$, catalog $\Delta k_2 = 0.198$!) |
| **`2D-DEM-042`** | 0.85 | 0.0287 | `(+1.9268, -1.9697)` · `0.0303` | `(+1.9999999 ± 3.8e-7, -1.9913952)` (`1.2e-9`) | **ULTRA-NEAR-PARABOLIC** ($|k_1-2| \lesssim 4\times 10^{-7}$; fails `> 1e-4` by $>250\times$, $\Delta k_1 = 0.073$) |
| `2D-DEM-039` | 1.50 | 0.0355 | `(+1.9575, -0.0610)` · `0.0425` | `(+1.9999871 ± 6.8e-7, -0.0509839)` (`3.3e-9`) | **Fails `> 1e-4` filter** (true margin `1.29e-5`, overstated **3,295×**!) |
| `2D-DEM-033` | 1.35 | 0.0696 | `(+1.9451, +0.3453)` · `0.0549` | `(+1.9999567 ± 1.1e-7, +0.3001128)` (`4.4e-10`) | **Fails `> 1e-4` filter** (true margin `4.33e-5`, overstated **1,269×**!) |
| `2D-DEM-035` | 1.10 | 0.0195 | `(+1.9498, +0.5286)` · `0.0502` | `(+1.9999524 ± 1.3e-7, +0.5204483)` (`9.3e-10`) | **Fails `> 1e-4` filter** (true margin `4.76e-5`, overstated **1,054×**!) |
| `2D-DEM-021` | 2.00 | 0.0793 | `(+1.7767, +1.5804)` · `0.2233` | `(+1.9999023 ± 8.1e-7, +1.3882430)` (`4.9e-9`) | **Fails `> 1e-4` filter** (true margin `9.77e-5`, overstated **2,286×**, $\Delta k_2 = 0.192$!) |

*(On the 4 ultra-near-parabolic orbits `030, 034, 036, 042`, $k_1 - 2$ fluctuates across zero at the $\pm(3\times 10^{-7}\text{–}10^{-6})$ close-encounter sensitivity floor of 80-bit DOPRI5, so determining the exact sign of $k_1 - 2$ requires quad-precision / arbitrary-precision Taylor integration; in either case, all four violate the catalog's `stab_margin > 1e-4` cutoff by more than two orders of magnitude. As a benchmark, the exact same half-period verifier recovers the Figure-8 `2D-DEM-067` to `retH = 4.3e-17`, `sympH = 9.7e-16`, `k1 - 2 = -2.800e-3 ± 1.2e-15`, `k2 = -0.5951933`.)*

*(Why did all 5 confirmed unstable orbits pass the repo's 100-period kick test? Because with $\lambda \in [1.0005, 1.0288]$, $(1.02875)^{100} \times 10^{-5} \approx 1.7\times 10^{-4} \ll 5.0$ — 100 periods is far too short to detect weak hyperbolic instabilities!)*

### 5.2 Catastrophic Free-Group Classifier Failure (`scripts/classify_topology.py`): 44 / 74 Words (59.5%) Are Wrong

Three separate bugs corrupt the reported Montgomery shape-sphere topology words and the claimed *"32 distinct free-group homotopy classes"*:

1. **Fixed-step RK4 blowup on close encounters (`classify_topology.py:108, 137-140`):**
   At $r_{\min} = 0.0078$ (`2D-DEM-017`), the encounter timescale is $r_{\min}^{3/2} \approx 6.9\times 10^{-4}$, comparable to `dt/32 = 2.6e-4`. Fixed-step RK4 loses the orbit and hallucinates short or trivial words:
   - `2D-DEM-022` (featured in `NEW_FAMILIES.md`): reported `AABABABabaB` (odd length 11!) $\to$ true word **`ABABABababab`** (length 12, `(AB)^3(ab)^3` = LL2018 `I.A^{i.c.}_2 (2)`).
   - `2D-DEM-016`: reported `1` (trivial!) $\to$ true word **`ABABaBABabABabaBababAbabABabABAb`** (length 32).
   - `2D-DEM-017`: reported `AAAbAb` (len 6) $\to$ true word **`ABABabABAbabABababABabaBABab`** (length 28 = LL2018 `I.A^{i.c.}_10 (0.75)`).
   - `2D-DEM-018` (`Braid-1.10`): reported `ABABab` (len 6) $\to$ true word **`ABAbaBAbaBABabABab`** (length 18, identical to `2D-DEM-020` and `023`).
   - `2D-DEM-043, 044`: reported `AA` and `1` $\to$ true word **`ABAbaBABab`** (Moth I, length 10).
   - In total, **41 of 48** `collinear_L0` words computed by `classify_topology.py` are wrong!

2. **Start-on-equator bug drops $t=0$ crossing (`classify_topology.py:128, 142`):**
   Every catalog orbit starts collinear ($y_i = 0 \implies n_z(0) = 0$). Checking `if prev[2] * cur[2] < 0.0:` skips $t = 0$ because `prev[2] == 0.0`. Whenever Body 3 is NOT between Bodies 1 and 2 at $t = 0$, the initial point lies on the `a`/`A` or `b`/`B` branch cut ($|u(0)| > 1/\sqrt{3}$), dropping a generator letter:
   - **`2D-DEM-006`** (featured in `README.md` & `NEW_FAMILIES.md` as `AABaaB`, $u(0) = +2.27 > 0.577$): `classify_topology.py` reports `AABaaB`, whose first-homology abelianization (taking $a, b = +1$ and $A, B = -1$) is $(0, -2)$ ($\gcd = 2$). The true exact word in the $r_{12}=\infty$ basis is `aBaBabb` (canonicalized to `AABaBaB`, abelianization $(0, -3)$, $\gcd = 3$). Because `aBaBabb` $= (a b^{-1})^3 b^3$, placing the third collision puncture $r_{23} = \infty$ (setting $c = a^{-1} b$ so $a b^{-1} = c^{-1}$) transforms it under $S_3$ relabeling directly to **`AAAbbb`** (abelianization $(-3, +3)$, $\gcd = 3$). Since the $\gcd$ of the abelianization is an $\mathrm{Aut}(F_2)$-invariant, **`AABaaB` is topologically impossible for this orbit!**
   - **`2D-DEM-048`** ($r_{\min} = 1.3907$, no close encounter): reported `A` instead of true `AA`.

3. **Missing $S_3$ body-relabeling quotient (`classify_topology.py:80-101`):**
   `canonical_word` mods out only `a ↔ b`, inversion, and cyclic shift (fixing puncture $r_{12} = \infty$), omitting the full $S_3$ outer automorphism group of the 3-punctured sphere. This artificially separates **`(Ab)^7`** (`2D-DEM-013`) from **`A^7`** (`2D-DEM-001`) and **`AbAb`** (`2D-DEM-071, 072`) from **`AA`** (`2D-DEM-056..060`).

### 5.3 Deduplication, Multiple-Cover, Continuation, and Vertical-Stability Issues

1. **Multiple covers evade deduplication (`src/campaign_2d.cpp:449`, `scripts/build_catalog.py:27-66`):**
   - `campaign_2d.cpp --collinear` only checks divisors `{2, 3, 4}` (missing 5× Figure-8 `2D-DEM-064`), and `--perp` checks **no** divisors (admitting 7 exact 2× covers in `2D-HIE`).
   - `build_catalog.py` deduplicates using $E T^{2/3}$, which scales by $k^{2/3}$ under $T \to kT$, so multiple covers never match their primitive orbit.
2. **Non-contiguous continuation windows & branch jumps (`scripts/build_catalog.py:97-112`, `data/continuations/`):**
   - `braid110_dn.dat` turns **unstable** for $m_3 \in [0.990, 1.040]$ and **jumps across 3 different orbit branches** ($T$ jumps $28.49 \to 30.15 \to 28.60 \to 25.69 \to 23.82 \to 28.72$, $k_1$ reaches `50.5`), yet `CATALOG.md` reports `st[0]["m3"]..st[-1]["m3"]` as a continuous *"Linearly stable window: $m_3 \in [0.9850, 1.1350]$"*.
   - `heavy300_up.dat` has an unstable point at $m_3 = 3.08$ ($k_2 = -2.00023$) hidden inside the reported $[2.84, 3.42]$ window.
3. **Degenerate vertical stability for $L = 0$ planar orbits (`src/vertical_bifurcation.cpp:83-96`):**
   For planar $L = 0$ orbits, the infinitesimal rotations $(R_x, R_y)$ form a Lagrangian subspace ($\{L_x, L_y\} = L_z = 0$) of the reduced 4D vertical phase space, forcing **all 6 eigenvalues of $M_z$ to $+1$ identically** ($\operatorname{tr}(M_z) - 4 \equiv 2$). Marking `stable_3d = 0` for all $L = 0$ orbits (including the Figure-8) reflects this symplectic degeneracy rather than exponential vertical instability.

---

## 6. Actionable Recommendations

To elevate this repository to publication-grade rigor:

1. **Replace finite-difference monodromy with analytical variational equations** (see `review/audit_orbits.cpp` and half-period involution reduction in `review/verify_collinear_ultra.cpp`), enforce canonical symplecticity $\|M_{\text{can}}^T J M_{\text{can}} - J\|_\infty < 10^{-7}$, compute the full eigenvalue spectrum via QR, and remove the **5 confirmed unstable**, **4 ultra-near-parabolic**, and **4 sub-threshold** `collinear_L0` orbits.
2. **Fix `scripts/classify_topology.py`** to use adaptive Dormand–Prince integration, start the crossing scan at $\arg\max_t |n_z(t)|$ (never on the equator $n_z = 0$), and canonicalize over all three $S_3$ collision projections and primitive roots.
3. **Add primitive-period and rotating-frame RPO reduction** before cataloging: strip all 8 exact inertial multiple covers (`2D-DEM-064` + 7 `2D-HIE`), group $q$-fold BHH inertial closures under their parent BHH relative periodic family (citing Broucke 1975, Hadjidemetriou 1975, Hénon 1976, and Janković & Dmitrašinović 2016), and replace the `swaps >= 2` democracy heuristic with a true bounded-pair / hierarchy check (requiring all 3 pairs to undergo genuine exchanges or testing $\max_t r_{ij} / \min_t r_{ij}$).
4. **Reframe `3D-Δ14` accurately** as the **5:2 vertical subharmonic bifurcation of the planar BHH retrograde family** at $m = [1, 1, 0.5]$ (recomputing exact $T_{\text{inert}}$ Floquet indices, $E, L$, and swap counts for `3D-002, 003, 005`, fixing the tangent bifurcation description at $v_{z2} \approx 0.678$, and fixing the scrambled IDs in `NEW_FAMILIES.md` and `scripts/generate_gifs.py`).
5. **Focus genuine novelty claims** on:
   - the `3D-Δ14` 5:2 BHH vertical subharmonic branch (and tracing it to equal masses / mapping neighboring $n:m$ vertical subharmonic branches),
   - the mass-stabilized **`AAAbbb` shuttle exchange family** (`2D-DEM-006`, `2D-DEM-024`), and
   - the **6 verified stable intermediate-$m_3$ $L=0$ orbits** not in Li & Liao (2017, 2018), while citing Li & Liao (2017, 2018) for the 38 matched $L=0$ families.
6. **Fix continuation window reporting** in `scripts/build_catalog.py` to detect branch jumps (continuity of $T, v_x, v_y$) and report contiguous stable sub-intervals rather than `[min_stable, max_stable]`.

---

## 7. References

- Broucke, R. (1975). *On relative periodic solutions of the planar general three-body problem.* Celestial Mechanics, 12, 439–462.
- Hadjidemetriou, J. D. (1975). *The continuation of periodic orbits from the restricted to the general three-body problem.* Celestial Mechanics, 12, 155–174.
- Hénon, M. (1976). *A family of periodic solutions of the planar three-body problem, and their stability.* Celestial Mechanics, 13, 267–285.
- Hénon, M. (1977). *Stability of interplay motions.* Celestial Mechanics, 15, 243–261.
- Ichtiaroglou, S., & Michalodimitrakis, M. (1980). *Three-body problem: the existence of families of three-dimensional periodic orbits which bifurcate from planar periodic orbits.* Astronomy & Astrophysics, 81, 30–32.
- Moore, C. (1993). *Braids in classical dynamics.* Phys. Rev. Lett., 70, 3675; Chenciner, A., & Montgomery, R. (2000). *A remarkable periodic solution of the three-body problem in the case of equal masses.* Annals of Mathematics, 152, 881–901.
- Šuvakov, M., & Dmitrašinović, V. (2013). *Three classes of Newtonian three-body planar periodic orbits.* Phys. Rev. Lett., 110, 114301.
- Janković, M. R., & Dmitrašinović, V. (2016). *Angular momentum and topological dependence of Kepler's third law in the Broucke–Hadjidemetriou–Hénon family of periodic three-body orbits.* Phys. Rev. Lett., 116, 064301.
- Li, X., & Liao, S. (2017). *More than six hundred new families of Newtonian periodic planar collisionless three-body orbits.* Sci. China Phys. Mech. Astron., 60, 129511.
- Li, X., Jing, Y., & Liao, S. (2018). *Over a thousand new periodic orbits of a planar three-body system with unequal masses.* Publ. Astron. Soc. Japan, 70, 64.
- Li, X., Li, X., & Liao, S. (2021). *One family of 13315 stable periodic orbits of non-hierarchical unequal-mass triple systems.* Sci. China Phys. Mech. Astron., 64, 219511.
- Li, X., & Liao, S. (2025). *Discovery of 10,059 new three-dimensional periodic orbits of general three-body problem.* (`sjtu-liao/three-body`).

---

## 8. Reproducing This Audit

Run the self-contained verification suite from the repository root (requires only `g++` with C++20 and `python3` stdlib):

```bash
bash review/reproduce.sh
```
