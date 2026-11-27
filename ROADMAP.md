# Next-Phase Discovery Roadmap

Planned directions for discovering genuinely new families of steady (linearly stable) $N$-body periodic orbits in this repository.

---

## 1. Free-Fall (Brake) Orbits ($\mathbf{v}_i(0) = \mathbf{0}$)
- **Subspace:** All three bodies start at rest on a triangle in the Agekyan–Anosova domain: $\mathbf{r}_1 = (-0.5, 0)$, $\mathbf{r}_2 = (0.5, 0)$, $\mathbf{r}_3 = (x_3, y_3)$ shifted to center-of-mass $= 0$, with $\mathbf{v}_1 = \mathbf{v}_2 = \mathbf{v}_3 = \mathbf{0}$.
- **Periodicity condition:** Simultaneous zero-velocity stop at $t = T/2$ ($|\mathbf{v}(T/2)| = 0$), which by time-reversal symmetry $t \to -t$ guarantees exact closure at $T$.
- **Goal:** Scan the 2D $(x_3, y_3)$ domain across equal and unequal masses $(m_1, m_2, m_3)$ to find **linearly stable** ($|k_1|, |k_2| < 2$) brake orbits (almost all published brake orbits are unstable).
- **Planned tool:** `src/search_brake.cpp`

---

## 2. Period-Doubled 3D Democratic Families via Vertical Bifurcation ($k_z \to \pm 2$)
- **Starting seeds:** Existing high-margin, wide-separation 2D democratic orbits already sitting near vertical bifurcations:
  - `2D-DEM-001` (`m = [1, 1, 0.5]`, $r_{\min} = 0.982$, $k_z = -1.916$)
  - `2D-DEM-002` (`m = [1, 1, 0.6]`, $r_{\min} = 0.991$, $k_z = -1.680$)
  - `3D-005` (`m = [1, 0.7, 0.4]`, all-distinct-mass 3D RPO)
- **Method:** Continue 2D orbits in $(m_3, x_2)$ until $k_z = -2$, then branch into the period-doubled ($2T$) 3D RPO manifold in the co-precessing frame and map the full 2-parameter $(v_{z2}, m_3)$ stability sheet + rational inertial closures ($q$-fold).
- **Planned tool:** `src/branch_3d_vertical.cpp`

---

## 3. 3-Body Rotating-Frame Choreographies ($L \neq 0$)
- **Concept:** Unlike the $L = 0$ Figure-8, all three equal-mass bodies trace the **exact same closed curve** in a uniformly rotating frame ($\omega \neq 0$) at phase offsets $0, T/3, 2T/3$.
- **Method:** Solve $\mathbf{r}_1(T/3) = R(\omega T/3)\,\mathbf{r}_2(0)$ with variational/Newton shooting in both 2D (trefoil / 5-petal choreographies) and 3D (helical choreographies), then verify $18\times 18$ rotating-frame Floquet stability and snap rational $\omega T / (2\pi) = p/q$ inertial closures.
- **Planned tool:** `src/search_choreo_rot.cpp`

---

## 4. Post-Newtonian (1PN) Relativistic Steady Families
- **Physics:** Add the Einstein–Infeld–Hoffmann (1PN) acceleration terms parameterized by $\beta = v^2 / c^2$.
- **Goal:** Trace how Newtonian steady families (`3D-Δ14`, `2D-DEM-001`, Figure-8) deform with $\beta$, measure whether 1PN pericenter precession **widens** stability islands or **stabilizes** marginally unstable Newtonian orbits, and snap new relativistic rational closures.
- **Planned tool:** `src/search_1pn.cpp`

---

## 5. 4-Body 3D Steady Extensions (Hip-Hop & 3+1 Resonances)
- **Subspaces:**
  - **3D Hip-Hop symmetry:** Two equal-mass pairs alternating above and below the $xy$-plane.
  - **3+1 Hybrid:** Place a 4th light or heavy body at the center or outer resonance of our deepest 3-body stability islands (`2D-DEM-001`, `3D-001`).
- **Verification:** Full $24\times 24$ (3D) or $16\times 16$ (2D) monodromy matrix stability + 100-period empirical perturbation survival.
- **Planned tool:** `src/search_4body.cpp`
