#!/usr/bin/env python3
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<title>3-Body &amp; 4-Body Steady Orbit Discovery Explorer</title>
<style>
  :root {
    --bg: #0b0e14;
    --panel: #121722;
    --card: #181f2e;
    --border: #263045;
    --fg: #e6edf7;
    --muted: #8b9bb4;
    --accent: #3b82f6;
    --green: #22c55e;
    --amber: #f59e0b;
    --pink: #ec4899;
    --emerald: #10b981;
    --cyan: #06b6d4;
    --purple: #a855f7;
  }
  * { box-sizing: border-box; }
  body {
    margin: 0; background: var(--bg); color: var(--fg);
    font-family: ui-sans-serif, system-ui, -apple-system, Segoe UI, Roboto, sans-serif;
    font-size: 13px; height: 100vh; display: flex; flex-direction: column; overflow: hidden;
  }
  header {
    display: flex; align-items: center; justify-content: space-between;
    padding: 10px 16px; background: var(--panel); border-bottom: 1px solid var(--border);
  }
  header h1 { margin: 0; font-size: 15px; font-weight: 600; letter-spacing: 0.01em; }
  .pill { display: inline-block; padding: 2px 8px; border-radius: 999px; font-size: 11px; font-weight: 600; background: var(--card); border: 1px solid var(--border); color: var(--muted); }
  .pill.green { color: var(--green); border-color: rgba(34,197,94,0.4); }
  .pill.cyan { color: var(--cyan); border-color: rgba(6,182,212,0.4); }
  .pill.purple { color: var(--purple); border-color: rgba(168,85,247,0.4); }
  .layout { display: grid; grid-template-columns: 325px 1fr 295px; flex: 1; min-height: 0; }
  .sidebar, .rightbar { background: var(--panel); border-right: 1px solid var(--border); display: flex; flex-direction: column; min-height: 0; }
  .rightbar { border-right: none; border-left: 1px solid var(--border); padding: 12px; gap: 12px; overflow-y: auto; }
  .filters { padding: 10px; border-bottom: 1px solid var(--border); display: flex; flex-direction: column; gap: 8px; }
  .tabs { display: flex; flex-wrap: wrap; gap: 4px; }
  .tab {
    padding: 4px 8px; border-radius: 5px; border: 1px solid var(--border);
    background: var(--card); color: var(--muted); cursor: pointer; font-size: 11px;
  }
  .tab.active { background: var(--accent); color: #fff; border-color: var(--accent); }
  select, button {
    background: var(--card); color: var(--fg); border: 1px solid var(--border);
    border-radius: 6px; padding: 5px 10px; font-size: 12px; cursor: pointer;
  }
  button:hover { border-color: var(--accent); }
  .orbit-list { flex: 1; overflow-y: auto; padding: 6px; display: flex; flex-direction: column; gap: 4px; }
  .orbit-item {
    padding: 8px 10px; border-radius: 6px; border: 1px solid transparent;
    background: var(--card); cursor: pointer; display: flex; flex-direction: column; gap: 3px;
  }
  .orbit-item:hover { border-color: var(--border); }
  .orbit-item.active { border-color: var(--accent); background: #1d273b; }
  .orbit-item .row { display: flex; justify-content: space-between; align-items: center; }
  .mono { font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-variant-numeric: tabular-nums; }
  .muted { color: var(--muted); }
  .main { display: flex; flex-direction: column; min-height: 0; position: relative; }
  .toolbar {
    display: flex; flex-wrap: wrap; align-items: center; gap: 8px;
    padding: 8px 14px; background: var(--panel); border-bottom: 1px solid var(--border);
  }
  .canvas-wrap { flex: 1; position: relative; min-height: 0; }
  #simCanvas { width: 100%; height: 100%; display: block; cursor: grab; }
  #simCanvas:active { cursor: grabbing; }
  .overlay-tl { position: absolute; top: 12px; left: 14px; pointer-events: none; display: flex; gap: 6px; align-items: center; flex-wrap: wrap; }
  .overlay-bl { position: absolute; bottom: 12px; left: 14px; pointer-events: none; display: flex; gap: 14px; font-size: 12px; }
  .overlay-br { position: absolute; bottom: 12px; right: 14px; pointer-events: none; font-size: 12px; color: var(--muted); }
  .dot { width: 9px; height: 9px; border-radius: 50%; display: inline-block; margin-right: 5px; }
  .stat-grid { display: grid; grid-template-columns: 1fr auto; gap: 5px 8px; font-size: 12px; }
  .card-box { background: var(--card); border: 1px solid var(--border); border-radius: 8px; padding: 10px; }
  .card-title { font-size: 11px; font-weight: 600; text-transform: uppercase; letter-spacing: 0.05em; color: var(--muted); margin-bottom: 8px; }
  #floquetCanvas { width: 100%; height: 220px; display: block; border-radius: 6px; background: #0e131d; cursor: crosshair; }
  .badge-kick { background: var(--amber); color: #000; padding: 2px 7px; border-radius: 4px; font-weight: 700; font-size: 11px; display: none; }
</style>
</head>
<body>
<header>
  <div style="display:flex;align-items:center;gap:10px;">
    <h1>3-Body &amp; 4-Body Steady Orbit Discovery Explorer</h1>
    <span id="hdrCounts" class="pill green"></span>
  </div>
  <div class="muted mono" style="font-size:11px;">DOPRI5 · Exact Floquet / Francis QR · Newtonian + 1PN GR + 4-Body</div>
</header>

<div class="layout">
  <aside class="sidebar">
    <div class="filters">
      <div class="tabs" id="filterTabs">
        <button class="tab active" data-f="all">All</button>
        <button class="tab" data-f="3d">3D Non-Planar</button>
        <button class="tab" data-f="dem">2D Democratic</button>
        <button class="tab" data-f="brake">Brake (v₀=0)</button>
        <button class="tab" data-f="choreo">Rot. Choreo</button>
        <button class="tab" data-f="1pn">1PN GR</button>
        <button class="tab" data-f="4b">4-Body</button>
        <button class="tab" data-f="3dok">3D-Stable</button>
        <button class="tab" data-f="col">Collinear L=0</button>
        <button class="tab" data-f="perp">Perp L≠0</button>
      </div>
      <div style="display:flex;align-items:center;gap:6px;">
        <span class="muted" style="font-size:11px;">Sort:</span>
        <select id="sortSelect" style="flex:1;">
          <option value="margin">Stability Margin (deepest)</option>
          <option value="swaps">Partner Swaps (most democratic)</option>
          <option value="rmin">Min Separation rₘᵢₙ (widest)</option>
          <option value="T">Period T (shortest)</option>
        </select>
      </div>
    </div>
    <div class="orbit-list" id="orbitList"></div>
  </aside>

  <main class="main">
    <div class="toolbar">
      <button id="btnPlay">Pause</button>
      <button id="btnKick" title="Inject a momentum-preserving +0.002 perturbation">Kick (+2e-3)</button>
      <button id="btnReset">Reset</button>
      <button id="btnCopy">Copy JSON</button>
      <label class="muted" style="display:flex;align-items:center;gap:6px;margin-left:auto;">
        Speed <input id="speedRange" type="range" min="0.2" max="5" step="0.1" value="1.2" />
      </label>
      <label class="muted" style="display:flex;align-items:center;gap:4px;">
        <input id="chkGhost" type="checkbox" checked /> Ghost Track
      </label>
    </div>
    <div class="canvas-wrap">
      <canvas id="simCanvas"></canvas>
      <div class="overlay-tl">
        <span id="badgeTitle" class="pill cyan mono"></span>
        <span id="badgeDim" class="pill"></span>
        <span id="badgePhys" class="pill purple mono" style="display:none;"></span>
        <span id="badgeKick" class="badge-kick mono">PERTURBED</span>
      </div>
      <div class="overlay-bl mono">
        <span><span class="dot" style="background:var(--accent)"></span>m₁=<span id="lblM1"></span></span>
        <span><span class="dot" style="background:var(--amber)"></span>m₂=<span id="lblM2"></span></span>
        <span><span class="dot" style="background:var(--pink)"></span>m₃=<span id="lblM3"></span></span>
        <span id="wrapM4" style="display:none;"><span class="dot" style="background:var(--emerald)"></span>m₄=<span id="lblM4"></span></span>
      </div>
      <div class="overlay-br mono" id="timeReadout">t = 0.00 (0.00 T)</div>
    </div>
  </main>

  <aside class="rightbar">
    <div class="card-box">
      <div class="card-title">Orbit Telemetry</div>
      <div class="stat-grid mono">
        <span class="muted">Subspace</span><span id="stSub"></span>
        <span class="muted">Physics</span><span id="stPhys"></span>
        <span class="muted">Period T</span><span id="stT"></span>
        <span class="muted">Energy E</span><span id="stE"></span>
        <span class="muted">Ang. Mom |L|</span><span id="stL"></span>
        <span class="muted">Min Sep rₘᵢₙ</span><span id="stRmin"></span>
        <span class="muted">Partner Swaps</span><span id="stSwaps"></span>
        <span class="muted">Return Error</span><span id="stRet"></span>
        <span class="muted">100T Drift</span><span id="stDrift" style="color:var(--green)"></span>
        <span class="muted">Stability Margin</span><span id="stMarg" style="color:var(--green);font-weight:700"></span>
        <span class="muted">Floquet k</span><span id="stK"></span>
        <span class="muted">3D Vertical k_z</span><span id="stKz"></span>
      </div>
    </div>

    <div class="card-box">
      <div class="card-title">Floquet Stability Square (k₁, k₂) ∈ [-2, 2]²</div>
      <canvas id="floquetCanvas" width="250" height="220"></canvas>
      <div class="muted" style="font-size:11px;margin-top:6px;">
        Click any dot to load that orbit. Closer to center (0,0) = deeper linear stability island.
        Drag canvas to rotate 3D orbits.
      </div>
    </div>
  </aside>
</div>

<script>
const CATALOG = __CATALOG_JSON__;
const ORBITS = [
  ...CATALOG.orbits.map(o => ({ ...o, _kind: "newt3", n_bodies: 3 })),
  ...(CATALOG.relativistic_1pn || []).map(o => ({ ...o, _kind: "1pn", n_bodies: 3 })),
  ...(CATALOG.four_body || []).map(o => ({ ...o, _kind: "4b", n_bodies: 4, k1: (o.k_indices||[0,0])[0], k2: (o.k_indices||[0,0])[1], stable_3d: o.dim === 3 })),
];

document.getElementById("hdrCounts").textContent =
  `${ORBITS.length} Verified Steady Orbits · ${CATALOG.counts["3d"]} 3D (3B) · ${CATALOG.counts["relativistic_1pn"]||0} 1PN GR · ${CATALOG.counts["four_body"]||0} 4-Body`;

let activeFilter = "all";
let activeSort = "margin";
let filtered = [...ORBITS];
let currentOrbit = ORBITS[0];
let state = null;
let simTime = 0;
let playing = true;
let trail = [];
let refTrack = [];
let bounds = { cx: 0, cy: 0, cz: 0, scale: 100 };
let camYaw = 0.65, camPitch = 0.45, dragging = false, lastMX = 0, lastMY = 0;

function derivNewt(s, m, nb, dim, out) {
  const half = nb * dim;
  for (let i = 0; i < half; i++) out[i] = s[half + i];
  for (let i = half; i < 2 * half; i++) out[i] = 0;
  for (let i = 0; i < nb; i++) {
    for (let j = i + 1; j < nb; j++) {
      const dx = s[dim*j] - s[dim*i];
      const dy = s[dim*j+1] - s[dim*i+1];
      const dz = dim === 3 ? (s[dim*j+2] - s[dim*i+2]) : 0;
      const r2 = dx*dx + dy*dy + dz*dz;
      const ir3 = 1 / (r2 * Math.sqrt(r2));
      out[half + dim*i]   += m[j] * dx * ir3;
      out[half + dim*i+1] += m[j] * dy * ir3;
      out[half + dim*j]   -= m[i] * dx * ir3;
      out[half + dim*j+1] -= m[i] * dy * ir3;
      if (dim === 3) {
        out[half + dim*i+2] += m[j] * dz * ir3;
        out[half + dim*j+2] -= m[i] * dz * ir3;
      }
    }
  }
}

function deriv1PN(s, m, dim, beta, out) {
  const half = 3 * dim;
  const rx = [s[0], s[dim], s[2*dim]];
  const ry = [s[1], s[dim+1], s[2*dim+1]];
  const rz = dim === 3 ? [s[2], s[5], s[8]] : [0, 0, 0];
  const px = [s[half], s[half+dim], s[half+2*dim]];
  const py = [s[half+1], s[half+dim+1], s[half+2*dim+1]];
  const pz = dim === 3 ? [s[half+2], s[half+5], s[half+8]] : [0, 0, 0];
  const drx = [0,0,0], dry = [0,0,0], drz = [0,0,0];
  const dpx = [0,0,0], dpy = [0,0,0], dpz = [0,0,0];
  const p2 = [0,0,0], phi = [0,0,0];
  for (let a = 0; a < 3; a++) {
    p2[a] = px[a]*px[a] + py[a]*py[a] + pz[a]*pz[a];
    const f0 = 1/m[a] - beta*0.5*p2[a]/(m[a]*m[a]*m[a]);
    drx[a] = f0*px[a]; dry[a] = f0*py[a]; drz[a] = f0*pz[a];
  }
  const pairs = [[0,1],[0,2],[1,2]];
  const dxA = [0,0,0], dyA = [0,0,0], dzA = [0,0,0], irA = [0,0,0], ir2A = [0,0,0], ir3A = [0,0,0];
  for (let q = 0; q < 3; q++) {
    const [a,b] = pairs[q];
    const dx = rx[a]-rx[b], dy = ry[a]-ry[b], dz = rz[a]-rz[b];
    const ir = 1 / Math.sqrt(dx*dx + dy*dy + dz*dz);
    const ir2 = ir*ir, ir3 = ir2*ir;
    dxA[q]=dx; dyA[q]=dy; dzA[q]=dz; irA[q]=ir; ir2A[q]=ir2; ir3A[q]=ir3;
    phi[a] += m[b]*ir; phi[b] += m[a]*ir;
    const fac = -m[a]*m[b]*ir3;
    dpx[a] += fac*dx; dpy[a] += fac*dy; dpz[a] += fac*dz;
    dpx[b] -= fac*dx; dpy[b] -= fac*dy; dpz[b] -= fac*dz;
  }
  for (let q = 0; q < 3; q++) {
    const [a,b] = pairs[q];
    const ir = irA[q], ir2 = ir2A[q], ir3 = ir3A[q];
    const nx = dxA[q]*ir, ny = dyA[q]*ir, nz = dzA[q]*ir;
    const papb = px[a]*px[b] + py[a]*py[b] + pz[a]*pz[b];
    const npa = nx*px[a] + ny*py[a] + nz*pz[a];
    const npb = nx*px[b] + ny*py[b] + nz*pz[b];
    const Kab = 3*((m[b]/m[a])*p2[a] + (m[a]/m[b])*p2[b]) - 7*papb - npa*npb;
    const c = -0.5*beta*ir;
    drx[a] += c*(6*(m[b]/m[a])*px[a] - 7*px[b] - npb*nx);
    dry[a] += c*(6*(m[b]/m[a])*py[a] - 7*py[b] - npb*ny);
    drz[a] += c*(6*(m[b]/m[a])*pz[a] - 7*pz[b] - npb*nz);
    drx[b] += c*(6*(m[a]/m[b])*px[b] - 7*px[a] - npa*nx);
    dry[b] += c*(6*(m[a]/m[b])*py[b] - 7*py[a] - npa*ny);
    drz[b] += c*(6*(m[a]/m[b])*pz[b] - 7*pz[a] - npa*nz);
    const potF = -m[a]*m[b]*(phi[a]+phi[b])*ir3;
    const dKx = -ir*((px[a]-npa*nx)*npb + (px[b]-npb*nx)*npa);
    const dKy = -ir*((py[a]-npa*ny)*npb + (py[b]-npb*ny)*npa);
    const dKz = -ir*((pz[a]-npa*nz)*npb + (pz[b]-npb*nz)*npa);
    const tx = beta*(0.5*Kab*ir2*nx - 0.5*ir*dKx + potF*dxA[q]);
    const ty = beta*(0.5*Kab*ir2*ny - 0.5*ir*dKy + potF*dyA[q]);
    const tz = beta*(0.5*Kab*ir2*nz - 0.5*ir*dKz + potF*dzA[q]);
    dpx[a] -= tx; dpy[a] -= ty; dpz[a] -= tz;
    dpx[b] += tx; dpy[b] += ty; dpz[b] += tz;
  }
  for (let a = 0; a < 3; a++) {
    out[dim*a] = drx[a]; out[dim*a+1] = dry[a];
    out[half+dim*a] = dpx[a]; out[half+dim*a+1] = dpy[a];
    if (dim === 3) { out[dim*a+2] = drz[a]; out[half+dim*a+2] = dpz[a]; }
  }
}

function advance(y, dt, o, tol = 1e-10) {
  const nb = o.n_bodies || o.masses.length;
  const dim = o.dim;
  const beta = o.beta_inv_c2 || 0;
  const m = o.masses;
  const n = 2 * nb * dim;
  const k1 = new Float64Array(n), k2 = new Float64Array(n), k3 = new Float64Array(n),
        k4 = new Float64Array(n), k5 = new Float64Array(n), k6 = new Float64Array(n),
        k7 = new Float64Array(n), yt = new Float64Array(n), yn = new Float64Array(n);
  const evalFn = (s, out) => beta > 0 ? deriv1PN(s, m, dim, beta, out) : derivNewt(s, m, nb, dim, out);
  let t = 0, h = Math.min(0.005, dt), steps = 0;
  while (t < dt - 1e-14 && steps < 6000) {
    if (t + h > dt) h = dt - t;
    evalFn(y, k1);
    for (let i=0;i<n;i++) yt[i]=y[i]+h*0.2*k1[i]; evalFn(yt, k2);
    for (let i=0;i<n;i++) yt[i]=y[i]+h*(3/40*k1[i]+9/40*k2[i]); evalFn(yt, k3);
    for (let i=0;i<n;i++) yt[i]=y[i]+h*(44/45*k1[i]-56/15*k2[i]+32/9*k3[i]); evalFn(yt, k4);
    for (let i=0;i<n;i++) yt[i]=y[i]+h*(19372/6561*k1[i]-25360/2187*k2[i]+64448/6561*k3[i]-212/729*k4[i]); evalFn(yt, k5);
    for (let i=0;i<n;i++) yt[i]=y[i]+h*(9017/3168*k1[i]-355/33*k2[i]+46732/5247*k3[i]+49/176*k4[i]-5103/18656*k5[i]); evalFn(yt, k6);
    for (let i=0;i<n;i++) yn[i]=y[i]+h*(35/384*k1[i]+500/1113*k3[i]+125/192*k4[i]-2187/6784*k5[i]+11/84*k6[i]); evalFn(yn, k7);
    let err = 0;
    for (let i=0;i<n;i++) {
      const ei = h*(71/57600*k1[i]-71/16695*k3[i]+71/1920*k4[i]-17253/339200*k5[i]+22/525*k6[i]-1/40*k7[i]);
      const sc = tol*(1 + Math.max(Math.abs(y[i]), Math.abs(yn[i])));
      const r = ei/sc; err += r*r;
    }
    err = Math.sqrt(err/n);
    if (err <= 1) {
      y.set(yn); t += h; steps++;
      h = Math.min(h * (err===0 ? 5 : Math.min(5, Math.max(0.2, 0.9*Math.pow(err,-0.2)))), 0.06);
    } else {
      h *= Math.max(0.1, 0.9*Math.pow(err,-0.25));
      if (h < 1e-12) break;
    }
  }
}

const simCanvas = document.getElementById("simCanvas");
const ctx = simCanvas.getContext("2d");
const fCanvas = document.getElementById("floquetCanvas");
const fctx = fCanvas.getContext("2d");

function resize() {
  const rect = simCanvas.parentElement.getBoundingClientRect();
  simCanvas.width = Math.max(300, Math.floor(rect.width));
  simCanvas.height = Math.max(300, Math.floor(rect.height));
  computeBounds();
}
window.addEventListener("resize", resize);

function project(x, y, z) {
  const w = simCanvas.width, h = simCanvas.height;
  if (currentOrbit.dim === 2) {
    return [w*0.5 + (x - bounds.cx)*bounds.scale, h*0.5 - (y - bounds.cy)*bounds.scale, 0];
  }
  const dx = x - bounds.cx, dy = y - bounds.cy, dz = z - bounds.cz;
  const cy = Math.cos(camYaw), sy = Math.sin(camYaw);
  const cp = Math.cos(camPitch), sp = Math.sin(camPitch);
  const rx = dx*cy - dy*sy;
  const ry = dx*sy + dy*cy;
  const px = rx;
  const py = ry*sp + dz*cp;
  const pz = ry*cp - dz*sp;
  return [w*0.5 + px*bounds.scale, h*0.5 - py*bounds.scale, pz];
}

function computeBounds() {
  if (!refTrack.length) return;
  const d = currentOrbit.dim, nb = currentOrbit.n_bodies;
  let maxR = 0.5;
  for (const pt of refTrack) {
    for (let b = 0; b < nb; b++) {
      const r = Math.hypot(pt[d*b], pt[d*b+1], d===3 ? pt[d*b+2] : 0);
      if (r > maxR) maxR = r;
    }
  }
  bounds.cx = 0; bounds.cy = 0; bounds.cz = 0;
  bounds.scale = Math.min(simCanvas.width, simCanvas.height) * 0.41 / maxR;
}

function initialStateOf(o) {
  return o.z0_canonical ? o.z0_canonical : o.y0;
}

function selectOrbit(o) {
  currentOrbit = o;
  const initVec = initialStateOf(o);
  state = new Float64Array(initVec);
  simTime = 0;
  trail = [];
  document.getElementById("badgeKick").style.display = "none";

  // Precompute reference period track (cap at 120 time units for ultra-long q-fold closures if needed)
  refTrack = [];
  const tmp = new Float64Array(initVec);
  const T_ref = Math.min(o.T, 380);
  const N = T_ref > 100 ? 720 : 360;
  const dt = T_ref / N, d = o.dim, nb = o.n_bodies;
  for (let i = 0; i <= N; i++) {
    refTrack.push(Array.from(tmp.slice(0, nb*d)));
    if (i < N) advance(tmp, dt, o, 1e-10);
  }
  computeBounds();

  document.getElementById("badgeTitle").textContent = o.id;
  document.getElementById("badgeDim").textContent = o.dim === 3
    ? `${nb}B 3D (h=${(o.thickness_3d||0).toFixed(2)})`
    : (o.democratic ? `${nb}B 2D Democratic (${o.partner_swaps} sw)` : `${nb}B 2D Hierarchical`);
  const bp = document.getElementById("badgePhys");
  if (o._kind === "1pn") {
    bp.style.display = "inline-block";
    bp.textContent = `1PN GR β=${o.beta_inv_c2.toFixed(5)} (v²/c²=${(o.beta_eff_v2c2*100).toFixed(2)}%)`;
  } else if (o._kind === "4b") {
    bp.style.display = "inline-block";
    bp.textContent = `4-Body (${o.dim === 3 ? "24D" : "16D"} QR)`;
  } else {
    bp.style.display = "none";
  }

  document.getElementById("lblM1").textContent = o.masses[0];
  document.getElementById("lblM2").textContent = o.masses[1];
  document.getElementById("lblM3").textContent = o.masses[2];
  const w4 = document.getElementById("wrapM4");
  if (nb === 4) { w4.style.display = "inline"; document.getElementById("lblM4").textContent = o.masses[3]; }
  else { w4.style.display = "none"; }

  document.getElementById("stSub").textContent = o.subspace;
  document.getElementById("stPhys").textContent = o._kind === "1pn" ? `1PN ADM (β=${o.beta_inv_c2.toFixed(4)})` : (nb === 4 ? "4-Body Newton" : "3-Body Newton");
  document.getElementById("stT").textContent = o.T.toFixed(4);
  document.getElementById("stE").textContent = o.E.toFixed(5);
  document.getElementById("stL").textContent = o.L.toFixed(4);
  document.getElementById("stRmin").textContent = o.rmin.toFixed(4);
  document.getElementById("stSwaps").textContent = o.partner_swaps;
  document.getElementById("stRet").textContent = o.ret_err.toExponential(1);
  document.getElementById("stDrift").textContent = o.drift_100T.toFixed(5);
  document.getElementById("stMarg").textContent = o.stab_margin.toFixed(5);
  if (o._kind === "4b") {
    document.getElementById("stK").textContent = (o.k_indices || []).slice(0, 3).map(v => v.toFixed(3)).join(", ") + "…";
    document.getElementById("stKz").textContent = o.dim === 3 ? `24D Stable (Δ|λ|=${o.max_eig_rad_err.toExponential(1)})` : `16D Stable`;
  } else {
    document.getElementById("stK").textContent = o.dim === 3
      ? `${o.k1.toFixed(2)}, ${o.k2.toFixed(2)}, ${o.k3.toFixed(2)}`
      : `${o.k1.toFixed(3)}, ${o.k2.toFixed(3)}`;
    document.getElementById("stKz").textContent = o.dim === 3
      ? "Full 3D Stable"
      : (o.stable_3d ? `YES (${o.kz.toFixed(3)})` : (o.kz !== undefined ? `L=0 (${o.kz.toFixed(2)})` : "—"));
  }

  renderList();
  drawFloquet();
}

function applyFilter() {
  filtered = ORBITS.filter(o => {
    if (activeFilter === "3d") return o.dim === 3 && o._kind === "newt3";
    if (activeFilter === "dem") return o.dim === 2 && o.democratic && o._kind === "newt3";
    if (activeFilter === "brake") return o.subspace === "free_fall_brake";
    if (activeFilter === "choreo") return o.subspace.includes("choreo");
    if (activeFilter === "1pn") return o._kind === "1pn";
    if (activeFilter === "4b") return o._kind === "4b";
    if (activeFilter === "3dok") return o.stable_3d;
    if (activeFilter === "col") return o.subspace === "collinear_L0";
    if (activeFilter === "perp") return o.subspace === "perp_L_nonzero";
    return true;
  });
  filtered.sort((a, b) => {
    if (activeSort === "margin") return b.stab_margin - a.stab_margin;
    if (activeSort === "swaps") return b.partner_swaps - a.partner_swaps;
    if (activeSort === "rmin") return b.rmin - a.rmin;
    if (activeSort === "T") return a.T - b.T;
    return 0;
  });
  renderList();
  drawFloquet();
}

function renderList() {
  const el = document.getElementById("orbitList");
  el.innerHTML = "";
  for (const o of filtered) {
    const div = document.createElement("div");
    div.className = "orbit-item" + (o.id === currentOrbit.id ? " active" : "");
    const mStr = o.masses.map(x => Number(x.toFixed(3))).join(",");
    div.innerHTML = `
      <div class="row">
        <span class="mono" style="font-weight:600;">${o.id} <span class="muted" style="font-weight:400;font-size:10px;">${o.subspace}</span></span>
        <span class="mono" style="color:var(--green);font-size:11px;">${o.stab_margin.toFixed(3)}</span>
      </div>
      <div class="row muted mono" style="font-size:11px;">
        <span>m=[${mStr}]</span>
        <span>T=${o.T.toFixed(2)} · sw=${o.partner_swaps}</span>
      </div>`;
    div.onclick = () => selectOrbit(o);
    el.appendChild(div);
  }
}

function drawFloquet() {
  const w = fCanvas.width, h = fCanvas.height;
  fctx.clearRect(0, 0, w, h);
  const pad = 16;
  const mx = k => pad + ((k + 2.2)/4.4)*(w - 2*pad);
  const my = k => h - pad - ((k + 2.2)/4.4)*(h - 2*pad);

  fctx.fillStyle = "rgba(34,197,94,0.08)";
  fctx.fillRect(mx(-2), my(2), mx(2)-mx(-2), my(-2)-my(2));
  fctx.strokeStyle = "rgba(34,197,94,0.5)";
  fctx.setLineDash([4,4]);
  fctx.strokeRect(mx(-2), my(2), mx(2)-mx(-2), my(-2)-my(2));
  fctx.setLineDash([]);

  fctx.strokeStyle = "rgba(139,155,180,0.2)";
  fctx.beginPath();
  fctx.moveTo(mx(-2.2), my(0)); fctx.lineTo(mx(2.2), my(0));
  fctx.moveTo(mx(0), my(-2.2)); fctx.lineTo(mx(0), my(2.2));
  fctx.stroke();

  for (const o of filtered) {
    if (o.k1 === undefined || o.k2 === undefined) continue;
    const px = mx(o.k1), py = my(o.k2);
    const isCur = o.id === currentOrbit.id;
    fctx.beginPath();
    fctx.arc(px, py, isCur ? 5 : 2.8, 0, Math.PI*2);
    fctx.fillStyle = isCur ? "#3b82f6" : (o._kind==="1pn" ? "#a855f7" : (o._kind==="4b" ? "#10b981" : (o.dim===3 ? "#06b6d4" : (o.democratic ? "#22c55e" : "#8b9bb4"))));
    fctx.fill();
    if (isCur) { fctx.strokeStyle = "#fff"; fctx.lineWidth = 1.5; fctx.stroke(); }
  }
}

fCanvas.addEventListener("click", (ev) => {
  const rect = fCanvas.getBoundingClientRect();
  const x = (ev.clientX - rect.left) * (fCanvas.width / rect.width);
  const y = (ev.clientY - rect.top) * (fCanvas.height / rect.height);
  const pad = 16, w = fCanvas.width, h = fCanvas.height;
  const mx = k => pad + ((k + 2.2)/4.4)*(w - 2*pad);
  const my = k => h - pad - ((k + 2.2)/4.4)*(h - 2*pad);
  let best = null, bestD = 1e9;
  for (const o of filtered) {
    if (o.k1 === undefined || o.k2 === undefined) continue;
    const d = Math.hypot(mx(o.k1)-x, my(o.k2)-y);
    if (d < bestD) { bestD = d; best = o; }
  }
  if (best && bestD < 18) selectOrbit(best);
});

const COLS = ["#3b82f6", "#f59e0b", "#ec4899", "#10b981"];

function drawSim() {
  const w = simCanvas.width, h = simCanvas.height;
  ctx.clearRect(0, 0, w, h);
  const d = currentOrbit.dim, nb = currentOrbit.n_bodies;

  const o0 = project(0,0,0);
  ctx.strokeStyle = "rgba(139,155,180,0.2)";
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(o0[0]-8, o0[1]); ctx.lineTo(o0[0]+8, o0[1]);
  ctx.moveTo(o0[0], o0[1]-8); ctx.lineTo(o0[0], o0[1]+8);
  ctx.stroke();

  if (document.getElementById("chkGhost").checked) {
    for (let b = 0; b < nb; b++) {
      ctx.beginPath();
      for (let i = 0; i < refTrack.length; i++) {
        const pt = refTrack[i];
        const p = project(pt[d*b], pt[d*b+1], d===3 ? pt[d*b+2] : 0);
        if (i === 0) ctx.moveTo(p[0], p[1]); else ctx.lineTo(p[0], p[1]);
      }
      ctx.strokeStyle = COLS[b];
      ctx.globalAlpha = 0.22;
      ctx.lineWidth = 1.4;
      ctx.stroke();
    }
    ctx.globalAlpha = 1;
  }

  if (trail.length > 1) {
    for (let b = 0; b < nb; b++) {
      ctx.beginPath();
      for (let i = 0; i < trail.length; i++) {
        const pt = trail[i];
        const p = project(pt[d*b], pt[d*b+1], d===3 ? pt[d*b+2] : 0);
        if (i === 0) ctx.moveTo(p[0], p[1]); else ctx.lineTo(p[0], p[1]);
      }
      ctx.strokeStyle = COLS[b];
      ctx.globalAlpha = 0.85;
      ctx.lineWidth = 2.2;
      ctx.stroke();
    }
    ctx.globalAlpha = 1;
  }

  const bodies = [];
  for (let b = 0; b < nb; b++) {
    const p = project(state[d*b], state[d*b+1], d===3 ? state[d*b+2] : 0);
    bodies.push({ b, x: p[0], y: p[1], z: p[2], r: 4.8 + 3.2 * Math.cbrt(Math.max(0.02, currentOrbit.masses[b])) });
  }
  bodies.sort((a,b) => a.z - b.z);

  for (const bd of bodies) {
    ctx.beginPath();
    ctx.arc(bd.x, bd.y, bd.r, 0, Math.PI*2);
    ctx.fillStyle = COLS[bd.b];
    ctx.fill();
    ctx.strokeStyle = "#fff";
    ctx.lineWidth = 1.3;
    ctx.stroke();
  }

  document.getElementById("timeReadout").textContent =
    `t = ${simTime.toFixed(2)} (${(simTime / currentOrbit.T).toFixed(2)} T)`;
}

let lastTs = 0;
function loop(ts) {
  if (!lastTs) lastTs = ts;
  const elapsed = Math.min(0.05, (ts - lastTs) / 1000);
  lastTs = ts;
  if (playing && state) {
    const spd = parseFloat(document.getElementById("speedRange").value);
    const T_vis = Math.min(currentOrbit.T, currentOrbit.T_rpo ? currentOrbit.T_rpo * 3 : currentOrbit.T);
    const dt = elapsed * spd * (T_vis / 8.0);
    advance(state, dt, currentOrbit, 1e-10);
    simTime += dt;
    const d = currentOrbit.dim, nb = currentOrbit.n_bodies;
    trail.push(Array.from(state.slice(0, nb*d)));
    if (trail.length > 240) trail.shift();
    drawSim();
  }
  requestAnimationFrame(loop);
}

simCanvas.addEventListener("mousedown", e => { dragging = true; lastMX = e.clientX; lastMY = e.clientY; });
window.addEventListener("mouseup", () => dragging = false);
window.addEventListener("mousemove", e => {
  if (!dragging) return;
  camYaw += (e.clientX - lastMX) * 0.01;
  camPitch = Math.max(-1.5, Math.min(1.5, camPitch + (e.clientY - lastMY) * 0.01));
  lastMX = e.clientX; lastMY = e.clientY;
  drawSim();
});

document.querySelectorAll("#filterTabs .tab").forEach(btn => {
  btn.onclick = () => {
    document.querySelectorAll("#filterTabs .tab").forEach(b => b.classList.remove("active"));
    btn.classList.add("active");
    activeFilter = btn.dataset.f;
    applyFilter();
  };
});
document.getElementById("sortSelect").onchange = e => { activeSort = e.target.value; applyFilter(); };
document.getElementById("btnPlay").onclick = e => { playing = !playing; e.target.textContent = playing ? "Pause" : "Play"; };
document.getElementById("btnReset").onclick = () => selectOrbit(currentOrbit);
document.getElementById("btnKick").onclick = () => {
  const m = currentOrbit.masses, d = currentOrbit.dim, nb = currentOrbit.n_bodies;
  const half = nb * d, k = 0.002;
  if (currentOrbit._kind === "1pn") {
    state[half] += k; state[half + d] -= k;
    if (d === 3) { state[half+2] += k; state[half+d+2] -= k; }
  } else {
    state[half] += k; state[half + d] -= k * m[0] / m[1];
    if (d === 3) { state[half+2] += k; state[half+d+2] -= k * m[0] / m[1]; }
  }
  document.getElementById("badgeKick").style.display = "inline-block";
};
document.getElementById("btnCopy").onclick = () => {
  navigator.clipboard?.writeText(JSON.stringify(currentOrbit, null, 2));
};

resize();
applyFilter();
if (ORBITS.length) selectOrbit(ORBITS[0]);
requestAnimationFrame(loop);
</script>
</body>
</html>
"""

def main():
    cat_path = os.path.join(ROOT, "catalog.json")
    with open(cat_path) as f:
        cat = json.load(f)
    html = HTML_TEMPLATE.replace("__CATALOG_JSON__", json.dumps(cat))
    out_path = os.path.join(ROOT, "explorer.html")
    with open(out_path, "w") as f:
        f.write(html)
    total = len(cat["orbits"]) + len(cat.get("relativistic_1pn", [])) + len(cat.get("four_body", []))
    print(f"Wrote {out_path} ({len(html)//1024} KB) with {total} orbits.")

if __name__ == "__main__":
    main()
