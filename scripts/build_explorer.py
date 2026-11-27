#!/usr/bin/env python3
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<title>3-Body Steady Orbit Discovery Explorer</title>
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
    --cyan: #06b6d4;
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
  .layout { display: grid; grid-template-columns: 310px 1fr 290px; flex: 1; min-height: 0; }
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
  .overlay-tl { position: absolute; top: 12px; left: 14px; pointer-events: none; display: flex; gap: 6px; align-items: center; }
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
    <h1>3-Body Steady Orbit Discovery Explorer</h1>
    <span id="hdrCounts" class="pill green"></span>
  </div>
  <div class="muted mono" style="font-size:11px;">Adaptive Dormand-Prince 5(4) · Exact Floquet Monodromy · 100T Verified</div>
</header>

<div class="layout">
  <aside class="sidebar">
    <div class="filters">
      <div class="tabs" id="filterTabs">
        <button class="tab active" data-f="all">All</button>
        <button class="tab" data-f="3d">3D Non-Planar</button>
        <button class="tab" data-f="dem">2D Democratic</button>
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
      <button id="btnKick" title="Inject a momentum-preserving +0.002 velocity perturbation">Kick (+2e-3)</button>
      <button id="btnReset">Reset</button>
      <button id="btnCopy">Copy y₀ JSON</button>
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
        <span id="badgeKick" class="badge-kick mono">PERTURBED</span>
      </div>
      <div class="overlay-bl mono">
        <span><span class="dot" style="background:var(--accent)"></span>m₁=<span id="lblM1"></span></span>
        <span><span class="dot" style="background:var(--amber)"></span>m₂=<span id="lblM2"></span></span>
        <span><span class="dot" style="background:var(--pink)"></span>m₃=<span id="lblM3"></span></span>
      </div>
      <div class="overlay-br mono" id="timeReadout">t = 0.00 (0.00 T)</div>
    </div>
  </main>

  <aside class="rightbar">
    <div class="card-box">
      <div class="card-title">Orbit Telemetry</div>
      <div class="stat-grid mono">
        <span class="muted">Subspace</span><span id="stSub"></span>
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
const ORBITS = CATALOG.orbits;

document.getElementById("hdrCounts").textContent =
  `${CATALOG.total_unique} Steady Orbits · ${CATALOG.counts["3d"]} 3D · ${CATALOG.counts["2d_democratic"]} 2D Democratic`;

// State & Filters
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

// Adaptive DOPRI5 for 2D (12) and 3D (18)
function deriv(s, m, dim, out) {
  const half = dim === 3 ? 9 : 6;
  for (let i = 0; i < half; i++) out[i] = s[half + i];
  for (let i = half; i < 2 * half; i++) out[i] = 0;
  for (let i = 0; i < 3; i++) {
    for (let j = i + 1; j < 3; j++) {
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

function advance(y, dt, m, dim, tol = 1e-10) {
  const n = dim === 3 ? 18 : 12;
  const k1 = new Float64Array(n), k2 = new Float64Array(n), k3 = new Float64Array(n),
        k4 = new Float64Array(n), k5 = new Float64Array(n), k6 = new Float64Array(n),
        k7 = new Float64Array(n), yt = new Float64Array(n), yn = new Float64Array(n);
  let t = 0, h = Math.min(0.005, dt), steps = 0;
  while (t < dt - 1e-14 && steps < 4000) {
    if (t + h > dt) h = dt - t;
    deriv(y, m, dim, k1);
    for (let i=0;i<n;i++) yt[i]=y[i]+h*0.2*k1[i]; deriv(yt, m, dim, k2);
    for (let i=0;i<n;i++) yt[i]=y[i]+h*(3/40*k1[i]+9/40*k2[i]); deriv(yt, m, dim, k3);
    for (let i=0;i<n;i++) yt[i]=y[i]+h*(44/45*k1[i]-56/15*k2[i]+32/9*k3[i]); deriv(yt, m, dim, k4);
    for (let i=0;i<n;i++) yt[i]=y[i]+h*(19372/6561*k1[i]-25360/2187*k2[i]+64448/6561*k3[i]-212/729*k4[i]); deriv(yt, m, dim, k5);
    for (let i=0;i<n;i++) yt[i]=y[i]+h*(9017/3168*k1[i]-355/33*k2[i]+46732/5247*k3[i]+49/176*k4[i]-5103/18656*k5[i]); deriv(yt, m, dim, k6);
    for (let i=0;i<n;i++) yn[i]=y[i]+h*(35/384*k1[i]+500/1113*k3[i]+125/192*k4[i]-2187/6784*k5[i]+11/84*k6[i]); deriv(yn, m, dim, k7);
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
  const d = currentOrbit.dim;
  let maxR = 0.5;
  for (const pt of refTrack) {
    for (let b = 0; b < 3; b++) {
      const r = Math.hypot(pt[d*b], pt[d*b+1], d===3 ? pt[d*b+2] : 0);
      if (r > maxR) maxR = r;
    }
  }
  bounds.cx = 0; bounds.cy = 0; bounds.cz = 0;
  bounds.scale = Math.min(simCanvas.width, simCanvas.height) * 0.41 / maxR;
}

function selectOrbit(o) {
  currentOrbit = o;
  state = new Float64Array(o.y0);
  simTime = 0;
  trail = [];
  document.getElementById("badgeKick").style.display = "none";

  // Precompute reference period track
  refTrack = [];
  const tmp = new Float64Array(o.y0);
  const N = 360, dt = o.T / N, d = o.dim;
  for (let i = 0; i <= N; i++) {
    refTrack.push(Array.from(tmp.slice(0, 3*d)));
    if (i < N) advance(tmp, dt, o.masses, d, 1e-10);
  }
  computeBounds();

  // Telemetry
  document.getElementById("badgeTitle").textContent = o.id;
  document.getElementById("badgeDim").textContent = o.dim === 3
    ? `3D Non-Planar (h=${(o.thickness_3d||0).toFixed(2)})`
    : (o.democratic ? `2D Democratic (${o.partner_swaps} swaps)` : "2D Hierarchical");
  document.getElementById("lblM1").textContent = o.masses[0];
  document.getElementById("lblM2").textContent = o.masses[1];
  document.getElementById("lblM3").textContent = o.masses[2];
  document.getElementById("stSub").textContent = o.subspace;
  document.getElementById("stT").textContent = o.T.toFixed(5);
  document.getElementById("stE").textContent = o.E.toFixed(5);
  document.getElementById("stL").textContent = o.L.toFixed(4);
  document.getElementById("stRmin").textContent = o.rmin.toFixed(4);
  document.getElementById("stSwaps").textContent = o.partner_swaps;
  document.getElementById("stRet").textContent = o.ret_err.toExponential(1);
  document.getElementById("stDrift").textContent = o.drift_100T.toFixed(5);
  document.getElementById("stMarg").textContent = o.stab_margin.toFixed(4);
  document.getElementById("stK").textContent = o.dim === 3
    ? `${o.k1.toFixed(2)}, ${o.k2.toFixed(2)}, ${o.k3.toFixed(2)}`
    : `${o.k1.toFixed(3)}, ${o.k2.toFixed(3)}`;
  document.getElementById("stKz").textContent = o.dim === 3
    ? "Full 3D Stable"
    : (o.stable_3d ? `YES (${o.kz.toFixed(3)})` : (o.kz !== undefined ? `L=0 (${o.kz.toFixed(2)})` : "—"));

  renderList();
  drawFloquet();
}

function applyFilter() {
  filtered = ORBITS.filter(o => {
    if (activeFilter === "3d") return o.dim === 3;
    if (activeFilter === "dem") return o.dim === 2 && o.democratic;
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
    div.innerHTML = `
      <div class="row">
        <span class="mono" style="font-weight:600;">${o.id}</span>
        <span class="mono" style="color:var(--green);font-size:11px;">marg ${o.stab_margin.toFixed(3)}</span>
      </div>
      <div class="row muted mono" style="font-size:11px;">
        <span>m=[${o.masses.join(",")}]</span>
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

  // Stability box [-2,2]^2
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
    const px = mx(o.k1), py = my(o.k2);
    const isCur = o.id === currentOrbit.id;
    fctx.beginPath();
    fctx.arc(px, py, isCur ? 5 : 2.8, 0, Math.PI*2);
    fctx.fillStyle = isCur ? "#3b82f6" : (o.dim===3 ? "#06b6d4" : (o.democratic ? "#22c55e" : "#8b9bb4"));
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
    const d = Math.hypot(mx(o.k1)-x, my(o.k2)-y);
    if (d < bestD) { bestD = d; best = o; }
  }
  if (best && bestD < 18) selectOrbit(best);
});

const COLS = ["#3b82f6", "#f59e0b", "#ec4899"];

function drawSim() {
  const w = simCanvas.width, h = simCanvas.height;
  ctx.clearRect(0, 0, w, h);
  const d = currentOrbit.dim;

  // Origin axes
  const o0 = project(0,0,0);
  ctx.strokeStyle = "rgba(139,155,180,0.2)";
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(o0[0]-8, o0[1]); ctx.lineTo(o0[0]+8, o0[1]);
  ctx.moveTo(o0[0], o0[1]-8); ctx.lineTo(o0[0], o0[1]+8);
  ctx.stroke();

  if (document.getElementById("chkGhost").checked) {
    for (let b = 0; b < 3; b++) {
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
    for (let b = 0; b < 3; b++) {
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

  // Draw bodies depth-sorted
  const bodies = [0,1,2].map(b => {
    const p = project(state[d*b], state[d*b+1], d===3 ? state[d*b+2] : 0);
    return { b, x: p[0], y: p[1], z: p[2], r: 5 + 3.2 * Math.cbrt(currentOrbit.masses[b]) };
  }).sort((a,b) => a.z - b.z);

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
    const dt = elapsed * spd * (currentOrbit.T / 8.0);
    advance(state, dt, currentOrbit.masses, currentOrbit.dim, 1e-10);
    simTime += dt;
    const d = currentOrbit.dim;
    trail.push(Array.from(state.slice(0, 3*d)));
    if (trail.length > 220) trail.shift();
    drawSim();
  }
  requestAnimationFrame(loop);
}

// Mouse drag for 3D rotation
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
  const m = currentOrbit.masses, d = currentOrbit.dim, half = d===3 ? 9 : 6, k = 0.002;
  state[half] += k; state[half + d] -= k * m[0] / m[1];
  if (d === 3) { state[half+2] += k; state[half+d+2] -= k * m[0] / m[1]; }
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
    print(f"Wrote {out_path} ({len(html)//1024} KB) with {len(cat['orbits'])} orbits.")

if __name__ == "__main__":
    main()
