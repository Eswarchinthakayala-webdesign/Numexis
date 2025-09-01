"use client";

/**
 * src/pages/Chapter30.jsx
 *
 * CHAPTER 30 — Finite Difference: Parabolic Equations
 * Sections:
 *  30.1 The Heat-Conduction Equation
 *  30.2 Explicit Methods (FTCS)
 *  30.3 A Simple Implicit Method (BTCS)
 *  30.4 The Crank–Nicolson Method
 *  30.5 Parabolic Equations in Two Spatial Dimensions
 *  Problems
 *
 * Design & Features:
 *  - Dark responsive layout (Tailwind utility classes)
 *  - React-Markdown for text & KaTeX for math (components prop)
 *  - Responsive Plotly charts wrapped to avoid overflow
 *  - Numerical solvers implemented in JavaScript:
 *      FTCS (explicit), BTCS (implicit via tridiagonal Thomas solver),
 *      Crank–Nicolson (CN), 2D explicit method
 *  - Interactive controls for Δx, Δt, α, ICs and boundary conditions
 *  - Residual/error & stability checks, comparative plots
 *
 * Notes:
 *  - Assumes UI components at "@/components/ui/*" or substitute with plain inputs/buttons
 *  - This file is long and split into two code blocks in the UI response; paste both into a file.
 */

/* =========================
   Imports & Theme + Helpers
   ========================= */
import React, { useEffect, useMemo, useState } from "react";
import Plot from "react-plotly.js";
import { motion } from "framer-motion";
import ReactMarkdown from "react-markdown";
import remarkMath from "remark-math";
import rehypeKatex from "rehype-katex";
import "katex/dist/katex.min.css";
import BottomBar from "../components/BottomBar";
import {
  Activity,
  Zap,
  CircleDot,
  ListChecks,
  BookOpen,
  Code,
  Grid,
  Layers,
} from "lucide-react";

import { Card, CardContent } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Button } from "@/components/ui/button";
import { Textarea } from "@/components/ui/textarea";

/* Theme constants */
const THEME = {
  BG: "bg-zinc-950",
  PANEL: "bg-zinc-900/60 border border-zinc-700",
  PRIMARY: "#22d3ee",
  ACCENT: "#34d399",
  WARN: "#f59e0b",
  DANGER: "#ef4444",
};

/* ReactMarkdown wrapper without className prop (v8+ compatibility) */
function MD({ children, className = "" }) {
  const components = {
    p: ({ node, ...props }) => <p className={`text-sm text-zinc-300 leading-relaxed ${className}`} {...props} />,
    strong: ({ node, ...props }) => <strong className="text-emerald-300" {...props} />,
    em: ({ node, ...props }) => <em className="text-zinc-200" {...props} />,
    code: ({ node, inline, ...props }) =>
      inline ? (
        <code className="bg-zinc-800/70 px-1 rounded text-xs text-cyan-200" {...props} />
      ) : (
        <pre className="bg-zinc-800/60 rounded p-2 text-xs overflow-auto" {...props} />
      ),
    li: ({ node, ...props }) => <li className="ml-4 text-sm text-zinc-300" {...props} />,
    h1: ({ node, ...props }) => <h1 className="text-xl text-zinc-100 font-semibold" {...props} />,
    h2: ({ node, ...props }) => <h2 className="text-lg text-zinc-100 font-semibold" {...props} />,
  };
  return (
    <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]} components={components}>
      {children}
    </ReactMarkdown>
  );
}

/* Small presentational header */
function SectionHeader({ icon: Icon, number, title, accent = THEME.PRIMARY }) {
  return (
    <div className="flex items-center gap-3">
      <div className="p-2 rounded-lg" style={{ background: "rgba(255,255,255,0.02)" }}>
        <Icon className="w-6 h-6" style={{ color: accent }} />
      </div>
      <div>
        <div className="text-xs text-zinc-400 uppercase tracking-wider">{number}</div>
        <div className="text-lg font-semibold text-zinc-100">{title}</div>
      </div>
    </div>
  );
}

/* Responsive plot wrapper to avoid overflow */
function ResponsivePlot({ children }) {
  return (
    <div className="w-full max-w-full overflow-hidden rounded-lg">
      <div style={{ minHeight: 120 }}>{children}</div>
    </div>
  );
}

/* =========================
   Numerical Utilities
   ========================= */

/* linspace */
function linspace(a, b, n) {
  if (n <= 1) return [a];
  const out = [];
  const dx = (b - a) / (n - 1);
  for (let i = 0; i < n; i++) out.push(a + i * dx);
  return out;
}

/* Create 1D vector filled with val */
function createVector(n, val = 0) {
  return new Array(n).fill(val);
}

/* Thomas algorithm (tridiagonal solver)
   a: sub-diagonal (length n-1)
   b: diagonal (length n)
   c: super-diagonal (length n-1)
   d: RHS (length n)
   returns solution x (length n)
*/
function thomasSolve(a, b, c, d) {
  const n = b.length;
  const cp = new Array(n - 1);
  const dp = new Array(n);
  const x = new Array(n);
  cp[0] = c[0] / b[0];
  dp[0] = d[0] / b[0];
  for (let i = 1; i < n - 1; i++) {
    const denom = b[i] - a[i - 1] * cp[i - 1];
    cp[i] = c[i] / denom;
    dp[i] = (d[i] - a[i - 1] * dp[i - 1]) / denom;
  }
  dp[n - 1] = (d[n - 1] - a[n - 2] * dp[n - 2]) / (b[n - 1] - a[n - 2] * cp[n - 2]);
  x[n - 1] = dp[n - 1];
  for (let i = n - 2; i >= 0; i--) {
    x[i] = dp[i] - cp[i] * x[i + 1];
  }
  return x;
}

/* Norm helpers */
function linftyNorm(vec) {
  return Math.max(...vec.map((v) => Math.abs(v)));
}

/* =========================
   PDE-specific helpers
   ========================= */

/* 1D grid builder for domain [0, L] with M interior intervals giving M+1 nodes */
function build1DGrid(M, L = 1) {
  const dx = L / M;
  const xs = linspace(0, L, M + 1);
  return { M, dx, xs };
}

/* initial condition examples for 1D problems */
function initialCondition1D(name, x, L = 1) {
  if (name === "gaussian") {
    const sigma = 0.08;
    const center = 0.5 * L;
    return Math.exp(-((x - center) * (x - center)) / (2 * sigma * sigma));
  }
  if (name === "sine") {
    return Math.sin(Math.PI * x);
  }
  if (name === "step") {
    return x < 0.5 ? 1.0 : 0.0;
  }
  return 0;
}

/* analytic 1D heat solution for comparison: for example initial sine mode u(x,0)=sin(pi x), Dirichlet 0 BCs
   u(x,t)=exp(-α π^2 t) sin(π x)
*/
function analyticHeatSine(x, t, alpha) {
  return Math.exp(-alpha * Math.PI * Math.PI * t) * Math.sin(Math.PI * x);
}

/* =========================
   FTCS: explicit solver for 1D heat equation
   u_t = α u_xx
   Forward time, central space:
   u_j^{n+1} = u_j^n + r (u_{j+1}^n - 2 u_j^n + u_{j-1}^n)
   r = α Δt / Δx^2  (stability requires r <= 1/2 for explicit scheme)
   ========================= */
function solveFTCS1D({ M = 100, L = 1, alpha = 0.01, dt = 1e-4, tmax = 0.01, ic = "gaussian", bc = { left: 0, right: 0 } }) {
  const { dx, xs } = build1DGrid(M, L);
  const steps = Math.max(1, Math.ceil(tmax / dt));
  const r = (alpha * dt) / (dx * dx);
  const u = new Array(M + 1);
  const Ustore = []; // store snapshots every few steps
  for (let j = 0; j <= M; j++) u[j] = initialCondition1D(ic, xs[j], L);

  // enforce Dirichlet BCs at boundaries (simple)
  u[0] = bc.left;
  u[M] = bc.right;
  Ustore.push({ t: 0, u: u.slice() });

  for (let n = 1; n <= steps; n++) {
    const uNew = u.slice();
    for (let j = 1; j < M; j++) {
      uNew[j] = u[j] + r * (u[j + 1] - 2 * u[j] + u[j - 1]);
    }
    // enforce BCs
    uNew[0] = bc.left;
    uNew[M] = bc.right;
    // swap
    for (let j = 0; j <= M; j++) u[j] = uNew[j];
    // store at reasonable intervals
    if (n % Math.max(1, Math.floor(steps / 40)) === 0) Ustore.push({ t: n * dt, u: u.slice() });
  }
  return { xs, Ustore, dx, dt, r };
}

/* =========================
   BTCS: backward Euler implicit method (1D)
   (I - r A) u^{n+1} = u^n
   produces tridiagonal linear system solved via Thomas algorithm
   ========================= */
function solveBTCS1D({ M = 100, L = 1, alpha = 0.01, dt = 1e-4, tmax = 0.01, ic = "gaussian", bc = { left: 0, right: 0 } }) {
  const { dx, xs } = build1DGrid(M, L);
  const steps = Math.max(1, Math.ceil(tmax / dt));
  const r = (alpha * dt) / (dx * dx);

  // build tridiagonal coefficients for (I - r A)
  // interior unknowns 1..M-1 (size N = M-1)
  const N = M - 1;
  const a = new Array(N - 1).fill(-r); // sub-diagonal
  const b = new Array(N).fill(1 + 2 * r); // diagonal
  const c = new Array(N - 1).fill(-r); // super-diagonal

  // initial u
  let u = new Array(M + 1);
  for (let j = 0; j <= M; j++) u[j] = initialCondition1D(ic, xs[j], L);
  u[0] = bc.left;
  u[M] = bc.right;
  const Ustore = [{ t: 0, u: u.slice() }];

  for (let n = 1; n <= steps; n++) {
    // build RHS vector d = u_interior^n + boundary contributions
    const d = new Array(N);
    for (let i = 1; i <= N; i++) {
      d[i - 1] = u[i];
    }
    // include boundary contributions into first and last RHS entries
    d[0] += r * bc.left;
    d[N - 1] += r * bc.right;

    // solve tridiagonal system
    const uInterior = thomasSolve(a, b, c, d);
    // assemble u^{n+1}
    const uNew = u.slice();
    for (let i = 1; i <= N; i++) uNew[i] = uInterior[i - 1];
    // enforce BCs
    uNew[0] = bc.left;
    uNew[M] = bc.right;
    u = uNew;
    if (n % Math.max(1, Math.floor(steps / 40)) === 0) Ustore.push({ t: n * dt, u: u.slice() });
  }
  return { xs, Ustore, dx, dt, r };
}

/* =========================
   Crank–Nicolson (CN): (I - 0.5 r A) u^{n+1} = (I + 0.5 r A) u^n
   also tridiagonal; stable and second-order in time
   ========================= */
function solveCrankNicolson1D({ M = 100, L = 1, alpha = 0.01, dt = 1e-4, tmax = 0.01, ic = "gaussian", bc = { left: 0, right: 0 } }) {
  const { dx, xs } = build1DGrid(M, L);
  const steps = Math.max(1, Math.ceil(tmax / dt));
  const r = (alpha * dt) / (dx * dx);

  const N = M - 1;
  // left side tridiagonal: (I - 0.5 r A)
  const aL = new Array(N - 1).fill(-0.5 * r);
  const bL = new Array(N).fill(1 + r);
  const cL = new Array(N - 1).fill(-0.5 * r);
  // right side operator applied to u^n: (I + 0.5 r A) -> used to build RHS each step
  const aR = new Array(N - 1).fill(0.5 * r);
  const bR = new Array(N).fill(1 - r);
  const cR = new Array(N - 1).fill(0.5 * r);

  // initial condition
  let u = new Array(M + 1);
  for (let j = 0; j <= M; j++) u[j] = initialCondition1D(ic, xs[j], L);
  u[0] = bc.left;
  u[M] = bc.right;
  const Ustore = [{ t: 0, u: u.slice() }];

  for (let n = 1; n <= steps; n++) {
    // build RHS d = (I + 0.5 r A) * u^n (only interior entries)
    const d = new Array(N);
    for (let i = 1; i <= N; i++) {
      const left = i - 1 >= 1 ? u[i - 1] : bc.left;
      const right = i + 1 <= M - 1 ? u[i + 1] : bc.right;
      d[i - 1] = (1 - r) * u[i] + 0.5 * r * (left + right);
    }
    // include boundary terms if necessary (already included using left/right above due to bc values)
    // solve (I - 0.5 r A) u^{n+1} = d
    const uInterior = thomasSolve(aL, bL, cL, d);
    const uNew = u.slice();
    for (let i = 1; i <= N; i++) uNew[i] = uInterior[i - 1];
    uNew[0] = bc.left;
    uNew[M] = bc.right;
    u = uNew;
    if (n % Math.max(1, Math.floor(steps / 40)) === 0) Ustore.push({ t: n * dt, u: u.slice() });
  }
  return { xs, Ustore, dx, dt, r };
}

/* =========================
   2D explicit FTCS (for demonstration)
   u_t = α (u_xx + u_yy)
   u_{i,j}^{n+1} = u_{i,j}^n + r (u_{i+1,j} + u_{i-1,j} + u_{i,j+1} + u_{i,j-1} - 4 u_{i,j})
   where r = α Δt / Δx^2 (for square grid Δx=Δy)
   Stability roughly: r <= 1/4
   ========================= */
function solveFTCS2D({ N = 50, L = 1, alpha = 0.1, dt = 1e-4, tmax = 0.01, ic = "gaussian2d", bc = { value: 0 } }) {
  // grid nodes: (N+1)x(N+1)
  const dx = L / N;
  const xs = linspace(0, L, N + 1);
  const steps = Math.max(1, Math.ceil(tmax / dt));
  const r = (alpha * dt) / (dx * dx);

  // initialize grid u[i][j] with i along x, j along y
  let u = new Array(N + 1);
  for (let i = 0; i <= N; i++) {
    u[i] = new Array(N + 1);
    for (let j = 0; j <= N; j++) {
      // default Gaussian centered
      const xc = 0.5 * L;
      const yc = 0.5 * L;
      const sigma = 0.08;
      const val = Math.exp(-(((xs[i] - xc) ** 2) + ((xs[j] - yc) ** 2)) / (2 * sigma * sigma));
      u[i][j] = val;
    }
  }
  // enforce Dirichlet BCs (constant value)
  for (let i = 0; i <= N; i++) {
    u[i][0] = bc.value;
    u[i][N] = bc.value;
    u[0][i] = bc.value;
    u[N][i] = bc.value;
  }

  const snapshots = [{ t: 0, u: u.map((row) => row.slice()) }];

  for (let n = 1; n <= steps; n++) {
    const uNew = new Array(N + 1);
    for (let i = 0; i <= N; i++) {
      uNew[i] = u[i].slice();
    }
    for (let i = 1; i < N; i++) {
      for (let j = 1; j < N; j++) {
        uNew[i][j] = u[i][j] + r * (u[i + 1][j] + u[i - 1][j] + u[i][j + 1] + u[i][j - 1] - 4 * u[i][j]);
      }
    }
    // apply BCs (Dirichlet)
    for (let i = 0; i <= N; i++) {
      uNew[i][0] = bc.value;
      uNew[i][N] = bc.value;
      uNew[0][i] = bc.value;
      uNew[N][i] = bc.value;
    }
    u = uNew;
    if (n % Math.max(1, Math.floor(steps / 20)) === 0) snapshots.push({ t: n * dt, u: u.map((row) => row.slice()) });
  }

  return { xs, snapshots, dx, dt, r };
}

/* =========================
   Plot helpers (convert 1D/2D data into Plotly arrays)
   ========================= */

/* 1D series snapshots -> { x, traces: [{x,y,t}...] } */
function plot1DSnapshots(Ustore) {
  const traces = Ustore.map((snap) => ({ x: snap.u.map((v, idx) => idx), y: snap.u, t: snap.t }));
  const x = Ustore[0].u.map((_, idx) => idx);
  return { x, traces };
}

/* 1D xs + snapshot with actual x coordinates */
function plot1DSnapshotsWithXs(xs, Ustore) {
  const traces = Ustore.map((snap) => ({ x: xs, y: snap.u, t: snap.t }));
  return traces;
}

/* 2D grid to Plotly-friendly */
function grid2Plotly(uGrid, xs) {
  // uGrid should be u[i][j] where i corresponds to x index and j to y index
  // Plotly surface expects Z as rows corresponding to y
  const N = xs.length - 1;
  const Z = new Array(N + 1);
  for (let j = 0; j <= N; j++) {
    Z[j] = new Array(N + 1);
    for (let i = 0; i <= N; i++) {
      Z[j][i] = uGrid[i][j];
    }
  }
  return { x: xs, y: xs, Z };
}

/* =========================
   UI Sections: 30.1 - 30.5 (part 1)
   ========================= */

/* 30.1 The Heat-Conduction Equation */
function Section301() {
  return (
    <motion.div className={`${THEME.PANEL} p-4 rounded-2xl`} initial={{ opacity: 0, y: 8 }} whileInView={{ opacity: 1, y: 0 }}>
      <div className="flex items-center justify-between">
        <SectionHeader icon={Activity} number="30.1" title="The Heat-Conduction Equation" />
        <div className="text-xs text-zinc-400">Parabolic PDE — time-dependent diffusion</div>
      </div>

      <div className="grid gap-4 md:grid-cols-3 mt-4">
        <div className="md:col-span-1 space-y-3">
          <MD>
  {`The heat equation in one dimension (constant diffusivity $\\alpha$):

$$\\frac{\\partial u}{\\partial t} = \\alpha \\frac{\\partial^2 u}{\\partial x^2}, \\quad x \\in (0,L),\\; t>0.$$

Initial condition: $u(x,0) = u_0(x)$. Boundary conditions: Dirichlet, Neumann, or mixed.`}
</MD>

<MD>
  {"Physical interpretation: diffusion of heat (or concentration). The PDE smooths initial data; high-frequency modes decay faster (see analytic sine-mode solution)."}
</MD>

<MD>
  {"We'll implement explicit and implicit finite-difference schemes and compare their stability and accuracy."}
</MD>

        </div>

        <div className="md:col-span-2 space-y-3">
          <div className="rounded-lg border border-zinc-700 p-3 bg-zinc-900/50">
           <MD>
  {`Analytic sine-mode solution (Dirichlet):

$u(x,t) = e^{-\\alpha \\pi^2 t} \\sin(\\pi x)$

This is a good manufactured solution for testing.`}
</MD>

          </div>
        </div>
      </div>
    </motion.div>
  );
}

/* 30.2 Explicit Methods (FTCS) UI and example runner */
function Section302() {
  const [M, setM] = useState(100);
  const [alpha, setAlpha] = useState(0.01);
  const [dt, setDt] = useState(0.00005);
  const [tmax, setTmax] = useState(0.01);
  const [ic, setIc] = useState("sine");
  const [run, setRun] = useState(false);
  const [result, setResult] = useState(null);

  useEffect(() => {
    if (!run) return;
    const res = solveFTCS1D({ M: Math.max(10, M), alpha: parseFloat(alpha), dt: parseFloat(dt), tmax: parseFloat(tmax), ic, bc: { left: 0, right: 0 } });
    setResult(res);
    setRun(false);
  }, [run]); // eslint-disable-line

  // compute stability indicator
  const stabilityIndicator = useMemo(() => {
    const { dx } = build1DGrid(Math.max(10, M));
    const r = (alpha * dt) / (dx * dx);
    return r;
  }, [M, alpha, dt]);

  return (
    <motion.div className={`${THEME.PANEL} p-4 rounded-2xl`} initial={{ opacity: 0, y: 8 }} whileInView={{ opacity: 1, y: 0 }}>
      <div className="flex items-center justify-between">
        <SectionHeader icon={Zap} number="30.2" title="Explicit Methods (FTCS)" accent={THEME.ACCENT} />
        <div className="text-xs text-zinc-400">CFL-like stability: r = α Δt / Δx² ≤ 1/2</div>
      </div>

      <div className="grid gap-4 md:grid-cols-4 mt-4">
        <div className="md:col-span-1 space-y-3">
          <MD>{`The Forward-Time Central-Space (FTCS) scheme is straightforward but conditionally stable. Choose Δt so that r ≤ 1/2 for stability (1D).`}</MD>
          <div className="grid gap-2">
            <label className="text-xs text-zinc-400">M (intervals)</label>
            <Input className="text-white" value={M} type="number" onChange={(e) => setM(Math.max(10, parseInt(e.target.value || 10, 10)))} />
            <label className="text-xs text-zinc-400">α (diffusivity)</label>
            <Input className="text-white" value={alpha} type="number" step="0.0001" onChange={(e) => setAlpha(parseFloat(e.target.value || 0.01))} />
            <label className="text-xs text-zinc-400">Δt</label>
            <Input className="text-white" value={dt} type="number" step="1e-6" onChange={(e) => setDt(parseFloat(e.target.value || 1e-5))} />
            <label className="text-xs text-zinc-400">tₘₐₓ</label>
            <Input className="text-white" value={tmax} type="number" step="0.001" onChange={(e) => setTmax(parseFloat(e.target.value || 0.01))} />
            <label className="text-xs text-zinc-400">Initial condition</label>
            <select value={ic} onChange={(e) => setIc(e.target.value)} className="bg-zinc-800/90 text-zinc-100 rounded border border-zinc-700 p-2 w-full">
              <option value="sine">sin(π x)</option>
              <option value="gaussian">Gaussian bump</option>
              <option value="step">Step</option>
            </select>
            <div className="flex gap-2 mt-2">
              <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={() => setRun(true)}>Run FTCS</Button>
              <Button className="bg-white cursor-pointer" variant="ghost" onClick={() => setResult(null)}>Clear</Button>
            </div>
          </div>

          <div className="text-xs text-zinc-400 mt-2">
            <div>Stability r = <span className="text-zinc-100 font-medium">{stabilityIndicator.toFixed(6)}</span></div>
            <div className="text-xs text-zinc-500 mt-1">{stabilityIndicator <= 0.5 ? "Likely stable (r ≤ 0.5)" : "Unstable region (r > 0.5)"}</div>
          </div>
        </div>

        <div className="md:col-span-3 space-y-3">
          <div className="rounded-lg border border-zinc-700 p-3 bg-zinc-900/50">
            {!result && <div className="text-sm text-zinc-400">Run a simulation to see snapshots of u(x,t).</div>}
            {result && (
              <>
                <div className="grid md:grid-cols-2 gap-3">
                  <div className="rounded-lg border border-zinc-700 p-2 bg-zinc-950/20">
                    <ResponsivePlot>
                      <Plot
                        data={result.Ustore.slice(0, 30).map((snap, idx) => ({
                          x: result.xs,
                          y: snap.u,
                          mode: "lines",
                          name: `t=${snap.t.toFixed(4)}`,
                          visible: idx === result.Ustore.length - 1 ? true : "legendonly",
                        }))}
                        layout={{
                          title: "FTCS: snapshots (select trace in legend)",
                          autosize: true,
                          paper_bgcolor: "rgba(0,0,0,0)",
                          plot_bgcolor: "rgba(0,0,0,0)",
                          margin: { t: 30, b: 40 },
                          legend: { orientation: "h", y: -0.2, xanchor: "center", x: 0.5 },
                        }}
                        style={{ width: "100%", height: 360 }}
                        useResizeHandler
                        config={{ responsive: true }}
                      />
                    </ResponsivePlot>
                  </div>

                  <div className="rounded-lg border border-zinc-700 p-2 bg-zinc-950/20">
                    <ResponsivePlot>
                      <Plot
                        data={[
                          {
                            z: result.Ustore[result.Ustore.length - 1].u,
                            x: result.xs,
                            y: [0], // placeholder – we'll plot surface as 2D line time series instead
                            type: "scatter",
                            mode: "lines",
                          },
                        ]}
                        layout={{
                          title: "Final profile (FTCS)",
                          autosize: true,
                          paper_bgcolor: "rgba(0,0,0,0)",
                          plot_bgcolor: "rgba(0,0,0,0)",
                        }}
                        style={{ width: "100%", height: 360 }}
                        useResizeHandler
                        config={{ responsive: true }}
                      />
                    </ResponsivePlot>
                  </div>
                </div>

                <div className="mt-3">
                  <h4 className="text-sm text-zinc-100 mb-2">Compare with analytic (sine mode)</h4>
                  <ResponsivePlot>
                    <Plot
                      data={[
                        {
                          x: result.xs,
                          y: result.Ustore[result.Ustore.length - 1].u,
                          mode: "lines",
                          name: "FTCS numeric",
                          line: { color: THEME.PRIMARY },
                        },
                        {
                          x: result.xs,
                          y: result.xs.map((x) => analyticHeatSine(x, result.Ustore[result.Ustore.length - 1].t, alpha)),
                          mode: "lines",
                          name: "analytic sine",
                          line: { dash: "dash", color: THEME.ACCENT },
                        },
                      ]}
                      layout={{
                        title: `Comparison at t=${result.Ustore[result.Ustore.length - 1].t.toFixed(4)}`,
                        autosize: true,
                        paper_bgcolor: "rgba(0,0,0,0)",
                        plot_bgcolor: "rgba(0,0,0,0)",
                      }}
                      style={{ width: "100%", height: 280 }}
                      useResizeHandler
                      config={{ responsive: true }}
                    />
                  </ResponsivePlot>
                </div>
              </>
            )}
          </div>
        </div>
      </div>
    </motion.div>
  );
}

/* 30.3 BTCS (implicit) UI and runner */
function Section303() {
  const [M, setM] = useState(100);
  const [alpha, setAlpha] = useState(0.01);
  const [dt, setDt] = useState(0.001);
  const [tmax, setTmax] = useState(0.1);
  const [ic, setIc] = useState("sine");
  const [run, setRun] = useState(false);
  const [result, setResult] = useState(null);

  useEffect(() => {
    if (!run) return;
    const res = solveBTCS1D({ M: Math.max(10, M), alpha: parseFloat(alpha), dt: parseFloat(dt), tmax: parseFloat(tmax), ic, bc: { left: 0, right: 0 } });
    setResult(res);
    setRun(false);
  }, [run]); // eslint-disable-line

  const stabilityNote = useMemo(() => {
    const { dx } = build1DGrid(Math.max(10, M));
    const r = (alpha * dt) / (dx * dx);
    return { r, stable: true }; // BTCS is unconditionally stable
  }, [M, alpha, dt]);

  return (
    <motion.div className={`${THEME.PANEL} p-4 rounded-2xl`} initial={{ opacity: 0, y: 8 }} whileInView={{ opacity: 1, y: 0 }}>
      <div className="flex items-center justify-between">
        <SectionHeader icon={Code} number="30.3" title="A Simple Implicit Method (BTCS)" accent={THEME.WARN} />
        <div className="text-xs text-zinc-400">Backward Euler — unconditionally stable</div>
      </div>

      <div className="grid gap-4 md:grid-cols-4 mt-4">
        <div className="md:col-span-1 space-y-3">
          <MD>{`BTCS solves a linear system each time-step; stable for any Δt but only first-order accurate in time.`}</MD>
          <div className="grid gap-2">
            <label className="text-xs text-zinc-400">M (intervals)</label>
            <Input className="text-white" type="number" value={M} onChange={(e) => setM(Math.max(10, parseInt(e.target.value || 10, 10)))} />
            <label className="text-xs text-zinc-400">α (diffusivity)</label>
            <Input className="text-white" type="number" value={alpha} onChange={(e) => setAlpha(parseFloat(e.target.value || 0.01))} />
            <label className="text-xs text-zinc-400">Δt</label>
            <Input className="text-white" type="number" value={dt} onChange={(e) => setDt(parseFloat(e.target.value || 0.0001))} />
            <label className="text-xs text-zinc-400">tₘₐₓ</label>
            <Input className="text-white" type="number" value={tmax} onChange={(e) => setTmax(parseFloat(e.target.value || 0.1))} />
            <label className="text-xs text-zinc-400">IC</label>
            <select value={ic} onChange={(e) => setIc(e.target.value)} className="bg-zinc-800/90 text-zinc-100 rounded border border-zinc-700 p-2 w-full">
              <option value="sine">sin(π x)</option>
              <option value="gaussian">Gaussian</option>
            </select>
            <div className="flex gap-2 mt-2">
              <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={() => setRun(true)}>Run BTCS</Button>
              <Button className="bg-white cursor-pointer" variant="ghost" onClick={() => setResult(null)}>Clear</Button>
            </div>
          </div>

          <div className="text-xs text-zinc-400 mt-2">
            <div>Stability: BTCS is unconditionally stable. Current r estimate: <span className="text-zinc-100">{stabilityNote.r.toFixed(5)}</span></div>
          </div>
        </div>

        <div className="md:col-span-3 space-y-3">
          <div className="rounded-lg border border-zinc-700 p-3 bg-zinc-900/50">
            {!result && <div className="text-sm text-zinc-400">Run the implicit solver to see profiles and compare with explicit methods.</div>}
            {result && (
              <>
                <div className="grid md:grid-cols-2 gap-3">
                  <div className="rounded-lg border border-zinc-700 p-2 bg-zinc-950/20">
                    <ResponsivePlot>
                      <Plot
                        data={result.Ustore.map((snap, idx) => ({
                          x: result.xs,
                          y: snap.u,
                          mode: "lines",
                          name: `t=${snap.t.toFixed(4)}`,
                          visible: idx === result.Ustore.length - 1 ? true : "legendonly",
                        }))}
                        layout={{
                          title: "BTCS: snapshots",
                          autosize: true,
                          paper_bgcolor: "rgba(0,0,0,0)",
                          plot_bgcolor: "rgba(0,0,0,0)",
                          legend: { orientation: "h", y: -0.2 },
                        }}
                        style={{ width: "100%", height: 360 }}
                        useResizeHandler
                        config={{ responsive: true }}
                      />
                    </ResponsivePlot>
                  </div>

                  <div className="rounded-lg border border-zinc-700 p-2 bg-zinc-950/20">
                    <ResponsivePlot>
                      <Plot
                        data={[
                          {
                            x: result.xs,
                            y: result.Ustore[result.Ustore.length - 1].u,
                            mode: "lines",
                            name: "BTCS final",
                            line: { color: THEME.PRIMARY },
                          },
                          {
                            x: result.xs,
                            y: result.xs.map((x) => analyticHeatSine(x, result.Ustore[result.Ustore.length - 1].t, alpha)),
                            mode: "lines",
                            name: "analytic",
                            line: { dash: "dash", color: THEME.ACCENT },
                          },
                        ]}
                        layout={{
                          title: "BTCS vs analytic",
                          autosize: true,
                          paper_bgcolor: "rgba(0,0,0,0)",
                          plot_bgcolor: "rgba(0,0,0,0)",
                        }}
                        style={{ width: "100%", height: 360 }}
                        useResizeHandler
                        config={{ responsive: true }}
                      />
                    </ResponsivePlot>
                  </div>
                </div>
              </>
            )}
          </div>
        </div>
      </div>
    </motion.div>
  );
}

/* ---- STOP BLOCK 1 ----
   Paste the next code block (part 2) immediately after this block to complete Chapter30.jsx
 */
/* =========================
   Part 2 of Chapter30.jsx
   (continues from previous block)
   ========================= */

/* 30.4 Crank–Nicolson UI & comparison */
function Section304() {
  const [M, setM] = useState(100);
  const [alpha, setAlpha] = useState(0.01);
  const [dt, setDt] = useState(0.001);
  const [tmax, setTmax] = useState(0.1);
  const [ic, setIc] = useState("sine");
  const [run, setRun] = useState(false);
  const [resFTCS, setResFTCS] = useState(null);
  const [resBTCS, setResBTCS] = useState(null);
  const [resCN, setResCN] = useState(null);

  useEffect(() => {
    if (!run) return;
    const opts = { M: Math.max(10, M), alpha: parseFloat(alpha), dt: parseFloat(dt), tmax: parseFloat(tmax), ic, bc: { left: 0, right: 0 } };
    // run all three solvers for comparison (note: FTCS may be unstable if r large)
    const ftcs = solveFTCS1D(opts);
    const btcs = solveBTCS1D(opts);
    const cn = solveCrankNicolson1D(opts);
    setResFTCS(ftcs);
    setResBTCS(btcs);
    setResCN(cn);
    setRun(false);
  }, [run]); // eslint-disable-line

  const stabilityR = useMemo(() => {
    const { dx } = build1DGrid(Math.max(10, M));
    return (alpha * dt) / (dx * dx);
  }, [M, alpha, dt]);

  return (
    <motion.div className={`${THEME.PANEL} p-4 rounded-2xl`} initial={{ opacity: 0, y: 8 }} whileInView={{ opacity: 1, y: 0 }}>
      <div className="flex items-center justify-between">
        <SectionHeader icon={Layers} number="30.4" title="The Crank–Nicolson Method" accent={THEME.ACCENT} />
        <div className="text-xs text-zinc-400">Second-order accurate in time, unconditionally stable</div>
      </div>

      <div className="grid gap-4 md:grid-cols-4 mt-4">
        <div className="md:col-span-1 space-y-3">
          <MD>{`Crank–Nicolson averages FTCS and BTCS for second-order temporal accuracy. It leads to a tridiagonal solve each time-step.`}</MD>
          <div className="grid gap-2">
            <label className="text-xs text-zinc-400">M</label>
            <Input className="text-white"
 type="number" value={M} onChange={(e) => setM(Math.max(10, parseInt(e.target.value || 10, 10)))} />
            <label className="text-xs text-zinc-400">α</label>
            <Input className="text-white"
 type="number" value={alpha} onChange={(e) => setAlpha(parseFloat(e.target.value || 0.01))} />
            <label className="text-xs text-zinc-400">Δt</label>
            <Input className="text-white"
 type="number" value={dt} onChange={(e) => setDt(parseFloat(e.target.value || 0.001))} />
            <label className="text-xs text-zinc-400">tₘₐₓ</label>
            <Input className="text-white"
 type="number" value={tmax} onChange={(e) => setTmax(parseFloat(e.target.value || 0.1))} />
            <div className="flex gap-2 mt-2">
              <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={() => setRun(true)}>Run FTCS/BTCS/CN</Button>
              <Button className="bg-white cursor-pointer" variant="ghost" onClick={() => { setResFTCS(null); setResBTCS(null); setResCN(null); }}>Clear</Button>
            </div>
          </div>

          <div className="text-xs text-zinc-400 mt-2">
            <div>Stencil stability r ≈ <span className="text-zinc-100">{stabilityR.toFixed(6)}</span></div>
            <div className="text-xs text-zinc-500 mt-1">FTCS may be unstable for large r; BTCS and CN are stable.</div>
          </div>
        </div>

        <div className="md:col-span-3 space-y-3">
          <div className="rounded-lg border border-zinc-700 p-3 bg-zinc-900/50">
            {!resCN && <div className="text-sm text-zinc-400">Run to compare FTCS, BTCS and Crank–Nicolson.</div>}
            {resCN && (
              <>
                <div className="grid md:grid-cols-3 gap-3">
                  <div className="rounded-lg border border-zinc-700 p-2 bg-zinc-950/20">
                    <h4 className="text-sm text-zinc-100 mb-2">FTCS final</h4>
                    <ResponsivePlot>
                      <Plot
                        data={[
                          { x: resFTCS.xs, y: resFTCS.Ustore[resFTCS.Ustore.length - 1].u, mode: "lines", name: "FTCS", line: { color: THEME.PRIMARY } },
                        ]}
                        layout={{
                          autosize: true,
                          paper_bgcolor: "rgba(0,0,0,0)",
                          plot_bgcolor: "rgba(0,0,0,0)",
                        }}
                        style={{ width: "100%", height: 220 }}
                        useResizeHandler
                        config={{ responsive: true }}
                      />
                    </ResponsivePlot>
                  </div>

                  <div className="rounded-lg border border-zinc-700 p-2 bg-zinc-950/20">
                    <h4 className="text-sm text-zinc-100 mb-2">BTCS final</h4>
                    <ResponsivePlot>
                      <Plot
                        data={[
                          { x: resBTCS.xs, y: resBTCS.Ustore[resBTCS.Ustore.length - 1].u, mode: "lines", name: "BTCS", line: { color: THEME.WARN } },
                        ]}
                        layout={{ autosize: true, paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)" }}
                        style={{ width: "100%", height: 220 }}
                        useResizeHandler
                        config={{ responsive: true }}
                      />
                    </ResponsivePlot>
                  </div>

                  <div className="rounded-lg border border-zinc-700 p-2 bg-zinc-950/20">
                    <h4 className="text-sm text-zinc-100 mb-2">Crank–Nicolson final</h4>
                    <ResponsivePlot>
                      <Plot
                        data={[
                          { x: resCN.xs, y: resCN.Ustore[resCN.Ustore.length - 1].u, mode: "lines", name: "CN", line: { color: THEME.ACCENT } },
                        ]}
                        layout={{ autosize: true, paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)" }}
                        style={{ width: "100%", height: 220 }}
                        useResizeHandler
                        config={{ responsive: true }}
                      />
                    </ResponsivePlot>
                  </div>
                </div>

                <div className="mt-3">
                  <h4 className="text-sm text-zinc-100 mb-2">Error vs analytic (sine mode)</h4>
                  <ResponsivePlot>
                    <Plot
                      data={[
                        {
                          x: resCN.xs,
                          y: resFTCS.Ustore[resFTCS.Ustore.length - 1].u.map((u, i) => Math.abs(u - analyticHeatSine(resFTCS.xs[i], resFTCS.Ustore[resFTCS.Ustore.length - 1].t, alpha))),
                          mode: "lines",
                          name: "FTCS error",
                          line: { color: THEME.PRIMARY },
                        },
                        {
                          x: resCN.xs,
                          y: resBTCS.Ustore[resBTCS.Ustore.length - 1].u.map((u, i) => Math.abs(u - analyticHeatSine(resBTCS.xs[i], resBTCS.Ustore[resBTCS.Ustore.length - 1].t, alpha))),
                          mode: "lines",
                          name: "BTCS error",
                          line: { color: THEME.WARN },
                        },
                        {
                          x: resCN.xs,
                          y: resCN.Ustore[resCN.Ustore.length - 1].u.map((u, i) => Math.abs(u - analyticHeatSine(resCN.xs[i], resCN.Ustore[resCN.Ustore.length - 1].t, alpha))),
                          mode: "lines",
                          name: "CN error",
                          line: { color: THEME.ACCENT },
                        },
                      ]}
                      layout={{
                        autosize: true,
                        paper_bgcolor: "rgba(0,0,0,0)",
                        plot_bgcolor: "rgba(0,0,0,0)",
                        yaxis: { type: "log", title: "|error|" },
                        margin: { t: 10, b: 30, l: 50, r: 10 },
                      }}
                      style={{ width: "100%", height: 300 }}
                      useResizeHandler
                      config={{ responsive: true }}
                    />
                  </ResponsivePlot>
                </div>
              </>
            )}
          </div>
        </div>
      </div>
    </motion.div>
  );
}

/* 30.5 Parabolic Equations in Two Spatial Dimensions (2D demo) */
function Section305() {
  const [N, setN] = useState(60);
  const [alpha, setAlpha] = useState(0.1);
  const [dt, setDt] = useState(1e-5);
  const [tmax, setTmax] = useState(0.005);
  const [run, setRun] = useState(false);
  const [result2D, setResult2D] = useState(null);

  useEffect(() => {
    if (!run) return;
    const res = solveFTCS2D({ N: Math.max(10, N), alpha: parseFloat(alpha), dt: parseFloat(dt), tmax: parseFloat(tmax), bc: { value: 0 } });
    setResult2D(res);
    setRun(false);
  }, [run]); // eslint-disable-line

  const stabilityInfo = useMemo(() => {
    const dx = 1 / Math.max(10, N);
    const r = (alpha * dt) / (dx * dx);
    return r;
  }, [N, alpha, dt]);

  return (
    <motion.div className={`${THEME.PANEL} p-4 rounded-2xl`} initial={{ opacity: 0, y: 8 }} whileInView={{ opacity: 1, y: 0 }}>
      <div className="flex items-center justify-between">
        <SectionHeader icon={Grid} number="30.5" title="Parabolic Equations in Two Spatial Dimensions" accent={THEME.PRIMARY} />
        <div className="text-xs text-zinc-400">Explicit 2D FTCS (demonstration)</div>
      </div>

      <div className="grid gap-4 md:grid-cols-4 mt-4">
        <div className="md:col-span-1 space-y-3">
          <MD>
  {`2D heat: $u_t = \\alpha (u_{xx} + u_{yy})$. Explicit scheme has severe stability constraint ($r \\le 1/4$ for square grids).`}
</MD>

          <div className="grid gap-2">
            <label className="text-xs text-zinc-400">N (intervals)</label>
            <Input className="text-white" type="number" value={N} onChange={(e) => setN(Math.max(10, parseInt(e.target.value || 10, 10)))} />
            <label className="text-xs text-zinc-400">α</label>
            <Input className="text-white" type="number" value={alpha} onChange={(e) => setAlpha(parseFloat(e.target.value || 0.1))} />
            <label className="text-xs text-zinc-400">Δt</label>
            <Input className="text-white" type="number" value={dt} onChange={(e) => setDt(parseFloat(e.target.value || 1e-5))} />
            <label className="text-xs text-zinc-400">tₘₐₓ</label>
            <Input className="text-white" type="number" value={tmax} onChange={(e) => setTmax(parseFloat(e.target.value || 0.005))} />
            <div className="flex gap-2 mt-2">
              <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={() => setRun(true)}>Run 2D FTCS</Button>
              <Button className="bg-white cursor-pointer" variant="ghost" onClick={() => setResult2D(null)}>Clear</Button>
            </div>
          </div>

          <div className="text-xs text-zinc-400 mt-2">
            <div>2D stability parameter r = α Δt / Δx² ≈ <span className="text-zinc-100">{stabilityInfo.toFixed(6)}</span></div>
            <div className="text-xs text-zinc-500">Stability target r ≤ 0.25 for square grids (explicit).</div>
          </div>
        </div>

        <div className="md:col-span-3 space-y-3">
          <div className="rounded-lg border border-zinc-700 p-3 bg-zinc-900/50">
            {!result2D && <div className="text-sm text-zinc-400">Run the 2D explicit simulation (may be slow for large N).</div>}
            {result2D && (
              <>
                <div className="grid md:grid-cols-2 gap-3">
                  <div className="rounded-lg border border-zinc-700 p-2 bg-zinc-950/20">
                    <ResponsivePlot>
                      <Plot
                        data={[
                          {
                            z: grid2Plotly(result2D.snapshots[result2D.snapshots.length - 1].u, result2D.xs).Z,
                            x: result2D.xs,
                            y: result2D.xs,
                            type: "surface",
                            colorscale: "Viridis",
                            showscale: true,
                          },
                        ]}
                        layout={{
                          title: `2D surface at t=${result2D.snapshots[result2D.snapshots.length - 1].t.toFixed(4)}`,
                          autosize: true,
                          paper_bgcolor: "rgba(0,0,0,0)",
                          plot_bgcolor: "rgba(0,0,0,0)",
                          margin: { t: 30 },
                        }}
                        style={{ width: "100%", height: 420 }}
                        useResizeHandler
                        config={{ responsive: true }}
                      />
                    </ResponsivePlot>
                  </div>

                  <div className="rounded-lg border border-zinc-700 p-2 bg-zinc-950/20">
                    <ResponsivePlot>
                      <Plot
                        data={[
                          {
                            z: grid2Plotly(result2D.snapshots[Math.max(0, result2D.snapshots.length - 1)].u, result2D.xs).Z,
                            x: result2D.xs,
                            y: result2D.xs,
                            type: "contour",
                            colorscale: "Viridis",
                            contours: { coloring: "fill" },
                          },
                        ]}
                        layout={{
                          title: "Contour plot",
                          autosize: true,
                          paper_bgcolor: "rgba(0,0,0,0)",
                          plot_bgcolor: "rgba(0,0,0,0)",
                        }}
                        style={{ width: "100%", height: 420 }}
                        useResizeHandler
                        config={{ responsive: true }}
                      />
                    </ResponsivePlot>
                  </div>
                </div>

                <div className="mt-3 text-sm text-zinc-300">
                  Snapshots stored: <span className="text-zinc-100">{result2D.snapshots.length}</span>. Stability r: <span className="text-zinc-100">{result2D.r.toFixed(6)}</span>.
                </div>
              </>
            )}
          </div>
        </div>
      </div>
    </motion.div>
  );
}

/* Problems section */
function Problems30() {
  const [open, setOpen] = useState(-1);
  const problems = [
    {
      id: 1,
      title: "Stability of FTCS",
      prompt: "Derive the stability condition r ≤ 1/2 for FTCS in 1D by von Neumann analysis.",
      solution: "Perform Fourier mode substitution u_j^n = ξ^n e^{ikjΔx}; solve for amplification factor |ξ| ≤ 1 leading to r ≤ 1/2.",
    },
    {
      id: 2,
      title: "Convergence order",
      prompt: "Verify temporal order of BTCS (first) and Crank–Nicolson (second) numerically with mesh refinement.",
      solution: "Compute solutions at a fixed time with refined Δt and compare L∞ errors; estimate slopes on log-log plot.",
    },
    {
      id: 3,
      title: "2D scheme stability",
      prompt: "Show that the explicit 2D FTCS stability requires r ≤ 1/4 for square grids and explain physical intuition.",
      solution: "Von Neumann analysis on 2D Fourier modes; four neighbors reduce allowable r by factor 2 relative to 1D.",
    },
    {
      id: 4,
      title: "Implicit solver with variable coefficients",
      prompt: "Extend BTCS to u_t = (α(x) u_x)_x and discuss how the tridiagonal system changes (non-constant coefficients).",
      solution: "Discrete operator yields nonuniform tridiagonal with entries dependent on α_{i+1/2}, use generalized Thomas solver or sparse direct solver.",
    },
  ];

  return (
    <motion.div className={`${THEME.PANEL} p-4 rounded-2xl`} initial={{ opacity: 0, y: 8 }} whileInView={{ opacity: 1, y: 0 }}>
      <div className="flex items-center justify-between">
        <SectionHeader icon={ListChecks} number="Problems" title="Exercises" />
        <div className="text-xs text-zinc-400">Try implementing & verifying these numerically</div>
      </div>

      <div className="mt-4 space-y-3">
        {problems.map((p, i) => (
          <div key={p.id} className="rounded-lg border border-zinc-700 p-3 bg-zinc-900/50">
            <div className="flex items-center justify-between">
              <div>
                <div className="text-sm text-zinc-100 font-medium">{p.title}</div>
                <div className="text-xs text-zinc-400">{p.prompt}</div>
              </div>
              <div>
                <Button className="bg-white cursor-pointer" variant="ghost" onClick={() => setOpen(open === i ? -1 : i)}>{open === i ? "Hide" : "Show"}</Button>
              </div>
            </div>
            {open === i && (
              <div className="mt-3 bg-zinc-950/30 border border-zinc-800 rounded p-3 text-sm text-zinc-200">
                <div className="font-semibold mb-1">Solution sketch</div>
                <div>{p.solution}</div>
              </div>
            )}
          </div>
        ))}
      </div>
    </motion.div>
  );
}

/* Documentation / notes panel */
function ChapterDocs30() {
  const [active, setActive] = useState(0);
  const sections = [
    {
      title: "30.1 Heat Equation",
      notes: ["Diffusion smooths and damps initial data; use analytic modes for verification.", "Careful with units: α has units length^2 / time."],
    },
    {
      title: "30.2 FTCS",
      notes: ["Simple and explicit; conditionally stable. Use for quick prototyping only.", "Stability limit: r ≤ 1/2 (1D)."],
    },
    {
      title: "30.3 BTCS",
      notes: ["Implicit and unconditionally stable, first-order in time; requires linear solve each time step.", "Use tridiagonal solver for efficiency."],
    },
    {
      title: "30.4 Crank–Nicolson",
      notes: ["Second-order accurate in time, unconditionally stable; oscillatory behavior possible for discontinuous ICs.", "Common choice for accuracy/stability tradeoff."],
    },
    {
      title: "30.5 2D Parabolic",
      notes: ["Explicit 2D FTCS stability is stricter (r ≤ 1/4).", "For large-scale 2D problems use implicit schemes, ADI, or multigrid in time/space."],
    },
  ];

  return (
    <Card className={`${THEME.PANEL} rounded-2xl`}>
      <CardContent className="p-4">
        <div className="flex items-center justify-between mb-3">
          <div className="flex items-center gap-2">
            <BookOpen className="w-5 h-5 text-cyan-400" />
            <div>
              <div className="text-sm text-zinc-100 font-semibold">Chapter 30 — Notes & References</div>
              <div className="text-xs text-zinc-400">Key takeaways & references</div>
            </div>
          </div>
          <div className="text-xs text-zinc-400">Strikwerda, Morton & Mayers, LeVeque</div>
        </div>

        <div className="grid md:grid-cols-4 gap-3">
          <div className="space-y-2">
            {sections.map((s, i) => (
              <button key={i} onClick={() => setActive(i)} className={`w-full cursor-pointer text-gray-200 text-left px-3 py-2 rounded border ${i === active ? "bg-cyan-900/30 border-cyan-600" : "bg-zinc-900/40 border-zinc-700"}`}>
                {s.title}
              </button>
            ))}
          </div>

          <div className="md:col-span-3 bg-zinc-950/30 border border-zinc-800 rounded p-4">
            <div className="text-sm font-semibold text-zinc-100 mb-2">{sections[active].title}</div>
            <ul className="text-sm text-zinc-300 list-disc ml-5 space-y-1">
              {sections[active].notes.map((n, idx) => <li key={idx}>{n}</li>)}
            </ul>
            <div className="text-xs text-zinc-500 mt-3">References: Strikwerda (Finite Difference Schemes), Morton & Mayers (Numerical Solution of PDEs), LeVeque (Finite Difference Methods for ODEs & PDEs).</div>
          </div>
        </div>
      </CardContent>
    </Card>
  );
}

/* Assemble Chapter30 page */
export default function Chapter30() {
  useEffect(() => {
    if (typeof window !== "undefined") window.scrollTo({ top: 0, behavior: "smooth" });
  }, []);

  return (
    <div className={`p-4 sm:p-6 lg:p-8 space-y-6 ${THEME.BG} min-h-screen`}>
      <motion.div initial={{ opacity: 0, y: -8 }} animate={{ opacity: 1, y: 0 }}>
        <div className="flex flex-col sm:flex-row sm:items-center sm:justify-between gap-4">
          <div>
           
            <h1 className="text-2xl sm:text-3xl font-bold text-emerald-400">Finite Difference: Parabolic Equations</h1>
            <div className="text-zinc-400 mt-1">30.1 • 30.2 • 30.3 • 30.4 • 30.5 — explicit & implicit schemes</div>
          </div>

        
        </div>
      </motion.div>

      <div className="space-y-6">
        <Section301 />
        <Section302 />
        <Section303 />
        <Section304 />
        <Section305 />
        <Problems30 />
        <ChapterDocs30 />
        <BottomBar/>
      </div>

      <footer className="text-xs text-zinc-500 mt-6">Generated Chapter 30 module — interactive finite-difference solvers and responsive visualizations.</footer>
    </div>
  );
}
