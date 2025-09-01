// src/pages/Chapter32.jsx
"use client";

/**
 * CHAPTER 32 — Case Studies: Partial Differential Equations
 *
 * Sections:
 *   32.1 One-Dimensional Mass Balance of a Reactor (Chemical/Bio Engineering)
 *   32.2 Deflections of a Plate (Civil/Environmental Engineering)
 *   32.3 Two-Dimensional Electrostatic Field Problems (Electrical Engineering)
 *   32.4 Finite-Element Solution of a Series of Springs (Mechanical/Aerospace Engineering)
 *   Problems
 *
 * Features (pro-level, responsive):
 *  - Dark Tailwind-like UI (utility classes assumed in host project)
 *  - KaTeX math rendering via react-markdown + remark-math + rehype-katex
 *  - Plotly.js responsive plots (useResizeHandler + config.responsive)
 *  - Framer Motion animations
 *  - Multiple specialized solvers:
 *      * 1D transient reactor (finite-difference RK4 and analytic where available)
 *      * 2D plate deflection via finite difference (biharmonic approx / simple thin-plate)
 *      * 2D Laplace relaxation for electrostatic potential and field lines
 *      * 1D spring chain solved with FEM-like stiffness assembly
 *  - Many explanatory comments, split JSX props for readability
 *  - Responsive layout: stacks on small screens, multi-column on larger
 *
 * NOTES:
 *  - Intended for educational/demo use; dense/sparse solvers are naive and for clarity.
 *  - Assumes dependencies: react, react-dom, react-markdown, remark-math, rehype-katex, katex, react-plotly.js, plotly.js, framer-motion, lucide-react
 *  - Replace UI primitives if your project provides different components.
 */

/* ------------------------------------------------------------------ */
/* Imports                                                           */
/* ------------------------------------------------------------------ */
import React, { useEffect, useMemo, useState } from "react";
import Plot from "react-plotly.js";
import { motion } from "framer-motion";
import ReactMarkdown from "react-markdown";
import remarkMath from "remark-math";
import rehypeKatex from "rehype-katex";
import "katex/dist/katex.min.css";
import BottomBar from "../components/BottomBar";
import {
  BookOpen,
  Activity,
  Layers,
  Zap,
  Box,
  ListChecks,
  BarChart3,
} from "lucide-react";

/* ------------------------------------------------------------------ */
/* UI primitive fallbacks (simple)                                   */
/* ------------------------------------------------------------------ */
let Input = (props) => (
  <input
    {...props}
    className={`bg-zinc-800/80 text-white rounded p-2 w-full ${props.className || ""}`}
  />
);
let Button = (props) => (
  <button
    {...props}
    className={`px-3 py-1 rounded bg-zinc-800/70 border border-zinc-700 text-sm text-white hover:bg-zinc-700 transition ${props.className || ""}`}
  >
    {props.children}
  </button>
);
let Textarea = (props) => (
  <textarea
    {...props}
    className={`bg-zinc-800/80 text-white rounded p-2 w-full ${props.className || ""}`}
  />
);

/* ------------------------------------------------------------------ */
/* Theme and helpers                                                  */
/* ------------------------------------------------------------------ */
const THEME = {
  BG: "bg-zinc-950",
  PANEL: "bg-zinc-900/60 border border-zinc-700",
  ACCENT: "#22d3ee",
  ACCENT2: "#34d399",
};

function MD({ children }) {
  const components = {
    p: ({ node, ...props }) => <p className="text-sm text-zinc-300 leading-relaxed" {...props} />,
    li: ({ node, ...props }) => <li className="ml-4 text-sm text-zinc-300" {...props} />,
    code: ({ node, inline, ...props }) =>
      inline ? (
        <code className="bg-zinc-800/70 px-1 rounded text-xs text-cyan-200" {...props} />
      ) : (
        <pre className="bg-zinc-800/60 rounded p-2 text-xs overflow-auto" {...props} />
      ),
    h1: ({ node, ...props }) => <h1 className="text-xl text-zinc-100 font-semibold" {...props} />,
    h2: ({ node, ...props }) => <h2 className="text-lg text-zinc-100 font-semibold" {...props} />,
  };
  return (
    <ReactMarkdown
      remarkPlugins={[remarkMath]}
      rehypePlugins={[rehypeKatex]}
      components={components}
    >
      {children}
    </ReactMarkdown>
  );
}

/* Section header */
function SectionHeader({ icon: Icon, number, title, accent = THEME.ACCENT }) {
  return (
    <div className="flex items-center gap-3">
      <div className="p-2 rounded-lg bg-zinc-800/20">
        <Icon className="w-6 h-6" style={{ color: accent }} />
      </div>
      <div>
        <div className="text-xs text-zinc-400 uppercase tracking-wider">{number}</div>
        <div className="text-lg font-semibold text-zinc-100">{title}</div>
      </div>
    </div>
  );
}

/* Responsive plot wrapper */
function ResponsivePlot({ children, height = 360 }) {
  return (
    <div className="w-full max-w-full overflow-hidden rounded-lg">
      <div style={{ minHeight: Math.max(120, height) }}>{children}</div>
    </div>
  );
}

/* ------------------------------------------------------------------ */
/* Numerical helpers (finite differences, RK4, small linear solves)   */
/* ------------------------------------------------------------------ */

/* linspace */
function linspace(a, b, n) {
  if (n <= 1) return [a];
  const out = [];
  const dx = (b - a) / (n - 1);
  for (let i = 0; i < n; i++) out.push(a + i * dx);
  return out;
}

/* Simple RK4 integrator for ODE y' = f(t,y) with fixed step */
function rk4Step(y, t, dt, f) {
  const k1 = f(t, y);
  const k2 = f(t + dt / 2, addVec(y, mulVec(k1, dt / 2)));
  const k3 = f(t + dt / 2, addVec(y, mulVec(k2, dt / 2)));
  const k4 = f(t + dt, addVec(y, mulVec(k3, dt)));
  const yNext = addVec(y, mulVec(addVec(addVec(k1, mulVec(k2, 2)), addVec(mulVec(k3, 2), k4)), dt / 6));
  return yNext;
}

/* vector helpers */
function addVec(a, b) {
  if (typeof a === "number") return a + b;
  return a.map((v, i) => v + b[i]);
}
function mulVec(a, scalar) {
  if (typeof a === "number") return a * scalar;
  return a.map((v) => v * scalar);
}

/* Simple tridiagonal solver for 1D implicit steps */
function solveTridiagonal(a, b, c, d) {
  // a: sub-diagonal (length n-1), b: diagonal (length n), c: super-diagonal (length n-1), d: RHS (length n)
  const n = b.length;
  const cp = new Array(n - 1);
  const dp = new Array(n);
  cp[0] = c[0] / b[0];
  dp[0] = d[0] / b[0];
  for (let i = 1; i < n - 1; i++) {
    const m = b[i] - a[i - 1] * cp[i - 1];
    cp[i] = c[i] / m;
    dp[i] = (d[i] - a[i - 1] * dp[i - 1]) / m;
  }
  dp[n - 1] = (d[n - 1] - a[n - 2] * dp[n - 2]) / (b[n - 1] - a[n - 2] * cp[n - 2]);
  const x = new Array(n);
  x[n - 1] = dp[n - 1];
  for (let i = n - 2; i >= 0; i--) {
    x[i] = dp[i] - cp[i] * x[i + 1];
  }
  return x;
}

/* ------------------------------------------------------------------ */
/* 32.1 — 1D Mass Balance of a Reactor                                */
/* ------------------------------------------------------------------ */
/**
 * Model: dC/dt + v dC/dx = -k1 C - k2 C^2 (advection + reactions)
 * For educational simplicity we will implement:
 *  - pure reaction (ODE) at a point with analytic solutions for 1st and 2nd order
 *  - plug-flow reactor using 1D finite differences for spatial discretization (explicit & implicit)
 *
 * Analytic:
 *  - First-order: C(t) = C0 * exp(-k1 t)
 *  - Second-order: C(t) = 1 / (1/C0 + k2 t)
 *
 * We'll compare RK4 time integration (explicit) vs analytic where applicable.
 */

/* analytic solutions */
function analyticFirstOrder(C0, k1, t) {
  return C0 * Math.exp(-k1 * t);
}
function analyticSecondOrder(C0, k2, t) {
  return 1 / (1 / C0 + k2 * t);
}

/* simple reactor transient simulator (well-mixed CSTR style) for demonstration */
function simulateCSTR({ C0, k1, k2, tFinal, dt }) {
  const steps = Math.max(1, Math.ceil(tFinal / dt));
  const t = linspace(0, tFinal, steps + 1);
  const C = new Array(steps + 1).fill(0);
  C[0] = C0;
  for (let n = 0; n < steps; n++) {
    // ODE: dC/dt = -k1 C - k2 C^2
    const f = (tt, y) => -k1 * y - k2 * y * y;
    const y = C[n];
    // scalar RK4
    const k1r = f(t[n], y);
    const k2r = f(t[n] + dt / 2, y + k1r * dt / 2);
    const k3r = f(t[n] + dt / 2, y + k2r * dt / 2);
    const k4r = f(t[n] + dt, y + k3r * dt);
    const yNext = y + (dt / 6) * (k1r + 2 * k2r + 2 * k3r + k4r);
    C[n + 1] = Math.max(0, yNext);
  }
  return { t, C };
}

/* 1D plug-flow reactor finite-difference (explicit upwind advection + reaction) */
function simulatePlugFlow({
  L = 1.0,
  Nx = 100,
  v = 1.0,
  k1 = 0.0,
  k2 = 0.0,
  C_in = 1.0,
  tFinal = 1.0,
  dt = 0.001,
}) {
  const dx = L / (Nx - 1);
  const x = linspace(0, L, Nx);
  const steps = Math.max(1, Math.ceil(tFinal / dt));
  // concentration array, initial zero except inlet
  let C = new Array(Nx).fill(0);
  // for periodic or continuous inlet, set C[0] = C_in always (Dirichlet at inlet)
  const recordT = [];
  const recordC = [];
  for (let n = 0; n <= steps; n++) {
    const t = n * dt;
    if (n % Math.max(1, Math.floor(steps / 50)) === 0) {
      recordT.push(t);
      recordC.push(C.slice());
    }
    // advance using explicit upwind
    const Cnew = C.slice();
    for (let i = 1; i < Nx; i++) {
      // advection: upwind difference
      const adv = -v * (C[i] - C[i - 1]) / dx;
      const reac = -k1 * C[i] - k2 * C[i] * C[i];
      Cnew[i] = C[i] + dt * (adv + reac);
      if (Cnew[i] < 0) Cnew[i] = 0;
    }
    // inlet Dirichlet
    Cnew[0] = C_in;
    C = Cnew;
  }
  return { x, recordT, recordC };
}

/* ------------------------------------------------------------------ */
/* 32.2 — Plate deflection (approximate via finite difference)        */
/* ------------------------------------------------------------------ */
/**
 * Thin plate equation (small deflection, simply supported approximation):
 *   D ∇^4 w = q(x,y)
 * where D is flexural rigidity and ∇^4 is biharmonic operator.
 *
 * For simplicity, implement a low-order finite-difference approximate using
 * successive application of Laplacian (discrete biharmonic = Laplacian(Laplacian(w))).
 *
 * We solve (with relaxation) for a plate on [0,1]x[0,1] under q(x,y).
 */

function solvePlateBiharmonic({
  Nx = 41,
  Ny = 41,
  qFunc = (x, y) => 1.0,
  D = 1.0,
  tol = 1e-4,
  maxIter = 2000,
}) {
  // grid
  const dx = 1 / (Nx - 1);
  const dy = 1 / (Ny - 1);
  const nodes = [];
  for (let j = 0; j < Ny; j++) {
    for (let i = 0; i < Nx; i++) {
      nodes.push({ x: i * dx, y: j * dy });
    }
  }
  // initialize w (deflection)
  const w = Array.from({ length: Ny }, () => Array(Nx).fill(0));
  // boundary conditions: simply supported w=0 on boundary and moment-free? We'll enforce w=0 Dirichlet for demo.
  // RHS f = q/D
  const f = (i, j) => qFunc(i * dx, j * dy) / D;

  // helper: discrete Laplacian (5-point)
  function laplacian(arr, i, j) {
    const xm = i === 0 ? arr[j][i] : arr[j][i - 1];
    const xp = i === Nx - 1 ? arr[j][i] : arr[j][i + 1];
    const ym = j === 0 ? arr[j][i] : arr[j - 1][i];
    const yp = j === Ny - 1 ? arr[j][i] : arr[j + 1][i];
    const center = arr[j][i];
    return (xp + xm - 2 * center) / (dx * dx) + (yp + ym - 2 * center) / (dy * dy);
  }

  // Relaxation solve for Δ^2 w = f  => Δ(Δ w) = f
  // Approach: introduce v = Δ w, solve Δ v = f, and enforce v = Δ w iteratively
  const v = Array.from({ length: Ny }, () => Array(Nx).fill(0));

  // Initialize v as laplacian of initial w (zero)
  for (let j = 0; j < Ny; j++) {
    for (let i = 0; i < Nx; i++) {
      v[j][i] = laplacian(w, i, j);
    }
  }

  // Iterative relaxation (Jacobi) for coupled system — naive but illustrative
  let iter = 0;
  let err = 1;
  while (iter < maxIter && err > tol) {
    err = 0;
    // update v by solving Δ v = f  (Jacobi)
    const vNew = v.map((row) => row.slice());
    for (let j = 1; j < Ny - 1; j++) {
      for (let i = 1; i < Nx - 1; i++) {
        const rhs = f(i, j);
        const vp = v[j][i + 1];
        const vm = v[j][i - 1];
        const up = v[j + 1][i];
        const um = v[j - 1][i];
        const newv = ((vp + vm) / (dx * dx) + (up + um) / (dy * dy) - rhs) / ( -2 / (dx*dx) - 2 / (dy*dy));
        const de = Math.abs(newv - v[j][i]);
        if (de > err) err = de;
        vNew[j][i] = newv;
      }
    }
    // update w by solving Δ w = v (Jacobi)
    const wNew = w.map((row) => row.slice());
    for (let j = 1; j < Ny - 1; j++) {
      for (let i = 1; i < Nx - 1; i++) {
        const rhs = vNew[j][i];
        const wp = w[j][i + 1];
        const wm = w[j][i - 1];
        const up = w[j + 1][i];
        const um = w[j - 1][i];
        const neww = ((wp + wm) / (dx * dx) + (up + um) / (dy * dy) - rhs) / ( -2 / (dx*dx) - 2 / (dy*dy));
        const de = Math.abs(neww - w[j][i]);
        if (de > err) err = de;
        wNew[j][i] = neww;
      }
    }
    // enforce Dirichlet 0 on boundaries
    for (let i = 0; i < Nx; i++) {
      wNew[0][i] = 0; wNew[Ny - 1][i] = 0;
      vNew[0][i] = 0; vNew[Ny - 1][i] = 0;
    }
    for (let j = 0; j < Ny; j++) {
      wNew[j][0] = 0; wNew[j][Nx - 1] = 0;
      vNew[j][0] = 0; vNew[j][Nx - 1] = 0;
    }
    // swap
    for (let j = 0; j < Ny; j++) {
      for (let i = 0; i < Nx; i++) {
        w[j][i] = wNew[j][i];
        v[j][i] = vNew[j][i];
      }
    }
    iter++;
  }

  // Prepare arrays for plotting
  const xs = linspace(0, 1, Nx);
  const ys = linspace(0, 1, Ny);
  const Z = w.map((row) => row.slice());

  return { xs, ys, Z, iterations: iter, residual: err };
}

/* ------------------------------------------------------------------ */
/* 32.3 — 2D Electrostatics (Laplace)                                 */
/* ------------------------------------------------------------------ */
/**
 * Solve Laplace's equation (∇^2 φ = 0) on a rectangular grid with Dirichlet BCs.
 * Use successive over-relaxation (SOR) for speed.
 *
 * Then compute electric field E = -∇φ and show field lines / equipotential contours.
 */

function solveLaplaceSOR({
  Nx = 80,
  Ny = 80,
  phiBoundaryFunc = (x, y) => 0, // default zero boundary
  tol = 1e-5,
  maxIter = 5000,
  omega = 1.7, // relaxation parameter
}) {
  const dx = 1 / (Nx - 1);
  const dy = 1 / (Ny - 1);
  const phi = Array.from({ length: Ny }, () => Array(Nx).fill(0));

  // apply boundary conditions from phiBoundaryFunc (Dirichlet)
  for (let j = 0; j < Ny; j++) {
    for (let i = 0; i < Nx; i++) {
      const x = i * dx;
      const y = j * dy;
      if (i === 0 || i === Nx - 1 || j === 0 || j === Ny - 1) {
        phi[j][i] = phiBoundaryFunc(x, y);
      }
    }
  }

  let iter = 0;
  let maxDiff = 1;
  while (iter < maxIter && maxDiff > tol) {
    maxDiff = 0;
    for (let j = 1; j < Ny - 1; j++) {
      for (let i = 1; i < Nx - 1; i++) {
        // skip boundary nodes
        const newPhi = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i] + phi[j - 1][i]) / 4;
        const updated = (1 - omega) * phi[j][i] + omega * newPhi;
        const diff = Math.abs(updated - phi[j][i]);
        if (diff > maxDiff) maxDiff = diff;
        phi[j][i] = updated;
      }
    }
    iter++;
  }

  // compute E field: Ex = -dφ/dx, Ey = -dφ/dy (central differences)
  const Ex = Array.from({ length: Ny }, () => Array(Nx).fill(0));
  const Ey = Array.from({ length: Ny }, () => Array(Nx).fill(0));
  for (let j = 1; j < Ny - 1; j++) {
    for (let i = 1; i < Nx - 1; i++) {
      Ex[j][i] = -(phi[j][i + 1] - phi[j][i - 1]) / (2 * dx);
      Ey[j][i] = -(phi[j + 1][i] - phi[j - 1][i]) / (2 * dy);
    }
  }

  return { phi, Ex, Ey, iterations: iter, residual: maxDiff };
}

/* ------------------------------------------------------------------ */
/* 32.4 — Series of Springs (1D Finite-Element like)                   */
/* ------------------------------------------------------------------ */
/**
 * A chain of springs and masses: K matrix assembly and solve for static displacements
 * For simplicity: masses not used in static case; nodes connected by linear springs with stiffness k_i
 *
 * Boundary: first node fixed (u0=0), last node free or fixed depending on BC
 */

function assembleSpringChain(kList, fList, bc = { left: 0, right: null }) {
  // number of springs = n; number of nodes = n+1
  const n = kList.length;
  const N = n + 1;
  const K = Array.from({ length: N }, () => Array(N).fill(0));
  // assemble
  for (let i = 0; i < n; i++) {
    const k = kList[i];
    K[i][i] += k;
    K[i][i + 1] += -k;
    K[i + 1][i] += -k;
    K[i + 1][i + 1] += k;
  }
  // apply BC: left fixed u0 = bc.left
  const F = fList.slice();
  // enforce left Dirichlet
  for (let j = 0; j < N; j++) K[0][j] = 0;
  K[0][0] = 1;
  F[0] = bc.left;
  // if right BC is provided
  if (bc.right !== null) {
    for (let j = 0; j < N; j++) K[N - 1][j] = 0;
    K[N - 1][N - 1] = 1;
    F[N - 1] = bc.right;
  }
  // solve (naive dense)
  const U = solveDense(K, F);
  return { K, F, U };
}

/* ------------------------------------------------------------------ */
/* UI & Interaction Components                                       */
/* ------------------------------------------------------------------ */

/* 32.1 UI Component */
function ReactorSection() {
  // parameters
  const [C0, setC0] = useState(1.0);
  const [k1, setK1] = useState(0.5);
  const [k2, setK2] = useState(0.0);
  const [tFinal, setTFinal] = useState(5.0);
  const [dt, setDt] = useState(0.01);
  const [run, setRun] = useState(false);
  const [result, setResult] = useState(null);

  useEffect(() => {
    if (!run) return;
    // simulate CSTR ODE
    const out = simulateCSTR({
      C0: Number(C0),
      k1: Number(k1),
      k2: Number(k2),
      tFinal: Number(tFinal),
      dt: Number(dt),
    });
    setResult(out);
    setRun(false);
  }, [run]);

  // analytic if purely first-order or purely second-order where applicable
  const analyticSeries = useMemo(() => {
    if (!result) return null;
    const t = result.t;
    const Cfirst = t.map((tt) => analyticFirstOrder(C0, k1, tt));
    const Csecond = t.map((tt) => analyticSecondOrder(C0, k2, tt));
    return { t, Cfirst, Csecond };
  }, [result, C0, k1, k2]);

  return (
    <motion.div
      initial={{ opacity: 0, y: 6 }}
      whileInView={{ opacity: 1, y: 0 }}
      className={`${THEME.PANEL} p-4 rounded-2xl`}
    >
      <div className="flex items-center justify-between">
        <SectionHeader
          icon={Activity}
          number="32.1"
          title="1D Mass Balance: Reactor (CSTR / Plug-Flow demos)"
          accent={THEME.ACCENT}
        />
        <div className="text-xs text-zinc-400">Transient reaction kinetics & transport</div>
      </div>

      <div className="mt-3 grid gap-4 md:grid-cols-3">
        <div className="md:col-span-1 space-y-3">
          <div>
            <label className="text-xs text-zinc-400">Initial concentration C₀</label>
            <Input
              type="number"
              step="0.1"
              value={C0}
              onChange={(e) => setC0(e.target.value)}
            />
          </div>

          <div>
            <label className="text-xs text-zinc-400">First-order rate k₁</label>
            <Input
              type="number"
              step="0.01"
              value={k1}
              onChange={(e) => setK1(e.target.value)}
            />
          </div>

          <div>
            <label className="text-xs text-zinc-400">Second-order rate k₂</label>
            <Input
              type="number"
              step="0.01"
              value={k2}
              onChange={(e) => setK2(e.target.value)}
            />
          </div>

          <div>
            <label className="text-xs text-zinc-400">t final (s)</label>
            <Input
              type="number"
              step="0.1"
              value={tFinal}
              onChange={(e) => setTFinal(e.target.value)}
            />
          </div>

          <div>
            <label className="text-xs text-zinc-400">time step Δt</label>
            <Input
              type="number"
              step="0.001"
              value={dt}
              onChange={(e) => setDt(e.target.value)}
            />
          </div>

          <div className="flex gap-2 mt-2">
            <Button onClick={() => setRun(true)}>Run RK4</Button>
            <Button onClick={() => { setResult(null); setRun(false); }}>Clear</Button>
          </div>

          <div className="text-xs text-zinc-400 mt-2">
            When k₂=0 analytic solution is exponential. When k₁=0 analytic follows second-order formula.
          </div>
        </div>

        <div className="md:col-span-2">
          <div className="rounded border border-zinc-700 p-3 bg-zinc-900/40">
            {!result && <div className="text-sm text-zinc-400">No simulation yet — run to see concentration vs time.</div>}
            {result && analyticSeries && (
              <>
                <ResponsivePlot height={360}>
                  <Plot
                    data={[
                      {
                        x: result.t,
                        y: result.C,
                        mode: "lines+markers",
                        name: "RK4 numeric",
                        line: { color: THEME.ACCENT },
                      },
                      k1 > 0 && k2 == 0 && {
                        x: analyticSeries.t,
                        y: analyticSeries.Cfirst,
                        mode: "lines",
                        name: "Analytic (1st order)",
                        line: { dash: "dash", color: THEME.ACCENT2 },
                      },
                      k2 > 0 && k1 == 0 && {
                        x: analyticSeries.t,
                        y: analyticSeries.Csecond,
                        mode: "lines",
                        name: "Analytic (2nd order)",
                        line: { dash: "dash", color: "#f97316" },
                      },
                    ].filter(Boolean)}
                    layout={{
                      title: "Concentration vs time",
                      autosize: true,
                      paper_bgcolor: "rgba(0,0,0,0)",
                      plot_bgcolor: "rgba(0,0,0,0)",
                      margin: { t: 30, b: 40, l: 50, r: 20 },
                      font: { color: "#cbd5e1" },
                      legend: { orientation: "h" },
                    }}
                    style={{ width: "100%", height: 360 }}
                    useResizeHandler
                    config={{ responsive: true }}
                  />
                </ResponsivePlot>

                <div className="mt-3 grid grid-cols-1 sm:grid-cols-2 gap-3">
                  <div className="text-sm text-zinc-300">
                    <div>Final concentration: <span className="text-zinc-100">{result.C[result.C.length - 1].toExponential(3)}</span></div>
                    <div className="text-xs text-zinc-400 mt-1">Time steps: {result.t.length}</div>
                  </div>

                  <div className="text-sm text-zinc-300">
                    <div>Compare analytic vs numeric:</div>
                    <div className="text-xs text-zinc-400 mt-1">If analytic exists we overlay and compute max error below.</div>
                    {(() => {
                      if (!analyticSeries) return null;
                      let maxErr = 0;
                      for (let i = 0; i < result.t.length; i++) {
                        const num = result.C[i];
                        const ana = (k1 > 0 && k2 == 0) ? analyticSeries.Cfirst[i] : ((k2 > 0 && k1 == 0) ? analyticSeries.Csecond[i] : null);
                        if (ana !== null) {
                          const err = Math.abs(num - ana);
                          if (err > maxErr) maxErr = err;
                        }
                      }
                      return (
                        <div className="mt-2">Max absolute error vs analytic: <span className="text-zinc-100">{maxErr ? maxErr.toExponential(3) : "N/A"}</span></div>
                      );
                    })()}
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

/* 32.2 UI Component */
function PlateSection() {
  const [Nx, setNx] = useState(41);
  const [Ny, setNy] = useState(41);
  const [qMag, setQMag] = useState(1.0);
  const [run, setRun] = useState(false);
  const [result, setResult] = useState(null);

  useEffect(() => {
    if (!run) return;
    const out = solvePlateBiharmonic({
      Nx: Number(Nx),
      Ny: Number(Ny),
      qFunc: (x, y) => Number(qMag) * Math.sin(Math.PI * x) * Math.sin(Math.PI * y),
      D: 1.0,
      tol: 1e-4,
      maxIter: 2000,
    });
    setResult(out);
    setRun(false);
  }, [run]);

  return (
    <motion.div
      initial={{ opacity: 0, y: 6 }}
      whileInView={{ opacity: 1, y: 0 }}
      className={`${THEME.PANEL} p-4 rounded-2xl`}
    >
      <div className="flex items-center justify-between">
        <SectionHeader
          icon={Layers}
          number="32.2"
          title="Deflections of a Plate"
          accent={THEME.ACCENT2}
        />
        <div className="text-xs text-zinc-400">Biharmonic approximate solver (demo)</div>
      </div>

      <div className="mt-3 grid gap-4 md:grid-cols-3">
        <div className="space-y-3">
          <div>
            <label className="text-xs text-zinc-400">Nx</label>
            <Input type="number" min={21} max={121} value={Nx} onChange={(e) => setNx(Math.max(21, Math.min(121, parseInt(e.target.value || 41, 10))))} />
          </div>
          <div>
            <label className="text-xs text-zinc-400">Ny</label>
            <Input type="number" min={21} max={121} value={Ny} onChange={(e) => setNy(Math.max(21, Math.min(121, parseInt(e.target.value || 41, 10))))} />
          </div>
          <div>
            <label className="text-xs text-zinc-400">Load magnitude</label>
            <Input type="number" step="0.1" value={qMag} onChange={(e) => setQMag(e.target.value)} />
          </div>
          <div className="flex gap-2">
            <Button onClick={() => setRun(true)}>Solve Plate</Button>
            <Button onClick={() => setResult(null)}>Clear</Button>
          </div>
          <div className="text-xs text-zinc-400 mt-2">This is a demo solver illustrating shape and relative magnitudes; for production use FEM libraries.</div>
        </div>

        <div className="md:col-span-2">
          <div className="rounded border border-zinc-700 p-3 bg-zinc-900/40">
            {!result && <div className="text-sm text-zinc-400">No solution — run to compute plate deflection surface.</div>}
            {result && (
              <>
                <ResponsivePlot height={420}>
                  <Plot
                    data={[
                      {
                        z: result.Z,
                        x: result.xs,
                        y: result.ys,
                        type: "surface",
                        colorscale: "Viridis",
                        showscale: true,
                      },
                    ]}
                    layout={{
                      title: "Plate deflection surface (w)",
                      autosize: true,
                      margin: { t: 30 },
                      paper_bgcolor: "rgba(0,0,0,0)",
                      plot_bgcolor: "rgba(0,0,0,0)",
                    }}
                    style={{ width: "100%", height: 420 }}
                    useResizeHandler
                    config={{ responsive: true }}
                  />
                </ResponsivePlot>

                <div className="mt-3 grid grid-cols-1 sm:grid-cols-2 gap-3">
                  <div className="text-sm text-zinc-300">Iterations: <span className="text-zinc-100">{result.iterations}</span></div>
                  <div className="text-sm text-zinc-300">Final residual: <span className="text-zinc-100">{result.residual.toExponential(2)}</span></div>
                </div>
              </>
            )}
          </div>
        </div>
      </div>
    </motion.div>
  );
}

/* 32.3 UI Component */
function ElectrostaticsSection() {
  const [Nx, setNx] = useState(80);
  const [Ny, setNy] = useState(80);
  const [run, setRun] = useState(false);
  const [result, setResult] = useState(null);

  useEffect(() => {
    if (!run) return;
    // Example boundary: left=+1, right=-1, top/bottom 0
    const out = solveLaplaceSOR({
      Nx: Number(Nx),
      Ny: Number(Ny),
      phiBoundaryFunc: (x, y) => {
        if (x === 0) return 1;
        if (x === 1) return -1;
        return 0;
      },
      tol: 1e-5,
      maxIter: 5000,
      omega: 1.7,
    });
    setResult(out);
    setRun(false);
  }, [run]);

  return (
    <motion.div
      initial={{ opacity: 0, y: 6 }}
      whileInView={{ opacity: 1, y: 0 }}
      className={`${THEME.PANEL} p-4 rounded-2xl`}
    >
      <div className="flex items-center justify-between">
        <SectionHeader
          icon={Zap}
          number="32.3"
          title="2D Electrostatic Field Problems"
        />
        <div className="text-xs text-zinc-400">Laplace solver + field visualization</div>
      </div>

      <div className="mt-3 grid md:grid-cols-3 gap-4">
        <div className="space-y-3">
          <div>
            <label className="text-xs text-zinc-400">Nx</label>
            <Input type="number" min={20} max={200} value={Nx} onChange={(e) => setNx(Math.max(20, Math.min(200, parseInt(e.target.value || 80, 10))))} />
          </div>
          <div>
            <label className="text-xs text-zinc-400">Ny</label>
            <Input type="number" min={20} max={200} value={Ny} onChange={(e) => setNy(Math.max(20, Math.min(200, parseInt(e.target.value || 80, 10))))} />
          </div>
          <div className="flex gap-2">
            <Button onClick={() => setRun(true)}>Solve Laplace</Button>
            <Button onClick={() => setResult(null)}>Clear</Button>
          </div>
          <div className="text-xs text-zinc-400 mt-2">Boundary: left=+1V, right=-1V, top/bottom=0V — modify phiBoundaryFunc in code to change problem.</div>
        </div>

        <div className="md:col-span-2">
          <div className="rounded border border-zinc-700 p-3 bg-zinc-900/40">
            {!result && <div className="text-sm text-zinc-400">No solution — run to compute potential & field.</div>}
            {result && (
              <>
                <ResponsivePlot height={360}>
                  <Plot
                    data={[
                      {
                        z: result.phi,
                        x: linspace(0, 1, result.phi[0].length),
                        y: linspace(0, 1, result.phi.length),
                        type: "heatmap",
                        colorscale: "RdBu",
                        reversescale: true,
                        showscale: true,
                        colorbar: { title: "Φ (V)" },
                      },
                    ]}
                    layout={{
                      title: "Potential φ (V) heatmap",
                      autosize: true,
                      margin: { t: 30, b: 40, l: 40, r: 20 },
                      paper_bgcolor: "rgba(0,0,0,0)",
                      plot_bgcolor: "rgba(0,0,0,0)",
                    }}
                    style={{ width: "100%", height: 360 }}
                    useResizeHandler
                    config={{ responsive: true }}
                  />
                </ResponsivePlot>

                <div className="mt-3 grid grid-cols-1 sm:grid-cols-2 gap-3">
                  <ResponsivePlot height={320}>
                    <Plot
                      data={[
                        {
                          z: result.phi,
                          x: linspace(0, 1, result.phi[0].length),
                          y: linspace(0, 1, result.phi.length),
                          type: "contour",
                          contours: { coloring: "lines" },
                          colorscale: "Turbo",
                          name: "equipotentials",
                        },
                      ]}
                      layout={{
                        title: "Equipotential contours",
                        autosize: true,
                        margin: { t: 30 },
                        paper_bgcolor: "transparent",
                        plot_bgcolor: "transparent",
                      }}
                      style={{ width: "100%", height: 320 }}
                      useResizeHandler
                      config={{ responsive: true }}
                    />
                  </ResponsivePlot>

                  <ResponsivePlot height={320}>
                    <Plot
                      data={[
                        {
                          // quiver plot via scatter with lines (approx)
                          x: (() => {
                            const xs = [];
                            for (let j = 1; j < result.Ex.length - 1; j += Math.max(1, Math.floor(result.Ex.length / 30))) {
                              for (let i = 1; i < result.Ex[0].length - 1; i += Math.max(1, Math.floor(result.Ex[0].length / 30))) {
                                xs.push(i / (result.Ex[0].length - 1));
                              }
                            }
                            return xs;
                          })(),
                          y: (() => {
                            const ys = [];
                            for (let j = 1; j < result.Ex.length - 1; j += Math.max(1, Math.floor(result.Ex.length / 30))) {
                              for (let i = 1; i < result.Ex[0].length - 1; i += Math.max(1, Math.floor(result.Ex[0].length / 30))) {
                                ys.push(j / (result.Ex.length - 1));
                              }
                            }
                            return ys;
                          })(),
                          // u/v components flattened
                          u: (() => {
                            const us = [];
                            for (let j = 1; j < result.Ex.length - 1; j += Math.max(1, Math.floor(result.Ex.length / 30))) {
                              for (let i = 1; i < result.Ex[0].length - 1; i += Math.max(1, Math.floor(result.Ex[0].length / 30))) {
                                us.push(result.Ex[j][i]);
                              }
                            }
                            return us;
                          })(),
                          v: (() => {
                            const vs = [];
                            for (let j = 1; j < result.Ey.length - 1; j += Math.max(1, Math.floor(result.Ey.length / 30))) {
                              for (let i = 1; i < result.Ey[0].length - 1; i += Math.max(1, Math.floor(result.Ey[0].length / 30))) {
                                vs.push(result.Ey[j][i]);
                              }
                            }
                            return vs;
                          })(),
                          mode: "markers",
                          marker: { size: 2, color: "#fff" },
                          type: "scatter",
                        },
                      ]}
                      layout={{
                        title: "Electric field sample points (vector plot overlay recommended)",
                        autosize: true,
                        margin: { t: 30 },
                        paper_bgcolor: "transparent",
                        plot_bgcolor: "transparent",
                      }}
                      style={{ width: "100%", height: 320 }}
                      useResizeHandler
                      config={{ responsive: true }}
                    />
                  </ResponsivePlot>
                </div>
                <div className="mt-3 text-sm text-zinc-300">Iterations: <span className="text-zinc-100">{result.iterations}</span>, Residual: <span className="text-zinc-100">{result.residual.toExponential(3)}</span></div>
              </>
            )}
          </div>
        </div>
      </div>
    </motion.div>
  );
}

/* 32.4 UI Component */
function SpringsSection() {
  const [nSprings, setNSprings] = useState(6);
  const [kBase, setKBase] = useState(10);
  const [forceMag, setForceMag] = useState(1.0);
  const [run, setRun] = useState(false);
  const [result, setResult] = useState(null);

  useEffect(() => {
    if (!run) return;
    // build k list and f list
    const kList = new Array(Number(nSprings)).fill(Number(kBase));
    const fList = new Array(Number(nSprings) + 1).fill(0);
    // apply force at last node
    fList[fList.length - 1] = -Number(forceMag);
    const out = assembleSpringChain(kList, fList, { left: 0, right: null });
    setResult(out);
    setRun(false);
  }, [run]);

  return (
    <motion.div
      initial={{ opacity: 0, y: 6 }}
      whileInView={{ opacity: 1, y: 0 }}
      className={`${THEME.PANEL} p-4 rounded-2xl`}
    >
      <div className="flex items-center justify-between">
        <SectionHeader icon={Box} number="32.4" title="Series of Springs (static)" />
        <div className="text-xs text-zinc-400">Stiffness assembly & solve</div>
      </div>

      <div className="mt-3 grid md:grid-cols-3 gap-4">
        <div className="space-y-3">
          <div>
            <label className="text-xs text-zinc-400">Number of springs</label>
            <Input type="number" min={1} max={40} value={nSprings} onChange={(e) => setNSprings(Math.max(1, Math.min(40, parseInt(e.target.value || 6, 10))))} />
          </div>
          <div>
            <label className="text-xs text-zinc-400">Spring stiffness k (uniform)</label>
            <Input type="number" step="0.1" value={kBase} onChange={(e) => setKBase(e.target.value)} />
          </div>
          <div>
            <label className="text-xs text-zinc-400">Applied force at end (negative pulls)</label>
            <Input type="number" step="0.1" value={forceMag} onChange={(e) => setForceMag(e.target.value)} />
          </div>
          <div className="flex gap-2">
            <Button onClick={() => setRun(true)}>Assemble & Solve</Button>
            <Button onClick={() => setResult(null)}>Clear</Button>
          </div>
        </div>

        <div className="md:col-span-2">
          <div className="rounded border border-zinc-700 p-3 bg-zinc-900/40">
            {!result && <div className="text-sm text-zinc-400">No assembly yet — run to build stiffness matrix and solve for nodal displacements.</div>}
            {result && (
              <>
                <div className="text-sm text-zinc-300 mb-2">Nodal displacements (u):</div>
                <ResponsivePlot height={300}>
                  <Plot
                    data={[
                      {
                        x: result.U.map((_, i) => i),
                        y: result.U,
                        mode: "lines+markers",
                        name: "u (displacements)",
                        line: { color: THEME.ACCENT },
                      },
                    ]}
                    layout={{
                      title: "Nodal displacement along spring chain",
                      autosize: true,
                      margin: { t: 30, b: 40 },
                      paper_bgcolor: "transparent",
                      plot_bgcolor: "transparent",
                      font: { color: "#cbd5e1" },
                    }}
                    style={{ width: "100%", height: 300 }}
                    useResizeHandler
                    config={{ responsive: true }}
                  />
                </ResponsivePlot>

                <div className="mt-3 text-sm text-zinc-300">Stiffness matrix (K) size: <span className="text-zinc-100">{result.K.length}×{result.K.length}</span></div>
              </>
            )}
          </div>
        </div>
      </div>
    </motion.div>
  );
}

/* Problems Section */
function Chapter32Problems() {
  const [open, setOpen] = useState(-1);
  const problems = [
    {
      id: 1,
      title: "Compare RK4 vs analytic for first-order kinetics",
      prompt: "Use the analytic C(t)=C0 e^{-k1 t} and compare to RK4 numeric across dt. Report max error vs dt and observe O(dt^4) local error behavior.",
      solution: "Run simulateCSTR with k2=0 and vary dt, plot max error vs dt on log-log; slope ~4 for RK4 integrated error regime until round-off/consistency limits.",
    },
    {
      id: 2,
      title: "Plate deflection convergence",
      prompt: "Refine Nx,Ny and compute central deflection for the plate. Observe convergence trend and discuss numerical damping due to relaxation method.",
      solution: "Use doubling grids and observe approximate convergence; for rigorous results use FEM.",
    },
    {
      id: 3,
      title: "Electrostatics: field line tracing",
      prompt: "Trace streamlines of E field from a set of seed points using RK4 on vector field computed from potential.",
      solution: "Interpolate Ex,Ey at arbitrary points and integrate dx/ds = Ex, dy/ds = Ey; implement simple RK4 for pathlines.",
    },
    {
      id: 4,
      title: "Spring chain with variable stiffness",
      prompt: "Set k_i varying (e.g., geometric series) and examine where most deformation localizes.",
      solution: "Assemble K with kList and inspect nodal u; softer springs show larger local extension.",
    },
  ];

  return (
    <motion.div
      initial={{ opacity: 0, y: 6 }}
      whileInView={{ opacity: 1, y: 0 }}
      className={`${THEME.PANEL} p-4 rounded-2xl`}
    >
      <div className="flex items-center justify-between">
        <SectionHeader icon={ListChecks} number="Problems" title="Exercises & Projects" />
        <div className="text-xs text-zinc-400">Ideas for exploration and extensions</div>
      </div>

      <div className="mt-3 space-y-3">
        {problems.map((p, idx) => (
          <div key={p.id} className="rounded border border-zinc-700 p-3 bg-zinc-900/40">
            <div className="flex items-center justify-between">
              <div>
                <div className="text-sm text-zinc-100 font-medium">{p.title}</div>
                <div className="text-xs text-zinc-400">{p.prompt}</div>
              </div>
              <div>
                <Button variant="ghost" onClick={() => setOpen(open === idx ? -1 : idx)}>{open === idx ? "Hide" : "Show"}</Button>
              </div>
            </div>
            {open === idx && (
              <div className="mt-3 text-sm text-zinc-200 bg-zinc-950/30 border border-zinc-800 p-3 rounded">
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
// --- Utility: Simple dense solver for Ku = F ---
function solveDense(K, F) {
  const n = K.length;
  const A = K.map(row => row.slice()); // copy of K
  const b = F.slice();                 // copy of F

  // Forward elimination
  for (let i = 0; i < n; i++) {
    // Pivot if diagonal is zero
    if (Math.abs(A[i][i]) < 1e-12) {
      for (let k = i + 1; k < n; k++) {
        if (Math.abs(A[k][i]) > Math.abs(A[i][i])) {
          [A[i], A[k]] = [A[k], A[i]];
          [b[i], b[k]] = [b[k], b[i]];
          break;
        }
      }
    }

    // Normalize pivot row
    const pivot = A[i][i];
    for (let j = i; j < n; j++) {
      A[i][j] /= pivot;
    }
    b[i] /= pivot;

    // Eliminate below
    for (let k = i + 1; k < n; k++) {
      const factor = A[k][i];
      for (let j = i; j < n; j++) {
        A[k][j] -= factor * A[i][j];
      }
      b[k] -= factor * b[i];
    }
  }

  // Back substitution
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    x[i] = b[i];
    for (let j = i + 1; j < n; j++) {
      x[i] -= A[i][j] * x[j];
    }
  }

  return x;
}


/* ------------------------------------------------------------------ */
/* Chapter 32 assembly (page)                                         */
/* ------------------------------------------------------------------ */

export default function Chapter32() {
  useEffect(() => {
    if (typeof window !== "undefined") window.scrollTo({ top: 0, behavior: "smooth" });
  }, []);

  return (
    <div className={`p-4 sm:p-6 lg:p-8 space-y-6 ${THEME.BG} min-h-screen`}>
      <motion.div initial={{ opacity: 0, y: -8 }} animate={{ opacity: 1, y: 0 }}>
        <div className="flex flex-col sm:flex-row sm:items-center sm:justify-between gap-4">
          <div>

            <h1 className="text-2xl sm:text-3xl font-bold text-emerald-400">Case Studies: Partial Differential Equations</h1>
            <div className="text-zinc-400 mt-1 text-sm">32.1 Reactor • 32.2 Plate Deflection • 32.3 Electrostatics • 32.4 Springs • Problems</div>
          </div>

          
        </div>
      </motion.div>

      <div className="space-y-6">
        <ReactorSection />
        <PlateSection />
        <ElectrostaticsSection />
        <SpringsSection />
        <Chapter32Problems />
        <BottomBar/>
      </div>

      <footer className="text-xs text-zinc-500 mt-6">
        Chapter 32 — Advanced PDE case studies with responsive visualizations. For production-quality FEM/PDE work use specialized libraries (FEniCS, PETSc, deal.II).
      </footer>
    </div>
  );
}
