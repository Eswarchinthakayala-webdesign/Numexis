// src/pages/Chapter25.jsx
// Chapter 25 — Runge-Kutta Methods (expanded, interactive, KaTeX-ready)
// Diverse examples: Euler, Improved Euler, RK4, Systems, Adaptive RK
// Responsive Plotly charts wrapped safely to prevent overflow

import React, { useMemo, useState } from "react";
import Plot from "react-plotly.js";
import { motion } from "framer-motion";
import ReactMarkdown from "react-markdown";
import remarkMath from "remark-math";
import rehypeKatex from "rehype-katex";
import "katex/dist/katex.min.css";

import { Card, CardHeader, CardContent, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Button } from "@/components/ui/button";

import {
  Calculator,
  Zap,
  Activity,
  Layers,
  Compass,
  BookOpen,
  List,
  ChevronRight,
  ChevronDown,
  RefreshCcw,
  FileText,
} from "lucide-react";
import BottomBar from "../components/BottomBar";

/* ===========================
   Theme + small helpers
   =========================== */

const theme = {
  bg: "bg-zinc-950",
  panel: "bg-zinc-900/60 border border-zinc-700 overflow-hidden",
  accent: "#22d3ee",
  accent2: "#34d399",
  warn: "#f59e0b",
};

const fadeUp = {
  initial: { opacity: 0, y: 10 },
  animate: { opacity: 1, y: 0 },
  transition: { duration: 0.28 },
};

/* ===========================
   Numeric helpers
   =========================== */

function linspace(a, b, n) {
  if (n <= 1) return [a];
  const dx = (b - a) / (n - 1);
  const out = new Array(n);
  for (let i = 0; i < n; i++) out[i] = a + i * dx;
  return out;
}

function trapezoidIntegrate(x, y) {
  let s = 0;
  for (let i = 0; i < x.length - 1; i++) {
    const dx = x[i + 1] - x[i];
    s += 0.5 * (y[i] + y[i + 1]) * dx;
  }
  return s;
}

function simpsonIntegrate(x, y) {
  const n = x.length;
  if (n < 3) return trapezoidIntegrate(x, y);
  const h = (x[n - 1] - x[0]) / (n - 1);
  let m = n;
  if ((m - 1) % 2 !== 0) m = n - 1; // make even number of intervals
  let s = y[0] + y[m - 1];
  for (let i = 1; i < m - 1; i++) {
    s += (i % 2 === 0 ? 2 : 4) * y[i];
  }
  s = (h / 3) * s;
  if (m !== n) {
    // handle last interval by trapezoid
    s += 0.5 * (y[m - 1] + y[m]) * (x[m] - x[m - 1]);
  }
  return s;
}

/* Small linear solver used by some filters (naive Gauss elimination) */
function solveLinearSystem(A, b) {
  const n = A.length;
  const M = A.map((row) => row.slice());
  const rhs = b.slice();
  for (let k = 0; k < n; k++) {
    let piv = k;
    for (let i = k + 1; i < n; i++) if (Math.abs(M[i][k]) > Math.abs(M[piv][k])) piv = i;
    if (Math.abs(M[piv][k]) < 1e-15) return null;
    [M[k], M[piv]] = [M[piv], M[k]];
    [rhs[k], rhs[piv]] = [rhs[piv], rhs[k]];
    for (let i = k + 1; i < n; i++) {
      const f = M[i][k] / M[k][k];
      for (let j = k; j < n; j++) M[i][j] -= f * M[k][j];
      rhs[i] -= f * rhs[k];
    }
  }
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    let s = rhs[i];
    for (let j = i + 1; j < n; j++) s -= M[i][j] * x[j];
    x[i] = s / M[i][i];
  }
  return x;
}

/* ===========================
   Reusable UI components
   =========================== */

function Formula({ children }) {
  return (
    <div className="my-2">
      <div className="px-3 py-2 rounded-lg text-gray-300 border border-zinc-700 bg-zinc-900/80 overflow-x-auto">
        <div className="text-xs uppercase tracking-widest mb-1" style={{ color: theme.accent }}>
          Formula
        </div>
        <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>
          {children}
        </ReactMarkdown>
      </div>
    </div>
  );
}

function Note({ children }) {
  return (
    <div className="p-3 bg-zinc-800/40 border-l-4 border-emerald-500 rounded-md text-zinc-300 text-sm">
      {children}
    </div>
  );
}

function Labeled({ label, children }) {
  return (
    <div>
      <div className="text-xs text-zinc-400 mb-1">{label}</div>
      {children}
    </div>
  );
}

function Summary({ items }) {
  return (
    <div className="grid grid-cols-1 sm:grid-cols-3 gap-3">
      {items.map((it, i) => (
        <div key={i} className="rounded-lg bg-zinc-900/70 border border-zinc-700 px-3 py-2">
          <div className="text-xs text-zinc-400">{it.label}</div>
          <div className="text-zinc-100 text-sm">{it.value}</div>
        </div>
      ))}
    </div>
  );
}

function CollapsibleExercise({ title, body, solution }) {
  const [open, setOpen] = useState(false);
  return (
    <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-3">
      <div className="flex items-center justify-between">
        <div className="text-zinc-100 font-medium">{title}</div>
        <Button variant="ghost" className="bg-white cursor-pointer hover:bg-gray-200" size="sm" onClick={() => setOpen((o) => !o)}>
          {open ? <ChevronDown className="w-4 h-4" /> : <ChevronRight className="w-4 h-4" />} 
        </Button>
      </div>
      <div className="text-zinc-300 mt-2">
        <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>
          {body}
        </ReactMarkdown>
      </div>
      {open && (
        <div className="mt-3 border-t border-zinc-800 pt-3 text-zinc-200">
          <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>
            {solution}
          </ReactMarkdown>
        </div>
      )}
    </div>
  );
}

/* ===========================
   ODE solvers (vector-friendly)
   =========================== */

function eulerSolver(f, x0, y0, h, nSteps) {
  const xs = [x0];
  const ys = [Array.isArray(y0) ? y0.slice() : y0];
  let x = x0;
  let y = Array.isArray(y0) ? y0.slice() : y0;
  for (let i = 0; i < nSteps; i++) {
    const dy = f(x, y);
    if (Array.isArray(y)) {
      y = y.map((val, idx) => val + h * dy[idx]);
    } else {
      y = y + h * dy;
    }
    x += h;
    xs.push(x);
    ys.push(Array.isArray(y) ? y.slice() : y);
  }
  return { xs, ys };
}

function rk4Solver(f, x0, y0, h, nSteps) {
  const xs = [x0];
  const ys = [Array.isArray(y0) ? y0.slice() : y0];
  let x = x0;
  let y = Array.isArray(y0) ? y0.slice() : y0;
  for (let i = 0; i < nSteps; i++) {
    const add = (a, b) => a.map((v, j) => v + b[j]);
    const scale = (a, s) => a.map((v) => v * s);
    if (Array.isArray(y)) {
      const k1 = f(x, y);
      const yk2 = add(y, scale(k1, h / 2));
      const k2 = f(x + h / 2, yk2);
      const yk3 = add(y, scale(k2, h / 2));
      const k3 = f(x + h / 2, yk3);
      const yk4 = add(y, scale(k3, h));
      const k4 = f(x + h, yk4);
      y = y.map((val, j) => val + (h / 6) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]));
    } else {
      const k1 = f(x, y);
      const k2 = f(x + h / 2, y + (h / 2) * k1);
      const k3 = f(x + h / 2, y + (h / 2) * k2);
      const k4 = f(x + h, y + h * k3);
      y = y + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
    }
    x += h;
    xs.push(x);
    ys.push(Array.isArray(y) ? y.slice() : y);
  }
  return { xs, ys };
}

/* Dormand-Prince / RKF45-style adaptive stepper (simple) */
function rkf45Adaptive(f, x0, y0, h0, xEnd, tol = 1e-6, hMin = 1e-8, hMax = 1.0) {
  // coefficients for Dormand-Prince (Butcher tableau)
  function add(a, b) {
    return a.map((v, i) => v + b[i]);
  }
  function addScaled(a, b, s) {
    return a.map((v, i) => v + s * b[i]);
  }
  function scale(a, s) {
    return a.map((v) => v * s);
  }
  let x = x0;
  let y = Array.isArray(y0) ? y0.slice() : y0;
  let h = h0;
  const xs = [x];
  const ys = [Array.isArray(y0) ? y0.slice() : y0];

  // coefficients (Dormand-Prince)
  const c = [0, 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1, 1];
  const a = [
    [],
    [1 / 5],
    [3 / 40, 9 / 40],
    [44 / 45, -56 / 15, 32 / 9],
    [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729],
    [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656],
    [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84],
  ];
  const b = [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84, 0];
  const bStar = [
    5179 / 57600,
    0,
    7571 / 16695,
    393 / 640,
    -92097 / 339200,
    187 / 2100,
    1 / 40,
  ];

  while (x < xEnd - 1e-12) {
    if (h < hMin) h = hMin;
    if (h > hMax) h = hMax;
    if (x + h > xEnd) h = xEnd - x;

    const ks = [];
    // compute k1..k7
    for (let stage = 0; stage < 7; stage++) {
      if (stage === 0) {
        ks[0] = f(x, y);
      } else {
        let ystage = y.slice ? y.slice() : y;
        for (let j = 0; j < stage; j++) {
          if (Array.isArray(y)) {
            ystage = add(ystage, scale(ks[j], a[stage][j] * h));
          } else {
            ystage = ystage + a[stage][j] * h * ks[j];
          }
        }
        ks[stage] = f(x + c[stage] * h, ystage);
      }
    }
    // compute 4th and 5th order estimates
    let y5, y4;
    if (Array.isArray(y)) {
      const y5temp = y.slice();
      const y4temp = y.slice();
      for (let j = 0; j < 7; j++) {
        y5temp.forEach((_, idx) => {
          y5temp[idx] = y5temp[idx] + h * b[j] * ks[j][idx];
        });
        y4temp.forEach((_, idx) => {
          y4temp[idx] = y4temp[idx] + h * bStar[j] * ks[j][idx];
        });
      }
      y5 = y5temp;
      y4 = y4temp;
    } else {
      let y5temp = y;
      let y4temp = y;
      for (let j = 0; j < 7; j++) {
        y5temp = y5temp + h * b[j] * ks[j];
        y4temp = y4temp + h * bStar[j] * ks[j];
      }
      y5 = y5temp;
      y4 = y4temp;
    }

    // error estimate
    let err;
    if (Array.isArray(y)) {
      let s = 0;
      for (let i = 0; i < y.length; i++) s += Math.pow(y5[i] - y4[i], 2);
      err = Math.sqrt(s / y.length);
    } else {
      err = Math.abs(y5 - y4);
    }

    // adapt step
    const safety = 0.9;
    const order = 5;
    if (err <= tol) {
      // accept
      x += h;
      y = Array.isArray(y5) ? y5.slice() : y5;
      xs.push(x);
      ys.push(Array.isArray(y) ? y.slice() : y);
      // increase
      let hNew = h * Math.min(5, Math.max(0.2, safety * Math.pow(tol / (err + 1e-15), 1 / (order + 1))));
      h = Math.min(hNew, hMax);
    } else {
      // reject and reduce
      let hNew = h * Math.max(0.1, safety * Math.pow(tol / (err + 1e-15), 1 / (order + 1)));
      h = Math.max(hNew, hMin);
    }
    if (!isFinite(h) || h <= 0) break;
  }

  return { xs, ys };
}

/* ===========================
   Section 25.1 — Euler's Method
   Example: dy/dx = -2y + x, y(0)=1 (analytical solution provided)
   =========================== */

function Section251() {
  // example ODE and analytic solution
  function f(x, y) {
    return -2 * y + x;
  }
  function exact(x) {
    // exact solution for dy/dx = -2y + x, y(0)=1
    // solution computed analytically for demonstration
    // solve linear ODE y' + 2y = x -> integrating factor e^{2x}
    // y*e^{2x} = ∫ x e^{2x} dx + C; ∫ x e^{2x} dx = (x - 1/2) e^{2x} / 2
    // after simplification we pick constant to match y(0)=1
    return 0.25 * x - 0.125 + 1.125 * Math.exp(-2 * x);
  }

  const [h, setH] = useState(0.2);
  const [x0] = useState(0);
  const [xEnd] = useState(2.0);
  const n = Math.max(2, Math.floor((xEnd - x0) / h));
  const xs = linspace(x0, xEnd, n + 1);
  const euler = useMemo(() => eulerSolver(f, x0, 1.0, h, n).ys.map((y) => Array.isArray(y) ? y[0] : y), [h, n]);
  const rk4 = useMemo(() => rk4Solver(f, x0, 1.0, h, n).ys.map((y) => Array.isArray(y) ? y[0] : y), [h, n]);
  const exactVals = useMemo(() => xs.map((xx) => exact(xx)), [xs]);

  // error vs step sizes
  const hs = [0.8, 0.4, 0.2, 0.1, 0.05, 0.025];
  const errs = hs.map((hh) => {
    const m = Math.floor((xEnd - x0) / hh);
    const e = eulerSolver(f, x0, 1.0, hh, m).ys;
    const last = e[e.length - 1];
    const approx = Array.isArray(last) ? last[0] : last;
    return { h: hh, err: Math.abs(approx - exact(xEnd)) };
  });
  const errsRK = hs.map((hh) => {
    const m = Math.floor((xEnd - x0) / hh);
    const r = rk4Solver(f, x0, 1.0, hh, m).ys;
    const last = r[r.length - 1];
    const approx = Array.isArray(last) ? last[0] : last;
    return { h: hh, err: Math.abs(approx - exact(xEnd)) };
  });

  return (
    <motion.div {...fadeUp}>
      <Card className={`${theme.panel} p-4 rounded-2xl`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Zap className="w-5 h-5" /> 25.1 Euler’s Method
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="grid md:grid-cols-2 gap-4">
            <div>
              <Formula>{`**Euler step:**\n\n$$y_{n+1} = y_n + h f(x_n,y_n).$$\n\nLocal truncation error is $\\mathcal{O}(h^2)$ and global error $\\mathcal{O}(h)$.`}</Formula>
              <Note>Euler is explicit and simple — good to demonstrate behavior but inefficient for high accuracy.</Note>
            </div>
            <div>
              <Labeled label="Step size h">
                <Input value={h} onChange={(e) => setH(Number(e.target.value))} className="bg-zinc-800 text-white" />
              </Labeled>
              <Summary items={[
                { label: "Interval", value: `x∈[${x0},${xEnd}]` },
                { label: "Steps (approx)", value: n },
                { label: "Initial y(0)", value: 1.0 },
              ]} />
            </div>
          </div>

          <div className="grid md:grid-cols-3 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: xs, y: exactVals, type: "scatter", mode: "lines", name: "Exact", line: { color: theme.accent2 } },
                    { x: xs, y: euler, type: "scatter", mode: "lines+markers", name: "Euler", line: { color: theme.accent } },
                    { x: xs, y: rk4, type: "scatter", mode: "lines+markers", name: "RK4(h)", line: { color: theme.warn } },
                  ]}
                  layout={{
                    title: "Euler vs Exact vs RK4",
                    autosize: true,
                    responsive: true,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    margin: { t: 36, l: 48, r: 20, b: 36 },
                    font: { color: "#e5e7eb" },
                  }}
                  useResizeHandler
                  style={{ width: "100%", height: "100%" }}
                  config={{ responsive: true, displayModeBar: false }}
                />
              </div>
            </div>

            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: hs.map((d) => d), y: errs.map((e) => e.err), type: "scatter", mode: "lines+markers", name: "Euler error", line: { color: theme.accent } },
                    { x: hs.map((d) => d), y: errsRK.map((e) => e.err), type: "scatter", mode: "lines+markers", name: "RK4 error", line: { color: theme.accent2 } },
                  ]}
                  layout={{
                    title: "Error at x=2 vs h (log-log)",
                    xaxis: { type: "log", title: "h" },
                    yaxis: { type: "log", title: "abs error" },
                    autosize: true,
                    responsive: true,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    margin: { t: 36, l: 48, r: 20, b: 36 },
                    font: { color: "#e5e7eb" },
                  }}
                  useResizeHandler
                  style={{ width: "100%", height: "100%" }}
                  config={{ responsive: true, displayModeBar: false }}
                />
              </div>
            </div>

            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: xs, y: exactVals.map((v, i) => Math.abs(v - euler[i])), type: "scatter", mode: "lines", name: "abs error (Euler)", line: { color: theme.accent } },
                    { x: xs, y: exactVals.map((v, i) => Math.abs(v - rk4[i])), type: "scatter", mode: "lines", name: "abs error (RK4)", line: { color: theme.accent2 } },
                  ]}
                  layout={{
                    title: "Error vs x (local/global)",
                    autosize: true,
                    responsive: true,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    margin: { t: 36, l: 48, r: 20, b: 36 },
                    font: { color: "#e5e7eb" },
                  }}
                  useResizeHandler
                  style={{ width: "100%", height: "100%" }}
                  config={{ responsive: true, displayModeBar: false }}
                />
              </div>
            </div>
          </div>

          <CollapsibleExercise
            title="Exercise 25.1.1 — Step halving"
            body={`Use Euler's method with h and h/2 to approximate y(2). Use Richardson extrapolation to estimate error and improved value. Explain why extrapolation improves accuracy.`}
            solution={`Richardson cancels leading-order error (global O(h)) when combining approximations, improving to a higher effective order.`}
          />
        </CardContent>
      </Card>
    </motion.div>
  );
}

/* ===========================
   Section 25.2 — Improvements of Euler (Heun, Midpoint)
   Example: dy/dx = y - x^2 + 1, y(0)=0.5 (known analytic)
   =========================== */

function Section252() {
  // Example ODE with known exact solution y(x) = (x+1)^2 - 0.5*e^{x}
  // We'll use a test function with known behavior (for demonstration).
  function f(x, y) {
    return y - x * x + 1;
  }
  function exact(x) {
    // Not the exact solution of that ODE unless specified; we'll use numeric reference
    // We'll build a dense numeric reference using RK4 with very small h to serve as "exact".
    return null;
  }

  const [h, setH] = useState(0.2);
  const x0 = 0;
  const xEnd = 2.0;

  const n = Math.max(2, Math.floor((xEnd - x0) / h));
  const xs = linspace(x0, xEnd, n + 1);

  // Methods: Euler, Heun (Improved Euler), Midpoint, RK4
  function heunSolver(f, x0, y0, h, nSteps) {
    const xs = [x0];
    const ys = [y0];
    let x = x0;
    let y = y0;
    for (let i = 0; i < nSteps; i++) {
      const k1 = f(x, y);
      const yPred = y + h * k1;
      const k2 = f(x + h, yPred);
      y = y + (h / 2) * (k1 + k2);
      x += h;
      xs.push(x);
      ys.push(y);
    }
    return { xs, ys };
  }

  function midpointSolver(f, x0, y0, h, nSteps) {
    const xs = [x0];
    const ys = [y0];
    let x = x0;
    let y = y0;
    for (let i = 0; i < nSteps; i++) {
      const k1 = f(x, y);
      const yMid = y + (h / 2) * k1;
      const k2 = f(x + h / 2, yMid);
      y = y + h * k2;
      x += h;
      xs.push(x);
      ys.push(y);
    }
    return { xs, ys };
  }

  // dense RK4 as reference
  const dense = useMemo(() => {
    const Nref = 4000;
    return rk4Solver(f, x0, 0.5, (xEnd - x0) / Nref, Nref).ys;
  }, []);

  const euler = useMemo(() => eulerSolver(f, x0, 0.5, h, n).ys, [h, n]);
  const heun = useMemo(() => heunSolver(f, x0, 0.5, h, n).ys, [h, n]);
  const mid = useMemo(() => midpointSolver(f, x0, 0.5, h, n).ys, [h, n]);
  const rk4vals = useMemo(() => rk4Solver(f, x0, 0.5, h, n).ys, [h, n]);

  // build reference values aligned with xs by interpolating dense
  const denseXs = useMemo(() => linspace(x0, xEnd, dense.length), [dense.length]);
  const denseInterp = (xq) => {
    // linear interpolation onto dense grid (simple)
    const N = denseXs.length;
    if (xq <= denseXs[0]) return dense[0];
    if (xq >= denseXs[N - 1]) return dense[N - 1];
    let i = 0;
    while (i < N - 1 && denseXs[i + 1] < xq) i++;
    const xL = denseXs[i], xR = denseXs[i + 1];
    const yL = dense[i], yR = dense[i + 1];
    const t = (xq - xL) / (xR - xL);
    if (Array.isArray(yL)) return yL[0] * (1 - t) + yR[0] * t;
    return yL * (1 - t) + yR * t;
  };

  const exactVals = xs.map((xx) => denseInterp(xx));

  // errors at xEnd for convergence
  const hs = [0.8, 0.4, 0.2, 0.1, 0.05];
  const conv = hs.map((hh) => {
    const m = Math.floor((xEnd - x0) / hh);
    const e = heunSolver(f, x0, 0.5, hh, m).ys;
    const last = e[e.length - 1];
    const approx = Array.isArray(last) ? last[0] : last;
    // reference from dense
    const ref = dense[dense.length - 1];
    const refVal = Array.isArray(ref) ? ref[0] : ref;
    return { h: hh, err: Math.abs(approx - refVal) };
  });

  return (
    <motion.div {...fadeUp}>
      <Card className={`${theme.panel} p-4 rounded-2xl`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Activity className="w-5 h-5" /> 25.2 Improvements of Euler’s Method (Heun / Midpoint)
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid md:grid-cols-2 gap-4">
            <div>
              <Formula>{`**Heun (improved Euler)**\n\n$$y_{n+1} = y_n + \\frac{h}{2}\\left[f(x_n,y_n)+f(x_{n+1}, y_n + h f(x_n,y_n))\\right]$$\n\n**Midpoint method** (second-order):\n\n$$k_1=f(x_n,y_n),\\quad k_2=f\\left(x_n+\\frac{h}{2}, y_n+\\frac{h}{2}k_1\\right),\\quad y_{n+1}=y_n+h k_2.$$`}</Formula>
              <Note>Both Heun and midpoint achieve second-order accuracy (global error O(h^2)) while staying explicit.</Note>
            </div>

            <div>
              <Labeled label="Step size h"><Input value={h} onChange={(e) => setH(Number(e.target.value))} className="bg-zinc-800 text-white" /></Labeled>
              <Summary items={[
                { label: "Method types", value: "Euler / Heun / Midpoint / RK4" },
                { label: "Initial y0", value: 0.5 },
                { label: "Reference", value: "Dense RK4" },
              ]} />
            </div>
          </div>

          <div className="grid md:grid-cols-3 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: xs, y: exactVals, type: "scatter", mode: "lines", name: "dense RK4 (ref)", line: { color: theme.accent2 } },
                    { x: xs, y: euler.map((v) => Array.isArray(v) ? v[0] : v), type: "scatter", mode: "lines+markers", name: "Euler", line: { color: theme.warn } },
                    { x: xs, y: heun.map((v) => Array.isArray(v) ? v[0] : v), type: "scatter", mode: "lines+markers", name: "Heun", line: { color: theme.accent } },
                    { x: xs, y: mid.map((v) => Array.isArray(v) ? v[0] : v), type: "scatter", mode: "lines+markers", name: "Midpoint", line: { color: theme.accent2 } },
                  ]}
                  layout={{
                    title: "Method comparison vs dense RK4",
                    autosize: true,
                    responsive: true,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    margin: { t: 36, l: 48, r: 20, b: 36 },
                    font: { color: "#e5e7eb" },
                  }}
                  useResizeHandler
                  style={{ width: "100%", height: "100%" }}
                  config={{ responsive: true, displayModeBar: false }}
                />
              </div>
            </div>

            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: xs, y: exactVals.map((v, i) => Math.abs(v - (euler[i] ? (Array.isArray(euler[i]) ? euler[i][0] : euler[i]) : 0))), type: "scatter", mode: "lines", name: "Euler err", line: { color: theme.warn } },
                    { x: xs, y: exactVals.map((v, i) => Math.abs(v - (heun[i] ? (Array.isArray(heun[i]) ? heun[i][0] : heun[i]) : 0))), type: "scatter", mode: "lines", name: "Heun err", line: { color: theme.accent } },
                  ]}
                  layout={{
                    title: "Pointwise Absolute Errors",
                    autosize: true,
                    responsive: true,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    margin: { t: 36, l: 48, r: 20, b: 36 },
                    font: { color: "#e5e7eb" },
                  }}
                  useResizeHandler
                  style={{ width: "100%", height: "100%" }}
                  config={{ responsive: true, displayModeBar: false }}
                />
              </div>
            </div>

            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: conv.map((c) => c.h), y: conv.map((c) => c.err), type: "scatter", mode: "lines+markers", name: "Heun error vs h", line: { color: theme.accent } },
                  ]}
                  layout={{
                    title: "Convergence (Heun) at x=2",
                    xaxis: { title: "h", type: "log" },
                    yaxis: { title: "abs error", type: "log" },
                    autosize: true,
                    responsive: true,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    margin: { t: 36, l: 48, r: 20, b: 36 },
                    font: { color: "#e5e7eb" },
                  }}
                  useResizeHandler
                  style={{ width: "100%", height: "100%" }}
                  config={{ responsive: true, displayModeBar: false }}
                />
              </div>
            </div>
          </div>

          <CollapsibleExercise
            title="Exercise 25.2.1 — Constructive comparison"
            body={`For the given ODE, compare global error of Euler, Heun and Midpoint for h = 0.2, 0.1, 0.05. Which is cheapest for a target error of 1e-3?`}
            solution={`Heun and Midpoint both reach O(h^2) global error — typically they are preferable to use over Euler for modest extra cost.`}
          />
        </CardContent>
      </Card>
    </motion.div>
  );
}

/* ===========================
   Section 25.3 — RK4 + error vs h
   Example: dy/dx = y * cos(x), y(0)=1 (nonlinear but simple)
   =========================== */

function Section253() {
  function f(x, y) {
    return y * Math.cos(x);
  }
  // analytic reference via dense RK4 or known solution (no simple closed form), so use dense RK4
  const [h, setH] = useState(0.1);
  const x0 = 0;
  const xEnd = 6.28; // ~2π
  const n = Math.max(2, Math.floor((xEnd - x0) / h));
  const xs = linspace(x0, xEnd, n + 1);
  const rk4vals = useMemo(() => rk4Solver(f, x0, 1.0, h, n).ys.map((v) => Array.isArray(v) ? v[0] : v), [h, n]);

  // dense reference:
  const dense = useMemo(() => {
    const Nref = 8000;
    return rk4Solver(f, x0, 1.0, (xEnd - x0) / Nref, Nref).ys.map((v) => Array.isArray(v) ? v[0] : v);
  }, []);

  // align dense with xs by interpolation (linear)
  const denseXs = useMemo(() => linspace(x0, xEnd, dense.length), [dense.length]);
  function interpDense(xq) {
    if (xq <= denseXs[0]) return dense[0];
    if (xq >= denseXs[denseXs.length - 1]) return dense[dense.length - 1];
    let i = 0;
    while (i < denseXs.length - 1 && denseXs[i + 1] < xq) i++;
    const xL = denseXs[i], xR = denseXs[i + 1];
    const yL = dense[i], yR = dense[i + 1];
    const t = (xq - xL) / (xR - xL);
    return yL * (1 - t) + yR * t;
  }
  const exactVals = xs.map((xx) => interpDense(xx));

  const hs = [0.5, 0.25, 0.125, 0.0625, 0.03125];
  const conv = hs.map((hh) => {
    const m = Math.floor((xEnd - x0) / hh);
    const arr = rk4Solver(f, x0, 1.0, hh, m).ys;
    const last = arr[arr.length - 1];
    const approx = Array.isArray(last) ? last[0] : last;
    const ref = dense[dense.length - 1];
    const refVal = Array.isArray(ref) ? ref[0] : ref;
    return { h: hh, err: Math.abs(approx - refVal) };
  });

  return (
    <motion.div {...fadeUp}>
      <Card className={`${theme.panel} p-4 rounded-2xl`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Calculator className="w-5 h-5" /> 25.3 Runge-Kutta Methods (RK4)
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="grid md:grid-cols-2 gap-4">
            <div>
             <Formula>{`
RK4 step (explicit):

$$
\\begin{align*}
k_1 &= f(x_n, y_n), \\\\
k_2 &= f\\left(x_n + \\frac{h}{2}, y_n + \\frac{h}{2} k_1\\right), \\\\
k_3 &= f\\left(x_n + \\frac{h}{2}, y_n + \\frac{h}{2} k_2\\right), \\\\
k_4 &= f(x_n + h, y_n + h k_3), \\\\
y_{n+1} &= y_n + \\frac{h}{6} (k_1 + 2 k_2 + 2 k_3 + k_4).
\\end{align*}
$$
`}</Formula>

              <Note>RK4 commonly gives global O(h^4) error for smooth problems and is a reliable general-purpose integrator.</Note>
            </div>

            <div>
              <Labeled label="Step size h"><Input value={h} onChange={(e) => setH(Number(e.target.value))} className="bg-zinc-800 text-white" /></Labeled>
              <Summary items={[
                { label: "Interval", value: `[${x0}, ${xEnd}]` },
                { label: "Steps", value: n },
                { label: "Method", value: "RK4" },
              ]} />
            </div>
          </div>

          <div className="grid md:grid-cols-3 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: xs, y: exactVals, type: "scatter", mode: "lines", name: "dense (ref)", line: { color: theme.accent2 } },
                    { x: xs, y: rk4vals, type: "scatter", mode: "lines+markers", name: "RK4 (h)", line: { color: theme.accent } },
                  ]}
                  layout={{
                    title: "RK4 solution vs dense reference",
                    autosize: true,
                    responsive: true,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    margin: { t: 36, l: 48, r: 20, b: 36 },
                    font: { color: "#e5e7eb" },
                  }}
                  useResizeHandler
                  style={{ width: "100%", height: "100%" }}
                  config={{ responsive: true, displayModeBar: false }}
                />
              </div>
            </div>

            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: xs, y: exactVals.map((v, i) => Math.abs(v - rk4vals[i])), type: "scatter", mode: "lines", name: "abs error", line: { color: theme.warn } },
                  ]}
                  layout={{
                    title: "Pointwise absolute error (RK4)",
                    autosize: true,
                    responsive: true,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    margin: { t: 36, l: 48, r: 20, b: 36 },
                    font: { color: "#e5e7eb" },
                  }}
                  useResizeHandler
                  style={{ width: "100%", height: "100%" }}
                  config={{ responsive: true, displayModeBar: false }}
                />
              </div>
            </div>

            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: conv.map((c) => c.h), y: conv.map((c) => c.err), type: "scatter", mode: "lines+markers", name: "RK4 err vs h", line: { color: theme.accent2 } },
                  ]}
                  layout={{
                    title: "Convergence (RK4) at x=end (log-log)",
                    xaxis: { title: "h", type: "log" },
                    yaxis: { title: "abs error", type: "log" },
                    autosize: true,
                    responsive: true,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    margin: { t: 36, l: 48, r: 20, b: 36 },
                    font: { color: "#e5e7eb" },
                  }}
                  useResizeHandler
                  style={{ width: "100%", height: "100%" }}
                  config={{ responsive: true, displayModeBar: false }}
                />
              </div>
            </div>
          </div>

          <CollapsibleExercise
            title="Exercise 25.3.1 — RK4 order verification"
            body={`Use the error vs h data to estimate the observed order of convergence of RK4 by fitting a line on the log-log plot (slope ~ 4).`}
            solution={`On log-log plot slope ~ 4 indicates global O(h^4) behavior for RK4 on smooth problems.`}
          />
        </CardContent>
      </Card>
    </motion.div>
  );
}

/* ===========================
   Section 25.4 — Systems of Equations
   Example: Harmonic oscillator (x'' + x = 0) rewritten as system
   dx/dt = v, dv/dt = -x  with x(0)=1, v(0)=0
   =========================== */

function Section254() {
  function f(t, y) {
    // y = [x, v]
    const x = y[0], v = y[1];
    return [v, -x];
  }

  const [h, setH] = useState(0.05);
  const t0 = 0;
  const tEnd = 20;
  const n = Math.max(2, Math.floor((tEnd - t0) / h));
  const solEuler = useMemo(() => eulerSolver(f, t0, [1, 0], h, n).ys, [h, n]);
  const solRK4 = useMemo(() => rk4Solver(f, t0, [1, 0], h, n).ys, [h, n]);
  const ts = useMemo(() => linspace(t0, tEnd, solRK4.length), [solRK4.length]);

  // energy = 0.5*(v^2 + x^2) should be constant for the ideal oscillator
  const energyEuler = solEuler.map((s) => 0.5 * (s[0] * s[0] + s[1] * s[1]));
  const energyRK4 = solRK4.map((s) => 0.5 * (s[0] * s[0] + s[1] * s[1]));

  return (
    <motion.div {...fadeUp}>
      <Card className={`${theme.panel} p-4 rounded-2xl`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Layers className="w-5 h-5" /> 25.4 Systems of Equations (Harmonic Oscillator)
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid md:grid-cols-2 gap-4">
            <div>
              <Formula>{`System form:\n\n$$\\mathbf{y}' = \\mathbf{f}(t,\\mathbf{y}), \\quad \\mathbf{y}\\in\\mathbb{R}^n.$$`}</Formula>
              <Note>Runge-Kutta methods apply component-wise to vector y. For conservative systems (like harmonic oscillator), monitor conserved quantities (energy) to assess long-term behavior.</Note>
            </div>
            <div>
              <Labeled label="Step size h"><Input value={h} onChange={(e) => setH(Number(e.target.value))} className="bg-zinc-800 text-white" /></Labeled>
              <Summary items={[
                { label: "Initial state", value: "[x0=1, v0=0]" },
                { label: "Integration steps", value: n },
                { label: "Time span", value: `${t0} to ${tEnd}` },
              ]} />
            </div>
          </div>

          <div className="grid md:grid-cols-3 gap-3">
            {/* x(t), v(t) */}
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: ts, y: solRK4.map((s) => s[0]), type: "scatter", mode: "lines", name: "x (RK4)", line: { color: theme.accent2 } },
                    { x: ts, y: solRK4.map((s) => s[1]), type: "scatter", mode: "lines", name: "v (RK4)", line: { color: theme.accent } },
                    { x: ts, y: solEuler.map((s) => s[0]), type: "scatter", mode: "lines", name: "x (Euler)", line: { color: theme.warn, dash: "dash" } },
                  ]}
                  layout={{
                    title: "State vs Time (x, v)",
                    autosize: true,
                    responsive: true,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    margin: { t: 36, l: 48, r: 20, b: 36 },
                    font: { color: "#e5e7eb" },
                  }}
                  useResizeHandler
                  style={{ width: "100%", height: "100%" }}
                  config={{ responsive: true, displayModeBar: false }}
                />
              </div>
            </div>

            {/* phase plot */}
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: solRK4.map((s) => s[0]), y: solRK4.map((s) => s[1]), type: "scatter", mode: "lines", name: "phase (RK4)", line: { color: theme.accent2 } },
                    { x: solEuler.map((s) => s[0]), y: solEuler.map((s) => s[1]), type: "scatter", mode: "lines", name: "phase (Euler)", line: { color: theme.warn, dash: "dash" } },
                  ]}
                  layout={{
                    title: "Phase plot (x vs v)",
                    xaxis: { title: "x" },
                    yaxis: { title: "v" },
                    autosize: true,
                    responsive: true,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    margin: { t: 36, l: 48, r: 20, b: 36 },
                    font: { color: "#e5e7eb" },
                  }}
                  useResizeHandler
                  style={{ width: "100%", height: "100%" }}
                  config={{ responsive: true, displayModeBar: false }}
                />
              </div>
            </div>

            {/* energy conservation */}
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: ts, y: energyRK4, type: "scatter", mode: "lines", name: "Energy RK4", line: { color: theme.accent2 } },
                    { x: ts, y: energyEuler, type: "scatter", mode: "lines", name: "Energy Euler", line: { color: theme.warn } },
                  ]}
                  layout={{
                    title: "Energy vs Time (conservation check)",
                    autosize: true,
                    responsive: true,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    margin: { t: 36, l: 48, r: 20, b: 36 },
                    font: { color: "#e5e7eb" },
                  }}
                  useResizeHandler
                  style={{ width: "100%", height: "100%" }}
                  config={{ responsive: true, displayModeBar: false }}
                />
              </div>
            </div>
          </div>

          <CollapsibleExercise
            title="Exercise 25.4.1 — Symplectic integrator"
            body={`Discuss how symplectic integrators (e.g., Stormer-Verlet) preserve energy properties better than explicit RK for long-time integration of conservative systems. Implement Verlet and compare energy drift.`}
            solution={`Symplectic methods preserve a modified energy and show bounded energy error over long time, making them preferred for Hamiltonian systems like orbital mechanics.`}
          />
        </CardContent>
      </Card>
    </motion.div>
  );
}

/* ===========================
   Section 25.5 — Adaptive Runge-Kutta (RKF45-style)
   Example: logistic growth dy/dt = r*y*(1 - y/K)
   =========================== */

function Section255() {
  function f(t, y) {
    const r = 1.0;
    const K = 10.0;
    return r * y * (1 - y / K);
  }

  const [tol, setTol] = useState(1e-5);
  const [h0, setH0] = useState(0.1);
  const t0 = 0;
  const tEnd = 15;
  // fixed-step RK4 baseline:
  const Nfixed = 1000;
  const tfixed = linspace(t0, tEnd, Nfixed);
  const yfixed = useMemo(() => rk4Solver(f, t0, 1.0, (tEnd - t0) / (Nfixed - 1), Nfixed - 1).ys.map((v) => Array.isArray(v) ? v[0] : v), []);

  const adaptive = useMemo(() => rkf45Adaptive((t, y) => typeof y === "number" ? f(t, y) : f(t, y[0]), t0, 1.0, h0, tEnd, tol, 1e-8, 1.0), [tol, h0]);
  const tAdaptive = adaptive.xs;
  const yAdaptive = adaptive.ys.map((v) => Array.isArray(v) ? v[0] : v);

  // interpolation of adaptive to fixed grid for comparison (linear)
  function interpAdaptive(tq) {
    if (tAdaptive.length === 0) return 0;
    if (tq <= tAdaptive[0]) return yAdaptive[0];
    if (tq >= tAdaptive[tAdaptive.length - 1]) return yAdaptive[yAdaptive.length - 1];
    let i = 0;
    while (i < tAdaptive.length - 1 && tAdaptive[i + 1] < tq) i++;
    const tL = tAdaptive[i], tR = tAdaptive[i + 1];
    const yL = yAdaptive[i], yR = yAdaptive[i + 1];
    const u = (tq - tL) / (tR - tL);
    return yL * (1 - u) + yR * u;
  }

  const yAdaptOnFixed = tfixed.map((tt) => interpAdaptive(tt));

  // step-size visualization
  const steps = tAdaptive.map((t, i) => {
    if (i === 0) return { t: tAdaptive[i], h: tAdaptive[i] - (tAdaptive[i - 1] || 0) };
    return { t: tAdaptive[i], h: tAdaptive[i] - tAdaptive[i - 1] };
  });

  return (
    <motion.div {...fadeUp}>
      <Card className={`${theme.panel} p-4 rounded-2xl`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Compass className="w-5 h-5" /> 25.5 Adaptive Runge-Kutta Methods (RKF45)
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid md:grid-cols-2 gap-4">
            <div>
              <Formula>{`Adaptive RK (Dormand-Prince / RKF45) uses embedded methods to estimate local error and adapt step size:\n\nIf err < tol accept step and increase h; if err > tol reduce h and retry.`}</Formula>
              <Note>Adaptive methods concentrate steps where the solution changes quickly and take larger steps in smooth regions—efficient for stiff or varying problems.</Note>
            </div>
            <div>
              <div className="grid grid-cols-2 gap-2">
                <Labeled label="Tolerance tol"><Input value={tol} onChange={(e) => setTol(Number(e.target.value))} className="bg-zinc-800 text-white" /></Labeled>
                <Labeled label="Initial h0"><Input value={h0} onChange={(e) => setH0(Number(e.target.value))} className="bg-zinc-800 text-white" /></Labeled>
              </div>
              <Summary items={[
                { label: "Fixed RK4 steps", value: Nfixed },
                { label: "Adaptive steps taken", value: tAdaptive.length },
                { label: "Tol", value: tol },
              ]} />
            </div>
          </div>

          <div className="grid md:grid-cols-3 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: tfixed, y: yfixed, type: "scatter", mode: "lines", name: "Fixed RK4", line: { color: theme.accent2 } },
                    { x: tAdaptive, y: yAdaptive, type: "scatter", mode: "markers+lines", name: "Adaptive RKF45", line: { color: theme.accent } },
                  ]}
                  layout={{
                    title: "Solution: adaptive vs fixed-step",
                    autosize: true,
                    responsive: true,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    margin: { t: 36, l: 48, r: 20, b: 36 },
                    font: { color: "#e5e7eb" },
                  }}
                  useResizeHandler
                  style={{ width: "100%", height: "100%" }}
                  config={{ responsive: true, displayModeBar: false }}
                />
              </div>
            </div>

            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: tfixed, y: yAdaptOnFixed.map((y, i) => Math.abs(y - yfixed[i])), type: "scatter", mode: "lines", name: "abs difference", line: { color: theme.warn } },
                  ]}
                  layout={{
                    title: "Difference (adaptive vs fixed reference)",
                    xaxis: { title: "t" },
                    autosize: true,
                    responsive: true,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    margin: { t: 36, l: 48, r: 20, b: 36 },
                    font: { color: "#e5e7eb" },
                  }}
                  useResizeHandler
                  style={{ width: "100%", height: "100%" }}
                  config={{ responsive: true, displayModeBar: false }}
                />
              </div>
            </div>

            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: steps.map((s) => s.t), y: steps.map((s) => s.h), type: "bar", name: "step size h", marker: { color: theme.accent } },
                  ]}
                  layout={{
                    title: "Adaptive step-size vs time",
                    xaxis: { title: "t" },
                    yaxis: { title: "h (step size)" },
                    autosize: true,
                    responsive: true,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    margin: { t: 36, l: 48, r: 20, b: 36 },
                    font: { color: "#e5e7eb" },
                  }}
                  useResizeHandler
                  style={{ width: "100%", height: "100%" }}
                  config={{ responsive: true, displayModeBar: false }}
                />
              </div>
            </div>
          </div>

          <CollapsibleExercise
            title="Exercise 25.5.1 — Tolerance effect"
            body={`Experiment with tolerance tol across orders of magnitude (1e-3 .. 1e-8). How does the number of steps and error change? Plot cost vs accuracy.`}
            solution={`Lower tol -> more steps -> lower error. Adaptive steppers concentrate work in difficult regions, often giving better cost/accuracy than fixed-step RK4.`}
          />
        </CardContent>
      </Card>
    </motion.div>
  );
}

/* ===========================
   Documentation / Problems panel (end)
   =========================== */

function ChapterDocs() {
  const [open, setOpen] = useState("25.1 Euler’s Method");
  const docs = {
  "25.1 Euler’s Method": [
    `Euler's method (explicit, first-order): $$y_{n+1} = y_n + h f(x_n, y_n).$$`,
    `Local truncation error: $$\\tau_n \\sim \\mathcal{O}(h^2),$$ global error: $$E_n \\sim \\mathcal{O}(h).$$`,
    `Simple and easy to implement, but low accuracy; step size must be small for stiff problems.`,
    `Useful for initial exploration of solutions or when high precision is not required.`,
    `Euler can be unstable for stiff equations; implicit versions (backward Euler) improve stability.`
  ],
  "25.2 Improvements of Euler’s Method": [
    `Heun's method (predictor-corrector) and the Midpoint method are explicit second-order methods:`,
    `Heun: $$y_{n+1} = y_n + \\frac{h}{2} [f(x_n, y_n) + f(x_n+h, y_n+h f(x_n, y_n))].$$`,
    `Midpoint: $$y_{n+1} = y_n + h f\\left(x_n + \\frac{h}{2}, y_n + \\frac{h}{2} f(x_n, y_n)\\right).$$`,
    `They improve accuracy over Euler's method with minimal additional cost.`,
    `Good choice when RK4 is overkill but higher accuracy than Euler is desired.`
  ],
"25.3 Runge-Kutta Methods": [
  `RK4 is the classical fourth-order Runge-Kutta method with excellent balance of accuracy and computational cost:`,
  `$$
  \\begin{aligned}
  k_1 &= f(x_n, y_n), \\\\
  k_2 &= f\\left(x_n + \\frac{h}{2}, y_n + \\frac{h}{2} k_1\\right), \\\\
  k_3 &= f\\left(x_n + \\frac{h}{2}, y_n + \\frac{h}{2} k_2\\right), \\\\
  k_4 &= f(x_n + h, y_n + h k_3), \\\\
  y_{n+1} &= y_n + \\frac{h}{6} (k_1 + 2 k_2 + 2 k_3 + k_4).
  \\end{aligned}
  $$`,
  `Local truncation error: $\\mathcal{O}(h^5)$, global error: $\\mathcal{O}(h^4)$.`,
  `Widely used for smooth problems; robust and stable for moderate step sizes.`
],

  "25.4 Systems of Equations": [
    `For systems of ODEs, apply numerical methods component-wise:`,
    `$$\\dot{x} = f(t, x, y), \\quad \\dot{y} = g(t, x, y).$$`,
    `Example RK4 for a system: apply $k_1, k_2, k_3, k_4$ for each variable simultaneously.`,
    `Monitor invariants (energy, momentum, total mass) for long-time integration to check stability.`,
    `Vectorization in software (MATLAB, NumPy) allows efficient computation of multi-dimensional systems.`
  ],
  "25.5 Adaptive Runge-Kutta Methods": [
    `Adaptive methods adjust step size to meet a specified error tolerance, improving efficiency and stability.`,
    `Example: Dormand-Prince (RKF45) uses an embedded fourth- and fifth-order pair to estimate error:`,
    `$$E_n = |y_{n+1}^{(5)} - y_{n+1}^{(4)}|.$$`,
    `Step size $h$ is updated as: $$h_{new} = h \\left( \\frac{\\text{tolerance}}{E_n} \\right)^{1/5}.$$`,
    `Useful for stiff or highly varying solutions; avoids unnecessary computations on smooth regions.`,
    `Widely implemented in numerical libraries: MATLAB \`ode45\`, Python \`scipy.integrate.solve_ivp(method='RK45')\`.`
  ],
};

  return (
    <div className="w-full mt-6 bg-zinc-900 border border-zinc-700 rounded-xl shadow-lg overflow-hidden">
      <div className="p-3 border-b border-zinc-700 flex items-center justify-between">
        <div className="flex items-center gap-3">
          <BookOpen className="w-5 h-5 text-cyan-400" />
          <div>
            <div className="text-zinc-100 font-semibold">Chapter 25 — Notes & Derivations</div>
            <div className="text-zinc-400 text-xs">KaTeX-rendered formulas and short references</div>
          </div>
        </div>
        <div className="text-zinc-400 text-xs">Interactive demos above</div>
      </div>

      <div className="grid grid-cols-1 md:grid-cols-4">
        <div className="col-span-1 border-r border-zinc-800">
          {Object.keys(docs).map((k) => (
            <button key={k} onClick={() => setOpen(k)} className={`w-full p-3 text-left hover:bg-zinc-800/40 flex items-center justify-between ${open === k ? "bg-zinc-800/20" : ""}`}>
              <div className="flex items-center gap-2">
                <List className="w-4 h-4 text-zinc-300" />
                <div className="text-zinc-100 text-sm">{k}</div>
              </div>
              {open === k ? <ChevronDown className="w-4 h-4 text-zinc-400" /> : <ChevronRight className="w-4 h-4 text-zinc-400" />}
            </button>
          ))}
        </div>

        <div className="col-span-3 p-4 overflow-auto" style={{ maxHeight: 460 }}>
          <motion.h3 initial={{ opacity: 0, y: -6 }} animate={{ opacity: 1, y: 0 }} className="text-xl font-semibold text-zinc-100 mb-3">{open}</motion.h3>
          <div className="text-zinc-300 space-y-3 text-sm leading-relaxed">
            {docs[open].map((p, i) => (
              <div key={i}>
                <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>
                  {p}
                </ReactMarkdown>
              </div>
            ))}
            <div className="mt-3 text-zinc-400 text-xs">References: Burden & Faires; Shampine & Gordon; Hairer, Norsett & Wanner.</div>
          </div>
        </div>
      </div>
    </div>
  );
}

/* ===========================
   Problems panel
   =========================== */

function ProblemsPanel() {
  return (
    <div className={`${theme.panel} p-4 rounded-2xl`}>
      <div className="flex items-center justify-between mb-3">
        <div className="text-xl text-cyan-300 font-semibold flex items-center gap-2"><FileText className="w-5 h-5" /> Problems</div>
        <Button variant="ghost" size="sm" onClick={() => window.scrollTo({ top: 0, behavior: "smooth" })}><RefreshCcw className="w-4 h-4" /> Top</Button>
      </div>

      <div className="grid md:grid-cols-2 gap-3">
        <CollapsibleExercise
          title="Problem 25.1 — Euler and RK4 comparison"
          body={`For dy/dx = -2y + x, y(0)=1 compute y(2) using Euler with h=0.2 and RK4 with h=0.2. Estimate absolute errors using a dense RK4 reference. Which method is more accurate per step?`}
          solution={`RK4 is significantly more accurate per step; error scaling demonstrates Euler O(h) vs RK4 O(h^4).`}
        />
        <CollapsibleExercise
          title="Problem 25.2 — Midpoint vs Heun"
          body={`Implement both midpoint and Heun on the sample ODE and measure global error vs h. Compare computational cost (function evaluations) to reach a target error of 1e-4.`}
          solution={`Midpoint and Heun are both second-order; cost tradeoffs depend on number of function evaluations per step (Midpoint uses 2 evals, Heun uses 2 evals as well).`}
        />
        <CollapsibleExercise
          title="Problem 25.3 — Long-term stability"
          body={`For the harmonic oscillator, compare Euler, RK4, and symplectic Verlet over long time T=1000. Plot energy drift. Which method is best for conserving energy?`}
          solution={`Symplectic Verlet shows bounded energy error; explicit RK4 shows slow drift while Euler quickly diverges in energy.`}
        />
        <CollapsibleExercise
          title="Problem 25.4 — Adaptive RK usage"
          body={`Use an adaptive RK (RKF45) to integrate a stiff-ish test equation and compare steps taken vs fixed-step RK4 achieving similar error. Which is more efficient?`}
          solution={`Adaptive RK focuses effort where needed, often outperforming fixed-step when solution has localized rapid changes.`}
        />
      </div>
    </div>
  );
}

/* ===========================
   Page assembly (export default)
   =========================== */

export default function Chapter25() {
  return (
    <div className={`p-6 space-y-6 ${theme.bg}`}>
      <motion.div initial={{ opacity: 0, y: -10 }} animate={{ opacity: 1, y: 0 }} className="space-y-1">
        <h1 className="text-3xl font-bold" style={{ color: theme.accent2 }}>Runge-Kutta Methods</h1>
        <div className="text-zinc-400">Euler, improved Euler, RK4, systems, and adaptive Runge-Kutta (interactive demos + derivations).</div>
      </motion.div>

      {/* Sections */}
      <Section251 />
      <Section252 />
      <Section253 />
      <Section254 />
      <Section255 />

      {/* Docs and Problems */}
      <ChapterDocs />
      <ProblemsPanel />
      <BottomBar/>
    </div>
  );
}
