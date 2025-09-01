// src/pages/Chapter23.jsx
// Chapter 23 — Numerical Differentiation
// Mirrors Chapter7.jsx design patterns
// All Plotly charts constrained in responsive containers to prevent overflow

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
  Ruler,
  Activity,
  Layers,
  Zap,
  BookOpen,
  List,
  ChevronRight,
  ChevronDown,
  RefreshCcw,
  FileText,
} from "lucide-react";
import BottomBar from "../components/BottomBar";

/* ============================
   Theme & small helpers
   ============================ */

const theme = {
  bg: "bg-zinc-950",
  panel: "bg-zinc-900/60 border border-zinc-700 overflow-hidden",
  text: "text-zinc-200",
  subtext: "text-zinc-400",
  accent: "#22d3ee",
  accent2: "#34d399",
  warn: "#f59e0b",
  danger: "#ef4444",
};

const fadeUp = {
  initial: { opacity: 0, y: 12 },
  animate: { opacity: 1, y: 0 },
  transition: { duration: 0.28 },
};

/* ============================
   Numeric utilities
   ============================ */

function linspace(a, b, n) {
  if (n <= 1) return [a];
  const dx = (b - a) / (n - 1);
  const arr = new Array(n);
  for (let i = 0; i < n; i++) arr[i] = a + i * dx;
  return arr;
}

function gaussianNoise(std = 1) {
  let u = 0,
    v = 0;
  while (u === 0) u = Math.random();
  while (v === 0) v = Math.random();
  return Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v) * std;
}

function solveLinearSystem(A, b) {
  const n = A.length;
  const M = new Array(n);
  for (let i = 0; i < n; i++) M[i] = A[i].slice();
  const rhs = b.slice();
  for (let k = 0; k < n; k++) {
    let piv = k;
    for (let i = k + 1; i < n; i++) if (Math.abs(M[i][k]) > Math.abs(M[piv][k])) piv = i;
    if (Math.abs(M[piv][k]) < 1e-14) return null;
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

/* Finite differences (several schemes) */

function forwardDiff(y, h) {
  const n = y.length;
  const out = new Array(n).fill(NaN);
  for (let i = 0; i < n - 1; i++) out[i] = (y[i + 1] - y[i]) / h;
  out[n - 1] = (y[n - 1] - y[n - 2]) / h;
  return out;
}

function backwardDiff(y, h) {
  const n = y.length;
  const out = new Array(n).fill(NaN);
  out[0] = (y[1] - y[0]) / h;
  for (let i = 1; i < n; i++) out[i] = (y[i] - y[i - 1]) / h;
  return out;
}

function centralDiff(y, h) {
  const n = y.length;
  const out = new Array(n).fill(NaN);
  out[0] = (y[1] - y[0]) / h;
  for (let i = 1; i < n - 1; i++) out[i] = (y[i + 1] - y[i - 1]) / (2 * h);
  out[n - 1] = (y[n - 1] - y[n - 2]) / h;
  return out;
}

function fivePointCentral(y, h) {
  const n = y.length;
  const out = new Array(n).fill(NaN);
  if (n < 5) return centralDiff(y, h);
  out[0] = (y[1] - y[0]) / h;
  out[1] = (-y[3] + 4 * y[2] - 3 * y[1]) / (2 * h);
  for (let i = 2; i < n - 2; i++) {
    out[i] = (y[i - 2] - 8 * y[i - 1] + 8 * y[i + 1] - y[i + 2]) / (12 * h);
  }
  out[n - 2] = (3 * y[n - 2] - 4 * y[n - 3] + y[n - 4]) / (2 * h);
  out[n - 1] = (y[n - 1] - y[n - 2]) / h;
  return out;
}

/* Savitzky-Golay derivative (local poly fit) */

function savitzkyGolayDerivative(t, y, window = 7, degree = 3) {
  const n = y.length;
  const half = Math.floor(window / 2);
  const deriv = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    const lo = Math.max(0, i - half);
    const hi = Math.min(n - 1, i + half);
    const m = hi - lo + 1;
    const A = new Array(m);
    const b = new Array(m);
    for (let k = 0; k < m; k++) {
      const idx = lo + k;
      const dt = t[idx] - t[i];
      A[k] = new Array(degree + 1);
      for (let p = 0; p <= degree; p++) A[k][p] = Math.pow(dt, p);
      b[k] = y[idx];
    }
    const ATA = new Array(degree + 1).fill(0).map(() => new Array(degree + 1).fill(0));
    const ATb = new Array(degree + 1).fill(0);
    for (let r = 0; r < degree + 1; r++) {
      for (let c = 0; c < degree + 1; c++) {
        let s = 0;
        for (let k = 0; k < m; k++) s += A[k][r] * A[k][c];
        ATA[r][c] = s;
      }
      let sb = 0;
      for (let k = 0; k < m; k++) sb += A[k][r] * b[k];
      ATb[r] = sb;
    }
    const coeffs = solveLinearSystem(ATA, ATb);
    deriv[i] = coeffs ? coeffs[1] ?? 0 : 0;
  }
  return deriv;
}

/* Unequal spacing derivative via local quadratic fit */

function derivativeUnequal(xarr, yarr) {
  const n = xarr.length;
  const dy = new Array(n).fill(NaN);
  for (let i = 0; i < n; i++) {
    // select up to 3 points around i
    let idxs = [i];
    let l = i - 1,
      r = i + 1;
    while (idxs.length < Math.min(3, n)) {
      if (l >= 0) {
        idxs.push(l);
        l--;
      }
      if (idxs.length >= Math.min(3, n)) break;
      if (r < n) {
        idxs.push(r);
        r++;
      }
    }
    idxs.sort((a, b) => a - b);
    const m = idxs.length;
    const A = new Array(m);
    const bb = new Array(m);
    for (let k = 0; k < m; k++) {
      const xi = xarr[idxs[k]];
      A[k] = [xi * xi, xi, 1];
      bb[k] = yarr[idxs[k]];
    }
    const ATA = new Array(3).fill(0).map(() => new Array(3).fill(0));
    const ATb = new Array(3).fill(0);
    for (let r2 = 0; r2 < 3; r2++) {
      for (let c = 0; c < 3; c++) {
        let s = 0;
        for (let k = 0; k < m; k++) s += (A[k][r2] ?? 0) * (A[k][c] ?? 0);
        ATA[r2][c] = s;
      }
      let sb = 0;
      for (let k = 0; k < m; k++) sb += (A[k][r2] ?? 0) * bb[k];
      ATb[r2] = sb;
    }
    const coeffs = solveLinearSystem(ATA, ATb);
    if (coeffs) {
      const [a, bcoef] = coeffs;
      dy[i] = 2 * a * xarr[i] + bcoef;
    } else dy[i] = NaN;
  }
  return dy;
}

/* Partial derivatives for 2D grid */

function partialDerivatives2D(Z, dx, dy) {
  const ny = Z.length;
  const nx = Z[0].length;
  const zx = new Array(ny).fill(0).map(() => new Array(nx).fill(0));
  const zy = new Array(ny).fill(0).map(() => new Array(nx).fill(0));
  for (let j = 0; j < ny; j++) {
    for (let i = 0; i < nx; i++) {
      if (i === 0) zx[j][i] = (Z[j][i + 1] - Z[j][i]) / dx;
      else if (i === nx - 1) zx[j][i] = (Z[j][i] - Z[j][i - 1]) / dx;
      else zx[j][i] = (Z[j][i + 1] - Z[j][i - 1]) / (2 * dx);
      if (j === 0) zy[j][i] = (Z[j + 1][i] - Z[j][i]) / dy;
      else if (j === ny - 1) zy[j][i] = (Z[j][i] - Z[j - 1][i]) / dy;
      else zy[j][i] = (Z[j + 1][i] - Z[j - 1][i]) / (2 * dy);
    }
  }
  return { zx, zy };
}

/* Richardson extrapolation for numerical derivative estimate function D(h) */

function richardsonExtrapolate(Dh_func, h) {
  const D_h = Dh_func(h);
  const D_h2 = Dh_func(h / 2);
  if (Array.isArray(D_h) && Array.isArray(D_h2)) {
    return D_h.map((_, i) => (4 * D_h2[i] - D_h[i]) / 3);
  } else if (!Array.isArray(D_h) && !Array.isArray(D_h2)) {
    return (4 * D_h2 - D_h) / 3;
  } else return null;
}

/* ============================
   Reusable UI components
   (match Chapter7 patterns)
   ============================ */

function Formula({ children, color = theme.accent }) {
  return (
    <div className="my-2">
      <div className="px-3 py-2 rounded-lg border border-zinc-700 bg-zinc-900/80 overflow-x-auto">
        <div className="text-xs uppercase tracking-widest mb-1" style={{ color }}>
          formula
        </div>
        <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>
          {children}
        </ReactMarkdown>
      </div>
    </div>
  );
}

function NoteBox({ children }) {
  return <div className="text-sm text-zinc-300">{children}</div>;
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
        <div key={i} className="rounded-lg bg-zinc-900/70 border border-zinc-700 px-3 py-2 flex flex-col">
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
      <div className="flex items-center justify-between gap-2">
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
        <div className="mt-3 border-t border-zinc-800 pt-3">
          <div className="text-xs text-zinc-400 mb-1">Solution</div>
          <div className="text-zinc-200">
            <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>
              {solution}
            </ReactMarkdown>
          </div>
        </div>
      )}
    </div>
  );
}

/* ============================
   Docs text (KaTeX-ready)
   ============================ */

const docsText = {
  "23.1 High-Accuracy Differentiation Formulas": [
    `Central difference (second-order):\n\n$$f'(x) \\approx \\frac{f(x+h)-f(x-h)}{2h} + \\mathcal{O}(h^2).$$`,
    `Five-point (fourth-order) formula:\n\n$$f'(x) \\approx \\frac{1}{12h} \\big(f(x-2h) - 8f(x-h) + 8f(x+h) - f(x+2h)\\big) + \\mathcal{O}(h^4).$$`,
    "Higher-order formulas reduce truncation error but may amplify round-off for very small $h$.",
    "Choose $h$ to balance truncation error ($\\sim h^2$ or $h^4$) and round-off ($\\sim 1/h$).",
    "Use these formulas on smooth functions; for noisy data, smoothing is recommended first."
  ],
  "23.2 Richardson Extrapolation": [
    `If $D(h) = f'(x) + C h^2 + \\mathcal{O}(h^4)$, Richardson extrapolation improves accuracy:\n\n$$D_{\\text{rich}} = \\frac{4 D(h/2) - D(h)}{3} = f'(x) + \\mathcal{O}(h^4).$$`,
    "Extrapolation can be applied iteratively for even higher-order accuracy.",
    "Use with central difference or other finite-difference estimates to reduce error without decreasing $h$ too much.",
    "Practical tip: Richardson reduces truncation error but still sensitive to noise in $f(x)$."
  ],
  "23.3 Derivatives of Unequally Spaced Data": [
    "When nodes are irregular, finite-difference coefficients can be derived from polynomial interpolation:",
    "Fit a polynomial of degree $n$ through $n+1$ nearby points and differentiate analytically at the point of interest.",
    "Local polynomial fitting (quadratic, cubic) improves robustness against uneven spacing.",
    "Weighted least squares can be used for noisy or scattered data.",
    "Ensure that the interpolation stencil includes points sufficiently close to avoid ill-conditioning."
  ],
  "23.4 Derivatives and Integrals for Data with Errors": [
    "Numerical differentiation amplifies high-frequency noise; direct finite differences may be unstable.",
    "Smoothing options include moving averages, Gaussian filters, or Savitzky–Golay filters before differentiation.",
    "Alternatively, fit a model (polynomial, spline, or regression) to the data and differentiate the fitted function.",
    "For integration, cumulative sums are more stable, and smoothing can reduce error propagation.",
    "Always visualize residuals and derivative estimates to detect artifacts caused by noise."
  ],
  "23.5 Partial Derivatives": [
    "For a 2D scalar field $f(x,y)$ on a grid:",
    "$$\\frac{\\partial f}{\\partial x}\\bigg|_{i,j} \\approx \\frac{f_{i+1,j} - f_{i-1,j}}{2\\Delta x}, \\quad \\frac{\\partial f}{\\partial y}\\bigg|_{i,j} \\approx \\frac{f_{i,j+1} - f_{i,j-1}}{2\\Delta y}.$$",
    "Higher-order approximations (fourth-order) can also be used for smoother results:",
    "$$\\frac{\\partial f}{\\partial x}\\bigg|_{i,j} \\approx \\frac{-f_{i+2,j} + 8 f_{i+1,j} - 8 f_{i-1,j} + f_{i-2,j}}{12 \\Delta x}.$$",
    "For mixed derivatives, use combinations like:",
    "$$\\frac{\\partial^2 f}{\\partial x \\partial y} \\approx \\frac{f_{i+1,j+1} - f_{i+1,j-1} - f_{i-1,j+1} + f_{i-1,j-1}}{4 \\Delta x \\Delta y}.$$",
    "Boundary points may require one-sided differences or extrapolation."
  ],
  "23.6 Numerical Integration/Differentiation with Software Packages": [
    "Python: `numpy.gradient`, `scipy.misc.derivative`, `scipy.integrate.simps`, `scipy.integrate.trapz`.",
    "MATLAB: `gradient`, `diff`, `integral`, `trapz`, `cumtrapz`.",
    "For noisy data, use smoothing functions: `scipy.signal.savgol_filter` or MATLAB's `sgolayfilt`.",
    "High-dimensional differentiation: `numpy.gradient` can handle n-dimensional arrays directly.",
    "Always check step sizes and spacing; functions assume uniform grids unless specified."
  ],
};


/* ============================
   Sections (23.1 - 23.6)
   Each section follows pattern: inputs, summary, responsive plots in wrapper divs
   ============================ */

/* -------- 23.1 High-Accuracy -------- */
function Section231() {
  const [funcExpr, setFuncExpr] = useState("Math.sin(x)");
  const [x0, setX0] = useState(1.0);
  const [h, setH] = useState(0.2);
  const [N, setN] = useState(201);

  const fn = useMemo(() => {
    try {
      // eslint-disable-next-line no-new-func
      return new Function("x", `return ${funcExpr};`);
    } catch {
      return (x) => 0;
    }
  }, [funcExpr]);

  const X = useMemo(() => linspace(x0 - 2, x0 + 2, Number(N)), [x0, N]);
  const Y = useMemo(() => X.map((xx) => {
    try { return fn(xx); } catch { return NaN; }
  }), [X, fn]);

  const central = useMemo(() => (fn(x0 + h) - fn(x0 - h)) / (2 * h), [fn, x0, h]);
  const fivePoint = useMemo(() => (fn(x0 - 2 * h) - 8 * fn(x0 - h) + 8 * fn(x0 + h) - fn(x0 + 2 * h)) / (12 * h), [fn, x0, h]);

  const tinyRef = useMemo(() => {
    const hh = 1e-8;
    return (fn(x0 + hh) - fn(x0 - hh)) / (2 * hh);
  }, [fn, x0]);

  const summaryItems = [
    { label: "Central (2nd order)", value: central.toPrecision(6) },
    { label: "Five-point (4th order)", value: fivePoint.toPrecision(6) },
    { label: "Tiny-h ref", value: tinyRef.toPrecision(8) },
  ];

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Ruler className="w-5 h-5" />
            23.1 High-Accuracy Differentiation Formulas
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid md:grid-cols-3 gap-3">
            <Labeled label="Function (JS expression)">
              <Input value={funcExpr} onChange={(e) => setFuncExpr(e.target.value)} className="bg-zinc-800 text-white" />
            </Labeled>
            <Labeled label="x₀">
              <Input value={x0} onChange={(e) => setX0(Number(e.target.value))} className="bg-zinc-800 text-white" />
            </Labeled>
            <Labeled label="Step h">
              <Input value={h} onChange={(e) => setH(Number(e.target.value))} className="bg-zinc-800 text-white" />
            </Labeled>
          </div>

          <div className="grid md:grid-cols-2 gap-3">
            {/* Left: function plot */}
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: X, y: Y, type: "scatter", mode: "lines", name: "f(x)", line: { color: theme.accent } },
                    { x: [x0], y: [fn(x0)], type: "scatter", mode: "markers", marker: { color: theme.accent2, size: 8 }, name: "x₀" },
                    { x: [x0 - h], y: [fn(x0 - h)], type: "scatter", mode: "markers", marker: { color: "#f59e0b", size: 6 }, name: "x₀-h" },
                    { x: [x0 + h], y: [fn(x0 + h)], type: "scatter", mode: "markers", marker: { color: "#f59e0b", size: 6 }, name: "x₀+h" },
                  ]}
                  layout={{
                    title: "Function around x₀",
                    autosize: true,
                    responsive: true,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    font: { color: "#e5e7eb" },
                    margin: { t: 36, l: 48, r: 20, b: 36 },
                  }}
                  useResizeHandler
                  style={{ width: "100%", height: "100%" }}
                  config={{ responsive: true, displayModeBar: false }}
                />
              </div>
            </div>

            {/* Right: summary + formula */}
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-3 overflow-hidden">
              <Summary items={summaryItems} />
              <Formula>{`**Five-point (4th order)**\n\n$$f'(x)\\approx\\frac{f(x-2h)-8f(x-h)+8f(x+h)-f(x+2h)}{12h}.$$`}</Formula>
              <NoteBox>For very small h, rounding error increases — pick h by balancing truncation & round-off (e.g., sqrt(epsilon) scale).</NoteBox>
            </div>
          </div>

          <CollapsibleExercise
        title="Exercise 23.1.1 — Convergence study"
        body={`For $f(x)=\\sin x$ at $x=1$, compute derivatives using central and five-point formulas for $h=10^{-k}$, $k=1\\ldots6$. Plot the absolute error vs $h$ (log-log).`}
        solution={`Central difference should display slope ≈ -2, five-point formula slope ≈ -4 (until round-off errors dominate).`}
        />

        </CardContent>
      </Card>
    </motion.div>
  );
}

/* -------- 23.2 Richardson -------- */

function Section232() {
  const [funcExpr, setFuncExpr] = useState("Math.exp(x)");
  const [x0, setX0] = useState(0.5);
  const [h, setH] = useState(0.1);

  const fn = useMemo(() => {
    try { return new Function("x", `return ${funcExpr};`); } catch { return (x) => 0; }
  }, [funcExpr]);

  const D = (hh) => (fn(x0 + hh) - fn(x0 - hh)) / (2 * hh);
  const D_rich = useMemo(() => richardsonExtrapolate(D, h), [D, h]);
  const ref = useMemo(() => {
    const eps = 1e-9;
    return (fn(x0 + eps) - fn(x0 - eps)) / (2 * eps);
  }, [fn, x0]);

  const hvals = useMemo(() => linspace(h / 16, h, 8).reverse(), [h]);

  const errs = useMemo(() => hvals.map((hv) => Math.abs(D(hv) - ref)), [hvals, D, ref]);
  const errsRich = useMemo(() => hvals.map((hv) => Math.abs(richardsonExtrapolate(D, hv) - ref)), [hvals, D, ref]);

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Activity className="w-5 h-5" />
            23.2 Richardson Extrapolation
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid md:grid-cols-3 gap-3">
            <Labeled label="Function (JS expr)">
              <Input value={funcExpr} onChange={(e) => setFuncExpr(e.target.value)} className="bg-zinc-800 text-white" />
            </Labeled>
            <Labeled label="x₀">
              <Input value={x0} onChange={(e) => setX0(Number(e.target.value))} className="bg-zinc-800 text-white" />
            </Labeled>
            <Labeled label="h">
              <Input value={h} onChange={(e) => setH(Number(e.target.value))} className="bg-zinc-800 text-white" />
            </Labeled>
          </div>

          <div className="grid md:grid-cols-2 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: hvals, y: errs, type: "scatter", mode: "lines+markers", name: "Central error", line: { color: theme.accent } },
                    { x: hvals, y: errsRich, type: "scatter", mode: "lines+markers", name: "Richardson error", line: { color: theme.accent2 } },
                  ]}
                  layout={{
                    title: "Error vs h (log-log)",
                    xaxis: { type: "log", title: "h" },
                    yaxis: { type: "log", title: "abs(error)" },
                    autosize: true,
                    responsive: true,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    font: { color: "#e5e7eb" },
                    margin: { t: 36, l: 48, r: 20, b: 36 },
                  }}
                  useResizeHandler
                  style={{ width: "100%", height: "100%" }}
                  config={{ responsive: true, displayModeBar: false }}
                />
              </div>
            </div>

            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-3 overflow-hidden">
              <Summary items={[
                { label: "D(h)", value: D(h).toPrecision(8) },
                { label: "D(h/2)", value: D(h / 2).toPrecision(8) },
                { label: "Richardson", value: (Array.isArray(D_rich) ? D_rich[0] : D_rich).toPrecision ? (Array.isArray(D_rich) ? D_rich[0].toPrecision(8) : D_rich.toPrecision(8)) : `${D_rich}` },
                { label: "Reference", value: ref.toPrecision(10) },
              ]} />
              <Formula>{`If $D(h)=f'(x)+C h^2+\\mathcal{O}(h^4)$ then:\n\n$$D_{rich} = \\frac{4 D(h/2) - D(h)}{3}$$`}</Formula>
              <NoteBox>Richardson is highly effective when the error expansion is dominated by a known power of h (e.g., h²). Beware of round-off when h is tiny.</NoteBox>
            </div>
          </div>

<CollapsibleExercise
  title="Exercise 23.2.1 — Richardson in practice"
  body={`Apply Richardson extrapolation for $f(x)=\\ln x$ at $x=2$ with $h=0.1$ and compare with the analytic derivative $1/2$.`}
  solution={`Compute $D(h)$, $D(h/2)$, and $D_{\\text{rich}}$. Expect improved accuracy approaching $0.5$.`}
/>

        </CardContent>
      </Card>
    </motion.div>
  );
}

/* -------- 23.3 Unequally Spaced Data -------- */

function Section233() {
  const [n, setN] = useState(60);
  const [noiseStd, setNoiseStd] = useState(0.0);
  const [func, setFunc] = useState("Math.sin(x)");

  const xarr = useMemo(() => {
    const base = linspace(0, Math.PI * 2, n);
    // small deterministic pseudo-random perturbation for reproducibility
    return base.map((v, i) => v + (Math.sin(i * 17.3) * 0.02));
  }, [n]);

  const fn = useMemo(() => {
    try { return new Function("x", `return ${func};`); } catch { return (x) => 0; }
  }, [func]);

  const yarr = useMemo(() => xarr.map((xi) => fn(xi) + gaussianNoise(noiseStd)), [xarr, fn, noiseStd]);

  const derivUneq = useMemo(() => derivativeUnequal(xarr, yarr), [xarr, yarr]);

  const uniformX = useMemo(() => linspace(xarr[0], xarr[xarr.length - 1], n), [xarr, n]);
  const uniformY = useMemo(() => uniformX.map((xu) => {
    // linear interpolation
    let i = 0;
    while (i < xarr.length - 1 && xarr[i + 1] < xu) i++;
    const x0 = xarr[i], x1 = xarr[Math.min(i + 1, xarr.length - 1)];
    const y0 = yarr[i], y1 = yarr[Math.min(i + 1, xarr.length - 1)];
    if (x1 === x0) return y0;
    return y0 + ((xu - x0) / (x1 - x0)) * (y1 - y0);
  }), [uniformX, xarr, yarr]);

  const centralUniform = useMemo(() => centralDiff(uniformY, uniformX[1] - uniformX[0]), [uniformY, uniformX]);

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Activity className="w-5 h-5" />
            23.3 Derivatives of Unequally Spaced Data
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid md:grid-cols-4 gap-3">
            <Labeled label="Points n">
              <Input value={n} onChange={(e) => setN(Number(e.target.value))} className="bg-zinc-800 text-white"/>
            </Labeled>
            <Labeled label="Noise std">
              <Input value={noiseStd} onChange={(e) => setNoiseStd(Number(e.target.value))} className="bg-zinc-800 text-white"/>
            </Labeled>
            <Labeled label="Function (expr)">
              <Input value={func} onChange={(e) => setFunc(e.target.value)} className="bg-zinc-800 text-white"/>
            </Labeled>
            <div className="flex items-end">
              <Button  
className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black"
onClick={() => window.location.reload()}>Regenerate</Button>
            </div>
          </div>

          <div className="grid md:grid-cols-2 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: xarr, y: yarr, type: "scatter", mode: "markers", name: "Unequal", marker: { color: theme.accent } },
                    { x: uniformX, y: uniformY, type: "scatter", mode: "lines", name: "Uniform interp", line: { color: theme.accent2 } },
                  ]}
                  layout={{
                    title: "Unequal samples vs interpolated uniform",
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
                    { x: xarr, y: derivUneq, type: "scatter", mode: "lines+markers", name: "Local poly derivative", line: { color: theme.accent } },
                    { x: uniformX, y: centralUniform, type: "scatter", mode: "lines", name: "Uniform central", line: { color: theme.accent2 } },
                  ]}
                  layout={{
                    title: "Derivative estimates",
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
  title="Exercise 23.3.1 — Local vs global"
  body={`Compare local quadratic derivative estimates with global cubic spline derivative (if available) on noisy data. Which is more robust?`}
  solution={`Local quadratic derivatives are less sensitive to outliers and local noise; global cubic splines are smooth but may introduce bias in the presence of outliers.`}
/>

        </CardContent>
      </Card>
    </motion.div>
  );
}

/* -------- 23.4 Derivatives & Integrals for Data with Errors -------- */

function Section234() {
  const [N, setN] = useState(220);
  const [noiseStd, setNoiseStd] = useState(0.2);
  const [method, setMethod] = useState("none");

  const t = useMemo(() => linspace(0, 4 * Math.PI, N), [N]);
  const clean = useMemo(() => t.map((tt) => Math.sin(2 * tt) + 0.5 * Math.cos(0.5 * tt)), [t]);
  const noisy = useMemo(() => clean.map((v) => v + gaussianNoise(noiseStd)), [clean, noiseStd]);

  const noisyMov = useMemo(() => {
    if (method === "movavg") {
      const w = 7;
      return movingAverage(noisy, w);
    }
    return noisy;
  }, [noisy, method]);

  function movingAverage(signal, window = 5) {
    const n = signal.length; const half = Math.floor(window / 2); const out = new Array(n);
    for (let i = 0; i < n; i++) {
      let s = 0; let cnt = 0;
      for (let j = i - half; j <= i + half; j++) {
        if (j >= 0 && j < n) { s += signal[j]; cnt++; }
      }
      out[i] = s / cnt;
    }
    return out;
  }

  const derivRaw = useMemo(() => centralDiff(noisy, t[1] - t[0]), [noisy, t]);
  const derivMov = useMemo(() => centralDiff(noisyMov, t[1] - t[0]), [noisyMov, t]);
  const derivSG = useMemo(() => savitzkyGolayDerivative(t, noisy, 11, 3), [t, noisy]);

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Zap className="w-5 h-5" />
            23.4 Derivatives & Integrals for Data with Errors
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid md:grid-cols-4 gap-3">
            <Labeled label="Samples (N)">
              <Input value={N} onChange={(e) => setN(Number(e.target.value))} className="bg-zinc-800 text-white" />
            </Labeled>
            <Labeled label="Noise std">
              <Input value={noiseStd} onChange={(e) => setNoiseStd(Number(e.target.value))} className="bg-zinc-800 text-white" />
            </Labeled>
            <Labeled label="Smoothing">
              <select value={method} onChange={(e) => setMethod(e.target.value)} className="bg-zinc-800 text-white rounded px-2 py-1">
                <option value="none">None</option>
                <option value="movavg">Moving average</option>
                <option value="savgol">Savitzky–Golay (derivative)</option>
              </select>
            </Labeled>
            <div className="flex items-end">
              <Button 
className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black"
 onClick={() => { setN(220); setNoiseStd(0.2); setMethod("none"); }}>Reset</Button>
            </div>
          </div>

          <div className="grid md:grid-cols-2 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: t, y: noisy, type: "scatter", mode: "lines", name: "Noisy", line: { color: "#f97316" } },
                    { x: t, y: clean, type: "scatter", mode: "lines", name: "Clean", line: { color: theme.accent2 } },
                  ]}
                  layout={{
                    title: "Noisy signal vs Clean",
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
                    { x: t, y: derivRaw, type: "scatter", mode: "lines", name: "deriv raw", line: { color: "#ef4444" } },
                    { x: t, y: derivMov, type: "scatter", mode: "lines", name: "deriv movavg", line: { color: theme.accent } },
                    { x: t, y: derivSG, type: "scatter", mode: "lines", name: "deriv savgol", line: { color: theme.accent2 } },
                  ]}
                  layout={{
                    title: "Derivative estimates (various smoothing)",
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
  title="Exercise 23.4.1 — Noise amplification"
  body={`Show that as noise increases, derivative noise increases roughly like $\\mathcal{O}(1/h)$ for finite differences. Use experiments varying noise and $h$.`}
  solution={`Differentiate noisy signal and plot derivative variance vs noise amplitude for fixed $h$. Observe amplification roughly proportional to $\\sim 1/h^2$ depending on the finite difference scheme.`}
/>

        </CardContent>
      </Card>
    </motion.div>
  );
}

/* -------- 23.5 Partial Derivatives -------- */

function Section235() {
  const [nx, setNx] = useState(80);
  const [ny, setNy] = useState(60);

  const xs = useMemo(() => linspace(-2, 2, nx), [nx]);
  const ys = useMemo(() => linspace(-2, 2, ny), [ny]);

  const Z = useMemo(() => {
    const out = new Array(ny);
    for (let j = 0; j < ny; j++) {
      out[j] = new Array(nx);
      for (let i = 0; i < nx; i++) {
        const xv = xs[i], yv = ys[j];
        out[j][i] = Math.sin(Math.PI * xv) * Math.cos(Math.PI * yv);
      }
    }
    return out;
  }, [nx, ny, xs, ys]);

  const { zx, zy } = useMemo(() => partialDerivatives2D(Z, xs[1] - xs[0], ys[1] - ys[0]), [Z, xs, ys]);

  const midY = Math.floor(ny / 2);
  const zxSlice = zx[midY];
  const zySlice = zy.map((r) => r[Math.floor(nx / 2)]);

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Layers className="w-5 h-5" />
            23.5 Partial Derivatives
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid md:grid-cols-3 gap-3">
            <Labeled label="nx">
              <Input value={nx} onChange={(e) => setNx(Number(e.target.value))} className="bg-zinc-800 text-white" />
            </Labeled>
            <Labeled label="ny">
              <Input value={ny} onChange={(e) => setNy(Number(e.target.value))} className="bg-zinc-800 text-white" />
            </Labeled>
            <div className="flex items-end">
              <Button 
className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black"
 onClick={() => { setNx(80); setNy(60); }}>Reset</Button>
            </div>
          </div>

          <div className="grid md:grid-cols-2 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { z: Z, x: xs, y: ys, type: "surface", colorscale: "Viridis", showscale: false },
                  ]}
                  layout={{
                    title: "Scalar field f(x,y)",
                    autosize: true,
                    responsive: true,
                    scene: { xaxis: { title: "x", color: "#e5e7eb" }, yaxis: { title: "y", color: "#e5e7eb" }, zaxis: { title: "f", color: "#e5e7eb" } },
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
                    { x: xs, y: zxSlice, type: "scatter", mode: "lines", name: "∂f/∂x slice", line: { color: theme.accent } },
                    { x: ys, y: zySlice, type: "scatter", mode: "lines", name: "∂f/∂y slice", line: { color: theme.accent2 } },
                  ]}
                  layout={{
                    title: "Partial derivative slices",
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
              <Formula>{`Central differences for partials:\n\n$$\\frac{\\partial f}{\\partial x}\\approx \\frac{f_{i+1,j} - f_{i-1,j}}{2\\Delta x},\\quad \\frac{\\partial f}{\\partial y}\\approx \\frac{f_{i,j+1} - f_{i,j-1}}{2\\Delta y}.$$`}</Formula>
            </div>
          </div>

          <CollapsibleExercise
            title="Exercise 23.5.1 — Boundary stencils"
            body={`Derive a second-order accurate one-sided finite-difference formula for the boundary: \\(f'(x_0)\\) using f_0,f_1,f_2.`}
            solution={`One-sided second-order: \\(f'(x_0) \\approx \\frac{-3f_0 + 4f_1 - f_2}{2\\Delta x}\\).`}
          />
        </CardContent>
      </Card>
    </motion.div>
  );
}

/* -------- 23.6 Integration/Differentiation with Packages (demo) -------- */

function Section236() {
  const [N, setN] = useState(160);
  const [noiseStd, setNoiseStd] = useState(0.05);
  const t = useMemo(() => linspace(0, 2 * Math.PI, N), [N]);
  const y = useMemo(() => t.map((tt) => Math.sin(3 * tt) + 0.3 * Math.cos(0.7 * tt) + gaussianNoise(noiseStd)), [t, noiseStd]);

  function polyFitAndDeriv(tarr, yarr, degree = 6) {
    const m = tarr.length;
    const A = new Array(m);
    for (let i = 0; i < m; i++) {
      A[i] = new Array(degree + 1);
      for (let p = 0; p <= degree; p++) A[i][p] = Math.pow(tarr[i], p);
    }
    const ATA = new Array(degree + 1).fill(0).map(() => new Array(degree + 1).fill(0));
    const ATb = new Array(degree + 1).fill(0);
    for (let r = 0; r <= degree; r++) {
      for (let c = 0; c <= degree; c++) {
        let s = 0;
        for (let i = 0; i < m; i++) s += A[i][r] * A[i][c];
        ATA[r][c] = s;
      }
      let sb = 0;
      for (let i = 0; i < m; i++) sb += A[i][r] * yarr[i];
      ATb[r] = sb;
    }
    const coeffs = solveLinearSystem(ATA, ATb);
    if (!coeffs) return { yfit: new Array(m).fill(0), yderiv: new Array(m).fill(0) };
    const yfit = tarr.map((tt) => {
      let s = 0;
      for (let p = 0; p <= degree; p++) s += coeffs[p] * Math.pow(tt, p);
      return s;
    });
    const yderiv = tarr.map((tt) => {
      let s = 0;
      for (let p = 1; p <= degree; p++) s += p * coeffs[p] * Math.pow(tt, p - 1);
      return s;
    });
    return { yfit, yderiv };
  }

  const { yfit, yderiv } = useMemo(() => polyFitAndDeriv(t, y, 8), [t, y]);

  function simpsonIntegrate(xarr, yarr) {
    const n = xarr.length;
    if (n < 3) return 0;
    const h = (xarr[n - 1] - xarr[0]) / (n - 1);
    let s = yarr[0] + yarr[n - 1];
    for (let i = 1; i < n - 1; i++) s += (i % 2 === 0 ? 2 : 4) * yarr[i];
    return (h / 3) * s;
  }

  const integralNumeric = useMemo(() => simpsonIntegrate(t, y), [t, y]);
  const integralFit = useMemo(() => simpsonIntegrate(t, yfit), [t, yfit]);

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <FileText className="w-5 h-5" />
            23.6 Numerical Integration / Differentiation (Package demo)
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid md:grid-cols-3 gap-3">
            <Labeled label="N">
              <Input value={N} onChange={(e) => setN(Number(e.target.value))} className="bg-zinc-800 text-white"/>
            </Labeled>
            <Labeled label="Noise std">
              <Input value={noiseStd} onChange={(e) => setNoiseStd(Number(e.target.value))} className="bg-zinc-800 text-white"/>
            </Labeled>
            <div className="flex items-end">
              <Button 
className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black"
 onClick={() => { setN(160); setNoiseStd(0.05); }}>Reset</Button>
            </div>
          </div>

          <div className="grid md:grid-cols-2 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: t, y: y, type: "scatter", mode: "markers", name: "samples", marker: { color: theme.warn } },
                    { x: t, y: yfit, type: "scatter", mode: "lines", name: "poly fit", line: { color: theme.accent2 } },
                  ]}
                  layout={{
                    title: "Fitting for smoothing",
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
                    { x: t, y: yderiv, type: "scatter", mode: "lines", name: "derivative from fit", line: { color: theme.accent } },
                  ]}
                  layout={{
                    title: "Derivative of fitted polynomial",
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

              <div className="mt-2">
                <Summary items={[
                  { label: "Integral (Simpson)", value: integralNumeric.toFixed(6) },
                  { label: "Integral from fit", value: integralFit.toFixed(6) },
                  { label: "Difference", value: (integralNumeric - integralFit).toExponential(3) },
                ]} />
              </div>
            </div>
          </div>

          <CollapsibleExercise
            title="Exercise 23.6.1 — Fit-based derivatives"
            body={`Compare derivative estimates from local finite differences and from differentiating a smoothing spline or polynomial fit. Test varying degrees and smoothing parameters.`}
            solution={`Fit-based derivatives are smoother and reduce noise amplification, but can introduce bias if model is not flexible enough.`}
          />
        </CardContent>
      </Card>
    </motion.div>
  );
}

/* ============================
   Chapter documentation panel
   (mirror of Chapter7 style)
   ============================ */

function ChapterDocs() {
  const [open, setOpen] = useState("23.1 High-Accuracy Differentiation Formulas");

  return (
    <div className="w-full mt-6 bg-zinc-900 border border-zinc-700 rounded-xl shadow-lg overflow-hidden">
      <div className="p-3 border-b border-zinc-700 flex items-center justify-between">
        <div className="flex items-center gap-3">
          <BookOpen className="w-5 h-5 text-cyan-400" />
          <div>
            <div className="text-zinc-100 font-semibold">Chapter 23 — Documentation</div>
            <div className="text-zinc-400 text-xs">Theory, formulas, and problems (KaTeX rendered)</div>
          </div>
        </div>
        <div className="text-zinc-400 text-xs">Interactive demos & notes</div>
      </div>

      <div className="grid grid-cols-1 md:grid-cols-4">
        <div className="col-span-1 border-r border-zinc-800">
          {Object.keys(docsText).map((k) => (
            <button key={k} onClick={() => setOpen(k)} className={`w-full p-3 text-left hover:bg-zinc-800/40 flex items-center justify-between ${open === k ? "bg-zinc-800/20" : ""}`}>
              <div className="flex items-center gap-2">
                <List className="w-4 h-4 text-zinc-300" />
                <div className="text-zinc-100">{k}</div>
              </div>
              {open === k ? <ChevronDown className="w-4 h-4 text-zinc-400" /> : <ChevronRight className="w-4 h-4 text-zinc-400" />}
            </button>
          ))}
        </div>

        <div className="col-span-3 p-4 overflow-auto" style={{ maxHeight: 460 }}>
          <motion.h3 initial={{ opacity: 0, y: -6 }} animate={{ opacity: 1, y: 0 }} className="text-xl font-semibold text-zinc-100 mb-3">
            {open}
          </motion.h3>
          <div className="text-zinc-300 space-y-3 text-sm leading-relaxed">
            {docsText[open].map((para, i) => (
              <div key={i}>
                <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>
                  {para}
                </ReactMarkdown>
              </div>
            ))}
            <div className="mt-3 text-zinc-400 text-xs">References: Burden & Faires; Stoer & Bulirsch; Press et al. (Numerical Recipes).</div>
          </div>
        </div>
      </div>
    </div>
  );
}

/* ============================
   Problems panel (end of page)
   ============================ */

function ProblemsPanel() {
  return (
    <div className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
      <div className="flex items-center justify-between mb-3">
        <div className="text-xl text-cyan-300 font-semibold flex items-center gap-2"><FileText className="w-5 h-5" /> Problems</div>
        <Button variant="ghost" size="sm" onClick={() => window.scrollTo({ top: 0, behavior: "smooth" })}><RefreshCcw className="w-4 h-4" /> Top</Button>
      </div>

      <div className="grid md:grid-cols-2 gap-3">
        <CollapsibleExercise
          title="Problem 23.1 — Convergence plots"
          body={`For \\(f(x)=\\sin x\\) compute derivative at x=1 using central and five-point formulas for h=10^{-k}, k=1..6. Plot error vs h on log-log scale and estimate slopes.`}
          solution={`Central difference slope ≈ -2; five-point slope ≈ -4 until rounding error.`}
        />
        <CollapsibleExercise
          title="Problem 23.2 — Richardson accumulation"
          body={`Implement Richardson iteratively to boost accuracy from O(h^2) to higher order; test on e^x at x=0.5. How many levels are practical?`}
          solution={`Usually 1–2 levels are practical; more levels may amplify round-off with tiny h.`}
        />
        <CollapsibleExercise
          title="Problem 23.3 — Unequally spaced real data"
          body={`Given irregularly sampled sensor data with noise, compute derivatives using local quadratic fits and quantify robustness with bootstrap.`}
          solution={`Local quadratic fits give robust estimates; bootstrap yields CI for derivative estimates.`}
        />
        <CollapsibleExercise
          title="Problem 23.4 — Partial derivatives"
          body={`Compute gradient magnitude on a sampled 2D temperature field using central differences and visualize heatmap.`}
          solution={`Compute ∂f/∂x, ∂f/∂y, then gradient magnitude sqrt(fx^2 + fy^2).`}
        />
      </div>
    </div>
  );
}

/* ============================
   Page assembly
   ============================ */

export default function Chapter23() {
  return (
    <div className={`p-6 space-y-6 ${theme.bg}`}>
      <motion.div initial={{ opacity: 0, y: -10 }} animate={{ opacity: 1, y: 0 }}>
        <h1 className="text-3xl font-bold" style={{ color: theme.accent2 }}>
          Numerical Differentiation
        </h1>
        <div className="text-zinc-400">High-accuracy formulas, Richardson extrapolation, unequal spacing, noise handling, partial derivatives, and software tips.</div>
      </motion.div>

      <Section231 />
      <Section232 />
      <Section233 />
      <Section234 />
      <Section235 />
      <Section236 />

      <ChapterDocs />

      <ProblemsPanel />
      <BottomBar/>
    </div>
  );
}
