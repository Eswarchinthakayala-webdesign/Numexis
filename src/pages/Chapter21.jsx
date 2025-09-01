// src/pages/Chapter21.jsx
// ======================================================================
// Chapter 21 — Numerical Integration (Responsive, full-featured)
// - Trapezoidal, Simpson (1/3 & 3/8), Midpoint (open), Unequal segments,
//   Open Newton–Cotes, Multiple integrals (double integrals)
// - Responsive Plotly plots (autosize + useResizeHandler + parent sizing)
// - CSV / PNG export for plots & tabular data
// - Robust safe function parsing (JS Math namespace) with fallbacks
// - Mobile-first layout, accessible controls, overflow protection
// - Detailed inline documentation & helpful UI components
// ======================================================================

import React, { useMemo, useState, useRef, useEffect } from "react";
import Plot from "react-plotly.js";
import { motion } from "framer-motion";
import ReactMarkdown from "react-markdown";
import remarkMath from "remark-math";
import rehypeKatex from "rehype-katex";
import "katex/dist/katex.min.css";

import { Card, CardHeader, CardContent, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Textarea } from "@/components/ui/textarea";
import { Button } from "@/components/ui/button";

import {
  FunctionSquare,
  Percent,
  Layers,
  Grid,
  ArrowRightSquare,
  Download,
  ImageIcon,
  Calculator,
  Zap,
  BookOpen,
  Activity,
  ChevronDown,
  ChevronRight,
} from "lucide-react";
import BottomBar from "../components/BottomBar";

// ======================================================================
// Theme & small helpers
// ======================================================================
const theme = {
  bg: "bg-zinc-950 text-zinc-200",
  panel: "bg-zinc-900/60 border border-zinc-700",
  accent: "#22d3ee",
  accent2: "#34d399",
  monospace: "font-mono",
};

const fadeUp = {
  initial: { opacity: 0, y: 8 },
  animate: { opacity: 1, y: 0 },
  transition: { duration: 0.24 },
};

function fmt(v, p = 6) {
  if (!Number.isFinite(v)) return "NaN";
  return Number(v).toFixed(p).replace(/(?:\.0+|(\.\d+?)0+)$/, "$1");
}

function clamp(v, min, max) {
  return Math.max(min, Math.min(max, v));
}

// ----------------------------------------------------------------------
// UI Microcomponents
// ----------------------------------------------------------------------
function Labeled({ label, children, desc }) {
  return (
    <div className="min-w-0">
      <div className="text-xs text-zinc-400 mb-1">{label}</div>
      {children}
      {desc && <div className="text-[11px] text-zinc-500 mt-1">{desc}</div>}
    </div>
  );
}

function Equation({ children, color = theme.accent }) {
  return (
    <div className="my-2">
      <div className="px-3 py-2 rounded-lg border text-gray-300 border-zinc-700 bg-zinc-900/80 overflow-x-auto">
        <div
          className="text-xs uppercase tracking-widest mb-1"
          style={{ color }}
        >
          formula
        </div>
        <div className="min-w-fit">
          <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>
            {children}
          </ReactMarkdown>
        </div>
      </div>
    </div>
  );
}


function Exercise({ title, body, solution }) {
  const [open, setOpen] = useState(false);
  return (
    <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-3">
      <div className="flex items-center justify-between gap-2">
        <div className="text-zinc-100 font-medium">{title}</div>
        <Button variant="ghost" className="bg-white cursor-pointer py-2 px-2 rounded text-black flex items-center hover:bg-gray-200" size="sm" onClick={() => setOpen((o) => !o)}>
          {open ? <ChevronDown className="w-4 h-4" /> : <ChevronRight className="w-4 h-4" /> } 
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

// ======================================================================
// Utility: Safe function creators (1D & 2D), CSV/PNG export, linspace
// ======================================================================

function makeFunctionSafe1D(expr) {
  // Create a safe numeric function from user input string using JS Math.
  // Replace ^ with ** to allow '^' shorthand for exponent.
  try {
    const clean = String(expr).replace(/\^/g, "**");
    // eslint-disable-next-line no-new-func
    const fn = new Function("x", `with (Math) { return (${clean}); }`);
    // return wrapper that always returns Number or NaN
    return (x) => {
      try {
        const v = fn(x);
        return Number(v);
      } catch {
        return NaN;
      }
    };
  } catch {
    return () => NaN;
  }
}

function makeFunctionSafe2D(expr) {
  try {
    const clean = String(expr).replace(/\^/g, "**");
    // eslint-disable-next-line no-new-func
    const fn = new Function("x", "y", `with (Math) { return (${clean}); }`);
    return (x, y) => {
      try {
        const v = fn(x, y);
        return Number(v);
      } catch {
        return NaN;
      }
    };
  } catch {
    return () => NaN;
  }
}

function downloadCSV(filename, rows) {
  const csv = rows.map((r) => r.map((c) => `"${String(c).replace(/"/g, '""')}"`).join(",")).join("\n");
  const blob = new Blob([csv], { type: "text/csv;charset=utf-8;" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  setTimeout(() => {
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  }, 150);
}

async function exportPlotPNG(divId, filename = "plot.png") {
  try {
    if (window.Plotly && document.getElementById(divId)) {
      // choose a reasonably large raster size for export
      const imgData = await window.Plotly.toImage(document.getElementById(divId), {
        format: "png",
        width: 1200,
        height: 800,
      });
      const a = document.createElement("a");
      a.href = imgData;
      a.download = filename;
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
    } else {
      alert("Export failed: Plotly not found or element missing.");
    }
  } catch (err) {
    console.error(err);
    alert("Export failed: " + String(err));
  }
}

function linspace(a, b, n) {
  n = Math.max(1, Math.floor(n));
  if (n === 1) return [a];
  const arr = new Array(n);
  const dx = (b - a) / (n - 1);
  for (let i = 0; i < n; i++) arr[i] = a + i * dx;
  return arr;
}

// ======================================================================
// Numerical algorithms (robust, well-tested implementations)
// ======================================================================

// Composite trapezoid
function compositeTrapezoid(f, a, b, m) {
  if (!isFinite(a) || !isFinite(b) || m < 1) return NaN;
  const h = (b - a) / m;
  let s = 0.5 * (f(a) + f(b));
  for (let i = 1; i < m; i++) s += f(a + i * h);
  return s * h;
}

// Composite Simpson 1/3 (m even)
function compositeSimpson13(f, a, b, m) {
  if (!isFinite(a) || !isFinite(b) || m < 2) return NaN;
  if (m % 2 === 1) m += 1;
  const h = (b - a) / m;
  let s = f(a) + f(b);
  for (let i = 1; i < m; i++) {
    s += f(a + i * h) * (i % 2 === 0 ? 2 : 4);
  }
  return (h / 3) * s;
}

// Composite Simpson 3/8 (m multiple of 3)
function compositeSimpson38(f, a, b, m) {
  if (!isFinite(a) || !isFinite(b) || m < 3) return NaN;
  while (m % 3 !== 0) m += 1;
  const h = (b - a) / m;
  let s = f(a) + f(b);
  for (let i = 1; i < m; i++) s += f(a + i * h) * (i % 3 === 0 ? 2 : 3);
  return (3 * h / 8) * s;
}

// Open Newton-Cotes (midpoint & simple open 2-point)
function compositeOpenNewtonCotes(f, a, b, n, order = 0) {
  if (!isFinite(a) || !isFinite(b) || n < 1) return NaN;
  const h = (b - a) / n;
  let sum = 0;
  if (order === 0) {
    // midpoint rule: sample midpoints
    for (let i = 0; i < n; i++) {
      const xm = a + (i + 0.5) * h;
      sum += f(xm);
    }
    return sum * h;
  } else if (order === 1) {
    // 2-point open: sample at 1/3 and 2/3 of each panel
    for (let i = 0; i < n; i++) {
      const x0 = a + i * h;
      sum += f(x0 + h / 3) + f(x0 + (2 * h) / 3);
    }
    return (h / 2) * sum;
  } else {
    // fallback to midpoint
    for (let i = 0; i < n; i++) {
      const xm = a + (i + 0.5) * h;
      sum += f(xm);
    }
    return sum * h;
  }
}

// Unequal trapezoid segments given nodes array
function trapezoidUnequal(f, xs) {
  if (!Array.isArray(xs) || xs.length < 2) return NaN;
  let s = 0;
  for (let i = 0; i < xs.length - 1; i++) {
    const x0 = xs[i];
    const x1 = xs[i + 1];
    s += 0.5 * (x1 - x0) * (f(x0) + f(x1));
  }
  return s;
}

// Romberg integration table
function rombergIntegration(f, a, b, maxK = 6) {
  if (!isFinite(a) || !isFinite(b)) return [];
  const R = Array.from({ length: maxK }, () => []);
  for (let k = 0; k < maxK; k++) {
    const n = 2 ** k;
    const T = compositeTrapezoid(f, a, b, n);
    R[k][0] = T;
    for (let j = 1; j <= k; j++) {
      R[k][j] = R[k][j - 1] + (R[k][j - 1] - R[k - 1][j - 1]) / (4 ** j - 1);
    }
  }
  return R;
}

// Adaptive Simpson's method (recursive)
function adaptiveSimpson(f, a, b, tol = 1e-8, maxDepth = 20) {
  function simpson(a0, b0) {
    const c0 = 0.5 * (a0 + b0);
    return ((b0 - a0) / 6) * (f(a0) + 4 * f(c0) + f(b0));
  }
  function recurse(a0, b0, eps, S, depth) {
    const c0 = 0.5 * (a0 + b0);
    const Sleft = simpson(a0, c0);
    const Sright = simpson(c0, b0);
    const err = Sleft + Sright - S;
    if (depth <= 0 || Math.abs(err) < 15 * eps) {
      return Sleft + Sright + err / 15;
    }
    return recurse(a0, c0, eps / 2, Sleft, depth - 1) + recurse(c0, b0, eps / 2, Sright, depth - 1);
  }
  const S = simpson(a, b);
  return recurse(a, b, tol, S, maxDepth);
}

// Gauss-Legendre nodes & weights (for general n)
function legendreP(n, x) {
  if (n === 0) return 1;
  if (n === 1) return x;
  let p0 = 1;
  let p1 = x;
  for (let k = 2; k <= n; k++) {
    const pk = ((2 * k - 1) * x * p1 - (k - 1) * p0) / k;
    p0 = p1;
    p1 = pk;
  }
  return p1;
}
function legendrePDeriv(n, x) {
  if (n === 0) return 0;
  return (n / (x * x - 1)) * (x * legendreP(n, x) - legendreP(n - 1, x));
}
function gaussNodesWeights(n, tol = 1e-14) {
  // Returns { nodes: [], weights: [] } on [-1,1]
  const nodes = new Array(n);
  const weights = new Array(n);
  const m = Math.floor((n + 1) / 2);
  for (let i = 0; i < m; i++) {
    // initial guess (Chebyshev)
    let x = Math.cos(Math.PI * (i + 0.75) / (n + 0.5));
    let dx = 1;
    let its = 0;
    while (Math.abs(dx) > tol && its < 100) {
      const p = legendreP(n, x);
      const dp = legendrePDeriv(n, x);
      dx = -p / dp;
      x += dx;
      its += 1;
    }
    nodes[i] = x;
    nodes[n - 1 - i] = -x;
    const dpFinal = legendrePDeriv(n, x);
    const w = 2 / ((1 - x * x) * dpFinal * dpFinal);
    weights[i] = w;
    weights[n - 1 - i] = w;
  }
  return { nodes, weights };
}
function gaussLegendreQuad(f, a, b, n) {
  if (!isFinite(a) || !isFinite(b) || n <= 0) return NaN;
  const { nodes, weights } = gaussNodesWeights(n);
  let sum = 0;
  const mid = 0.5 * (a + b);
  const half = 0.5 * (b - a);
  for (let i = 0; i < n; i++) {
    const xi = mid + half * nodes[i];
    sum += weights[i] * f(xi);
  }
  return half * sum;
}

// Improper integrals: transforms for infinite intervals
function improperIntegralInfinite(f, a, b, tol = 1e-8) {
  if (!isFinite(a) && !isFinite(b)) {
    // both infinite: x = tan(pi*(t-0.5))
    const transform = (t) => {
      const x = Math.tan(Math.PI * (t - 0.5));
      const dxdt = Math.PI / Math.cos(Math.PI * (t - 0.5)) ** 2;
      return f(x) * dxdt;
    };
    return adaptiveSimpson(transform, 0, 1, tol);
  }
  if (!isFinite(b)) {
    // upper infinite: x = a + t/(1-t)
    const transform = (t) => {
      const x = a + t / (1 - t);
      const dxdt = 1 / (1 - t) ** 2;
      return f(x) * dxdt;
    };
    return adaptiveSimpson(transform, 0, 1 - 1e-12, tol);
  }
  if (!isFinite(a)) {
    // lower infinite: x = b - t/(1-t)
    const transform = (t) => {
      const x = b - t / (1 - t);
      const dxdt = 1 / (1 - t) ** 2;
      return f(x) * dxdt;
    };
    return adaptiveSimpson(transform, 0, 1 - 1e-12, tol);
  }
  // finite endpoints
  return adaptiveSimpson(f, a, b, tol);
}

// Double integral tensor-product composite trapezoid
function doubleIntegralTensor(f, ax, bx, ay, by, nx = 40, ny = 40) {
  if (!isFinite(ax) || !isFinite(bx) || !isFinite(ay) || !isFinite(by)) return NaN;
  const hx = (bx - ax) / nx;
  const hy = (by - ay) / ny;
  let sum = 0;
  for (let i = 0; i <= nx; i++) {
    const x = ax + i * hx;
    const wx = (i === 0 || i === nx) ? 0.5 : 1;
    for (let j = 0; j <= ny; j++) {
      const y = ay + j * hy;
      const wy = (j === 0 || j === ny) ? 0.5 : 1;
      sum += wx * wy * f(x, y);
    }
  }
  return sum * hx * hy;
}

// ======================================================================
// Responsive Plot wrapper (ensures plots never overflow)
// ======================================================================
function ResponsivePlot({ id, data, layout = {}, config = {}, style = {}, title }) {
  // The outer div determines size; Plot uses useResizeHandler and autosize.
  return (
    <div id={id} className="rounded-xl border border-zinc-700 p-2 bg-zinc-900/40 overflow-hidden w-full" style={{ maxWidth: "100%" }}>
      <Plot
        data={data}
        layout={{
          title,
          paper_bgcolor: "rgba(0,0,0,0)",
          plot_bgcolor: "rgba(0,0,0,0)",
          font: { color: "#e5e7eb" },
          autosize: true,
          margin: { l: 50, r: 20, t: 50, b: 45 },
          legend: { orientation: "h", x: 0, y: -0.18 },
          ...layout,
        }}
        config={{ responsive: true, displayModeBar: true, ...config }}
        useResizeHandler
        style={{ width: "100%", height: "clamp(260px, 52vh, 620px)", ...style }}
      />
    </div>
  );
}

// ======================================================================
// Section implementations (21.1 - 21.5) — expanded with responsive UI
// ======================================================================

// ---------------- Section 21.1: Trapezoidal Rule -----------------------
function Section211() {
  const [expr, setExpr] = useState("Math.sin(x)");
  const [a, setA] = useState("0");
  const [b, setB] = useState(String(Math.PI));
  const [n, setN] = useState(8);
  const [showFill, setShowFill] = useState(true);
  const plotId = "plot-21-1";

  const f = useMemo(() => makeFunctionSafe1D(expr), [expr]);
  const aNum = Number(a);
  const bNum = Number(b);
  const nNum = Math.max(1, Math.floor(Number(n)));

  const approx = useMemo(() => compositeTrapezoid(f, aNum, bNum, nNum), [f, aNum, bNum, nNum]);

  const xsPlot = useMemo(() => {
    const L = 400;
    if (!isFinite(aNum) || !isFinite(bNum) || bNum === aNum) return [];
    return linspace(aNum, bNum, L);
  }, [aNum, bNum]);

  const ysPlot = useMemo(() => xsPlot.map((x) => f(x)), [xsPlot, f]);

  const subdivisions = useMemo(() => {
    const arr = [];
    for (let i = 0; i <= nNum; i++) arr.push(aNum + (i / nNum) * (bNum - aNum));
    return arr;
  }, [aNum, bNum, nNum]);

  const areaTraces = useMemo(() => {
    return subdivisions.slice(0, subdivisions.length - 1).map((x0, i) => {
      const x1 = subdivisions[i + 1];
      return {
        x: [x0, x0, x1, x1],
        y: [0, f(x0), f(x1), 0],
        fill: "toself",
        type: "scatter",
        mode: "lines",
        fillcolor: "rgba(34,211,238,0.12)",
        line: { color: "rgba(34,211,238,0.3)" },
        showlegend: i === 0,
        name: i === 0 ? "Trapezoids" : undefined,
      };
    });
  }, [subdivisions, f]);

  function handleExportCSV() {
    const rows = [["x", "f(x)"]];
    for (let i = 0; i < xsPlot.length; i++) rows.push([xsPlot[i], ysPlot[i]]);
    downloadCSV("21_1_trapezoid.csv", rows);
  }
  function handleExportPNG() {
    exportPlotPNG(plotId, "21_1_trapezoid.png");
  }

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <FunctionSquare className="w-5 h-5" />
            21.1 Trapezoidal Rule — Visualizer
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-4 gap-3">
            <Labeled label="f(x) (JavaScript Math namespace)">
              <Textarea value={expr} onChange={(e) => setExpr(e.target.value)} className="h-24 w-full text-white" />
            </Labeled>
            <Labeled  label="a"><Input className="text-white" value={a} onChange={(e) => setA(e.target.value)} /></Labeled>
            <Labeled label="b"><Input className="text-white" value={b} onChange={(e) => setB(e.target.value)} /></Labeled>
            <Labeled label="n (subintervals)"><Input className="text-white" value={n} onChange={(e) => setN(e.target.value)} /></Labeled>
          </div>

          <div className="flex items-center gap-3 flex-wrap">
            <div className="text-xs text-zinc-400">Show fills:</div>
            <Button className={`px-2 py-1 rounded ${showFill ? "bg-emerald-600 cursor-pointer hover:bg-emerald-500 text-white" : "bg-zinc-500 cursor-pointer hover:bg-zinc-600 text-zinc-200"}`} onClick={() => setShowFill((s) => !s)}>{showFill ? "On" : "Off"}</Button>
            <div className="text-xs text-zinc-400 ml-3">Interactive plot is responsive — resize your window to test.</div>
          </div>

         <Equation>{`Composite trapezoidal rule: $\\int_a^b f(x) \\, dx \\approx \\frac{b-a}{2n} \\Big[f(a) + 2\\sum_{i=1}^{n-1} f(x_i) + f(b)\\Big]$.`}</Equation>


          <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
            <ResponsivePlot
              id={plotId}
              data={[
                { x: xsPlot, y: ysPlot, mode: "lines", type: "scatter", name: "f(x)", line: { color: theme.accent } },
                ...(showFill ? areaTraces : []),
              ]}
              title={`Trapezoidal visualization (approx ≈ ${fmt(approx, 10)})`}
            />

            <div className="rounded-xl border border-zinc-700 p-3 bg-zinc-900/40">
              <div className="text-zinc-300">Composite trapezoidal approximation:</div>
              <div className="text-2xl mt-2 text-emerald-300 font-mono">{fmt(approx, 10)}</div>
              <div className="mt-3 flex gap-2 flex-wrap">
                <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={handleExportCSV} variant="ghost"><Download className="w-4 h-4" /> CSV</Button>
                <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={handleExportPNG} variant="ghost"><ImageIcon className="w-4 h-4" /> PNG</Button>
              </div>

              <div className="mt-4 text-xs text-zinc-400">
                <Exercise
                  title="Exercise 21.1.1"
                  body="Use the trapezoidal rule with n=8 to approximate ∫_0^π sin(x) dx. Compare to exact value 2."
                  solution="With n=8 the approximation will be close to 2; refine n to see error decay O(h^2)."
                />
              </div>
            </div>
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ---------------- Section 21.2: Simpson's Rules ------------------------
function Section212() {
  const [expr, setExpr] = useState("Math.exp(-x*x)");
  const [a, setA] = useState("0");
  const [b, setB] = useState("1");
  const [n, setN] = useState(12);
  const [mode, setMode] = useState("1/3"); // "1/3" or "3/8"
  const plotId = "plot-21-2";

  const f = useMemo(() => makeFunctionSafe1D(expr), [expr]);
  const aNum = Number(a);
  const bNum = Number(b);
  let nNum = Math.max(1, Math.floor(Number(n)));

  const approx = useMemo(() => {
    if (mode === "1/3") {
      if (nNum % 2 === 1) nNum += 1;
      return compositeSimpson13(f, aNum, bNum, nNum);
    } else {
      while (nNum % 3 !== 0) nNum += 1;
      return compositeSimpson38(f, aNum, bNum, nNum);
    }
  }, [f, aNum, bNum, nNum, mode]);

  const xsPlot = useMemo(() => {
    const L = 400;
    if (!isFinite(aNum) || !isFinite(bNum) || bNum === aNum) return [];
    return linspace(aNum, bNum, L);
  }, [aNum, bNum]);

  const ysPlot = useMemo(() => xsPlot.map((x) => f(x)), [xsPlot, f]);

  function handleExportCSV() {
    const rows = [["x", "f(x)"], ...xsPlot.map((x, i) => [x, ysPlot[i]])];
    downloadCSV(`21_2_simpson_${mode}.csv`, rows);
  }
  function handleExportPNG() {
    exportPlotPNG(plotId, `21_2_simpson_${mode}.png`);
  }

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Layers className="w-5 h-5" />
            21.2 Simpson's Rules — 1/3 & 3/8
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-4 gap-3">
            <Labeled label="f(x)"><Textarea value={expr} onChange={(e) => setExpr(e.target.value)} className="h-24 w-full text-white" /></Labeled>
            <Labeled label="a"><Input className="text-white" value={a} onChange={(e) => setA(e.target.value)} /></Labeled>
            <Labeled label="b"><Input className="text-white" value={b} onChange={(e) => setB(e.target.value)} /></Labeled>
            <Labeled label="n (subintervals)"><Input className="text-white" value={n} onChange={(e) => setN(e.target.value)} /></Labeled>
          </div>

          <div className="flex items-center gap-2 flex-wrap">
            <Button className={`px-3 py-1 rounded ${mode === "1/3" ? "bg-emerald-600 hover:bg-emerald-500 cursor-pointer text-white" : "bg-zinc-500 hover:bg-zinc-600 cursor-pointer text-zinc-200"}`} onClick={() => setMode("1/3")}>Simpson 1/3</Button>
            <Button className={`px-3 py-1 rounded ${mode === "3/8" ? "bg-emerald-600 hover:bg-emerald-500 cursor-pointer text-white" : "bg-zinc-500 hover:bg-zinc-600 cursor-pointer text-zinc-200"}`} onClick={() => setMode("3/8")}>Simpson 3/8</Button>
            <div className="text-xs text-zinc-400 ml-3">Note: 1/3 needs n even; 3/8 needs n multiple of 3. The UI adjusts automatically.</div>
          </div>

        <Equation>{`Simpson $1/3$: $\\int_a^b f(x) \\, dx \\approx \\frac{h}{3} \\Big[f(a) + 4 f(a+h) + 2 f(a+2h) + \\dots + f(b)\\Big]$`}</Equation>


          <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
            <ResponsivePlot id={plotId} data={[{ x: xsPlot, y: ysPlot, mode: "lines", type: "scatter", name: "f(x)", line: { color: theme.accent } }]} title={`Simpson ${mode} approximation (≈ ${fmt(approx, 10)})`} />

            <div className="rounded-xl border border-zinc-700 p-3 bg-zinc-900/40">
              <div className="text-zinc-300">Approximation:</div>
              <div className="text-2xl mt-2 text-emerald-300 font-mono">{fmt(approx, 10)}</div>
              <div className="mt-3 flex gap-2 flex-wrap">
                <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={handleExportCSV} variant="ghost"><Download className="w-4 h-4" /> CSV</Button>
                <Button  className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={handleExportPNG} variant="ghost"><ImageIcon className="w-4 h-4" /> PNG</Button>
              </div>

              <div className="mt-4 text-xs text-zinc-400">
                <Exercise title="Exercise 21.2.1" body="Compare Simpson 1/3 and 3/8 on ∫_0^2 e^{-x^2} dx with n=12. Which is more accurate?" solution="Both converge; 3/8 requires n multiple of 3. For smooth functions both give very similar results when using compatible n." />
              </div>
            </div>
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ---------------- Section 21.3: Unequal Segments ---------------------
function Section213() {
  const [expr, setExpr] = useState("Math.log(x+1)");
  const [xNodesStr, setXNodesStr] = useState("0 0.5 1.5 2.3 3");
  const plotId = "plot-21-3";

  const f = useMemo(() => makeFunctionSafe1D(expr), [expr]);

  const xsNodes = useMemo(() => {
    try {
      return xNodesStr
        .split(/[,\s]+/)
        .map(Number)
        .filter((v) => !Number.isNaN(v))
        .sort((a, b) => a - b);
    } catch {
      return [];
    }
  }, [xNodesStr]);

  const approx = useMemo(() => (xsNodes.length >= 2 ? trapezoidUnequal(f, xsNodes) : NaN), [f, xsNodes]);

  const xsPlot = useMemo(() => {
    if (xsNodes.length < 2) return [];
    const a = xsNodes[0];
    const b = xsNodes[xsNodes.length - 1];
    return linspace(a, b, 400);
  }, [xsNodes]);

  const ysPlot = useMemo(() => xsPlot.map((x) => f(x)), [xsPlot, f]);

  function handleExportCSV() {
    const rows = [["x", "f(x)"], ...xsPlot.map((x, i) => [x, ysPlot[i]])];
    downloadCSV("21_3_unequal.csv", rows);
  }
  function handleExportPNG() {
    exportPlotPNG(plotId, "21_3_unequal.png");
  }

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Grid className="w-5 h-5" />
            21.3 Integration with Unequal Segments
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 md:grid-cols-3 gap-3">
            <Labeled label="f(x)"><Input className="text-white" value={expr} onChange={(e) => setExpr(e.target.value)} /></Labeled>
            <Labeled label="x nodes (space or comma separated)">
              <Textarea value={xNodesStr} onChange={(e) => setXNodesStr(e.target.value)} className="h-24 w-full text-white" />
            </Labeled>
            <div>
              <div className="text-xs text-zinc-400 mb-1">Result</div>
              <div className="text-2xl text-emerald-300">{fmt(approx, 8)}</div>
            </div>
          </div>

          <Equation>{`Trapezoidal rule with unequal segments: $\\sum_{i} \\frac{x_{i+1}-x_i}{2} \\Big[f(x_i) + f(x_{i+1})\\Big]$.`}</Equation>


          <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
            <ResponsivePlot
              id={plotId}
              data={[
                { x: xsPlot, y: ysPlot, mode: "lines", name: "f(x)", line: { color: theme.accent } },
                { x: xsNodes, y: xsNodes.map((x) => f(x)), mode: "markers", name: "nodes", marker: { size: 8, color: theme.accent2 } },
                ...xsNodes.slice(0, xsNodes.length - 1).map((x0, i) => {
                  const x1 = xsNodes[i + 1];
                  return {
                    x: [x0, x0, x1, x1],
                    y: [0, f(x0), f(x1), 0],
                    fill: "toself",
                    type: "scatter",
                    mode: "lines",
                    fillcolor: "rgba(52,211,153,0.12)",
                    line: { color: "rgba(52,211,153,0.3)" },
                    showlegend: i === 0,
                    name: i === 0 ? "Unequal traps" : undefined,
                  };
                }),
              ]}
              title={`Unequal segments trapezoids (approx ≈ ${fmt(approx, 8)})`}
            />

            <div className="rounded-xl border border-zinc-700 p-3 bg-zinc-900/40 overflow-x-auto">
              <div className="flex gap-2">
                <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={handleExportCSV} variant="ghost"><Download className="w-4 h-4" /> CSV</Button>
                <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={handleExportPNG} variant="ghost"><ImageIcon className="w-4 h-4" /> PNG</Button>
              </div>

              <div className="mt-4 text-xs text-zinc-400">
                <Exercise
                  title="Exercise 21.3.1"
                  body="Integrate f(x)=e^{-x^2} on nonuniform partition [0,0.2,1,1.5,2] using trapezoid unequal segments. Compare to uniform trapezoid with same number of panels."
                  solution="Nonuniform partitions can concentrate nodes where f changes rapidly to improve accuracy; test by varying partitions."
                />
              </div>
            </div>
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ---------------- Section 21.4: Open formulas (Midpoint & open NC) --------------
function Section214() {
  const [expr, setExpr] = useState("Math.cos(x)");
  const [a, setA] = useState("0");
  const [b, setB] = useState(String(Math.PI / 2));
  const [n, setN] = useState(8);
  const [order, setOrder] = useState(0); // 0 midpoint, 1 two-point open
  const plotId = "plot-21-4";

  const f = useMemo(() => makeFunctionSafe1D(expr), [expr]);
  const aNum = Number(a);
  const bNum = Number(b);
  const nNum = Math.max(1, Math.floor(Number(n)));

  const approx = useMemo(() => compositeOpenNewtonCotes(f, aNum, bNum, nNum, Number(order)), [f, aNum, bNum, nNum, order]);

  const xsPlot = useMemo(() => {
    const L = 400;
    if (!isFinite(aNum) || !isFinite(bNum) || bNum === aNum) return [];
    return linspace(aNum, bNum, L);
  }, [aNum, bNum]);

  const ysPlot = useMemo(() => xsPlot.map((x) => f(x)), [xsPlot, f]);

  function handleExportCSV() {
    const rows = [["x", "f(x)"], ...xsPlot.map((x, i) => [x, ysPlot[i]])];
    downloadCSV("21_4_open_nc.csv", rows);
  }
  function handleExportPNG() {
    exportPlotPNG(plotId, "21_4_open_nc.png");
  }

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <ArrowRightSquare className="w-5 h-5" />
            21.4 Open Integration Formulas (Midpoint & Open NC)
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-4 gap-3">
            <Labeled label="f(x)"><Textarea value={expr} onChange={(e) => setExpr(e.target.value)} className="h-24 w-full text-white" /></Labeled>
            <Labeled label="a"><Input className="text-white" value={a} onChange={(e) => setA(e.target.value)} /></Labeled>
            <Labeled label="b"><Input className="text-white" value={b} onChange={(e) => setB(e.target.value)} /></Labeled>
            <Labeled label="n"><Input className="text-white" value={n} onChange={(e) => setN(e.target.value)} /></Labeled>
          </div>

          <div className="flex items-center gap-2">
            <button className={`px-3 py-1 rounded ${order === 0 ? "bg-emerald-600 cursor-pointer text-white" : "bg-zinc-800 cursor-pointer text-zinc-200"}`} onClick={() => setOrder(0)}>Midpoint</button>
            <button className={`px-3 py-1 rounded ${order === 1 ? "bg-emerald-600 cursor-pointer text-white" : "bg-zinc-800 cursor-pointer text-zinc-200"}`} onClick={() => setOrder(1)}>2-point open</button>
            <div className="text-xs text-zinc-400 ml-3">Open rules avoid endpoint evaluations (useful for endpoint singularities).</div>
          </div>

          <Equation>{`Open Newton-Cotes (midpoint): use interior points of each subinterval to avoid endpoints; midpoint approximates the integral by evaluating f at center of each panel.`}</Equation>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
            <ResponsivePlot id={plotId} data={[{ x: xsPlot, y: ysPlot, mode: "lines", name: "f(x)", line: { color: theme.accent } }]} title={`Open Newton-Cotes (order ${order}) approx ≈ ${fmt(approx, 8)}`} />

            <div className="rounded-xl border border-zinc-700 p-3 bg-zinc-900/40">
              <div className="text-zinc-300">Approximation:</div>
              <div className="text-2xl mt-2 text-emerald-300 font-mono">{fmt(approx, 10)}</div>
              <div className="mt-3 flex gap-2">
                <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={handleExportCSV} variant="ghost"><Download className="w-4 h-4" /> CSV</Button>
                <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={handleExportPNG} variant="ghost"><ImageIcon className="w-4 h-4" /> PNG</Button>
              </div>

              <div className="mt-4 text-xs text-zinc-400">
                <Exercise title="Exercise 21.4.1" body="Use the midpoint rule to approximate ∫_0^1 x^{-1/2} dx with n=10. Explain why midpoint is advantageous for endpoint singularity." solution="Midpoint avoids x=0 (singularity), sampling interior points yields finite evaluations; the overflow from the endpoint is avoided." />
              </div>
            </div>
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ---------------- Section 21.5: Multiple Integrals (Double integrals) ------------
function Section215() {
  const [expr, setExpr] = useState("Math.sin(x)*Math.cos(y)");
  const [ax, setAx] = useState("0");
  const [bx, setBx] = useState(String(Math.PI));
  const [ay, setAy] = useState("0");
  const [by, setBy] = useState(String(Math.PI / 2));
  const [nx, setNx] = useState(60);
  const [ny, setNy] = useState(60);
  const plotId = "plot-21-5";

  const f = useMemo(() => makeFunctionSafe2D(expr), [expr]);

  const axNum = Number(ax);
  const bxNum = Number(bx);
  const ayNum = Number(ay);
  const byNum = Number(by);
  const nxNum = Math.max(4, Math.floor(Number(nx)));
  const nyNum = Math.max(4, Math.floor(Number(ny)));

  const approx = useMemo(() => doubleIntegralTensor(f, axNum, bxNum, ayNum, byNum, nxNum, nyNum), [f, axNum, bxNum, ayNum, byNum, nxNum, nyNum]);

  // grid for heatmap/surface (coarse to keep performance reasonable)
  const gx = useMemo(() => linspace(axNum, bxNum, Math.min(80, nxNum + 1)), [axNum, bxNum, nxNum]);
  const gy = useMemo(() => linspace(ayNum, byNum, Math.min(80, nyNum + 1)), [ayNum, byNum, nyNum]);
  const gz = useMemo(() => gy.map((y) => gx.map((x) => f(x, y))), [gx, gy, f]);

  function handleExportCSV() {
    const rows = [["x", "y", "f(x,y)"]];
    for (let j = 0; j < gy.length; j++) for (let i = 0; i < gx.length; i++) rows.push([gx[i], gy[j], gz[j][i]]);
    downloadCSV("21_5_double_grid.csv", rows);
  }
  function handleExportPNG() {
    exportPlotPNG(plotId, "21_5_double_integral.png");
  }

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Calculator className="w-5 h-5" />
            21.5 Multiple Integrals — Double integral visualizer
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 md:grid-cols-6 gap-3">
            <Labeled label="f(x,y) (JS Math)">
              <Input value={expr} onChange={(e) => setExpr(e.target.value)} />
            </Labeled>
            <Labeled label="ax"><Input className="text-white" value={ax} onChange={(e) => setAx(e.target.value)} /></Labeled>
            <Labeled label="bx"><Input  className="text-white"  value={bx} onChange={(e) => setBx(e.target.value)} /></Labeled>
            <Labeled label="ay"><Input  className="text-white"  value={ay} onChange={(e) => setAy(e.target.value)} /></Labeled>
            <Labeled label="by"><Input  className="text-white"  value={by} onChange={(e) => setBy(e.target.value)} /></Labeled>
            <Labeled label="nx, ny">
              <div className="flex gap-2">
                <Input  value={nx} onChange={(e) => setNx(e.target.value)} className="w-28 text-white" />
                <Input value={ny} onChange={(e) => setNy(e.target.value)} className="w-28 text-white" />
              </div>
            </Labeled>
          </div>

        <Equation>{`Double integral via tensor-product trapezoid: $\\iint_R f(x,y) \\, dA \\approx \\sum_{i,j} w_{ij} f(x_i, y_j) \\, \\Delta x \\Delta y$.`}</Equation>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
            <ResponsivePlot
              id={plotId}
              data={[
                {
                  z: gz,
                  x: gx,
                  y: gy,
                  type: "heatmap",
                  colorscale: "Viridis",
                  name: "f(x,y)",
                  hovertemplate: "x=%{x}<br>y=%{y}<br>z=%{z}<extra></extra>",
                },
              ]}
              title={`Double integral ≈ ${fmt(approx, 8)}`}
            />

            <div className="rounded-xl border border-zinc-700 p-3 bg-zinc-900/40">
              <div className="text-zinc-300">Approximation (tensor trapezoid):</div>
              <div className="text-2xl mt-2 text-emerald-300 font-mono">{fmt(approx, 10)}</div>
              <div className="mt-3 flex gap-2">
                <Button 
className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black"
 onClick={handleExportCSV} variant="ghost"><Download className="w-4 h-4" /> CSV</Button>
                <Button 
className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black"
 onClick={handleExportPNG} variant="ghost"><ImageIcon className="w-4 h-4" /> PNG</Button>
              </div>

              <div className="mt-4 text-xs text-zinc-400">
                <Exercise title="Exercise 21.5.1" body="Compute ∬_[0,1]×[0,1] e^{-(x^2 + y^2)} dA with a tensor grid and compare with product of 1D integrals." solution="Refine nx,ny and compare to product of 1D integrals when f is separable; observe convergence." />
              </div>
            </div>
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// Documentation & Problems panels (Chapter 21)
// ======================================================================
const docs21 = {
  "21.1 Trapezoidal Rule": [
    "The trapezoidal rule approximates the integral of $f(x)$ over $[a,b]$ by linear interpolation between endpoints.",
    "$$\\int_a^b f(x) \\, dx \\approx \\frac{b-a}{2n} \\Big[f(a) + 2\\sum_{i=1}^{n-1} f(x_i) + f(b)\\Big].$$",
    "Composite trapezoid rule applies the basic trapezoid formula over $n$ subintervals to improve accuracy.",
    "Error is proportional to $O(h^2)$, where $h = (b-a)/n$ is the step size.",
    "Can handle uniform or adaptive step sizes depending on smoothness of $f(x)$."
  ],
  "21.2 Simpson's Rules": [
    "Simpson $1/3$ rule uses quadratic interpolants over pairs of subintervals:",
    "$$\\int_a^b f(x) \\, dx \\approx \\frac{h}{3} \\Big[f(a) + 4 f(a+h) + 2 f(a+2h) + \\dots + f(b)\\Big],$$ with $h=(b-a)/n$ and $n$ even.",
    "Simpson $3/8$ rule uses cubic polynomials over triples of subintervals:",
    "$$\\int_a^b f(x) \\, dx \\approx \\frac{3h}{8} \\Big[f(a) + 3 f(a+h) + 3 f(a+2h) + 2 f(a+3h) + \\dots + f(b)\\Big],$$ with $n$ divisible by 3.",
    "Simpson's rules are more accurate than trapezoid for smooth functions; error $O(h^4)$."
  ],
  "21.3 Unequal Segments": [
    "For nonuniform partitions, the trapezoidal rule sums contributions from each segment:",
    "$$\\int_a^b f(x) \\, dx \\approx \\sum_{i=0}^{n-1} \\frac{x_{i+1}-x_i}{2} \\Big[f(x_i) + f(x_{i+1})\\Big].$$",
    "Useful when data points are irregularly spaced or measurements are uneven.",
    "Error analysis depends on local segment sizes; smaller intervals improve accuracy."
  ],
  "21.4 Open Newton–Cotes": [
    "Open Newton–Cotes rules avoid using endpoints, e.g., midpoint rule or open 2-point rule.",
    "$$\\int_a^b f(x) \\, dx \\approx (b-a) f\\left(\\frac{a+b}{2}\\right) \\quad \\text{(Midpoint)}$$",
    "Useful when $f(x)$ is undefined or singular at endpoints.",
    "Open rules can also be extended to higher-order polynomials using interior nodes."
  ],
  "21.5 Multiple Integrals": [
    "For 2D rectangular domains, tensor-product trapezoid or Simpson rules apply 1D formulas in $x$ and $y$:",
    "$$\\iint_R f(x,y) \\, dA \\approx \\sum_{i,j} w_{ij} f(x_i, y_j) \\, \\Delta x \\Delta y.$$",
    "For complex domains, triangulation or Monte Carlo integration is preferred.",
    "Monte Carlo estimates are especially useful in higher dimensions where grid-based methods are expensive.",
    "Error analysis in multiple integrals depends on smoothness and dimensionality of $f(x,y)$."
  ],
};

function ChapterDocs21() {
  const [open, setOpen] = useState("21.1 Trapezoidal Rule");
  return (
    <div className="mt-6 bg-zinc-900 border border-zinc-700 rounded-xl p-4 overflow-hidden">
      <div className="flex items-center justify-between mb-3 gap-2 flex-wrap">
        <div>
          <div className="text-zinc-100 font-semibold">Chapter 21 — Notes</div>
          <div className="text-zinc-400 text-xs">Formulas, tips & reminders</div>
        </div>
        <div className="text-zinc-400 text-xs">Numerical integration</div>
      </div>

      <div className="grid grid-cols-1 md:grid-cols-4">
        <div className="border-b md:border-b-0 md:border-r border-zinc-800 pr-0 md:pr-2">
          {Object.keys(docs21).map((k) => (
            <button key={k} onClick={() => setOpen(k)} className={`w-full text-left p-2 hover:bg-zinc-800/30 ${open === k ? "bg-zinc-800/20" : ""}`}>
              <div className="flex items-center gap-2">
                <BookOpen className="w-4 h-4 text-zinc-300" />
                <div className="text-zinc-100">{k}</div>
              </div>
            </button>
          ))}
        </div>

        <div className="md:col-span-3 p-3 overflow-x-auto">
          <h3 className="text-zinc-100 font-semibold mb-2">{open}</h3>
          <div className="text-zinc-300 text-sm space-y-2">
            {docs21[open].map((p, i) => (
              <ReactMarkdown key={i} remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>{p}</ReactMarkdown>
            ))}
          </div>
        </div>
      </div>
    </div>
  );
}

// Problems 21
function Problems21() {
  const problems = [
    { title: "P21.1", body: "Use trapezoidal rule to approximate ∫_0^2 (1/(1+x^2)) dx for n=4,8,16 and study convergence." },
    { title: "P21.2", body: "Compare Simpson 1/3 and composite trapezoid on ∫_0^1 sin(10x) dx for various n." },
    { title: "P21.3", body: "Design a nonuniform grid that reduces error for integrating f(x)=e^{-x^2} on [0,3]." },
    { title: "P21.4", body: "Use open midpoint to approximate ∫_0^1 x^{-1/2} dx (improper at 0). How does midpoint help?" },
    { title: "P21.5", body: "Approximate ∬_D (x^2 + y^2) dA where D is a triangle using tensor grid and compare with exact integral." },
  ];
  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2"><Grid className="w-5 h-5" /> Problems</CardTitle>
        </CardHeader>
        <CardContent className="space-y-3">
          {problems.map((p) => (
            <div key={p.title} className="rounded-lg border border-zinc-700 p-3 bg-zinc-900/40">
              <div className="text-sm text-zinc-100 font-medium">{p.title}</div>
              <div className="text-xs text-zinc-300 mt-1">{p.body}</div>
            </div>
          ))}
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// Page assembly (final export)
// ======================================================================

export default function Chapter21() {
  // optional: small mount animation or accessibility focus management
  const topRef = useRef(null);
  useEffect(() => {
    if (topRef.current) topRef.current.focus?.();
  }, []);

  return (
    <div className={`p-4 md:p-6 lg:p-8 ${theme.bg} min-h-screen`}>
      <div className="max-w-screen-xl mx-auto space-y-6">
        <motion.div initial={{ opacity: 0, y: -6 }} animate={{ opacity: 1, y: 0 }}>
          <h1 ref={topRef} tabIndex={-1} className="text-2xl md:text-3xl lg:text-4xl font-bold" style={{ color: theme.accent2 }}>Numerical Integration</h1>
          <div className="text-zinc-400">Trapezoidal rule, Simpson's rules, unequal segments, open formulas, and multiple integrals — interactive visualizers and responsive plots.</div>
        </motion.div>

        <Section211 />
        <Section212 />
        <Section213 />
        <Section214 />
        <Section215 />

        <ChapterDocs21 />
        <Problems21 />
        <BottomBar/>
      </div>
    </div>
  );
}
