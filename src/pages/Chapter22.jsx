// src/pages/Chapter22.jsx
// ======================================================================
// Chapter 22 — Integration of Equations (responsive version)
// - All grids are mobile-first with safe wrapping
// - Plots are responsive, use clamp() height, no overflow
// - Tables and long content scroll horizontally on small screens
// - Kept your original functionality, exports, and styles
// ======================================================================

import React, { useMemo, useState } from "react";
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
  Grid,
  Layers,
  Activity,
  Download,
  ImageIcon,
  Zap,
  Percent,
  BookOpen,
  ChevronDown,
  ChevronRight,
} from "lucide-react";
import BottomBar from "../components/BottomBar";

// ---------- Theme & tiny helpers ----------
const theme = {
  bg: "bg-zinc-950 text-zinc-200",
  panel: "bg-zinc-900/60 border border-zinc-700",
  accent: "#22d3ee",
  accent2: "#34d399",
};

const fadeUp = {
  initial: { opacity: 0, y: 8 },
  animate: { opacity: 1, y: 0 },
  transition: { duration: 0.28 },
};

function fmt(v, p = 6) {
  if (!Number.isFinite(v)) return "NaN";
  return Number(v).toFixed(p).replace(/(?:\.0+|(\.\d+?)0+)$/, "$1");
}

function Labeled({ label, children }) {
  return (
    <div className="min-w-0">
      <div className="text-xs text-zinc-400 mb-1">{label}</div>
      {children}
    </div>
  );
}

function Equation({ children, color = theme.accent }) {
  return (
    <div className="my-2">
      <div className="px-3 py-2 rounded-lg text-gray-300  border border-zinc-700 bg-zinc-900/80 overflow-x-auto">
        <div className="text-xs uppercase tracking-widest mb-1" style={{ color }}>
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

function Note({ children }) {
  return <div className="text-sm text-zinc-300">{children}</div>;
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

// ---------- Utilities: safe function maker, CSV/PNG export ----------
function makeFunctionSafe1D(expr) {
  try {
    const clean = String(expr).replace(/\^/g, "**");
    // eslint-disable-next-line no-new-func
    const fn = new Function("x", `with (Math) { return (${clean}); }`);
    return (x) => {
      try {
        return Number(fn(x));
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
        return Number(fn(x, y));
      } catch {
        return NaN;
      }
    };
  } catch {
    return () => NaN;
  }
}

function downloadCSV(filename, rows) {
  const csv = rows
    .map((r) => r.map((c) => `"${String(c).replace(/"/g, '""')}"`).join(","))
    .join("\n");
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
      alert("Export failed: Plotly not available or element not found.");
    }
  } catch (err) {
    console.error(err);
    alert("Export failed: " + String(err));
  }
}

// ---------- Shared Plot wrapper (responsive, no overflow) ----------
function ResponsivePlot({ plotId, data, title, extraLayout = {} }) {
  return (
    <div
      id={plotId}
      className="rounded-xl border border-zinc-700 p-2 bg-zinc-900/40 overflow-hidden w-full"
      style={{ maxWidth: "100%" }}
    >
      <Plot
        data={data}
        layout={{
          title,
          paper_bgcolor: "rgba(0,0,0,0)",
          plot_bgcolor: "rgba(0,0,0,0)",
          font: { color: "#e5e7eb" },
          autosize: true,
          margin: { l: 50, r: 20, t: 50, b: 45 },
          legend: { orientation: "h", x: 0, y: -0.15 },
          ...extraLayout,
        }}
        config={{ responsive: true, displayModeBar: false }}
        useResizeHandler
        style={{ width: "100%", height: "clamp(260px, 52vh, 520px)" }}
      />
    </div>
  );
}

// ---------- Numerical routines ----------

// Composite trapezoid (m subintervals)
function compositeTrapezoid(f, a, b, m) {
  if (m < 1) m = 1;
  const h = (b - a) / m;
  let s = 0.5 * (f(a) + f(b));
  for (let i = 1; i < m; i++) s += f(a + i * h);
  return s * h;
}

// Composite Simpson 1/3 (m must be even)
function compositeSimpson13(f, a, b, m) {
  if (m % 2 === 1) m += 1;
  const h = (b - a) / m;
  let s = f(a) + f(b);
  for (let i = 1; i < m; i++) s += f(a + i * h) * (i % 2 === 0 ? 2 : 4);
  return (h / 3) * s;
}

// Composite Simpson 3/8 (m must be multiple of 3)
function compositeSimpson38(f, a, b, m) {
  while (m % 3 !== 0) m += 1;
  const h = (b - a) / m;
  let s = f(a) + f(b);
  for (let i = 1; i < m; i++) s += f(a + i * h) * (i % 3 === 0 ? 2 : 3);
  return (3 * h) / 8 * s;
}

// Newton-Cotes closed general (weights for small n can be hard-coded)
function newtonCotesClosed(f, a, b, n) {
  // n = degree: 1 -> trapezoid, 2 -> Simpson 1/3, 3 -> 3/8
  if (n === 1) return compositeTrapezoid(f, a, b, 1);
  if (n === 2) return compositeSimpson13(f, a, b, 2);
  if (n === 3) return compositeSimpson38(f, a, b, 3);
  // For higher degree, fallback to composite trapezoid with many panels
  return compositeTrapezoid(f, a, b, 200);
}

// Romberg integration table
function rombergIntegration(f, a, b, maxK = 6) {
  const R = Array.from({ length: maxK }, () => []);
  for (let k = 0; k < maxK; k++) {
    const n = 2 ** k;
    const T = compositeTrapezoid(f, a, b, n);
    R[k][0] = T;
    for (let j = 1; j <= k; j++) {
      R[k][j] =
        R[k][j - 1] + (R[k][j - 1] - R[k - 1][j - 1]) / (4 ** j - 1);
    }
  }
  return R;
}

// Adaptive Simpson's method (recursive)
function adaptiveSimpson(f, a, b, tol = 1e-8, maxDepth = 20) {
  function simpson(f, a, b) {
    const c = 0.5 * (a + b);
    return ((b - a) / 6) * (f(a) + 4 * f(c) + f(b));
  }
  function recurse(f, a, b, eps, S, depth) {
    const c = 0.5 * (a + b);
    const Sleft = simpson(f, a, c);
    const Sright = simpson(f, c, b);
    const err = Sleft + Sright - S;
    if (depth <= 0 || Math.abs(err) < 15 * eps) {
      return Sleft + Sright + err / 15;
    }
    return (
      recurse(f, a, c, eps / 2, Sleft, depth - 1) +
      recurse(f, c, b, eps / 2, Sright, depth - 1)
    );
  }
  const S = simpson(f, a, b);
  return recurse(f, a, b, tol, S, maxDepth);
}

// Legendre polynomials and root finding (for general Gauss n)
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
  const nodes = new Array(n);
  const weights = new Array(n);
  const m = Math.floor((n + 1) / 2);
  for (let i = 0; i < m; i++) {
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
  if (n <= 0) return 0;
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

// Improper integrals handling: infinite limits and endpoint singularities
function improperIntegralInfinite(f, a, b, tol = 1e-8) {
  if (!isFinite(a) && !isFinite(b)) {
    const transform = (t) => {
      const x = Math.tan(Math.PI * (t - 0.5));
      const dxdt = Math.PI / Math.cos(Math.PI * (t - 0.5)) ** 2;
      return f(x) * dxdt;
    };
    return adaptiveSimpson(transform, 0, 1, tol);
  }
  if (!isFinite(b)) {
    const transform = (t) => {
      const x = a + t / (1 - t);
      const dxdt = 1 / (1 - t) ** 2;
      return f(x) * dxdt;
    };
    return adaptiveSimpson(transform, 0, 1 - 1e-12, tol);
  }
  if (!isFinite(a)) {
    const transform = (t) => {
      const x = b - t / (1 - t);
      const dxdt = 1 / (1 - t) ** 2;
      return f(x) * dxdt;
    };
    return adaptiveSimpson(transform, 0, 1 - 1e-12, tol);
  }
  return adaptiveSimpson(f, a, b, tol);
}

// ---------- UI Sections ----------

// Section 22.1: Newton-Cotes algorithms visualizer
function Section221() {
  const [expr, setExpr] = useState("Math.sin(x) / x");
  const [a, setA] = useState("0.1");
  const [b, setB] = useState("10");
  const [panels, setPanels] = useState(10);
  const [degree, setDegree] = useState(2);
  const plotId = "plot-22-1";

  const f = useMemo(() => makeFunctionSafe1D(expr), [expr]);
  const aNum = Number(a);
  const bNum = Number(b);
  const m = Math.max(1, Math.floor(Number(panels)));
  const deg = Math.max(1, Math.min(3, Math.floor(Number(degree))));

  const approx = useMemo(() => {
    const h = (bNum - aNum) / m;
    let total = 0;
    for (let i = 0; i < m; i++) {
      const xa = aNum + i * h;
      const xb = xa + h;
      if (deg === 1) total += compositeTrapezoid(f, xa, xb, 1);
      else if (deg === 2) total += compositeSimpson13(f, xa, xb, 2);
      else if (deg === 3) total += compositeSimpson38(f, xa, xb, 3);
    }
    return total;
  }, [f, aNum, bNum, m, deg]);

  const xs = useMemo(() => {
    const L = 400;
    if (!isFinite(aNum) || !isFinite(bNum) || bNum === aNum) return [];
    return new Array(L)
      .fill(0)
      .map((_, i) => aNum + (i / (L - 1)) * (bNum - aNum));
  }, [aNum, bNum]);

  const ys = useMemo(() => xs.map((xx) => f(xx)), [xs, f]);

  function handleExportCSV() {
    const rows = [["x", "f(x)"]];
    for (let i = 0; i < xs.length; i++) rows.push([xs[i], ys[i]]);
    downloadCSV("22_1_newton_cotes.csv", rows);
  }
  function handleExportPNG() {
    exportPlotPNG(plotId, "22_1_newton_cotes.png");
  }

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <FunctionSquare className="w-5 h-5" />
            22.1 Newton–Cotes Algorithms for Equations
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-4 gap-3">
            <Labeled label="f(x) (JS Math namespace)">
              <Textarea
                value={expr}
                onChange={(e) => setExpr(e.target.value)}
                className="h-24 text-white"
              />
            </Labeled>
            <Labeled label="a">
              <Input className="text-white" value={a} onChange={(e) => setA(e.target.value)} />
            </Labeled>
            <Labeled label="b">
              <Input className="text-white"  value={b} onChange={(e) => setB(e.target.value)} />
            </Labeled>
            <Labeled label="Panels m">
              <Input
                value={panels}
                onChange={(e) => setPanels(e.target.value)}
                className="text-white" 
              />
            </Labeled>
          </div>

          <div className="flex items-center gap-2 flex-wrap">
            <div className="text-xs text-zinc-400">Choose NC degree:</div>
            <button
              className={`px-3 py-1 rounded ${
                degree === 1
                  ? "bg-emerald-600 cursor-pointer text-white"
                  : "bg-zinc-800 cursor-pointer text-zinc-200"
              }`}
              onClick={() => setDegree(1)}
            >
              Trapezoid
            </button>
            <button
              className={`px-3 py-1 rounded ${
                degree === 2
                  ? "bg-emerald-600 cursor-pointer text-white"
                  : "bg-zinc-800 cursor-pointer text-zinc-200"
              }`}
              onClick={() => setDegree(2)}
            >
              Simpson 1/3
            </button>
            <button
              className={`px-3 py-1 rounded ${
                degree === 3
                  ? "bg-emerald-600 cursor-pointer text-white"
                  : "bg-zinc-800 cursor-pointer text-zinc-200"
              }`}
              onClick={() => setDegree(3)}
            >
              Simpson 3/8
            </button>
          </div>

          <Equation>
  {`Composite Newton–Cotes: $\\int_a^b f(x) \\, dx \\approx \\sum_{i=0}^{n} w_i f(x_i)$. Closed Newton–Cotes use endpoint values.`}
</Equation>


          <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
            <ResponsivePlot
              plotId={plotId}
              data={[
                {
                  x: xs,
                  y: ys,
                  mode: "lines",
                  type: "scatter",
                  name: "f(x)",
                  line: { color: theme.accent },
                },
              ]}
              title={`Newton–Cotes composite (approx ≈ ${fmt(approx, 8)})`}
            />

            <div className="rounded-xl border border-zinc-700 p-3 bg-zinc-900/40 overflow-hidden">
              <div className="text-zinc-300">Approximation:</div>
              <div className="text-2xl text-emerald-300 mt-2 font-mono">
                {fmt(approx, 10)}
              </div>
              <div className="mt-3 flex gap-2 flex-wrap">
                <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={handleExportCSV} variant="ghost">
                  <Download className="w-4 h-4" /> CSV
                </Button>
                <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={handleExportPNG} variant="ghost">
                  <ImageIcon className="w-4 h-4" /> PNG
                </Button>
              </div>

              <div className="mt-4 text-xs text-zinc-400">
                <Exercise
                  title="Exercise 22.1.1"
                  body="Compare composite trapezoid and Simpson 1/3 on ∫_0^2 e^{-x^2} dx for m=16 and m=32."
                  solution="Simpson typically converges faster for smooth functions; compute numerical errors to observe rates."
                />
              </div>
            </div>
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// Section 22.2: Romberg Integration visualizer & table
function Section222() {
  const [expr, setExpr] = useState("Math.sin(x)/x");
  const [a, setA] = useState("0.1");
  const [b, setB] = useState("10");
  const [levels, setLevels] = useState(6);
  const plotId = "plot-22-2";

  const f = useMemo(() => makeFunctionSafe1D(expr), [expr]);
  const aNum = Number(a);
  const bNum = Number(b);
  const lev = Math.max(1, Math.min(10, Math.floor(Number(levels))));

  const R = useMemo(
    () => rombergIntegration(f, aNum, bNum, lev),
    [f, aNum, bNum, lev]
  );

  const best = R && R.length ? R[R.length - 1][R.length - 1] : NaN;

  const xs = useMemo(() => {
    const L = 400;
    if (!isFinite(aNum) || !isFinite(bNum) || bNum === aNum) return [];
    return new Array(L)
      .fill(0)
      .map((_, i) => aNum + (i / (L - 1)) * (bNum - aNum));
  }, [aNum, bNum]);

  const ys = useMemo(() => xs.map((xx) => f(xx)), [xs, f]);

  function handleExportCSV() {
    const rows = [["level", ...Array.from({ length: R.length }, (_, j) => `R_${j}`)]];
    for (let i = 0; i < R.length; i++) {
      rows.push([i + 1, ...R[i].map((v) => fmt(v, 12))]);
    }
    downloadCSV("22_2_romberg_table.csv", rows);
  }
  function handleExportPNG() {
    exportPlotPNG(plotId, "22_2_romberg.png");
  }

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Layers className="w-5 h-5" />
            22.2 Romberg Integration
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-4 gap-3">
            <Labeled label="f(x)">
              <Textarea
                value={expr}
                onChange={(e) => setExpr(e.target.value)}
                className="h-24 text-white"
              />
            </Labeled>
            <Labeled label="a">
              <Input className="text-white"  value={a} onChange={(e) => setA(e.target.value)} />
            </Labeled>
            <Labeled label="b">
              <Input className="text-white"  value={b} onChange={(e) => setB(e.target.value)} />
            </Labeled>
            <Labeled label="levels (max 10)">
              <Input
                value={levels}
                onChange={(e) => setLevels(e.target.value)}
                className="text-white" 
              />
            </Labeled>
          </div>

         <Equation>
  {`Romberg extrapolation: $R_{k,j} = R_{k,j-1} + \\frac{R_{k,j-1} - R_{k-1,j-1}}{4^{j}-1}$`}
</Equation>


          <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
            <ResponsivePlot
              plotId={plotId}
              data={[
                {
                  x: xs,
                  y: ys,
                  mode: "lines",
                  type: "scatter",
                  name: "f(x)",
                  line: { color: theme.accent },
                },
              ]}
              title={`Romberg (best ≈ ${fmt(best, 10)})`}
            />

            <div className="rounded-xl border border-zinc-700 p-3 bg-zinc-900/40 overflow-x-auto">
              <div className="text-zinc-300 mb-2">
                Romberg table (rows: k, columns: j)
              </div>
              <div className="text-xs text-zinc-100 font-mono min-w-fit">
                {R.map((row, i) => (
                  <div key={i} className="whitespace-nowrap">
                    {row.map((v, j) => (
                      <span
                        key={j}
                        style={{ display: "inline-block", minWidth: 220 }}
                      >{`R[${i},${j}] = ${fmt(v, 12)}`}</span>
                    ))}
                  </div>
                ))}
              </div>

              <div className="mt-3 flex gap-2 flex-wrap">
                <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={handleExportCSV} variant="ghost">
                  <Download className="w-4 h-4" /> CSV
                </Button>
                <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={handleExportPNG} variant="ghost">
                  <ImageIcon className="w-4 h-4" /> PNG
                </Button>
              </div>

              <div className="mt-4 text-xs text-zinc-400">
                <Exercise
                  title="Exercise 22.2.1"
                  body="Build Romberg table for ∫_0^1 e^{-x^2} dx and observe convergence as level increases."
                  solution="Romberg accelerates trap rule with Richardson extrapolation; expect rapid convergence for smooth f."
                />
              </div>
            </div>
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// Section 22.3: Adaptive Quadrature visualizer (Adaptive Simpson)
function Section223() {
  const [expr, setExpr] = useState("1/(1 + x*x)");
  const [a, setA] = useState("0");
  const [b, setB] = useState("10");
  const [tol, setTol] = useState("1e-6");
  const [maxDepth, setMaxDepth] = useState("20");
  const plotId = "plot-22-3";

  const f = useMemo(() => makeFunctionSafe1D(expr), [expr]);
  const aNum = Number(a);
  const bNum = Number(b);
  const tolNum = Number(tol);
  const maxD = Math.max(1, Math.floor(Number(maxDepth)));

  const approx = useMemo(
    () => adaptiveSimpson(f, aNum, bNum, tolNum, maxD),
    [f, aNum, bNum, tolNum, maxD]
  );

  const xs = useMemo(() => {
    const L = 400;
    if (!isFinite(aNum) || !isFinite(bNum) || bNum === aNum) return [];
    return new Array(L)
      .fill(0)
      .map((_, i) => aNum + (i / (L - 1)) * (bNum - aNum));
  }, [aNum, bNum]);

  const ys = useMemo(() => xs.map((xx) => f(xx)), [xs, f]);

  function handleExportCSV() {
    const rows = [["x", "f(x)"]];
    for (let i = 0; i < xs.length; i++) rows.push([xs[i], ys[i]]);
    downloadCSV("22_3_adaptive_simpson.csv", rows);
  }
  function handleExportPNG() {
    exportPlotPNG(plotId, "22_3_adaptive_simpson.png");
  }

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Activity className="w-5 h-5" />
            22.3 Adaptive Quadrature (Adaptive Simpson)
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-4 gap-3">
            <Labeled label="f(x)">
              <Textarea
                value={expr}
                onChange={(e) => setExpr(e.target.value)}
                className="h-24 text-white"
              />
            </Labeled>
            <Labeled label="a">
              <Input className="text-white"  value={a} onChange={(e) => setA(e.target.value)} />
            </Labeled>
            <Labeled label="b">
              <Input className="text-white"  value={b} onChange={(e) => setB(e.target.value)} />
            </Labeled>
            <Labeled label="tol">
              <Input className="text-white"  value={tol} onChange={(e) => setTol(e.target.value)} />
            </Labeled>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
            <ResponsivePlot
              plotId={plotId}
              data={[
                {
                  x: xs,
                  y: ys,
                  mode: "lines",
                  type: "scatter",
                  name: "f(x)",
                  line: { color: theme.accent },
                },
              ]}
              title={`Adaptive Simpson approx ≈ ${fmt(approx, 10)}`}
            />

            <div className="rounded-xl border border-zinc-700 p-3 bg-zinc-900/40 overflow-hidden">
              <div className="text-zinc-300">Approximation:</div>
              <div className="text-2xl mt-2 text-emerald-300 font-mono">
                {fmt(approx, 10)}
              </div>
              <div className="mt-3 flex gap-2 flex-wrap">
                <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={handleExportCSV} variant="ghost">
                  <Download className="w-4 h-4" /> CSV
                </Button>
                <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={handleExportPNG} variant="ghost">
                  <ImageIcon className="w-4 h-4" /> PNG
                </Button>
              </div>

              <div className="mt-4 text-xs text-zinc-400">
                <Exercise
                  title="Exercise 22.3.1"
                  body="Use adaptive Simpson to compute ∫_0^100 e^{-x} dx with varying tolerances. Observe required recursion depth."
                  solution="For rapidly decaying integrands adaptive Simpson concentrates nodes near origin; error tolerance controls subdivisions."
                />
              </div>
            </div>
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// Section 22.4: Gauss Quadrature (general n)
function Section224() {
  const [expr, setExpr] = useState("Math.exp(-x*x)");
  const [a, setA] = useState("-1");
  const [b, setB] = useState("1");
  const [n, setN] = useState(4);
  const plotId = "plot-22-4";

  const f = useMemo(() => makeFunctionSafe1D(expr), [expr]);
  const aNum = Number(a);
  const bNum = Number(b);
  const nNum = Math.max(1, Math.min(50, Math.floor(Number(n))));

  const approx = useMemo(
    () => gaussLegendreQuad(f, aNum, bNum, nNum),
    [f, aNum, bNum, nNum]
  );

  const gw = useMemo(() => gaussNodesWeights(nNum), [nNum]);

  const xs = useMemo(() => {
    const L = 300;
    if (!isFinite(aNum) || !isFinite(bNum) || bNum === aNum) return [];
    return new Array(L)
      .fill(0)
      .map((_, i) => aNum + (i / (L - 1)) * (bNum - aNum));
  }, [aNum, bNum]);

  const ys = useMemo(() => xs.map((xx) => f(xx)), [xs, f]);

  function handleExportCSV() {
    const rows = [["node(-1,1)", "weight(-1,1)"]];
    for (let i = 0; i < gw.nodes.length; i++)
      rows.push([gw.nodes[i], gw.weights[i]]);
    downloadCSV(`22_4_gauss_n${nNum}.csv`, rows);
  }
  function handleExportPNG() {
    exportPlotPNG(plotId, `22_4_gauss_n${nNum}.png`);
  }

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Zap className="w-5 h-5" />
            22.4 Gauss Quadrature (Legendre-based, general n)
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-4 gap-3">
            <Labeled label="f(x)">
              <Textarea
                value={expr}
                onChange={(e) => setExpr(e.target.value)}
                className="h-24 text-white"
              />
            </Labeled>
            <Labeled label="a">
              <Input className="text-white"  value={a} onChange={(e) => setA(e.target.value)} />
            </Labeled>
            <Labeled   label="b">
              <Input className="text-white"  value={b} onChange={(e) => setB(e.target.value)} />
            </Labeled>
            <Labeled label="n (points)">
              <Input className="text-white"  value={n} onChange={(e) => setN(e.target.value)} />
            </Labeled>
          </div>

          <Equation>
  {`Gauss–Legendre: $\\int_a^b f(x) \\, dx \\approx \\frac{b-a}{2} \\sum_{i=1}^n w_i f\\Big(\\frac{b-a}{2} x_i + \\frac{a+b}{2}\\Big), \\quad x_i, w_i \\text{ from roots of } P_n.$`}
</Equation>


          <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
            <ResponsivePlot
              plotId={plotId}
              data={[
                {
                  x: xs,
                  y: ys,
                  mode: "lines",
                  type: "scatter",
                  name: "f(x)",
                  line: { color: theme.accent },
                },
                {
                  x: gw.nodes.map(
                    (xi) => ((bNum - aNum) / 2) * xi + (aNum + bNum) / 2
                  ),
                  y: gw.nodes.map((xi) =>
                    f(((bNum - aNum) / 2) * xi + (aNum + bNum) / 2)
                  ),
                  mode: "markers",
                  marker: { size: 8, color: theme.accent2 },
                  name: "Gauss nodes",
                },
              ]}
              title={`Gauss n=${nNum} (approx ≈ ${fmt(approx, 10)})`}
            />

            <div className="rounded-xl border border-zinc-700 p-3 bg-zinc-900/40 overflow-x-auto">
              <div className="text-zinc-300 mb-2">
                Nodes & weights (on [-1,1])
              </div>
              <div className="text-xs text-zinc-100 font-mono min-w-fit">
                {gw.nodes.map((nd, i) => (
                  <div key={i} className="whitespace-nowrap">
                    x_{i + 1}={fmt(nd, 10)} &nbsp;&nbsp; w_{i + 1}=
                    {fmt(gw.weights[i], 12)}
                  </div>
                ))}
              </div>

              <div className="mt-3 flex gap-2 flex-wrap">
                <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={handleExportCSV} variant="ghost">
                  <Download className="w-4 h-4" /> CSV
                </Button>
                <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={handleExportPNG} variant="ghost">
                  <ImageIcon className="w-4 h-4" /> PNG
                </Button>
              </div>

              <div className="mt-4 text-xs text-zinc-400">
                <Exercise
                  title="Exercise 22.4.1"
                  body="Compare Gauss n=4 vs composite Simpson 1/3 with many panels on ∫_{-1}^{1} (1/(1+x^2)) dx."
                  solution="Gauss of moderate n often outperforms composite rules for smooth integrands since it has optimal degree of exactness."
                />
              </div>
            </div>
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// Section 22.5: Improper Integrals (infinite limits & endpoint singularities)
function Section225() {
  const [expr, setExpr] = useState("Math.exp(-x)");
  const [a, setA] = useState("0");
  const [b, setB] = useState("Infinity");
  const [tol, setTol] = useState("1e-8");
  const plotId = "plot-22-5";

  const f = useMemo(() => {
    try {
      const clean = String(expr).replace(/\^/g, "**");
      // eslint-disable-next-line no-new-func
      const fn = new Function("x", `with (Math) { return (${clean}); }`);
      return (x) => {
        try {
          return Number(fn(x));
        } catch {
          return NaN;
        }
      };
    } catch {
      return () => NaN;
    }
  }, [expr]);

  const parseLimit = (s) => {
    if (typeof s !== "string") return Number(s);
    const t = s.trim().toLowerCase();
    if (t === "infinity" || t === "inf" || t === "+infinity" || t === "+inf")
      return Infinity;
    if (t === "-infinity" || t === "-inf") return -Infinity;
    const v = Number(s);
    return Number.isFinite(v) ? v : NaN;
  };

  const aNum = useMemo(() => parseLimit(String(a)), [a]);
  const bNum = useMemo(() => parseLimit(String(b)), [b]);
  const tolNum = Number(tol);

  const approx = useMemo(() => {
    try {
      if (!isFinite(aNum) || !isFinite(bNum)) {
        return improperIntegralInfinite(f, aNum, bNum, tolNum);
      }
      return adaptiveSimpson(f, aNum, bNum, tolNum, 20);
    } catch {
      return NaN;
    }
  }, [f, aNum, bNum, tolNum]);

  const xs = useMemo(() => {
    let A = aNum,
      B = bNum;
    if (!isFinite(A) && isFinite(B)) A = B - 10;
    if (!isFinite(B) && isFinite(A)) B = A + 10;
    if (!isFinite(A) && !isFinite(B)) {
      A = -10;
      B = 10;
    }
    const L = 400;
    return new Array(L).fill(0).map((_, i) => A + (i / (L - 1)) * (B - A));
  }, [aNum, bNum]);

  const ys = useMemo(() => xs.map((xx) => f(xx)), [xs, f]);

  function handleExportCSV() {
    const rows = [["x", "f(x)"]];
    for (let i = 0; i < xs.length; i++) rows.push([xs[i], ys[i]]);
    downloadCSV("22_5_improper.csv", rows);
  }
  function handleExportPNG() {
    exportPlotPNG(plotId, "22_5_improper.png");
  }

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Percent className="w-5 h-5" />
            22.5 Improper Integrals — Infinite limits & singular endpoints
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-4 gap-3">
            <Labeled label="f(x)">
              <Textarea
                value={expr}
                onChange={(e) => setExpr(e.target.value)}
                className="h-24 text-white"
              />
            </Labeled>
            <Labeled label="a (use 'Infinity' or '-Infinity' for infinite)">
              <Input className="text-white"  value={a} onChange={(e) => setA(e.target.value)} />
            </Labeled>
            <Labeled label="b (use 'Infinity' or '-Infinity')">
              <Input className="text-white"  value={b} onChange={(e) => setB(e.target.value)} />
            </Labeled>
            <Labeled label="tol">
              <Input className="text-white"  value={tol} onChange={(e) => setTol(e.target.value)} />
            </Labeled>
          </div>

          <Equation>
  {`Transformations map infinite intervals to finite ones (e.g., $x = \\tan(\\pi(t-1/2))$ for $(-\\infty,\\infty)$) or use rational transforms for semi-infinite intervals.`}
</Equation>


          <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
            <ResponsivePlot
              plotId={plotId}
              data={[
                {
                  x: xs,
                  y: ys,
                  mode: "lines",
                  type: "scatter",
                  name: "f(x)",
                  line: { color: theme.accent },
                },
              ]}
              title={`Improper integral approx ≈ ${fmt(approx, 10)}`}
            />

            <div className="rounded-xl border border-zinc-700 p-3 bg-zinc-900/40 overflow-hidden">
              <div className="text-zinc-300">Approximation (numerical)</div>
              <div className="text-2xl mt-2 text-emerald-300 font-mono">
                {fmt(approx, 10)}
              </div>
              <div className="mt-3 flex gap-2 flex-wrap">
                <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={handleExportCSV} variant="ghost">
                  <Download className="w-4 h-4" /> CSV
                </Button>
                <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={handleExportPNG} variant="ghost">
                  <ImageIcon className="w-4 h-4" /> PNG
                </Button>
              </div>

              <div className="mt-4 text-xs text-zinc-400">
                <Exercise
                  title="Exercise 22.5.1"
                  body="Compute ∫_0^∞ e^{-x} dx by transformation and compare to exact value 1."
                  solution="Transform x=t/(1-t) with t∈[0,1), apply adaptive quadrature; result ≈ 1."
                />
              </div>
            </div>
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ---------- Documentation & Problems panels ----------
const docs22 = {
  "22.1 Newton–Cotes": [
    "Newton–Cotes formulas approximate integrals using weighted sums of function values at equally spaced nodes:",
    "$$\\int_a^b f(x) \\, dx \\approx \\sum_{i=0}^n w_i f(x_i),$$ where $w_i$ are weights determined by interpolating polynomials.",
    "Closed Newton–Cotes rules include the endpoints $a$ and $b$ (e.g., trapezoidal rule, Simpson's 1/3 rule).",
    "Open Newton–Cotes rules use only interior nodes, which is useful when $f(x)$ is undefined or singular at endpoints.",
    "Higher-order Newton–Cotes formulas can improve accuracy but may suffer from oscillations (Runge phenomenon) for large $n$."
  ],
  "22.2 Romberg": [
    "Romberg integration refines trapezoidal rule estimates using Richardson extrapolation.",
    "$$R_{k,j} = R_{k,j-1} + \\frac{R_{k,j-1} - R_{k-1,j-1}}{4^{j-1} - 1},$$ where $R_{k,0}$ are trapezoidal estimates with step size $h_k$.",
    "By building a Romberg table, one can systematically improve integral estimates.",
    "Error decreases rapidly with $j$, often achieving very high accuracy with few function evaluations."
  ],
  "22.3 Adaptive Quadrature": [
    "Adaptive quadrature recursively subdivides intervals based on estimated local error.",
    "For example, adaptive Simpson's rule checks:",
    "$$|S(a,c) + S(c,b) - S(a,b)| > \\text{tolerance},$$ where $S(a,b)$ is the Simpson estimate for interval $[a,b]$ and $c=(a+b)/2$.",
    "Intervals are split further only where the function varies rapidly, making computation efficient.",
    "Adaptive quadrature is especially useful for functions with sharp peaks or near singularities."
  ],
  "22.4 Gauss Quadrature": [
    "Gauss quadrature achieves exact integration for polynomials of degree up to $2n-1$ using $n$ points and weights.",
    "$$\\int_{-1}^1 f(x) \\, dx \\approx \\sum_{i=1}^n w_i f(x_i),$$ where $x_i$ are the roots of the Legendre polynomial $P_n(x)$ and $w_i$ are corresponding weights.",
    "Change of variables allows application on arbitrary intervals $[a,b]$:",
    "$$\\int_a^b f(x) \\, dx = \\frac{b-a}{2} \\int_{-1}^1 f\\Big(\\frac{b-a}{2} t + \\frac{a+b}{2}\\Big) \\, dt.$$",
    "Highly efficient for smooth integrands, requiring fewer function evaluations than Newton–Cotes for the same accuracy."
  ],
  "22.5 Improper Integrals": [
    "Improper integrals arise from infinite intervals or singularities at endpoints.",
    "Infinite intervals can be transformed using substitutions, e.g., $x = 1/t$ maps $[1,\\infty)$ to $[0,1]$.",
    "Near singularities, split the integral and apply quadrature separately on each piece.",
    "Example: $$\\int_0^1 \\frac{f(x)}{\\sqrt{x}} \\, dx = \\int_0^\\epsilon \\frac{f(x)}{\\sqrt{x}} \\, dx + \\int_\\epsilon^1 \\frac{f(x)}{\\sqrt{x}} \\, dx,$$ with small $\\epsilon$.",
    "Adaptive quadrature often handles improper integrals efficiently by focusing on difficult regions."
  ],
};


function ChapterDocs22() {
  const [open, setOpen] = useState("22.1 Newton–Cotes");
  return (
    <div className="mt-6 bg-zinc-900 border border-zinc-700 rounded-xl p-4 overflow-hidden">
      <div className="flex items-center justify-between mb-3 gap-2 flex-wrap">
        <div>
          <div className="text-zinc-100 font-semibold">Chapter 22 — Notes</div>
          <div className="text-zinc-400 text-xs">
            Integration of equations: overview & formulas
          </div>
        </div>
        <div className="text-zinc-400 text-xs">Numerical integration</div>
      </div>

      <div className="grid grid-cols-1 md:grid-cols-4">
        <div className="border-b md:border-b-0 md:border-r border-zinc-800 pr-0 md:pr-2">
          {Object.keys(docs22).map((k) => (
            <button
              key={k}
              onClick={() => setOpen(k)}
              className={`w-full text-left p-2 hover:bg-zinc-800/30 ${
                open === k ? "bg-zinc-800/20" : ""
              }`}
            >
              <div className="flex items-center gap-2">
                <BookOpen className="w-4 h-4 text-zinc-300" />
                <div className="text-zinc-100">{k}</div>
              </div>
            </button>
          ))}
        </div>

        <div className="md:col-span-3 p-3 overflow-x-auto">
          <h3 className="text-zinc-100 font-semibold mb-2 min-w-fit">{open}</h3>
          <div className="text-zinc-300 text-sm space-y-2 min-w-fit">
            {docs22[open].map((p, i) => (
              <ReactMarkdown
                key={i}
                remarkPlugins={[remarkMath]}
                rehypePlugins={[rehypeKatex]}
              >
                {p}
              </ReactMarkdown>
            ))}
          </div>
        </div>
      </div>
    </div>
  );
}

// Problems
function Problems22() {
  const problems = [
    {
      title: "P22.1",
      body:
        "Use Romberg integration to compute ∫_0^1 ln(1+x^2) dx and report convergence table.",
    },
    {
      title: "P22.2",
      body:
        "Implement adaptive Simpson for ∫_0^100 x^2 e^{-x} dx and determine tolerance needed for 6-digit accuracy.",
    },
    {
      title: "P22.3",
      body:
        "Compare Gauss n=5 with composite Simpson on ∫_{-1}^1 1/(1+x^2) dx for error and efficiency.",
    },
    {
      title: "P22.4",
      body:
        "Compute ∫_0^∞ x^{1/2} e^{-x} dx (Gamma function) using appropriate transform and compare with exact value Γ(3/2)=√π/2.",
    },
  ];
  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Grid className="w-5 h-5" /> Problems
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-3">
          {problems.map((p) => (
            <div
              key={p.title}
              className="rounded-lg border border-zinc-700 p-3 bg-zinc-900/40"
            >
              <div className="text-sm text-zinc-100 font-medium">
                {p.title}
              </div>
              <div className="text-xs text-zinc-300 mt-1">{p.body}</div>
            </div>
          ))}
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ---------- Page assembly ----------
export default function Chapter22() {
  return (
    <div className={`p-4 md:p-6 lg:p-8 ${theme.bg} min-h-screen`}>
      <div className="max-w-screen-xl mx-auto space-y-6">
        <motion.div initial={{ opacity: 0, y: -10 }} animate={{ opacity: 1, y: 0 }}>
          <h1
            className="text-2xl md:text-3xl lg:text-4xl font-bold"
            style={{ color: theme.accent2 }}
          >
            Integration of Equations
          </h1>
          <div className="text-zinc-400">
            Newton–Cotes, Romberg, Adaptive Quadrature, Gauss Quadrature, and
            Improper Integrals — interactive visualizers.
          </div>
        </motion.div>

        <Section221 />
        <Section222 />
        <Section223 />
        <Section224 />
        <Section225 />

        <ChapterDocs22 />
        <Problems22 />
        <BottomBar/>
      </div>
    </div>
  );
}
