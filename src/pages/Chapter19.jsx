// src/pages/Chapter19.jsx
// ======================================================================
// Chapter 19 — Numerical ODEs & Initial Value Problems (Responsive)
// - Explicit & implicit Euler, Midpoint, RK4, adaptive RK, stability notes
// - Phase portraits, solution overlays, error plots, CSV/PNG export
// - Fully responsive Plotly plots (autosize + useResizeHandler + parent sizing)
// - Mobile-first layout, overflow protection, helpful exercises & documentation
// - Expanded content and comments to create a substantial, study-ready file
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
  Zap,
  Layers,
  Activity,
  Download,
  ImageIcon,
  Grid,
  BookOpen,
  Calculator,
  FunctionSquare,
  Percent,
  ChevronDown,
  ChevronRight,
} from "lucide-react";
import BottomBar from "../components/BottomBar";

// ======================================================================
// Theme & tiny helpers
// ======================================================================
const theme = {
  bg: "bg-zinc-950 text-zinc-200",
  panel: "bg-zinc-900/60 border border-zinc-700",
  accent: "#22d3ee",
  accent2: "#34d399",
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

function linspace(a, b, n) {
  n = Math.max(1, Math.floor(n));
  if (n === 1) return [a];
  const out = new Array(n);
  const dx = (b - a) / (n - 1);
  for (let i = 0; i < n; i++) out[i] = a + i * dx;
  return out;
}

// ======================================================================
// UI Microcomponents
// ======================================================================
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

function Exercise({ title, body, solution }) {
  const [open, setOpen] = useState(false);
  return (
    <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-3">
      <div className="flex items-center justify-between">
        <div className="text-zinc-100 font-medium">{title}</div>
        <Button variant="ghost" className="bg-white text-black hover:bg-gray-200 cursor-pointer" size="sm" onClick={() => setOpen((o) => !o)}>
          {open ? <ChevronDown className="h-4 w-4"/> : <ChevronRight className="h-4 w-4"/>} 
        </Button>
      </div>
      <div className="text-zinc-300 mt-2">
        <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>{String(body)}</ReactMarkdown>
      </div>
      {open && (
        <div className="mt-3 border-t border-zinc-800 pt-3">
          <div className="text-xs text-zinc-400 mb-1">Solution</div>
          <div className="text-zinc-200">
            <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>{String(solution)}</ReactMarkdown>
          </div>
        </div>
      )}
    </div>
  );
}

// ======================================================================
// Utilities: safe function makers (1D & system), CSV/PNG export
// ======================================================================
function makeFunctionSafe1D(expr) {
  try {
    const clean = String(expr).replace(/\^/g, "**");
    // eslint-disable-next-line no-new-func
    const fn = new Function("t", "y", `with (Math) { return (${clean}); }`);
    return (t, y) => {
      try {
        return Number(fn(t, y));
      } catch {
        return NaN;
      }
    };
  } catch {
    return () => NaN;
  }
}

function makeSystemSafe(exprX, exprY) {
  // exprX, exprY are strings using x,y,t and Math
  try {
    const cX = String(exprX).replace(/\^/g, "**");
    const cY = String(exprY).replace(/\^/g, "**");
    // eslint-disable-next-line no-new-func
    const fx = new Function("t", "x", "y", `with (Math) { return (${cX}); }`);
    // eslint-disable-next-line no-new-func
    const fy = new Function("t", "x", "y", `with (Math) { return (${cY}); }`);
    return (t, x, y) => {
      try {
        return [Number(fx(t, x, y)), Number(fy(t, x, y))];
      } catch {
        return [NaN, NaN];
      }
    };
  } catch {
    return () => [NaN, NaN];
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
      const imgData = await window.Plotly.toImage(document.getElementById(divId), { format: "png", width: 1200, height: 800 });
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

// ======================================================================
// ODE Solvers (explicit/implicit integrators, RK methods, adaptive RK)
// ======================================================================

// Explicit Euler: y_{n+1} = y_n + h f(t_n, y_n)
function explicitEuler(f, t0, y0, h, nSteps) {
  const ts = new Array(nSteps + 1);
  const ys = new Array(nSteps + 1);
  ts[0] = t0; ys[0] = y0;
  for (let i = 0; i < nSteps; i++) {
    const t = ts[i];
    const y = ys[i];
    const k1 = f(t, y);
    const y1 = y + h * k1;
    ts[i + 1] = t + h;
    ys[i + 1] = y1;
  }
  return { t: ts, y: ys };
}

// Explicit Midpoint (2nd order)
function explicitMidpoint(f, t0, y0, h, nSteps) {
  const ts = new Array(nSteps + 1);
  const ys = new Array(nSteps + 1);
  ts[0] = t0; ys[0] = y0;
  for (let i = 0; i < nSteps; i++) {
    const t = ts[i];
    const y = ys[i];
    const k1 = f(t, y);
    const yMid = y + 0.5 * h * k1;
    const k2 = f(t + 0.5 * h, yMid);
    const y1 = y + h * k2;
    ts[i + 1] = t + h;
    ys[i + 1] = y1;
  }
  return { t: ts, y: ys };
}

// RK4 classical
function rk4(f, t0, y0, h, nSteps) {
  const ts = new Array(nSteps + 1);
  const ys = new Array(nSteps + 1);
  ts[0] = t0; ys[0] = y0;
  for (let i = 0; i < nSteps; i++) {
    const t = ts[i];
    const y = ys[i];
    const k1 = f(t, y);
    const k2 = f(t + 0.5 * h, y + 0.5 * h * k1);
    const k3 = f(t + 0.5 * h, y + 0.5 * h * k2);
    const k4 = f(t + h, y + h * k3);
    const y1 = y + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
    ts[i + 1] = t + h;
    ys[i + 1] = y1;
  }
  return { t: ts, y: ys };
}

// Simple adaptive RK45 (Dormand-Prince-ish like, but simplified)
// We'll implement a pair of 4th/5th order coefficients for error estimation (classic DP coefficients)
function rk45Adaptive(f, t0, y0, h0, tEnd, tol = 1e-6, hMin = 1e-8, hMax = 1.0) {
  const ts = [];
  const ys = [];
  let t = t0;
  let y = y0;
  let h = h0;
  ts.push(t);
  ys.push(y);

  // coefficients (Dormand-Prince-ish simplified)
  while (t < tEnd - 1e-14) {
    if (t + h > tEnd) h = tEnd - t;
    // use RK4 as base and embedded lower order for error (this is a simplified controller)
    // compute RK4 step
    const k1 = f(t, y);
    const k2 = f(t + 0.5 * h, y + 0.5 * h * k1);
    const k3 = f(t + 0.5 * h, y + 0.5 * h * k2);
    const k4 = f(t + h, y + h * k3);
    const y4 = y + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4);

    // take two half RK4 steps for error estimate
    const h2 = h / 2;
    const k1a = f(t, y);
    const k2a = f(t + 0.5 * h2, y + 0.5 * h2 * k1a);
    const k3a = f(t + 0.5 * h2, y + 0.5 * h2 * k2a);
    const k4a = f(t + h2, y + h2 * k3a);
    const yHalf = y + (h2 / 6) * (k1a + 2 * k2a + 2 * k3a + k4a);

    const k1b = f(t + h2, yHalf);
    const k2b = f(t + h2 + 0.5 * h2, yHalf + 0.5 * h2 * k1b);
    const k3b = f(t + h2 + 0.5 * h2, yHalf + 0.5 * h2 * k2b);
    const k4b = f(t + h, yHalf + h2 * k3b);
    const yFull = yHalf + (h2 / 6) * (k1b + 2 * k2b + 2 * k3b + k4b);

    const err = Math.abs(yFull - y4);
    const safety = 0.9;
    const eps = tol;

    // Accept if error small
    if (err <= eps || h <= hMin) {
      t = t + h;
      y = yFull;
      ts.push(t);
      ys.push(y);
      // adjust step
      const factor = err === 0 ? 2 : Math.min(5, Math.max(0.2, safety * Math.pow(eps / (err + 1e-16), 0.25)));
      h = Math.min(hMax, Math.max(hMin, h * factor));
    } else {
      // reduce h
      const factor = Math.max(0.1, safety * Math.pow(eps / (err + 1e-16), 0.25));
      h = Math.max(hMin, h * factor);
    }
  }

  return { t: ts, y: ys };
}

// Simple fixed-step solver for systems (explicit RK4)
function rk4System(sys, t0, y0, h, nSteps) {
  // y0 is array [x0, y0]
  const ts = new Array(nSteps + 1);
  const xs = new Array(nSteps + 1);
  const ys = new Array(nSteps + 1);
  ts[0] = t0; xs[0] = y0[0]; ys[0] = y0[1];
  for (let i = 0; i < nSteps; i++) {
    const t = ts[i];
    const x = xs[i];
    const y = ys[i];
    const k1 = sys(t, x, y);
    const k2 = sys(t + 0.5 * h, x + 0.5 * h * k1[0], y + 0.5 * h * k1[1]);
    const k3 = sys(t + 0.5 * h, x + 0.5 * h * k2[0], y + 0.5 * h * k2[1]);
    const k4 = sys(t + h, x + h * k3[0], y + h * k3[1]);
    const x1 = x + (h / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
    const y1 = y + (h / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
    ts[i + 1] = t + h;
    xs[i + 1] = x1;
    ys[i + 1] = y1;
  }
  return { t: ts, x: xs, y: ys };
}

// ======================================================================
// Responsive Plot wrapper
// Ensures plots adapt to container size and do not overflow
// ======================================================================
function ResponsivePlot({ plotId, data, title, layout = {}, config = {}, style = {} }) {
  return (
    <div id={plotId} className="rounded-xl border border-zinc-700 p-2 bg-zinc-900/40 overflow-hidden w-full" style={{ maxWidth: "100%" }}>
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
        style={{ width: "100%", height: "clamp(280px, 52vh, 540px)", ...style }}
      />
    </div>
  );
}

// ======================================================================
// Sections
// - 19.1 Single-step integrators (Euler, Midpoint)
// - 19.2 RK4 & fixed-step convergence
// - 19.3 Adaptive RK (rk45 simplified)
// - 19.4 Systems & phase portraits (2D)
// - 19.5 Stability & stiff example (implicit Euler simple solver demo)
// ======================================================================

// ---------------- Section 19.1: Euler & Midpoint -----------------------
function Section191() {
  const [expr, setExpr] = useState("y - t^2 + 1"); // dy/dt = y - t^2 + 1
  const [t0, setT0] = useState("0");
  const [y0, setY0] = useState("0.5");
  const [tEnd, setTEnd] = useState("2");
  const [n, setN] = useState(20);
  const [method, setMethod] = useState("euler"); // 'euler' or 'midpoint'
  const plotId = "plot-19-1";

  const f = useMemo(() => makeFunctionSafe1D(expr), [expr]);
  const t0Num = Number(t0);
  const y0Num = Number(y0);
  const tEndNum = Number(tEnd);
  const nNum = Math.max(1, Math.floor(Number(n)));
  const h = (tEndNum - t0Num) / nNum;

  // numerical solution
  const sol = useMemo(() => {
    if (method === "euler") return explicitEuler((tt, yy) => f(tt, yy), t0Num, y0Num, h, nNum);
    return explicitMidpoint((tt, yy) => f(tt, yy), t0Num, y0Num, h, nNum);
  }, [f, t0Num, y0Num, h, nNum, method]);

  // exact solution for test equation (if known). For y' = y - t^2 + 1 with y(0)=0.5 exact is:
  // y = (t + 1)^2 - 0.5 e^{t} (??) Actually known exact: y(t) = (t^2 + 2t + 2 - 0.5 * e^{t})? But avoid assuming.
  // We'll optionally allow user to provide exact solution string.
  const [exactExpr, setExactExpr] = useState("");
  const exactFn = useMemo(() => {
    if (!exactExpr || exactExpr.trim() === "") return null;
    try {
      // eslint-disable-next-line no-new-func
      const fn = new Function("t", `with (Math) { return (${String(exactExpr).replace(/\^/g, "**")}); }`);
      return (t) => {
        try { return Number(fn(t)); } catch { return NaN; }
      };
    } catch { return null; }
  }, [exactExpr]);

  const xsPlot = useMemo(() => linspace(t0Num, tEndNum, 400), [t0Num, tEndNum]);
  const ysPlot = useMemo(() => xsPlot.map((tt) => {
    const val = f(tt, undefined); // if user uses t only
    return Number.isFinite(val) ? val : NaN;
  }), [xsPlot, f]);

  function handleExportCSV() {
    const rows = [["t", "numerical y", "maybe exact"]];
    for (let i = 0; i < sol.t.length; i++) {
      const t = sol.t[i];
      const y = sol.y[i];
      const ye = exactFn ? exactFn(t) : "";
      rows.push([t, y, ye]);
    }
    downloadCSV("19_1_euler_midpoint.csv", rows);
  }
  function handleExportPNG() {
    exportPlotPNG(plotId, "19_1_euler_midpoint.png");
  }

  // Prepare Plot traces: numeric and optionally exact (if provided)
  const traces = useMemo(() => {
    const numericTrace = { x: sol.t, y: sol.y, mode: "lines+markers", name: `${method} (h=${fmt(h,6)})`, line: { color: theme.accent } };
    const exactTrace = exactFn ? { x: linspace(t0Num, tEndNum, 400), y: linspace(t0Num, tEndNum, 400).map((tt) => exactFn(tt)), mode: "lines", name: "Exact", line: { dash: "dash", color: theme.accent2 } } : null;
    return exactTrace ? [exactTrace, numericTrace] : [numericTrace];
  }, [sol, exactFn, method, h, t0Num, tEndNum]);

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader><CardTitle className="text-cyan-300 flex items-center gap-2"><Activity className="w-5 h-5" />19.1 Single-step Integrators — Euler & Midpoint</CardTitle></CardHeader>
        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-4 gap-3">
            <Labeled label="dy/dt = f(t,y)"><Textarea value={expr} onChange={(e) => setExpr(e.target.value)} className="h-24 text-white" /></Labeled>
            <Labeled label="t0"><Input  className="text-white" value={t0} onChange={(e) => setT0(e.target.value)} /></Labeled>
            <Labeled label="y0"><Input  className="text-white" value={y0} onChange={(e) => setY0(e.target.value)} /></Labeled>
            <Labeled label="t end / n"><div className="flex gap-2"><Input  className="text-white" value={tEnd} onChange={(e) => setTEnd(e.target.value)} /><Input  value={n} onChange={(e) => setN(e.target.value)} className="w-28 text-white" /></div></Labeled>
          </div>

          <div className="flex items-center gap-3 flex-wrap">
            <div className="text-xs text-zinc-400">Method:</div>
            <button className={`px-3 py-1 rounded cursor-pointer ${method === "euler" ? "bg-emerald-600 text-white" : "bg-zinc-800 text-zinc-200"}`} onClick={() => setMethod("euler")}>Explicit Euler</button>
            <button className={`px-3 py-1 rounded cursor-pointer ${method === "midpoint" ? "bg-emerald-600 text-white" : "bg-zinc-800 text-zinc-200"}`} onClick={() => setMethod("midpoint")}>Midpoint</button>
            <div className="text-xs text-zinc-400 ml-3">Optional: paste exact y(t) if known to overlay.</div>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
            <ResponsivePlot plotId={plotId} data={traces} title={`Numerical solution (${method}) — approx final y ≈ ${fmt(sol.y[sol.y.length - 1], 8)}`} />

            <div className="rounded-xl border border-zinc-700 p-3 bg-zinc-900/40 overflow-auto">
              <div className="text-zinc-300">Numerical summary</div>
              <div className="text-2xl mt-2 text-emerald-300 font-mono">{fmt(sol.y[sol.y.length - 1], 8)}</div>

              <div className="mt-3 text-xs text-zinc-400">
                <div className="mb-2">Exact solution (optional):</div>
                <Textarea value={exactExpr} onChange={(e) => setExactExpr(e.target.value)} className="h-20 text-white" />
              </div>

              <div className="mt-3 flex gap-2">
                <Button  className="bg-white text-black hover:bg-gray-200 cursor-pointer" onClick={handleExportCSV} variant="ghost"><Download className="w-4 h-4" /> CSV</Button>
                <Button className="bg-white text-black hover:bg-gray-200 cursor-pointer" onClick={handleExportPNG} variant="ghost"><ImageIcon className="w-4 h-4" /> PNG</Button>
              </div>

              <div className="mt-4 text-xs text-zinc-400">
                <Exercise title="Exercise 19.1.1" body="Solve y' = y - t^2 + 1, y(0)=0.5 using Euler with n=20 and compare with RK4 with same step. Which is more accurate?" solution="RK4 will be significantly more accurate than Euler for moderate n because of higher order local truncation error." />
              </div>
            </div>
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ---------------- Section 19.2: RK4 & convergence -----------------------
function Section192() {
  const [expr, setExpr] = useState("-2*t*y"); // sample stiffish-like but simple
  const [t0, setT0] = useState("0");
  const [y0, setY0] = useState("1");
  const [tEnd, setTEnd] = useState("2");
  const [nBase, setNBase] = useState(40);
  const plotId = "plot-19-2";

  const f = useMemo(() => makeFunctionSafe1D(expr), [expr]);
  const t0Num = Number(t0);
  const y0Num = Number(y0);
  const tEndNum = Number(tEnd);
  const nBaseNum = Math.max(2, Math.floor(Number(nBase)));

  // compute RK4 for multiple step refinements to illustrate convergence
  const results = useMemo(() => {
    const cases = [];
    for (let k = 0; k < 6; k++) {
      const n = nBaseNum * (2 ** k);
      const h = (tEndNum - t0Num) / n;
      const sol = rk4((tt, yy) => f(tt, yy), t0Num, y0Num, h, n);
      cases.push({ n, h, t: sol.t, y: sol.y });
    }
    return cases;
  }, [f, t0Num, y0Num, tEndNum, nBaseNum]);

  // pick finest as "reference"
  const ref = results[results.length - 1];

  // compute error at tEnd for each
  const errors = useMemo(() => results.map((r) => {
    const yEnd = r.y[r.y.length - 1];
    const yRef = ref.y[ref.y.length - 1];
    return Math.abs(yEnd - yRef);
  }), [results, ref]);

  // prepare plot: overlay solutions
  const traces = useMemo(() => {
    const out = [];
    for (let i = 0; i < results.length; i++) {
      const r = results[i];
      out.push({ x: r.t, y: r.y, name: `RK4 n=${r.n}`, mode: "lines", line: { width: i === results.length - 1 ? 3 : 1, dash: i === results.length - 1 ? "solid" : "dot" } });
    }
    return out;
  }, [results]);

  function handleExportCSV() {
    // export a table of t,y for the finest simulation
    const rows = [["t", ...results.map((r) => `y_n${r.n}`)]];
    const maxLen = Math.max(...results.map((r) => r.t.length));
    for (let i = 0; i < maxLen; i++) {
      const row = [results[results.length - 1].t[i] || ""];
      for (let rIdx = 0; rIdx < results.length; rIdx++) {
        row.push(results[rIdx].y[i] || "");
      }
      rows.push(row);
    }
    downloadCSV("19_2_rk4_convergence.csv", rows);
  }
  function handleExportPNG() {
    exportPlotPNG(plotId, "19_2_rk4_convergence.png");
  }

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader><CardTitle className="text-cyan-300 flex items-center gap-2"><Layers className="w-5 h-5" />19.2 RK4 & Convergence</CardTitle></CardHeader>
        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-4 gap-3">
            <Labeled label="dy/dt = f(t,y)"><Textarea value={expr} onChange={(e) => setExpr(e.target.value)} className="h-24 text-white" /></Labeled>
            <Labeled label="t0"><Input className="text-white" value={t0} onChange={(e) => setT0(e.target.value)} /></Labeled>
            <Labeled label="y0"><Input className="text-white" value={y0} onChange={(e) => setY0(e.target.value)} /></Labeled>
            <Labeled label="t end / n base"><div className="flex gap-2"><Input className="text-white" value={tEnd} onChange={(e) => setTEnd(e.target.value)} /><Input value={nBase} onChange={(e) => setNBase(e.target.value)} className="w-28 text-white" /></div></Labeled>
          </div>

         <Equation>{`$$y_{n+1} = y_n + \\dfrac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4)$$`}</Equation>


          <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
            <ResponsivePlot plotId={plotId} data={traces} title={`RK4 solutions (ref n=${ref.n})`} />
            <div className="rounded-xl border border-zinc-700 p-3 bg-zinc-900/40 overflow-auto">
              <div className="text-zinc-300">Convergence summary (error at t_end vs finest)</div>
              <div className="mt-2 text-xs text-zinc-100 font-mono">
                {results.map((r, i) => (
                  <div key={r.n} className="flex items-center justify-between whitespace-nowrap">
                    <div>n = {r.n}, h = {fmt((tEndNum - t0Num) / r.n, 8)}</div>
                    <div className="text-emerald-300 font-mono">error ≈ {fmt(errors[i], 12)}</div>
                  </div>
                ))}
              </div>

              <div className="mt-3 flex gap-2">
                <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={handleExportCSV} variant="ghost"><Download className="w-4 h-4" /> CSV</Button>
                <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={handleExportPNG} variant="ghost"><ImageIcon className="w-4 h-4" /> PNG</Button>
              </div>

              <div className="mt-4 text-xs text-zinc-400">
                <Exercise title="Exercise 19.2.1" body="Demonstrate 4th-order convergence by halving h and tracking error decrease; compute observed order." solution="For RK4, error should reduce roughly by factor 16 when halving h (error ~ C h^4)." />
              </div>
            </div>
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ---------------- Section 19.3: Adaptive RK (RK45 simplified) -----------------------
function Section193() {
  const [expr, setExpr] = useState("t*Math.sin(y) - 0.5*y");
  const [t0, setT0] = useState("0");
  const [y0, setY0] = useState("1");
  const [tEnd, setTEnd] = useState("10");
  const [h0, setH0] = useState("0.1");
  const [tol, setTol] = useState("1e-5");
  const plotId = "plot-19-3";

  const f = useMemo(() => makeFunctionSafe1D(expr), [expr]);
  const t0Num = Number(t0), y0Num = Number(y0), tEndNum = Number(tEnd), h0Num = Math.max(1e-8, Number(h0));
  const tolNum = Math.max(1e-12, Number(tol));

  const sol = useMemo(() => rk45Adaptive((tt, yy) => f(tt, yy), t0Num, y0Num, h0Num, tEndNum, tolNum, 1e-8, 1.0), [f, t0Num, y0Num, h0Num, tEndNum, tolNum]);

  function handleExportCSV() {
    const rows = [["t", "y"]];
    for (let i = 0; i < sol.t.length; i++) rows.push([sol.t[i], sol.y[i]]);
    downloadCSV("19_3_adaptive_rk.csv", rows);
  }
  function handleExportPNG() {
    exportPlotPNG(plotId, "19_3_adaptive_rk.png");
  }

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader><CardTitle className="text-cyan-300 flex items-center gap-2"><Zap className="w-5 h-5" />19.3 Adaptive RK (RK45 simplified)</CardTitle></CardHeader>
        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 md:grid-cols-4 gap-3">
            <Labeled label="dy/dt = f(t,y)"><Textarea value={expr} onChange={(e) => setExpr(e.target.value)} className="h-24 text-white" /></Labeled>
            <Labeled label="t0"><Input className="text-white" value={t0} onChange={(e) => setT0(e.target.value)} /></Labeled>
            <Labeled label="y0"><Input  className="text-white" value={y0} onChange={(e) => setY0(e.target.value)} /></Labeled>
            <Labeled label="t end / h0 / tol"><div className="flex gap-2"><Input  className="text-white" value={tEnd} onChange={(e) => setTEnd(e.target.value)} /><Input value={h0} onChange={(e) => setH0(e.target.value)} className="w-28 text-white" /><Input value={tol} onChange={(e) => setTol(e.target.value)} className="w-28 text-white" /></div></Labeled>
          </div>

          <Equation>{`Adaptive RK uses embedded lower-order method to estimate local error and adjust step size to meet tolerance.`}</Equation>

          <div className="grid md:grid-cols-2 gap-3">
            <ResponsivePlot plotId={plotId} data={[{ x: sol.t, y: sol.y, mode: "lines+markers", name: "y(t)", line: { color: theme.accent } }]} title={`Adaptive RK (steps: ${sol.t.length})`} />

            <div className="rounded-xl border border-zinc-700 p-3 bg-zinc-900/40 overflow-auto">
              <div className="text-zinc-300">Adaptive integration summary</div>
              <div className="text-2xl mt-2 text-emerald-300 font-mono">{fmt(sol.y[sol.y.length - 1], 8)}</div>

              <div className="mt-3 flex gap-2">
                <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={handleExportCSV} variant="ghost"><Download className="w-4 h-4" /> CSV</Button>
                <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={handleExportPNG} variant="ghost"><ImageIcon className="w-4 h-4" /> PNG</Button>
              </div>

              <div className="mt-4 text-xs text-zinc-400">
                <Exercise title="Exercise 19.3.1" body="Use adaptive RK to integrate a problem with rapid transients (e.g., y' = -1000(y - cos(t))). Observe step size adaptation near transients." solution="Adaptive scheme reduces step size near sharp features and increases it elsewhere to maintain tol with fewer function evaluations overall." />
              </div>
            </div>
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ---------------- Section 19.4: Systems & Phase Portraits -----------------------
function Section194() {
  const [exprX, setExprX] = useState("y"); // dx/dt = y
  const [exprY, setExprY] = useState("-x"); // dy/dt = -x  --> harmonic oscillator
  const [t0, setT0] = useState("0");
  const [x0, setX0] = useState("1");
  const [y0, setY0] = useState("0");
  const [tEnd, setTEnd] = useState("12.566"); // 4*pi
  const [n, setN] = useState(400);
  const plotIdPhase = "plot-19-phase";
  const plotIdXY = "plot-19-xy";

  const sys = useMemo(() => makeSystemSafe(exprX, exprY), [exprX, exprY]);
  const t0Num = Number(t0), x0Num = Number(x0), y0Num = Number(y0), tEndNum = Number(tEnd), nNum = Math.max(4, Math.floor(Number(n)));
  const h = (tEndNum - t0Num) / nNum;

  const sol = useMemo(() => rk4System(sys, t0Num, [x0Num, y0Num], h, nNum), [sys, t0Num, x0Num, y0Num, h, nNum]);

  // Phase portrait vector field (grid)
  const vf = useMemo(() => {
    // create grid around trajectory extents
    const xs = sol.x;
    const ys = sol.y;
    const xmin = Math.min(...xs) - 0.5, xmax = Math.max(...xs) + 0.5;
    const ymin = Math.min(...ys) - 0.5, ymax = Math.max(...ys) + 0.5;
    const gx = linspace(xmin, xmax, 20);
    const gy = linspace(ymin, ymax, 20);
    const vecs = [];
    for (let i = 0; i < gx.length; i++) {
      for (let j = 0; j < gy.length; j++) {
        const x = gx[i], y = gy[j];
        const [dx, dy] = sys(0, x, y);
        vecs.push({ x, y, dx, dy });
      }
    }
    return { gx, gy, vecs, xmin, xmax, ymin, ymax };
  }, [sol, sys]);

  function handleExportCSV() {
    const rows = [["t", "x", "y"]];
    for (let i = 0; i < sol.t.length; i++) rows.push([sol.t[i], sol.x[i], sol.y[i]]);
    downloadCSV("19_4_phase_xy.csv", rows);
  }
  function handleExportPNGPhase() { exportPlotPNG(plotIdPhase, "19_4_phase.png"); }
  function handleExportPNGXY() { exportPlotPNG(plotIdXY, "19_4_xy.png"); }

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader><CardTitle className="text-cyan-300 flex items-center gap-2"><FunctionSquare className="w-5 h-5" />19.4 Systems & Phase Portraits</CardTitle></CardHeader>
        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 md:grid-cols-6 gap-3">
            <Labeled label="dx/dt"><Input className="text-white" value={exprX} onChange={(e) => setExprX(e.target.value)} /></Labeled>
            <Labeled label="dy/dt"><Input className="text-white" value={exprY} onChange={(e) => setExprY(e.target.value)} /></Labeled>
            <Labeled label="t0"><Input className="text-white" value={t0} onChange={(e) => setT0(e.target.value)} /></Labeled>
            <Labeled label="x0"><Input className="text-white" value={x0} onChange={(e) => setX0(e.target.value)} /></Labeled>
            <Labeled label="y0"><Input className="text-white" value={y0} onChange={(e) => setY0(e.target.value)} /></Labeled>
            <Labeled label="tend / n"><div className="flex gap-2"><Input className="text-white" value={tEnd} onChange={(e) => setTEnd(e.target.value)} /><Input value={n} onChange={(e) => setN(e.target.value)} className="w-28 text-white" /></div></Labeled>
          </div>

          <Equation>{`System: $\\dot{x} = f(t,x,y), \\; \\dot{y} = g(t,x,y)$. Phase portrait shows vector field and trajectories.`}</Equation>


          <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
            <ResponsivePlot plotId={plotIdPhase} data={[
              { x: vf.vecs.map((v) => v.x), y: vf.vecs.map((v) => v.y), mode: "markers", marker: { size: 1 }, name: "grid points", showlegend: false },
              { x: sol.x, y: sol.y, mode: "lines", name: "trajectory", line: { color: theme.accent, width: 2 } },
              // overlay arrows as quiver-like using segments
              ...vf.vecs.slice(0, 400).map((v, i) => ({ x: [v.x, v.x + 0.15 * v.dx], y: [v.y, v.y + 0.15 * v.dy], mode: "lines", line: { width: 1, color: "rgba(255,255,255,0.4)" }, showlegend: false }))
            ]} title="Phase portrait & trajectory" />

            <ResponsivePlot plotId={plotIdXY} data={[
              { x: sol.t, y: sol.x, mode: "lines", name: "x(t)", line: { color: theme.accent } },
              { x: sol.t, y: sol.y, mode: "lines", name: "y(t)", line: { color: theme.accent2 } },
            ]} title="State components over time" />
          </div>

          <div className="mt-2 flex gap-2">
            <Button onClick={handleExportCSV} variant="ghost"><Download className="w-4 h-4" /> CSV</Button>
            <Button onClick={handleExportPNGPhase} variant="ghost"><ImageIcon className="w-4 h-4" /> PNG Phase</Button>
            <Button onClick={handleExportPNGXY} variant="ghost"><ImageIcon className="w-4 h-4" /> PNG x,y</Button>
          </div>

          <div className="mt-4 text-xs text-zinc-400">
            <Exercise title="Exercise 19.4.1" body="Explore trajectories for dx/dt = y, dy/dt = -x + 0.1(1 - x^2) y (Van der Pol like). Plot phase portrait and comment on limit cycles." solution="Adjust parameters and observe emergence of limit cycles when nonlinearity is present; RK4 with sufficient resolution visualizes the orbit." />
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ---------------- Section 19.5: Stability & implicit Euler demo -----------------------
// We'll implement a simple implicit Euler using fixed-point iteration (simple demo)
function Section195() {
  const [expr, setExpr] = useState("-50 * (y - Math.cos(t))"); // stiff-like example
  const [t0, setT0] = useState("0");
  const [y0, setY0] = useState("0");
  const [tEnd, setTEnd] = useState("1");
  const [n, setN] = useState(40);
  const plotId = "plot-19-5";

  const f = useMemo(() => makeFunctionSafe1D(expr), [expr]);
  const t0Num = Number(t0), y0Num = Number(y0), tEndNum = Number(tEnd), nNum = Math.max(1, Math.floor(Number(n)));
  const h = (tEndNum - t0Num) / nNum;

  // explicit Euler (for comparison)
  const expl = useMemo(() => explicitEuler((t, y) => f(t, y), t0Num, y0Num, h, nNum), [f, t0Num, y0Num, h, nNum]);

  // implicit Euler using simple fixed point: y_{n+1} = y_n + h f(t_{n+1}, y_{n+1}) solve by iteration y_{k+1} = y_n + h f(t_{n+1}, y_k)
  function implicitEulerSolve(sysf, t0_, y0_, h_, nSteps, maxIter = 50, tolIter = 1e-10) {
    const ts = new Array(nSteps + 1);
    const ys = new Array(nSteps + 1);
    ts[0] = t0_; ys[0] = y0_;
    for (let i = 0; i < nSteps; i++) {
      const tnp1 = ts[i] + h_;
      const yPrev = ys[i];
      let yGuess = yPrev; // initial guess
      for (let k = 0; k < maxIter; k++) {
        const yNew = yPrev + h_ * sysf(tnp1, yGuess); // fixed-point
        if (Math.abs(yNew - yGuess) < tolIter) { yGuess = yNew; break; }
        yGuess = yNew;
      }
      ts[i + 1] = tnp1;
      ys[i + 1] = yGuess;
    }
    return { t: ts, y: ys };
  }

  const impl = useMemo(() => implicitEulerSolve((t, y) => f(t, y), t0Num, y0Num, h, nNum, 200, 1e-12), [f, t0Num, y0Num, h, nNum]);

  function handleExportCSV() {
    const rows = [["t", "explicit y", "implicit y"]];
    const L = Math.max(expl.t.length, impl.t.length);
    for (let i = 0; i < L; i++) {
      rows.push([expl.t[i] || "", expl.y[i] || "", impl.y[i] || ""]);
    }
    downloadCSV("19_5_implicit_explicit.csv", rows);
  }
  function handleExportPNG() { exportPlotPNG(plotId, "19_5_stability.png"); }

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader><CardTitle className="text-cyan-300 flex items-center gap-2"><Percent className="w-5 h-5" />19.5 Stability & Implicit Euler</CardTitle></CardHeader>
        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-4 gap-3">
            <Labeled label="dy/dt"><Textarea value={expr} onChange={(e) => setExpr(e.target.value)} className="h-24 text-white" /></Labeled>
            <Labeled label="t0"><Input className="text-white" value={t0} onChange={(e) => setT0(e.target.value)} /></Labeled>
            <Labeled label="y0"><Input className="text-white" value={y0} onChange={(e) => setY0(e.target.value)} /></Labeled>
            <Labeled label="tend / n"><div className="flex gap-2"><Input className="text-white" value={tEnd} onChange={(e) => setTEnd(e.target.value)} /><Input value={n} onChange={(e) => setN(e.target.value)} className="w-28 text-white" /></div></Labeled>
          </div>

         <Equation>{`Implicit Euler: $y_{n+1} = y_n + h f(t_{n+1}, y_{n+1})$`}</Equation>


          <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
            <ResponsivePlot plotId={plotId} data={[
              { x: expl.t, y: expl.y, mode: "lines+markers", name: "Explicit Euler", line: { color: "#f97316" } },
              { x: impl.t, y: impl.y, mode: "lines+markers", name: "Implicit Euler", line: { color: theme.accent2 } }
            ]} title="Explicit vs Implicit Euler — stability demonstration" />

            <div className="rounded-xl border border-zinc-700 p-3 bg-zinc-900/40 overflow-auto">
              <div className="text-zinc-300">Stability notes</div>
              <div className="mt-2 text-xs text-zinc-100">
                Implicit methods (like backward Euler) are A-stable for many linear problems and handle stiff ODEs with larger steps without instability, though they require solving nonlinear equations each step (iteration or Newton).
              </div>

              <div className="mt-3 flex gap-2">
                <Button 
className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black"
 onClick={handleExportCSV} variant="ghost"><Download className="w-4 h-4" /> CSV</Button>
                <Button 
className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black"
 onClick={handleExportPNG} variant="ghost"><ImageIcon className="w-4 h-4" /> PNG</Button>
              </div>

              <div className="mt-4 text-xs text-zinc-400">
                <Exercise title="Exercise 19.5.1" body="Test stiffness by integrating y' = -1000(y - cos(t)) with explicit and implicit Euler for h=0.01 and h=0.001. Observe instabilities with explicit method." solution="Explicit method will blow up unless h is very small; implicit remains stable with larger h though slightly dissipative." />
              </div>
            </div>
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// Documentation & problems for Chapter 19
// ======================================================================
const docs19 = {
  "19.1 Single-step methods": [
    "Single-step methods compute $y_{n+1}$ using only information from the current step $y_n$.",
    "$$y_{n+1} = y_n + h f(t_n, y_n) \\quad \\text{(Explicit Euler)}$$",
    "$$\\text{Midpoint: } y_{n+1} = y_n + h f\\left(t_n + \\frac{h}{2}, y_n + \\frac{h}{2} f(t_n, y_n)\\right)$$",
    "Heun's method (improved Euler) uses an average of slopes:",
    "$$y_{n+1} = y_n + \\frac{h}{2} \\left(f(t_n, y_n) + f(t_{n+1}, y_n + h f(t_n, y_n))\\right)$$",
    "These methods are easy to implement but may require small step sizes for stability."
  ],
  "19.2 RK4 & convergence": [
    "Classical RK4 is a 4th-order single-step method with good accuracy for smooth problems.",
    "$$\\begin{aligned} k_1 &= f(t_n, y_n), \\\\ k_2 &= f(t_n + \\frac{h}{2}, y_n + \\frac{h}{2}k_1), \\\\ k_3 &= f(t_n + \\frac{h}{2}, y_n + \\frac{h}{2}k_2), \\\\ k_4 &= f(t_n + h, y_n + hk_3), \\\\ y_{n+1} &= y_n + \\frac{h}{6} (k_1 + 2k_2 + 2k_3 + k_4) \\end{aligned}$$",
    "RK4 achieves global error $O(h^4)$ and is widely used in engineering and physics.",
    "Convergence analysis ensures that as $h \\to 0$, the numerical solution approaches the exact solution."
  ],
  "19.3 Adaptive RK": [
    "Adaptive Runge-Kutta methods adjust step size $h$ to maintain a desired error tolerance.",
    "Embedded RK pairs (e.g., Dormand-Prince, Fehlberg) compute two solutions of different orders:",
    "$$y_{n+1}^{(p)} \\text{ and } y_{n+1}^{(p-1)}$$",
    "The local error estimate is:",
    "$$\\text{Error} = \\lVert y_{n+1}^{(p)} - y_{n+1}^{(p-1)} \\rVert$$",
    "If the error exceeds tolerance, $h$ is reduced; if the error is small, $h$ can be increased.",
    "Adaptive RK is especially useful for stiff or highly varying problems."
  ],
  "19.4 Systems & phase portraits": [
    "For systems of ODEs $\\dot{x} = f(t,x,y), \\; \\dot{y} = g(t,x,y)$, the phase plane shows trajectories.",
    "$$\\text{System: } \\dot{x} = f(t,x,y), \\; \\dot{y} = g(t,x,y)$$",
    "Phase portraits reveal fixed points, limit cycles, and stability regions.",
    "Vector fields indicate direction of motion; trajectories show solution paths.",
    "Linearization near equilibrium points helps classify stability (node, saddle, spiral)."
  ],
  "19.5 Stability": [
    "Stability refers to how numerical solutions behave for long time integration.",
    "Explicit methods (Euler, RK) may be unstable for stiff ODEs.",
    "Implicit methods (backward Euler, BDF) are often A-stable:",
    "$$y_{n+1} = y_n + h f(t_{n+1}, y_{n+1})$$",
    "A-stable methods remain stable for all step sizes for linear test equation $y' = \\lambda y, \\; \\text{Re}(\\lambda) < 0$.",
    "Stability analysis is crucial for stiff systems, chemical kinetics, and control problems."
  ],
};


function ChapterDocs19() {
  const [open, setOpen] = useState("19.1 Single-step methods");
  return (
    <div className="mt-6 bg-zinc-900 border border-zinc-700 rounded-xl p-4 overflow-hidden">
      <div className="flex items-center justify-between mb-3 gap-2 flex-wrap">
        <div>
          <div className="text-zinc-100 font-semibold">Chapter 19 — Notes</div>
          <div className="text-zinc-400 text-xs">Methods, stability, and tips</div>
        </div>
        <div className="text-zinc-400 text-xs">Initial value problems</div>
      </div>

      <div className="grid grid-cols-1 md:grid-cols-4">
        <div className="border-b md:border-b-0 md:border-r border-zinc-800 pr-0 md:pr-2">
          {Object.keys(docs19).map((k) => (
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
            {docs19[open].map((p, i) => (
              <ReactMarkdown key={i} remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>{p}</ReactMarkdown>
            ))}
          </div>
        </div>
      </div>
    </div>
  );
}

function Problems19() {
  const problems = [
    { title: "P19.1", body: "Use explicit Euler and RK4 on y' = y - t^2 + 1 with y(0)=0.5, compare global errors for n=10,20,40." },
    { title: "P19.2", body: "Implement adaptive RK and integrate a stiff test problem: y' = -1000(y - cos(t)). Compare step sizes." },
    { title: "P19.3", body: "Plot phase portrait for Van der Pol oscillator for mu=0.5 and mu=5 and discuss behavior." },
    { title: "P19.4", body: "Compare explicit and implicit Euler on y' = -1000 y with various step sizes to study stability." },
  ];
  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader><CardTitle className="text-cyan-300 flex items-center gap-2"><Grid className="w-5 h-5" />Problems</CardTitle></CardHeader>
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
// Page assembly
// ======================================================================
export default function Chapter19() {
  const topRef = useRef(null);
  useEffect(() => { if (topRef.current) topRef.current.focus?.(); }, []);
  return (
    <div className={`p-4 md:p-6 lg:p-8 ${theme.bg} min-h-screen`}>
      <div className="max-w-screen-xl mx-auto space-y-6">
        <motion.div initial={{ opacity: 0, y: -6 }} animate={{ opacity: 1, y: 0 }}>
          <h1 ref={topRef} tabIndex={-1} className="text-2xl md:text-3xl lg:text-4xl font-bold" style={{ color: theme.accent2 }}>Numerical ODEs & Initial Value Problems</h1>
          <div className="text-zinc-400">Explicit/implicit integrators, Runge-Kutta methods, adaptive control, phase portraits, and stability — interactive visualizers.</div>
        </motion.div>

        <Section191 />
        <Section192 />
        <Section193 />
        <Section194 />
        <Section195 />

        <ChapterDocs19 />
        <Problems19 />
         <BottomBar/>
      </div>
    </div>
  );
}
