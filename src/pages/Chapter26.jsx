// src/pages/Chapter26.jsx
// ======================================================================
// Chapter 26 — Stiffness and Multistep Methods (responsive, expanded)
// - Sections:
//   26.1 Stiffness (test equations, stability regions, explicit vs implicit)
//   26.2 Multistep Methods (AB2, AB3, AM2 predictor-corrector, comparisons)
//   Problems
// - Fully responsive layout (mobile-first), multiple plots per section
// - KaTeX rendering for formulas via local Equation component
// - Plotly plots use useResizeHandler + responsive config
// - Drop into src/pages/Chapter26.jsx (adjust imports if needed)
// ======================================================================

import React, { useEffect, useMemo, useRef, useState } from "react";
import Plot from "react-plotly.js";
import { motion } from "framer-motion";
import katex from "katex";
import "katex/dist/katex.min.css";
import ReactMarkdown from "react-markdown";
import remarkMath from "remark-math";
import rehypeKatex from "rehype-katex";

import { Card, CardHeader, CardContent, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Textarea } from "@/components/ui/textarea";
import { Button } from "@/components/ui/button";

import {
  Zap,
  Layers,
  Activity,
  BookOpen,
  List,
  ChevronRight,
  ChevronDown,
  FileText,
  RefreshCcw,
  Download,
  ImageIcon,
  Grid,
} from "lucide-react";
import BottomBar from "../components/BottomBar";

// ======================================================================
// Theme & helpers
// ======================================================================

const theme = {
  bg: "bg-zinc-950 text-zinc-200",
  panel: "bg-zinc-900/60 border border-zinc-700",
  accent: "#22d3ee",
  accent2: "#34d399",
  warn: "#f59e0b",
};

const fadeUp = {
  initial: { opacity: 0, y: 10 },
  animate: { opacity: 1, y: 0 },
  transition: { duration: 0.28 },
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
      const imgData = await window.Plotly.toImage(document.getElementById(divId), { format: "png", width: 1400, height: 900 });
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
// KaTeX equation renderer
// ======================================================================
function Equation({ latex, display = true, className = "" }) {
  const ref = useRef(null);
  useEffect(() => {
    if (!ref.current) return;
    try {
      katex.render(latex, ref.current, { throwOnError: false, displayMode: display });
    } catch (err) {
      ref.current.innerText = latex;
    }
  }, [latex, display]);
  return <div ref={ref} className={className} />;
}

// Provide a small Formula wrapper to avoid errors if used elsewhere
function Formula({ children }) {
  // children may be a string or JSX; we want a LaTeX string
  const latex = typeof children === "string" ? children : String(children);
  return <Equation latex={latex} display={true} />;
}

// ======================================================================
// Numerical routines (kept and adapted for clarity)
// ======================================================================

// scalar/vector-friendly euler
function eulerSolver(f, t0, y0, h, nSteps) {
  const ts = [t0];
  const ys = [Array.isArray(y0) ? [...y0] : y0];
  let t = t0;
  let y = Array.isArray(y0) ? [...y0] : y0;
  for (let i = 0; i < nSteps; i++) {
    const fy = f(t, y);
    if (Array.isArray(y)) {
      y = y.map((v, j) => v + h * fy[j]);
    } else {
      y = y + h * fy;
    }
    t = t + h;
    ts.push(t);
    ys.push(Array.isArray(y) ? [...y] : y);
  }
  return { t: ts, y: ys };
}

// RK4 scalar
function rk4Solver(f, t0, y0, h, nSteps) {
  const ts = [t0];
  const ys = [Array.isArray(y0) ? [...y0] : y0];
  let t = t0;
  let y = Array.isArray(y0) ? [...y0] : y0;
  for (let i = 0; i < nSteps; i++) {
    if (Array.isArray(y)) {
      const add = (a, b) => a.map((v, j) => v + b[j]);
      const scale = (a, s) => a.map((v) => v * s);
      const k1 = f(t, y);
      const k2 = f(t + h / 2, add(y, scale(k1, h / 2)));
      const k3 = f(t + h / 2, add(y, scale(k2, h / 2)));
      const k4 = f(t + h, add(y, scale(k3, h)));
      y = y.map((val, j) => val + (h / 6) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]));
    } else {
      const k1 = f(t, y);
      const k2 = f(t + h / 2, y + (h / 2) * k1);
      const k3 = f(t + h / 2, y + (h / 2) * k2);
      const k4 = f(t + h, y + h * k3);
      y = y + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
    }
    t = t + h;
    ts.push(t);
    ys.push(Array.isArray(y) ? [...y] : y);
  }
  return { t: ts, y: ys };
}

// Adams-Bashforth 2-step
function adamsBashforth2(f, t0, y0, h, nSteps) {
  const ts = [t0];
  const ys = [Array.isArray(y0) ? [...y0] : y0];
  // initialize with one RK4 step
  const s = rk4Solver(f, t0, y0, h, 1);
  ts.push(s.t[1]);
  ys.push(s.y[1]);
  for (let i = 1; i < nSteps; i++) {
    const tn = ts[i];
    const tn1 = ts[i - 1];
    const yn = ys[i];
    const yn1 = ys[i - 1];
    const fn = f(tn, yn);
    const fn1 = f(tn1, yn1);
    if (Array.isArray(yn)) {
      const next = yn.map((v, j) => v + (h / 2) * (3 * fn[j] - fn1[j]));
      ts.push(tn + h);
      ys.push(next);
    } else {
      const next = yn + (h / 2) * (3 * fn - fn1);
      ts.push(tn + h);
      ys.push(next);
    }
  }
  return { t: ts, y: ys };
}

// Adams-Bashforth 3-step
function adamsBashforth3(f, t0, y0, h, nSteps) {
  const start = rk4Solver(f, t0, y0, h, 2);
  const ts = start.t.slice();
  const ys = start.y.slice();
  for (let i = 2; i < nSteps; i++) {
    const tn = ts[i], yn = ys[i];
    const tn1 = ts[i - 1], yn1 = ys[i - 1];
    const tn2 = ts[i - 2], yn2 = ys[i - 2];
    const fn = f(tn, yn), fn1 = f(tn1, yn1), fn2 = f(tn2, yn2);
    if (Array.isArray(yn)) {
      const next = yn.map((v, j) => v + (h / 12) * (23 * fn[j] - 16 * fn1[j] + 5 * fn2[j]));
      ts.push(tn + h);
      ys.push(next);
    } else {
      const next = yn + (h / 12) * (23 * fn - 16 * fn1 + 5 * fn2);
      ts.push(tn + h);
      ys.push(next);
    }
  }
  return { t: ts, y: ys };
}

// Adams-Moulton 2 predictor-corrector
function adamsMoulton2_pc(f, t0, y0, h, nSteps, iters = 3) {
  const ts = [t0];
  const ys = [Array.isArray(y0) ? [...y0] : y0];
  const s = rk4Solver(f, t0, y0, h, 1);
  ts.push(s.t[1]);
  ys.push(s.y[1]);
  for (let i = 1; i < nSteps; i++) {
    const tn = ts[i], yn = ys[i];
    const tn1 = ts[i - 1], yn1 = ys[i - 1];
    const fn = f(tn, yn), fn1 = f(tn1, yn1);
    // predictor AB2
    let ypred = Array.isArray(yn) ? yn.map((v, j) => v + (h / 2) * (3 * fn[j] - fn1[j])) : yn + (h / 2) * (3 * fn - fn1);
    for (let k = 0; k < iters; k++) {
      const fpred = f(tn + h, ypred);
      if (Array.isArray(yn)) {
        ypred = yn.map((v, j) => v + (h / 2) * (fn[j] + fpred[j]));
      } else {
        ypred = yn + (h / 2) * (fn + fpred);
      }
    }
    ts.push(tn + h);
    ys.push(Array.isArray(ypred) ? [...ypred] : ypred);
  }
  return { t: ts, y: ys };
}

// ======================================================================
// Reusable UI bits
// ======================================================================
function SmallLabel({ label, children }) {
  return (
    <div className="min-w-0">
      <div className="text-xs text-zinc-400 mb-1">{label}</div>
      {children}
    </div>
  );
}
function CollapsibleExercise({ title, body, solution }) {
  const [open, setOpen] = useState(false);
  return (
    <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-3">
      <div className="flex items-center justify-between">
        <div className="text-zinc-100  font-medium">{title}</div>
        <Button variant="ghost" className="bg-white cursor-pointer text-black hover:bg-gray-200" size="sm" onClick={() => setOpen((o) => !o)}>
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

function Note({ children }) {
  return <div className="p-3 bg-zinc-800/40 rounded-md text-sm text-zinc-300 border-l-4 border-emerald-500">{children}</div>;
}

function Exercise({ title, body, solution }) {
  const [open, setOpen] = useState(false);
  return (
    <div className="rounded-lg border border-zinc-700 bg-zinc-900/40 p-3">
      <div className="flex items-center justify-between">
        <div className="text-sm font-medium text-zinc-100">{title}</div>
        <Button variant="ghost" size="sm" onClick={() => setOpen((o) => !o)}>{open ? "Hide" : "Show"} solution</Button>
      </div>
      <div className="text-xs text-zinc-300 mt-2"><Equation latex={body} display={false} /></div>
      {open && <div className="mt-3 border-t border-zinc-800 pt-3 text-zinc-200"><Equation latex={solution} display={true} /></div>}
    </div>
  );
}

// ======================================================================
// Responsive Plot wrapper
// ======================================================================
function ResponsivePlot({ id, data, title, layout = {}, style = {}, config = {} }) {
  return (
    <div id={id} className="rounded-xl border border-zinc-700 bg-zinc-900/40 p-2 overflow-hidden w-full" style={{ maxWidth: "100%" }}>
      <Plot
        data={data}
        layout={{
          title,
          paper_bgcolor: "rgba(0,0,0,0)",
          plot_bgcolor: "rgba(0,0,0,0)",
          font: { color: "#e5e7eb" },
          autosize: true,
          margin: { l: 48, r: 20, t: 48, b: 48 },
          legend: { orientation: "h", x: 0, y: -0.22 },
          ...layout,
        }}
        config={{ responsive: true, ...config }}
        useResizeHandler
        style={{ width: "100%", height: "clamp(240px, 52vh, 680px)", ...style }}
      />
    </div>
  );
}

// ======================================================================
// Section 26.1 — Stiffness (expanded with extra plots)
// ======================================================================
function Section261() {
  // test equation y' = lambda * y with y(0)=1
  const [lambda, setLambda] = useState(-50.0);
  const [hExp, setHExp] = useState(0.04);
  const [hImp, setHImp] = useState(0.2);
  const t0 = 0;
  const tEnd = 1.0;

  const lambdaNum = Number(lambda) || -50;
  const hExpNum = Math.max(1e-8, Number(hExp) || 0.04);
  const hImpNum = Math.max(1e-8, Number(hImp) || 0.2);

  // exact
  const exact = useMemo(() => {
    return (t) => Math.exp(lambdaNum * t);
  }, [lambdaNum]);

  // samples for plotting exact solution
  const tExact = linspace(t0, tEnd, 400);
  const yExact = tExact.map((tt) => exact(tt));

  // explicit Euler
  const nExp = Math.max(1, Math.floor((tEnd - t0) / hExpNum));
  const eulerExp = useMemo(() => eulerSolver((t, y) => lambdaNum * y, t0, 1.0, hExpNum, nExp), [lambdaNum, hExpNum, nExp]);

  // implicit Euler analytic update for scalar: y_{n+1} = y_n / (1 - h lambda)
  const implicitEuler = useMemo(() => {
    const N = Math.max(1, Math.floor((tEnd - t0) / hImpNum));
    const ts = [t0];
    const ys = [1.0];
    let y = 1.0;
    let t = t0;
    for (let i = 0; i < N; i++) {
      y = y / (1 - hImpNum * lambdaNum);
      t = t + hImpNum;
      ts.push(t);
      ys.push(y);
    }
    return { t: ts, y: ys };
  }, [lambdaNum, hImpNum]);

  // stability region circle for explicit Euler: |1+z| < 1 -> disk centered at -1 radius 1
  const theta = linspace(0, 2 * Math.PI, 240);
  const circleX = theta.map((th) => -1 + Math.cos(th));
  const circleY = theta.map((th) => Math.sin(th));

  // sample lambda*h across range to show stability markers (real axis)
  const markerHs = linspace(0.01, Math.max(0.01, Math.abs(2 / lambdaNum) * 0.9), 8);
  const markerZ = markerHs.map((hh) => ({ h: hh, z: lambdaNum * hh }));

  // convergence: error vs h for both explicit and implicit
  const hs = [0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625];
  const conv = hs.map((hh) => {
    const N = Math.max(1, Math.floor((tEnd - t0) / hh));
    const e = eulerSolver((t, y) => lambdaNum * y, t0, 1.0, hh, N).y;
    const lastE = e[e.length - 1];
    const approxE = Array.isArray(lastE) ? lastE[0] : lastE;
    const errE = Math.abs(approxE - exact(tEnd));
    // implicit analytic
    let yB = 1.0;
    for (let i = 0; i < N; i++) yB = yB / (1 - hh * lambdaNum);
    const errB = Math.abs(yB - exact(tEnd));
    return { h: hh, errE, errB };
  });

  // additional plot: spectral-like view (real axis) showing allowed h for stability vs lambda
  const lambdaRange = linspace(lambdaNum * 0.2, lambdaNum * 1.2, 200);

  const idSol = "c26-1-solution";
  const idErr = "c26-1-error";
  const idStab = "c26-1-stability";
  const idSpec = "c26-1-spectrum";

  return (
    <motion.div {...fadeUp}>
      <Card className={`${theme.panel} p-4 rounded-2xl`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2"><Zap className="w-5 h-5" /> 26.1 Stiffness — test problems & stability</CardTitle>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 md:grid-cols-3 gap-3">
            <div className="md:col-span-2">
              <div className="grid grid-cols-1 sm:grid-cols-2 gap-3">
                <SmallLabel label="λ (stiffness parameter)">
                  <Input className="text-white" value={lambdaNum} onChange={(e) => setLambda(Number(e.target.value))} />
                </SmallLabel>
                <SmallLabel label="Explicit Euler step h">
                  <Input className="text-white" value={hExpNum} onChange={(e) => setHExp(Number(e.target.value))} />
                </SmallLabel>
              </div>

              <div className="grid grid-cols-1 md:grid-cols-2 gap-3 mt-3">
                <ResponsivePlot
                  id={idSol}
                  data={[
                    { x: tExact, y: yExact, type: "scatter", mode: "lines", name: "Exact", line: { color: theme.accent2 } },
                    { x: eulerExp.t, y: eulerExp.y.map((v) => (Array.isArray(v) ? v[0] : v)), type: "scatter", mode: "lines+markers", name: `Explicit Euler h=${fmt(hExpNum,6)}`, line: { color: theme.warn } },
                    { x: implicitEuler.t, y: implicitEuler.y, type: "scatter", mode: "lines+markers", name: `Implicit Euler h=${fmt(hImpNum,6)}`, line: { color: theme.accent } },
                  ]}
                  title="Test equation: explicit vs implicit"
                />

                <ResponsivePlot
                  id={idErr}
                  data={[
                    { x: conv.map((d) => d.h), y: conv.map((d) => d.errE), mode: "lines+markers", name: "Explicit Euler", type: "scatter", line: { color: theme.warn }, marker: { size: 6 } },
                    { x: conv.map((d) => d.h), y: conv.map((d) => d.errB), mode: "lines+markers", name: "Implicit Euler", type: "scatter", line: { color: theme.accent2 }, marker: { size: 6 } },
                  ]}
                  layout={{ xaxis: { title: "h", type: "log" }, yaxis: { title: "abs error", type: "log" } }}
                  title="Error at t=end vs h (log-log)"
                />
              </div>
            </div>

            <div>
              <div className="space-y-3">
                <div className="text-sm text-emerald-500 font-mono whitespace-pre-wrap break-words">
                  {"Test:  y' = λ y,   y(0) = 1"}
                </div>
                <div className="text-sm text-emerald-500 font-mono whitespace-pre-wrap break-words">
                  {"Explicit Euler:  y_{n+1} = (1 + hλ) y_n"}
                </div>
                <div className="text-sm text-emerald-500 font-mono whitespace-pre-wrap break-words">
                  {"Backward Euler:  y_{n+1} = y_n / (1 - hλ)"}
                </div>


                <Note>
                  Stability requirement for explicit Euler on linear test problem: <span className="font-mono">|1 + h λ| &lt; 1</span>. This leads to restrictive h when λ has large negative real part.
                </Note>

                <div className="mt-2 grid gap-2">
                  <div className="text-xs text-zinc-400">Quick numbers</div>
                  <div className="rounded-md text-gray-200 bg-zinc-900/50 p-2 border border-zinc-700 text-sm">
                    <div>Exact y(1) = <span className="font-mono">{exact(1).toExponential(6)}</span></div>
                    <div>Exp Euler final (h={fmt(hExpNum,6)}): <span className="font-mono">{fmt((Array.isArray(eulerExp.y[eulerExp.y.length-1]) ? eulerExp.y[eulerExp.y.length-1][0] : eulerExp.y[eulerExp.y.length-1]),6)}</span></div>
                    <div>Implicit Euler final (h={fmt(hImpNum,6)}): <span className="font-mono">{fmt(implicitEuler.y[implicitEuler.y.length-1],6)}</span></div>
                  </div>
                </div>

                <div className="mt-2 flex gap-2 flex-wrap">
                  <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" variant="ghost" onClick={() => exportPlotPNG(idSol, "26_1_solution.png")}><ImageIcon className="w-4 h-4" /> Export solution</Button>
                  <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" variant="ghost" onClick={() => exportPlotPNG(idErr, "26_1_error.png")}><Download className="w-4 h-4" /> Export error</Button>
                </div>
              </div>
            </div>
          </div>

          {/* stability region and spectrum */}
          <div className="grid grid-cols-1 md:grid-cols-3 gap-3">
            <ResponsivePlot
              id={idStab}
              data={[
                { x: circleX, y: circleY, mode: "lines", type: "scatter", fill: "toself", name: "Stability disk |1+z|<1", line: { color: theme.accent2 } },
                { x: markerZ.map((m) => m.z), y: markerZ.map(() => 0), mode: "markers+text", name: "λh (real)", marker: { color: theme.warn, size: 8 }, text: markerZ.map((m) => `h=${m.h.toFixed(3)}`), textposition: "top center" },
              ]}
              layout={{ xaxis: { title: "Re(z)" }, yaxis: { title: "Im(z)" } }}
              title="Explicit Euler stability region (z-plane)"
            />

            <ResponsivePlot
              id={idSpec}
              data={[
                { x: lambdaRange, y: lambdaRange.map((l) => 0), mode: "lines", type: "scatter", name: "λ (real axis)", line: { color: theme.accent } },
                { x: markerZ.map((m) => m.z), y: markerZ.map(() => 0), mode: "markers+text", marker: { size: 8, color: theme.warn }, text: markerZ.map((m) => `h=${m.h.toFixed(3)}`), textposition: "top center", name: "λh points" },
              ]}
              layout={{ xaxis: { title: "λ (real)" }, yaxis: { visible: false } }}
              title="Spectrum view (real axis)"
            />

            <ResponsivePlot
              id={"c26-1-margin"}
              data={[
                { x: linspace(0.001, Math.abs(2 / lambdaNum) * 2 + 1e-6, 120), y: linspace(0.001, Math.abs(2 / lambdaNum) * 2 + 1e-6, 120).map((h) => Math.abs(1 + lambdaNum * h)), mode: "lines", type: "scatter", name: "|1 + λh|" , line: { color: theme.accent2 } },
                { x: [hExpNum, hImpNum], y: [Math.abs(1 + lambdaNum * hExpNum), Math.abs(1 + lambdaNum * hImpNum)], mode: "markers+text", marker: { size: 8, color: theme.warn }, text: ["h_exp", "h_imp"], textposition: "top center" },
              ]}
              layout={{ xaxis: { title: "h" }, yaxis: { title: "|1+λh|" } }}
              title="Stability margin |1 + λh| vs h"
            />
          </div>

          <CollapsibleExercise
          title="Exercise 26.1.1"
          body={`Find the maximum stable step size h for the explicit Euler method applied to y' = λy. The stability condition is |1 + hλ| < 1.`}
          solution={`The inequality -1 < 1 + hλ < 1 gives -2 < hλ < 0. For λ < 0, this means the stable step size is 0 < h < -2/λ.`}
        />

        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// Section 26.2 — Multistep Methods (expanded with extra plots)
// ======================================================================
function Section262() {
  // sample non-stiff test ODE: y' = y - t^2 + 1, y(0)=0.5
  // exact solution: y = (t + 1)^2 - 0.5 e^t  (for verification)
  function f(t, y) {
    return y - t * t + 1;
  }
  function exact(t) {
    return (t + 1) * (t + 1) - 0.5 * Math.exp(t);
  }

  const [h, setH] = useState(0.1);
  const [pcIters, setPcIters] = useState(3);
  const t0 = 0;
  const tEnd = 2;
  const hNum = Math.max(1e-8, Number(h) || 0.1);
  const nSteps = Math.max(2, Math.floor((tEnd - t0) / hNum));

  const rk4 = useMemo(() => rk4Solver(f, t0, 0.5, hNum, nSteps), [hNum, nSteps]);
  const ab2 = useMemo(() => adamsBashforth2(f, t0, 0.5, hNum, nSteps), [hNum, nSteps]);
  const ab3 = useMemo(() => adamsBashforth3(f, t0, 0.5, hNum, nSteps), [hNum, nSteps]);
  const am2 = useMemo(() => adamsMoulton2_pc(f, t0, 0.5, hNum, nSteps, pcIters), [hNum, nSteps, pcIters]);

  const tExact = linspace(t0, tEnd, 400);
  const yExact = tExact.map((tt) => exact(tt));

  const xRK = rk4.t;
  const yRK = rk4.y.map((v) => (Array.isArray(v) ? v[0] : v));
  const xAB2 = ab2.t;
  const yAB2 = ab2.y.map((v) => (Array.isArray(v) ? v[0] : v));
  const xAB3 = ab3.t;
  const yAB3 = ab3.y.map((v) => (Array.isArray(v) ? v[0] : v));
  const xAM2 = am2.t;
  const yAM2 = am2.y.map((v) => (Array.isArray(v) ? v[0] : v));

  const errRK = Math.abs(yRK[yRK.length - 1] - exact(tEnd));
  const errAB2 = Math.abs(yAB2[yAB2.length - 1] - exact(tEnd));
  const errAB3 = Math.abs(yAB3[yAB3.length - 1] - exact(tEnd));
  const errAM2 = Math.abs(yAM2[yAM2.length - 1] - exact(tEnd));

  // convergence study across hs
  const hs = [0.5, 0.25, 0.125, 0.0625, 0.03125];
  const conv = hs.map((hh) => {
    const N = Math.max(2, Math.floor((tEnd - t0) / hh));
    const r = rk4Solver(f, t0, 0.5, hh, N).y;
    const a2 = adamsBashforth2(f, t0, 0.5, hh, N).y;
    const a3 = adamsBashforth3(f, t0, 0.5, hh, N).y;
    const am = adamsMoulton2_pc(f, t0, 0.5, hh, N, pcIters).y;
    const rlast = r[r.length - 1];
    const a2last = a2[a2.length - 1];
    const a3last = a3[a3.length - 1];
    const amlast = am[am.length - 1];
    const ref = exact(tEnd);
    return {
      h: hh,
      rkErr: Math.abs((Array.isArray(rlast) ? rlast[0] : rlast) - ref),
      ab2Err: Math.abs((Array.isArray(a2last) ? a2last[0] : a2last) - ref),
      ab3Err: Math.abs((Array.isArray(a3last) ? a3last[0] : a3last) - ref),
      am2Err: Math.abs((Array.isArray(amlast) ? amlast[0] : amlast) - ref),
    };
  });

  const idSol = "c26-2-solution";
  const idErrPoint = "c26-2-pointwise";
  const idConv = "c26-2-conv";
  const idHist = "c26-2-hist";

  // histogram of local errors for each method (approximate)
  const localErrors = {
    rk4: xRK.map((t, i) => Math.abs(yRK[i] - exact(t))),
    ab2: xAB2.map((t, i) => Math.abs(yAB2[i] - exact(t))),
    ab3: xAB3.map((t, i) => Math.abs(yAB3[i] - exact(t))),
    am2: xAM2.map((t, i) => Math.abs(yAM2[i] - exact(t))),
  };

  return (
    <motion.div {...fadeUp}>
      <Card className={`${theme.panel} p-4 rounded-2xl`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2"><Layers className="w-5 h-5" /> 26.2 Multistep Methods — AB / AM</CardTitle>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 md:grid-cols-3 gap-3">
            <div className="md:col-span-2">
              <div className="grid grid-cols-1 sm:grid-cols-2 gap-3">
                <SmallLabel label="Step size h">
                  <Input className="text-white" value={hNum} onChange={(e) => setH(Number(e.target.value))} />
                </SmallLabel>
                <SmallLabel label="PC iterations (AM2)">
                  <Input className="text-white" value={pcIters} onChange={(e) => setPcIters(Number(e.target.value))} />
                </SmallLabel>
              </div>

              <div className="grid grid-cols-1 md:grid-cols-2 gap-3 mt-3">
                <ResponsivePlot
                  id={idSol}
                  data={[
                    { x: tExact, y: yExact, mode: "lines", type: "scatter", name: "Exact", line: { color: theme.accent2 } },
                    { x: xRK, y: yRK, mode: "lines+markers", name: `RK4 h=${hNum}`, line: { color: theme.accent } },
                    { x: xAB2, y: yAB2, mode: "lines+markers", name: `AB2 h=${hNum}`, line: { color: theme.warn } },
                    { x: xAB3, y: yAB3, mode: "lines+markers", name: `AB3 h=${hNum}`, line: { color: "#f59e0b" } },
                    { x: xAM2, y: yAM2, mode: "lines+markers", name: `AM2 (PC) h=${hNum}`, line: { color: theme.accent2 } },
                  ]}
                  title="Solutions: multistep vs RK4"
                />

                <ResponsivePlot
                  id={idErrPoint}
                  data={[
                    { x: xRK, y: localErrors.rk4, mode: "lines", name: "RK4 local err", line: { color: theme.accent } },
                    { x: xAB2, y: localErrors.ab2, mode: "lines", name: "AB2 local err", line: { color: theme.warn } },
                    { x: xAB3, y: localErrors.ab3, mode: "lines", name: "AB3 local err", line: { color: "#f59e0b" } },
                    { x: xAM2, y: localErrors.am2, mode: "lines", name: "AM2 local err", line: { color: theme.accent2 } },
                  ]}
                  title="Pointwise absolute errors"
                />
              </div>
            </div>

            <div>
              <div className="space-y-3">
                <div className="text-sm text-emerald-500 font-mono whitespace-pre-wrap break-words">
                  {"AB2: y_{n+1} = y_n + (h/2)(3f_n - f_{n-1})"}
                </div>
                <div className="text-sm text-emerald-500 font-mono whitespace-pre-wrap break-words">
                  {"AM2 (PC): y_{n+1} = y_n + (h/2)(f_{n+1} + f_n)"}
                </div>
     <Note>Use RK4 for start-up. Predictor-corrector reduces cost of implicit methods by avoiding full Newton solves for many problems.</Note>

                <div className="mt-2 text-xs text-zinc-400">
                  <div>Errors at t = {tEnd}:</div>
                  <div className="text-sm font-mono mt-1">
                    RK4: {fmt(errRK, 8)} • AB2: {fmt(errAB2, 8)} • AB3: {fmt(errAB3, 8)} • AM2: {fmt(errAM2, 8)}
                  </div>
                </div>

                <div className="mt-3 flex gap-2">
                  <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" variant="ghost" onClick={() => exportPlotPNG(idSol, "26_2_solution.png")}><ImageIcon className="w-4 h-4" /> Export</Button>
                  <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" variant="ghost" onClick={() => exportPlotPNG(idErrPoint, "26_2_pointwise.png")}><Download className="w-4 h-4" /> Export</Button>
                </div>
              </div>
            </div>
          </div>

          {/* convergence and histograms */}
          <div className="grid grid-cols-1 md:grid-cols-3 gap-3">
            <ResponsivePlot
              id={idConv}
              data={[
                { x: conv.map((d) => d.h), y: conv.map((d) => d.rkErr), mode: "lines+markers", name: "RK4", line: { color: theme.accent } },
                { x: conv.map((d) => d.h), y: conv.map((d) => d.ab2Err), mode: "lines+markers", name: "AB2", line: { color: theme.warn } },
                { x: conv.map((d) => d.h), y: conv.map((d) => d.ab3Err), mode: "lines+markers", name: "AB3", line: { color: "#f59e0b" } },
                { x: conv.map((d) => d.h), y: conv.map((d) => d.am2Err), mode: "lines+markers", name: "AM2", line: { color: theme.accent2 } },
              ]}
              layout={{ xaxis: { type: "log", title: "h" }, yaxis: { type: "log", title: "abs error" } }}
              title="Convergence at final time (log-log)"
            />

            <ResponsivePlot
              id={idHist}
              data={[
                { x: localErrors.rk4, type: "histogram", name: "RK4 local error", opacity: 0.7 },
                { x: localErrors.ab2, type: "histogram", name: "AB2 local error", opacity: 0.7 },
                { x: localErrors.ab3, type: "histogram", name: "AB3 local error", opacity: 0.7 },
                { x: localErrors.am2, type: "histogram", name: "AM2 local error", opacity: 0.7 },
              ]}
              layout={{ barmode: "overlay", xaxis: { title: "local abs error" } }}
              title="Distribution of local errors (histogram)"
              style={{ height: "clamp(240px, 40vh, 520px)" }}
            />

            <div className="rounded-xl border border-zinc-700 bg-zinc-900/40 p-3">
              <div className="text-xs text-zinc-400">Cost comparison</div>
              <div className="mt-2 text-sm text-zinc-100">
                <div>RK4: 4 f-evals per step</div>
                <div>AB2/AB3: 1 f-eval per step (after start-up)</div>
                <div>AM2 (PC): predictor (AB2) + k corrector iterations (each requires f-eval at predicted time) — trade-off between cost and accuracy</div>
              </div>

              <div className="mt-3 text-xs text-zinc-400">Export</div>
              <div className="mt-2 flex gap-2">
                <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" variant="ghost" onClick={() => downloadCSV("26_2_multistep_summary.csv", [["method", "error"], ["RK4", errRK], ["AB2", errAB2], ["AB3", errAB3], ["AM2", errAM2]])}><Download className="w-4 h-4" /> CSV</Button>
              </div>
            </div>
          </div>

          <CollapsibleExercise
            title="Exercise 26.2.1"
            body={`Implement the Adams–Bashforth 3-step (AB3) method and the classical Runge–Kutta method (RK4) for y' = y - t^2 + 1 on the interval [0, 2] with step size h = 0.05. Compare the absolute error at t = 2 and count the total number of function evaluations.`}
            solution={`RK4 requires 4 function evaluations per step, while AB3 requires only 1 after startup. For moderate step sizes, RK4 usually produces smaller error but at a higher cost. Compare cost-normalized error (error per function evaluation) to determine efficiency.`}
          />

        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// Chapter notes and problems panels
// ======================================================================
function ChapterDocs() {
  const [open, setOpen] = useState("26.1 Stiffness");
  const docs = {
    "26.1 Stiffness": [
      `Stiffness refers to systems where there are rapid transients together with slow dynamics. For the linear test problem y' = λy, stiffness occurs when Re(λ) is very negative and large step sizes make explicit methods unstable. Implicit methods such as Backward Euler or BDF have larger stability regions (often A-stable) and are preferred for stiff problems.`,
    ],
    "26.2 Multistep Methods": [
      `Adams–Bashforth methods are explicit multistep integrators that reuse previous function evaluations. Adams–Moulton methods are implicit and often used as correctors in predictor–corrector schemes. Multistep methods require a start-up phase (for example, RK4 or repeated Euler) to generate the first few values.`,
    ],
  };

  return (
    <div className="mt-6 bg-zinc-900 border border-zinc-700 rounded-xl p-3 overflow-hidden">
      {/* Header */}
      <div className="flex items-center justify-between mb-3">
        <div className="flex items-center gap-3">
          <BookOpen className="w-5 h-5 text-cyan-400" />
          <div>
            <div className="text-zinc-100 font-semibold">Chapter 26 — Notes</div>
            <div className="text-zinc-400 text-xs">
              Stiffness, stability, multistep methods
            </div>
          </div>
        </div>
        <div className="text-zinc-400 text-xs">Interactive examples above</div>
      </div>

      {/* Content */}
      <div className="grid grid-cols-1 md:grid-cols-4">
        {/* Sidebar */}
        <div className="col-span-1 border-r border-zinc-800">
          {Object.keys(docs).map((k) => (
            <button
              key={k}
              onClick={() => setOpen(k)}
              className={`w-full text-left p-3 hover:bg-zinc-800/30 ${
                open === k ? "bg-zinc-800/20" : ""
              }`}
            >
              <div className="flex items-center gap-2">
                <List className="w-4 h-4 text-zinc-300" />
                <div className="text-zinc-100 text-sm">{k}</div>
              </div>
            </button>
          ))}
        </div>

        {/* Main content */}
        <div
          className="col-span-3 p-4 overflow-auto"
          style={{ maxHeight: 460 }}
        >
          <h3 className="text-lg text-zinc-100 font-semibold mb-2">{open}</h3>
          <div className="text-zinc-300 text-sm space-y-3">
            {docs[open].map((p, i) => (
              <p key={i} className="whitespace-pre-wrap break-words">
                {p}
              </p>
            ))}
            <div className="text-xs text-zinc-400 mt-2">
              References: Hairer & Wanner — *Solving Ordinary Differential
              Equations II*; Butcher — *Numerical Methods for ODEs*.
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}


function ProblemsPanel() {
  return (
    <div className={`${theme.panel} p-4 rounded-2xl`}>
      {/* Header */}
      <div className="flex items-center justify-between mb-3">
        <div className="text-xl text-cyan-300 font-semibold flex items-center gap-2">
          <FileText className="w-5 h-5" /> Problems
        </div>
        <Button
          variant="ghost"
          size="sm"
          onClick={() => window.scrollTo({ top: 0, behavior: "smooth" })}
        >
          <RefreshCcw className="w-4 h-4" /> Top
        </Button>
      </div>

      {/* Exercises Grid */}
      <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
        <CollapsibleExercise
          title="Problem 25.1 — Euler and RK4 comparison"
          body={`For dy/dx = -2y + x, with y(0) = 1, compute y(2) using Euler with h = 0.2 and RK4 with h = 0.2. Estimate absolute errors using a dense RK4 reference. Which method is more accurate per step?`}
          solution={`RK4 is significantly more accurate per step; the error scaling demonstrates Euler is O(h), while RK4 is O(h^4).`}
        />

        <CollapsibleExercise
          title="Problem 25.2 — Midpoint vs Heun"
          body={`Implement both midpoint and Heun methods on the sample ODE and measure the global error versus step size h. Compare the computational cost (function evaluations) to reach a target error of 1e-4.`}
          solution={`Both Midpoint and Heun are second-order methods. The tradeoff depends on the number of function evaluations per step (both use 2 evals).`}
        />

        <CollapsibleExercise
          title="Problem 25.3 — Long-term stability"
          body={`For the harmonic oscillator, compare Euler, RK4, and symplectic Verlet over a long time horizon T = 1000. Plot energy drift. Which method best conserves energy?`}
          solution={`Symplectic Verlet conserves energy with bounded error. RK4 shows slow drift, while Euler quickly diverges.`}
        />

        <CollapsibleExercise
          title="Problem 25.4 — Adaptive RK usage"
          body={`Use an adaptive Runge–Kutta method (RKF45) to integrate a stiff test equation. Compare the steps taken versus a fixed-step RK4 achieving similar error. Which approach is more efficient?`}
          solution={`Adaptive RK adjusts effort where needed, often outperforming fixed-step methods when the solution has localized rapid changes.`}
        />
      </div>
    </div>
  );
}

// ======================================================================
// Page assembly and export
// ======================================================================

export default function Chapter26() {
  // focus handling for accessibility
  const topRef = useRef(null);
  useEffect(() => { topRef.current?.focus?.(); }, []);
  return (
    <div className={`p-4 md:p-6 lg:p-8 ${theme.bg} min-h-screen`}>
      <div className="max-w-screen-xl mx-auto space-y-6">
        <motion.div initial={{ opacity: 0, y: -8 }} animate={{ opacity: 1, y: 0 }}>
          <h1 ref={topRef} tabIndex={-1} className="text-2xl md:text-3xl lg:text-4xl font-bold" style={{ color: theme.accent2 }}>Stiffness & Multistep Methods</h1>
          <div className="text-zinc-400 mt-1">Stiffness • Multistep Methods • Problems — interactive demos and plots below.</div>
        </motion.div>

        <Section261 />
        <Section262 />

        <ChapterDocs />
        <ProblemsPanel />
        <BottomBar/>
      </div>
    </div>
  );
}
