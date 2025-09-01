// src/pages/Chapter14.jsx
// ======================================================================
// Chapter 14 — Multidimensional Unconstrained Optimization
// Sections: 14.1 Direct Methods, 14.2 Gradient Methods, Problems
// Dark-gray theme, emerald/amber accents, stacked section design
// ======================================================================

import React, { useMemo, useState } from "react";
import Plot from "react-plotly.js";
import { motion } from "framer-motion";
import ReactMarkdown from "react-markdown";
import remarkMath from "remark-math";
import rehypeKatex from "rehype-katex";
import "katex/dist/katex.min.css";

import { Grid, ArrowRightCircle, Zap, BookOpen, Code, ChevronDown, ChevronRight } from "lucide-react";
import BottomBar from "../components/BottomBar";

// UI fallbacks if your project doesn't have custom components
let Input, Button;
try {
  // eslint-disable-next-line
  Input = require("@/components/ui/input").Input;
  // eslint-disable-next-line
  Button = require("@/components/ui/button").Button;
} catch (e) {
  Input = (props) => <input {...props} />;
  Button = (props) => <button {...props} />;
}

// Theme and animation
const THEME = {
  panelBg: "#0b0f20",
  text: "#e6eef6",
  muted: "#9ca3af",
  accent: "#16a34a",
  accent2: "#d97706",
};

const fadeUp = { initial: { opacity: 0, y: 8 }, animate: { opacity: 1, y: 0 }, transition: { duration: 0.28 } };

// ------------------------ Utility: safe function parser ------------------------
function makeFun(expr) {
  try {
    // Provide Math.* helpers
    const wrapped = `
      const { sin, cos, tan, asin, acos, atan, exp, log, sqrt, pow, abs, PI, E } = Math;
      return (x) => { const [x0,x1,x2,x3] = Array.isArray(x) ? x : [x]; return (${expr}); }
    `;
    return new Function(wrapped)();
  } catch (e) {
    return null;
  }
}

// gradient approximations (central differences)
function gradApprox(f, x, h = 1e-6) {
  const n = x.length;
  const g = Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    const xp = x.slice(); xp[i] += h;
    const xm = x.slice(); xm[i] -= h;
    g[i] = (f(xp) - f(xm)) / (2 * h);
  }
  return g;
}
function normVec(v) { return Math.sqrt(v.reduce((s, a) => s + a * a, 0)); }
function addVec(a, b) { return a.map((x, i) => x + b[i]); }
function subVec(a, b) { return a.map((x, i) => x - b[i]); }
function scaleVec(a, s) { return a.map((x) => x * s); }
function dot(a, b) { return a.reduce((s, x, i) => s + x * b[i], 0); }

// ------------------------ Small UI helpers ------------------------
function SectionHeader({ Icon, title, subtitle }) {
  return (
    <div className="flex items-center gap-3 mb-3">
      <div style={{ width: 44, height: 44, display: "flex", alignItems: "center", justifyContent: "center", background: "#071018", borderRadius: 8 }}>
        <Icon className="w-5 h-5" />
      </div>
      <div>
        <div style={{ color: THEME.text, fontSize: 18, fontWeight: 700 }}>{title}</div>
        <div style={{ color: THEME.muted, fontSize: 13 }}>{subtitle}</div>
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

// ------------------------ Documentation strings ------------------------
const docs = {
  "14.1 Direct Methods": [
    `Direct (derivative-free) methods explore the domain without explicit gradient information. Examples: grid search, coordinate search, Nelder–Mead simplex.`,
    `Nelder–Mead (simplex) adapts a simplex through reflection, expansion, contraction and shrink steps to seek minima without derivatives.`,
  ],
  "14.2 Gradient Methods": [
    `Gradient-based methods (steepest descent, conjugate gradient, quasi-Newton) leverage gradient information for faster convergence.`,
    `Backtracking line-search (Armijo) is commonly used to choose a step size satisfying sufficient decrease.`,
  ],
};

// ------------------------ SECTION: 14.1 Direct Methods ------------------------
function Section141() {
  // 2D examples for visualization & methods
  const [expr, setExpr] = useState(" (x0-1)*(x0-1) + 2*(x1+0.5)*(x1+0.5) + 0.3*Math.cos(3*x0+2*x1) ");
  const [gridRange, setGridRange] = useState("-3,3");
  const [gridPts, setGridPts] = useState(41);

  const f = useMemo(() => makeFun(expr), [expr]);

  // grid-search (coarse)
  function gridSearch(f, xlo, xhi, ylo, yhi, nx) {
    const xs = [];
    let best = { x: [xlo, ylo], fx: Infinity };
    for (let i = 0; i < nx; i++) {
      const xi = xlo + (xhi - xlo) * (i / (nx - 1));
      for (let j = 0; j < nx; j++) {
        const yj = ylo + (yhi - ylo) * (j / (nx - 1));
        const fx = f([xi, yj]);
        xs.push({ x: xi, y: yj, fx });
        if (fx < best.fx) best = { x: [xi, yj], fx };
      }
    }
    return { samples: xs, best };
  }

  // simple coordinate search (pattern search)
  function coordinateSearch(f, x0, step0 = 0.5, tol = 1e-6, maxIter = 200) {
    let x = x0.slice();
    let step = step0;
    const trace = [{ k: 0, x: x.slice(), fx: f(x) }];
    for (let k = 0; k < maxIter; k++) {
      let improved = false;
      for (let d = 0; d < x.length; d++) {
        const xp = x.slice(); xp[d] += step;
        const xm = x.slice(); xm[d] -= step;
        const fp = f(xp), fm = f(xm), f0 = f(x);
        if (fp < f0) { x = xp; improved = true; }
        else if (fm < f0) { x = xm; improved = true; }
      }
      trace.push({ k: k + 1, x: x.slice(), fx: f(x) });
      if (!improved) step *= 0.5;
      if (step < tol) break;
    }
    return { x, fx: f(x), trace };
  }

  // Nelder-Mead (simple educational implementation)
  function nelderMead(f, x0, step = 0.5, tol = 1e-6, maxIter = 200) {
    const n = x0.length;
    // initialize simplex: x0 and x0 + step * ei
    let simplex = [x0.slice()];
    for (let i = 0; i < n; i++) {
      const v = x0.slice(); v[i] += step; simplex.push(v);
    }
    let fs = simplex.map((s) => f(s));
    const trace = [{ iter: 0, simplex: simplex.map(s => s.slice()), fs: fs.slice() }];
    for (let iter = 0; iter < maxIter; iter++) {
      // order
      const idx = fs.map((val, i) => ({ val, i })).sort((a, b) => a.val - b.val).map(o => o.i);
      const best = simplex[idx[0]];
      const worst = simplex[idx[n]];
      const secondWorst = simplex[idx[n - 1]];
      const centroid = Array(n).fill(0);
      for (let k = 0; k <= n - 1; k++) {
        const s = simplex[idx[k]];
        for (let j = 0; j < n; j++) centroid[j] += s[j] / n;
      }
      // reflection
      const alpha = 1.0, gamma = 2.0, rho = 0.5, sigma = 0.5;
      const xr = centroid.map((c, i) => c + alpha * (c - worst[i]));
      const fxr = f(xr);
      if (fxr < fs[idx[0]]) {
        // expansion
        const xe = centroid.map((c, i) => c + gamma * (xr[i] - c));
        const fxe = f(xe);
        if (fxe < fxr) {
          simplex[idx[n]] = xe; fs[idx[n]] = fxe;
        } else {
          simplex[idx[n]] = xr; fs[idx[n]] = fxr;
        }
      } else if (fxr < fs[idx[n - 1]]) {
        simplex[idx[n]] = xr; fs[idx[n]] = fxr;
      } else {
        // contraction
        let xc;
        if (fxr < fs[idx[n]]) {
          // outside contraction
          xc = centroid.map((c, i) => c + rho * (xr[i] - c));
        } else {
          // inside contraction
          xc = centroid.map((c, i) => c + rho * (worst[i] - c));
        }
        const fxc = f(xc);
        if (fxc < fs[idx[n]]) {
          simplex[idx[n]] = xc; fs[idx[n]] = fxc;
        } else {
          // shrink
          for (let k = 1; k < simplex.length; k++) {
            const si0 = simplex[idx[0]];
            simplex[idx[k]] = simplex[idx[k]].map((v, i) => si0[i] + sigma * (v - si0[i]));
            fs[idx[k]] = f(simplex[idx[k]]);
          }
        }
      }
      trace.push({ iter: iter + 1, simplex: simplex.map(s => s.slice()), fs: fs.slice() });
      const fvals = fs.slice().sort((a, b) => a - b);
      if (Math.abs(fvals[0] - fvals[fvals.length - 1]) < tol) break;
    }
    // return best
    const bestIdx = fs.map((v, i) => ({ v, i })).sort((a, b) => a.v - b.v)[0].i;
    return { x: simplex[bestIdx], fx: fs[bestIdx], trace };
  }

  const parsedRange = useMemo(() => {
    const parts = gridRange.split(",").map((s) => Number(s.trim()));
    if (parts.length >= 2) {
      return { lo: parts[0], hi: parts[1] };
    }
    return { lo: -3, hi: 3 };
  }, [gridRange]);

  const gridResult = useMemo(() => {
    if (!f) return null;
    try {
      return gridSearch(f, parsedRange.lo, parsedRange.hi, parsedRange.lo, parsedRange.hi, Number(gridPts));
    } catch (e) {
      return null;
    }
  }, [f, parsedRange, gridPts]);

  const coordResult = useMemo(() => {
    if (!f) return null;
    try {
      return coordinateSearch(f, [0, 0], 0.8, 1e-6, 200);
    } catch (e) {
      return null;
    }
  }, [f]);

  const nmResult = useMemo(() => {
    if (!f) return null;
    try {
      return nelderMead(f, [0, 0], 0.8, 1e-6, 300);
    } catch (e) {
      return null;
    }
  }, [f]);

  // prepare plot grid samples for heatmap/contour
  const surface = useMemo(() => {
    if (!f) return null;
    const lo = parsedRange.lo, hi = parsedRange.hi, N = 80;
    const xs = Array.from({ length: N }, (_, i) => lo + (hi - lo) * (i / (N - 1)));
    const ys = xs.slice();
    const z = Array.from({ length: N }, () => Array(N).fill(0));
    for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) z[i][j] = f([xs[j], ys[i]]);
    return { xs, ys, z };
  }, [f, parsedRange]);

  return (
    <section className="rounded-lg p-4" style={{ background: THEME.panelBg }}>
      <motion.div {...fadeUp}>
        <SectionHeader Icon={Grid} title="14.1 Direct Methods" subtitle="Grid search, coordinate search, Nelder–Mead simplex (derivative-free)" />

        <div className="rounded bg-zinc-800 border border-zinc-700 p-3 mb-3">
          <div className="grid grid-cols-1 md:grid-cols-3 gap-3">
            <div className="md:col-span-2">
              <label className="text-xs text-zinc-400">Multivariate function f([x0,x1])</label>
              <Input value={expr} onChange={(e) => setExpr(e.target.value)} className="w-full bg-zinc-900 p-2 rounded text-zinc-100" />
              <div className="text-xs text-zinc-400 mt-1">Use x0,x1 variables (JS Math allowed)</div>
            </div>
            <div>
              <label className="text-xs text-zinc-400">Plot range (lo,hi)</label>
              <Input value={gridRange} onChange={(e) => setGridRange(e.target.value)} className="w-full bg-zinc-900 p-2 rounded text-zinc-100" />
              <label className="text-xs text-zinc-400 mt-2">Grid pts (N x N)</label>
              <Input value={gridPts} onChange={(e) => setGridPts(e.target.value)} className="w-full bg-zinc-900 p-2 rounded text-zinc-100" />
            </div>
          </div>
        </div>

        <div className="grid grid-cols-1 lg:grid-cols-3 gap-3 mb-3">
          <div className="lg:col-span-2 rounded border border-zinc-700 bg-zinc-900 p-3">
            <div className="text-xs text-zinc-400 mb-2">Contour & samples</div>
            {surface ? (
              <Plot
                data={[
                  { x: surface.xs, y: surface.ys, z: surface.z, type: "contour", contours: { coloring: "heatmap" } },
                  ...(gridResult ? [{ x: gridResult.samples.map(s => s.x), y: gridResult.samples.map(s => s.y), mode: "markers", marker: { size: 4 }, name: "grid" }] : []),
                  ...(nmResult ? [{ x: nmResult.trace.map(t => t.simplex.map(s=>s[0])).flat(), y: nmResult.trace.map(t => t.simplex.map(s=>s[1])).flat(), mode: "markers", name: "nm simplex verts", marker: { size: 4 } }] : []),
                ]}
                layout={{ autosize: true, height: 420, paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)" }}
                useResizeHandler
                style={{ width: "100%" }}
              />
            ) : <div className="text-xs text-zinc-400">Invalid function</div>}
          </div>

          <div className="rounded border border-zinc-700 bg-zinc-800 p-3">
            <div className="text-xs text-zinc-400 mb-2">Direct method results</div>
            <div className="text-sm text-zinc-100">Grid best: {gridResult ? prettyPoint(gridResult.best.x) + " f=" + fmt(gridResult.best.fx) : "—"}</div>
            <div className="text-sm text-zinc-100 mt-1">Coord search: {coordResult ? prettyPoint(coordResult.x) + " f=" + fmt(coordResult.fx) : "—"}</div>
            <div className="text-sm text-zinc-100 mt-1">Nelder-Mead: {nmResult ? prettyPoint(nmResult.x) + " f=" + fmt(nmResult.fx) : "—"}</div>
          </div>
        </div>

        <div>
          <Exercise
            title="Exercise 14.1.1 — Nelder–Mead steps"
            body={`Run Nelder–Mead starting from different initial points. Observe how the simplex evolves (inspect trace). What behavior do you see near flat regions or ridges?`}
            solution={`Near flat regions Nelder–Mead can shrink the simplex and make slow progress; on ridges it may adapt shape but can stall unless contraction/expansion parameters are well chosen.`}
          />
        </div>
      </motion.div>
    </section>
  );
}

// ------------------------ SECTION: 14.2 Gradient Methods ------------------------
function Section142() {
  // example 2D quadratic + nonquadratic functions
  const [expr, setExpr] = useState(" 4*(x0-1)*(x0-1) + 2*(x1+0.5)*(x1+0.5) + 0.5*(x0-0.5)*(x1+0.7) ");
  const [x0Text, setX0Text] = useState(" -2, 2 "); // starting point
  const [alpha, setAlpha] = useState(1.0); // initial step size for backtracking
  const [c, setC] = useState(1e-4);
  const [rho, setRho] = useState(0.5);
  const [tol, setTol] = useState(1e-6);
  const [maxIter, setMaxIter] = useState(200);

  const f = useMemo(() => makeFun(expr), [expr]);

  function parsePoint(text) {
    const parts = text.split(/[\s,]+/).map(s => Number(s.trim())).filter(n => !Number.isNaN(n));
    if (parts.length < 2) return [0, 0];
    return [parts[0], parts[1]];
  }

  // backtracking line search (Armijo)
  function backtracking(f, x, p, g, alpha0 = 1.0, rho = 0.5, c = 1e-4) {
    let alpha = alpha0;
    const fx = f(x);
    const dotpg = dot(p, g);
    for (let iter = 0; iter < 50; iter++) {
      const xnew = addVec(x, scaleVec(p, alpha));
      const fxnew = f(xnew);
      if (fxnew <= fx + c * alpha * dotpg) return alpha;
      alpha *= rho;
    }
    return alpha;
  }

  function gradientDescent(f, x0, opts = {}) {
    const { alpha0 = 1.0, rho = 0.5, c = 1e-4, tol = 1e-6, maxIter = 200 } = opts;
    let x = x0.slice();
    const trace = [{ k: 0, x: x.slice(), fx: f(x), g: gradApprox(f, x) }];
    for (let k = 0; k < maxIter; k++) {
      const g = gradApprox(f, x);
      const p = scaleVec(g, -1);
      if (normVec(g) < tol) break;
      const alpha = backtracking(f, x, p, g, alpha0, rho, c);
      x = addVec(x, scaleVec(p, alpha));
      trace.push({ k: k + 1, x: x.slice(), fx: f(x), g });
    }
    return { x, fx: f(x), trace };
  }

  // Conjugate Gradient for quadratic Ax=b problems (we'll detect if user provided quadratic by analytic Hessian)
  // Here, for education, we'll implement CG for symmetric positive-definite A when f(x)=0.5 x^T A x - b^T x + const
  function conjugateGradient(A, b, x0 = null, tol = 1e-8, maxIter = 1000) {
    const n = b.length;
    let x = x0 ? x0.slice() : Array(n).fill(0);
    let r = subVec(b, matVec(A, x));
    let p = r.slice();
    let rsold = dot(r, r);
    const trace = [{ k: 0, x: x.slice(), resid: Math.sqrt(rsold) }];
    for (let k = 0; k < Math.min(maxIter, n); k++) {
      const Ap = matVec(A, p);
      const alpha = rsold / Math.max(dot(p, Ap), 1e-18);
      x = addVec(x, scaleVec(p, alpha));
      r = subVec(r, scaleVec(Ap, alpha));
      const rsnew = dot(r, r);
      trace.push({ k: k + 1, x: x.slice(), resid: Math.sqrt(rsnew) });
      if (Math.sqrt(rsnew) < tol) break;
      p = addVec(r, scaleVec(p, rsnew / rsold));
      rsold = rsnew;
    }
    return { x, trace };
  }

  // plotting surface for visualization (project along region around starting point)
  const point = useMemo(() => parsePoint(x0Text), [x0Text]);

  const gradResult = useMemo(() => {
    if (!f) return null;
    try {
      return gradientDescent(f, point.slice(), { alpha0: Number(alpha), rho: Number(rho), c: Number(c), tol: Number(tol), maxIter: Number(maxIter) });
    } catch (e) {
      return null;
    }
  }, [f, point, alpha, rho, c, tol, maxIter]);

  // For CG example, attempt to build A,b if expr looks quadratic: heuristically skip for now and provide demo A
  const demoA = useMemo(() => [[4, 1], [1, 3]], []);
  const demob = useMemo(() => [1, 2], []);
  const cgResult = useMemo(() => conjugateGradient(demoA, demob, [0, 0], 1e-8, 50), [demoA, demob]);

  // prepare 2D contour for plotting around start
  const surface = useMemo(() => {
    if (!f) return null;
    const left = point[0] - 3, right = point[0] + 3;
    const bottom = point[1] - 3, top = point[1] + 3;
    const N = 120;
    const xs = Array.from({ length: N }, (_, i) => left + (right - left) * (i / (N - 1)));
    const ys = Array.from({ length: N }, (_, j) => bottom + (top - bottom) * (j / (N - 1)));
    const z = Array.from({ length: N }, () => Array(N).fill(0));
    for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) z[i][j] = f([xs[j], ys[i]]);
    return { xs, ys, z };
  }, [f, point]);

  return (
    <section className="rounded-lg p-4" style={{ background: THEME.panelBg }}>
      <motion.div {...fadeUp}>
        <SectionHeader Icon={Zap} title="14.2 Gradient Methods" subtitle="Steepest descent, line search, conjugate gradient (quadratic demo)" />

        <div className="rounded bg-zinc-800 border border-zinc-700 p-3 mb-3">
          <div className="grid grid-cols-1 md:grid-cols-3 gap-3">
            <div className="md:col-span-2">
              <label className="text-xs text-zinc-400">Function f([x0,x1])</label>
              <Input value={expr} onChange={(e) => setExpr(e.target.value)} className="w-full bg-zinc-900 p-2 rounded text-zinc-100" />
              <div className="text-xs text-zinc-400 mt-1">Use x0,x1 variables (JS Math allowed)</div>
            </div>

            <div>
              <label className="text-xs text-zinc-400">Start point (x0,y0)</label>
              <Input value={x0Text} onChange={(e) => setX0Text(e.target.value)} className="w-full bg-zinc-900 p-2 rounded text-zinc-100" />
              <div className="grid grid-cols-2 gap-2 mt-2">
                <div>
                  <label className="text-xs text-zinc-400">alpha0</label>
                  <Input value={alpha} onChange={(e) => setAlpha(e.target.value)} className="bg-zinc-900 p-2 rounded text-zinc-100" />
                </div>
                <div>
                  <label className="text-xs text-zinc-400">rho</label>
                  <Input value={rho} onChange={(e) => setRho(e.target.value)} className="bg-zinc-900 p-2 rounded text-zinc-100" />
                </div>
              </div>
              <div className="grid grid-cols-2 gap-2 mt-2">
                <div>
                  <label className="text-xs text-zinc-400">c (Armijo)</label>
                  <Input value={c} onChange={(e) => setC(e.target.value)} className="bg-zinc-900 p-2 rounded text-zinc-100" />
                </div>
                <div>
                  <label className="text-xs text-zinc-400">tol</label>
                  <Input value={tol} onChange={(e) => setTol(e.target.value)} className="bg-zinc-900 p-2 rounded text-zinc-100" />
                </div>
              </div>
            </div>
          </div>
        </div>

        <div className="grid grid-cols-1 lg:grid-cols-3 gap-3 mb-3">
          <div className="lg:col-span-2 rounded border border-zinc-700 bg-zinc-900 p-3">
            <div className="text-xs text-zinc-400 mb-2">Contour & gradient path</div>
            {surface ? (
              <Plot
                data={[
                  { x: surface.xs, y: surface.ys, z: surface.z, type: "contour", contours: { coloring: "heatmap" } },
                  ...(gradResult ? [{ x: gradResult.trace.map(t => t.x[0]), y: gradResult.trace.map(t => t.x[1]), mode: "markers+lines", marker: { size: 6 }, name: "GD path" }] : []),
                  ...(cgResult ? [{ x: cgResult.trace.map(t => t.x[0]), y: cgResult.trace.map(t => t.x[1]), mode: "markers+lines", name: "CG path (demo A)", marker: { size: 6 } }] : []),
                ]}
                layout={{ autosize: true, height: 420, paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)" }}
                useResizeHandler
                style={{ width: "100%" }}
              />
            ) : <div className="text-xs text-zinc-400">Invalid function</div>}
          </div>

          <div className="rounded border border-zinc-700 bg-zinc-800 p-3">
            <div className="text-xs text-zinc-400 mb-2">Results</div>
            <div className="text-sm text-zinc-100">Gradient descent final: {gradResult ? prettyPoint(gradResult.x) + " f=" + fmt(gradResult.fx) : "—"}</div>
            <div className="text-sm text-zinc-100 mt-1">Iterations: {gradResult ? gradResult.trace.length : "—"}</div>
            <div className="text-sm text-zinc-100 mt-3">CG (demo A): final {cgResult ? prettyPoint(cgResult.x) : "—"}</div>
          </div>
        </div>

        <div>
          <Exercise
            title="Exercise 14.2.1 — Armijo tuning"
            body={`Experiment with different rho and c in the backtracking line search. How do these parameters affect step sizes and convergence?`}
            solution={`Smaller rho reduces step sizes faster when sufficient decrease fails. Larger c requires more decrease, potentially reducing step sizes; tune for balance.`}
          />
        </div>
      </motion.div>
    </section>
  );
}

// ------------------------ Problems section ------------------------
function SectionProblems() {
  const problems = [
    {
      title: "Problem 14.1",
      body: `Implement Nelder–Mead and compare its performance to coordinate search on the Rosenbrock function in 2D. Report iterations and final function value.`,
      solution: `Nelder–Mead often outperforms coordinate search on smooth nonconvex problems but can stall; Rosenbrock is challenging and benefits from gradient-based methods.`,
    },
    {
      title: "Problem 14.2",
      body: `Implement gradient descent with Armijo line search and test on a quadratic with ill-conditioned Hessian. Measure convergence vs condition number.`,
      solution: `Convergence rate is affected by condition number; steepest descent slows when condition number is large; CG or preconditioning helps.`,
    },
  ];
  return (
    <section className="rounded-lg p-4" style={{ background: THEME.panelBg }}>
      <motion.div {...fadeUp}>
        <SectionHeader Icon={Code} title="Problems (Chapter 14)" subtitle="Practice tasks" />
        <div className="grid grid-cols-1 gap-3">
          {problems.map((p, i) => <Exercise key={i} title={p.title} body={p.body} solution={p.solution} />)}
        </div>
      </motion.div>
    </section>
  );
}

// ------------------------ Utilities ------------------------
function fmt(x) {
  if (Array.isArray(x)) return `[${x.map(v => Number(v).toPrecision(6)).join(", ")}]`;
  return Number(x).toPrecision(6);
}
function prettyPoint(x) {
  return `(${x.map(v => Number(v).toFixed(4)).join(", ")})`;
}
// small mat helpers for CG demo
function matVec(A, x) { return A.map(row => row.reduce((s, v, i) => s + v * x[i], 0)); }

// ------------------------ Page assembly ------------------------
export default function Chapter14() {
  return (
    <div className="min-h-screen p-6 bg-zinc-950 text-zinc-100">
      <header className="mb-6">
        <div className="flex items-center justify-between">
          <div>
            <h1 className="text-3xl font-bold" style={{ color: THEME.accent }}>Multidimensional Unconstrained Optimization</h1>
            <div className="text-sm text-zinc-400 mt-1">Direct (derivative-free) and gradient-based methods with visualization</div>
          </div>
        </div>
      </header>

      <main className="space-y-6">
        <Section141 />
        <Section142 />
        <SectionProblems />
        <BottomBar/>
      </main>

      <footer className="mt-8 text-xs text-zinc-500">Reminder: for production use highly optimized libraries (e.g., NLopt, SciPy's optimize). This page is for teaching and visualization.</footer>
    </div>
  );
}
