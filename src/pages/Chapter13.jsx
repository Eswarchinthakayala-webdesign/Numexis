// src/pages/Chapter13.jsx
// ======================================================================
// Chapter 13 — One-Dimensional Unconstrained Optimization
// Sections: 13.1 Golden-Section Search, 13.2 Parabolic Interpolation,
// 13.3 Newton's Method (optim.), 13.4 Brent-like hybrid, Problems
// Dark-gray theme, emerald/amber accents, stacked layout
// ======================================================================

import React, { useMemo, useState, useEffect, useRef } from "react";
import Plot from "react-plotly.js";
import { motion } from "framer-motion";
import ReactMarkdown from "react-markdown";
import remarkMath from "remark-math";
import rehypeKatex from "rehype-katex";
import "katex/dist/katex.min.css";

import { Triangle, Circle, Square, Repeat, Zap, BookOpen, Code, ChevronDown, ChevronRight } from "lucide-react";
import BottomBar from "../components/BottomBar";

// UI primitive fallbacks (try to load your project's components; fall back to simple ones)
let Card, CardHeader, CardContent, CardTitle, Input, Button;
try {
  // eslint-disable-next-line
  Card = require("@/components/ui/card").Card;
  // eslint-disable-next-line
  CardHeader = require("@/components/ui/card").CardHeader;
  // eslint-disable-next-line
  CardContent = require("@/components/ui/card").CardContent;
  // eslint-disable-next-line
  CardTitle = require("@/components/ui/card").CardTitle;
  // eslint-disable-next-line
  Input = require("@/components/ui/input").Input;
  // eslint-disable-next-line
  Button = require("@/components/ui/button").Button;
} catch (e) {
  // simple fallbacks
  Card = ({ children, className, style }) => <div className={className} style={style}>{children}</div>;
  CardHeader = ({ children }) => <div>{children}</div>;
  CardContent = ({ children }) => <div>{children}</div>;
  CardTitle = ({ children }) => <div>{children}</div>;
  Input = (props) => <input {...props} />;
  Button = (props) => <button {...props} />;
}

// ---------------------------------------------------------------------
// THEME & ANIMATION
// ---------------------------------------------------------------------
const THEME = {
  bg: "bg-zinc-900",
  panelBg: "#0b0f12",
  text: "#e6eef6",
  muted: "#9ca3af",
  accent: "#16a34a", // emerald
  accent2: "#d97706", // amber
};

const fadeUp = {
  initial: { opacity: 0, y: 8 },
  animate: { opacity: 1, y: 0 },
  transition: { duration: 0.28 },
};

// ---------------------------------------------------------------------
// HELPERS: safe evaluator for user-entered functions
// ---------------------------------------------------------------------
// We provide a lightweight, best-effort evaluator using Function.
// This is not sandboxed. If you need a sandbox, replace with a safe parser (mathjs).
function makeFun(expr, fallback = null) {
  try {
    // Allow Math.* functions without prefix
    const wrapped = `
      const { sin, cos, tan, asin, acos, atan, exp, log, sqrt, pow, abs, PI, E } = Math;
      return (x) => { return ${expr}; }
    `;
    return new Function(wrapped)();
  } catch (e) {
    return fallback;
  }
}

// numerical differentiation (central finite difference)
function derivative(f, x, h = 1e-6) {
  return (f(x + h) - f(x - h)) / (2 * h);
}
function secondDerivative(f, x, h = 1e-4) {
  return (f(x + h) - 2 * f(x) + f(x - h)) / (h * h);
}

// clamp helper
function clamp(v, a, b) {
  return Math.max(a, Math.min(b, v));
}

// format vector nicely
function fmt(x) {
  if (Array.isArray(x)) return `[${x.map((v) => Number(v).toPrecision(6)).join(", ")}]`;
  return Number(x).toPrecision(6);
}

// ---------------------------------------------------------------------
// DOCUMENTATION STRINGS (KaTeX-ready)
const docs = {
  "13.1 Golden-Section Search": [
    `Golden-section search is a derivative-free bracket-shrinking method for unimodal functions on an interval \\([a,b]\\).`,
    `Using the golden ratio \\(\\varphi = \\frac{1+\\sqrt{5}}{2}\\), choose interior points\n\n$$x_1 = b - \\frac{b-a}{\\varphi}, \\qquad x_2 = a + \\frac{b-a}{\\varphi}.$$`,
    `At each iteration remove the half-interval that cannot contain the minimizer and reuse one interior sample to avoid re-evaluation.`,
  ],
  "13.2 Parabolic Interpolation": [
    `Fit a parabola through three points \\((x_0,f_0),(x_1,f_1),(x_2,f_2)\\) and take the vertex as a candidate minimizer. Parabolic interpolation is fast for smooth functions but may be unstable if points are ill-conditioned.`,
  ],
  "13.3 Newton's Method (optimization)": [
    `Newton minimization uses derivatives:\n\n$$x_{k+1} = x_k - \\frac{f'(x_k)}{f''(x_k)}.$$`,
    `Quadratic convergence near a non-degenerate minimum when \\(f''(r)>0\\). Requires first and second derivatives; may diverge if \\(f''\\) is small or negative.`,
  ],
  "13.4 Brent's Method (sketch)": [
    `Brent's method combines golden-section search with parabolic interpolation, using the latter when it appears reliable and falling back to the golden step otherwise. It is robust and efficient in practice.`,
  ],
};

// ---------------------------------------------------------------------
// SECTION: Golden-Section Search implementation & UI
// ---------------------------------------------------------------------
function Section131() {
  // default example: a multimodal-ish but unimodal inside [0,4]
  const [expr, setExpr] = useState(" (x-2)*(x-2) + Math.sin(3*x)*0.03 ");
  const [a, setA] = useState(0);
  const [b, setB] = useState(4);
  const [tol, setTol] = useState(1e-5);
  const [maxIter, setMaxIter] = useState(60);

  const f = useMemo(() => makeFun(expr), [expr]);

  // golden-section implementation returning trace
  function goldenSearch(f, a0, b0, tol = 1e-6, maxIter = 100) {
    const phi = (1 + Math.sqrt(5)) / 2;
    let a = a0, b = b0;
    let x1 = b - (b - a) / phi;
    let x2 = a + (b - a) / phi;
    let f1 = f(x1), f2 = f(x2);
    const trace = [{ a, b, x1, x2, f1, f2 }];
    for (let k = 0; k < maxIter && (b - a) > tol; k++) {
      if (f1 > f2) {
        a = x1;
        x1 = x2;
        f1 = f2;
        x2 = a + (b - a) / phi;
        f2 = f(x2);
      } else {
        b = x2;
        x2 = x1;
        f2 = f1;
        x1 = b - (b - a) / phi;
        f1 = f(x1);
      }
      trace.push({ a, b, x1, x2, f1, f2 });
    }
    const xmin = (a + b) / 2;
    return { a, b, xmin, fmin: f(xmin), trace };
  }

  const result = useMemo(() => {
    if (!f) return null;
    try {
      return goldenSearch(f, Number(a), Number(b), Number(tol), Number(maxIter));
    } catch (e) {
      return null;
    }
  }, [f, a, b, tol, maxIter]);

  // prepare plot data: function samples + trace intervals markers
  const plotData = useMemo(() => {
    if (!f) return null;
    const left = Number(a), right = Number(b);
    const N = 600;
    const xs = [];
    const ys = [];
    for (let i = 0; i <= N; i++) {
      const x = left + (right - left) * (i / N);
      xs.push(x);
      try { ys.push(f(x)); } catch { ys.push(NaN); }
    }
    const intervals = result ? result.trace.map((t) => ({ a: t.a, b: t.b })) : [];
    const points = result ? result.trace.map((t) => ({ x: (t.x1 + t.x2) / 2, y: (t.f1 + t.f2) / 2 })) : [];
    return { xs, ys, intervals, points };
  }, [f, a, b, result]);

  return (
    <section className="rounded-lg p-4" style={{ background: THEME.panelBg }}>
      <motion.div {...fadeUp}>
        <div className="flex items-center gap-3 mb-3">
          <Triangle className="w-6 h-6" />
          <div>
            <h3 className="text-lg font-semibold" style={{ color: THEME.text }}>13.1 Golden-Section Search</h3>
            <div className="text-sm" style={{ color: THEME.muted }}>Derivative-free bracket-shrinking method, ideal for unimodal functions.</div>
          </div>
        </div>

        <div className="rounded bg-zinc-800 border border-zinc-700 p-3 mb-3">
          <div className="grid grid-cols-1 overflow-auto md:grid-cols-3 gap-3">
            <div className="md:col-span-2">
              <label className="text-xs text-zinc-400">Function f(x)</label>
              <Input value={expr} onChange={(e) => setExpr(e.target.value)} className="w-full bg-zinc-900 p-2 rounded text-zinc-100" />
              <div className="text-xs text-zinc-400 mt-1">Use JS Math functions (sin, cos, exp, log, etc.). Example: <code>(x-2)*(x-2)+Math.sin(3*x)*0.03</code></div>
            </div>

            <div>
              <div className="grid grid-cols-2 gap-2">
                <div>
                  <label className="text-xs text-zinc-400">a</label>
                  <Input value={a} onChange={(e) => setA(e.target.value)} className="bg-zinc-900 p-2 rounded text-zinc-100" />
                </div>
                <div>
                  <label className="text-xs text-zinc-400">b</label>
                  <Input value={b} onChange={(e) => setB(e.target.value)} className="bg-zinc-900 p-2 rounded text-zinc-100" />
                </div>
              </div>

              <div className="grid grid-cols-2 overflow-auto gap-2 mt-2">
                <div>
                  <label className="text-xs text-zinc-400">tol</label>
                  <Input value={tol} onChange={(e) => setTol(Number(e.target.value))} className="bg-zinc-900 p-2 rounded text-zinc-100" />
                </div>
                <div>
                  <label className="text-xs text-zinc-400">maxIter</label>
                  <Input value={maxIter} onChange={(e) => setMaxIter(Number(e.target.value))} className="bg-zinc-900 p-2 rounded text-zinc-100" />
                </div>
              </div>
            </div>
          </div>
        </div>

        <div className="rounded border border-zinc-700 bg-zinc-900 p-3 mb-3">
          <div className="text-sm text-zinc-200 mb-2">Visualization</div>
          {plotData ? (
            <Plot
              data={[
                { x: plotData.xs, y: plotData.ys, mode: "lines", name: "f(x)" },
                // plot interval midpoints as markers along iterations
                ...(plotData.points.length ? [{ x: plotData.points.map(p => p.x), y: plotData.points.map(p => p.y), mode: "markers+lines", name: "interval centers", marker: { size: 8 } }] : []),
                // show final xmin
                ...(result ? [{ x: [result.xmin], y: [result.fmin], mode: "markers", name: "xmin", marker: { size: 12, symbol: "diamond" } }] : []),
              ]}
              layout={{
                autosize: true,
                height: 360,
                margin: { l: 40, r: 10, t: 30, b: 40 },
                paper_bgcolor: "rgba(0,0,0,0)",
                plot_bgcolor: "rgba(0,0,0,0)",
                xaxis: { title: "x" },
                yaxis: { title: "f(x)" },
              }}
              useResizeHandler
              style={{ width: "100%" }}
            />
          ) : (
            <div className="text-sm text-zinc-400">Invalid function or parameters — check expression and interval.</div>
          )}
        </div>

        <div className="rounded bg-zinc-800 border border-zinc-700 p-3">
          <div className="text-sm text-zinc-200 mb-2">Result</div>
          {result ? (
            <div className="grid grid-cols-1 md:grid-cols-3 gap-3">
              <div className="rounded p-2 bg-zinc-900 border border-zinc-700">
                <div className="text-xs text-zinc-400">xmin</div>
                <div className="text-sm text-zinc-100">{fmt(result.xmin)}</div>
              </div>
              <div className="rounded p-2 bg-zinc-900 border border-zinc-700">
                <div className="text-xs text-zinc-400">f(xmin)</div>
                <div className="text-sm text-zinc-100">{fmt(result.fmin)}</div>
              </div>
              <div className="rounded p-2 bg-zinc-900 border border-zinc-700">
                <div className="text-xs text-zinc-400">iterations</div>
                <div className="text-sm text-zinc-100">{result.trace.length}</div>
              </div>
            </div>
          ) : (
            <div className="text-sm text-zinc-400">No result (check function or interval)</div>
          )}
        </div>

        <div className="mt-3">
          <Exercise
            title="Exercise 13.1.1 — Golden steps"
            body={`Run golden-section on the provided function and inspect how the interval [a,b] shrinks. Why do we reuse one interior point at each step?`}
            solution={`Reusing one interior point avoids re-evaluating f at two new points each iteration — only one new evaluation is required, saving function calls.`}
          />
        </div>
      </motion.div>
    </section>
  );
}

// ---------------------------------------------------------------------
// SECTION: Parabolic Interpolation
// ---------------------------------------------------------------------
function Section132() {
  const [expr, setExpr] = useState(" (x-2)*(x-2) + 0.05*Math.cos(5*x) ");
  const [x0, setX0] = useState(0.5);
  const [x1, setX1] = useState(2.0);
  const [x2, setX2] = useState(3.5);
  const [maxIter, setMaxIter] = useState(30);
  const [tol, setTol] = useState(1e-6);

  const f = useMemo(() => makeFun(expr), [expr]);

  // parabolic interpolation step given three points returns vertex
  function parabolicVertex(x0, f0, x1, f1, x2, f2) {
    // explicit formula for parabola vertex:
    // numerator = ( (x1-x0)^2*(f2-f1) - (x2-x1)^2*(f1-f0) )
    // denom = ( (x1-x0)*(f2-f1) - (x2-x1)*(f1-f0) )
    // but safer: compute using Lagrange interpolation and derivative = 0
    const A = [
      [x0 * x0, x0, 1],
      [x1 * x1, x1, 1],
      [x2 * x2, x2, 1],
    ];
    const b = [f0, f1, f2];
    // Solve 3x3 linear system for parabola coefficients [a,b,c] where p(x)=ax^2+bx+c
    // Use Cramer's rule or simple elimination (explicit small system)
    const det = (m) =>
      m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
      m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
      m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    const detA = det(A);
    if (Math.abs(detA) < 1e-18) return null;
    const replace = (col, vec) => A.map((row, i) => row.map((val, j) => (j === col ? vec[i] : val)));
    const a = det(replace(0, b)) / detA;
    const bb = det(replace(1, b)) / detA;
    const c = det(replace(2, b)) / detA;
    // vertex x = -b/(2a)
    if (Math.abs(a) < 1e-18) return null;
    return -bb / (2 * a);
  }

  function parabolicInterpolation(f, x0, x1, x2, tol = 1e-6, maxIter = 30) {
    let a = x0, b = x2;
    let xi0 = x0, xi1 = x1, xi2 = x2;
    let f0 = f(xi0), f1 = f(xi1), f2 = f(xi2);
    const trace = [{ xi0, xi1, xi2, f0, f1, f2 }];
    for (let k = 0; k < maxIter; k++) {
      const xv = parabolicVertex(xi0, f0, xi1, f1, xi2, f2);
      if (xv === null || !isFinite(xv)) break;
      const fv = f(xv);
      // replace the worst of the three points
      const worstIndex = ([f0, f1, f2].map((v, i) => ({ v, i })).sort((a, b) => b.v - a.v))[0].i;
      if (worstIndex === 0) {
        xi0 = xv; f0 = fv;
      } else if (worstIndex === 1) {
        xi1 = xv; f1 = fv;
      } else {
        xi2 = xv; f2 = fv;
      }
      trace.push({ xi0, xi1, xi2, f0, f1, f2, xv, fv });
      const interval = Math.max(Math.abs(xi2 - xi1), Math.abs(xi1 - xi0));
      if (interval < tol) break;
    }
    // choose best point
    const pts = [[xi0, f0], [xi1, f1], [xi2, f2]];
    const best = pts.reduce((p, c) => (c[1] < p[1] ? c : p), pts[0]);
    return { xmin: best[0], fmin: best[1], trace };
  }

  const result = useMemo(() => {
    if (!f) return null;
    try {
      return parabolicInterpolation(f, Number(x0), Number(x1), Number(x2), Number(tol), Number(maxIter));
    } catch (e) {
      return null;
    }
  }, [f, x0, x1, x2, tol, maxIter]);

  const plotData = useMemo(() => {
    if (!f) return null;
    // choose plotting interval around provided x0..x2
    const left = Math.min(Number(x0), Number(x1), Number(x2)) - 1;
    const right = Math.max(Number(x0), Number(x1), Number(x2)) + 1;
    const N = 600;
    const xs = [], ys = [];
    for (let i = 0; i <= N; i++) {
      const x = left + (right - left) * (i / N);
      xs.push(x);
      try { ys.push(f(x)); } catch { ys.push(NaN); }
    }
    return { xs, ys };
  }, [f, x0, x1, x2]);

  return (
    <section className="rounded-lg p-4" style={{ background: THEME.panelBg }}>
      <motion.div {...fadeUp}>
        <div className="flex items-center gap-3 mb-3">
          <Circle className="w-6 h-6" />
          <div>
            <h3 className="text-lg font-semibold" style={{ color: THEME.text }}>13.2 Parabolic Interpolation</h3>
            <div className="text-sm" style={{ color: THEME.muted }}>Fit a parabola through three points and use its vertex as a candidate minimizer.</div>
          </div>
        </div>

        <div className="rounded bg-zinc-800 border border-zinc-700 p-3 mb-3">
          <div className="grid grid-cols-1 overflow-auto md:grid-cols-3 gap-3">
            <div className="md:col-span-2">
              <label className="text-xs text-zinc-400">Function f(x)</label>
              <Input value={expr} onChange={(e) => setExpr(e.target.value)} className="w-full bg-zinc-900 p-2 rounded text-zinc-100" />
            </div>

            <div>
              <div className="grid grid-cols-3 overflow-auto gap-2">
                <div>
                  <label className="text-xs text-zinc-400">x0</label>
                  <Input value={x0} onChange={(e) => setX0(e.target.value)} className="bg-zinc-900 p-2 rounded text-zinc-100" />
                </div>
                <div>
                  <label className="text-xs text-zinc-400">x1</label>
                  <Input value={x1} onChange={(e) => setX1(e.target.value)} className="bg-zinc-900 p-2 rounded text-zinc-100" />
                </div>
                <div>
                  <label className="text-xs text-zinc-400">x2</label>
                  <Input value={x2} onChange={(e) => setX2(e.target.value)} className="bg-zinc-900 p-2 rounded text-zinc-100" />
                </div>
              </div>

              <div className="grid grid-cols-2 overflow-auto gap-2 mt-2">
                <div>
                  <label className="text-xs text-zinc-400">tol</label>
                  <Input value={tol} onChange={(e) => setTol(Number(e.target.value))} className="bg-zinc-900 p-2 rounded text-zinc-100" />
                </div>
                <div>
                  <label className="text-xs text-zinc-400">maxIter</label>
                  <Input value={maxIter} onChange={(e) => setMaxIter(Number(e.target.value))} className="bg-zinc-900 p-2 rounded text-zinc-100" />
                </div>
              </div>
            </div>
          </div>
        </div>

        <div className="rounded border border-zinc-700 bg-zinc-900 p-3 mb-3">
          <div className="text-sm text-zinc-200 mb-2">Visualization</div>
          {plotData ? (
            <Plot
              data={[
                { x: plotData.xs, y: plotData.ys, mode: "lines", name: "f(x)" },
                ...(result && result.trace.length ? [{ x: result.trace.map((t, i) => t.xv ?? ((t.xi0 + t.xi1 + t.xi2)/3)).filter(Boolean), y: result.trace.map((t) => t.fv ?? Math.min(t.f0 ?? Infinity, t.f1 ?? Infinity, t.f2 ?? Infinity)).filter(Boolean), mode: "markers+lines", name: "parabolic candidates", marker: { size: 8 } }] : []),
                ...(result ? [{ x: [result.xmin], y: [result.fmin], mode: "markers", marker: { size: 10, symbol: "diamond" }, name: "xmin" }] : []),
              ]}
              layout={{
                autosize: true,
                height: 360,
                margin: { l: 40, r: 10, t: 30, b: 40 },
                paper_bgcolor: "rgba(0,0,0,0)",
                plot_bgcolor: "rgba(0,0,0,0)",
              }}
              useResizeHandler
              style={{ width: "100%" }}
            />
          ) : (
            <div className="text-sm text-zinc-400">Invalid function or inputs.</div>
          )}
        </div>

        <div className="rounded bg-zinc-800 border border-zinc-700 p-3">
          {result ? (
            <div className="grid md:grid-cols-3 gap-3">
              <div className="rounded p-2 bg-zinc-900 border border-zinc-700">
                <div className="text-xs text-zinc-400">xmin</div>
                <div className="text-sm text-zinc-100">{fmt(result.xmin)}</div>
              </div>
              <div className="rounded p-2 bg-zinc-900 border border-zinc-700">
                <div className="text-xs text-zinc-400">f(xmin)</div>
                <div className="text-sm text-zinc-100">{fmt(result.fmin)}</div>
              </div>
              <div className="rounded p-2 bg-zinc-900 border border-zinc-700">
                <div className="text-xs text-zinc-400">iterations</div>
                <div className="text-sm text-zinc-100">{result.trace.length}</div>
              </div>
            </div>
          ) : (
            <div className="text-sm text-zinc-400">No result.</div>
          )}
        </div>

        <div className="mt-3">
          <Exercise
            title="Exercise 13.2.1 — Stability of parabola fit"
            body={`Why might parabolic interpolation produce erratic steps on noisy functions? How can you guard against ill-conditioned fits?`}
            solution={`Ill-conditioning occurs when sample x are nearly collinear in the Vandermonde system; use safeguards: bracket checks, fallback to golden step, or regularization of fit.`}
          />
        </div>
      </motion.div>
    </section>
  );
}

// ---------------------------------------------------------------------
// SECTION: Newton's method for optimization
// ---------------------------------------------------------------------
function Section133() {
  const [expr, setExpr] = useState(" (x-1.5)*(x-1.5) + 0.1*Math.cos(6*x) ");
  const [x0, setX0] = useState(0.5);
  const [tol, setTol] = useState(1e-8);
  const [maxIter, setMaxIter] = useState(50);
  const [useAnalytic, setUseAnalytic] = useState(false);
  const [dExpr, setDExpr] = useState(""); // optional analytic derivative
  const [ddExpr, setDdExpr] = useState("");

  const f = useMemo(() => makeFun(expr), [expr]);
  const df = useMemo(() => (dExpr.trim() ? makeFun(dExpr) : null), [dExpr]);
  const ddf = useMemo(() => (ddExpr.trim() ? makeFun(ddExpr) : null), [ddExpr]);

  function newtonOptimize(f, x0, tol = 1e-8, maxIter = 50, dfAnalytic = null, ddfAnalytic = null) {
    let x = Number(x0);
    const trace = [{ k: 0, x, fx: f(x) }];
    for (let k = 0; k < maxIter; k++) {
      const f1 = dfAnalytic ? dfAnalytic(x) : derivative(f, x, 1e-6);
      const f2 = ddfAnalytic ? ddfAnalytic(x) : secondDerivative(f, x, 1e-4);
      if (!isFinite(f1) || !isFinite(f2)) break;
      // if second derivative near zero, break (can't trust Newton)
      if (Math.abs(f2) < 1e-12) break;
      const dx = f1 / f2;
      const xnew = x - dx;
      trace.push({ k: k + 1, x: xnew, fx: f(xnew), dx });
      if (Math.abs(xnew - x) < tol) {
        x = xnew;
        break;
      }
      x = xnew;
    }
    const xmin = x;
    return { xmin, fmin: f(xmin), trace };
  }

  const result = useMemo(() => {
    if (!f) return null;
    try {
      return newtonOptimize(f, Number(x0), Number(tol), Number(maxIter), useAnalytic ? df : null, useAnalytic ? ddf : null);
    } catch (e) {
      return null;
    }
  }, [f, x0, tol, maxIter, useAnalytic, df, ddf]);

  const plotData = useMemo(() => {
    if (!f) return null;
    const left = Number(x0) - 5;
    const right = Number(x0) + 5;
    const N = 600;
    const xs = [], ys = [];
    for (let i = 0; i <= N; i++) {
      const x = left + (right - left) * (i / N);
      xs.push(x);
      try { ys.push(f(x)); } catch { ys.push(NaN); }
    }
    return { xs, ys };
  }, [f, x0]);

  return (
    <section className="rounded-lg p-4" style={{ background: THEME.panelBg }}>
      <motion.div {...fadeUp}>
        <div className="flex items-center gap-3 mb-3">
          <Zap className="w-6 h-6" />
          <div>
            <h3 className="text-lg font-semibold" style={{ color: THEME.text }}>13.3 Newton's Method (optimization)</h3>
            <div className="text-sm" style={{ color: THEME.muted }}>Uses f' and f'' to jump to minima; quadratic convergence near well-behaved minima.</div>
          </div>
        </div>

        <div className="rounded bg-zinc-800 border border-zinc-700 p-3 mb-3">
          <div className="grid grid-cols-1 md:grid-cols-3 overflow-auto gap-3">
            <div className="md:col-span-2">
              <label className="text-xs text-zinc-400">Function f(x)</label>
              <Input value={expr} onChange={(e) => setExpr(e.target.value)} className="w-full bg-zinc-900 p-2 rounded text-zinc-100" />
              <div className="text-xs text-zinc-400 mt-1">Provide analytic derivatives if available (optional) to improve accuracy.</div>
              <div className="mt-2">
                <label className="text-xs text-zinc-400"><input type="checkbox" checked={useAnalytic} onChange={(e) => setUseAnalytic(e.target.checked)} /> use analytic derivatives</label>
              </div>
              {useAnalytic && (
                <div className="grid grid-cols-1 overflow-auto md:grid-cols-2 gap-2 mt-2">
                  <div>
                    <label className="text-xs text-zinc-400">f'(x)</label>
                    <Input value={dExpr} onChange={(e) => setDExpr(e.target.value)} className="w-full bg-zinc-900 p-2 rounded text-zinc-100" />
                  </div>
                  <div>
                    <label className="text-xs text-zinc-400">f''(x)</label>
                    <Input value={ddExpr} onChange={(e) => setDdExpr(e.target.value)} className="w-full bg-zinc-900 p-2 rounded text-zinc-100" />
                  </div>
                </div>
              )}
            </div>

            <div>
              <label className="text-xs text-zinc-400">x0</label>
              <Input value={x0} onChange={(e) => setX0(e.target.value)} className="w-full bg-zinc-900 p-2 rounded text-zinc-100" />
              <div className="grid grid-cols-2 gap-2 mt-2">
                <div>
                  <label className="text-xs text-zinc-400">tol</label>
                  <Input value={tol} onChange={(e) => setTol(Number(e.target.value))} className="bg-zinc-900 p-2 rounded text-zinc-100" />
                </div>
                <div>
                  <label className="text-xs text-zinc-400">maxIter</label>
                  <Input value={maxIter} onChange={(e) => setMaxIter(Number(e.target.value))} className="bg-zinc-900 p-2 rounded text-zinc-100" />
                </div>
              </div>
            </div>
          </div>
        </div>

        <div className="rounded border border-zinc-700 bg-zinc-900 p-3 mb-3">
          <div className="text-sm text-zinc-200 mb-2">Visualization (neighborhood of x0)</div>
          {plotData ? (
            <Plot
              data={[
                { x: plotData.xs, y: plotData.ys, mode: "lines", name: "f(x)" },
                ...(result ? [{ x: result.trace.map(t => t.x), y: result.trace.map(t => t.fx), mode: "markers+lines", name: "iteration", marker: { size: 8 } }] : []),
                ...(result ? [{ x: [result.xmin], y: [result.fmin], mode: "markers", name: "xmin", marker: { size: 12, symbol: "diamond" } }] : []),
              ]}
              layout={{
                autosize: true,
                height: 360,
                margin: { l: 40, r: 10, t: 30, b: 40 },
                paper_bgcolor: "rgba(0,0,0,0)",
                plot_bgcolor: "rgba(0,0,0,0)",
              }}
              useResizeHandler
              style={{ width: "100%" }}
            />
          ) : (
            <div className="text-sm text-zinc-400">Invalid inputs or function.</div>
          )}
        </div>

        <div className="rounded bg-zinc-800 border border-zinc-700 p-3">
          {result ? (
            <div className="grid grid-cols-1 md:grid-cols-3 gap-3">
              <div className="rounded p-2 bg-zinc-900 border border-zinc-700">
                <div className="text-xs text-zinc-400">xmin</div>
                <div className="text-sm text-zinc-100">{fmt(result.xmin)}</div>
              </div>
              <div className="rounded p-2 bg-zinc-900 border border-zinc-700">
                <div className="text-xs text-zinc-400">f(xmin)</div>
                <div className="text-sm text-zinc-100">{fmt(result.fmin)}</div>
              </div>
              <div className="rounded p-2 bg-zinc-900 border border-zinc-700">
                <div className="text-xs text-zinc-400">iterations</div>
                <div className="text-sm text-zinc-100">{result.trace.length}</div>
              </div>
            </div>
          ) : (
            <div className="text-sm text-zinc-400">No result.</div>
          )}
        </div>

        <div className="mt-3">
          <Exercise
            title="Exercise 13.3.1 — Newton pitfalls"
            body={`Try Newton's method on functions with flat curvature (small f'') or near inflection points. What behaviour do you observe? How can damping help?`}
            solution={`When f'' is tiny or negative, Newton step may be huge or go uphill. Damping (line search or step-size control) restricts step size and provides stability.`}
          />
        </div>
      </motion.div>
    </section>
  );
}

// ---------------------------------------------------------------------
// SECTION: Brent-like hybrid (practical)
 // (A simplified, educational hybrid inspired by Brent's approach)
 // We'll implement a simple hybrid: try parabolic step if within bracket and decreasing; otherwise do golden step.
//
function Section134() {
  const [expr, setExpr] = useState(" (x-1.6)*(x-1.6) + 0.05*Math.cos(8*x) ");
  const [a, setA] = useState(-1.0);
  const [b, setB] = useState(4.0);
  const [tol, setTol] = useState(1e-6);
  const [maxIter, setMaxIter] = useState(60);

  const f = useMemo(() => makeFun(expr), [expr]);

  function goldenStep(a, b) {
    const phi = (1 + Math.sqrt(5)) / 2;
    const x1 = b - (b - a) / phi;
    const x2 = a + (b - a) / phi;
    return { x1, x2 };
  }

  function parabolicVertexSimple(x0, f0, x1, f1, x2, f2) {
    // safe vertex formula from three points (see section 13.2)
    const denom = (x1 - x0) * (x2 - x0) * (x2 - x1);
    if (Math.abs(denom) < 1e-18) return null;
    const A = ( (f2 - f0) * (x1 - x0) - (f1 - f0) * (x2 - x0) ) / denom;
    if (Math.abs(A) < 1e-18) return null;
    const B = (f1 - f0 - A * (x1 * x1 - x0 * x0)) / (x1 - x0);
    const xv = -B / (2 * A);
    return xv;
  }

  function brentLike(f, a0, b0, tol = 1e-6, maxIter = 100) {
    let a = a0, b = b0;
    let x = a + (b - a) / 2;
    let w = x, v = x;
    let fx = f(x), fw = fx, fv = fx;
    const phi = (1 + Math.sqrt(5)) / 2;
    const trace = [{ a, b, x, fx }];
    for (let k = 0; k < maxIter && (b - a) > tol; k++) {
      // attempt parabolic interpolation using x,w,v
      let triedParabolic = false;
      let u = null;
      const xv = parabolicVertexSimple(v, fv, w, fw, x, fx);
      if (xv !== null && xv > a && xv < b) {
        triedParabolic = true;
        u = xv;
      }
      if (!triedParabolic) {
        // golden step
        const gs = goldenStep(a, b);
        u = (Math.abs(x - gs.x1) > Math.abs(x - gs.x2)) ? gs.x1 : gs.x2;
      }
      const fu = f(u);
      // update bracket
      if (fu <= fx) {
        if (u < x) b = x; else a = x;
        v = w; fv = fw;
        w = x; fw = fx;
        x = u; fx = fu;
      } else {
        if (u < x) a = u; else b = u;
        if (fu <= fw || w === x) {
          v = w; fv = fw;
          w = u; fw = fu;
        } else if (fu <= fv || v === x || v === w) {
          v = u; fv = fu;
        }
      }
      trace.push({ a, b, x, fx, u, fu, triedParabolic });
    }
    return { xmin: x, fmin: fx, trace };
  }

  const result = useMemo(() => {
    if (!f) return null;
    try {
      return brentLike(f, Number(a), Number(b), Number(tol), Number(maxIter));
    } catch (e) {
      return null;
    }
  }, [f, a, b, tol, maxIter]);

  const plotData = useMemo(() => {
    if (!f) return null;
    const left = Number(a), right = Number(b);
    const N = 700;
    const xs = [], ys = [];
    for (let i = 0; i <= N; i++) {
      const x = left + (right - left) * (i / N);
      xs.push(x);
      try { ys.push(f(x)); } catch { ys.push(NaN); }
    }
    return { xs, ys };
  }, [f, a, b]);

  return (
    <section className="rounded-lg p-4" style={{ background: THEME.panelBg }}>
      <motion.div {...fadeUp}>
        <div className="flex items-center gap-3 mb-3">
          <Repeat className="w-6 h-6" />
          <div>
            <h3 className="text-lg font-semibold" style={{ color: THEME.text }}>13.4 Brent-like Hybrid</h3>
            <div className="text-sm" style={{ color: THEME.muted }}>Attempt parabolic interpolation when safe; otherwise use golden steps (educational simplified variant)</div>
          </div>
        </div>

        <div className="rounded bg-zinc-800 border border-zinc-700 p-3 mb-3">
          <div className="grid grid-cols-1 overflow-auto md:grid-cols-3 gap-3">
            <div className="md:col-span-2">
              <label className="text-xs text-zinc-400">Function f(x)</label>
              <Input value={expr} onChange={(e) => setExpr(e.target.value)} className="w-full bg-zinc-900 p-2 rounded text-zinc-100" />
            </div>
            <div>
              <label className="text-xs text-zinc-400">a</label>
              <Input value={a} onChange={(e) => setA(e.target.value)} className="w-full bg-zinc-900 p-2 rounded text-zinc-100" />
              <label className="text-xs text-zinc-400 mt-2">b</label>
              <Input value={b} onChange={(e) => setB(e.target.value)} className="w-full bg-zinc-900 p-2 rounded text-zinc-100" />
              <div className="grid grid-cols-2 overflow-auto gap-2 mt-2">
                <Input value={tol} onChange={(e) => setTol(Number(e.target.value))} className="bg-zinc-900 p-2 rounded text-zinc-100" />
                <Input value={maxIter} onChange={(e) => setMaxIter(Number(e.target.value))} className="bg-zinc-900 p-2 rounded text-zinc-100" />
              </div>
            </div>
          </div>
        </div>

        <div className="rounded border border-zinc-700 bg-zinc-900 p-3 mb-3">
          <div className="text-sm text-zinc-200 mb-2">Visualization</div>
          {plotData ? (
            <Plot
              data={[
                { x: plotData.xs, y: plotData.ys, mode: "lines", name: "f(x)" },
                ...(result ? [{ x: result.trace.map(t => t.x), y: result.trace.map(t => t.fx), mode: "markers+lines", name: "current best", marker: { size: 8 } }] : []),
                ...(result ? [{ x: [result.xmin], y: [result.fmin], mode: "markers", name: "xmin", marker: { size: 12, symbol: "diamond" } }] : []),
              ]}
              layout={{
                autosize: true,
                height: 360,
                margin: { l: 40, r: 10, t: 30, b: 40 },
                paper_bgcolor: "rgba(0,0,0,0)",
                plot_bgcolor: "rgba(0,0,0,0)",
              }}
              useResizeHandler
              style={{ width: "100%" }}
            />
          ) : <div className="text-sm text-zinc-400">Invalid input.</div>}
        </div>

        <div className="rounded bg-zinc-800 border border-zinc-700 p-3">
          {result ? (
            <div className="grid md:grid-cols-3 gap-3">
              <div className="rounded p-2 bg-zinc-900 border border-zinc-700">
                <div className="text-xs text-zinc-400">xmin</div>
                <div className="text-sm text-zinc-100">{fmt(result.xmin)}</div>
              </div>
              <div className="rounded p-2 bg-zinc-900 border border-zinc-700">
                <div className="text-xs text-zinc-400">f(xmin)</div>
                <div className="text-sm text-zinc-100">{fmt(result.fmin)}</div>
              </div>
              <div className="rounded p-2 bg-zinc-900 border border-zinc-700">
                <div className="text-xs text-zinc-400">iterations</div>
                <div className="text-sm text-zinc-100">{result.trace.length}</div>
              </div>
            </div>
          ) : <div className="text-sm text-zinc-400">No result.</div>}
        </div>

        <div className="mt-3">
          <Exercise
            title="Exercise 13.4.1 — Hybrid robustness"
            body={`Compare the hybrid to golden-section on functions where parabolic steps are unreliable (e.g., noisy or flat regions). Which method is more robust? Why?`}
            solution={`Brent-like hybrids are robust: they accept parabolic steps when reliable and fallback to safe bracket-halving (golden) otherwise. This balances speed and safety.`}
          />
        </div>
      </motion.div>
    </section>
  );
}

// ---------------------------------------------------------------------
// Problems: Chapter 13
// ---------------------------------------------------------------------
function SectionProblems13() {
  const problems = [
  {
    title: "Problem 13.1",
    body: `Use golden-section search to minimize $f(x)= (x-2)^2 + 0.03\\cos(8x)$ on $[0,4]$. Report $x_{min}$ and $f(x_{min})$.`,
    solution: `Answer: approx $x_{min} \\approx 2$ (mod small perturbation); exact value depends on the small cosine term. Use the provided golden-section UI to verify numerically.`,
  },
  {
    title: "Problem 13.2",
    body: `Implement parabolic interpolation and compare the convergence rate with golden-section on a smooth quadratic and on a noisy function. Explain differences.`,
    solution: `Parabolic interpolation often converges superlinearly on smooth functions; on noisy data it may be erratic — golden is stable but slower.`,
  },
  {
    title: "Problem 13.3",
    body: `Show a case where Newton minimization without damping diverges. Implement damping $x_{k+1}=x_k - \\lambda \\tfrac{f'(x_k)}{f''(x_k)}$ and show how choosing $\\lambda<1$ stabilizes the iteration.`,
    solution: `Divergence occurs when $f''$ is near zero or negative; damping reduces step size and prevents overshoot.`,
  },
  {
    title: "Problem 13.4",
    body: `Compare Brent-like hybrid with golden-section in terms of function evaluations on smooth unimodal functions. Which is faster? Why?`,
    solution: `Brent hybrid typically uses fewer function evaluations by using parabolic steps where appropriate, achieving near-superlinear performance while retaining safety.`,
  },
];


  return (
    <section className="rounded-lg p-4" style={{ background: THEME.panelBg }}>
      <motion.div {...fadeUp}>
        <div className="flex items-center gap-3 mb-3">
          <Code className="w-6 h-6" />
          <div>
            <h3 className="text-lg font-semibold" style={{ color: THEME.text }}>Problems (Chapter 13)</h3>
            <div className="text-sm" style={{ color: THEME.muted }}>Practice tasks and applied exploration</div>
          </div>
        </div>

        <div className="grid grid-cols-1 gap-3">
          {problems.map((p, i) => (
            <Exercise key={i} title={p.title} body={p.body} solution={p.solution} />
          ))}
        </div>
      </motion.div>
    </section>
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


// ---------------------------------------------------------------------
// Page: assemble all sections
// ---------------------------------------------------------------------
export default function Chapter13() {
  return (
    <div className={`min-h-screen p-6 ${THEME.bg} text-zinc-100`}>
      <header className="mb-6">
        <div className="flex items-center justify-between">
          <div>
            <h1 className="text-3xl font-bold" style={{ color: THEME.accent }}>One-Dimensional Unconstrained Optimization</h1>
            <div className="text-sm text-zinc-400 mt-1">Golden-section, parabolic interpolation, Newton, and Brent-like hybrids — visualization-first.</div>
          </div>
        </div>
      </header>

      <main className="space-y-6">
        <Section131 />
        <Section132 />
        <Section133 />
        <Section134 />
        <SectionProblems13 />
        <BottomBar/>
      </main>

      <footer className="mt-8 text-xs text-zinc-500">Tip: For robust optimization in production, use tested libraries (e.g., Brent's method in SciPy). This page is for exploration and learning.</footer>
    </div>
  );
}
