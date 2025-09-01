// src/pages/Chapter17.jsx
// ======================================================================
// Chapter 17 — Least-Squares Regression (responsive, one-file)
// - Linear / Polynomial / Multiple / General / Nonlinear
// - Responsive layout: fluid widths, min-w-0, plot heights use viewport units
// - No canvas/draft mode used — full code provided below
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
import { Button } from "@/components/ui/button";
import { Textarea } from "@/components/ui/textarea";

import {
  BookOpen,
  List,
  ChevronRight,
  ChevronDown,
  Play,
  Calculator,
  TrendingUp,
  BarChart2,
  Grid,
} from "lucide-react";
import BottomBar from "../components/BottomBar";

// ======================================================================
// Theme helpers (keeps same look as other chapters)
// ======================================================================
const theme = {
  bg: "bg-zinc-950",
  panel: "bg-zinc-900/60 border border-zinc-700",
  text: "text-zinc-200",
  subtext: "text-zinc-400",
  accent: "#22d3ee", // cyan-400
  accent2: "#34d399", // emerald-400
  warn: "#f59e0b",
  danger: "#ef4444",
};

const fadeUp = {
  initial: { opacity: 0, y: 12 },
  animate: { opacity: 1, y: 0 },
  transition: { duration: 0.28 },
};

// ======================================================================
// Utilities: parsing, numerics, QR, LM etc.
// (kept intentionally small / self-contained for portability)
// ======================================================================
function parsePairs(s) {
  const rows = String(s)
    .trim()
    .split(/[\n;]+/g)
    .map((r) => r.trim())
    .filter(Boolean);
  const xs = [];
  const ys = [];
  for (const row of rows) {
    const parts = row.split(/[, \t]+/g).filter(Boolean).map(Number);
    if (parts.length >= 2) {
      xs.push(parts[0]);
      ys.push(parts[1]);
    }
  }
  return { x: xs, y: ys };
}

function linspace(a, b, n) {
  if (n <= 1) return [a];
  const dx = (b - a) / (n - 1);
  const arr = new Array(n);
  for (let i = 0; i < n; i++) arr[i] = a + i * dx;
  return arr;
}

function zeros(n, m) {
  if (m == null) return Array.from({ length: n }, () => 0);
  return Array.from({ length: n }, () => Array.from({ length: m }, () => 0));
}

function transpose(A) {
  if (!A || A.length === 0) return [];
  const m = A.length;
  const n = A[0].length;
  const B = Array.from({ length: n }, () => Array.from({ length: m }, () => 0));
  for (let i = 0; i < m; i++) for (let j = 0; j < n; j++) B[j][i] = A[i][j];
  return B;
}

function matMul(A, B) {
  const m = A.length;
  const p = A[0].length;
  const n = B[0].length;
  const C = zeros(m, n);
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      let s = 0;
      for (let k = 0; k < p; k++) s += A[i][k] * B[k][j];
      C[i][j] = s;
    }
  }
  return C;
}

function matVec(A, x) {
  const m = A.length;
  const n = A[0].length;
  const y = new Array(m).fill(0);
  for (let i = 0; i < m; i++) {
    let s = 0;
    for (let j = 0; j < n; j++) s += A[i][j] * x[j];
    y[i] = s;
  }
  return y;
}

function vecSub(a, b) {
  return a.map((v, i) => v - b[i]);
}
function vecAdd(a, b) {
  return a.map((v, i) => v + b[i]);
}
function dot(a, b) {
  let s = 0;
  for (let i = 0; i < a.length; i++) s += a[i] * b[i];
  return s;
}
function norm2(a) {
  return Math.sqrt(dot(a, a));
}

function solveGEPP(Ain, bin) {
  const A = Ain.map((r) => r.slice());
  const b = bin.slice();
  const n = A.length;
  for (let k = 0; k < n - 1; k++) {
    // pivot selection
    let piv = k;
    let maxv = Math.abs(A[k][k]);
    for (let i = k + 1; i < n; i++) {
      const v = Math.abs(A[i][k]);
      if (v > maxv) {
        maxv = v;
        piv = i;
      }
    }
    if (piv !== k) {
      const t = A[k]; A[k] = A[piv]; A[piv] = t;
      const tb = b[k]; b[k] = b[piv]; b[piv] = tb;
    }
    if (!Number.isFinite(A[k][k]) || Math.abs(A[k][k]) < 1e-16) return { x: Array(n).fill(NaN), ok: false };
    for (let i = k + 1; i < n; i++) {
      const f = A[i][k] / A[k][k];
      A[i][k] = 0;
      for (let j = k + 1; j < n; j++) A[i][j] -= f * A[k][j];
      b[i] -= f * b[k];
    }
  }
  if (Math.abs(A[n - 1][n - 1]) < 1e-18) return { x: Array(n).fill(NaN), ok: false };
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    let s = b[i];
    for (let j = i + 1; j < n; j++) s -= A[i][j] * x[j];
    x[i] = s / A[i][i];
  }
  return { x, ok: true };
}

// Simple QR (thin) via Householder (keeps code self-contained)
// Returns {Q, R} where Q is m x n and R is n x n
function qrHouseholder(Ain) {
  const m = Ain.length;
  const n = Ain[0].length;
  let A = Ain.map((r) => r.slice());
  const Q = zeros(m, n);
  // initialize Q as first n columns of identity
  for (let i = 0; i < m; i++) for (let j = 0; j < n; j++) Q[i][j] = i === j ? 1 : 0;

  for (let k = 0; k < n; k++) {
    // compute norm of column below k
    let normx = 0;
    for (let i = k; i < m; i++) normx += A[i][k] * A[i][k];
    normx = Math.sqrt(normx);
    if (normx === 0) continue;
    const alpha = A[k][k] >= 0 ? -normx : normx;
    const v = new Array(m).fill(0);
    const denom = Math.sqrt(0.5 * (alpha * alpha - A[k][k] * alpha));
    if (denom === 0) continue;
    v[k] = (A[k][k] - alpha) / (2 * denom);
    for (let i = k + 1; i < m; i++) v[i] = A[i][k] / (2 * denom);

    // apply to A
    for (let j = k; j < n; j++) {
      let s = 0;
      for (let i = k; i < m; i++) s += v[i] * A[i][j];
      for (let i = k; i < m; i++) A[i][j] -= 2 * v[i] * s;
    }
    // apply to Q
    for (let j = 0; j < n; j++) {
      let s = 0;
      for (let i = k; i < m; i++) s += v[i] * Q[i][j];
      for (let i = k; i < m; i++) Q[i][j] -= 2 * v[i] * s;
    }
  }
  const R = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i <= j ? A[i][j] : 0)));
  return { Q, R };
}

function buildPolynomialFeatures(x, degree) {
  const n = x.length;
  const X = zeros(n, degree + 1);
  for (let i = 0; i < n; i++) {
    for (let d = 0; d <= degree; d++) X[i][d] = Math.pow(x[i], d);
  }
  return X;
}

function predict(X, beta) {
  return matVec(X, beta);
}

function rSquared(y, yhat) {
  if (!y || y.length === 0) return NaN;
  const ybar = y.reduce((a, b) => a + b, 0) / y.length;
  let ssr = 0, sst = 0;
  for (let i = 0; i < y.length; i++) {
    ssr += (yhat[i] - ybar) * (yhat[i] - ybar);
    sst += (y[i] - ybar) * (y[i] - ybar);
  }
  return sst === 0 ? 1 : ssr / sst;
}

// Small Levenberg-Marquardt solver (vectorized residuals)
function levenbergMarquardt(x, y, modelFn, jacobianFn, beta0, opts = {}) {
  const maxIt = opts.maxIt || 80;
  let lambda = opts.lambda || 1e-3;
  const tol = opts.tol || 1e-8;
  let beta = beta0.slice();
  const m = y.length;
  for (let it = 0; it < maxIt; it++) {
    const J = jacobianFn(x, beta); // m x p
    const f = y.map((yi, i) => yi - modelFn(x[i], beta)); // residual vector
    const Jt = transpose(J); // p x m
    const JTJ = matMul(Jt, J); // p x p
    const JTf = matVec(Jt, f); // p
    const p = beta.length;
    const A = JTJ.map((row, i) => row.map((v, j) => (i === j ? v + lambda : v)));
    const sol = solveGEPP(A, JTf);
    if (!sol.ok) break;
    const delta = sol.x;
    const betaNew = beta.map((b, i) => b + delta[i]);
    // evaluate error
    const yhat = x.map((xi) => modelFn(xi, beta));
    const yhatNew = x.map((xi) => modelFn(xi, betaNew));
    const err = norm2(vecSub(y, yhat));
    const errNew = norm2(vecSub(y, yhatNew));
    if (errNew < err) {
      beta = betaNew;
      lambda *= 0.1;
    } else {
      lambda *= 10;
    }
    if (norm2(delta) < tol) break;
  }
  return beta;
}

function fmt(v, p = 6) {
  if (!Number.isFinite(v)) return "NaN";
  return Number(v).toFixed(p).replace(/(?:\.0+|(\.\d+?)0+)$/, "$1");
}

function latexVec(v) {
  return `\\begin{bmatrix} ${v.map((x) => fmt(x, 3)).join(" \\\\ ")} \\end{bmatrix}`;
}

// ======================================================================
// Small shared UI components (responsive-friendly)
// ======================================================================
function Equation({ children, color = theme.accent }) {
  return (
    <div className="my-2">
      <div className="px-3 py-2 rounded-lg border border-zinc-700 bg-zinc-900/80">
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

function Note({ children }) {
  return <div className="text-sm text-zinc-300">{children}</div>;
}

function Labeled({ label, children }) {
  return (
    <div className="min-w-0">
      <div className="text-xs text-zinc-400 mb-1">{label}</div>
      {children}
    </div>
  );
}

function Summary({ items }) {
  return (
    <div className="grid grid-cols-1 sm:grid-cols-3 gap-3">
      {items.map((it, i) => (
        <div key={i} className="rounded-lg bg-zinc-900/70 border border-zinc-700 px-3 py-2 flex flex-col min-w-0">
          <div className="text-xs text-zinc-400 truncate">{it.label}</div>
          <div className="text-zinc-100 text-sm truncate">{it.value}</div>
        </div>
      ))}
    </div>
  );
}

function Exercise({ title, body, solution }) {
  const [open, setOpen] = useState(false);
  return (
    <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-3">
      <div className="flex items-center justify-between gap-2">
        <div className="text-zinc-100 font-medium truncate">{title}</div>
        <Button variant="ghost" className="bg-white cursor-pointer hover:bg-gray-200" size="sm" onClick={() => setOpen((o) => !o)}>
          {open ? <ChevronDown className="w-4 h-4" /> : <ChevronRight className="w-4 h-4" />}
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
// 17.1 Linear Regression — OLS
// ======================================================================
function Section171() {
  const defaultData = "0 1\n1 2\n2 2\n3 3\n4 5";
  const [pairs, setPairs] = useState(defaultData);
  const { x, y } = useMemo(() => parsePairs(pairs), [pairs]);

  const X = useMemo(() => buildPolynomialFeatures(x, 1), [x]);
  const sol = useMemo(() => {
    if (!X || !y) return { x: [], ok: false };
    return (X.length && y.length) ? ( (function(){ const s = (function(){ const Xt = transpose(X); const XtX = matMul(Xt, X); const Xty = matVec(Xt, y); return solveGEPP(XtX, Xty); })(); return s; })() ) : { x: [], ok: false };
  }, [X, y]);

  const beta = sol.ok ? sol.x : Array.from({ length: (X[0] ? X[0].length : 2) }, () => NaN);
  const yhat = (X && sol.ok) ? predict(X, beta) : Array.from({ length: y.length }, () => NaN);
  const rsq = y.length ? rSquared(y, yhat) : NaN;

  const lineX = x.length ? linspace(Math.min(...x) - 1, Math.max(...x) + 1, 120) : linspace(-1, 1, 50);
  const lineXmat = buildPolynomialFeatures(lineX, 1);
  const lineY = predict(lineXmat, beta);

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <TrendingUp className="w-5 h-5" />
            17.1 Linear Regression — OLS
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 md:grid-cols-3 gap-3 items-start min-w-0">
            <Labeled label="Data (x y per line)">
              <Textarea value={pairs} onChange={(e) => setPairs(e.target.value)} className="bg-zinc-800 text-white h-36 w-full min-w-0" />
            </Labeled>
            <Summary items={[{ label: "n (points)", value: x.length }, { label: "Parameters", value: beta.length }, { label: "R²", value: fmt(rsq, 3) }]} />
            <div className="text-sm text-zinc-300">
              <div className="mb-1">Model</div>
              <Equation>{"Ordinary least squares: $\\hat{\\beta}=(X^T X)^{-1}X^T y$"}</Equation>
            </div>
          </div>

          <div className="grid grid-cols-1 lg:grid-cols-2 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 min-w-0 overflow-hidden">
              <Plot
                data={[
                  { x: x, y: y, mode: "markers", type: "scatter", name: "data", marker: { size: 6 } },
                  { x: lineX, y: lineY, mode: "lines", type: "scatter", name: "fit", line: { color: theme.accent } },
                ]}
                layout={{
                  title: "Linear fit",
                  autosize: true,
                  paper_bgcolor: "rgba(0,0,0,0)",
                  plot_bgcolor: "rgba(0,0,0,0)",
                  font: { color: "#e5e7eb" },
                  margin: { l: 40, r: 10, b: 40, t: 40 },
                }}
                useResizeHandler
                style={{ width: "100%", height: "min(55vh,420px)" }}
              />
            </div>

            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-3 min-w-0">
              <div className="text-sm text-zinc-300 mb-2">Estimated coefficients</div>
              <div className="text-zinc-100 text-sm truncate">
                {beta.map((b, i) => (
                  <div key={i}>β_{i} = <span style={{ color: theme.accent2 }}>{fmt(b, 6)}</span></div>
                ))}
              </div>
              <div className="mt-3">
                <Summary items={[{ label: "Residual norm", value: fmt(norm2(vecSub(y, yhat))) }, { label: "R²", value: fmt(rsq, 4) }]} />
              </div>
            </div>
          </div>

          <div className="space-y-3">
            <Exercise
              title="Exercise 17.1.1 — Simple fit"
              body="Use the default dataset. Compare OLS line to eyeballed slope. Compute residuals and comment on fit quality."
              solution="Compute β via normal equations; residuals show deviation at some points. R² quantifies fit quality."
            />
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// 17.2 Polynomial Regression
// ======================================================================
function Section172() {
  const defaultPairs = "-3 -15\n-2 -7\n-1 -1\n0 1\n1 3\n2 9\n3 15";
  const [pairs, setPairs] = useState(defaultPairs);
  const [degree, setDegree] = useState(3);
  const { x, y } = useMemo(() => parsePairs(pairs), [pairs]);

  const X = useMemo(() => buildPolynomialFeatures(x, Number(degree)), [x, degree]);
  const sol = useMemo(() => {
    if (!X || !y) return { x: [], ok: false };
    const Xt = transpose(X);
    const XtX = matMul(Xt, X);
    const Xty = matVec(Xt, y);
    return solveGEPP(XtX, Xty);
  }, [X, y]);

  const beta = sol.ok ? sol.x : Array.from({ length: (X[0] ? X[0].length : (Number(degree) + 1)) }, () => NaN);
  const yhat = sol.ok ? predict(X, beta) : Array.from({ length: y.length }, () => NaN);
  const rsq = y.length ? rSquared(y, yhat) : NaN;

  const lineX = x.length ? linspace(Math.min(...x) - 1, Math.max(...x) + 1, 240) : linspace(-1, 1, 50);
  const lineXmat = buildPolynomialFeatures(lineX, Number(degree));
  const lineY = predict(lineXmat, beta);

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <BarChart2 className="w-5 h-5" />
            17.2 Polynomial Regression
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 md:grid-cols-4 gap-3 items-start min-w-0">
            <Labeled label="Data (x y per line)">
              <Textarea value={pairs} onChange={(e) => setPairs(e.target.value)} className="bg-zinc-800 text-white h-36 w-full min-w-0" />
            </Labeled>
            <Labeled label="Polynomial degree">
              <Input value={degree} onChange={(e) => setDegree(Math.max(0, Math.min(14, Number(e.target.value))))} className="bg-zinc-800 text-white" />
            </Labeled>
            <Summary items={[{ label: "n", value: x.length }, { label: "deg", value: degree }, { label: "R²", value: fmt(rsq, 3) }]} />
            <div className="text-sm text-zinc-300">Watch out for overfitting as degree approaches data count.</div>
          </div>

          <div className="grid grid-cols-1 lg:grid-cols-2 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 min-w-0 overflow-hidden">
              <Plot
                data={[
                  { x: x, y: y, mode: "markers", type: "scatter", name: "data", marker: { size: 6 } },
                  { x: lineX, y: lineY, mode: "lines", type: "scatter", name: `poly deg ${degree}`, line: { color: theme.accent } },
                ]}
                layout={{
                  title: "Polynomial fit",
                  autosize: true,
                  paper_bgcolor: "rgba(0,0,0,0)",
                  plot_bgcolor: "rgba(0,0,0,0)",
                  font: { color: "#e5e7eb" },
                  margin: { l: 40, r: 10, b: 40, t: 40 },
                }}
                useResizeHandler
                style={{ width: "100%", height: "min(60vh,520px)" }}
              />
            </div>

            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-3 min-w-0">
              <div className="text-sm text-zinc-300 mb-2">Coefficients (lowest→highest)</div>
              <div className="text-zinc-100 text-xs space-y-1">
                {beta.map((b, i) => (
                  <div key={i} className="truncate">
                    a_{i} = <span style={{ color: theme.accent2 }}>{fmt(b, 6)}</span>
                  </div>
                ))}
              </div>
              <div className="mt-3"><Summary items={[{ label: "Residual norm", value: fmt(norm2(vecSub(y, yhat))) }, { label: "R²", value: fmt(rsq, 4) }]} /></div>
            </div>
          </div>

          <div className="space-y-3">
            <Exercise title="Exercise 17.2.1 — Overfit demo" body="Increase degree to n-1 and inspect coefficients. Why do coefficients explode?" solution="High-degree polynomials can produce large oscillatory coefficients; basis orthogonalization or regularization helps." />
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// 17.3 Multiple Linear Regression
// ======================================================================
function Section173() {
  const defaultPairs = "1 2 3\n2 3 5\n3 5 7\n4 7 9\n5 11 14";
  const [pairs, setPairs] = useState(defaultPairs);
  const parsed = useMemo(() => {
    const rows = String(pairs)
      .trim()
      .split(/[\n;]+/g)
      .map((r) => r.trim())
      .filter(Boolean)
      .map((r) => r.split(/[, \t]+/g).map(Number));
    if (!rows.length) return { X: [], y: [] };
    const X = rows.map((r) => [1, ...r.slice(0, r.length - 1)]);
    const y = rows.map((r) => r[r.length - 1]);
    return { X, y };
  }, [pairs]);

  const sol = useMemo(() => {
    if (!parsed.X.length) return { x: [], ok: false };
    const Xt = transpose(parsed.X);
    const XtX = matMul(Xt, parsed.X);
    const Xty = matVec(Xt, parsed.y);
    return solveGEPP(XtX, Xty);
  }, [parsed]);

  const beta = sol.ok ? sol.x : Array.from({ length: (parsed.X[0] ? parsed.X[0].length : 1) }, () => NaN);
  const yhat = sol.ok ? predict(parsed.X, beta) : Array.from({ length: parsed.y.length }, () => NaN);
  const rsq = parsed.y.length ? rSquared(parsed.y, yhat) : NaN;

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Grid className="w-5 h-5" />
            17.3 Multiple Linear Regression
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 md:grid-cols-3 gap-3 items-start min-w-0">
            <Labeled label="Data (x1 x2 ... y per line)">
              <Textarea value={pairs} onChange={(e) => setPairs(e.target.value)} className="bg-zinc-800 text-white h-36 w-full min-w-0" />
            </Labeled>
            <Summary items={[{ label: "n rows", value: parsed.X.length }, { label: "features", value: (parsed.X[0] ? parsed.X[0].length - 1 : 0) }, { label: "R²", value: fmt(rsq, 3) }]} />
            <div className="text-sm text-zinc-300">Intercept included by default (column of ones).</div>
          </div>

          <div className="grid grid-cols-1 lg:grid-cols-2 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-3 min-w-0">
              <div className="text-sm text-zinc-300 mb-2">Estimated β (intercept first)</div>
              <div className="text-zinc-100 text-xs space-y-1">
                {beta.map((b, i) => (
                  <div key={i} className="truncate">β_{i} = <span style={{ color: theme.accent2 }}>{fmt(b, 6)}</span></div>
                ))}
              </div>
            </div>

            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-3 min-w-0">
              <Summary items={[{ label: "Residual norm", value: fmt(norm2(vecSub(parsed.y, yhat))) }, { label: "R²", value: fmt(rsq, 4) }]} />
              <Note className="mt-2">Use feature scaling if coefficients differ in scale significantly to improve numerical stability.</Note>
            </div>
          </div>

          <div className="space-y-3">
            <Exercise title="Exercise 17.3.1 — Multicollinearity" body="Construct correlated features (x2 ≈ 2 x1) and fit. Inspect coefficients and condition of X^T X." solution="Collinear features produce ill-conditioned normal equations; use QR or regularization." />
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// 17.4 General Linear Least Squares — QR option
// ======================================================================
function Section174() {
  const defaultPairs = "0 1\n1 2.1\n2 1.8\n3 3.3\n4 3.9";
  const [pairs, setPairs] = useState(defaultPairs);
  const [method, setMethod] = useState("normal"); // normal | qr
  const { x, y } = useMemo(() => parsePairs(pairs), [pairs]);
  const X = useMemo(() => buildPolynomialFeatures(x, 1), [x]);

  const sol = useMemo(() => {
    if (!X || !y) return { x: [], ok: false };
    if (method === "normal") {
      const Xt = transpose(X);
      const XtX = matMul(Xt, X);
      const Xty = matVec(Xt, y);
      return solveGEPP(XtX, Xty);
    } else {
      const qq = qrHouseholder(X);
      const Q = qq.Q; // m x n
      const R = qq.R; // n x n
      const Qt = transpose(Q); // n x m
      const qty = matVec(Qt, y); // n
      return solveGEPP(R, qty);
    }
  }, [X, y, method]);

  const beta = sol.ok ? sol.x : Array.from({ length: (X[0] ? X[0].length : 2) }, () => NaN);
  const yhat = sol.ok ? predict(X, beta) : Array.from({ length: (y ? y.length : 0) }, () => NaN);
  const rsq = y.length ? rSquared(y, yhat) : NaN;

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Calculator className="w-5 h-5" />
            17.4 General Linear Least Squares
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 md:grid-cols-3 gap-3 items-start min-w-0">
            <Labeled label="Data (x y per line)">
              <Textarea value={pairs} onChange={(e) => setPairs(e.target.value)} className="bg-zinc-800 text-white h-36 w-full min-w-0" />
            </Labeled>
            <Labeled label="Method">
              <div className="flex gap-2">
                <button className={`px-2 py-1 rounded ${method === "normal" ? "bg-emerald-600 text-white" : "bg-zinc-800 text-zinc-200"}`} onClick={() => setMethod("normal")}>Normal</button>
                <button className={`px-2 py-1 rounded ${method === "qr" ? "bg-emerald-600 text-white" : "bg-zinc-800 text-zinc-200"}`} onClick={() => setMethod("qr")}>QR</button>
              </div>
            </Labeled>
            <Summary items={[{ label: "R²", value: fmt(rsq, 3) }, { label: "Residual", value: fmt(norm2(vecSub(y, yhat))) }, { label: "Method", value: method }]} />
          </div>

          <div className="grid grid-cols-1 lg:grid-cols-2 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 min-w-0 overflow-hidden">
              <Plot
                data={[
                  { x: x, y: y, mode: "markers", type: "scatter", name: "data", marker: { size: 6 } },
                  { x: x, y: yhat, mode: "lines", type: "scatter", name: "fit", line: { color: theme.accent } },
                ]}
                layout={{
                  title: "Fit (compare methods)",
                  autosize: true,
                  paper_bgcolor: "rgba(0,0,0,0)",
                  plot_bgcolor: "rgba(0,0,0,0)",
                  font: { color: "#e5e7eb" },
                  margin: { l: 40, r: 10, b: 40, t: 40 },
                }}
                useResizeHandler
                style={{ width: "100%", height: "min(55vh,460px)" }}
              />
            </div>

            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-3 min-w-0">
              <div className="text-sm text-zinc-300 mb-2">Notes</div>
              <div className="text-zinc-200 text-sm">QR avoids forming X^T X and is numerically more stable when features are correlated. For severely ill-conditioned problems, SVD-based solvers or regularization (Ridge) are recommended.</div>
            </div>
          </div>

          <div className="space-y-3">
            <Exercise title="Exercise 17.4.1 — Compare normal vs QR" body="Construct a dataset with near-collinear columns and compare solutions from normal equations and QR. Which is more stable?" solution="QR tends to be more stable; normal may give large errors due to X^T X conditioning." />
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// 17.5 Nonlinear Regression — Exponential example using LM
// ======================================================================
function Section175() {
  const defaultPairs = "0 2.1\n1 4.3\n2 9.1\n3 20.2\n4 41.3";
  const [pairs, setPairs] = useState(defaultPairs);
  const { x, y } = useMemo(() => parsePairs(pairs), [pairs]);

  function modelExp(xi, beta) {
    const a = beta[0], b = beta[1];
    return a * Math.exp(b * xi);
  }
  function jacExp(xs, beta) {
    const a = beta[0], b = beta[1];
    return xs.map((xi) => [Math.exp(b * xi), a * xi * Math.exp(b * xi)]);
  }

  const betaInit = [1, 0.5];
  const beta = useMemo(() => {
    if (!x || !y || !x.length) return betaInit;
    return levenbergMarquardt(x, y, modelExp, jacExp, betaInit, { maxIt: 80 });
  }, [x, y]);

  const yhat = x.map((xi) => modelExp(xi, beta));
  const rsq = y.length ? rSquared(y, yhat) : NaN;

  const lineX = x.length ? linspace(Math.min(...x) - 0.5, Math.max(...x) + 0.5, 240) : linspace(-1, 1, 50);
  const lineY = lineX.map((xi) => modelExp(xi, beta));

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Play className="w-5 h-5" />
            17.5 Nonlinear Regression — Exponential Example
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 md:grid-cols-3 gap-3 items-start min-w-0">
            <Labeled label="Data (x y per line)">
              <Textarea value={pairs} onChange={(e) => setPairs(e.target.value)} className="bg-zinc-800 text-white h-36 w-full min-w-0" />
            </Labeled>
            <Summary items={[{ label: "Parameters", value: beta.length }, { label: "R²", value: fmt(rsq, 3) }, { label: "Residual", value: fmt(norm2(vecSub(y, yhat))) }]} />
            <div className="text-sm text-zinc-300">Levenberg–Marquardt used for optimization; jacobian provided analytically where possible.</div>
          </div>

          <div className="grid grid-cols-1 lg:grid-cols-2 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 min-w-0 overflow-hidden">
              <Plot
                data={[
                  { x: x, y: y, mode: "markers", type: "scatter", name: "data", marker: { size: 6 } },
                  { x: lineX, y: lineY, mode: "lines", type: "scatter", name: "exp fit", line: { color: theme.accent } },
                ]}
                layout={{
                  title: "Exponential fit",
                  autosize: true,
                  paper_bgcolor: "rgba(0,0,0,0)",
                  plot_bgcolor: "rgba(0,0,0,0)",
                  font: { color: "#e5e7eb" },
                  margin: { l: 40, r: 10, b: 40, t: 40 },
                }}
                useResizeHandler
                style={{ width: "100%", height: "min(55vh,480px)" }}
              />
            </div>

            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-3 min-w-0">
              <div className="text-sm text-zinc-300 mb-2">Estimated params</div>
              <div className="text-zinc-100 text-sm">
                a = <span style={{ color: theme.accent2 }}>{fmt(beta[0], 6)}</span><br />
                b = <span style={{ color: theme.accent2 }}>{fmt(beta[1], 6)}</span>
              </div>
              <div className="mt-3"><Summary items={[{ label: "R²", value: fmt(rsq, 4) }, { label: "Residual", value: fmt(norm2(vecSub(y, yhat))) }]} /></div>
            </div>
          </div>

          <div className="space-y-3">
            <Exercise title="Exercise 17.5.1 — Fit logistic" body="Try logistic model y = L / (1 + exp(-k(x-x0))). Provide jacobian and fit with LM." solution="Logistic can be fit with LM; provide good initial guesses for stability." />
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// Documentation & Problems panels
// ======================================================================
const docs17 = {
  "17.1 Linear Regression": [
    "Ordinary least squares minimizes the sum of squared residuals: $\\min_\\beta \\|y - X\\beta\\|_2^2$.",
    "Closed-form solution via the normal equations: $\\hat{\\beta}=(X^T X)^{-1}X^T y$ when $X^T X$ is invertible.",
    "Geometric view: regression projects $y$ onto the column space of $X$."
  ],
  "17.2 Polynomial Regression": [
    "Polynomial basis expansion turns linear regression into polynomial fitting: $x \\mapsto [1, x, x^2, \\dots, x^d]$.",
    "Higher degrees increase flexibility but may cause overfitting and oscillations (Runge’s phenomenon).",
    "Regularization (Ridge, Lasso) mitigates variance from large coefficients."
  ],
  "17.3 Multiple Linear Regression": [
    "Model: $y = X\\beta + \\varepsilon$, where $X$ includes multiple predictors and possibly an intercept column.",
    "Check for multicollinearity: large variance inflation factors (VIF) indicate near-linear dependence among predictors.",
    "Standardize features if scales differ significantly: $x_j^{(scaled)} = \\frac{x_j - \\mu_j}{\\sigma_j}$."
  ],
  "17.4 General Least Squares": [
    "When errors have covariance matrix $\\Sigma$, use generalized least squares (GLS): $\\hat{\\beta}=(X^T \\Sigma^{-1} X)^{-1}X^T \\Sigma^{-1} y$.",
    "Numerical stability: prefer QR decomposition for solving $X^T X \\beta = X^T y$ instead of forming $(X^T X)^{-1}$ directly.",
    "For rank-deficient $X$, use singular value decomposition (SVD) or apply regularization (Ridge regression)."
  ],
  "17.5 Nonlinear Regression": [
    "Model: $y = f(x,\\theta) + \\varepsilon$, with nonlinear dependence on parameters $\\theta$.",
    "Solve iteratively using Gauss–Newton or Levenberg–Marquardt algorithms.",
    "Provide analytical Jacobian $J = \\frac{\\partial f}{\\partial \\theta}$ when possible for faster convergence and stability.",
    "Initialization matters: poor starting guesses may lead to local minima."
  ],
  "17.6 Regularization Methods": [
    "Ridge regression (L2): $\\min_\\beta \\|y - X\\beta\\|_2^2 + \\lambda \\|\\beta\\|_2^2$ shrinks coefficients toward zero, reducing variance.",
    "Lasso (L1): $\\min_\\beta \\|y - X\\beta\\|_2^2 + \\lambda \\|\\beta\\|_1$ promotes sparsity in coefficients.",
    "Elastic Net: combination of Ridge and Lasso, useful when predictors are correlated."
  ],
  "17.7 Model Evaluation": [
    "Coefficient of determination: $R^2 = 1 - \\frac{\\|y - \\hat{y}\\|_2^2}{\\|y - \\bar{y}\\|_2^2}$.",
    "Use cross-validation to estimate generalization error.",
    "Information criteria: AIC, BIC balance fit quality with model complexity."
  ]
};


function ChapterDocs17() {
  const [open, setOpen] = useState("17.1 Linear Regression");
  return (
    <div className="w-full mt-6 bg-zinc-900 border border-zinc-700 rounded-xl shadow-lg overflow-hidden">
      <div className="p-3 border-b border-zinc-700 flex items-center justify-between">
        <div className="flex items-center gap-3">
          <BookOpen className="w-5 h-5 text-cyan-400" />
          <div>
            <div className="text-zinc-100 font-semibold">Chapter 17 — Documentation</div>
            <div className="text-zinc-400 text-xs">Theory, notes, and references</div>
          </div>
        </div>
        <div className="text-zinc-400 text-xs">Least-squares methods & diagnostics</div>
      </div>

      <div className="grid grid-cols-1 md:grid-cols-4">
        <div className="col-span-1 border-r border-zinc-800">
          {Object.keys(docs17).map((k) => (
            <button key={k} onClick={() => setOpen(k)} className={`w-full p-3 text-left hover:bg-zinc-800/40 flex items-center justify-between ${open === k ? "bg-zinc-800/20" : ""}`}>
              <div className="flex items-center gap-2">
                <List className="w-4 h-4 text-zinc-300" />
                <div className="text-zinc-100 truncate">{k}</div>
              </div>
              {open === k ? <ChevronDown className="w-4 h-4 text-zinc-400" /> : <ChevronRight className="w-4 h-4 text-zinc-400" />}
            </button>
          ))}
        </div>

        <div className="col-span-3 p-4 overflow-auto" style={{ maxHeight: 460 }}>
          <motion.h3 initial={{ opacity: 0, y: -6 }} animate={{ opacity: 1, y: 0 }} className="text-xl font-semibold text-zinc-100 mb-3">{open}</motion.h3>
          <div className="text-zinc-300 space-y-3 text-sm leading-relaxed">
            {docs17[open].map((p, i) => (
              <div key={i}><ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>{p}</ReactMarkdown></div>
            ))}
            <div className="mt-3 text-zinc-400 text-xs">References: Hastie, Tibshirani & Friedman; Seber & Lee; Björck; Golub & Van Loan.</div>
          </div>
        </div>
      </div>
    </div>
  );
}

function Problems17() {
  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Grid className="w-5 h-5" /> Problems
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-3">
          <Exercise title="P17.1 — Linear model diagnostics" body="Given a dataset, compute residuals, plot residuals vs fitted, and assess homoscedasticity." solution="Check residual scatter, QQ plot for normality, and influence points (Cook's distance)." />
          <Exercise title="P17.2 — Polynomial basis" body="Fit polynomial degrees 1..6 to a noisy dataset; select degree by cross-validation." solution="Cross-validated MSE typically gives optimal degree balancing bias-variance." />
          <Exercise title="P17.3 — Ridge regression" body="Implement ridge via (X^T X + λI)^{-1} X^T y and explore effect of λ." solution="Ridge shrinks coefficients and stabilizes ill-conditioned problems; λ chosen by CV." />
          <Exercise title="P17.4 — Nonlinear fit" body="Fit logistic growth to population data using LM and report parameter uncertainties (approx from Jacobian)." solution="Estimate covariance ≈ σ^2 (J^T J)^{-1} at solution; σ^2 from residual variance." />
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// Page assembly (responsive container)
// ======================================================================
export default function Chapter17() {
  return (
    <div className={`p-4 md:p-6 lg:p-8 space-y-6 ${theme.bg} min-h-screen`}>
      <div className="max-w-screen-xl mx-auto">
        <motion.div initial={{ opacity: 0, y: -10 }} animate={{ opacity: 1, y: 0 }}>
          <h1 className="text-2xl md:text-3xl lg:text-4xl font-bold" style={{ color: theme.accent2 }}>Least-Squares Regression</h1>
          <div className="text-zinc-400">Linear, polynomial, multivariate, and nonlinear regression — interactive demos. Responsive across devices.</div>
        </motion.div>

        <div className="space-y-6 mt-4">
          <Section171 />
          <Section172 />
          <Section173 />
          <Section174 />
          <Section175 />
          <ChapterDocs17 />
          <Problems17 />
          <BottomBar/>
        </div>
      </div>
    </div>
  );
}
