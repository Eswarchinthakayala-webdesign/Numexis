// src/pages/Chapter20.jsx
// ======================================================================
// Chapter 20 — Case Studies: Curve Fitting (Interactive, pro-level, responsive)
// Sections:
// 20.1 Linear Regression and Population Models
// 20.2 Use of Splines to Estimate Heat Transfer
// 20.3 Fourier Analysis
// 20.4 Analysis of Experimental Data
// Problems
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
  TrendingUp,
  Thermometer,
  Zap,
  Activity,
  Grid2X2,
} from "lucide-react";
import BottomBar from "../components/BottomBar";

// ======================================================================
// Theme & helpers (matches earlier chapters' style)
// ======================================================================
const theme = {
  bg: "bg-zinc-950",
  panel: "bg-zinc-900/60 border border-zinc-700",
  accent: "#22d3ee",
  accent2: "#34d399",
  text: "text-zinc-200",
  subtext: "text-zinc-400",
};

const fadeUp = {
  initial: { opacity: 0, y: 12 },
  animate: { opacity: 1, y: 0 },
  transition: { duration: 0.28 },
};

function fmt(v, p = 6) {
  if (!Number.isFinite(v)) return "NaN";
  return Number(v).toFixed(p).replace(/(?:\.0+|(\.\d+?)0+)$/, "$1");
}

// ======================================================================
// Small numeric utilities
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

function dot(a, b) {
  let s = 0;
  for (let i = 0; i < a.length; i++) s += a[i] * b[i];
  return s;
}
function vecSub(a, b) {
  return a.map((v, i) => v - b[i]);
}
function norm2(a) {
  return Math.sqrt(dot(a, a));
}

// Simple GEPP solver (for normal equations and small linear systems)
function solveGEPP(Ain, bin) {
  const A = Ain.map((r) => r.slice());
  const b = bin.slice();
  const n = A.length;
  for (let k = 0; k < n - 1; k++) {
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
      const t = A[k];
      A[k] = A[piv];
      A[piv] = t;
      const tb = b[k];
      b[k] = b[piv];
      b[piv] = tb;
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

// Linear regression (ordinary least squares) for design matrix X (m x p) and y (m)
function leastSquares(X, y) {
  const Xt = transpose(X);
  const XtX = matMul(Xt, X);
  const Xty = matVec(Xt, y);
  return solveGEPP(XtX, Xty);
}

// Build Vandermonde-like polynomial features up to degree d
function polyFeatures(x, degree) {
  const n = x.length;
  const X = zeros(n, degree + 1);
  for (let i = 0; i < n; i++) {
    for (let d = 0; d <= degree; d++) X[i][d] = Math.pow(x[i], d);
  }
  return X;
}

// Natural cubic spline building (reused concept from earlier chapters)
function naturalCubicSpline(xs, ys) {
  const n = xs.length;
  if (n < 2) return null;
  const h = new Array(n - 1);
  for (let i = 0; i < n - 1; i++) h[i] = xs[i + 1] - xs[i];
  if (n === 2) return { xs, ys, m: [0, 0], h };
  const a = [];
  const b = [];
  const c = [];
  const d = [];
  for (let i = 1; i < n - 1; i++) {
    a.push(h[i - 1]);
    b.push(2 * (h[i - 1] + h[i]));
    c.push(h[i]);
    d.push(6 * ((ys[i + 1] - ys[i]) / h[i] - (ys[i] - ys[i - 1]) / h[i - 1]));
  }
  // Solve tridiagonal using simple Thomas algorithm
  const n2 = b.length;
  if (n2 === 0) return { xs, ys, m: [0, 0], h };
  const cp = c.slice();
  const bp = b.slice();
  const dp = d.slice();
  for (let i = 1; i < n2; i++) {
    const m = a[i - 1] / bp[i - 1];
    bp[i] = bp[i] - m * cp[i - 1];
    dp[i] = dp[i] - m * dp[i - 1];
  }
  const mInterior = Array(n2).fill(0);
  mInterior[n2 - 1] = dp[n2 - 1] / bp[n2 - 1];
  for (let i = n2 - 2; i >= 0; i--) mInterior[i] = (dp[i] - cp[i] * mInterior[i + 1]) / bp[i];
  const m = [0, ...mInterior, 0];
  return { xs, ys, m, h };
}

function evalNaturalCubic(spline, x) {
  const { xs, ys, m, h } = spline;
  const n = xs.length;
  if (x <= xs[0]) return ys[0];
  if (x >= xs[n - 1]) return ys[n - 1];
  let i = 0;
  for (let k = 0; k < n - 1; k++) if (x >= xs[k] && x <= xs[k + 1]) { i = k; break; }
  const xi = xs[i];
  const xi1 = xs[i + 1];
  const hi = h[i];
  const A = (xi1 - x) / hi;
  const B = (x - xi) / hi;
  const Si = A * ys[i] + B * ys[i + 1] + ((A * (A * A - 1) * m[i] + B * (B * B - 1) * m[i + 1]) * (hi * hi)) / 6;
  return Si;
}

// Discrete Fourier transform (simple O(N^2), fine for moderate N)
function dftReal(x) {
  const N = x.length;
  const re = new Array(N).fill(0);
  const im = new Array(N).fill(0);
  for (let k = 0; k < N; k++) {
    for (let n = 0; n < N; n++) {
      const angle = (-2 * Math.PI * k * n) / N;
      re[k] += x[n] * Math.cos(angle);
      im[k] += x[n] * Math.sin(angle);
    }
  }
  return { re, im };
}

// Inverse DFT real reconstruction
function idftReconstruct(re, im) {
  const N = re.length;
  const x = new Array(N).fill(0);
  for (let n = 0; n < N; n++) {
    let s = 0;
    for (let k = 0; k < N; k++) {
      const angle = (2 * Math.PI * k * n) / N;
      s += re[k] * Math.cos(angle) - im[k] * Math.sin(angle);
    }
    x[n] = s / N;
  }
  return x;
}

// ======================================================================
// Small UI components
// ======================================================================
function Equation({ children, color = theme.accent }) {
  return (
    <div className="my-2">
      <div className="px-3 py-2 rounded-lg border text-gray-300 overflow-auto border-zinc-700 bg-zinc-900/80">
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
// Section 20.1 — Linear Regression and Population Models
// ======================================================================
function Section201() {
  const defaultData = "0 2.5\n1 3.0\n2 3.8\n3 5.1\n4 6.0\n5 7.9";
  const [pairs, setPairs] = useState(defaultData);
  const [degree, setDegree] = useState(1); // linear by default
  const [predictX, setPredictX] = useState(6);
  const { x, y } = useMemo(() => parsePairs(pairs), [pairs]);

  const X = useMemo(() => polyFeatures(x, Number(degree)), [x, degree]);
  const sol = useMemo(() => (X.length && y.length ? leastSquares(X, y) : { x: [], ok: false }), [X, y]);
  const beta = sol.ok ? sol.x : Array.from({ length: (X[0] ? X[0].length : 1) }, () => NaN);
  const yhat = (X && sol.ok) ? matVec(X, beta) : Array.from({ length: y.length }, () => NaN);
  const rss = norm2(vecSub(y, yhat));
  const lineX = x.length ? linspace(Math.min(...x) - 0.5, Math.max(...x) + 1.5, 200) : [];
  const lineXmat = polyFeatures(lineX, Number(degree));
  const lineY = lineXmat.map((row) => dot(row, beta));

  const predictVal = useMemo(() => {
    if (!beta.length) return NaN;
    const row = polyFeatures([Number(predictX)], Number(degree))[0];
    return dot(row, beta);
  }, [predictX, beta, degree]);

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <TrendingUp className="w-5 h-5" />
            20.1 Linear Regression & Population Models
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 md:grid-cols-3 gap-3 items-start min-w-0">
            <Labeled label="Data: time (t) population / measurement">
              <Textarea value={pairs} onChange={(e) => setPairs(e.target.value)} className="bg-zinc-800 text-white h-36 w-full min-w-0" />
            </Labeled>

            <div>
              <Labeled label="Model degree (0 = constant, 1 = linear, 2 = quadratic)">
                <Input value={degree} onChange={(e) => setDegree(Math.max(0, Math.min(6, Number(e.target.value))))} className="bg-zinc-800 text-white w-28" />
              </Labeled>
              <div className="mt-2">
                <Labeled label="Predict at t">
                  <Input value={predictX} onChange={(e) => setPredictX(e.target.value)} className="bg-zinc-800 text-white w-28" />
                </Labeled>
              </div>
            </div>

            <Summary items={[
              { label: "Parameters", value: beta.length },
              { label: "Residual norm (RSS)", value: fmt(rss, 6) },
              { label: "Prediction", value: fmt(predictVal, 6) }
            ]} />
          </div>

          <div className="grid grid-cols-1 lg:grid-cols-2 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 min-w-0 overflow-hidden">
              <Plot
                data={[
                  { x: x, y: y, mode: "markers", type: "scatter", name: "data", marker: { size: 7 } },
                  { x: lineX, y: lineY, mode: "lines", type: "scatter", name: `poly deg ${degree}`, line: { color: theme.accent } },
                  { x: [Number(predictX)], y: [predictVal], mode: "markers", type: "scatter", name: "prediction", marker: { size: 10, color: theme.accent2 } }
                ]}
                layout={{
                  title: "Regression fit (least squares)",
                  autosize: true,
                  paper_bgcolor: "rgba(0,0,0,0)",
                  plot_bgcolor: "rgba(0,0,0,0)",
                  font: { color: "#e5e7eb" },
                  margin: { l: 40, r: 10, b: 40, t: 40 },
                }}
                useResizeHandler
                style={{ width: "100%", height: "min(55vh,520px)" }}
              />
            </div>

            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-3 min-w-0 overflow-auto">
              <div className="text-sm text-zinc-300 mb-2">Estimated coefficients (lowest → highest)</div>
              <div className="text-zinc-100 text-xs space-y-1">
                {beta.map((b, i) => (<div key={i}>β_{i} = <span style={{ color: theme.accent2 }}>{fmt(b, 6)}</span></div>))}
              </div>
              <Equation>{`Least squares: $$\\hat{\\beta}=\\arg\\min_\\beta \\|y-X\\beta\\|_2^2,\\\\ \\text{(closed form via normal equations) }(X^TX)\\hat{\\beta}=X^Ty.$$`}</Equation>
              <Note>For population models, linear fits may be applied to transformed models (e.g., log-transform for exponential growth). See exercises.</Note>
            </div>
          </div>

          <div className="space-y-3">
            <Exercise
              title="Exercise 20.1.1 — Exponential model via transform"
              body="Fit an exponential population model y = A e^{rt} by applying a log-transform and performing linear regression on ln(y). Compare residuals and back-transform predictions."
              solution="Transform ln(y)=ln(A)+r t, fit slope r, intercept ln(A), then exponentiate to get A and predictions. Check residuals on original scale."
            />
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// Section 20.2 — Use of Splines to Estimate Heat Transfer
// ======================================================================
function Section202() {
  // Example: depth vs temperature measurement along wall/ground
  const defaultData = "0 15.2\n1 15.6\n2 16.5\n3 17.8\n4 18.2\n5 18.0\n6 17.5\n7 17.0";
  const [pairs, setPairs] = useState(defaultData);
  const [estimateAt, setEstimateAt] = useState(3.5);
  const { x, y } = useMemo(() => parsePairs(pairs), [pairs]);

  const spline = useMemo(() => (x.length ? naturalCubicSpline(x, y) : null), [x, y]);
  const lineX = x.length ? linspace(Math.min(...x), Math.max(...x), 400) : [];
  const splineY = (spline && spline.m) ? lineX.map((xx) => evalNaturalCubic(spline, xx)) : [];

  const estimateVal = useMemo(() => (spline ? evalNaturalCubic(spline, Number(estimateAt)) : NaN), [spline, estimateAt]);

  // Simple finite difference estimate of heat flux ~ -k dT/dx; approximate dT/dx via central diff on spline
  function splineDerivative(spline, xq) {
    const { xs, ys, m, h } = spline;
    const n = xs.length;
    if (xq <= xs[0]) {
      const i = 0; const hi = h[i];
      const t = (xq - xs[i]) / hi;
      // derivative of cubic piece S_i' formula
      const dS = (ys[i + 1] - ys[i]) / hi - (hi / 6) * (m[i + 1] + 2 * m[i]) + (xq - xs[i]) * 0; // approximate rough
      return dS;
    }
    if (xq >= xs[n - 1]) {
      const i = n - 2; const hi = h[i];
      const dS = (ys[i + 1] - ys[i]) / hi + (hi / 6) * (2 * m[i + 1] + m[i]);
      return dS;
    }
    let i = 0;
    for (let k = 0; k < n - 1; k++) if (xq >= xs[k] && xq <= xs[k + 1]) { i = k; break; }
    const xi = xs[i];
    const xi1 = xs[i + 1];
    const hi = h[i];
    const A = (xi1 - xq) / hi;
    const B = (xq - xi) / hi;
    // derivative formula for cubic spline S'(x)
    const term1 = (ys[i + 1] - ys[i]) / hi;
    const term2 = (m[i + 1] - m[i]) * ( ( ( ( -3 * A * A + 1 ) * (hi) ) / 6 ) ); // approximate factor, simplified
    // To keep things robust we compute derivative numerically via small h
    const delta = 1e-4;
    const f1 = evalNaturalCubic(spline, xq + delta);
    const f0 = evalNaturalCubic(spline, xq - delta);
    return (f1 - f0) / (2 * delta);
  }

  const dTdx = useMemo(() => (spline ? splineDerivative(spline, Number(estimateAt)) : NaN), [spline, estimateAt]);

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Thermometer className="w-5 h-5" />
            20.2 Use of Splines to Estimate Heat Transfer
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 md:grid-cols-4 gap-3 items-start min-w-0">
            <Labeled label="Depth (m) vs Temperature (°C)">
              <Textarea value={pairs} onChange={(e) => setPairs(e.target.value)} className="bg-zinc-800 text-white h-36 w-full min-w-0" />
            </Labeled>

            <div>
              <Labeled label="Estimate temperature at depth">
                <Input value={estimateAt} onChange={(e) => setEstimateAt(e.target.value)} className="bg-zinc-800 text-white w-28" />
              </Labeled>
              <div className="mt-2 text-sm text-zinc-100">Estimated T: <span style={{ color: theme.accent2 }}>{fmt(estimateVal, 4)} °C</span></div>
            </div>

            <Summary items={[
              { label: "Nodes", value: x.length },
              { label: "Spline type", value: "Natural cubic" },
              { label: "dT/dx at query", value: fmt(dTdx, 6) }
            ]} />

            <div>
             <Equation>{`For heat conduction normal to a surface, heat flux $q(x) \\sim -k \\frac{dT}{dx}$. Use splines to estimate $\\frac{dT}{dx}$ robustly from noisy measurements.`}</Equation>

            </div>
          </div>

          <div className="grid lg:grid-cols-2 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 min-w-0 overflow-hidden">
              <Plot
                data={[
                  { x: x, y: y, mode: "markers", type: "scatter", name: "measurements", marker: { size: 7 } },
                  { x: lineX, y: splineY, mode: "lines", type: "scatter", name: "natural cubic spline", line: { color: theme.accent } },
                  { x: [Number(estimateAt)], y: [estimateVal], mode: "markers", type: "scatter", name: "estimate", marker: { size: 10, color: theme.accent2 } }
                ]}
                layout={{
                  title: "Temperature profile and spline interpolation",
                  xaxis: { title: "Depth (m)" },
                  yaxis: { title: "Temperature (°C)" },
                  autosize: true,
                  paper_bgcolor: "rgba(0,0,0,0)",
                  plot_bgcolor: "rgba(0,0,0,0)",
                  font: { color: "#e5e7eb" },
                  margin: { l: 50, r: 10, b: 40, t: 40 },
                }}
                useResizeHandler
                style={{ width: "100%", height: "min(55vh,520px)" }}
              />
            </div>

            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-3 min-w-0 overflow-auto">
              <div className="text-sm text-zinc-300 mb-2">Interpretation</div>
              <div className="text-zinc-100 text-xs space-y-2">
                <div>Estimated gradient dT/dx ≈ <span style={{ color: theme.accent2 }}>{fmt(dTdx, 6)}</span></div>
                <Note>Multiply by thermal conductivity k (material property) and apply sign to get heat flux: q = -k dT/dx.</Note>
                <Note>Spline smoothing is local and helps reduce amplification of noise in derivative estimates compared to polynomial fits.</Note>
              </div>
            </div>
          </div>

          <div className="space-y-3">
            <Exercise title="Exercise 20.2.1 — Flux estimation" body="Given thermal conductivity k=0.8 W/m·K, compute heat flux at measured point using spline derivative. Compare with finite-difference derivative." solution="Multiply -k * dT/dx from spline and central finite differences; compare percent difference." />
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// Section 20.3 — Fourier Analysis (signal fitting, harmonics)
// ======================================================================
function Section203() {
  // sample signal with noise
  const defaultSignal = "0 0.12\n1 0.84\n2 0.91\n3 0.12\n4 -0.70\n5 -0.96\n6 -0.23\n7 0.66";
  const [pairs, setPairs] = useState(defaultSignal);
  const [numHarmonics, setNumHarmonics] = useState(3);
  const { x: tRaw, y: sig } = useMemo(() => parsePairs(pairs), [pairs]);
  // assume uniform sampling interval 1 for simplicity when computing DFT
  const N = sig.length;
  const dft = useMemo(() => (N ? dftReal(sig) : { re: [], im: [] }), [sig]);
  const re = dft.re || [];
  const im = dft.im || [];

  // Reconstruct using first M harmonics (low-pass)
  const recon = useMemo(() => {
    if (!N) return [];
    const M = Math.min(Number(numHarmonics), Math.floor(N / 2));
    const re2 = new Array(N).fill(0);
    const im2 = new Array(N).fill(0);
    // keep DC and first M harmonics (positive and negative symmetric)
    for (let k = 0; k <= M; k++) {
      re2[k] = re[k];
      im2[k] = im[k];
      // mirror
      const km = (N - k) % N;
      re2[km] = re[km];
      im2[km] = im[km];
    }
    return idftReconstruct(re2, im2);
  }, [re, im, N, numHarmonics]);

  const time = tRaw.length ? tRaw : Array.from({ length: sig.length }, (_, i) => i);
  const reconTime = time;

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Zap className="w-5 h-5" />
            20.3 Fourier Analysis — Harmonic decomposition
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 md:grid-cols-4 gap-3 items-start min-w-0">
            <Labeled label="Time (n) vs signal value">
              <Textarea value={pairs} onChange={(e) => setPairs(e.target.value)} className="bg-zinc-800 text-white h-36 w-full min-w-0" />
            </Labeled>

            <div>
              <Labeled label="Number of harmonics to keep (low-pass)">
                <Input value={numHarmonics} onChange={(e) => setNumHarmonics(Math.max(0, Math.min(50, Number(e.target.value))))} className="bg-zinc-800 text-white w-28" />
              </Labeled>
              <div className="mt-2 text-sm text-zinc-100">Samples: {sig.length}</div>
            </div>

            <Summary items={[
              { label: "DFT bins", value: N },
              { label: "DC magnitude", value: fmt(Math.hypot(re[0] || 0, im[0] || 0), 6) },
              { label: "Reconstruction error", value: fmt(norm2(vecSub(sig, recon)), 6) }
            ]} />

            <div>
              <Equation>{`Discrete Fourier transform: $$X_k=\\sum_{n=0}^{N-1} x_n e^{-i2\\pi kn/N},\\quad x_n=\\frac{1}{N}\\sum_{k=0}^{N-1} X_k e^{i2\\pi kn/N}.$$`}</Equation>
            </div>
          </div>

          <div className="grid lg:grid-cols-2 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 min-w-0 overflow-hidden">
              <Plot
                data={[
                  { x: time, y: sig, mode: "markers+lines", type: "scatter", name: "original", marker: { size: 6 } },
                  { x: reconTime, y: recon, mode: "lines", type: "scatter", name: "reconstruction", line: { color: theme.accent } },
                ]}
                layout={{
                  title: "Signal and truncated Fourier reconstruction",
                  autosize: true,
                  paper_bgcolor: "rgba(0,0,0,0)",
                  plot_bgcolor: "rgba(0,0,0,0)",
                  font: { color: "#e5e7eb" },
                  margin: { l: 40, r: 10, b: 40, t: 40 },
                }}
                useResizeHandler
                style={{ width: "100%", height: "min(55vh,520px)" }}
              />
            </div>

            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-3 min-w-0 overflow-auto">
              <div className="text-sm text-zinc-300 mb-2">Spectral magnitudes</div>
              <div className="text-zinc-100 text-xs space-y-1">
                {re.map((rv, k) => (<div key={k}>k={k}: |X_k| = {fmt(Math.hypot(re[k] || 0, im[k] || 0), 6)}</div>))}
              </div>
              <Note className="mt-2">Truncating harmonics is equivalent to low-pass filtering. For nonuniform sampling or continuous signals, consider Lomb-Scargle or windowing techniques.</Note>
            </div>
          </div>

          <div className="space-y-3">
            <Exercise title="Exercise 20.3.1 — Windowing & leakage" body="Apply different windows (Hann, Hamming) to the signal before DFT and inspect spectral leakage. Compare main-lobe width & side-lobe attenuation." solution="Windows trade main-lobe width vs side-lobe level; Hann reduces leakage compared to rectangular window." />
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// Section 20.4 — Analysis of Experimental Data (smoothing, outlier handling)
// ======================================================================
function Section204() {
  const defaultData = "0 0.1\n1 0.9\n2 2.0\n3 2.8\n4 3.9\n5 5.0\n6 4.8\n7 6.1\n8 6.9"; // noisy rising experimental set
  const [pairs, setPairs] = useState(defaultData);
  const [window, setWindow] = useState(3);
  const [robust, setRobust] = useState(true);
  const { x, y } = useMemo(() => parsePairs(pairs), [pairs]);

  // moving average smoothing
  const smoothed = useMemo(() => {
    if (!y.length) return [];
    const w = Math.max(1, Math.min(25, Number(window)));
    const out = new Array(y.length).fill(0);
    for (let i = 0; i < y.length; i++) {
      const half = Math.floor(w / 2);
      const lo = Math.max(0, i - half);
      const hi = Math.min(y.length - 1, i + half);
      const vals = y.slice(lo, hi + 1);
      if (robust) {
        // median filter
        const sorted = vals.slice().sort((a, b) => a - b);
        out[i] = sorted[Math.floor(sorted.length / 2)];
      } else {
        const s = vals.reduce((a, b) => a + b, 0);
        out[i] = s / vals.length;
      }
    }
    return out;
  }, [y, window, robust]);

  // polynomial fit (degree 2) for trend
  const deg = 2;
  const X = useMemo(() => polyFeatures(x, deg), [x]);
  const sol = useMemo(() => (X.length && y.length ? leastSquares(X, y) : { x: [], ok: false }), [X, y]);
  const beta = sol.ok ? sol.x : [];
  const trend = (X && beta.length) ? matVec(X, beta) : Array.from({ length: y.length }, () => NaN);

  const residuals = y.length ? vecSub(y, trend) : [];

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Activity className="w-5 h-5" />
            20.4 Analysis of Experimental Data
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 md:grid-cols-4 gap-3 items-start min-w-0">
            <Labeled label="Experimental data (t value per line)">
              <Textarea value={pairs} onChange={(e) => setPairs(e.target.value)} className="bg-zinc-800 text-white h-36 w-full min-w-0" />
            </Labeled>

            <div>
              <Labeled label="Smoother window size (odd)">
                <Input value={window} onChange={(e) => setWindow(Math.max(1, Math.min(25, Number(e.target.value))))} className="bg-zinc-800 text-white w-28" />
              </Labeled>
              <div className="mt-2 flex items-center gap-2">
                <label className="text-sm text-zinc-300 flex items-center gap-2">
                  <input type="checkbox" checked={robust} onChange={(e) => setRobust(e.target.checked)} />
                  Robust (median)
                </label>
              </div>
            </div>

            <Summary items={[
              { label: "Points", value: x.length },
              { label: "Trend degree", value: deg },
              { label: "Residual norm", value: fmt(norm2(residuals), 6) }
            ]} />

            <div>
              <Equation>{`Smoothing and trend extraction. Residuals $r_i = y_i - \\hat{y}_i$ help identify outliers and model misspecification.`}</Equation>

            </div>
          </div>

          <div className="grid lg:grid-cols-2 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 min-w-0 overflow-hidden">
              <Plot
                data={[
                  { x: x, y: y, mode: "markers", type: "scatter", name: "raw data", marker: { size: 7 } },
                  { x: x, y: smoothed, mode: "lines", type: "scatter", name: "smoothed", line: { color: theme.accent } },
                  { x: x, y: trend, mode: "lines", type: "scatter", name: "polynomial trend", line: { dash: "dash", color: theme.accent2 } },
                ]}
                layout={{
                  title: "Experimental data: smoothing & trend",
                  autosize: true,
                  paper_bgcolor: "rgba(0,0,0,0)",
                  plot_bgcolor: "rgba(0,0,0,0)",
                  font: { color: "#e5e7eb" },
                  margin: { l: 40, r: 10, b: 40, t: 40 },
                }}
                useResizeHandler
                style={{ width: "100%", height: "min(55vh,520px)" }}
              />
            </div>

            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-3 min-w-0 overflow-auto">
              <div className="text-sm text-zinc-300 mb-2">Residual analysis</div>
              <div className="text-zinc-100 text-xs">
                {residuals.map((r, i) => (<div key={i}>t={x[i]}: residual = {fmt(r, 6)}</div>))}
              </div>
              <Note className="mt-2">Inspect residuals for patterns (nonrandom structure indicates model misspecification). Use robust regression or transform variables if outliers dominate.</Note>
            </div>
          </div>

          <div className="space-y-3">
            <Exercise title="Exercise 20.4.1 — Outlier removal" body="Apply an iterative outlier removal (points whose residual > 3σ) and refit trend. Compare residual norms before and after." solution="Iteratively remove outliers and refit; robust regression (e.g., Huber loss) is preferable to hard removal." />
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// Documentation & Problems (sidebar + final problems panel)
// ======================================================================
const docs20 = {
  "20.1 Linear Regression": [
    "Linear regression models the relationship between a dependent variable $y$ and independent variable(s) $x_i$ using a linear function of parameters.",
    "Ordinary least squares (OLS) minimizes the sum of squared residuals:",
    "$$\\text{Minimize } S = \\sum_i (y_i - \\hat{y}_i)^2,$$ where $\\hat{y}_i$ is the predicted value.",
    "Transformations (log, square root) allow modeling nonlinear relationships as linear in parameters.",
    "$$\\text{Exponential model: } y = A e^{r t} \\Rightarrow \\ln y = \\ln A + r t,$$ which can be fit using linear regression on $\\ln y$.",
    "Assess model quality using $R^2$, adjusted $R^2$, and residual plots to check assumptions of linearity and homoscedasticity."
  ],
  "20.2 Splines & Heat Transfer": [
    "Splines provide smooth interpolation of data points and allow derivative estimation for physical applications.",
    "Cubic splines minimize bending energy and ensure continuity of first and second derivatives.",
    "Natural splines set second derivatives to zero at endpoints; clamped splines fix first derivatives at endpoints.",
    "$$q = -k \\frac{dT}{dx},$$ the Fourier law of heat conduction, requires accurate derivative estimation from discrete temperature data.",
    "Spline-based derivative estimates help compute fluxes in 1D and multi-dimensional heat transfer problems."
  ],
  "20.3 Fourier Analysis": [
    "Discrete Fourier Transform (DFT) expresses a signal $x_n$ as a sum of harmonics:",
    "$$X_k = \\sum_{n=0}^{N-1} x_n e^{-i 2 \\pi k n / N}, \\quad k = 0, 1, \\dots, N-1.$$",
    "Inverse DFT reconstructs the original signal:",
    "$$x_n = \\frac{1}{N} \\sum_{k=0}^{N-1} X_k e^{i 2 \\pi k n / N}.$$",
    "Truncating high-frequency components acts as a low-pass filter to smooth data.",
    "Windowing (Hanning, Hamming) reduces spectral leakage when analyzing finite signals."
  ],
  "20.4 Experimental Data": [
    "Residuals help assess the goodness-of-fit and reveal outliers or model misspecification:",
    "$$r_i = y_i - \\hat{y}_i.$$",
    "Smoothing techniques (moving average, LOESS) reduce noise to reveal underlying trends.",
    "Robust filters (e.g., Tukey's biweight) mitigate the influence of outliers on trend estimates.",
    "Plot residuals against predictors or fitted values to check assumptions like homoscedasticity and independence."
  ],
  "Problems Section": [
    "1. Fit an exponential population model $y = A e^{r t}$ using linear regression on $\ln y$ and plot residuals.",
    "2. Given temperature measurements at discrete points, construct a cubic spline and compute flux $q = -k dT/dx$ at midpoints.",
    "3. Apply DFT to a noisy signal, filter high-frequency components, and reconstruct the smoothed signal using inverse DFT.",
    "4. Perform robust regression on a dataset with outliers and compare residuals with ordinary least squares.",
    "5. Analyze experimental data for trend, smooth using LOESS, and visualize residuals to detect anomalies."
  ],
};


function ChapterDocs20() {
  const [open, setOpen] = useState("20.1 Linear Regression");
  return (
    <div className="w-full mt-6 bg-zinc-900 border border-zinc-700 rounded-xl shadow-lg overflow-hidden">
      <div className="p-3 border-b border-zinc-700 flex items-center justify-between">
        <div className="flex items-center gap-3">
          <BookOpen className="w-5 h-5 text-cyan-400" />
          <div>
            <div className="text-zinc-100 font-semibold">Chapter 20 — Documentation</div>
            <div className="text-zinc-400 text-xs">Case studies & practical formulas</div>
          </div>
        </div>
        <div className="text-zinc-400 text-xs">Applied curve fitting</div>
      </div>

      <div className="grid grid-cols-1 md:grid-cols-4">
        <div className="col-span-1 border-r border-zinc-800">
          {Object.keys(docs20).map((k) => (
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
            {docs20[open].map((p, i) => (<div key={i}><ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>{p}</ReactMarkdown></div>))}
            <div className="mt-3 text-zinc-400 text-xs">References: Applied regression textbooks, standard heat transfer texts, signal processing references.</div>
          </div>
        </div>
      </div>
    </div>
  );
}

function Problems20() {
  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2"><Grid2X2 className="w-5 h-5" /> Problems</CardTitle>
        </CardHeader>
        <CardContent className="space-y-3">
          <Exercise title="P20.1 — Exponential population fit" body="Fit exponential model to provided population data. Provide confidence on fitted growth rate estimated via transformed residuals." solution="Transform to linear using ln(y) and fit with OLS. Estimate residual variance on original scale via delta method or bootstrap." />
          <Exercise title="P20.2 — Spline flux estimation" body="Given temperature profile along depth, estimate heat flux at multiple depths and compare with finite-difference estimates." solution="Spline derivatives provide smoother flux estimates and are less sensitive to noise." />
          <Exercise title="P20.3 — Fourier denoising" body="Use truncated DFT to denoise a sampled noisy signal. Compare RMSE for different cutoff frequencies." solution="Select cutoff based on energy (retain X% of total power) or via cross-validation on held-out samples." />
          <Exercise title="P20.4 — Robust experimental fitting" body="Implement Huber regression to fit trend to noisy experimental data and compare with OLS and median smoothing." solution="Huber mixes L2 and L1 behaviors; use IRLS or gradient-based optimization to estimate robust coefficients." />
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// Page assembly (responsive container)
// ======================================================================
export default function Chapter20() {
  return (
    <div className={`p-4 md:p-6 lg:p-8 space-y-6 ${theme.bg} min-h-screen`}>
      <div className="max-w-screen-xl mx-auto">
        <motion.div initial={{ opacity: 0, y: -10 }} animate={{ opacity: 1, y: 0 }}>
          <h1 className="text-2xl md:text-3xl lg:text-4xl font-bold" style={{ color: theme.accent2 }}>Case Studies: Curve Fitting</h1>
          <div className="text-zinc-400">Applied regression, splines, Fourier methods, and experimental data analysis — interactive case studies.</div>
        </motion.div>

        <div className="space-y-6 mt-4">
          <Section201 />
          <Section202 />
          <Section203 />
          <Section204 />
          <ChapterDocs20 />
          <Problems20 />
          <BottomBar/>
        </div>
      </div>
    </div>
  );
}
