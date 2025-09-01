// src/pages/Chapter18.jsx
// ======================================================================
// Chapter 18 — Interpolation (Newton, Lagrange, Splines, Multidimensional)
// Design: professional, responsive, KaTeX-friendly, Plotly visualizations
// Consistent with Chapter17/12 look-and-feel; single-file React page
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
  Grid,
  Layers,
  Ruler,
  Hash,
  Grid2X2,
  GitBranch,
} from "lucide-react";
import BottomBar from "../components/BottomBar";

// ======================================================================
// Theme & small helpers
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
// Numeric utilities (self-contained)
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

function solveTridiagonal(a, b, c, d) {
  // a: sub-diagonal (n-1), b: diag (n), c: super-diag (n-1), d: rhs (n)
  // returns x (n)
  const n = b.length;
  const cp = c.slice();
  const dp = d.slice();
  // forward
  for (let i = 1; i < n; i++) {
    const m = a[i - 1] / b[i - 1];
    b[i] = b[i] - m * cp[i - 1];
    dp[i] = dp[i] - m * dp[i - 1];
  }
  // back substitution
  const x = Array(n).fill(0);
  x[n - 1] = dp[n - 1] / b[n - 1];
  for (let i = n - 2; i >= 0; i--) {
    x[i] = (dp[i] - cp[i] * x[i + 1]) / b[i];
  }
  return x;
}

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

// ======================================================================
// Interpolation algorithms
// ======================================================================

/* Newton's divided differences:
   Given nodes x0..xn and values y0..yn,
   compute table of divided differences and evaluate Newton form.
*/
function dividedDifferences(xs, ys) {
  const n = xs.length;
  const table = Array.from({ length: n }, (_, i) => Array.from({ length: n - i }, () => 0));
  for (let i = 0; i < n; i++) table[0][i] = ys[i]; // level 0 stored in table[0][i]
  // We'll instead build a 2D triangular table differently: d[i][j] = f[x_i ... x_{i+j}]
  const d = Array.from({ length: n }, () => Array.from({ length: n }, () => 0));
  for (let i = 0; i < n; i++) d[i][0] = ys[i];
  for (let j = 1; j < n; j++) {
    for (let i = 0; i < n - j; i++) {
      const denom = xs[i + j] - xs[i];
      d[i][j] = denom === 0 ? NaN : (d[i + 1][j - 1] - d[i][j - 1]) / denom;
    }
  }
  // coefficients are d[0][0], d[0][1], ..., d[0][n-1]
  const coeffs = new Array(n).fill(0);
  for (let j = 0; j < n; j++) coeffs[j] = d[0][j];
  return { table: d, coeffs };
}

function evalNewton(xs, coeffs, nodes, x) {
  // Horner-like evaluation in Newton basis
  let n = coeffs.length;
  let val = coeffs[n - 1];
  for (let k = n - 2; k >= 0; k--) {
    val = val * (x - nodes[k]) + coeffs[k];
  }
  return val;
}

/* Lagrange interpolation:
   L(x) = sum_{i} y_i * l_i(x)
   where l_i(x) = product_{j != i} (x - x_j)/(x_i - x_j)
*/
function evalLagrange(xs, ys, x) {
  const n = xs.length;
  let s = 0;
  for (let i = 0; i < n; i++) {
    let li = 1;
    for (let j = 0; j < n; j++) if (j !== i) li *= (x - xs[j]) / (xs[i] - xs[j]);
    s += ys[i] * li;
  }
  return s;
}

// Vandermonde coefficients (explicit polynomial coefficients)
function polyCoeffsFromVandermonde(xs, ys) {
  const n = xs.length;
  // Build Vandermonde V(ij) = xs[i]^j (0..n-1)
  const V = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => Math.pow(xs[i], j)));
  // Solve V c = y for polynomial coeffs c_0..c_{n-1}
  const sol = solveGEPP(V, ys);
  return sol;
}

// Natural cubic spline (1D) construction:
// Solve for second derivatives m_i (or spline coefficients) then evaluate piecewise cubic.
function naturalCubicSpline(xs, ys) {
  const n = xs.length;
  if (n < 2) return null;
  const h = new Array(n - 1);
  for (let i = 0; i < n - 1; i++) h[i] = xs[i + 1] - xs[i];
  // Build tridiagonal system for m (second derivatives)
  const a = []; // sub-diagonal
  const b = []; // diag
  const c = []; // super-diagonal
  const d = []; // rhs
  // Natural boundary -> m0 = mn-1 = 0, so we have (n-2)x(n-2) system for interior points
  if (n === 2) {
    // linear between two points, second derivatives 0
    return { xs, ys, m: [0, 0], h };
  }
  for (let i = 1; i < n - 1; i++) {
    a.push(h[i - 1]);
    b.push(2 * (h[i - 1] + h[i]));
    c.push(h[i]);
    d.push(6 * ((ys[i + 1] - ys[i]) / h[i] - (ys[i] - ys[i - 1]) / h[i - 1]));
  }
  // Solve tridiagonal
  const mInterior = solveTridiagonal(a.slice(), b.slice(), c.slice(), d.slice());
  const m = [0, ...mInterior, 0];
  return { xs, ys, m, h };
}

function evalNaturalCubic(spline, x) {
  const { xs, ys, m, h } = spline;
  const n = xs.length;
  if (x <= xs[0]) return ys[0];
  if (x >= xs[n - 1]) return ys[n - 1];
  // find interval i such that xs[i] <= x < xs[i+1]
  let i = 0;
  for (let k = 0; k < n - 1; k++) {
    if (x >= xs[k] && x <= xs[k + 1]) {
      i = k;
      break;
    }
  }
  const xi = xs[i];
  const xi1 = xs[i + 1];
  const hi = h[i];
  const Ai = (xi1 - x) / hi;
  const Bi = (x - xi) / hi;
  const Si = Ai * ys[i] + Bi * ys[i + 1] + ((Ai * (Ai * Ai - 1) * m[i] + Bi * (Bi * Bi - 1) * m[i + 1]) * (hi * hi)) / 6;
  return Si;
}

// Bilinear interpolation on regular grid
function bilinearInterpolate(xGrid, yGrid, zGrid, x, y) {
  // xGrid: sorted unique x nodes (length nx)
  // yGrid: sorted unique y nodes (length ny)
  // zGrid: ny x nx matrix (rows correspond to y)
  const nx = xGrid.length;
  const ny = yGrid.length;
  // find ix, iy such that xGrid[ix] <= x <= xGrid[ix+1]
  let ix = 0;
  for (let i = 0; i < nx - 1; i++) {
    if (x >= xGrid[i] && x <= xGrid[i + 1]) {
      ix = i; break;
    }
  }
  let iy = 0;
  for (let j = 0; j < ny - 1; j++) {
    if (y >= yGrid[j] && y <= yGrid[j + 1]) {
      iy = j; break;
    }
  }
  const x1 = xGrid[ix], x2 = xGrid[ix + 1];
  const y1 = yGrid[iy], y2 = yGrid[iy + 1];
  const Q11 = zGrid[iy][ix], Q21 = zGrid[iy][ix + 1], Q12 = zGrid[iy + 1][ix], Q22 = zGrid[iy + 1][ix + 1];
  const denom = (x2 - x1) * (y2 - y1);
  if (denom === 0) return NaN;
  const val = (Q11 * (x2 - x) * (y2 - y) + Q21 * (x - x1) * (y2 - y) + Q12 * (x2 - x) * (y - y1) + Q22 * (x - x1) * (y - y1)) / denom;
  return val;
}

// ======================================================================
// Small UI building blocks (Equation, Note, Summary, Exercise)
// ======================================================================
function Equation({ children, color = theme.accent }) {
  return (
    <div className="my-2">
      <div className="px-3 py-2 overflow-auto rounded-lg border text-gray-300 border-zinc-700 bg-zinc-900/80">
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
        <Button variant="ghost" className="cursor-pointer bg-white hover:bg-gray-200" size="sm" onClick={() => setOpen((o) => !o)}>
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
// Section 18.1 — Newton's Divided-Difference Interpolating Polynomials
// ======================================================================
function Section181() {
  const defaultPairs = "0 1\n1 2\n2 0\n3 5";
  const [pairs, setPairs] = useState(defaultPairs);
  const { x, y } = useMemo(() => parsePairs(pairs), [pairs]);

  const dd = useMemo(() => (x.length ? dividedDifferences(x, y) : { coeffs: [], table: [] }), [x, y]);
  const coeffs = dd.coeffs;
  const lineX = x.length ? linspace(Math.min(...x) - 0.5, Math.max(...x) + 0.5, 300) : [];
  const newtonY = lineX.map((xx) => (coeffs.length ? evalNewton(x, coeffs, x, xx) : NaN));
  // compute Lagrange for comparison
  const lagrangeY = lineX.map((xx) => (x.length ? evalLagrange(x, y, xx) : NaN));

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Hash className="w-5 h-5" />
            18.1 Newton's Divided-Difference Interpolating Polynomials
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 md:grid-cols-3 gap-3 items-start min-w-0">
            <Labeled label="Data points (x y per line)">
              <Textarea value={pairs} onChange={(e) => setPairs(e.target.value)} className="bg-zinc-800 text-white h-36 w-full min-w-0" />
            </Labeled>

            <Summary items={[
              { label: "Nodes", value: x.length },
              { label: "Newton coeffs", value: coeffs.length },
              { label: "Degree", value: Math.max(0, coeffs.length - 1) }
            ]} />

            <div>
              <div className="text-sm text-zinc-300 mb-2">Newton form</div>
              <Equation>{`Newton's form:\n\n$$P_n(x)=c_0 + c_1(x-x_0) + c_2(x-x_0)(x-x_1)+\\dots + c_n\\prod_{j=0}^{n-1}(x-x_j)$$\n\nCoefficients: $c_k = f[x_0,\\dots,x_k]$ (divided differences).`}</Equation>
            </div>
          </div>

          <div className="grid grid-cols-1 lg:grid-cols-2 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 min-w-0 overflow-hidden">
              <Plot
                data={[
                  { x: x, y: y, mode: "markers", type: "scatter", name: "nodes", marker: { size: 8 } },
                  { x: lineX, y: newtonY, mode: "lines", type: "scatter", name: "Newton", line: { color: theme.accent } },
                  { x: lineX, y: lagrangeY, mode: "lines", type: "scatter", name: "Lagrange", line: { dash: "dash", color: theme.accent2 } },
                ]}
                layout={{
                  title: "Newton vs Lagrange interpolation",
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
              <div className="text-sm text-zinc-300 mb-2">Divided differences (top row = coefficients)</div>
              <div className="text-zinc-100 text-xs">
                {coeffs.map((c, i) => (<div key={i}>c_{i} = <span style={{ color: theme.accent2 }}>{fmt(c, 8)}</span></div>))}
              </div>

              <div className="mt-3">
               <Equation>{`Divided difference definition:

$$
\\begin{aligned}
f[x_i] &= f(x_i), \\\\
f[x_i, x_{i+1}] &= \\frac{f[x_{i+1}] - f[x_i]}{x_{i+1} - x_i}, \\\\
f[x_i, \\dots, x_{i+k}] &= \\frac{f[x_{i+1}, \\dots, x_{i+k}] - f[x_i, \\dots, x_{i+k-1}]}{x_{i+k} - x_i}.
\\end{aligned}
$$
`}</Equation>

              </div>
            </div>
          </div>

          <div className="space-y-3">
            <Exercise
              title="Exercise 18.1.1 — Stability of Newton form"
              body="Change the node ordering (permutations of x) and observe how the Newton coefficients and evaluation behave. Which ordering gives better numerical stability?"
              solution="Often ordering nodes near the evaluation region first improves stability. Chebyshev-like spacing reduces oscillation; using barycentric/Lagrange forms can be more stable in practice."
            />
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// Section 18.2 — Lagrange Interpolating Polynomials
// ======================================================================
function Section182() {
  const defaultPairs = " -2 4\n -1 1\n 0 0\n 1 1\n 2 4";
  const [pairs, setPairs] = useState(defaultPairs);
  const { x, y } = useMemo(() => parsePairs(pairs), [pairs]);

  const lineX = x.length ? linspace(Math.min(...x) - 0.5, Math.max(...x) + 0.5, 300) : [];
  const lagY = lineX.map((xx) => (x.length ? evalLagrange(x, y, xx) : NaN));

  // Precompute Lagrange basis polynomials at nodes for display (symbolic-like)
  const basisDescriptors = useMemo(() => {
    if (!x.length) return [];
    const n = x.length;
    const basis = [];
    for (let i = 0; i < n; i++) {
      const denomParts = [];
      for (let j = 0; j < n; j++) if (j !== i) denomParts.push(`(${x[i]} - ${x[j]})`);
      basis.push({ i, denom: denomParts.join("*") });
    }
    return basis;
  }, [x]);

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Layers className="w-5 h-5" />
            18.2 Lagrange Interpolating Polynomials
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 md:grid-cols-3 gap-3 items-start min-w-0">
            <Labeled label="Data (x y per line)">
              <Textarea value={pairs} onChange={(e) => setPairs(e.target.value)} className="bg-zinc-800 text-white h-36 w-full min-w-0" />
            </Labeled>

            <Summary items={[
              { label: "Nodes", value: x.length },
              { label: "Degree", value: Math.max(0, x.length - 1) },
              { label: "Formula", value: "Lagrange basis" }
            ]} />

            <div>
              <Equation>{`Lagrange form:\n\n$$P_n(x)=\\sum_{i=0}^n y_i\\,\\ell_i(x),\\\\ \\ell_i(x)=\\prod_{j\\ne i}\\frac{x-x_j}{x_i-x_j}.$$`}</Equation>
            </div>
          </div>

          <div className="grid grid-cols-1 lg:grid-cols-2 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 min-w-0 overflow-hidden">
              <Plot
                data={[
                  { x: x, y: y, mode: "markers", type: "scatter", name: "nodes", marker: { size: 8 } },
                  { x: lineX, y: lagY, mode: "lines", type: "scatter", name: "Lagrange", line: { color: theme.accent } },
                ]}
                layout={{
                  title: "Lagrange interpolation curve",
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
              <div className="text-sm text-zinc-300 mb-2">Lagrange basis details</div>
              <div className="text-zinc-100 text-xs space-y-1">
                {basisDescriptors.map((b) => (
                  <div key={b.i}>ℓ_{b.i}(x) denominator: <span style={{ color: theme.accent2 }}>{b.denom}</span></div>
                ))}
              </div>

              <div className="mt-3">
                <Note>Direct evaluation of Lagrange basis is O(n^2) per x; use barycentric form for better performance and stability.</Note>
                <Equation>{`Barycentric form (first form):\n\n$$P_n(x)=\\frac{\\sum_{j=0}^n\\frac{w_j}{x-x_j}y_j}{\\sum_{j=0}^n\\frac{w_j}{x-x_j}},\\\\ w_j=\\frac{1}{\\prod_{k\\ne j}(x_j-x_k)}.$$`}</Equation>
              </div>
            </div>
          </div>

          <div className="space-y-3">
            <Exercise
              title="Exercise 18.2.1 — Barycentric advantage"
              body="Implement barycentric interpolation (weights w_j) and compare evaluation speed & accuracy with naive Lagrange for n=50 nodes."
              solution="Barycentric reduces repeated products; numerically stable and allows O(n) evaluation with precomputed weights. Use Chebyshev nodes to reduce Runge oscillations."
            />
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// Section 18.3 — Coefficients of an Interpolating Polynomial (Vandermonde)
// ======================================================================
function Section183() {
  const defaultPairs = "0 1\n1 2\n2 0\n3 5";
  const [pairs, setPairs] = useState(defaultPairs);
  const { x, y } = useMemo(() => parsePairs(pairs), [pairs]);

  const vand = useMemo(() => {
    if (!x.length) return { coeffs: [], ok: false };
    const sol = polyCoeffsFromVandermonde(x, y);
    return sol;
  }, [x, y]);

  const polyCoeffs = vand.ok ? vand.x : [];

  // Evaluate polynomial (power basis)
  function evalPoly(coeffs, xv) {
    let s = 0;
    for (let j = coeffs.length - 1; j >= 0; j--) s = s * xv + coeffs[j];
    return s;
  }

  const lineX = x.length ? linspace(Math.min(...x) - 0.5, Math.max(...x) + 0.5, 300) : [];
  const polyY = lineX.map((xx) => (polyCoeffs.length ? evalPoly(polyCoeffs, xx) : NaN));

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Ruler className="w-5 h-5" />
            18.3 Coefficients of an Interpolating Polynomial
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 md:grid-cols-3 gap-3 items-start min-w-0">
            <Labeled label="Data (x y per line)">
              <Textarea value={pairs} onChange={(e) => setPairs(e.target.value)} className="bg-zinc-800 text-white h-36 w-full min-w-0" />
            </Labeled>

            <Summary items={[
              { label: "Nodes", value: x.length },
              { label: "Power-basis coeffs", value: polyCoeffs.length },
              { label: "Condition notes", value: "Vandermonde often ill-conditioned" }
            ]} />

            <div>
              <Equation>{`Vandermonde approach: solve \\(V c = y\\) where V_{ij}=x_i^{j}\\). Power-basis coefficients can be large and unstable for high degree.`}</Equation>
            </div>
          </div>

          <div className="grid lg:grid-cols-2 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 min-w-0 overflow-hidden">
              <Plot
                data={[
                  { x: x, y: y, mode: "markers", type: "scatter", name: "nodes", marker: { size: 8 } },
                  { x: lineX, y: polyY, mode: "lines", type: "scatter", name: "power-basis poly", line: { color: theme.accent } },
                ]}
                layout={{
                  title: "Interpolating polynomial (power basis)",
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
              <div className="text-sm text-zinc-300 mb-2">Power-basis coefficients (c_0 + c_1 x + ...)</div>
              <div className="text-zinc-100 text-xs">
                {polyCoeffs.map((c, i) => (<div key={i}>c_{i} = <span style={{ color: theme.accent2 }}>{fmt(c, 8)}</span></div>))}
              </div>
              <div className="mt-3">
                <Note>Because the Vandermonde matrix is often ill-conditioned, prefer Newton/barycentric forms or orthogonal bases (Chebyshev) for stability.</Note>
              </div>
            </div>
          </div>

          <div className="space-y-3">
            <Exercise
              title="Exercise 18.3.1 — Vandermonde Ill-conditioning"
              body="Construct equally spaced nodes x_i in [-1,1] for n=12 and examine condition of Vandermonde and resulting coefficients' magnitude. Compare with Chebyshev nodes."
              solution="Vandermonde with equally spaced nodes will be severely ill-conditioned; coefficients explode. Chebyshev spacing alleviates ill-conditioning and Runge effect."
            />
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// Section 18.4 — Inverse Interpolation
// ======================================================================
function Section184() {
  const defaultPairs = "0 1\n1 3\n2 7\n3 13";
  const [pairs, setPairs] = useState(defaultPairs);
  const [targetY, setTargetY] = useState(7.0);
  const { x, y } = useMemo(() => parsePairs(pairs), [pairs]);

  // For inverse interpolation, we interpolate x as a function of y (swap roles)
  const invX = useMemo(() => {
    if (!y.length) return [];
    // Build Newton on (y, x)
    const dd = dividedDifferences(y, x);
    return { coeffs: dd.coeffs, nodes: y };
  }, [x, y]);

  const invEstimate = useMemo(() => {
    if (!invX.coeffs) return NaN;
    return evalNewton(invX.nodes, invX.coeffs, invX.nodes, Number(targetY));
  }, [invX, targetY]);

  // Also do simple root-finding for f(x) - targetY using interpolant P(x) - use Newton interpolant from 18.1
  const polyFromXY = useMemo(() => {
    if (!x.length) return { eval: () => NaN };
    const dd = dividedDifferences(x, y);
    const coeffs = dd.coeffs;
    return { eval: (xx) => evalNewton(x, coeffs, x, xx) };
  }, [x, y]);

  // Simple bracketed search on P(x) - targetY (bisection-like)
  function findRootOnInterval(a, b, tol = 1e-6, maxIt = 60) {
    let fa = polyFromXY.eval(a) - targetY;
    let fb = polyFromXY.eval(b) - targetY;
    if (isNaN(fa) || isNaN(fb)) return NaN;
    if (fa * fb > 0) {
      // try sampling
      const N = 200;
      let lastX = a, lastF = fa;
      const xs = linspace(a, b, N);
      for (let xi of xs) {
        const fi = polyFromXY.eval(xi) - targetY;
        if (lastF * fi <= 0) {
          a = lastX; b = xi; fa = lastF; fb = fi; break;
        }
        lastX = xi; lastF = fi;
      }
      if (fa * fb > 0) return NaN;
    }
    let lo = a, hi = b, flo = fa, fhi = fb;
    for (let it = 0; it < maxIt; it++) {
      const mid = 0.5 * (lo + hi);
      const fm = polyFromXY.eval(mid) - targetY;
      if (Math.abs(fm) < tol) return mid;
      if (flo * fm <= 0) { hi = mid; fhi = fm; } else { lo = mid; flo = fm; }
    }
    return 0.5 * (lo + hi);
  }

  const rootEstimate = useMemo(() => {
    if (!x.length) return NaN;
    const a = Math.min(...x) - 0.5, b = Math.max(...x) + 0.5;
    return findRootOnInterval(a, b);
  }, [x, y, targetY, polyFromXY]);

  const lineX = x.length ? linspace(Math.min(...x) - 0.5, Math.max(...x) + 0.5, 300) : [];
  const polyY = lineX.map((xx) => (x.length ? polyFromXY.eval(xx) : NaN));

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <GitBranch className="w-5 h-5" />
            18.4 Inverse Interpolation
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 md:grid-cols-4 gap-3 items-start min-w-0">
            <Labeled label="Data (x y per line)">
              <Textarea value={pairs} onChange={(e) => setPairs(e.target.value)} className="bg-zinc-800 text-white h-36 w-full min-w-0" />
            </Labeled>

            <Labeled label="Target y for inverse interp">
              <Input value={targetY} onChange={(e) => setTargetY(Number(e.target.value))} className="bg-zinc-800 text-white" />
            </Labeled>

            <Summary items={[
              { label: "Newton inv. estimate x", value: fmt(invEstimate, 6) },
              { label: "Root (P(x)=y) estimate", value: fmt(rootEstimate, 6) },
              { label: "Method notes", value: "Swap x,y for Newton or solve P(x)-y=0" }
            ]} />

            <div>
              <Equation>{`Inverse interpolation: if y_i = f(x_i) is monotone around target y, interpolate x as a function of y (swap roles). Alternatively, solve $$P(x)-y_{\\text{target}}=0$$ for x using interpolation-based root-finding.`}</Equation>
            </div>
          </div>

          <div className="grid lg:grid-cols-2 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 min-w-0 overflow-hidden">
              <Plot
                data={[
                  { x: lineX, y: polyY, mode: "lines", type: "scatter", name: "Interpolant P(x)", line: { color: theme.accent } },
                  { x: [rootEstimate], y: [targetY], mode: "markers", type: "scatter", name: "root", marker: { size: 10, color: theme.accent2 } },
                  { x: x, y: y, mode: "markers", type: "scatter", name: "nodes", marker: { size: 8 } },
                ]}
                layout={{
                  title: "Interpolant and inverse target",
                  xaxis: { title: "x" },
                  yaxis: { title: "y" },
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
              <div className="text-sm text-zinc-300 mb-2">Inverse interpolation notes</div>
              <div className="text-zinc-100 text-xs">
                Newton inverse (interpolating x as function of y) works when mapping is locally invertible. If y is not monotone near target, solving P(x)-y=0 on an interval is safer.
              </div>
            </div>
          </div>

          <div className="space-y-3">
            <Exercise
              title="Exercise 18.4.1 — Inverse interpolation test"
              body="Generate a monotone quadratic dataset and compare inverse Newton vs root-finding for several target y values. Report errors vs true inverse."
              solution="Inverse Newton works well for monotone data with small degree; root-finding on polynomial interpolant is more general but needs bracket."
            />
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// Section 18.5 — Additional Comments (Runge phenomenon, node selection)
// ======================================================================
function Section185() {
  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Grid2X2 className="w-5 h-5" />
            18.5 Additional Comments: Runge, Node Selection, Barycentric
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-4">
          <Equation>{`**Runge phenomenon:** High-degree interpolation on equally spaced nodes can oscillate strongly near boundaries. Use Chebyshev nodes or piecewise (spline) interpolation to mitigate.`}</Equation>
          
          <Equation>{`**Chebyshev nodes (first kind) on [-1,1]:**\n\n$$x_j = \\cos\\left( \\frac{(2j+1)\\pi}{2(n+1)} \\right), \\quad j = 0, \\dots, n$$`}</Equation>
          
          <Note>
            Practical guidance: for smooth globally-behaved functions, Chebyshev nodes + barycentric
            interpolation or orthogonal polynomial bases perform well. For local features / large n,
            splines are preferred.
          </Note>

          <div className="space-y-3">
            <Exercise
              title="Exercise 18.5.1 — Runge demo"
              body="Interpolate f(x)=1/(1+25 x^2) on [-1,1] with equally spaced nodes and Chebyshev nodes for n=10, 20. Plot and compare maximum error."
              solution="Equally spaced nodes show large boundary oscillations (Runge phenomenon). Chebyshev nodes drastically reduce max error."
            />
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}


// ======================================================================
// Section 18.6 — Spline Interpolation (natural cubic and clamped)
// ======================================================================
function Section186() {
  const defaultPairs = "-2 4\n-1 1\n0 0\n1 1\n2 4";
  const [pairs, setPairs] = useState(defaultPairs);
  const [clampedLeft, setClampedLeft] = useState(false);
  const [clampedRight, setClampedRight] = useState(false);
  const [leftDeriv, setLeftDeriv] = useState(0);
  const [rightDeriv, setRightDeriv] = useState(0);
  const { x, y } = useMemo(() => parsePairs(pairs), [pairs]);

  const naturalSpline = useMemo(() => (x.length ? naturalCubicSpline(x, y) : null), [x, y]);

  // For clamped spline (specified end derivatives), we can modify system: we will implement clamped cubic by building full system
  function clampedCubicSpline(xs, ys, d0, dn) {
    const n = xs.length;
    if (n < 2) return null;
    const h = new Array(n - 1);
    for (let i = 0; i < n - 1; i++) h[i] = xs[i + 1] - xs[i];
    const A = zeros(n, n);
    const rhs = Array(n).fill(0);
    // clamped boundary:
    A[0][0] = 2 * h[0];
    A[0][1] = h[0];
    rhs[0] = 6 * ((ys[1] - ys[0]) / h[0] - d0);
    A[n - 1][n - 2] = h[n - 2];
    A[n - 1][n - 1] = 2 * h[n - 2];
    rhs[n - 1] = 6 * (dn - (ys[n - 1] - ys[n - 2]) / h[n - 2]);
    for (let i = 1; i < n - 1; i++) {
      A[i][i - 1] = h[i - 1];
      A[i][i] = 2 * (h[i - 1] + h[i]);
      A[i][i + 1] = h[i];
      rhs[i] = 6 * ((ys[i + 1] - ys[i]) / h[i] - (ys[i] - ys[i - 1]) / h[i - 1]);
    }
    const sol = solveGEPP(A, rhs);
    return { xs, ys, m: sol.ok ? sol.x : null, h };
  }

  const clampedSpline = useMemo(() => {
    if (!clampedLeft && !clampedRight) return null;
    if (!x.length) return null;
    const d0 = clampedLeft ? Number(leftDeriv) : (3 * (y[1] - y[0]) / (x[1] - x[0])); // fallback
    const dn = clampedRight ? Number(rightDeriv) : (3 * (y[y.length - 1] - y[y.length - 2]) / (x[x.length - 1] - x[x.length - 2]));
    return clampedCubicSpline(x, y, d0, dn);
  }, [x, y, clampedLeft, clampedRight, leftDeriv, rightDeriv]);

  const lineX = x.length ? linspace(Math.min(...x), Math.max(...x), 400) : [];
  const naturalY = (naturalSpline && naturalSpline.m) ? lineX.map((xx) => evalNaturalCubic(naturalSpline, xx)) : [];
  const clampedY = (clampedSpline && clampedSpline.m) ? lineX.map((xx) => evalNaturalCubic(clampedSpline, xx)) : [];

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Ruler className="w-5 h-5" />
            18.6 Spline Interpolation (Cubic)
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 md:grid-cols-4 gap-3 items-start min-w-0">
            <Labeled label="Data (x y per line)">
              <Textarea value={pairs} onChange={(e) => setPairs(e.target.value)} className="bg-zinc-800 text-white h-36 w-full min-w-0" />
            </Labeled>

            <div>
              <div className="text-xs text-zinc-400 mb-1">Boundary</div>
              <div className="flex gap-2 items-center">
                <label className="text-sm text-zinc-300 flex items-center gap-2">
                  <input type="checkbox" checked={clampedLeft} onChange={(e) => setClampedLeft(e.target.checked)} />
                  Left clamped
                </label>
                <Input value={leftDeriv} onChange={(e) => setLeftDeriv(e.target.value)} className="bg-zinc-800 text-white w-28" />
              </div>
              <div className="flex gap-2 items-center mt-2">
                <label className="text-sm text-zinc-300 flex items-center gap-2">
                  <input type="checkbox" checked={clampedRight} onChange={(e) => setClampedRight(e.target.checked)} />
                  Right clamped
                </label>
                <Input value={rightDeriv} onChange={(e) => setRightDeriv(e.target.value)} className="bg-zinc-800 text-white w-28" />
              </div>
            </div>

            <Summary items={[
              { label: "Nodes", value: x.length },
              { label: "Natural / Clamped", value: clampedSpline ? "Clamped" : "Natural" },
              { label: "Smoothness", value: "C2 continuity (second derivative continuous)" }
            ]} />

            <div>
              <Equation>{`Natural cubic spline: piecewise cubics S_i(x) on [x_i,x_{i+1}] with S_i''(x_0)=S_n''(x_n)=0 (natural). Solve tridiagonal system for second derivatives m_i.`}</Equation>
            </div>
          </div>

          <div className="grid lg:grid-cols-2 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 min-w-0 overflow-hidden">
              <Plot
                data={[
                  { x: x, y: y, mode: "markers", type: "scatter", name: "nodes", marker: { size: 8 } },
                  { x: lineX, y: naturalY, mode: "lines", type: "scatter", name: "natural spline", line: { color: theme.accent } },
                  ...(clampedSpline && clampedSpline.m ? [{ x: lineX, y: clampedY, mode: "lines", type: "scatter", name: "clamped spline", line: { dash: "dash", color: theme.accent2 } }] : []),
                ]}
                layout={{
                  title: "Spline interpolation (natural vs clamped)",
                  autosize: true,
                  paper_bgcolor: "rgba(0,0,0,0)",
                  plot_bgcolor: "rgba(0,0,0,0)",
                  font: { color: "#e5e7eb" },
                  margin: { l: 40, r: 10, b: 40, t: 40 },
                }}
                useResizeHandler
                style={{ width: "100%", height: "min(60vh,560px)" }}
              />
            </div>

            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-3 min-w-0 overflow-auto">
              <div className="text-sm text-zinc-300 mb-2">Spline details</div>
              <div className="text-zinc-100 text-xs space-y-1">
                <div>Natural spline second derivatives: {naturalSpline && naturalSpline.m ? naturalSpline.m.map((m,i) => (<div key={i}>m_{i} = {fmt(m,6)}</div>)) : "—"}</div>
                {clampedSpline && clampedSpline.m && (<div className="mt-2">Clamped spline second derivatives: {clampedSpline.m.map((m,i) => (<div key={i}>m_{i} = {fmt(m,6)}</div>))}</div>)}
              </div>
              <div className="mt-3">
                <Note>Spline interpolation avoids Runge oscillations and is local: changing one data point only affects nearby intervals (local support).</Note>
              </div>
            </div>
          </div>

          <div className="space-y-3">
            <Exercise title="Exercise 18.6.1 — Spline vs high-degree poly" body="Interpolate noisy sample of a smooth function with high-degree polynomial and with cubic spline. Compare L∞ error on a fine grid." solution="Splines typically give lower boundary oscillations and more stable interpolation for moderate node counts." />
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// Section 18.7 — Multidimensional Interpolation (bilinear example)
// ======================================================================
function Section187() {
  // default grid 5x5
  const defaultXs = "0 1 2 3 4";
  const defaultYs = "0 1 2 3 4";
  const defaultZ = `0 1 2 3 4
1 2 3 4 5
2 3 4 5 6
3 4 5 6 7
4 5 6 7 8`;
  const [xGridStr, setXGridStr] = useState(defaultXs);
  const [yGridStr, setYGridStr] = useState(defaultYs);
  const [zGridStr, setZGridStr] = useState(defaultZ);
  const [queryX, setQueryX] = useState(1.3);
  const [queryY, setQueryY] = useState(2.7);

  const xGrid = useMemo(() => parsePairs(xGridStr + "\n0 0").x.concat ? xGridStr.split(/[,\s]+/).filter(Boolean).map(Number) : xGridStr.split(/[,\s]+/).filter(Boolean).map(Number), [xGridStr]);
  const yGrid = useMemo(() => yGridStr.split(/[,\s]+/).filter(Boolean).map(Number), [yGridStr]);

  const zGrid = useMemo(() => {
    const rows = zGridStr.trim().split(/[\n;]+/g).map((r) => r.trim().split(/[,\s]+/g).filter(Boolean).map(Number));
    return rows;
  }, [zGridStr]);

  const interpVal = useMemo(() => {
    if (!xGrid.length || !yGrid.length || !zGrid.length) return NaN;
    return bilinearInterpolate(xGrid, yGrid, zGrid, Number(queryX), Number(queryY));
  }, [xGrid, yGrid, zGrid, queryX, queryY]);

  // build heatmap for display
  const heatX = xGrid;
  const heatY = yGrid;
  const zForPlot = zGrid; // assume rows correspond to y

  // contour grid fine for visualization
  const fineX = linspace(Math.min(...xGrid), Math.max(...xGrid), 80);
  const fineY = linspace(Math.min(...yGrid), Math.max(...yGrid), 80);
  const fineZ = fineY.map((yy) => fineX.map((xx) => bilinearInterpolate(xGrid, yGrid, zGrid, xx, yy)));

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Grid className="w-5 h-5" />
            18.7 Multidimensional Interpolation (Bilinear demo)
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 md:grid-cols-4 gap-3 items-start min-w-0">
            <Labeled label="x grid (space-separated)">
              <Input value={xGridStr} onChange={(e) => setXGridStr(e.target.value)} className="bg-zinc-800 text-white" />
            </Labeled>
            <Labeled label="y grid (space-separated)">
              <Input value={yGridStr} onChange={(e) => setYGridStr(e.target.value)} className="bg-zinc-800 text-white" />
            </Labeled>
            <Labeled label="z grid (rows per y line)">
              <Textarea value={zGridStr} onChange={(e) => setZGridStr(e.target.value)} className="bg-zinc-800 text-white h-28 w-full min-w-0" />
            </Labeled>
            <div>
              <div className="text-xs text-zinc-400 mb-1">Query point</div>
              <div className="flex gap-2">
                <Input value={queryX} onChange={(e) => setQueryX(e.target.value)} className="bg-zinc-800 text-white w-28" />
                <Input value={queryY} onChange={(e) => setQueryY(e.target.value)} className="bg-zinc-800 text-white w-28" />
              </div>
              <div className="mt-2 text-sm text-zinc-100">Interpolated value: <span style={{ color: theme.accent2 }}>{fmt(interpVal, 6)}</span></div>
            </div>
          </div>

          <div className="grid lg:grid-cols-2 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 min-w-0 overflow-hidden">
              <Plot
                data={[
                  { z: zForPlot, x: heatX, y: heatY, type: "heatmap", name: "grid", colorscale: "Viridis" },
                  { x: [queryX], y: [queryY], z: [[interpVal]], mode: "markers", type: "scatter", name: "query", marker: { size: 12, color: theme.accent2 } },
                ]}
                layout={{
                  title: "Grid (heatmap) and query point",
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

            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-3 min-w-0 overflow-auto">
              <div className="text-sm text-zinc-300 mb-2">Bilinear interpolation formula</div>
             <Equation>{`Given grid $(x_1,x_2)$ and $(y_1,y_2)$ with values $Q_{11}, Q_{21}, Q_{12}, Q_{22}$:

$$
f(x,y) = \\frac{1}{(x_2 - x_1)(y_2 - y_1)} \\Big(
    Q_{11}(x_2 - x)(y_2 - y) \\\\
    + Q_{21}(x - x_1)(y_2 - y) \\\\
    + Q_{12}(x_2 - x)(y - y_1) \\\\
    + Q_{22}(x - x_1)(y - y_1)
\\Big).
$$`}</Equation>

              <Note>Bilinear is tensor product of 1D linear interpolations. For smoother surfaces use bicubic or radial basis functions for scattered data.</Note>
            </div>
          </div>

          <div className="space-y-3">
           <Exercise
  title="Exercise 18.7.1 — Bilinear accuracy"
  body={`Compare bilinear interpolation vs bicubic on a smooth test function (e.g., $\\sin(x) \\cos(y)$) on a coarse grid. Report RMSE on a fine grid.`}
  solution={`Bicubic will typically give lower error for smooth functions; bilinear is faster and local.`}
/>

          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// Documentation panel (left nav + KaTeX content) and Problems
// ======================================================================
const docs18 = {
  "18.1 Newton's Divided-Difference": [
    "\\textbf{Newton's divided difference} constructs an interpolant incrementally in Newton basis, where each new node adds a higher-order term.",
    "The interpolant has the form $$P_n(x) = f[x_0] + f[x_0,x_1](x-x_0) + f[x_0,x_1,x_2](x-x_0)(x-x_1)+\\cdots+f[x_0,\\ldots,x_n](x-x_0)\\cdots(x-x_{n-1}).$$",
    "Here, the coefficients are \\emph{divided differences}, recursively defined by $$f[x_i,\\ldots,x_{i+k}] = \\frac{f[x_{i+1},\\ldots,x_{i+k}] - f[x_i,\\ldots,x_{i+k-1}]}{x_{i+k}-x_i}.$$",
    "This form allows efficient updating when new nodes are added without recomputing the entire interpolant.",
    "Numerical stability improves if nodes are ordered carefully; Chebyshev nodes mitigate large oscillations."
  ],

  "18.2 Lagrange Interpolating Polynomials": [
    "The \\textbf{Lagrange form} expresses the interpolant as a linear combination of basis polynomials:",
    "$$P_n(x) = \\sum_{i=0}^n f(x_i) \\ell_i(x), \\quad \\ell_i(x) = \\prod_{j=0, j\\neq i}^n \\frac{x-x_j}{x_i-x_j}.$$",
    "The basis \\(\\ell_i(x)\\) satisfies the cardinal property \\(\\ell_i(x_j)=\\delta_{ij}\\).",
    "Direct evaluation is inefficient; instead the \\textbf{barycentric form} is preferred:",
    "$$P_n(x) = \\frac{\\sum_{i=0}^n \\frac{w_i}{x-x_i} f(x_i)}{\\sum_{i=0}^n \\frac{w_i}{x-x_i}}, \\quad w_i = \\frac{1}{\\prod_{j\\neq i}(x_i-x_j)}.$$",
    "Weights \\(w_i\\) can be precomputed, making this representation both numerically stable and efficient."
  ],

  "18.3 Polynomial Coefficients": [
    "Interpolants can be represented in \\textbf{power basis}: $$P_n(x) = a_0 + a_1x + a_2x^2 + \\cdots + a_nx^n.$$",
    "Coefficients \\(a_k\\) are found by solving a \\textbf{Vandermonde system}, which is notoriously ill-conditioned for large \\(n\\).",
    "For large degree, the errors grow rapidly due to sensitivity of the Vandermonde matrix.",
    "To mitigate instability, one prefers orthogonal bases such as Chebyshev or Legendre polynomials.",
    "In practice, the Newton or barycentric form is far superior for stability and computational cost."
  ],

  "18.4 Inverse Interpolation": [
    "\\textbf{Inverse interpolation} aims to find \\(x\\) given \\(y\\), when monotonicity allows swapping roles of dependent/independent variables.",
    "If \\(P(x)\\) interpolates data \\((x_i,f(x_i))\\), then inverse interpolation builds an interpolant for \\(x\\) as a function of \\(y=f(x)\\).",
    "Alternatively, one can solve the equation $$P(x) - y^* = 0$$ for a target value \\(y^*\\), which reduces to interpolation-based root finding.",
    "Methods such as the secant method can be interpreted as a two-point inverse interpolation.",
    "Accuracy depends strongly on monotonicity of data; otherwise the inverse is multi-valued."
  ],

  "18.5 Additional Comments": [
    "\\textbf{Runge phenomenon:} Equally spaced interpolation on large intervals leads to oscillations near boundaries.",
    "Chebyshev nodes (roots or extrema of Chebyshev polynomials) minimize the maximum error.",
    "For practical use, piecewise interpolation (splines) avoids Runge-type oscillations.",
    "Barycentric interpolation is recommended for efficiency: \\(O(n)\\) evaluation time vs naive \\(O(n^2)\\).",
    "Stability improves when using orthogonal polynomial bases or carefully chosen interpolation nodes."
  ],

  "18.6 Spline Interpolation": [
    "\\textbf{Cubic splines} are piecewise polynomials of degree 3, ensuring \\(C^2\\) continuity across intervals.",
    "The natural spline boundary condition imposes zero second derivative at endpoints.",
    "The clamped spline specifies first derivatives at the boundaries, often from physical slope conditions.",
    "The spline system reduces to solving a tridiagonal linear system for the unknown second derivatives.",
    "Splines achieve small oscillations and good approximation properties compared to high-degree global polynomials."
  ],

  "18.7 Multidimensional Interpolation": [
    "For data on tensor grids, interpolation generalizes to tensor-product forms: bilinear (2D), trilinear (3D), bicubic, etc.",
    "Example: \\textbf{Bilinear interpolation} on rectangle uses a weighted average of 4 corner values.",
    "Bicubic interpolation incorporates derivatives, providing smoother results for smooth data.",
    "For scattered data, popular methods include radial basis functions (RBF), Shepard's method, and Delaunay triangulation (linear on simplices).",
    "Efficiency considerations: kd-trees or nearest-neighbor search can accelerate local interpolation for large data sets."
  ],

  "Problems Section": [
    "Problems include constructing Newton and Lagrange interpolants for small data sets.",
    "Experiments with Runge function $$f(x)=\\frac{1}{1+25x^2}$$ illustrate oscillations on equispaced vs Chebyshev nodes.",
    "Implement barycentric form and compare evaluation stability with direct Lagrange form.",
    "Spline problems: natural vs clamped splines on smooth test functions; examine second derivatives.",
    "Inverse interpolation tasks: find roots of nonlinear functions by constructing inverse interpolants.",
    "Multidimensional tasks: compare bilinear vs bicubic accuracy on test functions; compute RMSE."
  ]
};


function ChapterDocs18() {
  const [open, setOpen] = useState("18.1 Newton's Divided-Difference");
  return (
    <div className="w-full mt-6 bg-zinc-900 border border-zinc-700 rounded-xl shadow-lg overflow-hidden">
      <div className="p-3 border-b border-zinc-700 flex items-center justify-between">
        <div className="flex items-center gap-3">
          <BookOpen className="w-5 h-5 text-cyan-400" />
          <div>
            <div className="text-zinc-100 font-semibold">Chapter 18 — Documentation</div>
            <div className="text-zinc-400 text-xs">Theory, formulas, and notes</div>
          </div>
        </div>
        <div className="text-zinc-400 text-xs">Interpolation basics & advanced topics</div>
      </div>

      <div className="grid grid-cols-1 md:grid-cols-4">
        <div className="col-span-1 border-r border-zinc-800">
          {Object.keys(docs18).map((k) => (
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
            {docs18[open].map((p, i) => (
              <div key={i}><ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>{p}</ReactMarkdown></div>
            ))}
            <div className="mt-3 text-zinc-400 text-xs">References: Press et al. (Numerical Recipes), Trefethen (Approximation Theory), de Boor (Splines), Berrut & Trefethen (Barycentric Lagrange).</div>
          </div>
        </div>
      </div>
    </div>
  );
}

// ======================================================================
// Problems panel
// ======================================================================
function Problems18() {
  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Grid2X2 className="w-5 h-5" />
            Problems
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-3">
          <Exercise title="P18.1 — Implement barycentric interpolation" body="Implement first-and second-form barycentric interpolation and compare with naive Lagrange on 50 random nodes." solution="Barycentric should match and be numerically better; first form directly uses precomputed weights." />
          <Exercise title="P18.2 — Chebyshev vs equally-spaced" body="Interpolate Runge function f(x)=1/(1+25x^2) with n=10,20 equally spaced vs Chebyshev nodes; plot and compute max error." solution="Chebyshev nodes dramatically reduce boundary oscillation and max error." />
          <Exercise title="P18.3 — Spline derivatives" body="Compute cubic spline derivatives at nodes and verify C2 continuity numerically." solution="Check left/right derivatives at internal nodes; they should match within numerical tolerance." />
          <Exercise title="P18.4 — Inverse interpolation application" body="For monotone data, use inverse interpolation to estimate the inverse function at sample y values; compare with direct root-finding." solution="Inverse interpolation is efficient when monotonic; root-finding is more general." />
          <Exercise title="P18.5 — Multidimensional interpolation" body="Compare bilinear vs bicubic vs RBF on scattered samples from f(x,y)=sin(x)cos(y) and report RMSE on a fine grid." solution="RBF/bicubic typically better for smooth data; bilinear is simplest and fastest on grid." />
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// Page assembly
// ======================================================================
export default function Chapter18() {
  return (
    <div className={`p-4 md:p-6 lg:p-8 space-y-6 ${theme.bg} min-h-screen`}>
      <div className="max-w-screen-xl mx-auto">
        <motion.div initial={{ opacity: 0, y: -10 }} animate={{ opacity: 1, y: 0 }}>
          <h1 className="text-2xl md:text-3xl lg:text-4xl font-bold" style={{ color: theme.accent2 }}>Interpolation</h1>
          <div className="text-zinc-400">Newton, Lagrange, coefficients, inverse interpolation, splines, and multidimensional methods.</div>
        </motion.div>

        <div className="space-y-6 mt-4">
          <Section181 />
          <Section182 />
          <Section183 />
          <Section184 />
          <Section185 />
          <Section186 />
          <Section187 />
          <ChapterDocs18 />
          <Problems18 />
          <BottomBar/>
        </div>
      </div>
    </div>
  );
}
