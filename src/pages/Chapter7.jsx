// src/pages/Chapter7.jsx
// ======================================================================
// Chapter 7 — Roots of Polynomials
// Visualization-first page: every subsection has 2D + 3D interactive plots
// UI: dark emerald + neon cyan theme, framer-motion micro-animations
// Math: KaTeX via react-markdown + remark-math + rehype-katex
// Exercises: Embedded per subsection (with toggles & solutions)
// ======================================================================

import React, { useMemo, useState } from "react";
import Plot from "react-plotly.js";
import { motion } from "framer-motion";
import ReactMarkdown from "react-markdown";
import remarkMath from "remark-math";
import rehypeKatex from "rehype-katex";
import "katex/dist/katex.min.css";
import BottomBar from "../components/BottomBar";
import { Card, CardHeader, CardContent, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Button } from "@/components/ui/button";

import {
  BookOpen,
  List,
  ChevronRight,
  ChevronDown,
  Play,
  RefreshCcw,
  Target,
  Shapes,
  FunctionSquare,
  Ruler,
  Braces,
  PlusCircle,
  Calculator,
} from "lucide-react";

// ======================================================================
// Theme helpers
// ======================================================================
const theme = {
  bg: "bg-zinc-950",
  panel: "bg-zinc-900/60 border border-zinc-700",
  text: "text-zinc-200",
  subtext: "text-zinc-400",
  accent: "#22d3ee", // cyan-400
  accent2: "#34d399", // emerald-400
  warn: "#f59e0b", // amber-500
  danger: "#ef4444", // red-500
};

// motion presets
const fadeUp = {
  initial: { opacity: 0, y: 16 },
  animate: { opacity: 1, y: 0 },
  transition: { duration: 0.3 },
};

// ======================================================================
/** Utilities: polynomial parsing, evaluation, derivatives, deflation */
// ======================================================================

/** parseCoeffs: string -> number[]
 * Accepts space/comma separated, highest degree first. Example: "1 -6 11 -6"
 */
function parseCoeffs(s) {
  return s
    .split(/[,\s]+/g)
    .filter(Boolean)
    .map((v) => Number(v));
}

/** horner evaluation */
function horner(coeffs, x) {
  let y = 0;
  for (let i = 0; i < coeffs.length; i++) y = y * x + coeffs[i];
  return y;
}

/** polynomial derivative coefficients */
function polyDeriv(coeffs) {
  const n = coeffs.length - 1;
  if (n <= 0) return [0];
  const out = new Array(n);
  for (let i = 0; i < n; i++) out[i] = coeffs[i] * (n - i);
  return out;
}

/** synthetic division by (x - r) => {q, remainder} */
function syntheticDiv(coeffs, r) {
  const n = coeffs.length;
  const q = new Array(n - 1);
  let b = coeffs[0];
  q[0] = b;
  for (let i = 1; i < n - 1; i++) {
    b = coeffs[i] + b * r;
    q[i] = b;
  }
  const remainder = coeffs[n - 1] + b * r;
  return { q, remainder };
}

/** polynomial evaluate on array xs */
function polyEvalArray(coeffs, xs) {
  return xs.map((x) => horner(coeffs, x));
}

/** safe numeric */
function safeN(n, fallback = 0) {
  const v = Number(n);
  return Number.isFinite(v) ? v : fallback;
}

/** linspace helper */
function linspace(a, b, n) {
  if (n <= 1) return [a];
  const dx = (b - a) / (n - 1);
  const arr = new Array(n);
  for (let i = 0; i < n; i++) arr[i] = a + i * dx;
  return arr;
}

/** meshgrid-like helper for 3D surfaces */
function grid2D(xmin, xmax, xN, ymin, ymax, yN) {
  const xs = linspace(xmin, xmax, xN);
  const ys = linspace(ymin, ymax, yN);
  const Z = new Array(yN);
  for (let j = 0; j < yN; j++) Z[j] = new Array(xN).fill(0);
  return { xs, ys, Z };
}

/** |p(x+iy)| surface */
function surfaceModulus(coeffs, xmin, xmax, xN, ymin, ymax, yN) {
  const xs = linspace(xmin, xmax, xN);
  const ys = linspace(ymin, ymax, yN);
  const Z = ys.map((y) =>
    xs.map((x) => {
      // evaluate complex via Horner: real+imag
      let ar = 0,
        ai = 0;
      for (let i = 0; i < coeffs.length; i++) {
        // multiply (ar+ i ai) * (x+i y) + c
        const rr = ar * x - ai * y;
        const ii = ar * y + ai * x;
        ar = rr + coeffs[i];
        ai = ii;
      }
      return Math.sqrt(ar * ar + ai * ai);
    })
  );
  return { xs, ys, Z };
}

/** Newton for polynomial root (real) with damping */
function newtonPoly(coeffs, x0, maxIt = 30, tol = 1e-10, lambda = 1) {
  const d = polyDeriv(coeffs);
  let x = x0;
  const trace = [];
  for (let k = 0; k < maxIt; k++) {
    const fx = horner(coeffs, x);
    const dfx = horner(d, x);
    trace.push({ k, x, fx, dfx });
    if (!Number.isFinite(fx) || !Number.isFinite(dfx)) break;
    if (Math.abs(fx) < tol) break;
    if (Math.abs(dfx) < 1e-14) break;
    x = x - lambda * fx / dfx;
  }
  return trace;
}

/** Secant for polynomial root (real) */
function secantPoly(coeffs, x0, x1, maxIt = 30, tol = 1e-10) {
  let a = x0,
    b = x1;
  const trace = [];
  let fa = horner(coeffs, a);
  let fb = horner(coeffs, b);
  for (let k = 0; k < maxIt; k++) {
    trace.push({ k, a, b, fa, fb });
    if (Math.abs(fb) < tol) break;
    const denom = fb - fa;
    if (Math.abs(denom) < 1e-14) break;
    const c = b - fb * (b - a) / denom;
    a = b;
    fa = fb;
    b = c;
    fb = horner(coeffs, b);
  }
  return trace;
}

/** Muller's method (real/complex in R2 form) simplified */
function mullerStep(p, x0, x1, x2) {
  const f0 = horner(p, x0);
  const f1 = horner(p, x1);
  const f2 = horner(p, x2);
  const h0 = x1 - x0;
  const h1 = x2 - x1;
  const d0 = (f1 - f0) / h0;
  const d1 = (f2 - f1) / h1;
  const a = (d1 - d0) / (h1 + h0);
  const b = a * h1 + d1;
  const c = f2;
  const disc = b * b - 4 * a * c;
  const denom = Math.abs(b + Math.sign(b) * Math.sqrt(Math.max(0, disc))) > Math.abs(b - Math.sign(b) * Math.sqrt(Math.max(0, disc)))
    ? b + Math.sign(b) * Math.sqrt(Math.max(0, disc))
    : b - Math.sign(b) * Math.sqrt(Math.max(0, disc));
  const dx = denom !== 0 ? (-2 * c) / denom : 0;
  return x2 + dx;
}

/** Bairstow: find r,s (quadratic factor x^2 + r x + s)
 * Returns {r,s, iters, ok}
 */
function bairstow(coeffs, r0 = -1, s0 = -1, maxIt = 50, tol = 1e-10) {
  // coeffs: [a_n, ..., a_0]
  const n = coeffs.length - 1;
  if (n < 2) return { r: 0, s: 0, iters: [], ok: false };

  let r = r0,
    s = s0;
  const iters = [];

  const b = new Array(n + 1).fill(0);
  const c = new Array(n + 1).fill(0);

  for (let k = 0; k < maxIt; k++) {
    // synthetic for quadratic
    b[n] = coeffs[n];
    b[n - 1] = coeffs[n - 1] + r * b[n];

    for (let i = n - 2; i >= 0; i--) {
      b[i] = coeffs[i] + r * b[i + 1] + s * b[i + 2];
    }

    c[n] = b[n];
    c[n - 1] = b[n - 1] + r * c[n];
    for (let i = n - 2; i >= 0; i--) {
      c[i] = b[i + 1] + r * c[i + 1] + s * c[i + 2];
    }

    const F = b[1];
    const G = b[0];

    const J11 = c[2];
    const J12 = c[3];
    const J21 = c[1];
    const J22 = c[2];

    const det = J11 * J22 - J12 * J21;
    if (Math.abs(det) < 1e-16) break;

    const dr = (-G * J22 + F * J12) / det;
    const ds = (-F * J11 + G * J21) / det;

    r += dr;
    s += ds;

    iters.push({ k, r, s, dr, ds, F, G });

    if (Math.max(Math.abs(dr), Math.abs(ds)) < tol && Math.max(Math.abs(F), Math.abs(G)) < tol) {
      return { r, s, iters, ok: true };
    }
  }
  return { r, s, iters, ok: false };
}

/** Quadratic roots given r,s for factor x^2 + r x + s */
function rootsFromRS(r, s) {
  const a = 1, b = r, c = s;
  const D = b * b - 4 * a * c;
  if (D >= 0) {
    const sd = Math.sqrt(D);
    return [(-b + sd) / (2 * a), (-b - sd) / (2 * a)];
  }
  const real = -b / (2 * a);
  const imag = Math.sqrt(-D) / (2 * a);
  return [`${real}+${imag}i`, `${real}-${imag}i`];
}

/** Descartes sign changes for positive axis */
function descartsSignChanges(coeffs) {
  const signs = coeffs.filter((c) => c !== 0).map((c) => (c > 0 ? 1 : -1));
  let changes = 0;
  for (let i = 1; i < signs.length; i++) if (signs[i] !== signs[i - 1]) changes++;
  return changes;
}

/** Cauchy bound: 1 + max |a_i / a_n| for i=0..n-1 */
function cauchyBound(coeffs) {
  const an = coeffs[0];
  const rest = coeffs.slice(1).map((c) => Math.abs(c / an));
  return 1 + Math.max(...rest);
}

// ======================================================================
// Shared small components
// ======================================================================

function Equation({ children, color = theme.accent }) {
  // render KaTeX via Markdown; wrap in a colored chip bar
  return (
    <div className="my-2">
      <div className="px-3 py-2 text-white  rounded-lg border border-zinc-700 bg-zinc-900/80">
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
        <div
          key={i}
          className="rounded-lg bg-zinc-900/70 border border-zinc-700 px-3 py-2 flex flex-col"
        >
          <div className="text-xs text-zinc-400">{it.label}</div>
          <div className="text-zinc-100 text-sm">{it.value}</div>
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
        <Button variant="ghost" size="sm" className="cursor-pointer bg-white hover:bg-gray-200" onClick={() => setOpen((o) => !o)}>
          {open ? <ChevronDown className="w-4 h-4" /> : <ChevronRight className="w-4 h-4" />} Solution
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
// Docs text (KaTeX-ready) with extra details & theory
// ======================================================================
const docsText = {
  "7.1 Polynomials in Engineering and Science": [
    `Polynomials appear in **vibration analysis**, **control characteristic equations**, **signal processing** (filter design), **robot kinematics**, and **computational geometry**.`,
    `**General form**:\n\n$$p(x)=a_n x^n + a_{n-1}x^{n-1}+\\cdots + a_1 x + a_0, \\quad a_n\\neq 0.$$\n\nThe **roots** of \\(p\\) encode system behaviours such as natural frequencies, equilibria, or stability margins.`,
    `**Continuity & Differentiability:** Polynomials are entire functions; derivatives reduce degree: $$p'(x)=n a_n x^{n-1} + \\cdots + a_1.$$`,
    `**Factorization:** Over \\(\\mathbb{C}\\), $$p(x)=a_n\\prod_{k=1}^n (x-r_k).$$ Multiplicities and clustering matter for numerical conditioning.`,
  ],
  "7.2 Computing with Polynomials": [
    `**Horner's scheme** evaluates in \\(\\mathcal{O}(n)\\) operations with good numerical stability:\n\n$$p(x) = (\\cdots((a_n x + a_{n-1})x + a_{n-2})x + \\cdots + a_0).$$`,
    `**Synthetic division** enables fast deflation: dividing by \\(x-r\\) is linear-time.`,
    `**Conditioning:** Coefficient scaling and variable shift \\(x=y+c\\) can improve stability.`,
  ],
  "7.3 Conventional Methods": [
    `Classical root-finding applies **Newton's method** to \\(p\\), often with **deflation** as roots are found.`,
    `Newton update: $$x_{k+1}=x_k - \\frac{p(x_k)}{p'(x_k)}.$$ Quadratic convergence near simple roots.`,
    `**Secant** (derivative-free) and **Aitken-accelerated** schemes are common alternatives.`,
  ],
  "7.4 Müller’s Method": [
    `Müller fits a parabola through \\((x_{k-2},f_{k-2}), (x_{k-1},f_{k-1}), (x_k,f_k)\\) and takes a root:\n\n$$x_{k+1}=x_k-\\frac{2c}{b\\pm\\sqrt{b^2-4ac}},$$\nwhere \\(a,b,c\\) are quadratic coefficients in the local fit. It can converge to **complex roots** without complex arithmetic in the function evaluations.`,
  ],
  "7.5 Bairstow’s Method": [
    `Bairstow iteratively extracts a quadratic factor \\(x^2+rx+s\\): solve \\(F(r,s)=0\\), \\(G(r,s)=0\\) from synthetic division remainders. Works well for real-coefficient polynomials to obtain complex-conjugate pairs.`,
  ],
  "7.6 Other Methods": [
    `**Companion matrix:** roots of \\(p\\) are eigenvalues of \\(C\\). Eigen-solvers are robust and widely used:\n\nFor $$p(x)=x^n+a_{n-1}x^{n-1}+\\dots+a_0,$$ the companion matrix is\n\n$$C=\\begin{bmatrix}0&1&0&\\cdots&0\\\\0&0&1&\\cdots&0\\\\\\vdots&\\vdots&\\vdots&\\ddots&\\vdots\\\\0&0&0&\\cdots&1\\\\-a_0&-a_1&-a_2&\\cdots&-a_{n-1}\\end{bmatrix}.$$`,
    `**Complex-plane visualization:** The surface \\(z\\mapsto |p(z)|\\) highlights valleys at roots.`,
  ],
  "7.7 Root Location Theorems": [
    `**Descartes' Rule of Signs:** Number of positive real roots \\(\\le\\) sign changes in \\((a_n,\\dots,a_0)\\), difference by an even integer.`,
    `**Cauchy Bound:** All roots lie in $$|z|\\le 1+\\max_{0\\le i<n}\\left|\\frac{a_i}{a_n}\\right|.$$`,
    `**Rouché & Gershgorin (matrix view):** Offer alternative complex-plane counts/bounds.`,
  ],
};

// ======================================================================
// 7.1 Panel — Intro visualization (shape explorer)
// ======================================================================
function Section71() {
  const [coeffStr, setCoeffStr] = useState("1 -6 11 -6"); // (x-1)(x-2)(x-3)
  const coeffs = useMemo(() => parseCoeffs(coeffStr), [coeffStr]);

  const X = useMemo(() => linspace(-4, 6, 400), []);
  const Y = useMemo(() => polyEvalArray(coeffs, X), [coeffs, X]);

  // 3D |p(x+iy)| surface (shows minima at roots)
  const surface = useMemo(() => surfaceModulus(coeffs, -4, 4, 50, -4, 4, 50), [coeffs]);

  const posChanges = useMemo(() => descartsSignChanges(coeffs), [coeffs]);
  const cbound = useMemo(() => cauchyBound(coeffs), [coeffs]);

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Shapes className="w-5 h-5" />
            7.1 Polynomials in Engineering & Science — Shape Explorer
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-3">
          <div className="grid md:grid-cols-3 gap-3">
            <Labeled label="Coefficients (highest degree → lowest)">
              <Input
                value={coeffStr}
                onChange={(e) => setCoeffStr(e.target.value)}
                className="bg-zinc-800 text-white"
                placeholder="e.g. 1 -6 11 -6"
              />
            </Labeled>
            <Summary
              items={[
                { label: "Degree", value: Math.max(0, coeffs.length - 1) },
                { label: "Descartes sign changes (+)", value: posChanges },
                { label: "Cauchy bound", value: cbound.toFixed(3) },
              ]}
            />
          </div>

          <div className="grid md:grid-cols-2 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2">
              <Plot
                data={[
                  {
                    x: X,
                    y: Y,
                    type: "scatter",
                    mode: "lines",
                    line: { width: 2, color: theme.accent },
                    name: "p(x)",
                  },
                  {
                    x: X,
                    y: X.map(() => 0),
                    type: "scatter",
                    mode: "lines",
                    line: { width: 1, dash: "dot", color: "#71717a" },
                    name: "x-axis",
                  },
                ]}
                layout={{
                  title: "p(x) (2D)",
                  paper_bgcolor: "rgba(0,0,0,0)",
                  plot_bgcolor: "rgba(0,0,0,0)",
                  font: { color: "#e5e7eb" },
                  margin: { l: 40, r: 10, b: 40, t: 40 },
                }}
                style={{ width: "100%", height: 320 }}
                useResizeHandler
              />
            </div>
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2">
              <Plot
                data={[
                  {
                    x: surface.xs,
                    y: surface.ys,
                    z: surface.Z,
                    type: "surface",
                    colorscale: "Viridis",
                    showscale: false,
                  },
                ]}
                layout={{
                  title: "|p(x+iy)| surface (3D)",
                  paper_bgcolor: "rgba(0,0,0,0)",
                  plot_bgcolor: "rgba(0,0,0,0)",
                  font: { color: "#e5e7eb" },
                  scene: {
                    xaxis: { title: "Re", color: "#e5e7eb" },
                    yaxis: { title: "Im", color: "#e5e7eb" },
                    zaxis: { title: "|p|", color: "#e5e7eb" },
                  },
                  margin: { l: 0, r: 0, b: 0, t: 40 },
                }}
                style={{ width: "100%", height: 320 }}
                useResizeHandler
              />
            </div>
          </div>

          <div className="grid md:grid-cols-2 gap-3">
            <div>
              <Equation>
                {`**General form:**\n\n$$p(x)=a_n x^n + a_{n-1}x^{n-1}+\\cdots + a_1 x + a_0,\\quad a_n\\neq0.$$`}
              </Equation>
              <Note>
                Use the 3D surface to spot valleys: deep valleys indicate **near roots** in the complex plane.
              </Note>
            </div>
            <div>
              <Equation color={theme.accent2}>
                {`**Factorization:**\n\n$$p(x)=a_n\\prod_{k=1}^n (x-r_k).$$`}
              </Equation>
              <Note>Clustered roots can degrade conditioning — small coefficient perturbations may move them significantly.</Note>
            </div>
          </div>

          {/* Exercises */}
          <div className="mt-2 space-y-3">
            <Exercise
              title="Exercise 7.1.1 — Visual identification of roots"
              body={`Use \\(p(x)=(x-1)(x-2)(x-3)\\).  
Change coefficients to \\(1\\; -6\\; 11\\; -6\\) and confirm the **3D surface** shows minima near \\(1,2,3\\) on the real axis.`}
              solution={`For \\(p(x)=x^3-6x^2+11x-6\\), the minima of \\(|p(x+iy)|\\) occur near \\(x\\in\\{1,2,3\\}\\) and \\(y\\approx 0\\), matching the real simple roots.`}
            />
            <Exercise
              title="Exercise 7.1.2 — Cauchy bound"
              body={`For \\(p(x)=2x^4-3x^3+5x-9\\), compute the Cauchy bound \\(R\\).`}
              solution={`Here \\(a_4=2\\), other coefficients: \\([-3, 0, 5, -9]\\).  
\\( \\max |a_i/a_4| = \\max(1.5, 0, 2.5, 4.5) = 4.5\\).  
Thus \\(R=1+4.5=5.5\\). All roots lie in \\(|z|\\le 5.5\\).`}
            />
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// 7.2 Panel — Horner, synthetic division, derivatives
// ======================================================================
function Section72() {
  const [coeffStr, setCoeffStr] = useState("1 -3 2"); // x^2 - 3x + 2
  const [xval, setXval] = useState(1.5);

  const coeffs = useMemo(() => parseCoeffs(coeffStr), [coeffStr]);
  const d = useMemo(() => polyDeriv(coeffs), [coeffs]);
  const X = useMemo(() => linspace(-6, 6, 360), []);
  const Y = useMemo(() => polyEvalArray(coeffs, X), [coeffs, X]);

  const pAtX = useMemo(() => horner(coeffs, Number(xval)), [coeffs, xval]);
  const dpAtX = useMemo(() => horner(d, Number(xval)), [d, xval]);

  // 3D surface of |p| for additional view
  const surface = useMemo(() => surfaceModulus(coeffs, -4, 4, 50, -4, 4, 50), [coeffs]);

  const [deflateAt, setDeflateAt] = useState(1);
  const defl = useMemo(() => syntheticDiv(coeffs, Number(deflateAt)), [coeffs, deflateAt]);

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <FunctionSquare className="w-5 h-5" />
            7.2 Computing with Polynomials — Horner & Division
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="grid md:grid-cols-3 gap-3">
            <Labeled label="Coefficients">
              <Input
                value={coeffStr}
                onChange={(e) => setCoeffStr(e.target.value)}
                className="bg-zinc-800 text-white"
              />
            </Labeled>
            <Labeled label="Evaluate at x =">
              <Input
                value={xval}
                onChange={(e) => setXval(e.target.value)}
                className="bg-zinc-800 text-white"
              />
            </Labeled>
            <Summary
              items={[
                { label: "p(x)", value: pAtX.toPrecision(6) },
                { label: "p'(x)", value: dpAtX.toPrecision(6) },
                { label: "Degree", value: Math.max(0, coeffs.length - 1) },
              ]}
            />
          </div>

          <div className="grid md:grid-cols-2 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2">
              <Plot
                data={[
                  {
                    x: X,
                    y: Y,
                    type: "scatter",
                    mode: "lines",
                    line: { width: 2, color: theme.accent },
                    name: "p(x)",
                  },
                  {
                    x: [Number(xval)],
                    y: [pAtX],
                    type: "scatter",
                    mode: "markers",
                    marker: { size: 8, color: theme.accent2 },
                    name: "p(x0)",
                  },
                ]}
                layout={{
                  title: "Horner evaluation & mark",
                  paper_bgcolor: "rgba(0,0,0,0)",
                  plot_bgcolor: "rgba(0,0,0,0)",
                  font: { color: "#e5e7eb" },
                  margin: { l: 40, r: 10, b: 40, t: 40 },
                }}
                style={{ width: "100%", height: 300 }}
                useResizeHandler
              />
            </div>
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2">
              <Plot
                data={[
                  {
                    x: surface.xs,
                    y: surface.ys,
                    z: surface.Z,
                    type: "surface",
                    colorscale: "Viridis",
                    showscale: false,
                  },
                ]}
                layout={{
                  title: "|p(x+iy)| surface (3D)",
                  paper_bgcolor: "rgba(0,0,0,0)",
                  plot_bgcolor: "rgba(0,0,0,0)",
                  font: { color: "#e5e7eb" },
                  scene: {
                    xaxis: { title: "Re", color: "#e5e7eb" },
                    yaxis: { title: "Im", color: "#e5e7eb" },
                    zaxis: { title: "|p|", color: "#e5e7eb" },
                  },
                  margin: { l: 0, r: 0, b: 0, t: 40 },
                }}
                style={{ width: "100%", height: 300 }}
                useResizeHandler
              />
            </div>
          </div>

          <div className="grid md:grid-cols-3 gap-3">
            <div className="md:col-span-2">
              <Equation>{`**Horner**\n\n$$p(x)=(\\cdots((a_nx+a_{n-1})x+a_{n-2})x+\\cdots)+a_0.$$`}</Equation>
              <Note>Use synthetic division to deflate after you find a root (or near-root) to simplify further searches.</Note>
            </div>
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-3">
              <div className="text-sm text-zinc-300 mb-2">Synthetic division by (x − r)</div>
              <div className="flex gap-2 items-center">
                <Input
                  value={deflateAt}
                  onChange={(e) => setDeflateAt(e.target.value)}
                  className="bg-zinc-800 text-white w-28"
                />
                <div className="text-xs text-zinc-400">Try r = 1 or r = 2 for x² − 3x + 2</div>
              </div>
              <div className="text-sm text-zinc-300 mt-2">
                <div>
                  <span className="text-zinc-400">Quotient:</span>{" "}
                  [{defl.q.map((v) => Number.isFinite(v) ? v.toPrecision(4) : "NaN").join(", ")}]
                </div>
                <div>
                  <span className="text-zinc-400">Remainder:</span>{" "}
                  <span className="text-emerald-300">{Number.isFinite(defl.remainder) ? defl.remainder.toPrecision(6) : "NaN"}</span>
                </div>
              </div>
            </div>
          </div>

          {/* Exercises */}
          <div className="space-y-3">
            <Exercise
              title="Exercise 7.2.1 — Evaluate & differentiate"
              body={`For \\(p(x)=x^2-3x+2\\), compute \\(p(1.5)\\) and \\(p'(1.5)\\). Compare with the app output.`}
              solution={`\\(p(1.5)=1.5^2-3(1.5)+2=2.25-4.5+2=-0.25\\).  
\\(p'(x)=2x-3\\Rightarrow p'(1.5)=0\\).`}
            />
            <Exercise
              title="Exercise 7.2.2 — Synthetic division"
              body={`Divide \\(x^3 - 6x^2 + 11x - 6\\) by \\((x-2)\\) using synthetic division.`}
              solution={`Quotient: \\(x^2-4x+3\\), Remainder: 0.`}
            />
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// 7.3 Panel — Conventional methods (Newton & Secant) + 3D error surface
// ======================================================================
function Section73() {
  const [coeffStr, setCoeffStr] = useState("1 -6 11 -6");
  const [x0, setX0] = useState(3.5);
  const [x1, setX1] = useState(2.3);
  const [lambda, setLambda] = useState(1);

  const coeffs = useMemo(() => parseCoeffs(coeffStr), [coeffStr]);

  const newtonTrace = useMemo(() => newtonPoly(coeffs, Number(x0), 30, 1e-10, Number(lambda)), [coeffs, x0, lambda]);
  const secantTrace = useMemo(() => secantPoly(coeffs, Number(x0), Number(x1), 30, 1e-10), [coeffs, x0, x1]);

  const X = useMemo(() => linspace(-1, 5, 400), []);
  const Y = useMemo(() => polyEvalArray(coeffs, X), [coeffs, X]);

  // 3D error surface: (x, iter) -> |p(x)|, just a toy visual by scanning x for first few iterations
  const ErrSurf = useMemo(() => {
    const xs = linspace(-1, 5, 50);
    const iters = linspace(0, 10, 11); // 0..10
    const Z = iters.map((it) => xs.map((x) => Math.log10(1e-12 + Math.abs(horner(coeffs, x)))));
    return { xs, iters, Z };
  }, [coeffs]);

  const newtonXs = newtonTrace.map((r) => r.x);
  const newtonYs = newtonTrace.map((r) => horner(coeffs, r.x));

  const secantXs = secantTrace.map((r) => r.b);
  const secantYs = secantTrace.map((r) => horner(coeffs, r.b));

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Target className="w-5 h-5" />
            7.3 Conventional Methods — Newton & Secant
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="grid md:grid-cols-4 gap-3">
            <Labeled label="Coefficients">
              <Input value={coeffStr} onChange={(e) => setCoeffStr(e.target.value)} className="bg-zinc-800 text-white" />
            </Labeled>
            <Labeled label="Newton initial x₀">
              <Input value={x0} onChange={(e) => setX0(e.target.value)} className="bg-zinc-800 text-white" />
            </Labeled>
            <Labeled label="Secant x₀">
              <Input value={x1} onChange={(e) => setX1(e.target.value)} className="bg-zinc-800 text-white" />
            </Labeled>
            <Labeled label="Damping λ (Newton)">
              <Input value={lambda} onChange={(e) => setLambda(e.target.value)} className="bg-zinc-800 text-white" />
            </Labeled>
          </div>

          <div className="grid md:grid-cols-2 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2">
              <Plot
                data={[
                  {
                    x: X,
                    y: Y,
                    type: "scatter",
                    mode: "lines",
                    line: { color: theme.accent },
                    name: "p(x)",
                  },
                  {
                    x: newtonXs,
                    y: newtonYs,
                    type: "scatter",
                    mode: "markers+lines",
                    marker: { size: 6, symbol: "circle" },
                    line: { dash: "dot" },
                    name: "Newton path",
                  },
                  {
                    x: secantXs,
                    y: secantYs,
                    type: "scatter",
                    mode: "markers+lines",
                    marker: { size: 6, symbol: "square" },
                    line: { dash: "dash" },
                    name: "Secant path",
                  },
                ]}
                layout={{
                  title: "Iteration paths on p(x)",
                  paper_bgcolor: "rgba(0,0,0,0)",
                  plot_bgcolor: "rgba(0,0,0,0)",
                  font: { color: "#e5e7eb" },
                  margin: { l: 40, r: 10, b: 40, t: 40 },
                }}
                style={{ width: "100%", height: 320 }}
                useResizeHandler
              />
            </div>
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2">
              <Plot
                data={[
                  {
                    x: ErrSurf.xs,
                    y: ErrSurf.iters,
                    z: ErrSurf.Z,
                    type: "surface",
                    colorscale: "Viridis",
                    showscale: false,
                  },
                ]}
                layout={{
                  title: "log₁₀( |p(x)| ) surface",
                  paper_bgcolor: "rgba(0,0,0,0)",
                  plot_bgcolor: "rgba(0,0,0,0)",
                  font: { color: "#e5e7eb" },
                  scene: {
                    xaxis: { title: "x", color: "#e5e7eb" },
                    yaxis: { title: "iteration (index)", color: "#e5e7eb" },
                    zaxis: { title: "log₁₀|p|", color: "#e5e7eb" },
                  },
                }}
                style={{ width: "100%", height: 320 }}
                useResizeHandler
              />
            </div>
          </div>

          <div className="grid md:grid-cols-2 gap-3">
            <Equation>{`**Newton update**\n\n$$x_{k+1}=x_k-\\frac{p(x_k)}{p'(x_k)}.$$`}</Equation>
            <Equation color={theme.accent2}>{`**Secant update**\n\n$$x_{k+1}=x_k- p(x_k)\\,\\frac{x_k-x_{k-1}}{p(x_k)-p(x_{k-1})}.$$`}</Equation>
          </div>

          {/* Exercises */}
          <div className="space-y-3">
            <Exercise
              title="Exercise 7.3.1 — Damping and stability"
              body={`Try \\(x_0=3.5\\) with \\(\\lambda=0.5\\) and \\(1\\). Compare iterations to the nearest root. Why might damping help?`}
              solution={`Damping reduces overshoot when \\(|p'(x_k)|\\) is small or when starting far from a root; convergence may become more stable though slower.`}
            />
            <Exercise
              title="Exercise 7.3.2 — Secant vs Newton"
              body={`With \\(x_0=3.5\\), \\(x_1=2.3\\), compare iteration counts. When derivative is cheap, Newton is often faster; when not, Secant avoids derivative computation.`}
              solution={`Newton is locally quadratic; Secant is superlinear (≈1.618). Results depend on initial guesses and conditioning.`}
            />
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// 7.4 Panel — Müller’s method visualization (2D + 3D)
// ======================================================================
function Section74() {
  const [coeffStr, setCoeffStr] = useState("1 0 0 1"); // x^3 + 1 (has complex roots)
  const [x0, setX0] = useState(-1.5);
  const [x1, setX1] = useState(0.2);
  const [x2, setX2] = useState(1.0);
  const [steps, setSteps] = useState(6);

  const coeffs = useMemo(() => parseCoeffs(coeffStr), [coeffStr]);

  const trace = useMemo(() => {
    const xs = [Number(x0), Number(x1), Number(x2)];
    const path = [{ k: 0, x: xs[2], fx: horner(coeffs, xs[2]) }];
    for (let k = 0; k < Number(steps); k++) {
      const nxt = mullerStep(coeffs, xs[0], xs[1], xs[2]);
      xs[0] = xs[1];
      xs[1] = xs[2];
      xs[2] = nxt;
      path.push({ k: k + 1, x: xs[2], fx: horner(coeffs, xs[2]) });
    }
    return path;
  }, [coeffs, x0, x1, x2, steps]);

  const X = useMemo(() => linspace(-3, 3, 400), []);
  const Y = useMemo(() => polyEvalArray(coeffs, X), [coeffs, X]);
  const surface = useMemo(() => surfaceModulus(coeffs, -3, 3, 60, -3, 3, 60), [coeffs]);

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Ruler className="w-5 h-5" />
            7.4 Müller’s Method — Quadratic Interpolation Steps
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="grid md:grid-cols-5 gap-3">
            <Labeled label="Coefficients">
              <Input value={coeffStr} onChange={(e) => setCoeffStr(e.target.value)} className="bg-zinc-800 text-white" />
            </Labeled>
            <Labeled label="x₀">
              <Input value={x0} onChange={(e) => setX0(e.target.value)} className="bg-zinc-800 text-white" />
            </Labeled>
            <Labeled label="x₁">
              <Input value={x1} onChange={(e) => setX1(e.target.value)} className="bg-zinc-800 text-white" />
            </Labeled>
            <Labeled label="x₂">
              <Input value={x2} onChange={(e) => setX2(e.target.value)} className="bg-zinc-800 text-white" />
            </Labeled>
            <Labeled label="Steps">
              <Input value={steps} onChange={(e) => setSteps(e.target.value)} className="bg-zinc-800 text-white" />
            </Labeled>
          </div>

          <div className="grid md:grid-cols-2 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2">
              <Plot
                data={[
                  {
                    x: X,
                    y: Y,
                    type: "scatter",
                    mode: "lines",
                    line: { color: theme.accent },
                    name: "p(x)",
                  },
                  {
                    x: trace.map((t) => t.x),
                    y: trace.map((t) => t.fx),
                    type: "scatter",
                    mode: "markers+lines",
                    marker: { size: 6, symbol: "diamond-open" },
                    line: { dash: "dot" },
                    name: "Müller path",
                  },
                ]}
                layout={{
                  title: "Trajectory of Müller updates",
                  paper_bgcolor: "rgba(0,0,0,0)",
                  plot_bgcolor: "rgba(0,0,0,0)",
                  font: { color: "#e5e7eb" },
                }}
                style={{ width: "100%", height: 320 }}
                useResizeHandler
              />
            </div>
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2">
              <Plot
                data={[
                  { x: surface.xs, y: surface.ys, z: surface.Z, type: "surface", colorscale: "Viridis", showscale: false },
                ]}
                layout={{
                  title: "|p(x+iy)| surface (3D)",
                  paper_bgcolor: "rgba(0,0,0,0)",
                  plot_bgcolor: "rgba(0,0,0,0)",
                  font: { color: "#e5e7eb" },
                }}
                style={{ width: "100%", height: 320 }}
                useResizeHandler
              />
            </div>
          </div>

          <div className="grid md:grid-cols-2 gap-3">
            <Equation>{`**Update**\n\n$$x_{k+1}=x_k-\\frac{2c}{b\\pm\\sqrt{b^2-4ac}}.$$`}</Equation>
            <Note>Müller often converges to complex roots and may be less sensitive to derivative issues than Newton.</Note>
          </div>

          {/* Exercises */}
          <div className="space-y-3">
            <Exercise
              title="Exercise 7.4.1 — Complex roots of x³+1"
              body={`Run Müller with initial guesses around \\(-1, 0.2, 1\\). Which root do you approach? Compare with the known complex cube roots of \\(-1\\).`}
              solution={`\\(x^3+1=0\\Rightarrow x=-1,\\; \\frac{1}{2}\\pm\\frac{\\sqrt{3}}{2}i\\).  
Depending on initial seeds, Müller may converge to any of these; paths can cross the complex plane (here we visualize with real x only).`}
            />
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// 7.5 Panel — Bairstow’s method (quadratic factor) + residual surface
// ======================================================================
function Section75() {
  const [coeffStr, setCoeffStr] = useState("1 -6 11 -6");
  const [r, setR] = useState(-1);
  const [s, setS] = useState(1);
  const [maxIt, setMaxIt] = useState(30);

  const coeffs = useMemo(() => parseCoeffs(coeffStr), [coeffStr]);

  const result = useMemo(() => bairstow(coeffs, Number(r), Number(s), Number(maxIt), 1e-10), [coeffs, r, s, maxIt]);

  // residual surface over r,s grid: log10(|F|+|G|)
  const RSsurf = useMemo(() => {
    const rVals = linspace(-4, 4, 45);
    const sVals = linspace(-4, 4, 45);
    // compute residuals cheaply by synthetic division recurrence
    const Z = sVals.map((sv) =>
      rVals.map((rv) => {
        // compute remainders F=b[1], G=b[0]
        const n = coeffs.length - 1;
        const b = new Array(n + 1).fill(0);
        b[n] = coeffs[n];
        b[n - 1] = coeffs[n - 1] + rv * b[n];
        for (let i = n - 2; i >= 0; i--) b[i] = coeffs[i] + rv * b[i + 1] + sv * (b[i + 2] ?? 0);
        const F = b[1];
        const G = b[0];
        return Math.log10(1e-12 + Math.abs(F) + Math.abs(G));
      })
    );
    return { rVals, sVals, Z };
  }, [coeffs]);

  const roots = useMemo(() => {
    if (!result) return [];
    const rs = rootsFromRS(result.r, result.s);
    return rs;
  }, [result]);

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Braces className="w-5 h-5" />
            7.5 Bairstow’s Method — Quadratic Factor Extraction
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="grid md:grid-cols-5 gap-3">
            <Labeled label="Coefficients">
              <Input value={coeffStr} onChange={(e) => setCoeffStr(e.target.value)} className="bg-zinc-800 text-white" />
            </Labeled>
            <Labeled label="Initial r">
              <Input value={r} onChange={(e) => setR(e.target.value)} className="bg-zinc-800 text-white" />
            </Labeled>
            <Labeled label="Initial s">
              <Input value={s} onChange={(e) => setS(e.target.value)} className="bg-zinc-800 text-white" />
            </Labeled>
            <Labeled label="Max iterations">
              <Input value={maxIt} onChange={(e) => setMaxIt(e.target.value)} className="bg-zinc-800 text-white" />
            </Labeled>
            <Summary
              items={[
                { label: "Converged", value: result.ok ? "Yes" : "No" },
                { label: "r, s", value: `${result.r.toFixed(4)}, ${result.s.toFixed(4)}` },
                { label: "Roots of factor", value: `${roots[0]} , ${roots[1]}` },
              ]}
            />
          </div>

          <div className="grid md:grid-cols-2 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2">
              <Plot
                data={[
                  {
                    x: result.iters.map((t) => t.r),
                    y: result.iters.map((t) => t.s),
                    type: "scatter",
                    mode: "markers+lines",
                    marker: { size: 6, symbol: "circle" },
                    name: "Iteration in (r,s)",
                  },
                ]}
                layout={{
                  title: "Bairstow path in (r,s)",
                  xaxis: { title: "r", color: "#e5e7eb" },
                  yaxis: { title: "s", color: "#e5e7eb" },
                  paper_bgcolor: "rgba(0,0,0,0)",
                  plot_bgcolor: "rgba(0,0,0,0)",
                  font: { color: "#e5e7eb" },
                }}
                style={{ width: "100%", height: 320 }}
                useResizeHandler
              />
            </div>
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2">
              <Plot
                data={[
                  {
                    x: RSsurf.rVals,
                    y: RSsurf.sVals,
                    z: RSsurf.Z,
                    type: "surface",
                    colorscale: "Viridis",
                    showscale: false,
                  },
                ]}
                layout={{
                  title: "log₁₀(|F|+|G|) over (r,s)",
                  scene: {
                    xaxis: { title: "r", color: "#e5e7eb" },
                    yaxis: { title: "s", color: "#e5e7eb" },
                    zaxis: { title: "residual", color: "#e5e7eb" },
                  },
                  paper_bgcolor: "rgba(0,0,0,0)",
                  plot_bgcolor: "rgba(0,0,0,0)",
                  font: { color: "#e5e7eb" },
                }}
                style={{ width: "100%", height: 320 }}
                useResizeHandler
              />
            </div>
          </div>

          <div className="grid md:grid-cols-2 gap-3">
            <Equation>{`**Quadratic factor**\n\nFind \\(r,s\\) s.t. dividing by \\(x^2+rx+s\\) leaves zero remainder.`}</Equation>
            <Note>The residual surface reveals basins of attraction in the (r,s)-plane.</Note>
          </div>

          {/* Exercises */}
          <div className="space-y-3">
            <Exercise
              title="Exercise 7.5.1 — Factor (x−1)(x−2)(x−3)"
              body={`Run Bairstow on \\(x^3-6x^2+11x-6\\) with initial \\(r=-3, s=2\\).  
Extract a quadratic and then find its roots; the remaining root follows by deflation.`}
              solution={`One quadratic factor is \\(x^2-3x+2\\) (roots 1 and 2). The deflated linear factor gives the remaining root 3.`}
            />
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// 7.6 Panel — Other methods: Companion matrix via |p(z)| landscape
// ======================================================================
function Section76() {
  const [coeffStr, setCoeffStr] = useState("1 0 -7 6"); // x^3 -7x + 6
  const coeffs = useMemo(() => parseCoeffs(coeffStr), [coeffStr]);

  const surf = useMemo(() => surfaceModulus(coeffs, -4, 4, 64, -4, 4, 64), [coeffs]);

  // 2D contour of |p(z)| over complex plane (real-imag heatmap)
  const contour = useMemo(() => {
    return {
      x: surf.xs,
      y: surf.ys,
      z: surf.Z,
    };
  }, [surf]);

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Calculator className="w-5 h-5" />
            7.6 Other Methods — Companion Matrix & Complex Landscape
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="grid md:grid-cols-3 gap-3">
            <Labeled label="Coefficients">
              <Input value={coeffStr} onChange={(e) => setCoeffStr(e.target.value)} className="bg-zinc-800 text-white" />
            </Labeled>
            <div className="md:col-span-2">
              <Equation>{`**Companion matrix**\n\nFor $$p(x)=x^n+a_{n-1}x^{n-1}+\\cdots+a_0,$$ eigenvalues of $$C$$ equal roots of \\(p\\).`}</Equation>
            </div>
          </div>

          <div className="grid md:grid-cols-2 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2">
              <Plot
                data={[
                  {
                    x: contour.x,
                    y: contour.y,
                    z: contour.z,
                    type: "contour",
                    colorscale: "Portland",
                    contours: { coloring: "heatmap" },
                  },
                ]}
                layout={{
                  title: "Heatmap/Contours of |p(x+iy)|",
                  xaxis: { title: "Re", color: "#e5e7eb" },
                  yaxis: { title: "Im", color: "#e5e7eb" },
                  paper_bgcolor: "rgba(0,0,0,0)",
                  plot_bgcolor: "rgba(0,0,0,0)",
                  font: { color: "#e5e7eb" },
                }}
                style={{ width: "100%", height: 320 }}
                useResizeHandler
              />
            </div>
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2">
              <Plot
                data={[
                  {
                    x: surf.xs,
                    y: surf.ys,
                    z: surf.Z,
                    type: "surface",
                    colorscale: "Viridis",
                    showscale: false,
                  },
                ]}
                layout={{
                  title: "|p| surface (3D) — valleys ⇢ roots",
                  scene: {
                    xaxis: { title: "Re", color: "#e5e7eb" },
                    yaxis: { title: "Im", color: "#e5e7eb" },
                    zaxis: { title: "|p|", color: "#e5e7eb" },
                  },
                  paper_bgcolor: "rgba(0,0,0,0)",
                  plot_bgcolor: "rgba(0,0,0,0)",
                  font: { color: "#e5e7eb" },
                }}
                style={{ width: "100%", height: 320 }}
                useResizeHandler
              />
            </div>
          </div>

          <div className="space-y-3">
            <Exercise
              title="Exercise 7.6.1 — Landscapes and minima"
              body={`For \\(p(x)=x^3-7x+6\\), use the heatmap and 3D surface to estimate where roots lie in \\(\\mathbb{C}\\).  
Do you locate minima near real axis?`}
              solution={`This polynomial factors as \\((x-1)(x-2)(x+3)\\) so roots at \\(1,2,-3\\). The landscape shows three valleys on Re-axis.`}
            />
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// 7.7 Panel — Root Location Theorems (Descartes, Cauchy) + complex disk
// ======================================================================
function Section77() {
  const [coeffStr, setCoeffStr] = useState("1 -1 -7 1 6"); // x^4 - x^3 - 7x^2 + x + 6
  const coeffs = useMemo(() => parseCoeffs(coeffStr), [coeffStr]);

  const changes = useMemo(() => descartsSignChanges(coeffs), [coeffs]);
  const bound = useMemo(() => cauchyBound(coeffs), [coeffs]);

  // Simple complex-disk visualization: draw |z|=R circle and heatmap
  const surf = useMemo(() => surfaceModulus(coeffs, -bound, bound, 60, -bound, bound, 60), [coeffs, bound]);

  // circle param
  const circ = useMemo(() => {
    const t = linspace(0, 2 * Math.PI, 200);
    return {
      x: t.map((th) => bound * Math.cos(th)),
      y: t.map((th) => bound * Math.sin(th)),
    };
  }, [bound]);

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <PlusCircle className="w-5 h-5" />
            7.7 Root Location Theorems — Bounds & Counts
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="grid md:grid-cols-3 gap-3">
            <Labeled label="Coefficients">
              <Input value={coeffStr} onChange={(e) => setCoeffStr(e.target.value)} className="bg-zinc-800 text-white" />
            </Labeled>
            <Summary
              items={[
                { label: "Descartes sign changes (+)", value: changes },
                { label: "Cauchy bound R", value: bound.toFixed(4) },
                { label: "Degree", value: Math.max(0, coeffs.length - 1) },
              ]}
            />
            <div className="text-sm text-zinc-300">
              <Equation>{`**Cauchy bound** \\(R=1+\\max |a_i/a_n|\\).  
**Descartes:** # of positive roots ≤ sign changes (difference by even integer).`}</Equation>
            </div>
          </div>

          <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2">
            <Plot
              data={[
                {
                  x: surf.xs,
                  y: surf.ys,
                  z: surf.Z,
                  type: "contour",
                  colorscale: "Portland",
                  contours: { coloring: "heatmap" },
                  name: "|p| heat",
                },
                {
                  x: circ.x,
                  y: circ.y,
                  type: "scatter",
                  mode: "lines",
                  line: { color: theme.accent, width: 2 },
                  name: `|z|=R (Cauchy)`,
                },
              ]}
              layout={{
                title: "Root bounds in complex plane",
                xaxis: { title: "Re", color: "#e5e7eb", scaleanchor: "y", scaleratio: 1 },
                yaxis: { title: "Im", color: "#e5e7eb" },
                paper_bgcolor: "rgba(0,0,0,0)",
                plot_bgcolor: "rgba(0,0,0,0)",
                font: { color: "#e5e7eb" },
              }}
              style={{ width: "100%", height: 360 }}
              useResizeHandler
            />
          </div>

          {/* Exercises */}
          <div className="space-y-3">
            <Exercise
              title="Exercise 7.7.1 — Count possible positive roots"
              body={`For \\(p(x)=x^4-x^3-7x^2+x+6\\), compute sign changes and list possible numbers of positive roots.`}
              solution={`Coefficients: \\([1,-1,-7,1,6]\\) → signs \\([+,-,-,+,+]\\).  
Sign changes: from \\(+\\to-\\) (1), \\(-\\to-\\) (0), \\(-\\to+\\) (2), \\(+\\to+\\) (2). So **2** changes. Possible positive roots: 2 or 0.`}
            />
            <Exercise
              title="Exercise 7.7.2 — Cauchy radius"
              body={`Compute \\(R\\) for the same polynomial and verify that the valley of \\(|p|\\) lies inside \\(|z|\\le R\\).`}
              solution={`\\(a_4=1\\); maxima of \\(|a_i/a_n|=\\max(1,7,1,6)=7\\). Thus \\(R=8\\). The heatmap shows low \\(|p|\\) inside this disk.`}
            />
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// Documentation Panel — left nav + content
// ======================================================================
function ChapterDocs() {
  const [open, setOpen] = useState("7.1 Polynomials in Engineering and Science");

  return (
    <div className="w-full mt-6 bg-zinc-900 border border-zinc-700 rounded-xl shadow-lg overflow-hidden">
      <div className="p-3 border-b border-zinc-700 flex items-center justify-between">
        <div className="flex items-center gap-3">
          <BookOpen className="w-5 h-5 text-cyan-400" />
          <div>
            <div className="text-zinc-100 font-semibold">Chapter 7 — Documentation</div>
            <div className="text-zinc-400 text-xs">Theory, notes, and references</div>
          </div>
        </div>
        
      </div>

      <div className="grid grid-cols-1 md:grid-cols-4">
        <div className="col-span-1 border-r border-zinc-800">
          {Object.keys(docsText).map((k) => (
            <button
              key={k}
              onClick={() => setOpen(k)}
              className={`w-full p-3 text-left hover:bg-zinc-800/40 flex items-center justify-between ${open === k ? "bg-zinc-800/20" : ""}`}
            >
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
            <div className="mt-3  text-zinc-400 text-xs">
              References: Atkinson, Burden & Faires, Trefethen & Bau, Press et al. (Numerical Recipes).
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}

// ======================================================================
// Page assembly
// ======================================================================
export default function Chapter7() {
  return (
    <div className={`p-6 space-y-6 bg-zinc-950`}>
      <motion.div initial={{ opacity: 0, y: -10 }} animate={{ opacity: 1, y: 0 }}>
        <h1 className="text-3xl font-bold" style={{ color: theme.accent2 }}>
          Roots of Polynomials
        </h1>
        <div className="text-zinc-400">Visualization-first learning: every concept, plotted.</div>
      </motion.div>

      <Section71 />
      <Section72 />
      <Section73 />
      <Section74 />
      <Section75 />
      <Section76 />
      <Section77 />

      <ChapterDocs />
      <BottomBar/>
    </div>
  );
}
