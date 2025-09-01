// src/pages/Chapter27.jsx
/* 
  CHAPTER 27 — Boundary-Value and Eigenvalue Problems
  Sections:
    27.1 General Methods for Boundary-Value Problems
      - Shooting method (with RK4 solver)
      - Finite-difference method (Thomas algorithm)
    27.2 Eigenvalue Problems
      - Sturm–Liouville finite-difference discretization
      - Inverse iteration with Rayleigh quotient
      - Compare to analytic λ_n for y'' + λy = 0 with y(0)=0, y(π)=0
    27.3 ODEs and Eigenvalues with Software Packages
      - Guidance & reproducible parameter widgets
    Problems
      - Interactive exercises with togglable solutions
*/

import React, { useMemo, useState, useEffect, useRef, Fragment } from "react";
import Plot from "react-plotly.js";
import { motion } from "framer-motion";
import {
  BookOpen,
  Gauge,
  Route,
  Sigma,
  Binary,
  Settings,
  ChevronRight,
  ChevronDown,
  Play,
  GitBranch,
  Landmark,
  ListChecks,
  LineChart,
  Dot,
  Info,
} from "lucide-react";
import ReactMarkdown from "react-markdown";
import remarkMath from "remark-math";
import rehypeKatex from "rehype-katex";
import "katex/dist/katex.min.css";
import BottomBar from "../components/BottomBar";
// If you have shadcn/ui or your own UI kit, align these imports:
import { Card, CardContent } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Textarea } from "@/components/ui/textarea";

// ------------------------------
// Small presentational helpers
// ------------------------------
const SectionTitle = ({ icon: Icon, title, k }) => (
  <div className="flex items-center gap-2">
    <div className="p-2 rounded-xl bg-zinc-800/70 border border-zinc-700">
      <Icon className="h-5 w-5 text-cyan-300" />
    </div>
    <h2 className="text-lg sm:text-xl font-semibold tracking-wide">
      <span className="text-cyan-300 mr-2">{k}</span>
      <span className="text-zinc-100">{title}</span>
    </h2>
  </div>
);

const Chip = ({ children, tone = "cyan" }) => {
  const colors =
    tone === "cyan"
      ? "bg-cyan-950/40 text-cyan-200 border-cyan-700"
      : tone === "emerald"
      ? "bg-emerald-950/40 text-emerald-200 border-emerald-700"
      : tone === "amber"
      ? "bg-amber-950/40 text-amber-200 border-amber-700"
      : "bg-zinc-900/60 text-zinc-300 border-zinc-700";
  return (
    <span className={`px-2 py-1 rounded-lg text-xs border ${colors}`}>{children}</span>
  );
};

function Equation({ children }) {
  return (
    <div className="my-2">
      <div className="px-3 py-2 rounded-lg border text-gray-300 border-zinc-700 bg-zinc-900/80 overflow-x-auto">
        <div className="text-xs uppercase tracking-widest mb-1">
          Formula
        </div>
        <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>
          {children}
        </ReactMarkdown>
      </div>
    </div>
  );
}

const Note = ({ children }) => (
  <div className="text-sm text-zinc-400 leading-relaxed">{children}</div>
);

const Field = ({ label, right, children }) => (
  <label className="text-sm text-zinc-300">
    <div className="flex items-center justify-between mb-1">
      <span className="text-zinc-300">{label}</span>
      {right ? <span className="text-zinc-500 text-xs">{right}</span> : null}
    </div>
    {children}
  </label>
);

const NumberInput = (props) => (
  <Input
    {...props}
    className="bg-zinc-800/90 text-zinc-100 border-zinc-700 focus-visible:ring-cyan-600"
    type="number"
    step={props.step ?? "any"}
  />
);

const IntegerInput = (props) => (
  <Input
    {...props}
    className="bg-zinc-800/90 text-zinc-100 border-zinc-700 focus-visible:ring-cyan-600"
    type="number"
  />
);

const Toggle = ({ open, onClick, children }) => (
  <button
    onClick={onClick}
    className="flex items-center gap-2 px-2 py-1 rounded-lg bg-zinc-900/60 border border-zinc-700 hover:border-cyan-700 text-zinc-200"
  >
    {open ? <ChevronDown className="h-4 w-4" /> : <ChevronRight className="h-4 w-4" />}
    <span className="text-sm">{children}</span>
  </button>
);

// ------------------------------
// Math / Numerics utilities
// ------------------------------

// RK4 for vector ODE y' = F(x, y), y is array
function rk4Step(F, x, y, h) {
  const k1 = F(x, y);
  const yk2 = y.map((yi, i) => yi + 0.5 * h * k1[i]);
  const k2 = F(x + 0.5 * h, yk2);
  const yk3 = y.map((yi, i) => yi + 0.5 * h * k2[i]);
  const k3 = F(x + 0.5 * h, yk3);
  const yk4 = y.map((yi, i) => yi + h * k3[i]);
  const k4 = F(x + h, yk4);
  const yNext = y.map((yi, i) => yi + (h / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]));
  return yNext;
}

// Solve y'' = f(x,y,y') with shooting method by transforming to first-order system:
// y1 = y, y2 = y'  => y1' = y2, y2' = f(x, y1, y2)
function shootSolve({ f, a, b, alpha, beta, s, n }) {
  const h = (b - a) / n;
  let x = a;
  let y1 = alpha;
  let y2 = s; // initial slope guess
  const xs = [x];
  const ys = [y1];

  const F = (xx, yy) => {
    const y1 = yy[0];
    const y2 = yy[1];
    const y2p = f(xx, y1, y2);
    return [y2, y2p];
  };

  for (let i = 0; i < n; i++) {
    const y = [y1, y2];
    const yNext = rk4Step(F, x, y, h);
    x += h;
    y1 = yNext[0];
    y2 = yNext[1];
    xs.push(x);
    ys.push(y1);
  }

  const endError = y1 - beta;
  return { xs, ys, endError, slope: s };
}

// Simple secant-based shooting: iterate slope guesses to hit y(b) = beta
function shootingSecant({ f, a, b, alpha, beta, n, s0, s1, maxIter = 12, tol = 1e-6 }) {
  let sol0 = shootSolve({ f, a, b, alpha, beta, s: s0, n });
  let sol1 = shootSolve({ f, a, b, alpha, beta, s: s1, n });

  for (let k = 0; k < maxIter; k++) {
    const e0 = sol0.endError;
    const e1 = sol1.endError;
    const ds = s1 - s0;
    const de = e1 - e0;
    const s2 = s1 - e1 * (ds / (de === 0 ? 1e-14 : de));
    const sol2 = shootSolve({ f, a, b, alpha, beta, s: s2, n });

    if (Math.abs(sol2.endError) < tol) return sol2;

    s0 = s1;
    s1 = s2;
    sol0 = sol1;
    sol1 = sol2;
  }
  return sol1; // best effort
}

// Finite difference for y'' = f(x,y,y') with linearization assumption around y (nonlinear needs Newton);
// here we implement linear second-order BVP: y'' = p(x)*y' + q(x)*y + r(x), with Dirichlet BC.
function fdBVP({ a, b, alpha, beta, N, p, q, r }) {
  // interior points: i=1..N-1, grid h = (b-a)/N
  const h = (b - a) / N;
  const x = Array.from({ length: N + 1 }, (_, i) => a + i * h);

  // Tridiagonal A*y = d
  // For y'' = p y' + q y + r
  // Discretize: (y_{i-1} - 2 y_i + y_{i+1})/h^2 = p_i * (y_{i+1} - y_{i-1})/(2h) + q_i y_i + r_i
  // => (-1/h^2 - p_i/(2h)) y_{i-1} + (2/h^2 - q_i) y_i + (-1/h^2 + p_i/(2h)) y_{i+1} = r_i
  const aDiag = new Array(N - 1).fill(0); // sub
  const bDiag = new Array(N - 1).fill(0); // main
  const cDiag = new Array(N - 1).fill(0); // super
  const d = new Array(N - 1).fill(0);

  for (let i = 1; i <= N - 1; i++) {
    const xi = x[i];
    const pi = p(xi);
    const qi = q(xi);
    const ri = r(xi);

    const Ai = -1 / (h * h) - pi / (2 * h);
    const Bi = 2 / (h * h) - qi;
    const Ci = -1 / (h * h) + pi / (2 * h);
    aDiag[i - 1] = Ai;
    bDiag[i - 1] = Bi;
    cDiag[i - 1] = Ci;

    let di = ri;
    if (i === 1) di -= Ai * alpha;
    if (i === N - 1) di -= Ci * beta;
    d[i - 1] = di;
  }

  // Thomas algorithm
  const { xSol: yInterior } = thomasSolve({ a: aDiag, b: bDiag, c: cDiag, d });
  const y = [alpha, ...yInterior, beta];
  return { x, y };
}

function thomasSolve({ a, b, c, d }) {
  const n = b.length;
  const cPrime = new Array(n);
  const dPrime = new Array(n);

  cPrime[0] = c[0] / b[0];
  dPrime[0] = d[0] / b[0];

  for (let i = 1; i < n; i++) {
    const denom = b[i] - a[i] * cPrime[i - 1];
    cPrime[i] = i < n - 1 ? c[i] / denom : 0;
    dPrime[i] = (d[i] - a[i] * dPrime[i - 1]) / denom;
  }

  const xSol = new Array(n);
  xSol[n - 1] = dPrime[n - 1];
  for (let i = n - 2; i >= 0; i--) {
    xSol[i] = dPrime[i] - cPrime[i] * xSol[i + 1];
  }
  return { xSol };
}

// Build Sturm–Liouville tridiagonal for -(p y')' + q y = λ w y on [a,b] with Dirichlet y(a)=y(b)=0
// Using uniform grid and centered differences for interior nodes.
function buildSLMatrix({ a, b, N, p, q, w }) {
  // N is number of subintervals; we have interior points M = N-1
  const h = (b - a) / N;
  const xs = Array.from({ length: N + 1 }, (_, i) => a + i * h);
  const M = N - 1;

  const A = new Array(M).fill(0); // sub
  const B = new Array(M).fill(0); // main
  const C = new Array(M).fill(0); // super
  const W = new Array(M).fill(0); // weight diagonal
  for (let i = 1; i <= N - 1; i++) {
    const xm = xs[i - 1];
    const xi = xs[i];
    const xp = xs[i + 1];

    const pL = p((xi + xm) / 2);
    const pR = p((xp + xi) / 2);

    const a_i = pL / (h * h);
    const c_i = pR / (h * h);
    const b_i = (-(pL + pR)) / (h * h) + q(xi);
    A[i - 1] = a_i;
    B[i - 1] = -b_i; // Move to the left-hand standard form A y = λ W y, we’ll keep signs consistent below
    C[i - 1] = c_i;
    W[i - 1] = w(xi);
  }
  // We are effectively returning (K y) = λ W y where K is SPD-like with entries:
  // K = tridiag(-A, B + A + C, -C) after consistent sign handling in solve step.
  // To avoid confusion, below we’ll assemble action Ky for iteration directly.
  return { xs, h, A, B, C, W };
}

// Apply K (stiffness-like) to vector y (interior only)
function applyK({ A, B, C, y }) {
  const n = y.length;
  const Ky = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    const left = i > 0 ? -A[i] * y[i - 1] : 0;
    const main = (A[i] + C[i] + B[i]) * y[i];
    const right = i < n - 1 ? -C[i] * y[i + 1] : 0;
    Ky[i] = left + main + right;
  }
  return Ky;
}

function dot(u, v) {
  let s = 0;
  for (let i = 0; i < u.length; i++) s += u[i] * v[i];
  return s;
}

function norm2(v) {
  return Math.sqrt(dot(v, v));
}

// Inverse iteration with shift σ, solving (K - σ W) z = W y_k, then normalize.
// We solve tridiagonal system via Thomas by explicitly forming (K - σ W).
function inverseIterationSL({ A, B, C, W, sigma = 0, iters = 30, tol = 1e-10, y0 = null }) {
  const n = W.length;
  let y = y0 ?? Array.from({ length: n }, () => Math.random() - 0.5);
  const ynorm = norm2(y);
  if (ynorm === 0) y[0] = 1;
  else y = y.map((v) => v / ynorm);

  let lambda = 0;

  for (let k = 0; k < iters; k++) {
    // Build tridiagonal (K - σ W)
    const a = new Array(n).fill(0);
    const b = new Array(n).fill(0);
    const c = new Array(n).fill(0);
    for (let i = 0; i < n; i++) {
      a[i] = i > 0 ? -A[i] : 0;
      b[i] = A[i] + C[i] + B[i] - sigma * W[i];
      c[i] = i < n - 1 ? -C[i] : 0;
    }

    // RHS = W y
    const d = y.map((vi, i) => W[i] * vi);

    const { xSol: z } = thomasSolve({ a, b, c, d });
    const zn = norm2(z);
    if (zn === 0) break;
    const yNext = z.map((zi) => zi / zn);

    // Rayleigh quotient λ ≈ (yᵀ K y) / (yᵀ W y)
    const Ky = applyK({ A, B, C, y: yNext });
    const num = dot(yNext, Ky);
    const den = dot(yNext, W.map((wi, i) => wi * yNext[i]));
    const lambdaNext = num / (den === 0 ? 1e-16 : den);

    if (Math.abs(lambdaNext - lambda) < tol) {
      lambda = lambdaNext;
      y = yNext;
      break;
    }
    lambda = lambdaNext;
    y = yNext;
  }

  return { lambda, y }; // interior y
}

// ------------------------------
// Section 27.1 — BVPs
// ------------------------------
const Section271 = () => {
  const [a, setA] = useState(0);
  const [b, setB] = useState(1);
  const [alpha, setAlpha] = useState(0);
  const [beta, setBeta] = useState(0.5);

  // Problem: y'' = -π^2 y + sin(π x), (simple linear)
  // Exact-ish reference (for demonstration we won't solve exact; focus on numerics)
  const [Nfd, setNfd] = useState(40);

  // p, q, r for linear BVP form y'' = p y' + q y + r
  const p = (x) => 0;
  const q = (x) => -Math.PI * Math.PI; // y'' = -π^2 y + sin(π x)
  const r = (x) => Math.sin(Math.PI * x);

  const fdSol = useMemo(() => {
    try {
      return fdBVP({ a, b, alpha, beta, N: Nfd, p, q, r });
    } catch (e) {
      return { x: [], y: [] };
    }
  }, [a, b, alpha, beta, Nfd]);

  // Shooting method for same problem by re-writing f(x,y,y')
  const f = (x, y, yp) => -Math.PI * Math.PI * y + Math.sin(Math.PI * x);

  const [Ns, setNs] = useState(200);
  const [s0, setS0] = useState(0.0);
  const [s1, setS1] = useState(1.0);

  const shoot = useMemo(() => {
    try {
      return shootingSecant({
        f,
        a,
        b,
        alpha,
        beta,
        n: Ns,
        s0,
        s1,
        maxIter: 15,
        tol: 1e-8,
      });
    } catch (e) {
      return { xs: [], ys: [], endError: NaN, slope: NaN };
    }
  }, [a, b, alpha, beta, Ns, s0, s1]);

  const layout1 = useMemo(
    () => ({
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(0,0,0,0)",
      margin: { l: 40, r: 10, t: 10, b: 40 },
      xaxis: { title: "x", color: "#a1a1aa", gridcolor: "#27272a" },
      yaxis: { title: "y", color: "#a1a1aa", gridcolor: "#27272a" },
      showlegend: false,
    }),
    []
  );

  return (
    <motion.div
      initial={{ opacity: 0, y: 8 }}
      whileInView={{ opacity: 1, y: 0 }}
      viewport={{ once: true }}
      transition={{ duration: 0.5 }}
    >
      <Card className="bg-zinc-900/60 border-zinc-700 rounded-2xl">
        <CardContent className="p-5 space-y-4">
          <SectionTitle icon={Route} k="27.1" title="General Methods for Boundary-Value Problems" />

          <div className="grid md:grid-cols-2 gap-6">
            <div className="space-y-3">
              <Chip tone="cyan">Model BVP</Chip>
             <Equation>
              {`Solve $y'' = -\\pi^2 y + \\sin(\\pi x)$, $x \\in [a,b]$, $y(a)=\\alpha$, $y(b)=\\beta$.`}
            </Equation>

              <Note>
                We demo two classical techniques: <strong>finite differences</strong> and the{" "}
                <strong>shooting method</strong>. Adjust the grid size and initial slopes to see
                convergence behavior. The example is deliberately simple but representative.
              </Note>

              <div className="grid sm:grid-cols-2 gap-3">
                <Field label="a">
                  <NumberInput value={a} onChange={(e) => setA(parseFloat(e.target.value))} />
                </Field>
                <Field label="b">
                  <NumberInput value={b} onChange={(e) => setB(parseFloat(e.target.value))} />
                </Field>
                <Field label="α = y(a)">
                  <NumberInput
                    value={alpha}
                    onChange={(e) => setAlpha(parseFloat(e.target.value))}
                  />
                </Field>
                <Field label="β = y(b)">
                  <NumberInput
                    value={beta}
                    onChange={(e) => setBeta(parseFloat(e.target.value))}
                  />
                </Field>
              </div>

              <div className="grid sm:grid-cols-2 gap-3">
                <Field label="Finite difference N (subintervals)" right="N ≥ 10">
                  <IntegerInput
                    min={10}
                    value={Nfd}
                    onChange={(e) => setNfd(parseInt(e.target.value || "10", 10))}
                  />
                </Field>
                <Field label="Shooting steps (n)">
                  <IntegerInput
                    min={20}
                    value={Ns}
                    onChange={(e) => setNs(parseInt(e.target.value || "200", 10))}
                  />
                </Field>
                <Field label="Shooting slope s₀">
                  <NumberInput value={s0} onChange={(e) => setS0(parseFloat(e.target.value))} />
                </Field>
                <Field label="Shooting slope s₁">
                  <NumberInput value={s1} onChange={(e) => setS1(parseFloat(e.target.value))} />
                </Field>
              </div>
            </div>

            <div className="space-y-3">
              <Chip tone="emerald">Methods</Chip>
              <Equation>
              {`Finite differences: $\\text{tridiagonal solve (Thomas)}$`}
            </Equation>
            <Equation>
              {`Shooting method: $\\text{secant on initial slope} + \\text{RK4 stepping}$`}
            </Equation>
       <div className="text-xs text-zinc-400">
                Shooting terminal error:{" "}
                <span className="text-zinc-100">
                  {isFinite(shoot.endError) ? shoot.endError.toExponential(3) : "—"}
                </span>
                , slope ≈{" "}
                <span className="text-zinc-100">
                  {isFinite(shoot.slope) ? shoot.slope.toFixed(6) : "—"}
                </span>
              </div>
            </div>
          </div>

          <div className="grid lg:grid-cols-2 gap-6">
            <Card className="bg-zinc-950/50 border-zinc-800">
              <CardContent className="p-4">
                <div className="flex items-center gap-2 mb-2">
                  <LineChart className="h-4 w-4 text-cyan-300" />
                  <div className="text-sm text-zinc-200">Finite Difference Solution</div>
                </div>
                <Plot
                  data={[
                    {
                      x: fdSol.x,
                      y: fdSol.y,
                      mode: "lines+markers",
                      marker: { size: 4 },
                      line: { width: 2 },
                      name: "FD",
                    },
                  ]}
                  layout={layout1}
                  config={{ displayModeBar: false, responsive: true }}
                  style={{ width: "100%", height: 320 }}
                />
              </CardContent>
            </Card>

            <Card className="bg-zinc-950/50 border-zinc-800">
              <CardContent className="p-4">
                <div className="flex items-center gap-2 mb-2">
                  <GitBranch className="h-4 w-4 text-emerald-300" />
                  <div className="text-sm text-zinc-200">Shooting Method Solution</div>
                </div>
                <Plot
                  data={[
                    {
                      x: shoot.xs,
                      y: shoot.ys,
                      mode: "lines",
                      line: { width: 2 },
                      name: "Shooting",
                    },
                  ]}
                  layout={layout1}
                  config={{ displayModeBar: false, responsive: true }}
                  style={{ width: "100%", height: 320 }}
                />
              </CardContent>
            </Card>
          </div>

          <div className="space-y-3">
            <Chip tone="amber">Observations</Chip>
            <Note>
              With a sufficiently fine grid (large <em>N</em>) the finite difference solution
              stabilizes. The shooting method is sensitive to initial slope guesses; the secant
              refinement typically converges quickly for well-posed linear BVPs.
            </Note>
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
};

// ------------------------------
// Section 27.2 — Eigenvalue Problems
// ------------------------------
const Section272 = () => {
  // Default classic: y'' + λ y = 0 on [0, π], y(0)=y(π)=0
  // Sturm–Liouville form: -(1*y')' + 0*y = λ*1*y
  const [a, setA] = useState(0);
  const [b, setB] = useState(Math.PI);
  const [N, setN] = useState(200); // fine grid
  const p = (x) => 1.0;
  const q = (x) => 0.0;
  const w = (x) => 1.0;

  const [sigma, setSigma] = useState(1.0); // shift near first eigenvalue (≈1^2 = 1)
  const [iters, setIters] = useState(30);

  const { xs, h, A, B, C, W } = useMemo(
    () => buildSLMatrix({ a, b, N, p, q, w }),
    [a, b, N]
  );

  const eig = useMemo(() => {
    try {
      return inverseIterationSL({ A, B, C, W, sigma, iters, tol: 1e-12 });
    } catch (e) {
      return { lambda: NaN, y: [] };
    }
  }, [A, B, C, W, sigma, iters]);

  const yFull = useMemo(() => [0, ...(eig.y ?? []), 0], [eig.y]);
  const xFull = xs;

  // Compare numerics to analytic λ_n = n^2
  const [kmax, setKmax] = useState(5);
  const spectrum = useMemo(() => {
    const arr = [];
    for (let k = 1; k <= kmax; k++) arr.push({ k, lambda: k * k });
    return arr;
  }, [kmax]);

  // Normalize eigenfunction to max |y| = 1 for display
  const yScale = useMemo(() => {
    const m = Math.max(...yFull.map((v) => Math.abs(v)));
    return m > 0 ? yFull.map((v) => v / m) : yFull;
  }, [yFull]);

  const layout = useMemo(
    () => ({
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(0,0,0,0)",
      margin: { l: 40, r: 10, t: 10, b: 40 },
      xaxis: { title: "x", color: "#a1a1aa", gridcolor: "#27272a" },
      yaxis: { title: "y", color: "#a1a1aa", gridcolor: "#27272a" },
      showlegend: false,
    }),
    []
  );

  return (
    <motion.div
      initial={{ opacity: 0, y: 8 }}
      whileInView={{ opacity: 1, y: 0 }}
      viewport={{ once: true }}
      transition={{ duration: 0.5 }}
    >
      <Card className="bg-zinc-900/60 border-zinc-700 rounded-2xl">
        <CardContent className="p-5 space-y-4">
          <SectionTitle icon={Sigma} k="27.2" title="Eigenvalue Problems" />

          <div className="grid md:grid-cols-2 gap-6">
            <div className="space-y-3">
              <Chip tone="cyan">Sturm–Liouville</Chip>
              <Equation>
  {`$-(p(x) y')' + q(x) y = \\lambda w(x) y$, $y(a)=0$, $y(b)=0$ (Dirichlet)`}
</Equation>
<Equation>
  {`Default demo: $p=1$, $q=0$, $w=1$; $[a,b]=[0,\\pi]$`}
</Equation>

              <Note>
                Use <strong>inverse iteration</strong> with shift σ to extract the eigenpair closest
                to σ. The classic problem yields eigenvalues λ<sub>n</sub>=n² with eigenfunctions
                sin(n x).
              </Note>
              <div className="grid sm:grid-cols-2 gap-3">
                <Field label="a">
                  <NumberInput value={a} onChange={(e) => setA(parseFloat(e.target.value))} />
                </Field>
                <Field label="b">
                  <NumberInput value={b} onChange={(e) => setB(parseFloat(e.target.value))} />
                </Field>
                <Field label="N (subintervals)">
                  <IntegerInput
                    min={20}
                    value={N}
                    onChange={(e) => setN(parseInt(e.target.value || "200", 10))}
                  />
                </Field>
                <Field label="Shift σ (target λ)">
                  <NumberInput
                    value={sigma}
                    onChange={(e) => setSigma(parseFloat(e.target.value))}
                  />
                </Field>
                <Field label="Inverse iteration steps">
                  <IntegerInput
                    min={5}
                    value={iters}
                    onChange={(e) => setIters(parseInt(e.target.value || "30", 10))}
                  />
                </Field>
                <Field label="Show analytic λ up to k=">
                  <IntegerInput
                    min={1}
                    value={kmax}
                    onChange={(e) => setKmax(parseInt(e.target.value || "5", 10))}
                  />
                </Field>
              </div>
            </div>

            <div className="space-y-3">
              <Chip tone="emerald">Results</Chip>
              <Equation
                latex={`Computed \\; \\lambda \\;\\approx\\; ${
                  isFinite(eig.lambda) ? eig.lambda.toFixed(6) : "—"
                }`}
              />
              <Note>
                Choose σ near the desired eigenvalue, e.g. <em>σ ≈ 1</em> for the first mode,{" "}
                <em>σ ≈ 4</em> for the second, etc. Increase N for better accuracy.
              </Note>
              <div className="text-xs text-zinc-400">
                Grid: <span className="text-zinc-100">{N}</span> subintervals, h ={" "}
                <span className="text-zinc-100">{((b - a) / N).toExponential(3)}</span>
              </div>
            </div>
          </div>

          <div className="grid lg:grid-cols-2 gap-6">
            <Card className="bg-zinc-950/50 border-zinc-800">
              <CardContent className="p-4">
                <div className="flex items-center gap-2 mb-2">
                  <LineChart className="h-4 w-4 text-cyan-300" />
                  <div className="text-sm text-zinc-200">Eigenfunction (normalized)</div>
                </div>
                <Plot
                  data={[
                    {
                      x: xFull,
                      y: yScale,
                      mode: "lines",
                      line: { width: 2 },
                      name: "φ(x)",
                    },
                  ]}
                  layout={layout}
                  config={{ displayModeBar: false, responsive: true }}
                  style={{ width: "100%", height: 320 }}
                />
              </CardContent>
            </Card>

            <Card className="bg-zinc-950/50 border-zinc-800">
              <CardContent className="p-4">
                <div className="flex items-center gap-2 mb-2">
                  <Landmark className="h-4 w-4 text-emerald-300" />
                  <div className="text-sm text-zinc-200">Analytic Spectrum λ<sub>k</sub>=k²</div>
                </div>
                <Plot
                  data={[
                    {
                      x: spectrum.map((s) => s.k),
                      y: spectrum.map((s) => s.lambda),
                      mode: "lines+markers",
                      marker: { size: 6 },
                      line: { width: 2 },
                      name: "k ↦ k²",
                    },
                    {
                      x: [1],
                      y: [eig.lambda],
                      mode: "markers",
                      marker: { size: 10, symbol: "star" },
                      name: "Computed λ (near σ)",
                    },
                  ]}
                  layout={{
                    ...layout,
                    xaxis: { title: "mode index k", color: "#a1a1aa", gridcolor: "#27272a" },
                    yaxis: { title: "λ", color: "#a1a1aa", gridcolor: "#27272a" },
                  }}
                  config={{ displayModeBar: false, responsive: true }}
                  style={{ width: "100%", height: 320 }}
                />
              </CardContent>
            </Card>
          </div>

          <div className="space-y-3">
            <Chip tone="amber">Notes</Chip>
            <Note>
              The finite-difference discretization transforms the continuous eigenvalue problem into
              a generalized matrix eigenproblem <em>K y = λ W y</em>. Our inverse iteration solves
              shifted linear systems to refine an eigenvector, with the Rayleigh quotient providing
              a sharp eigenvalue estimate.
            </Note>
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
};

// ------------------------------
// Section 27.3 — ODEs & Eigenvalues with Software Packages
// ------------------------------
const Section273 = () => {
  const [pkg, setPkg] = useState("Python / SciPy");
  const [text, setText] = useState(
    [
      "Tips when using software packages:",
      "• Prefer boundary-value solvers (e.g., collocation) over naive shooting for stiff problems.",
      "• Validate with mesh refinement studies (double N, check convergence rates).",
      "• For eigenvalues, compare against symmetry, orthogonality, and known asymptotics.",
      "• Exploit tridiagonal structure to reduce complexity from O(n^3) to O(n).",
      "• Always nondimensionalize and scale to improve conditioning.",
    ].join("\n")
  );

  return (
    <motion.div
      initial={{ opacity: 0, y: 8 }}
      whileInView={{ opacity: 1, y: 0 }}
      viewport={{ once: true }}
      transition={{ duration: 0.5 }}
    >
      <Card className="bg-zinc-900/60 border-zinc-700 rounded-2xl">
        <CardContent className="p-5 space-y-4">
          <SectionTitle icon={Settings} k="27.3" title="ODEs and Eigenvalues with Software Packages" />

          <div className="grid lg:grid-cols-2 gap-6">
            <div className="space-y-3">
              <Chip tone="cyan">Package Hints</Chip>
          <Equation>
  {`Prefer collocation-based BVP solvers (e.g., $\\text{solve\\_bvp}$) for robust convergence.`}
</Equation>



              <Note>
                In practice, large-scale BVPs and eigenproblems benefit from <em>preconditioned</em>{" "}
                Krylov methods and mesh adaptivity. Packages provide hardened linear algebra and
                robust step-size control; your job is to nondimensionalize, scale, and validate.
              </Note>
              <div className="grid sm:grid-cols-2 gap-3">
                <Field label="Ecosystem / Package">
                  <Input
                    value={pkg}
                    onChange={(e) => setPkg(e.target.value)}
                    className="bg-zinc-800/90 text-zinc-100 border-zinc-700 focus-visible:ring-cyan-600"
                  />
                </Field>
              </div>
            </div>

            <div className="space-y-3">
              <Chip tone="emerald">Playbook</Chip>
              <Textarea
                rows={10}
                value={text}
                onChange={(e) => setText(e.target.value)}
                className="bg-zinc-800/90 text-zinc-100 border-zinc-700 focus-visible:ring-cyan-600"
              />
            </div>
          </div>

          <div className="space-y-3">
            <Chip>Checks</Chip>
            <Note>
              <ul className="list-disc ml-5 space-y-1">
                <li>
                  <strong>Mesh refinement:</strong> verify error decreases at expected order.
                </li>
                <li>
                  <strong>Conservation / invariants:</strong> for Hamiltonian or self-adjoint forms.
                </li>
                <li>
                  <strong>Orthogonality:</strong> eigenfunctions should be W-orthogonal.
                </li>
                <li>
                  <strong>Sensitivity:</strong> perturb data/BCs to assess conditioning.
                </li>
              </ul>
            </Note>
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
};

// ------------------------------
// Problems
// ------------------------------
const Problems27 = () => {
  const [open, setOpen] = useState([false, false, false, false, false]);

  return (
    <motion.div
      initial={{ opacity: 0, y: 8 }}
      whileInView={{ opacity: 1, y: 0 }}
      viewport={{ once: true }}
      transition={{ duration: 0.5 }}
    >
      <Card className="bg-zinc-900/60 border-zinc-700 rounded-2xl">
        <CardContent className="p-5 space-y-4">
          <SectionTitle icon={ListChecks} k="Problems" title="Practice & Exploration" />

          <div className="space-y-6">
            {/* P1 */}
            <div className="space-y-2">
              <div className="flex items-center gap-2">
                <Dot className="h-5 w-5 text-cyan-300" />
                <div className="font-semibold text-zinc-100">P1. Linear BVP (FD)</div>
              </div>
              <Note>
                Solve y'' = 4y - 3x with y(0)=0 and y(1)=1 using finite differences. Perform a mesh
                refinement study with N=20, 40, 80.
              </Note>
              <Toggle
                open={open[0]}
                onClick={() =>
                  setOpen((s) => {
                    const t = [...s];
                    t[0] = !t[0];
                    return t;
                  })
                }
              >
                Show solution sketch
              </Toggle>
              {open[0] && (
                <div className="text-sm text-zinc-300 bg-zinc-950/40 border border-zinc-800 rounded-xl p-3">
                  Discretize <em>y'' = p y' + q y + r</em> with p=0, q=-4, r=3x and Dirichlet BCs.
                  Assemble tridiagonal and solve via Thomas. Observe second-order convergence.
                </div>
              )}
            </div>

            {/* P2 */}
            <div className="space-y-2">
              <div className="flex items-center gap-2">
                <Dot className="h-5 w-5 text-emerald-300" />
                <div className="font-semibold text-zinc-100">P2. Shooting Stability</div>
              </div>
              <Note>
                For y'' = -100 y + e<sup>-x</sup>, y(0)=0, y(1)=0, test shooting with different
                initial slopes and discuss stiffness.
              </Note>
              <Toggle
                open={open[1]}
                onClick={() =>
                  setOpen((s) => {
                    const t = [...s];
                    t[1] = !t[1];
                    return t;
                  })
                }
              >
                Show solution sketch
              </Toggle>
              {open[1] && (
                <div className="text-sm text-zinc-300 bg-zinc-950/40 border border-zinc-800 rounded-xl p-3">
                  The problem is moderately stiff; naive shooting can be ill-conditioned. Use finite
                  differences or collocation. If shooting, apply robust root-finding and step-size
                  control with RK methods adapted to stiffness.
                </div>
              )}
            </div>

            {/* P3 */}
            <div className="space-y-2">
              <div className="flex items-center gap-2">
                <Dot className="h-5 w-5 text-amber-300" />
                <div className="font-semibold text-zinc-100">P3. Sturm–Liouville Weight</div>
              </div>
              <Note>
                Consider <em>-(p y')' + q y = λ w y</em> on [0,1] with p=1, q=0 and w(x)=1+x. Use
                the eigen-solver to approximate the first three eigenvalues by changing the shift.
              </Note>
              <Toggle
                open={open[2]}
                onClick={() =>
                  setOpen((s) => {
                    const t = [...s];
                    t[2] = !t[2];
                    return t;
                  })
                }
              >
                Show solution sketch
              </Toggle>
              {open[2] && (
                <div className="text-sm text-zinc-300 bg-zinc-950/40 border border-zinc-800 rounded-xl p-3">
                  Modify the weight vector W accordingly. Sweep σ≈λ<sub>1</sub>, then σ≈λ
                  <sub>2</sub>, etc. Normalize eigenvectors and check W-orthogonality via inner
                  product &Sigma; w<sub>i</sub> y<sub>i</sub> z<sub>i</sub>.
                </div>
              )}
            </div>

            {/* P4 */}
            <div className="space-y-2">
              <div className="flex items-center gap-2">
                <Dot className="h-5 w-5 text-cyan-300" />
                <div className="font-semibold text-zinc-100">P4. Robin BCs</div>
              </div>
              <Note>
                Adapt the finite-difference formulation for mixed (Robin) boundary conditions. Test
                with α y(0) + β y'(0) = γ and y(1)=0.
              </Note>
              <Toggle
                open={open[3]}
                onClick={() =>
                  setOpen((s) => {
                    const t = [...s];
                    t[3] = !t[3];
                    return t;
                  })
                }
              >
                Show implementation hint
              </Toggle>
              {open[3] && (
                <div className="text-sm text-zinc-300 bg-zinc-950/40 border border-zinc-800 rounded-xl p-3">
                  Replace the first interior equation with a ghost-point relation to encode the
                  Robin condition: y'(0) ≈ (y<sub>1</sub> - y<sub>-1</sub>)/(2h). Eliminate the
                  ghost via the BC to modify the first row coefficients.
                </div>
              )}
            </div>

            {/* P5 */}
            <div className="space-y-2">
              <div className="flex items-center gap-2">
                <Dot className="h-5 w-5 text-emerald-300" />
                <div className="font-semibold text-zinc-100">P5. Conditioning Study</div>
              </div>
              <Note>
                Examine how the condition number of the tridiagonal system scales with N for a
                second-order BVP. Discuss implications for floating-point error.
              </Note>
              <Toggle
                open={open[4]}
                onClick={() =>
                  setOpen((s) => {
                    const t = [...s];
                    t[4] = !t[4];
                    return t;
                  })
                }
              >
                Show discussion point
              </Toggle>
              {open[4] && (
                <div className="text-sm text-zinc-300 bg-zinc-950/40 border border-zinc-800 rounded-xl p-3">
                  As N↑, the discretized operator broadens the spectral span. Conditioning worsens,
                  amplifying roundoff. Preconditioning and scaling (and using SPD-aware solvers)
                  mitigate error growth.
                </div>
              )}
            </div>
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
};

// ------------------------------
// Doc Panel (like your ChapterDocs)
// ------------------------------
const ChapterDocs = () => {
  const [active, setActive] = useState(0);
  const sections = [
    {
      h: "27.1 General Methods for Boundary-Value Problems",
      body: [
        "• Shooting method transforms BVP to IVP; root-find on the missing slope.",
        "• Finite differences yield tridiagonal linear systems solved efficiently by Thomas.",
        "• For nonlinear right-hand sides, use Newton linearization on the discrete system.",
      ],
    },
    {
      h: "27.2 Eigenvalue Problems",
      body: [
        "• Sturm–Liouville form K y = λ W y on a grid with Dirichlet/Neumann/Robin BCs.",
        "• Use inverse iteration or Lanczos/LOBPCG for large-scale problems.",
        "• Check orthogonality in the W-inner-product and monitor Rayleigh quotient residual.",
      ],
    },
    {
      h: "27.3 Packages & Practice",
      body: [
        "• Prefer collocation BVP solvers and sparse eigensolvers (ARPACK/LOBPCG).",
        "• Validate with convergence rates and symmetry/energy checks.",
        "• Nondimensionalize to tame disparate scales.",
      ],
    },
  ];

  return (
    <Card className="bg-zinc-900/60 border-zinc-700 rounded-2xl">
      <CardContent className="p-5">
        <div className="flex items-center gap-2 mb-4">
          <BookOpen className="h-5 w-5 text-cyan-300" />
          <h3 className="text-zinc-100 font-semibold">Chapter Notes</h3>
        </div>

        <div className="grid md:grid-cols-3 gap-4">
          <div className="md:col-span-1">
            <div className="space-y-2">
              {sections.map((s, i) => (
                <button
                  key={i}
                  onClick={() => setActive(i)}
                  className={`w-full text-left px-3 py-2 rounded-xl border ${
                    i === active
                      ? "bg-cyan-950/50 cursor-pointer border-cyan-700 text-cyan-200"
                      : "bg-zinc-900/60 cursor-pointer border-zinc-700 text-zinc-300 hover:border-cyan-700"
                  }`}
                >
                  <div className="text-sm">{s.h}</div>
                </button>
              ))}
            </div>
          </div>

          <div className="md:col-span-2">
            <div className="bg-zinc-950/50 border border-zinc-800 rounded-2xl p-4 space-y-2">
              <div className="text-zinc-100 font-medium">{sections[active].h}</div>
              <ul className="list-disc ml-5 text-sm text-zinc-300 space-y-1">
                {sections[active].body.map((line, idx) => (
                  <li key={idx}>{line}</li>
                ))}
              </ul>
              <div className="text-xs text-zinc-500 mt-2">
                Tip: click other headings to switch notes.
              </div>
            </div>
          </div>
        </div>
      </CardContent>
    </Card>
  );
};

// ------------------------------
// Chapter header
// ------------------------------
const ChapterHeader = () => {
  return (
    <div className="rounded-2xl bg-gradient-to-b from-zinc-900/70 to-zinc-950/50 p-6">
      <div className="flex flex-col sm:flex-row sm:items-center sm:justify-between gap-4">
        <div>
      
          <h1 className="text-2xl sm:text-3xl font-bold text-emerald-400">
            Boundary-Value and Eigenvalue Problems
          </h1>
          <div className="text-zinc-400 mt-1">27.1 • 27.2 • 27.3 • Problems</div>
        </div>
       
      </div>
    </div>
  );
};

// ------------------------------
// Main Chapter component
// ------------------------------
const Chapter27 = () => {
  useEffect(() => {
    // Scroll to top when mounted (optional)
    window.scrollTo({ top: 0, behavior: "smooth" });
  }, []);

  return (
    <div className="p-6 space-y-6 bg-zinc-950">
      <ChapterHeader />
      <Section271 />
      <Section272 />
      <Section273 />
      <Problems27 />
      <ChapterDocs />
      <BottomBar/>
    </div>
  );
};

export default Chapter27;
