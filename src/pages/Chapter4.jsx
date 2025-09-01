// src/pages/Chapter4_TruncationAndTaylor.jsx
import React, { useState, useMemo, useEffect, useRef } from "react";
import { motion } from "framer-motion";
import Plot from "react-plotly.js";
import { Canvas, useFrame } from "@react-three/fiber";
import { OrbitControls } from "@react-three/drei";
import * as THREE from "three";
import {
  BookOpen,
  ChevronRight,
  ChevronDown,
  Triangle,
  Zap,
  RotateCw,
  LineChart,
  AlertTriangle,
  List,
  Layers,
  Info,
  GitMerge,
} from "lucide-react";

import { Card, CardHeader, CardTitle, CardContent } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import BottomBar from "../components/BottomBar";

// ---------------- Theme & motion ----------------
const theme = {
  background: "bg-zinc-900",
  panelBg: "#0f1720",
  text: "#e6eef6",
  accent1: "#60a5fa",
  accent2: "#34d399",
  accent3: "#f59e0b",
};

const cardMotion = {
  initial: { opacity: 0, y: 10, scale: 0.995 },
  enter: { opacity: 1, y: 0, scale: 1, transition: { duration: 0.36 } },
  hover: { scale: 1.02, y: -6, transition: { duration: 0.16 } },
};

// ---------------- Numerics helpers ----------------
function factorial(n) {
  let r = 1;
  for (let i = 2; i <= n; i++) r *= i;
  return r;
}

function binomialFraction(alpha, k) {
  if (k === 0) return 1;
  let num = 1;
  for (let i = 0; i < k; i++) num *= alpha - i;
  return num / factorial(k);
}

function taylorTerm(func, a, k, x) {
  const h = x - a;
  switch (func) {
    case "sin":
      return (Math.pow(-1, k) * Math.pow(h, 2 * k + 1)) / factorial(2 * k + 1);
    case "cos":
      return (Math.pow(-1, k) * Math.pow(h, 2 * k)) / factorial(2 * k);
    case "exp":
      return Math.pow(h, k) / factorial(k);
    case "ln1p":
      if (k === 0) return 0;
      return (Math.pow(-1, k + 1) * Math.pow(h, k)) / k;
    case "sqrt1p":
      return binomialFraction(0.5, k) * Math.pow(h, k);
    default:
      return 0;
  }
}

function taylorPoly(func, a, nTerms, x) {
  let s = 0;
  for (let k = 0; k < nTerms; k++) s += taylorTerm(func, a, k, x);
  return s;
}

function finiteDiffCentral(f, x, h) {
  return (f(x + h) - f(x - h)) / (2 * h);
}

function finiteDiffForward(f, x, h) {
  return (f(x + h) - f(x)) / h;
}

function kahanSum(arr) {
  let sum = 0.0,
    c = 0.0;
  for (let i = 0; i < arr.length; i++) {
    let y = arr[i] - c;
    let t = sum + y;
    c = (t - sum) - y;
    sum = t;
  }
  return sum;
}

// ---------------- Small 3D visual element ----------------
function SmallSurface({ color = theme.accent1 }) {
  const ref = useRef();
  useFrame(({ clock }) => {
    if (!ref.current) return;
    ref.current.rotation.y = clock.getElapsedTime() * 0.12;
    ref.current.rotation.x = Math.sin(clock.getElapsedTime() * 0.06) * 0.08;
  });

  const pts = useMemo(() => {
    const n = 20;
    const arr = [];
    for (let i = 0; i < n; i++) {
      const u = (i / (n - 1)) * Math.PI * 2;
      for (let j = 0; j < n; j++) {
        const v = (j / (n - 1)) * 2 - 1;
        const x = Math.cos(u) * (1 + 0.2 * v);
        const y = Math.sin(u) * (1 + 0.2 * v);
        const z = 0.4 * Math.sin(3 * u) * v;
        arr.push(x, y, z);
      }
    }
    return new Float32Array(arr);
  }, []);

  return (
    <points ref={ref}>
      <bufferGeometry>
        <bufferAttribute attach="attributes-position" array={pts} count={pts.length / 3} itemSize={3} />
      </bufferGeometry>
      <pointsMaterial size={0.05} sizeAttenuation color={color} />
    </points>
  );
}

// ---------------- Expanded documentation content ----------------
const docsText = {
  "4.1 The Taylor Series": [
    `The Taylor series expansion of a sufficiently smooth function f about the point a represents f locally as a polynomial plus a remainder:
      
f(x) = Σ_{n=0}^∞ (f^{(n)}(a)/n!) (x - a)^n.

The practical interest is the truncated polynomial P_N(x) = Σ_{n=0}^N (f^{(n)}(a)/n!) (x - a)^n and the remainder R_N(x) = f(x) - P_N(x). 

Lagrange form of the remainder:
R_N(x) = f^{(N+1)}(ξ) / (N+1)! * (x - a)^{N+1} for some ξ between a and x.

Key points:
- The region where the Taylor series converges equals the distance to the nearest singularity (complex analysis perspective) — for practical computation, this controls where Taylor approximations are accurate.
- When using Taylor expansions to derive finite-difference formulas, quadrature rules, or local approximations, the order of the leading neglected term determines the local truncation order.
- For global polynomial interpolation, plain Taylor polynomials can perform poorly away from a: use Chebyshev / minimax or piecewise polynomials instead.`,
    `Analytic examples and special cases:
- e^x = Σ x^n / n! (entire function; converges for all x).
- sin x = Σ (-1)^k x^{2k+1} / (2k+1)!, cos x = Σ (-1)^k x^{2k} / (2k)!.
- Binomial: (1 + x)^α = Σ C(α, k) x^k for |x| < 1; used for fractional powers and perturbation expansions.
- Practical usage: derive error bounds using the Lagrange remainder or alternating-series theorem (when applicable).`,
    `Practical guidance:
- Always compute or bound the remainder when reporting truncated series results.
- For moderate intervals, prefer polynomial approximations optimized on the interval (Chebyshev) rather than centered Taylor polynomials.
- Use automatic differentiation or symbolic differentiation for reliable high-order derivative evaluation.`,
  ],
  "4.2 Error Propagation": [
    `First-order linear propagation:
For y = f(x) with small input perturbation Δx, linearization gives:
Δy ≈ f'(x) Δx.

In relative terms (if x and y nonzero):
|Δy / y| ≈ |(x f'(x) / f(x))| * |Δx / x| which motivates the condition number:
κ(x) = |x f'(x) / f(x)|.

Multivariate propagation:
For y = f(x1,...,xn), small input covariance Σ_x propagates as:
Cov(y) ≈ J Σ_x J^T where J is the Jacobian matrix. This generalizes the scalar derivative to linear mappings and is essential for uncertainty quantification in engineering workflows.`,
    `Monte Carlo vs linearization:
- Linearization is fast and gives a first-order estimate; it fails when Δx are not small or f is strongly nonlinear within the uncertainty region.
- Monte Carlo sampling is robust and captures nonlinear and non-Gaussian effects; it is more expensive but offers empirical distributions of outputs (mean, variance, quantiles).`,
    `Stability & conditioning:
- Conditioning: sensitivity of the exact mathematical problem to perturbations in data (intrinsic).
- Stability: behavior of the numerical algorithm when executed in finite precision.
A well-conditioned problem with an unstable algorithm can produce inaccurate results; conversely, a stable algorithm may still be limited by a badly conditioned problem.`,
  ],
  "4.3 Total Numerical Error": [
    `Total numerical error typically decomposes into:
E_total = E_truncation + E_roundoff + E_model + E_other.

Truncation (discretization) error:
- Due to replacing infinite processes with finite approximations (Taylor truncation, finite difference discretization, time-stepping).
- Often scales like a power of the discretization parameter (h^p), where p is the method order.

Round-off error:
- Each floating-point operation introduces relative error O(ε) (machine epsilon; ~2.22e-16 for double).
- Accumulation, catastrophic cancellation, and subtractions of nearly equal numbers amplify round-off.

Balancing errors:
- Example: derivative estimates.
  - Forward diff: truncation ∼ C1 h, round-off ∼ C2 ε / h ⇒ total ∼ C1 h + C2 ε / h → h_opt ∝ sqrt(ε).
  - Central diff: truncation ∼ C3 h^2, round-off ∼ C2 ε / h ⇒ total ∼ C3 h^2 + C2 ε / h → h_opt ∝ ε^{1/3}.
Selecting h in practice benefits from running an error vs h curve (log-log) and finding the U-shaped minimum.`,
    `Mitigation strategies:
- Use compensated summation (Kahan) for long sums.
- Use higher-order methods to reduce truncation while keeping step moderate.
- Use extended precision for critical computations.
- Design algorithms that minimize subtractive cancellation (algebraic rearrangement).`,
  ],
  "4.4 Blunders & Uncertainty": [
    `Non-numerical errors often dominate a project:
- Blunders: typos, wrong sign, unit mistakes.
- Formulation errors: missing physics, wrong boundary conditions, incorrect assumptions.
- Data uncertainty: sensor noise, bias, poor calibration.

Detection & mitigation:
- Sanity checks: dimension analysis, monotonicity, conservation laws.
- Unit tests & regression tests: verify algorithms against known solutions.
- Validation: compare with experimental or benchmark data.
- Report uncertainty: provide error bars, confidence intervals, and clearly state assumptions.`,
    `Organizational practices:
- Keep small reproducible examples for each algorithm.
- Automate unit checks in CI to catch regressions.
- Use reproducible pipelines and document data provenance.`,
  ],
};

// ---------------- Reusable small card ----------------
function SummaryCard({ title, value, subtitle }) {
  return (
    <div className="bg-zinc-800 overflow-auto p-3 rounded min-w-[160px]">
      <div className="text-zinc-400 text-xs">{title}</div>
      <div className="text-zinc-100 font-mono text-lg mt-1">{value}</div>
      {subtitle ? <div className="text-zinc-500 text-xs mt-1">{subtitle}</div> : null}
    </div>
  );
}

// ---------------- TABS (simple local implementation) ----------------
function SimpleTabs({ tabs, value, onChange }) {
  return (
    <div>
      <div className="flex gap-2 mb-3">
        {tabs.map((t) => (
          <button
            key={t.key}
            onClick={() => onChange(t.key)}
            className={`px-3 py-1 rounded ${value === t.key ? "bg-zinc-700 text-white" : "bg-zinc-800 text-zinc-300"}`}
          >
            {t.label}
          </button>
        ))}
      </div>
    </div>
  );
}

// ---------------- Taylor Series Panel ----------------
function TaylorSeriesPanel() {
  const [func, setFunc] = useState("sin");
  const [a, setA] = useState(0);
  const [n, setN] = useState(5);
  const [xVal, setXVal] = useState(Math.PI / 3);
  const [range, setRange] = useState(Math.PI);
  const [autoRange, setAutoRange] = useState(true);

  const plot = useMemo(() => {
    const xs = [];
    const yExact = [];
    const yApprox = [];
    const steps = 240;
    const half = range;
    const left = a - half;
    const right = a + half;
    const dx = (right - left) / steps;
    for (let i = 0; i <= steps; i++) {
      const x = left + i * dx;
      xs.push(x);
      let exact;
      switch (func) {
        case "sin":
          exact = Math.sin(x);
          break;
        case "cos":
          exact = Math.cos(x);
          break;
        case "exp":
          exact = Math.exp(x);
          break;
        case "ln1p":
          exact = Math.log(1 + x);
          break;
        case "sqrt1p":
          exact = Math.sqrt(1 + x);
          break;
        default:
          exact = Math.sin(x);
      }
      yExact.push(exact);
      yApprox.push(taylorPoly(func, a, n, x));
    }
    const exactAtX = (() => {
      try {
        switch (func) {
          case "sin":
            return Math.sin(xVal);
          case "cos":
            return Math.cos(xVal);
          case "exp":
            return Math.exp(xVal);
          case "ln1p":
            return Math.log(1 + xVal);
          case "sqrt1p":
            return Math.sqrt(1 + xVal);
          default:
            return Math.sin(xVal);
        }
      } catch {
        return NaN;
      }
    })();
    const approxAtX = taylorPoly(func, a, n, xVal);
    const truncError = exactAtX - approxAtX;
    return { xs, yExact, yApprox, exactAtX, approxAtX, truncError, left, right };
  }, [func, a, n, xVal, range]);

  return (
    <motion.div variants={cardMotion} initial="initial" animate="enter" whileHover="hover">
      <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700 overflow-hidden`}>
        <CardHeader>
          <CardTitle className="flex items-center gap-2 text-blue-300">
            <Triangle className="w-5 h-5" /> 4.1 The Taylor Series — interactive
          </CardTitle>
        </CardHeader>
        <CardContent>
          {/* Controls */}
          <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-3 mb-3">
            <div>
              <label className="text-zinc-200 text-xs">Function</label>
              <select value={func} onChange={(e) => setFunc(e.target.value)} className="w-full bg-zinc-700 text-white p-2 rounded">
                <option value="sin">sin(x)</option>
                <option value="cos">cos(x)</option>
                <option value="exp">exp(x)</option>
                <option value="ln1p">ln(1 + x)</option>
                <option value="sqrt1p">sqrt(1 + x)</option>
              </select>
            </div>

            <div>
              <label className="text-zinc-200 text-xs">Expansion point a</label>
              <Input value={String(a)} onChange={(e) => setA(Number(e.target.value))} className="bg-zinc-700 text-white" />
            </div>

            <div>
              <label className="text-zinc-200 text-xs">Terms N</label>
              <div className="flex items-center gap-3">
                <input
                  type="range"
                  min={1}
                  max={50}
                  value={n}
                  onChange={(e) => setN(Number(e.target.value))}
                  className="w-full"
                />
                <div className="text-zinc-100 font-mono w-12 text-right">{n}</div>
              </div>
            </div>

            <div>
              <label className="text-zinc-200 text-xs">Evaluate at x</label>
              <Input value={String(xVal)} onChange={(e) => setXVal(Number(e.target.value))} className="bg-zinc-700 text-white" />
            </div>
          </div>

          <div className="grid grid-cols-1 sm:grid-cols-3 gap-3 items-center mb-4">
            <div className="flex items-center gap-2">
              <input id="autoRange" type="checkbox" checked={autoRange} onChange={(e) => setAutoRange(e.target.checked)} className="accent-sky-400" />
              <label htmlFor="autoRange" className="text-zinc-300 text-sm">Auto-range</label>
            </div>

            <div className="sm:col-span-2">
              <label className="text-zinc-200 text-xs mb-1">Plot half-width (±)</label>
              <div className="flex gap-2 items-center">
                <input type="range" min={0.1} max={Math.PI * 4} step={0.01} value={range} onChange={(e) => setRange(Number(e.target.value))} className="w-full" />
                <div className="w-20 text-right text-zinc-100 font-mono">{range.toFixed(2)}</div>
              </div>
            </div>
          </div>

          {/* Plot area (constrained) */}
          <div className="mb-3 overflow-hidden rounded">
            <div className="w-full" style={{ minHeight: 300 }}>
              <Plot
                data={[
                  { x: plot.xs, y: plot.yExact, mode: "lines", name: "Exact", line: { color: theme.accent2 } },
                  { x: plot.xs, y: plot.yApprox, mode: "lines", name: `Taylor (N=${n})`, line: { color: theme.accent1, dash: "dash" } },
                  { x: [xVal], y: [plot.exactAtX], mode: "markers", name: "Exact @ x", marker: { color: theme.accent2, size: 8 } },
                  { x: [xVal], y: [plot.approxAtX], mode: "markers", name: "Approx @ x", marker: { color: theme.accent1, size: 8 } },
                ]}
                layout={{
                  title: `Taylor approx of ${func} about a=${a}`,
                  paper_bgcolor: theme.panelBg,
                  plot_bgcolor: theme.panelBg,
                  font: { color: theme.text },
                  height: 360,
                  margin: { t: 40, l: 50, r: 20, b: 50 },
                  xaxis: { title: "x", range: [plot.left, plot.right] },
                  yaxis: { title: "f(x)" },
                  legend: { orientation: "h", y: -0.2 },
                }}
                useResizeHandler={true}
                style={{ width: "100%" }}
              />
            </div>
          </div>

          <div className="grid grid-cols-1 sm:grid-cols-3 gap-3">
            <SummaryCard title="Exact at x" value={String(plot.exactAtX)} />
            <SummaryCard title={`Approx at x (N=${n})`} value={String(plot.approxAtX)} />
            <SummaryCard title="Truncation error" value={String(plot.truncError)} subtitle="exact - approx" />
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ---------------- Error Propagation Panel (native selects and range inputs) ----------------
function ErrorPropagationPanel() {
  const [tab, setTab] = useState("linear");
  const [expr, setExpr] = useState("Math.sin(x)");
  const [x0, setX0] = useState(0.5);
  const [dx, setDx] = useState(0.01);
  const [samples, setSamples] = useState(500);

  const safeEval = (fStr) => {
    try {
      // eslint-disable-next-line no-new-func
      return new Function("x", `with (Math) { return (${fStr}); }`);
    } catch {
      return null;
    }
  };

  const f = useMemo(() => safeEval(expr), [expr]);

  const derivative = useMemo(() => {
    if (!f) return NaN;
    const h = 1e-6;
    try {
      return (f(x0 + h) - f(x0 - h)) / (2 * h);
    } catch {
      return NaN;
    }
  }, [f, x0]);

  const monte = useMemo(() => {
    if (!f) return { xs: [], ys: [] };
    const xs = [];
    const ys = [];
    for (let i = 0; i < samples; i++) {
      const eps = (Math.random() * 2 - 1) * dx;
      const xv = x0 + eps;
      try {
        xs.push(xv);
        ys.push(f(xv));
      } catch {
        xs.push(xv);
        ys.push(NaN);
      }
    }
    return { xs, ys };
  }, [f, x0, dx, samples]);

  const stats = useMemo(() => {
    const arr = (monte.ys || []).filter((v) => Number.isFinite(v));
    if (!arr.length) return null;
    const mean = arr.reduce((a, b) => a + b, 0) / arr.length;
    const variance = arr.reduce((s, v) => s + (v - mean) ** 2, 0) / (arr.length - 1 || 1);
    return { mean, std: Math.sqrt(variance), n: arr.length };
  }, [monte]);

  const tabs = [
    { key: "linear", label: "Linearization" },
    { key: "montecarlo", label: "Monte Carlo" },
  ];

  return (
    <motion.div variants={cardMotion} initial="initial" animate="enter" whileHover="hover">
      <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700 overflow-hidden`}>
        <CardHeader>
          <CardTitle className="flex items-center gap-2 text-emerald-300">
            <RotateCw className="w-5 h-5" /> 4.2 Error Propagation — Linearization & Monte Carlo
          </CardTitle>
        </CardHeader>
        <CardContent>
          <SimpleTabs tabs={tabs} value={tab} onChange={setTab} />

          {tab === "linear" && (
            <div>
              <div className="grid grid-cols-1 md:grid-cols-3 gap-3 mb-3">
                <div className="md:col-span-2">
                  <label className="text-zinc-200 text-xs">Function (JS expression)</label>
                  <Input value={expr} onChange={(e) => setExpr(e.target.value)} className="bg-zinc-700 text-white" />
                  <div className="text-zinc-400 text-xs mt-1">Example: Math.sin(x) or Math.exp(x) * (1 + x)</div>
                </div>

                <div className="flex flex-col gap-2">
                  <label className="text-zinc-200 text-xs">Point x₀</label>
                  <Input type="number" value={x0} onChange={(e) => setX0(Number(e.target.value))} className="bg-zinc-700 text-white" />
                  <label className="text-zinc-200 text-xs">Δx (uncertainty)</label>
                  <div className="flex items-center gap-2">
                    <input type="range" min={1e-6} max={1} step={1e-4} value={dx} onChange={(e) => setDx(Number(e.target.value))} className="w-full" />
                    <div className="w-24 text-right text-zinc-100 font-mono">{dx.toExponential(2)}</div>
                  </div>
                </div>
              </div>

              <div className="grid grid-cols-1 md:grid-cols-3 gap-3 mb-3">
                <SummaryCard title="Numeric f'(x₀)" value={Number(derivative).toExponential(3)} />
                <SummaryCard title="Linear Δy estimate" value={Number(derivative * dx).toExponential(3)} />
                <SummaryCard title="Notes" value="Linearization uses first derivative only" subtitle="Small Δx recommended" />
              </div>

              <div className="mt-3 text-zinc-400 text-sm">
                Linearization provides a first-order estimate. When f is nonlinear across the uncertainty region, this estimate can be optimistic. Use Monte Carlo to validate.
              </div>
            </div>
          )}

          {tab === "montecarlo" && (
            <div>
              <div className="grid grid-cols-1 md:grid-cols-3 gap-3 mb-3">
                <div className="md:col-span-2">
                  <label className="text-zinc-200 text-xs">Samples</label>
                  <div className="flex items-center gap-2">
                    <input type="range" min={10} max={2000} step={10} value={samples} onChange={(e) => setSamples(Number(e.target.value))} className="w-full" />
                    <div className="w-20 text-right text-zinc-100 font-mono">{samples}</div>
                  </div>
                </div>

                <div>
                  <label className="text-zinc-200 text-xs">Δx preview</label>
                  <div className="text-zinc-100 font-mono">{dx.toExponential(2)}</div>
                </div>
              </div>

              <div className="grid grid-cols-1 md:grid-cols-2 gap-3 mb-3">
                <div className="bg-zinc-800 p-3 rounded overflow-auto">
                  <div className="text-zinc-400 text-xs">Monte Carlo samples (x, f(x))</div>
                  <div className="text-zinc-100 font-mono text-sm mt-2 max-h-40 overflow-auto">
                    {monte.xs.slice(0, 20).map((xv, i) => (
                      <div key={i}>{xv.toFixed(6)} → {String(monte.ys[i]).slice(0, 12)}</div>
                    ))}
                    {monte.xs.length > 20 ? <div className="text-zinc-500 text-xs">... ({monte.xs.length} samples)</div> : null}
                  </div>
                </div>

                <div className="bg-zinc-800 p-3 rounded">
                  <div className="text-zinc-400 text-xs">Empirical stats</div>
                  <div className="text-zinc-100 font-mono text-lg mt-2">{stats ? `μ=${stats.mean.toFixed(6)}, σ=${stats.std.toExponential(2)}` : "—"}</div>
                </div>
              </div>

              <div className="mb-3">
                <div className="w-full overflow-auto rounded">
                  <Plot
                    data={[
                      { x: monte.xs, y: monte.ys, mode: "markers", name: "samples", marker: { color: theme.accent1, size: 6 } },
                      { x: [x0], y: [f ? f(x0) : NaN], mode: "markers", name: "f(x0)", marker: { color: theme.accent2, size: 10 } },
                    ]}
                    layout={{
                      title: "Monte Carlo propagation (samples of x0 ± Δx)",
                      paper_bgcolor: theme.panelBg,
                      plot_bgcolor: theme.panelBg,
                      font: { color: theme.text },
                      height: 300,
                      margin: { t: 40, b: 50 },
                      xaxis: { title: "x" },
                      yaxis: { title: "f(x)" },
                    }}
                    useResizeHandler={true}
                    style={{ width: "100%" }}
                  />
                </div>
              </div>

              <div className="text-zinc-400 text-sm">Monte Carlo sampling captures nonlinearity and distributional features; use it when linearization assumptions are doubtful.</div>
            </div>
          )}
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ---------------- Total Error Panel ----------------
function TotalErrorPanel() {
  const [func, setFunc] = useState("exp");
  const [x0, setX0] = useState(1.0);
  const [method, setMethod] = useState("both"); // 'forward', 'central', 'both'
  const eps = Number.EPSILON || 2.220446049250313e-16;

  const errorData = useMemo(() => {
    const hs = [];
    const fwdErr = [];
    const cenErr = [];
    const f = (x) => {
      switch (func) {
        case "sin":
          return Math.sin(x);
        case "cos":
          return Math.cos(x);
        case "exp":
          return Math.exp(x);
        default:
          return Math.exp(x);
      }
    };
    const trueDer = (() => {
      switch (func) {
        case "sin":
          return Math.cos(x0);
        case "cos":
          return -Math.sin(x0);
        case "exp":
          return Math.exp(x0);
        default:
          return Math.exp(x0);
      }
    })();

    for (let k = -18; k <= 0; k += 1) {
      const h = Math.pow(10, k / 3);
      hs.push(h);
      const fd = finiteDiffForward(f, x0, h);
      const cd = finiteDiffCentral(f, x0, h);
      fwdErr.push(Math.abs(fd - trueDer));
      cenErr.push(Math.abs(cd - trueDer));
    }

    const minFwd = fwdErr.reduce((m, v, i) => (v < m.val ? { val: v, h: hs[i] } : m), { val: Infinity, h: 0 });
    const minCen = cenErr.reduce((m, v, i) => (v < m.val ? { val: v, h: hs[i] } : m), { val: Infinity, h: 0 });

    return { hs, fwdErr, cenErr, minFwd, minCen, trueDer };
  }, [func, x0]);

  return (
    <motion.div variants={cardMotion} initial="initial" animate="enter" whileHover="hover">
      <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700 overflow-hidden`}>
        <CardHeader>
          <CardTitle className="flex items-center gap-2 text-amber-300">
            <LineChart className="w-5 h-5" /> 4.3 Total Numerical Error — forward vs central
          </CardTitle>
        </CardHeader>
        <CardContent>
          <div className="grid grid-cols-1 md:grid-cols-4 gap-3 mb-3">
            <div>
              <label className="text-zinc-200 text-xs">Function</label>
              <select value={func} onChange={(e) => setFunc(e.target.value)} className="w-full bg-zinc-700 text-white p-2 rounded">
                <option value="exp">exp(x)</option>
                <option value="sin">sin(x)</option>
                <option value="cos">cos(x)</option>
              </select>
            </div>

            <div>
              <label className="text-zinc-200 text-xs">Point x₀</label>
              <Input value={String(x0)} onChange={(e) => setX0(Number(e.target.value))} className="bg-zinc-700 text-white" />
            </div>

            <div>
              <label className="text-zinc-200 text-xs">Method</label>
              <select value={method} onChange={(e) => setMethod(e.target.value)} className="w-full bg-zinc-700 text-white p-2 rounded">
                <option value="both">both</option>
                <option value="forward">forward</option>
                <option value="central">central</option>
              </select>
            </div>

            <div>
              <label className="text-zinc-200 text-xs">Machine ε</label>
              <div className="text-zinc-100 font-mono">{eps.toExponential(2)}</div>
            </div>
          </div>

          <div className="mb-3">
            <div className="w-full overflow-auto rounded">
              <Plot
                data={[
                  method !== "central" ? { x: errorData.hs, y: errorData.fwdErr, mode: "lines+markers", name: "Forward", marker: { size: 6 } } : null,
                  method !== "forward" ? { x: errorData.hs, y: errorData.cenErr, mode: "lines+markers", name: "Central", marker: { size: 6 } } : null,
                ].filter(Boolean)}
                layout={{
                  title: "Absolute error vs h (log-log)",
                  paper_bgcolor: theme.panelBg,
                  plot_bgcolor: theme.panelBg,
                  font: { color: theme.text },
                  xaxis: { type: "log", title: "h" },
                  yaxis: { type: "log", title: "absolute error" },
                  height: 380,
                }}
                useResizeHandler={true}
                style={{ width: "100%" }}
              />
            </div>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-3 gap-3">
            <div className="bg-zinc-800 p-3 rounded">
              <div className="text-zinc-400 text-sm">True derivative</div>
              <div className="text-zinc-100 font-mono">{errorData.trueDer.toExponential(6)}</div>
            </div>
            <div className="bg-zinc-800 p-3 rounded">
              <div className="text-zinc-400 text-sm">Central best</div>
              <div className="text-zinc-100 font-mono">h ≈ {errorData.minCen.h.toExponential(3)}, err ≈ {errorData.minCen.val.toExponential(3)}</div>
            </div>
            <div className="bg-zinc-800 p-3 rounded">
              <div className="text-zinc-400 text-sm">Forward best</div>
              <div className="text-zinc-100 font-mono">h ≈ {errorData.minFwd.h.toExponential(3)}, err ≈ {errorData.minFwd.val.toExponential(3)}</div>
            </div>
          </div>

          <div className="mt-3 text-zinc-400 text-sm">
            As h decreases truncation error drops, but round-off increases due to cancellation; central difference attains smaller minimal error and larger optimal h (∝ ε^{1/3}) compared to forward (∝ ε^{1/2}).
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ---------------- Blunders Panel ----------------
function BlundersPanel() {
  const [example, setExample] = useState("typo");
  const [measurement, setMeasurement] = useState(12.34);
  const [uncert, setUncert] = useState(0.05);

  const result = useMemo(() => {
    const a = 2.5,
      b = -1.2;
    const y = a * measurement + b;
    const dy = Math.abs(a) * uncert;
    return { y, dy, rel: Math.abs(dy / (y || 1)) };
  }, [measurement, uncert]);

  return (
    <motion.div variants={cardMotion} initial="initial" animate="enter" whileHover="hover">
      <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700 overflow-hidden`}>
        <CardHeader>
          <CardTitle className="flex items-center gap-2 text-indigo-300">
            <AlertTriangle className="w-5 h-5" /> 4.4 Blunders, Formulation Errors & Data Uncertainty
          </CardTitle>
        </CardHeader>
        <CardContent>
          <div className="grid grid-cols-1 sm:grid-cols-3 gap-3 mb-3">
            <div>
              <label className="text-zinc-200 text-xs">Example</label>
              <select value={example} onChange={(e) => setExample(e.target.value)} className="w-full bg-zinc-700 text-white p-2 rounded">
                <option value="typo">Typo / sign flip</option>
                <option value="model">Wrong model</option>
                <option value="data">Data uncertainty</option>
              </select>
            </div>

            <div>
              <label className="text-zinc-200 text-xs">Measurement</label>
              <Input value={String(measurement)} onChange={(e) => setMeasurement(Number(e.target.value))} className="bg-zinc-700 text-white" />
            </div>

            <div>
              <label className="text-zinc-200 text-xs">Uncertainty ±</label>
              <div className="flex items-center gap-2">
                <input type="range" min={0} max={1} step={0.001} value={uncert} onChange={(e) => setUncert(Number(e.target.value))} className="w-full" />
                <div className="w-20 text-right text-zinc-100 font-mono">{uncert.toFixed(3)}</div>
              </div>
            </div>
          </div>

          <div className="grid grid-cols-1 sm:grid-cols-3 gap-3 mb-3">
            <SummaryCard title="Result y" value={result.y.toFixed(6)} />
            <SummaryCard title="Propagated Δy" value={result.dy.toExponential(3)} />
            <SummaryCard title="Relative Δy" value={result.rel.toExponential(3)} />
          </div>

          <div className="text-zinc-400 text-sm">
            Blunders and formulation errors frequently dominate careful numerical considerations. Practices to reduce risk:
            <ul className="list-disc list-inside mt-2 text-zinc-300">
              <li>Unit tests and regression suites</li>
              <li>Dimensional analysis and sanity checks</li>
              <li>Independent validation datasets and uncertainty reporting</li>
            </ul>
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ---------------- Documentation Panel ----------------
function ChapterDocs() {
  const [open, setOpen] = useState("4.1 The Taylor Series");

  return (
    <div className="w-full mt-6 bg-zinc-900 border border-zinc-700 rounded-xl shadow-lg overflow-hidden">
      <div className="p-3 border-b border-zinc-700 flex items-center justify-between">
        <div className="flex items-center gap-3">
          <BookOpen className="w-5 h-5 text-cyan-400" />
          <div>
            <div className="text-zinc-100 font-semibold">Chapter 4 — Documentation</div>
            <div className="text-zinc-400 text-xs">Theory, practical notes, and examples for each subsection</div>
          </div>
        </div>
        <div className="text-zinc-400 text-sm">Interactive examples above — experiment and read the notes</div>
      </div>

      <div className="grid grid-cols-1 md:grid-cols-4">
        <div className="col-span-1 border-r border-zinc-800">
          {Object.keys(docsText).map((k) => (
            <button key={k} onClick={() => setOpen(k)} className={`w-full p-3 text-left hover:bg-zinc-800/40 flex items-center justify-between ${open === k ? "bg-zinc-800/20" : ""}`}>
              <div className="flex items-center gap-2">
                <List className="w-4 h-4 text-zinc-300" />
                <div className="text-zinc-100">{k}</div>
              </div>
              {open === k ? <ChevronDown className="w-4 h-4 text-zinc-400" /> : <ChevronRight className="w-4 h-4 text-zinc-400" />}
            </button>
          ))}
        </div>

        <div className="col-span-3 p-4 overflow-auto" style={{ maxHeight: 420 }}>
          <motion.h3 initial={{ opacity: 0, y: -6 }} animate={{ opacity: 1, y: 0 }} className="text-xl font-semibold text-zinc-100 mb-3">
            {open}
          </motion.h3>
          <div className="text-zinc-300 space-y-3">
            {docsText[open].map((para, i) => (
              <div key={i} className="whitespace-pre-line text-sm">{para}</div>
            ))}
            <div className="mt-3 text-zinc-400 text-sm">
              For more advanced reading, consider textbooks and resources by Dahlquist & Björck, Atkinson, and other numerical analysis references on stability, conditioning, and error analysis.
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}

// ---------------- MAIN PAGE ----------------
export default function Chapter4() {
  return (
    <div className={`min-h-screen px-4 md:px-6 py-6 bg-zinc-950 text-zinc-100`}>
      <motion.header initial={{ opacity: 0, y: -8 }} animate={{ opacity: 1, y: 0 }} className="max-w-7xl mx-auto mb-6">
        <div className="flex flex-col md:flex-row items-start justify-between gap-4">
          <div className="max-w-3xl">
            <h1 className="text-2xl text-emerald-400 md:text-3xl font-bold">Truncation Errors & the Taylor Series</h1>
            <p className="text-zinc-300 mt-2">
              Truncation error, Taylor expansions, error propagation, and how to balance truncation and round-off in practice. Interactive visualizations are included — controls collapse and stack responsively for small screens.
            </p>
          </div>

         
        </div>
      </motion.header>

      <main className="max-w-7xl mx-auto grid grid-cols-1 gap-6">
        {/* Top overview cards */}
        <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
          <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700`}>
            <CardHeader>
              <CardTitle className="text-cyan-300 flex items-center gap-2"><BookOpen className="w-5 h-5" /> Overview</CardTitle>
            </CardHeader>
            <CardContent>
              <div className="text-zinc-300 text-sm">
                This chapter covers Taylor series, truncation remainder terms, error propagation, the interplay of truncation and round-off, and non-numerical modelling errors. Use the interactive panels to experiment and learn.
              </div>
            </CardContent>
          </Card>

          <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700`}>
            <CardHeader>
              <CardTitle className="text-emerald-300 flex items-center gap-2"><Layers className="w-5 h-5" /> Visual intuition</CardTitle>
            </CardHeader>
            <CardContent>
              <div className="text-zinc-300 text-sm">Graphs and numeric experiments reveal how errors behave with respect to parameters (order N, step h, input uncertainty). Try changing values on mobile — controls stack gracefully.</div>
            </CardContent>
          </Card>

          <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700`}>
            <CardHeader>
              <CardTitle className="text-amber-300 flex items-center gap-2"><AlertTriangle className="w-5 h-5" /> Practical notes</CardTitle>
            </CardHeader>
            <CardContent>
              <ul className="text-zinc-300 text-sm list-disc list-inside">
                <li>Prefer algebraic rearrangements when possible to reduce cancellation.</li>
                <li>Balance step sizes to minimize total error.</li>
                <li>Validate with exact values where possible and report uncertainty.</li>
              </ul>
            </CardContent>
          </Card>
        </div>

        {/* Main content: two-column responsive layout */}
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          <div className="space-y-6">
            <TaylorSeriesPanel />
            <ErrorPropagationPanel />
          </div>

          <div className="space-y-6">
            <TotalErrorPanel />
            <BlundersPanel />
          </div>
        </div>

        <ChapterDocs />

        <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700`}>
          <CardHeader>
            <CardTitle className="text-emerald-300 flex items-center gap-2"><List className="w-5 h-5" /> Problems & Exercises</CardTitle>
          </CardHeader>
          <CardContent>
            <ul className="list-decimal list-inside text-zinc-300 text-sm space-y-2">
              <li>Derive the degree-4 Taylor polynomial for <code className="font-mono">e^x</code> about <code className="font-mono">a=0</code> and compute truncation error at <code className="font-mono">x=1</code>. Include a Lagrange remainder bound and numeric evaluation.</li>
              <li>For <code className="font-mono">ln(1+x)</code>, use the alternating-series bound to find the smallest <code className="font-mono">N</code> such that the remainder &lt; <code className="font-mono">10^{-6}</code> at <code className="font-mono">x=0.5</code>. Show work.</li>
              <li>Compare forward vs central difference truncation orders and simulate round-off: produce a log-log plot of error vs <code className="font-mono">h</code> and identify approximate optimal <code className="font-mono">h</code> for each method. Explain the scaling with machine ε.</li>
              <li>Give three examples of modelling errors and mitigation strategies (unit testing, validation, uncertainty quantification).</li>
            </ul>
          </CardContent>
        </Card>

        <BottomBar/>
      </main>
    </div>
  );
}
