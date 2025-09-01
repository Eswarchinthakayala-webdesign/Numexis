// src/pages/Chapter24.jsx
// Chapter 24 — Case Studies: Numerical Integration and Differentiation
// Expanded, interactive, responsive, KaTeX-ready
// Mirrors the design patterns and structure from Chapter7.jsx and Chapter23.jsx

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
  Calculator,
  Flame,
  Anchor,
  Zap,
  Activity,
  BookOpen,
  List,
  ChevronRight,
  ChevronDown,
  RefreshCcw,
  FileText,
} from "lucide-react";

/* ==========================================================================
   Theme + small helpers (consistent with Chapter7/Chapter23)
   ========================================================================== */

const theme = {
  bg: "bg-zinc-950",
  panel: "bg-zinc-900/60 border border-zinc-700 overflow-hidden",
  text: "text-zinc-200",
  subtext: "text-zinc-400",
  accent: "#22d3ee",
  accent2: "#34d399",
  warn: "#f59e0b",
};

const fadeUp = {
  initial: { opacity: 0, y: 12 },
  animate: { opacity: 1, y: 0 },
  transition: { duration: 0.28 },
};

/* ==========================================================================
   Numerical utilities used across sections
   ========================================================================== */

function linspace(a, b, n) {
  if (n <= 1) return [a];
  const dx = (b - a) / (n - 1);
  const arr = new Array(n);
  for (let i = 0; i < n; i++) arr[i] = a + i * dx;
  return arr;
}

function trapezoidIntegrate(x, y) {
  const n = x.length;
  let s = 0;
  for (let i = 0; i < n - 1; i++) {
    const dx = x[i + 1] - x[i];
    s += 0.5 * (y[i] + y[i + 1]) * dx;
  }
  return s;
}

function simpsonIntegrate(x, y) {
  const n = x.length;
  if (n < 3) return trapezoidIntegrate(x, y);
  let s = 0;
  const h = (x[n - 1] - x[0]) / (n - 1);
  // require odd number of intervals (even number of segments -> n odd). We'll handle general n by trimming last point if needed.
  let m = n;
  if ((m - 1) % 2 === 1) {
    // odd number of segments -> trim last point
    m = n - 1;
  }
  s = y[0] + y[m - 1];
  for (let i = 1; i < m - 1; i++) {
    s += (i % 2 === 0 ? 2 : 4) * y[i];
  }
  s = (h / 3) * s;
  // if trimmed, add trapezoid over last interval
  if (m !== n) {
    const last = 0.5 * (y[m - 1] + y[m]) * (x[m] - x[m - 1]);
    return s + last;
  }
  return s;
}

function cumulativeTrapezoid(x, y) {
  const n = x.length;
  const out = new Array(n).fill(0);
  for (let i = 1; i < n; i++) {
    out[i] = out[i - 1] + 0.5 * (y[i - 1] + y[i]) * (x[i] - x[i - 1]);
  }
  return out;
}

function movingAverage(arr, window = 5) {
  const n = arr.length;
  const half = Math.floor(window / 2);
  const out = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    let sum = 0;
    let cnt = 0;
    for (let j = i - half; j <= i + half; j++) {
      if (j >= 0 && j < n) {
        sum += arr[j];
        cnt++;
      }
    }
    out[i] = sum / Math.max(1, cnt);
  }
  return out;
}

/* Small linear solver for SavGol and local poly fits (naive, OK for demos) */
function solveLinearSystem(A, b) {
  const n = A.length;
  const M = A.map((row) => row.slice());
  const rhs = b.slice();
  for (let k = 0; k < n; k++) {
    let piv = k;
    for (let i = k + 1; i < n; i++) if (Math.abs(M[i][k]) > Math.abs(M[piv][k])) piv = i;
    if (Math.abs(M[piv][k]) < 1e-14) return null;
    [M[k], M[piv]] = [M[piv], M[k]];
    [rhs[k], rhs[piv]] = [rhs[piv], rhs[k]];
    for (let i = k + 1; i < n; i++) {
      const f = M[i][k] / M[k][k];
      for (let j = k; j < n; j++) M[i][j] -= f * M[k][j];
      rhs[i] -= f * rhs[k];
    }
  }
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    let s = rhs[i];
    for (let j = i + 1; j < n; j++) s -= M[i][j] * x[j];
    x[i] = s / M[i][i];
  }
  return x;
}

/* Savitzky-Golay derivative (local least squares) */
function savitzkyGolayDerivative(t, y, window = 11, degree = 3) {
  const n = y.length;
  const half = Math.floor(window / 2);
  const out = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    const lo = Math.max(0, i - half);
    const hi = Math.min(n - 1, i + half);
    const m = hi - lo + 1;
    const A = new Array(m);
    const b = new Array(m);
    for (let k = 0; k < m; k++) {
      const idx = lo + k;
      const dt = t[idx] - t[i];
      A[k] = new Array(degree + 1);
      for (let p = 0; p <= degree; p++) A[k][p] = Math.pow(dt, p);
      b[k] = y[idx];
    }
    const ATA = new Array(degree + 1).fill(0).map(() => new Array(degree + 1).fill(0));
    const ATb = new Array(degree + 1).fill(0);
    for (let r = 0; r < degree + 1; r++) {
      for (let c = 0; c < degree + 1; c++) {
        let s = 0;
        for (let k = 0; k < m; k++) s += A[k][r] * A[k][c];
        ATA[r][c] = s;
      }
      let sb = 0;
      for (let k = 0; k < m; k++) sb += A[k][r] * b[k];
      ATb[r] = sb;
    }
    const coeffs = solveLinearSystem(ATA, ATb);
    out[i] = coeffs ? (coeffs[1] ?? 0) : 0;
  }
  return out;
}

/* ==========================================================================
   Reusable UI pieces (match Chapter7 patterns)
   ========================================================================== */

function Formula({ children }) {
  return (
    <div className="my-2">
      <div className="px-3 py-2 rounded-lg border text-gray-300 border-zinc-700 bg-zinc-900/80 overflow-x-auto">
        <div className="text-xs uppercase tracking-widest mb-1" style={{ color: theme.accent }}>
          Formula
        </div>
        <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>
          {children}
        </ReactMarkdown>
      </div>
    </div>
  );
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
        <div key={i} className="rounded-lg bg-zinc-900/70 border border-zinc-700 px-3 py-2">
          <div className="text-xs text-zinc-400">{it.label}</div>
          <div className="text-zinc-100 text-sm">{it.value}</div>
        </div>
      ))}
    </div>
  );
}

function CollapsibleExercise({ title, body, solution }) {
  const [open, setOpen] = useState(false);
  return (
    <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-3">
      <div className="flex items-center justify-between">
        <div className="text-zinc-100 font-medium">{title}</div>
        <Button variant="ghost" className="bg-white  hover:bg-gray-200 cursor-pointer" size="sm" onClick={() => setOpen((o) => !o)}>
          {open ? <ChevronDown className="w-4 h-4" /> : <ChevronRight className="w-4 h-4" />}
        </Button>
      </div>
      <div className="text-zinc-300 mt-2">
        <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>
          {body}
        </ReactMarkdown>
      </div>
      {open && (
        <div className="mt-3 border-t  border-zinc-800 pt-3 text-zinc-200">
          <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>
            {solution}
          </ReactMarkdown>
        </div>
      )}
    </div>
  );
}

/* ==========================================================================
   Section 24.1 — Integration to Determine Total Heat (Chemical/Bio Eng)
   Features:
   - interactive heat flux q(t) parameters
   - multi-plot: flux, cumulative heat (trapezoid & simpson), convergence error
   - derivations and formulas
   ========================================================================== */
const Note = ({ children }) => (
  <div className="p-3 bg-zinc-800/40 border-l-4 border-emerald-500 rounded-md text-zinc-300 text-sm">
    {children}
  </div>
);

function Section241() {
  // interactive parameters for q(t)
  const [A, setA] = useState(1.0); // amplitude
  const [decay, setDecay] = useState(0.12);
  const [omega, setOmega] = useState(1.0);
  const [t0, setT0] = useState(0.0);
  const [tf, setTf] = useState(10.0);
  const [N, setN] = useState(501);

  // define q(t) = A * exp(-decay * t) * sin(omega * t)
  const t = useMemo(() => linspace(Number(t0), Number(tf), Number(N)), [t0, tf, N]);
  const q = useMemo(() => t.map((tt) => A * Math.exp(-decay * tt) * Math.sin(omega * tt)), [t, A, decay, omega]);

  const Q_trap = useMemo(() => trapezoidIntegrate(t, q), [t, q]);
  const Q_simp = useMemo(() => simpsonIntegrate(t, q), [t, q]);

  // cumulative via trapezoid for plotting
  const Qcum = useMemo(() => cumulativeTrapezoid(t, q), [t, q]);

  // convergence study: vary N and compute error relative to a fine reference
  const refN = 2001;
  const tref = useMemo(() => linspace(t0, tf, refN), [t0, tf]);
  const qref = useMemo(() => tref.map((tt) => A * Math.exp(-decay * tt) * Math.sin(omega * tt)), [tref, A, decay, omega]);
  const Qref = useMemo(() => simpsonIntegrate(tref, qref), [tref, qref]);

  const Ns = [51, 101, 201, 401, 801, 1601].map((n) => Math.min(n, 2001));
  const conv = Ns.map((nn) => {
    const tt = linspace(t0, tf, nn);
    const y = tt.map((ttt) => A * Math.exp(-decay * ttt) * Math.sin(omega * ttt));
    const Qs = simpsonIntegrate(tt, y);
    return { N: nn, err: Math.abs(Qs - Qref) };
  });

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Flame className="w-5 h-5" />
            24.1 Integration to Determine Total Heat (Chemical/Bio Eng.)
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          {/* formulas / derivation */}
          <div className="grid md:grid-cols-2 gap-4">
            <div>
              <Formula>{`**Total heat**\n\nThe total heat delivered over time interval $[t_0,t_f]$ by heat flux $\\dot{q}(t)$ is\n\n$$Q=\\int_{t_0}^{t_f} \\dot{q}(t)\\, dt.$$`}</Formula>

              <Formula>{`**Numerical rules**\n\nTrapezoidal rule:\n\n$$Q\\approx\\sum_{i=0}^{n-1} \\frac{\\dot{q}_i+\\dot{q}_{i+1}}{2} (t_{i+1}-t_i).$$\n\nSimpson's rule (even segments):\n\n$$Q\\approx \\frac{h}{3}\\left(f_0+f_n + 4\\sum_{odd}f_i + 2\\sum_{even\\neq 0,n}f_i\\right).$$`}</Formula>

              <Note>
                Choose step size to balance truncation error and rounding. For decaying oscillatory heat flux, Simpson often converges faster than trapezoid for smooth data.
              </Note>
            </div>

            <div>
              <div className="grid grid-cols-2 gap-2">
                <Labeled label="Amplitude A">
                  <Input value={A} onChange={(e) => setA(Number(e.target.value))} className="bg-zinc-800 text-white" />
                </Labeled>
                <Labeled label="Decay">
                  <Input value={decay} onChange={(e) => setDecay(Number(e.target.value))} className="bg-zinc-800 text-white" />
                </Labeled>
                <Labeled label="Omega (rad/s)">
                  <Input value={omega} onChange={(e) => setOmega(Number(e.target.value))} className="bg-zinc-800 text-white" />
                </Labeled>
                <Labeled label="N samples">
                  <Input value={N} onChange={(e) => setN(Number(e.target.value))} className="bg-zinc-800 text-white" />
                </Labeled>
                <Labeled label="t0">
                  <Input value={t0} onChange={(e) => setT0(Number(e.target.value))} className="bg-zinc-800 text-white" />
                </Labeled>
                <Labeled label="tf">
                  <Input value={tf} onChange={(e) => setTf(Number(e.target.value))} className="bg-zinc-800 text-white" />
                </Labeled>
              </div>

              <div className="mt-3">
                <Summary items={[
                  { label: "Trapezoid Q", value: Q_trap.toFixed(6) },
                  { label: "Simpson Q", value: Q_simp.toFixed(6) },
                  { label: "Reference (fine Simpson)", value: Qref.toFixed(6) },
                ]} />
              </div>
            </div>
          </div>

          {/* plots: flux, cumulative, convergence */}
          <div className="grid md:grid-cols-3 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: t, y: q, type: "scatter", mode: "lines", name: "q(t)", line: { color: theme.accent } },
                  ]}
                  layout={{
                    title: "Heat Flux q(t)",
                    autosize: true,
                    responsive: true,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    margin: { t: 36, l: 48, r: 24, b: 36 },
                    font: { color: "#e5e7eb" },
                  }}
                  useResizeHandler
                  style={{ width: "100%", height: "100%" }}
                  config={{ responsive: true, displayModeBar: false }}
                />
              </div>
            </div>

            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: t, y: Qcum, type: "scatter", mode: "lines", name: "Q_cum (trap)", line: { color: theme.accent2 } },
                  ]}
                  layout={{
                    title: "Cumulative Heat (Trapezoid)",
                    autosize: true,
                    responsive: true,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    margin: { t: 36, l: 48, r: 24, b: 36 },
                    font: { color: "#e5e7eb" },
                  }}
                  useResizeHandler
                  style={{ width: "100%", height: "100%" }}
                  config={{ responsive: true, displayModeBar: false }}
                />
              </div>
            </div>

            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: Ns.map((d) => d), y: conv.map((c) => c.err), type: "scatter", mode: "lines+markers", name: "abs error", line: { color: theme.warn } },
                  ]}
                  layout={{
                    title: "Convergence: error vs N (Simpson)",
                    xaxis: { title: "N (samples)" },
                    yaxis: { title: "abs(error)", type: "log" },
                    autosize: true,
                    responsive: true,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    margin: { t: 36, l: 48, r: 24, b: 36 },
                    font: { color: "#e5e7eb" },
                  }}
                  useResizeHandler
                  style={{ width: "100%", height: "100%" }}
                  config={{ responsive: true, displayModeBar: false }}
                />
              </div>
            </div>
          </div>

          <CollapsibleExercise
            title="Exercise 24.1.1 — Adaptive integration"
            body={`Implement a simple adaptive Simpson's routine to compute $Q$. Compare performance (function evaluations) vs fixed-step Simpson for target tolerance $10^{-6}$.`}
            solution={`Adaptive Simpson reduces function evaluations for smooth regions and concentrates points where integrand varies rapidly. Implementation: recursively split intervals until estimated error < tol.`}
          />
        </CardContent>
      </Card>
    </motion.div>
  );
}

/* ==========================================================================
   Section 24.2 — Effective Force on Mast of a Racing Sailboat
   - pressure distribution p(z)
   - integrate p(z) over mast height to get total force and center of pressure
   - interactive mast height H, wind coeffs
   - plots: p(z), cumulative force, center of pressure marker, convergence
   ========================================================================== */

function Section242() {
  const [H, setH] = useState(20); // mast height meters
  const [p0, setP0] = useState(120); // near-base reference pressure (N/m^2) scale
  const [decay, setDecay] = useState(0.08); // how pressure decays with height
  const [N, setN] = useState(401);

  const z = useMemo(() => linspace(0, H, N), [H, N]);
  // build pressure profile: p(z) = p0 * exp(-decay*z) * (1 + 0.2 sin(2*pi*z/H))
  const p = useMemo(() => z.map((zz) => p0 * Math.exp(-decay * zz) * (1 + 0.2 * Math.sin((2 * Math.PI * zz) / H))), [z, p0, decay, H]);

  const F_total = useMemo(() => trapezoidIntegrate(z, p), [z, p]);
  // center of pressure z_cp = (1/F) * integral z * p(z) dz
  const zp = useMemo(() => z.map((zz, i) => zz * p[i]), [z, p]);
  const M = useMemo(() => trapezoidIntegrate(z, zp), [z, zp]);
  const z_cp = useMemo(() => (F_total !== 0 ? M / F_total : 0), [M, F_total]);

  // cumulative force
  const Fcum = useMemo(() => cumulativeTrapezoid(z, p), [z, p]);

  // convergence for different N relative to fine reference
  const refN = 2001;
  const zref = useMemo(() => linspace(0, H, refN), [H]);
  const pref = useMemo(() => zref.map((zz) => p0 * Math.exp(-decay * zz) * (1 + 0.2 * Math.sin((2 * Math.PI * zz) / H))), [zref, p0, decay, H]);
  const Fref = useMemo(() => simpsonIntegrate(zref, pref), [zref, pref]);

  const Ns = [51, 101, 201, 401, 801, 1601].map((n) => Math.min(n, refN));
  const conv = Ns.map((nn) => {
    const zz = linspace(0, H, nn);
    const yy = zz.map((zzz) => p0 * Math.exp(-decay * zzz) * (1 + 0.2 * Math.sin((2 * Math.PI * zzz) / H)));
    const Ft = simpsonIntegrate(zz, yy);
    return { N: nn, err: Math.abs(Ft - Fref) };
  });

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Anchor className="w-5 h-5" />
            24.2 Effective Force on the Mast of a Racing Sailboat
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid md:grid-cols-2 gap-4">
            <div>
              <Formula>{`Total force:\n\n$$F = \\int_0^{H} p(z)\\, dz.$$`}</Formula>
              <Formula>{`Center of pressure:\n\n$$z_{cp} = \\frac{1}{F} \\int_0^{H} z \\, p(z)\\, dz.$$`}</Formula>
              <Note>
                For structural design, both total force and center of pressure are required to compute bending moments at the base of the mast.
              </Note>
            </div>

            <div>
              <div className="grid grid-cols-2 gap-2">
                <Labeled label="Mast height H (m)">
                  <Input value={H} onChange={(e) => setH(Number(e.target.value))} className="bg-zinc-800 text-white" />
                </Labeled>

                <Labeled label="p0 (pressure scale)">
                  <Input value={p0} onChange={(e) => setP0(Number(e.target.value))} className="bg-zinc-800 text-white" />
                </Labeled>

                <Labeled label="Decay rate">
                  <Input value={decay} onChange={(e) => setDecay(Number(e.target.value))} className="bg-zinc-800 text-white" />
                </Labeled>

                <Labeled label="N samples">
                  <Input value={N} onChange={(e) => setN(Number(e.target.value))} className="bg-zinc-800 text-white" />
                </Labeled>
              </div>

              <div className="mt-3">
                <Summary items={[
                  { label: "Total Force F (N)", value: F_total.toFixed(2) },
                  { label: "Center of Pressure z_cp (m)", value: z_cp.toFixed(3) },
                  { label: "Reference (fine simpson)", value: Fref.toFixed(2) },
                ]} />
              </div>
            </div>
          </div>

          <div className="grid md:grid-cols-3 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: z, y: p, mode: "lines", type: "scatter", name: "p(z)", line: { color: theme.accent } },
                    { x: [z_cp], y: [p[Math.max(0, Math.floor((z_cp / H) * (z.length - 1)))]], mode: "markers", type: "scatter", name: "z_cp", marker: { color: theme.accent2, size: 10 } },
                  ]}
                  layout={{
                    title: "Pressure Distribution p(z)",
                    xaxis: { title: "z (m)" },
                    autosize: true,
                    responsive: true,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    margin: { t: 36, l: 48, r: 20, b: 36 },
                    font: { color: "#e5e7eb" },
                  }}
                  useResizeHandler
                  style={{ width: "100%", height: "100%" }}
                  config={{ responsive: true, displayModeBar: false }}
                />
              </div>
            </div>

            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: z, y: Fcum, mode: "lines", type: "scatter", name: "F_cum(z)", line: { color: theme.accent2 } },
                    { x: [z_cp], y: [M / F_total], mode: "markers", type: "scatter", name: "z_cp", marker: { color: theme.accent, size: 8 } },
                  ]}
                  layout={{
                    title: "Cumulative Force vs Height",
                    autosize: true,
                    responsive: true,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    margin: { t: 36, l: 48, r: 20, b: 36 },
                    font: { color: "#e5e7eb" },
                  }}
                  useResizeHandler
                  style={{ width: "100%", height: "100%" }}
                  config={{ responsive: true, displayModeBar: false }}
                />
              </div>
            </div>

            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: Ns.map((d) => d), y: conv.map((c) => c.err), type: "scatter", mode: "lines+markers", name: "abs error", line: { color: theme.warn } },
                  ]}
                  layout={{
                    title: "Convergence: abs error vs N",
                    xaxis: { title: "N (samples)" },
                    yaxis: { title: "abs(error)", type: "log" },
                    autosize: true,
                    responsive: true,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    margin: { t: 36, l: 48, r: 20, b: 36 },
                    font: { color: "#e5e7eb" },
                  }}
                  useResizeHandler
                  style={{ width: "100%", height: "100%" }}
                  config={{ responsive: true, displayModeBar: false }}
                />
              </div>
            </div>
          </div>

          <CollapsibleExercise
            title="Exercise 24.2.1 — Bending moment at base"
            body={`Given $p(z)$ and mast radius $r$, compute bending moment at base: \\(M = \\int_0^H z\, p(z)\\, dz \\times r\\). Numerically compute M and discuss how center of pressure affects required base design.`}
            solution={`Bending moment depends on both total force and center of pressure. Increasing z_cp moves resultant farther from the base increasing moment linearly. Compute integral numerically and multiply by lever arm r.`}
          />
        </CardContent>
      </Card>
    </motion.div>
  );
}

/* ==========================================================================
   Section 24.3 — Root-Mean-Square Current by Numerical Integration
   - compute current waveform i(t)
   - compute Irms by integral of i^2 over period
   - interactive harmonic content and windowing
   - plots: waveform, squared waveform, integration convergence
   ========================================================================== */

function Section243() {
  const [f0, setF0] = useState(1.0); // base frequency (Hz)
  const [A1, setA1] = useState(1.0);
  const [A3, setA3] = useState(0.3);
  const [phi3, setPhi3] = useState(0.3);
  const [cycles, setCycles] = useState(3);
  const [N, setN] = useState(2001);

  const T = 1 / f0;
  const tStart = 0;
  const tEnd = cycles * T;
  const t = useMemo(() => linspace(tStart, tEnd, N), [tStart, tEnd, N]);

  // waveform: i(t) = A1*sin(2πf0 t) + A3*sin(2π*3f0 t + phi3)
  const i = useMemo(() => t.map((tt) => A1 * Math.sin(2 * Math.PI * f0 * tt) + A3 * Math.sin(2 * Math.PI * 3 * f0 * tt + phi3)), [t, A1, A3, phi3, f0]);
  const i2 = useMemo(() => i.map((v) => v * v), [i]);

  // compute Irms over one period (take last full period)
  const periodIndices = useMemo(() => {
    const start = Math.max(0, Math.floor((cycles - 1) * N / cycles));
    const end = N - 1;
    return { start, end };
  }, [N, cycles]);

  const tPeriod = useMemo(() => t.slice(periodIndices.start, periodIndices.end + 1), [t, periodIndices]);
  const iPeriod = useMemo(() => i.slice(periodIndices.start, periodIndices.end + 1), [i, periodIndices]);

  const I2integral = useMemo(() => simpsonIntegrate(tPeriod, iPeriod.map((v) => v * v)), [tPeriod, iPeriod]);
  const Irms = useMemo(() => Math.sqrt(I2integral / (tPeriod[tPeriod.length - 1] - tPeriod[0])), [I2integral, tPeriod]);

  // convergence: vary N per period
  const Ns = [101, 201, 401, 801, 1601, 3201].map((nn) => Math.min(nn, 4001));
  const refN = 8001;
  const tref = linspace(tStart, tEnd, refN);
  const iref = tref.map((tt) => A1 * Math.sin(2 * Math.PI * f0 * tt) + A3 * Math.sin(2 * Math.PI * 3 * f0 * tt + phi3));
  const Iref = Math.sqrt(simpsonIntegrate(tref.slice(refN - Math.floor(refN / cycles), refN).map((_, idx) => tref[refN - Math.floor(refN / cycles) + idx]), iref.slice(iref.length - Math.floor(refN / cycles)).map((v) => v * v)) / T);

  const conv = Ns.map((nn) => {
    const tt = linspace(tStart, tEnd, nn);
    const ii = tt.map((ttt) => A1 * Math.sin(2 * Math.PI * f0 * ttt) + A3 * Math.sin(2 * Math.PI * 3 * f0 * ttt + phi3));
    const iper = ii.slice(Math.floor((cycles - 1) * (nn / cycles)));
    const tp = tt.slice(Math.floor((cycles - 1) * (nn / cycles)));
    const I2 = simpsonIntegrate(tp, iper.map((v) => v * v));
    const Ir = Math.sqrt(I2 / (tp[tp.length - 1] - tp[0]));
    return { N: nn, err: Math.abs(Ir - Iref) };
  });

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Zap className="w-5 h-5" />
            24.3 Root-Mean-Square Current by Numerical Integration
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid md:grid-cols-2 gap-4">
            <div>
              <Formula>{`RMS current (one period):\n\n$$I_{rms}=\\sqrt{\\frac{1}{T}\\int_0^{T} i^2(t)\\, dt}.$$`}</Formula>
              <Note>
                For periodic signals containing harmonics, ensure integration spans an integer number of periods. Windowing non-integer cycles creates bias.
              </Note>
              <div className="mt-3">
                <Summary items={[
                  { label: "Computed I_rms", value: Irms.toFixed(6) },
                  { label: "Period T", value: T.toFixed(4) },
                  { label: "Samples N", value: N },
                ]} />
              </div>
            </div>

            <div>
              <div className="grid grid-cols-2 gap-2">
                <Labeled label="Base freq f0 (Hz)">
                  <Input value={f0} onChange={(e) => setF0(Number(e.target.value))} className="bg-zinc-800 text-white" />
                </Labeled>
                <Labeled label="Cycles">
                  <Input value={cycles} onChange={(e) => setCycles(Number(e.target.value))} className="bg-zinc-800 text-white" />
                </Labeled>

                <Labeled label="A1 (fundamental)">
                  <Input value={A1} onChange={(e) => setA1(Number(e.target.value))} className="bg-zinc-800 text-white" />
                </Labeled>
                <Labeled label="A3 (3rd harmonic)">
                  <Input value={A3} onChange={(e) => setA3(Number(e.target.value))} className="bg-zinc-800 text-white" />
                </Labeled>

                <Labeled label="Phi3 (rad)">
                  <Input value={phi3} onChange={(e) => setPhi3(Number(e.target.value))} className="bg-zinc-800 text-white" />
                </Labeled>
                <Labeled label="Samples N">
                  <Input value={N} onChange={(e) => setN(Number(e.target.value))} className="bg-zinc-800 text-white" />
                </Labeled>
              </div>
            </div>
          </div>

          <div className="grid md:grid-cols-3 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: t, y: i, type: "scatter", mode: "lines", name: "i(t)", line: { color: theme.accent } },
                  ]}
                  layout={{
                    title: "Current waveform i(t)",
                    xaxis: { title: "t (s)" },
                    autosize: true,
                    responsive: true,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    margin: { t: 36, l: 48, r: 20, b: 36 },
                    font: { color: "#e5e7eb" },
                  }}
                  useResizeHandler
                  style={{ width: "100%", height: "100%" }}
                  config={{ responsive: true, displayModeBar: false }}
                />
              </div>
            </div>

            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: t, y: i2, type: "scatter", mode: "lines", name: "i^2(t)", line: { color: theme.accent2 } },
                  ]}
                  layout={{
                    title: "Squared waveform i^2(t)",
                    autosize: true,
                    responsive: true,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    margin: { t: 36, l: 48, r: 20, b: 36 },
                    font: { color: "#e5e7eb" },
                  }}
                  useResizeHandler
                  style={{ width: "100%", height: "100%" }}
                  config={{ responsive: true, displayModeBar: false }}
                />
              </div>
            </div>

            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: Ns, y: conv.map((c) => c.err), type: "scatter", mode: "lines+markers", name: "abs error", line: { color: theme.warn } },
                  ]}
                  layout={{
                    title: "Convergence: Irms error vs N",
                    xaxis: { title: "N (samples)" },
                    yaxis: { title: "abs error", type: "log" },
                    autosize: true,
                    responsive: true,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    margin: { t: 36, l: 48, r: 20, b: 36 },
                    font: { color: "#e5e7eb" },
                  }}
                  useResizeHandler
                  style={{ width: "100%", height: "100%" }}
                  config={{ responsive: true, displayModeBar: false }}
                />
              </div>
            </div>
          </div>

          <CollapsibleExercise
            title="Exercise 24.3.1 — Windowing and aliasing"
            body={`Investigate how sampling at low N causes aliasing in computed Irms for a signal with high harmonic content (increase A3). Discuss the effect of using anti-aliasing filtering or higher sampling rate before integration.`}
            solution={`Aliasing biases the integral and computed Irms. Filtering or oversampling resolves aliasing; use Nyquist criterion (fs > 2*f_max) to choose N per period.`}
          />
        </CardContent>
      </Card>
    </motion.div>
  );
}

/* ==========================================================================
   Section 24.4 — Numerical Integration to Compute Work (Mechanical/Aerospace)
   - F(x) functions and parameters
   - integrate F(x) dx to get work
   - show force-displacement curve, cumulative work, Simpson vs trapezoid, error, energy interpretation
   ========================================================================== */

function Section244() {
  const [F0, setF0] = useState(50); // baseline force
  const [amp, setAmp] = useState(20); // amplitude of variation
  const [L, setL] = useState(6.28); // displacement range ( ~ 2π )
  const [N, setN] = useState(401);

  const x = useMemo(() => linspace(0, L, N), [L, N]);
  const F = useMemo(() => x.map((xx) => F0 + amp * Math.sin(2 * Math.PI * xx / L)), [x, F0, amp, L]);

  const W_trap = useMemo(() => trapezoidIntegrate(x, F), [x, F]);
  const W_simp = useMemo(() => simpsonIntegrate(x, F), [x, F]);

  const Wcum = useMemo(() => cumulativeTrapezoid(x, F), [x, F]);

  // analytic integral if F(x) = F0 + amp * sin(2π x / L): integral = F0*L - (amp * (L/(2π)) * (cos(2π x / L) - 1)) from 0..L -> F0*L
  const W_exact = useMemo(() => F0 * L, [F0, L]);

  // convergence vs N
  const Ns = [41, 81, 161, 321, 641, 1281].map((n) => Math.min(n, 2001));
  const conv = Ns.map((nn) => {
    const xx = linspace(0, L, nn);
    const FF = xx.map((xxx) => F0 + amp * Math.sin((2 * Math.PI * xxx) / L));
    const W = simpsonIntegrate(xx, FF);
    return { N: nn, err: Math.abs(W - W_exact) };
  });

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Activity className="w-5 h-5" />
            24.4 Numerical Integration to Compute Work (Mechanical/Aerospace)
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid md:grid-cols-2 gap-4">
            <div>
              <Formula>{`Work done over [0,L]:\n\n$$W = \\int_{0}^{L} F(x)\\, dx.$$`}</Formula>

              <Formula>{`If \\(F(x)=F_0 + A\\sin\\left(\\frac{2\\pi x}{L}\\right)\\) then analytic result:\n\n$$W = F_0 L$$ (since integral of the sine over full period is zero).`}</Formula>

              <Note>
                In many engineering tasks the force varies cyclically; numerical integration confirms whether net work over a cycle is zero or yields expected baseline contribution.
              </Note>
            </div>

            <div>
              <div className="grid grid-cols-2 gap-2">
                <Labeled label="F0 (N)">
                  <Input value={F0} onChange={(e) => setF0(Number(e.target.value))} className="bg-zinc-800 text-white" />
                </Labeled>
                <Labeled label="Amplitude A (N)">
                  <Input value={amp} onChange={(e) => setAmp(Number(e.target.value))} className="bg-zinc-800 text-white" />
                </Labeled>

                <Labeled label="Displacement L">
                  <Input value={L} onChange={(e) => setL(Number(e.target.value))} className="bg-zinc-800 text-white" />
                </Labeled>

                <Labeled label="Samples N">
                  <Input value={N} onChange={(e) => setN(Number(e.target.value))} className="bg-zinc-800 text-white" />
                </Labeled>
              </div>

              <div className="mt-3">
                <Summary items={[
                  { label: "Trap W", value: W_trap.toFixed(4) },
                  { label: "Simpson W", value: W_simp.toFixed(6) },
                  { label: "Analytic W", value: W_exact.toFixed(6) },
                ]} />
              </div>
            </div>
          </div>

          <div className="grid md:grid-cols-3 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: x, y: F, type: "scatter", mode: "lines", name: "F(x)", line: { color: theme.accent } },
                  ]}
                  layout={{
                    title: "Force vs Displacement",
                    xaxis: { title: "x" },
                    autosize: true,
                    responsive: true,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    margin: { t: 36, l: 48, r: 20, b: 36 },
                    font: { color: "#e5e7eb" },
                  }}
                  useResizeHandler
                  style={{ width: "100%", height: "100%" }}
                  config={{ responsive: true, displayModeBar: false }}
                />
              </div>
            </div>

            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: x, y: Wcum, type: "scatter", mode: "lines", name: "W_cum", line: { color: theme.accent2 } },
                  ]}
                  layout={{
                    title: "Cumulative Work (Trapezoid)",
                    autosize: true,
                    responsive: true,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    margin: { t: 36, l: 48, r: 20, b: 36 },
                    font: { color: "#e5e7eb" },
                  }}
                  useResizeHandler
                  style={{ width: "100%", height: "100%" }}
                  config={{ responsive: true, displayModeBar: false }}
                />
              </div>
            </div>

            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2 overflow-hidden">
              <div className="w-full h-[360px]">
                <Plot
                  data={[
                    { x: Ns, y: conv.map((c) => c.err), type: "scatter", mode: "lines+markers", name: "abs error", line: { color: theme.warn } },
                  ]}
                  layout={{
                    title: "Convergence: abs error vs N",
                    xaxis: { title: "N (samples)" },
                    yaxis: { title: "abs(error)", type: "log" },
                    autosize: true,
                    responsive: true,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    margin: { t: 36, l: 48, r: 20, b: 36 },
                    font: { color: "#e5e7eb" },
                  }}
                  useResizeHandler
                  style={{ width: "100%", height: "100%" }}
                  config={{ responsive: true, displayModeBar: false }}
                />
              </div>
            </div>
          </div>

          <CollapsibleExercise
            title="Exercise 24.4.1 — Non-conservative forces"
            body={`If F(x) includes non-conservative elements (like drag dependent on velocity), numerical integration over a path may require coupled ODE solving for x(t) and F(x, v). Sketch an approach where you estimate work numerically when force depends on both position and velocity.`}
            solution={`Use time-stepping ODE integrator (e.g., RK4) to obtain x(t) and v(t), then form instantaneous power P = F(x,v)*v and integrate over time for work: W = \\int P dt.`}
          />
        </CardContent>
      </Card>
    </motion.div>
  );
}

/* ==========================================================================
   Documentation panel (contents / notes) - mirrors Chapter7 pattern
   ========================================================================== */

function ChapterDocs() {
  const [open, setOpen] = useState("24.1 Integration to Determine the Total Quantity of Heat (Chemical/Bio Engineering)");

  const docsMap = {
    "24.1 Integration to Determine the Total Quantity of Heat (Chemical/Bio Engineering)": [
      `Total heat over time interval:  $$Q=\\int_{t_0}^{t_f} \\dot{q}(t)\\, dt$$`,
      `Use adaptive integration for sharply varying flux; Simpson/trapezoid converge with rates O(h^4)/O(h^2) respectively.`,
    ],
    "24.2 Effective Force on the Mast of a Racing Sailboat (Civil/Environmental Engineering)": [
      `Total force: $$F=\\int_0^H p(z)\\, dz$$`,
      `Center of pressure: $$z_{cp}=\\frac{1}{F}\\int_0^H z p(z)\\, dz$$`,
    ],
    "24.3 Root-Mean-Square Current by Numerical Integration (Electrical Engineering)": [
      `RMS current: $$I_{rms}=\\sqrt{\\frac{1}{T}\\int_0^T i^2(t)\\, dt}$$`,
      `Make sure the integration window covers an integer number of periods to avoid bias.`,
    ],
    "24.4 Numerical Integration to Compute Work (Mechanical/Aerospace Engineering)": [
      `Work: $$W=\\int_{x_0}^{x_f} F(x)\\, dx$$`,
      `For analytic forms, compare numerical integration to closed-form results for verification.`,
    ],
  };

  return (
    <div className="w-full mt-6 bg-zinc-900 border border-zinc-700 rounded-xl shadow-lg overflow-hidden">
      <div className="p-3 border-b border-zinc-700 flex items-center justify-between">
        <div className="flex items-center gap-3">
          <BookOpen className="w-5 h-5 text-cyan-400" />
          <div>
            <div className="text-zinc-100 font-semibold">Chapter 24 — Documentation & Notes</div>
            <div className="text-zinc-400 text-xs">Step-by-step derivations, formulas, and references</div>
          </div>
        </div>
        <div className="text-zinc-400 text-xs">KaTeX-rendered</div>
      </div>

      <div className="grid grid-cols-1 md:grid-cols-4">
        <div className="col-span-1 border-r border-zinc-800">
          {Object.keys(docsMap).map((k) => (
            <button
              key={k}
              onClick={() => setOpen(k)}
              className={`w-full p-3 text-left hover:bg-zinc-800/40 flex items-center justify-between ${open === k ? "bg-zinc-800/20" : ""}`}
            >
              <div className="flex items-center gap-2">
                <List className="w-4 h-4 text-zinc-300" />
                <div className="text-zinc-100 text-sm">{k}</div>
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
            {docsMap[open].map((para, i) => (
              <div key={i}>
                <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>
                  {para}
                </ReactMarkdown>
              </div>
            ))}
            <div className="mt-3 text-zinc-400 text-xs">
              References: standard texts on numerical methods (Atkinson, Burden & Faires, Press et al. Numerical Recipes).
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}

/* ==========================================================================
   Problems panel (end-of-chapter)
   ========================================================================== */

function ProblemsPanel() {
  return (
    <div className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
      <div className="flex items-center justify-between mb-3">
        <div className="text-xl text-cyan-300 font-semibold flex items-center gap-2"><FileText className="w-5 h-5" /> Problems</div>
        <Button variant="ghost" size="sm" onClick={() => window.scrollTo({ top: 0, behavior: "smooth" })}><RefreshCcw className="w-4 h-4" /> Top</Button>
      </div>

      <div className="grid md:grid-cols-2 gap-3">
        <CollapsibleExercise
          title="Problem 24.1 — Heat accumulation in a batch reactor"
          body={`Given \\(\\dot{q}(t)=A e^{-\\alpha t}\\sin(\\omega t)\\), compute total heat from t=0 to t=10 numerically. Compare trap and Simpson for N=50,100,400 and discuss which is preferable for a given cost budget (function evaluations).`}
          solution={`Simpson converges faster for smooth integrands; for given budget choose Simpson with coarser N. If integrand has sharp spikes use adaptive methods.`}
        />
        <CollapsibleExercise
          title="Problem 24.2 — Mast: distributed loading"
          body={`Given measured pressure samples at irregular heights, compute total force and center of pressure using local polynomial fits and compare with trapezoid on interpolated uniform grid. Discuss sensitivity to measurement noise.`}
          solution={`Local polynomial fits reduce sensitivity to irregular spacing; smoothing before differentiation/integration can reduce noise amplification.`}
        />
        <CollapsibleExercise
          title="Problem 24.3 — RMS for distorted waveforms"
          body={`Construct a waveform with harmonics up to 7f0, compute Irms numerically and compare with spectral computation (Parseval) if FFT available. Which is more efficient for long signals?`}
          solution={`FFT + Parseval can compute power efficiently for long periodic signals; direct time-domain integration is simpler for short windows.`}
        />
        <CollapsibleExercise
          title="Problem 24.4 — Work done by aerodynamic forces"
          body={`Given a variable drag profile F(x,v)=C_d v^2(x) along a path x(t), outline a numerical approach to compute total work combining ODE integration and power integration.`}
          solution={`Integrate equations of motion for x(t) and v(t) using RK4, compute instantaneous power P=F(x,v)*v and integrate P over time using Simpson/trapezoid.`}
        />
      </div>
    </div>
  );
}

/* ==========================================================================
   Page assembly
   ========================================================================== */

export default function Chapter24() {
  return (
    <div className={`p-6 space-y-6 ${theme.bg}`}>
      <motion.div initial={{ opacity: 0, y: -10 }} animate={{ opacity: 1, y: 0 }}>
        <h1 className="text-3xl font-bold" style={{ color: theme.accent2 }}>
          Case Studies: Numerical Integration & Differentiation
        </h1>
        <div className="text-zinc-400">Applied examples from engineering domains with interactive numerical experiments.</div>
      </motion.div>

      <Section241 />
      <Section242 />
      <Section243 />
      <Section244 />

      <ChapterDocs />

      <ProblemsPanel />
      <BottomBar/>
    </div>
  );
}
