// src/pages/Chapter5_Bracketing.jsx
import React, { useState, useMemo, useRef } from "react";
import { motion } from "framer-motion";
import Plot from "react-plotly.js";
import { Canvas, useFrame } from "@react-three/fiber";
import { OrbitControls } from "@react-three/drei";
import * as THREE from "three";
import {
  Search as SearchIcon,
  SquarePlus,
  ChevronRight,
  ChevronDown,
  BookOpen,
  Circle,
  Zap,
  List,
  Layers,
  Info,
  GitMerge,
  AlertTriangle,
  Slash,
} from "lucide-react";

import { Card, CardHeader, CardTitle, CardContent } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import BottomBar from "../components/BottomBar";

// ---------------- Theme & motion (match Chapter4 style) ----------------
const theme = {
  background: "bg-zinc-900",
  panelBg: "#0f1720",
  text: "#e6eef6",
  accent1: "#60a5fa", // blue
  accent2: "#34d399", // emerald
  accent3: "#f59e0b", // amber
};

const cardMotion = {
  initial: { opacity: 0, y: 10, scale: 0.995 },
  enter: { opacity: 1, y: 0, scale: 1, transition: { duration: 0.36 } },
  hover: { scale: 1.02, y: -6, transition: { duration: 0.16 } },
};

// ---------------- Numeric helpers for root finding ----------------
function safeEvalFn(expr) {
  try {
    // eslint-disable-next-line no-new-func
    return new Function("x", `with (Math) { return (${expr}); }`);
  } catch (e) {
    return null;
  }
}

function bisectionStep(f, a, b) {
  const c = (a + b) / 2;
  const fa = f(a);
  const fc = f(c);
  if (Number.isNaN(fa) || Number.isNaN(fc)) return null;
  if (fa * fc <= 0) return { a, b: c, c, fa, fc };
  return { a: c, b, c, fa, fc };
}

function falsePositionStep(f, a, b) {
  const fa = f(a);
  const fb = f(b);
  const c = (a * fb - b * fa) / (fb - fa);
  const fc = f(c);
  if (Number.isNaN(fa) || Number.isNaN(fb) || Number.isNaN(fc)) return null;
  if (fa * fc <= 0) return { a, b: c, c, fa, fb, fc };
  return { a: c, b, c, fa, fb, fc };
}

function findSignChanges(f, left, right, nSteps) {
  const xs = [];
  const step = (right - left) / nSteps;
  let prev = left;
  let prevVal = f(prev);
  for (let i = 1; i <= nSteps; i++) {
    const x = left + i * step;
    const val = f(x);
    if (Number.isFinite(prevVal) && Number.isFinite(val) && prevVal * val <= 0) {
      xs.push([prev, x]);
    }
    prev = x;
    prevVal = val;
  }
  return xs;
}

// ---------------- Small 3D visual root cloud ----------------
function RootCloud({ roots = [], color = theme.accent2 }) {
  const ref = useRef();
  useFrame(({ clock }) => {
    if (!ref.current) return;
    ref.current.rotation.y = clock.getElapsedTime() * 0.08;
  });

  const pts = useMemo(() => {
    const arr = [];
    roots.forEach((r, i) => {
      const angle = (i / Math.max(1, roots.length)) * Math.PI * 2;
      arr.push(r, Math.sin(angle) * 0.4, Math.cos(angle) * 0.4);
    });
    return new Float32Array(arr);
  }, [roots]);

  if (!roots.length) return null;

  return (
    <points ref={ref} position={[0, 0, 0]}>
      <bufferGeometry>
        <bufferAttribute attach="attributes-position" array={pts} count={pts.length / 3} itemSize={3} />
      </bufferGeometry>
      <pointsMaterial size={0.08} sizeAttenuation color={color} />
    </points>
  );
}

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

// ---------------- Panels ----------------
function GraphicalPanel() {
  const [expr, setExpr] = useState("Math.sin(x) - 0.5");
  const [left, setLeft] = useState(-6.28);
  const [right, setRight] = useState(6.28);
  const [samples, setSamples] = useState(400);

  const f = useMemo(() => safeEvalFn(expr), [expr]);

  const plot = useMemo(() => {
    const xs = [];
    const ys = [];
    const step = (right - left) / samples;
    for (let i = 0; i <= samples; i++) {
      const x = left + i * step;
      xs.push(x);
      try {
        ys.push(f ? f(x) : NaN);
      } catch {
        ys.push(NaN);
      }
    }
    return { xs, ys };
  }, [f, left, right, samples]);

  const signIntervals = useMemo(() => {
    if (!f) return [];
    return findSignChanges(f, left, right, 200);
  }, [f, left, right]);

  return (
    <motion.div variants={cardMotion} initial="initial" animate="enter" whileHover="hover">
      <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700 overflow-hidden`}>
        <CardHeader>
          <CardTitle className="flex items-center gap-2 text-cyan-300"><SearchIcon className="w-5 h-5" /> 5.1 Graphical Methods — function plot & sign changes</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-3 mb-3">
            <div className="md:col-span-2">
              <label className="text-zinc-200 text-xs">Function (JS expression)</label>
              <Input value={expr} onChange={(e) => setExpr(e.target.value)} className="bg-zinc-700 text-white" />
              <div className="text-zinc-400 text-xs mt-1">Example: Math.sin(x) - 0.5 or x * Math.cos(x) - 1</div>
            </div>

            <div>
              <label className="text-zinc-200 text-xs">Range</label>
              <div className="flex gap-2">
                <Input value={String(left)} onChange={(e) => setLeft(Number(e.target.value))} className="bg-zinc-700 text-white" />
                <Input value={String(right)} onChange={(e) => setRight(Number(e.target.value))} className="bg-zinc-700 text-white" />
              </div>
            </div>
          </div>

          <div className="mb-3">
            <div className="w-full overflow-hidden rounded" style={{ minHeight: 260 }}>
              <Plot
                data={[
                  { x: plot.xs, y: plot.ys, mode: "lines", name: "f(x)", line: { color: theme.accent1 } },
                  { x: [left, right], y: [0, 0], mode: "lines", name: "x-axis", line: { dash: "dot" } },
                ]}
                layout={{
                  title: `Plot of f(x) = ${expr}`,
                  paper_bgcolor: theme.panelBg,
                  plot_bgcolor: theme.panelBg,
                  font: { color: theme.text },
                  height: 320,
                  margin: { t: 40, l: 50, r: 20, b: 50 },
                }}
                useResizeHandler={true}
                style={{ width: "100%" }}
              />
            </div>
          </div>

          <div className="grid grid-cols-1 sm:grid-cols-3 gap-3">
            <SummaryCard title="Range" value={`${left} → ${right}`} />
            <SummaryCard title="Samples" value={`${samples}`} />
            <SummaryCard title="Sign intervals" value={`${signIntervals.length}`} subtitle="potential brackets where f changes sign" />
          </div>

          <div className="mt-3 text-zinc-400 text-sm">
            Graphical inspection and sign-change scanning are useful to identify initial brackets for bracketing methods (bisection & false-position). Use fine sampling to avoid missing narrow roots.
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

function BisectionPanel() {
  const [expr, setExpr] = useState("Math.sin(x) - 0.5");
  const [a, setA] = useState(0.5);
  const [b, setB] = useState(2.5);
  const [tol, setTol] = useState(1e-6);
  const [maxIter, setMaxIter] = useState(40);

  const f = useMemo(() => safeEvalFn(expr), [expr]);

  const history = useMemo(() => {
    if (!f) return [];
    const arr = [];
    let A = a;
    let B = b;
    try {
      let FA = f(A);
      let FB = f(B);
      if (!Number.isFinite(FA) || !Number.isFinite(FB) || FA * FB > 0) return [];
      arr.push({ A, B, FA, FB });
      for (let i = 0; i < maxIter; i++) {
        const C = (A + B) / 2;
        const FC = f(C);
        arr.push({ A, B, C, FC });
        if (Math.abs(FC) < tol || Math.abs(B - A) / 2 < tol) break;
        if (FA * FC <= 0) {
          B = C;
          FB = FC;
        } else {
          A = C;
          FA = FC;
        }
      }
    } catch (e) {
      return [];
    }
    return arr;
  }, [f, a, b, tol, maxIter]);

  const last = history.length ? history[history.length - 1] : null;

  const xs = useMemo(() => {
    if (!f) return { xs: [], ys: [] };
    const left = Math.min(a, b) - 0.2 * Math.abs(b - a);
    const right = Math.max(a, b) + 0.2 * Math.abs(b - a);
    const steps = 240;
    const dx = (right - left) / steps;
    const xs = [];
    const ys = [];
    for (let i = 0; i <= steps; i++) {
      const x = left + i * dx;
      xs.push(x);
      try {
        ys.push(f(x));
      } catch {
        ys.push(NaN);
      }
    }
    return { xs, ys };
  }, [f, a, b]);

  return (
    <motion.div variants={cardMotion} initial="initial" animate="enter" whileHover="hover">
      <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700 overflow-hidden`}>
        <CardHeader>
          <CardTitle className="flex items-center gap-2 text-emerald-300"><Circle className="w-5 h-5" /> 5.2 Bisection Method — interval halving steps</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-3 mb-3">
            <div className="md:col-span-2">
              <label className="text-zinc-200 text-xs">Function (JS)</label>
              <Input value={expr} onChange={(e) => setExpr(e.target.value)} className="bg-zinc-700 text-white" />
            </div>

            <div>
              <label className="text-zinc-200 text-xs">Bracket [a, b]</label>
              <div className="flex gap-2">
                <Input value={String(a)} onChange={(e) => setA(Number(e.target.value))} className="bg-zinc-700 text-white" />
                <Input value={String(b)} onChange={(e) => setB(Number(e.target.value))} className="bg-zinc-700 text-white" />
              </div>
            </div>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-3 gap-3 mb-3">
            <div>
              <label className="text-zinc-200 text-xs">Tolerance</label>
              <Input value={String(tol)} onChange={(e) => setTol(Number(e.target.value))} className="bg-zinc-700 text-white" />
            </div>
            <div>
              <label className="text-zinc-200 text-xs">Max iter</label>
              <Input value={String(maxIter)} onChange={(e) => setMaxIter(Number(e.target.value))} className="bg-zinc-700 text-white" />
            </div>
            <div className="flex items-end">
              <div className="text-zinc-400 text-sm">Last step: {last ? `c ≈ ${last.c ?? ((last.A + last.B) / 2)}` : '—'}</div>
            </div>
          </div>

          <div className="mb-3">
            <div className="w-full overflow-hidden rounded" style={{ minHeight: 260 }}>
              <Plot
                data={[
                  { x: xs.xs, y: xs.ys, mode: "lines", name: "f(x)", line: { color: theme.accent1 } },
                  ...(history.length
                    ? [
                        { x: [history[0].A, history[0].B], y: [0, 0], mode: "lines", name: "initial bracket", line: { dash: "dot" } },
                      ]
                    : []),
                  ...history
                    .slice(0, 20)
                    .map((h, i) => ({ x: [h.A ?? h.a ?? h.A, h.B ?? h.b ?? h.B], y: [0, 0], mode: "lines", name: `step ${i}`, line: { width: 2, dash: "dash" } })),
                ]}
                layout={{
                  title: `Bisection steps for ${expr}`,
                  paper_bgcolor: theme.panelBg,
                  plot_bgcolor: theme.panelBg,
                  font: { color: theme.text },
                  height: 320,
                }}
                useResizeHandler={true}
                style={{ width: "100%" }}
              />
            </div>
          </div>

          <div className="grid grid-cols-1 sm:grid-cols-3 gap-3">
            <SummaryCard title="Iterations" value={`${history.length}`} />
            <SummaryCard title="Approx root" value={last ? String(last.c ?? ((last.A + last.B) / 2)) : '—'} />
            <SummaryCard title="Final interval" value={last ? `${last.A} → ${last.B}` : '—'} />
          </div>

          <div className="mt-3 text-zinc-400 text-sm">Bisection is robust and guarantees convergence for continuous functions with sign change in the bracket. Convergence is linear (factor 1/2 per iteration).</div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

function FalsePositionPanel() {
  const [expr, setExpr] = useState("Math.sin(x) - 0.5");
  const [a, setA] = useState(0.5);
  const [b, setB] = useState(2.5);
  const [tol, setTol] = useState(1e-6);
  const [maxIter, setMaxIter] = useState(50);

  const f = useMemo(() => safeEvalFn(expr), [expr]);

  const history = useMemo(() => {
    if (!f) return [];
    const arr = [];
    let A = a;
    let B = b;
    try {
      let FA = f(A);
      let FB = f(B);
      if (!Number.isFinite(FA) || !Number.isFinite(FB) || FA * FB > 0) return [];
      arr.push({ A, B, FA, FB });
      for (let i = 0; i < maxIter; i++) {
        const C = (A * FB - B * FA) / (FB - FA);
        const FC = f(C);
        arr.push({ A, B, C, FA, FB, FC });
        if (Math.abs(FC) < tol) break;
        if (FA * FC <= 0) {
          B = C;
          FB = FC;
        } else {
          A = C;
          FA = FC;
        }
      }
    } catch (e) {
      return [];
    }
    return arr;
  }, [f, a, b, tol, maxIter]);

  const last = history.length ? history[history.length - 1] : null;

  const xs = useMemo(() => {
    if (!f) return { xs: [], ys: [] };
    const left = Math.min(a, b) - 0.2 * Math.abs(b - a);
    const right = Math.max(a, b) + 0.2 * Math.abs(b - a);
    const steps = 240;
    const dx = (right - left) / steps;
    const xs = [];
    const ys = [];
    for (let i = 0; i <= steps; i++) {
      const x = left + i * dx;
      xs.push(x);
      try {
        ys.push(f(x));
      } catch {
        ys.push(NaN);
      }
    }
    return { xs, ys };
  }, [f, a, b]);

  return (
    <motion.div variants={cardMotion} initial="initial" animate="enter" whileHover="hover">
      <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700 overflow-hidden`}>
        <CardHeader>
          <CardTitle className="flex items-center gap-2 text-amber-300"><Zap className="w-5 h-5" /> 5.3 False-Position Method — secant update (regula falsi)</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-3 mb-3">
            <div className="md:col-span-2">
              <label className="text-zinc-200 text-xs">Function (JS)</label>
              <Input value={expr} onChange={(e) => setExpr(e.target.value)} className="bg-zinc-700 text-white" />
            </div>

            <div>
              <label className="text-zinc-200 text-xs">Bracket [a, b]</label>
              <div className="flex gap-2">
                <Input value={String(a)} onChange={(e) => setA(Number(e.target.value))} className="bg-zinc-700 text-white" />
                <Input value={String(b)} onChange={(e) => setB(Number(e.target.value))} className="bg-zinc-700 text-white" />
              </div>
            </div>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-3 gap-3 mb-3">
            <div>
              <label className="text-zinc-200 text-xs">Tolerance</label>
              <Input value={String(tol)} onChange={(e) => setTol(Number(e.target.value))} className="bg-zinc-700 text-white" />
            </div>
            <div>
              <label className="text-zinc-200 text-xs">Max iter</label>
              <Input value={String(maxIter)} onChange={(e) => setMaxIter(Number(e.target.value))} className="bg-zinc-700 text-white" />
            </div>
            <div className="flex items-end">
              <div className="text-zinc-400 text-sm">Last estimate: {last ? `c ≈ ${last.c ?? last.C}` : '—'}</div>
            </div>
          </div>

          <div className="mb-3">
            <div className="w-full overflow-hidden rounded" style={{ minHeight: 260 }}>
              <Plot
                data={[
                  { x: xs.xs, y: xs.ys, mode: "lines", name: "f(x)", line: { color: theme.accent1 } },
                  ...history
                    .slice(0, 20)
                    .map((h, i) => ({
                      x: [h.A ?? h.a ?? h.A, h.B ?? h.b ?? h.B],
                      y: [0, 0],
                      mode: "lines",
                      name: `step ${i}`,
                      line: { width: 2, dash: "dash" },
                    })),
                  ...(history.length
                    ? [
                        // add secant line for last step
                        (() => {
                          const h = history[history.length - 1];
                          if (!h || !h.A || !h.B) return null;
                          const x1 = h.A;
                          const y1 = h.FA ?? f(h.A);
                          const x2 = h.B;
                          const y2 = h.FB ?? f(h.B);
                          const slope = (y2 - y1) / (x2 - x1);
                          const left = Math.min(a, b) - 0.2 * Math.abs(b - a);
                          const right = Math.max(a, b) + 0.2 * Math.abs(b - a);
                          return {
                            x: [left, right],
                            y: [y1 + slope * (left - x1), y1 + slope * (right - x1)],
                            mode: "lines",
                            name: "secant",
                            line: { dash: "dot" },
                          };
                        })(),
                      ]
                    : []),
                ].filter(Boolean)}
                layout={{
                  title: `False-position steps for ${expr}`,
                  paper_bgcolor: theme.panelBg,
                  plot_bgcolor: theme.panelBg,
                  font: { color: theme.text },
                  height: 320,
                }}
                useResizeHandler={true}
                style={{ width: "100%" }}
              />
            </div>
          </div>

          <div className="grid grid-cols-1 sm:grid-cols-3 gap-3">
            <SummaryCard title="Iterations" value={`${history.length}`} />
            <SummaryCard title="Approx root" value={last ? String(last.c ?? last.C) : '—'} />
           <SummaryCard 
  title="Converged?" 
  value={
    last 
      ? (Number.isFinite(last.FC ?? last.fc) && Math.abs(last.FC ?? last.fc) < tol ? 'yes' : 'maybe') 
      : '—'
  } 
/>

          </div>

          <div className="mt-3 text-zinc-400 text-sm">False-position can converge faster than bisection for some functions, but may stagnate (one endpoint remains fixed). Modified regula falsi (Illinois, Pegasus) addresses that by adjusting stale endpoints.</div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

function IncrementalSearchPanel() {
  const [expr, setExpr] = useState("Math.sin(x) - 0.5");
  const [left, setLeft] = useState(-6.28);
  const [right, setRight] = useState(6.28);
  const [stepSize, setStepSize] = useState(0.5);

  const f = useMemo(() => safeEvalFn(expr), [expr]);

  // safe step size to avoid 0/NaN
  const safeStep = stepSize > 0 ? stepSize : 0.1;

  const intervals = useMemo(() => {
    if (!f) return [];
    try {
      const steps = Math.max(10, Math.ceil((right - left) / safeStep));
      if (!Number.isFinite(steps) || steps > 5000) return []; // safety guard
      return findSignChanges(f, left, right, steps);
    } catch (e) {
      return [];
    }
  }, [f, left, right, safeStep]);

  return (
    <motion.div variants={cardMotion} initial="initial" animate="enter" whileHover="hover">
      <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700 overflow-hidden`}>
        <CardHeader>
          <CardTitle className="flex items-center gap-2 text-indigo-300">
            <SquarePlus className="w-5 h-5" /> 5.4 Incremental Searches & initial guesses
          </CardTitle>
        </CardHeader>
        <CardContent>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-3 mb-3">
            <div className="md:col-span-2">
              <label className="text-zinc-200 text-xs">Function</label>
              <Input
                value={expr}
                onChange={(e) => setExpr(e.target.value)}
                className="bg-zinc-700 text-white"
              />
            </div>
            <div>
              <label className="text-zinc-200 text-xs">Step size</label>
              <Input
                type="number"
                step="0.1"
                value={String(stepSize)}
                onChange={(e) => {
                  let val = Number(e.target.value);
                  if (!Number.isFinite(val) || val <= 0) val = 0.1;
                  setStepSize(val);
                }}
                className="bg-zinc-700 text-white"
              />
            </div>
          </div>

          <div className="mb-3">
            <div className="w-full overflow-hidden rounded" style={{ minHeight: 220 }}>
              <Plot
                data={[
                  // function plot
                  (() => {
                    if (!f) return null;
                    const xs = [];
                    const ys = [];
                    const steps = 360;
                    const leftR = left;
                    const rightR = right;
                    const dx = (rightR - leftR) / steps;
                    for (let i = 0; i <= steps; i++) {
                      const x = leftR + i * dx;
                      xs.push(x);
                      try {
                        ys.push(f(x));
                      } catch {
                        ys.push(NaN);
                      }
                    }
                    return { x: xs, y: ys, mode: "lines", name: "f(x)", line: { color: theme.accent1 } };
                  })(),
                  // sign-change markers
                  intervals.length
                    ? {
                        x: intervals.map((iv) => (iv[0] + iv[1]) / 2),
                        y: intervals.map(() => 0),
                        mode: "markers",
                        name: "brackets",
                        marker: { size: 8 },
                      }
                    : null,
                ].filter(Boolean)}
                layout={{
                  title: `Incremental search (step=${safeStep})`,
                  paper_bgcolor: theme.panelBg,
                  plot_bgcolor: theme.panelBg,
                  font: { color: theme.text },
                  height: 320,
                }}
                useResizeHandler={true}
                style={{ width: "100%" }}
              />
            </div>
          </div>

          <div className="grid grid-cols-1 sm:grid-cols-3 gap-3">
            <SummaryCard title="Scanned range" value={`${left} → ${right}`} />
            <SummaryCard title="Step size" value={`${safeStep}`} />
            <SummaryCard title="Found brackets" value={`${intervals.length}`} subtitle="use these as initial guesses" />
          </div>

          <div className="mt-3 text-zinc-400 text-sm">
            Incremental searches are simple and reliable to locate sign changes; choose step-size considering
            function oscillation. Combine with plotting to avoid missing narrow-root intervals.
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}


// ---------------- Documentation Panel ----------------
const docsText = {
  "5.1 Graphical Methods": [
    `Graphical methods use function plots and sign-change scans to identify intervals where roots may exist. While visual inspection is never a proof, it is invaluable for initial guess selection and to notice multiple roots, nearby singularities, or discontinuities.`,
    `Practical tips: plot over a sufficiently wide range; use denser sampling where high-frequency oscillation is suspected; mark candidate intervals with sign changes for use with bracketing methods.`,
    `Graphical root finding is often the very first step in numerical analysis because it provides intuition about the behavior of the function. Engineers and scientists frequently overlay experimental data with a theoretical model to visually check root locations.`,
    `Advantages: intuitive, immediate feedback, and helps to avoid blindly applying algorithms to intervals with no roots.`,
    `Disadvantages: resolution of the graph may hide multiple close roots, numerical precision is limited by the plotting software, and it provides no rigorous stopping criterion.`,
    `Common usage: combined with incremental searches to select better initial guesses for bisection, regula falsi, or Newton-Raphson methods.`,
  ],

  "5.2 The Bisection Method": [
    `Bisection is a robust root-finding method that requires a continuous function and an interval [a,b] such that f(a) and f(b) have opposite signs. The method halves the interval each iteration, guaranteeing convergence to a root (or to a point where the sign changes).`,
    `Convergence is linear: the error reduces by a factor of 1/2 per iteration. Stopping criteria typically use absolute interval width |b-a| < tol or |f(c)| < tol.`,
    `The algorithm is guaranteed to converge if the assumptions hold, making it one of the safest methods in numerical analysis.`,
    `The main drawback is its slow convergence compared to methods like Newton-Raphson or secant methods.`,
    `Typical pseudocode: 
       1. Check that f(a) and f(b) have opposite signs.
       2. Compute c = (a+b)/2.
       3. If f(c) = 0 or |b-a| < tol, stop.
       4. Otherwise replace a or b with c depending on the sign, and repeat.`,
    `Applications: used where guaranteed convergence is more important than speed, e.g., safety-critical systems, control algorithms, or as a reliable "first stage" before switching to a faster method.`,
    `Practical tip: always test for discontinuities (like division by zero) in the interval, since they may also produce sign changes without an actual root.`,
  ],

  "5.3 The False-Position Method (Regula Falsi)": [
    `False-position (regula falsi) replaces the midpoint with the x-intercept of the secant line through (a,f(a)) and (b,f(b)). This can accelerate convergence for nearly-linear functions but may suffer from endpoint stagnation (one endpoint remaining fixed).`,
    `Variants (Illinois, Pegasus) modify the weight of the endpoint that does not change to avoid stagnation.`,
    `Convergence is usually faster than bisection, especially when the function is close to linear in the bracketing interval. However, it is not guaranteed to converge faster in all cases.`,
    `The formula for the new estimate is: 
       c = b - f(b)(b-a) / (f(b) - f(a))`,
    `If f(c) = 0 or the interval width is smaller than tolerance, the process stops.`,
    `Problem of stagnation: if the function is strongly curved, one endpoint may stay unchanged through many iterations. In such cases, the algorithm effectively reduces to bisection in performance.`,
    `Illinois method: reduces the function value at the stagnant endpoint by half each time it is reused, forcing interval update.`,
    `Pegasus method: a refined variant that rescales weights more smoothly, improving robustness.`,
    `Usage: suitable for engineering calculations where bisection is too slow but Newton’s method may fail due to poor derivatives.`,
  ],

  "5.4 Incremental Searches and Initial Guesses": [
    `An incremental search (scanning) partitions an interval and checks for sign changes between adjacent sample points to form candidate brackets. The step size controls sensitivity: too large and roots may be missed; too small and computation becomes expensive.`,
    `Use adaptive refinement: when a candidate interval is detected, subdivide it to isolate narrow roots before applying bracketing methods.`,
    `This technique is simple to implement and requires no derivatives, making it useful as a preprocessing stage before more sophisticated methods.`,
    `Advantages: guarantees that no root is overlooked in the scanned region, provided the step size is fine enough.`,
    `Disadvantages: computationally expensive for large ranges or highly oscillatory functions, and cannot directly estimate the root location with high accuracy.`,
    `Heuristics: 
       - Start with a coarse step size to cover the whole domain. 
       - Refine step size only in intervals where sign changes occur. 
       - Combine with plotting for additional confirmation.`,
    `In engineering, incremental searches are often applied to detect natural frequencies, resonance peaks, or load thresholds, where multiple roots may exist.`,
    `Mathematical note: a function with many oscillations may require logarithmic or adaptive step selection to avoid missing closely spaced roots.`,
  ],
};


function ChapterDocs() {
  const [open, setOpen] = useState("5.1 Graphical Methods");

  return (
    <div className="w-full mt-6 bg-zinc-900 border border-zinc-700 rounded-xl shadow-lg overflow-hidden">
      <div className="p-3 border-b border-zinc-700 flex items-center justify-between">
        <div className="flex items-center gap-3">
          <BookOpen className="w-5 h-5 text-cyan-400" />
          <div>
            <div className="text-zinc-100 font-semibold">Chapter 5 — Documentation</div>
            <div className="text-zinc-400 text-xs">Theory, practical notes, pseudocode, and usage tips for each subsection</div>
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
            <div className="mt-3 text-zinc-400 text-sm">References: Burden & Faires, Atkinson, and practical numerical methods notes about robustness and error behaviour in root-finding algorithms.</div>
          </div>
        </div>
      </div>
    </div>
  );
}

// ---------------- Problems ----------------
function ProblemsCard() {
  return (
    <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700`}>
      <CardHeader>
        <CardTitle className="text-emerald-300 flex items-center gap-2"><List className="w-5 h-5" /> Problems & Exercises</CardTitle>
      </CardHeader>
      <CardContent>
        <ul className="list-decimal list-inside text-zinc-300 text-sm space-y-2">
          <li>Use the bisection method to locate a root of <code className="font-mono">f(x)=x^3 - 2x - 5</code> in [1, 3]. Show iterations until |b-a| &lt; 10^{-6} and report the approximate root.</li>
          <li>For <code className="font-mono">f(x)=cos(x)-x</code>, apply false-position and bisection from the same initial bracket and compare the number of iterations required to reach <code className="font-mono">|f(x)| &lt; 10^{-8}</code>.</li>
          <li>Demonstrate how an incremental search with step 0.1 on [-10,10] finds candidate brackets for <code className="font-mono">f(x)=Math.sin(5*x)/(1+x^2)</code>. Discuss missed roots if step is doubled.</li>
          <li>Implement the Illinois or Pegasus modification to regula falsi and show numerical evidence that it reduces endpoint stagnation for <code className="font-mono">f(x)=x^3 - 2x + 2</code>.</li>
        </ul>
      </CardContent>
    </Card>
  );
}

// ---------------- MAIN PAGE ----------------
export default function Chapter5() {
  return (
    <div className={`min-h-screen px-4 md:px-6 py-6 bg-zinc-950 text-zinc-100`}>
      <motion.header initial={{ opacity: 0, y: -8 }} animate={{ opacity: 1, y: 0 }} className="max-w-7xl mx-auto mb-6">
        <div className="flex flex-col md:flex-row items-start justify-between gap-4">
          <div className="max-w-3xl">
            <h1 className="text-2xl text-emerald-400 md:text-3xl font-bold">Bracketing Methods — Bisection, False-Position, and Incremental Search</h1>
            <p className="text-zinc-300 mt-2">Interactive visualizations for root finding: use graphical inspection and incremental scanning to form initial brackets, then apply bisection or false-position. Experiment with tolerances, step sizes, and visualize iteration history.</p>
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
              <div className="text-zinc-300 text-sm">This chapter covers graphical root inspection, bisection, false-position (regula falsi), and incremental search for initial guesses. Use the interactive panels to experiment and internalize the methods.</div>
            </CardContent>
          </Card>

          <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700`}>
            <CardHeader>
              <CardTitle className="text-emerald-300 flex items-center gap-2"><Layers className="w-5 h-5" /> Visual intuition</CardTitle>
            </CardHeader>
            <CardContent>
              <div className="text-zinc-300 text-sm">Plots show function geometry, candidate brackets, and iteration progress. Use the 3D root cloud to visualise multiple roots in a compact form.</div>
            </CardContent>
          </Card>

          <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700`}>
            <CardHeader>
              <CardTitle className="text-amber-300 flex items-center gap-2"><AlertTriangle className="w-5 h-5" /> Practical notes</CardTitle>
            </CardHeader>
            <CardContent>
              <ul className="text-zinc-300 text-sm list-disc list-inside">
                <li>Always start with a graphical scan when the function is cheap to evaluate.</li>
                <li>Prefer bisection when robustness is critical; use false-position for speed if stagnation is unlikely.</li>
                <li>Check endpoints for multiple roots or derivative near-zero (flat regions) which can slow convergence.</li>
              </ul>
            </CardContent>
          </Card>
        </div>

        {/* Main content: two-column responsive layout */}
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          <div className="space-y-6">
            <GraphicalPanel />
            <BisectionPanel />
          </div>

          <div className="space-y-6">
            <FalsePositionPanel />
            <IncrementalSearchPanel />
          </div>
        </div>

        <ChapterDocs />

        <ProblemsCard />

        {/* 3D visual small canvas */}
        <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700`}>
          <CardHeader>
            <CardTitle className="text-cyan-300 flex items-center gap-2"><Layers className="w-5 h-5" /> Root cloud</CardTitle>
          </CardHeader>
          <CardContent>
            <div className="w-full h-48">
              <Canvas camera={{ position: [0, 0, 4] }}>
                <ambientLight intensity={0.6} />
                <pointLight position={[10, 10, 10]} />
                <RootCloud roots={[ -1.0, 0.5, 1.8 ]} />
                <OrbitControls />
              </Canvas>
            </div>
            <div className="text-zinc-400 text-sm mt-3">A compact 3D hinting widget — not for computation, but to add visual depth to the page.</div>
          </CardContent>
        </Card>
   <BottomBar/>

   
      </main>
    </div>
  );
}
