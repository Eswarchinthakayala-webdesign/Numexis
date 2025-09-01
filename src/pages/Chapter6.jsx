// src/pages/Chapter6_OpenMethods.jsx
// Professional Chapter 6 page — Open Methods (Fixed-Point, Newton, Secant, Brent, Multiple Roots, Systems)
// Uses: react, framer-motion, react-plotly.js, @react-three/fiber, @react-three/drei, lucide-react
// Tailwind + shadcn-style components assumed (Card, CardHeader, CardTitle, CardContent, Input, Button)

import React, { useState, useMemo, useRef } from "react";
import { motion } from "framer-motion";
import Plot from "react-plotly.js";
import { Canvas, useFrame } from "@react-three/fiber";
import { OrbitControls } from "@react-three/drei";
import ReactMarkdown from "react-markdown";
import remarkMath from "remark-math";
import rehypeKatex from "rehype-katex";
import "katex/dist/katex.min.css";   // <- import KaTeX CSS once globally

import {
  Target,
  Circle,
  Slash,
  RotateCw,
  Zap,
  BookOpen,
  List,
  ChevronDown,
  ChevronRight,
  Layers,
} from "lucide-react";

import { Card, CardHeader, CardTitle, CardContent } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Button } from "@/components/ui/button";
import BottomBar from "../components/BottomBar";

// ---------------- Theme & motion ----------------
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

// ---------------- Utilities ----------------
function safeEvalFn(expr, args = ["x"]) {
  // returns a function or null; sandboxed-ish (still use with care)
  try {
    // eslint-disable-next-line no-new-func
    if (args.length === 1) return new Function(args[0], `with (Math) { return (${expr}); }`);
    // for systems: x,y
    return new Function(...args, `with (Math) { return (${expr}); }`);
  } catch {
    return null;
  }
}

function approxDerivative(f, x, h = 1e-6) {
  try {
    return (f(x + h) - f(x - h)) / (2 * h);
  } catch {
    return NaN;
  }
}

const clamp = (v, lo, hi) => Math.min(Math.max(v, lo), hi);

// ---------------- Small 3D particle for polish ----------------
function ParticleCloud({ points = 48, size = 0.05 }) {
  const ref = useRef();
  useFrame(({ clock }) => {
    if (!ref.current) return;
    ref.current.rotation.y = clock.getElapsedTime() * 0.08;
    ref.current.rotation.x = Math.sin(clock.getElapsedTime() * 0.03) * 0.12;
  });

  const verts = useMemo(() => {
    const arr = [];
    for (let i = 0; i < points; i++) {
      const t = i / points;
      const phi = t * Math.PI * 2;
      const r = 1 + 0.25 * Math.sin(4 * phi);
      arr.push(r * Math.cos(phi), r * Math.sin(phi), 0.15 * Math.cos(3 * phi));
    }
    return new Float32Array(arr);
  }, [points]);

  return (
    <points ref={ref}>
      <bufferGeometry>
        <bufferAttribute attach="attributes-position" array={verts} count={verts.length / 3} itemSize={3} />
      </bufferGeometry>
      <pointsMaterial size={size} sizeAttenuation color={theme.accent2} />
    </points>
  );
}

// ---------------- Small summary card ----------------
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

/* 6.1 Fixed-Point Iteration */
function FixedPointPanel() {
  const [gExpr, setGExpr] = useState("(Math.cos(x) + x) / 2");
  const [x0, setX0] = useState(0.5);
  const [maxIter, setMaxIter] = useState(50);
  const [tol, setTol] = useState(1e-8);

  const g = useMemo(() => safeEvalFn(gExpr, ["x"]), [gExpr]);

  const history = useMemo(() => {
    if (!g) return [];
    const arr = [];
    let x = x0;
    arr.push(x);
    for (let i = 0; i < maxIter; i++) {
      try {
        const nx = g(x);
        if (!Number.isFinite(nx)) break;
        arr.push(nx);
        if (Math.abs(nx - x) < tol) break;
        x = nx;
      } catch {
        break;
      }
    }
    return arr;
  }, [g, x0, maxIter, tol]);

  const xsPlot = useMemo(() => {
    // plot g(x) and y=x between -6..6 but trimmed if g is invalid
    const xs = [];
    const ys = [];
    for (let i = -6; i <= 6; i += 0.05) {
      xs.push(i);
      try {
        ys.push(g ? g(i) : NaN);
      } catch {
        ys.push(NaN);
      }
    }
    return { xs, ys };
  }, [g]);

  return (
    <motion.div variants={cardMotion} initial="initial" animate="enter" whileHover="hover">
      <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700 overflow-hidden`}>
        <CardHeader>
          <CardTitle className="flex items-center gap-2 text-cyan-300">
            <Target className="w-5 h-5" /> 6.1 Fixed-Point Iteration
          </CardTitle>
        </CardHeader>

        <CardContent>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-3 mb-3">
            <div className="md:col-span-2">
              <label className="text-zinc-200 text-xs">g(x) (JS expression)</label>
              <Input value={gExpr} onChange={(e) => setGExpr(e.target.value)} className="bg-zinc-700 text-white" />
              <div className="text-zinc-400 text-xs mt-1">Example: (Math.cos(x)+x)/2</div>
            </div>

            <div>
              <label className="text-zinc-200 text-xs">x₀</label>
              <Input value={String(x0)} onChange={(e) => setX0(Number(e.target.value))} className="bg-zinc-700 text-white" />
            </div>
          </div>

          <div className="mb-3">
            <div className="w-full overflow-hidden rounded" style={{ minHeight: 260 }}>
              <Plot
                data={[
                  { x: xsPlot.xs, y: xsPlot.ys, mode: "lines", name: "g(x)", line: { color: theme.accent1 } },
                  { x: xsPlot.xs, y: xsPlot.xs, mode: "lines", name: "y=x", line: { dash: "dot" } },
                  history.length > 0 ? { x: history, y: history, mode: "markers+lines", name: "iterates", marker: { size: 6 } } : null,
                ].filter(Boolean)}
                layout={{
                  title: `Fixed-point iterates`,
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
            <SummaryCard title="Last iterate" value={`${history.length ? history[history.length - 1] : "—"}`} />
            <SummaryCard
              title="Converged?"
              value={history.length > 1 && Math.abs(history[history.length - 1] - history[history.length - 2]) < tol ? "yes" : "maybe"}
            />
          </div>

          <div className="mt-3 text-zinc-400 text-sm">
            Fixed-point iteration requires selecting g(x) to be contractive near the root (|g'(r)| &lt; 1). If slow, consider Aitken acceleration or rearranging the equation.
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

/* 6.2 Newton–Raphson */
function NewtonPanel() {
  const [expr, setExpr] = useState("Math.cos(x) - x");
  const [x0, setX0] = useState(0.7);
  const [maxIter, setMaxIter] = useState(40);
  const [tol, setTol] = useState(1e-10);

  const f = useMemo(() => safeEvalFn(expr, ["x"]), [expr]);

  const history = useMemo(() => {
    if (!f) return [];
    const arr = [];
    let x = x0;
    for (let i = 0; i < maxIter; i++) {
      try {
        const fx = f(x);
        const dfx = approxDerivative(f, x);
        if (!Number.isFinite(dfx) || Math.abs(dfx) < 1e-14) break;
        const nx = x - fx / dfx;
        arr.push({ x: nx, fx: f(nx), df: approxDerivative(f, nx) });
        if (Math.abs(nx - x) < tol) break;
        x = nx;
      } catch {
        break;
      }
    }
    return arr;
  }, [f, x0, maxIter, tol]);

  const xsPlot = useMemo(() => {
    const xs = [];
    const ys = [];
    for (let i = -5; i <= 5; i += 0.02) {
      xs.push(i);
      try {
        ys.push(f ? f(i) : NaN);
      } catch {
        ys.push(NaN);
      }
    }
    return { xs, ys };
  }, [f]);

  return (
    <motion.div variants={cardMotion} initial="initial" animate="enter" whileHover="hover">
      <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700 overflow-hidden`}>
        <CardHeader>
          <CardTitle className="flex items-center gap-2 text-amber-300">
            <Circle className="w-5 h-5" /> 6.2 Newton–Raphson Method
          </CardTitle>
        </CardHeader>

        <CardContent>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-3 mb-3">
            <div className="md:col-span-2">
              <label className="text-zinc-200 text-xs">f(x) (JS expression)</label>
              <Input value={expr} onChange={(e) => setExpr(e.target.value)} className="bg-zinc-700 text-white" />
              <div className="text-zinc-400 text-xs mt-1">Example: Math.cos(x) - x</div>
            </div>
            <div>
              <label className="text-zinc-200 text-xs">x₀</label>
              <Input value={String(x0)} onChange={(e) => setX0(Number(e.target.value))} className="bg-zinc-700 text-white" />
            </div>
          </div>

          <div className="mb-3">
            <div className="w-full overflow-hidden rounded" style={{ minHeight: 320 }}>
              <Plot
                data={[
                  { x: xsPlot.xs, y: xsPlot.ys, mode: "lines", name: "f(x)", line: { color: theme.accent1 } },
                  ...history.slice(0, 6).map((h, i) => {
                    const x = h.x;
                    const slope = h.df || approxDerivative(f, x);
                    const left = x - 1;
                    const right = x + 1;
                    return {
                      x: [left, right],
                      y: [h.fx + slope * (left - x), h.fx + slope * (right - x)],
                      mode: "lines",
                      name: `tangent ${i}`,
                      line: { dash: "dash" },
                    };
                  }),
                  history.length > 0 ? { x: history.map((h) => h.x), y: history.map((h) => h.fx), mode: "markers+lines", name: "iterates", marker: { size: 6 } } : null,
                ].filter(Boolean)}
                layout={{
                  title: `Newton iteration`,
                  paper_bgcolor: theme.panelBg,
                  plot_bgcolor: theme.panelBg,
                  font: { color: theme.text },
                  height: 380,
                }}
                useResizeHandler={true}
                style={{ width: "100%" }}
              />
            </div>
          </div>

          <div className="grid grid-cols-1 sm:grid-cols-3 gap-3">
            <SummaryCard title="Iterations" value={`${history.length}`} />
            <SummaryCard title="Approx root" value={`${history.length ? history[history.length - 1].x : "—"}`} />
            <SummaryCard
              title="Last |f(x)|"
              value={`${history.length ? Math.abs(history[history.length - 1].fx).toExponential(3) : "—"}`}
            />
          </div>

          <div className="mt-3 text-zinc-400 text-sm">
            Newton's method converges quadratically for simple roots but requires reliable derivative information. Use damping or fallback if divergence is observed.
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

/* 6.3 Secant Method */
function SecantPanel() {
  const [expr, setExpr] = useState("Math.cos(x) - x");
  const [x0, setX0] = useState(0.5);
  const [x1, setX1] = useState(1.0);
  const [maxIter, setMaxIter] = useState(50);
  const [tol, setTol] = useState(1e-8);

  const f = useMemo(() => safeEvalFn(expr, ["x"]), [expr]);

  const history = useMemo(() => {
    if (!f) return [];
    const arr = [];
    let a = x0;
    let b = x1;
    try {
      arr.push({ x: a, fx: f(a) }, { x: b, fx: f(b) });
    } catch {
      return arr;
    }
    for (let i = 0; i < maxIter; i++) {
      try {
        const fa = f(a);
        const fb = f(b);
        const denom = fb - fa;
        if (!Number.isFinite(denom) || Math.abs(denom) < 1e-14) break;
        const c = b - fb * (b - a) / denom;
        const fc = f(c);
        arr.push({ x: c, fx: fc });
        if (Math.abs(c - b) < tol) break;
        a = b;
        b = c;
      } catch {
        break;
      }
    }
    return arr;
  }, [f, x0, x1, maxIter, tol]);

  const xsPlot = useMemo(() => {
    const xs = [];
    const ys = [];
    for (let i = -5; i <= 5; i += 0.02) {
      xs.push(i);
      try {
        ys.push(f ? f(i) : NaN);
      } catch {
        ys.push(NaN);
      }
    }
    return { xs, ys };
  }, [f]);

  return (
    <motion.div variants={cardMotion} initial="initial" animate="enter" whileHover="hover">
      <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700 overflow-hidden`}>
        <CardHeader>
          <CardTitle className="flex items-center gap-2 text-emerald-300">
            <Slash className="w-5 h-5" /> 6.3 Secant Method
          </CardTitle>
        </CardHeader>

        <CardContent>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-3 mb-3">
            <div className="md:col-span-2">
              <label className="text-zinc-200 text-xs">f(x) (JS expression)</label>
              <Input value={expr} onChange={(e) => setExpr(e.target.value)} className="bg-zinc-700 text-white" />
            </div>
            <div>
              <label className="text-zinc-200 text-xs">x0, x1</label>
              <div className="flex gap-2">
                <Input value={String(x0)} onChange={(e) => setX0(Number(e.target.value))} className="bg-zinc-700 text-white" />
                <Input value={String(x1)} onChange={(e) => setX1(Number(e.target.value))} className="bg-zinc-700 text-white" />
              </div>
            </div>
          </div>

          <div className="mb-3">
            <div className="w-full overflow-hidden rounded" style={{ minHeight: 320 }}>
              <Plot
                data={[
                  { x: xsPlot.xs, y: xsPlot.ys, mode: "lines", name: "f(x)", line: { color: theme.accent1 } },
                  // secant lines (first few)
                  ...(() => {
                    const secants = [];
                    for (let i = 0; i < Math.min(5, history.length - 1); i++) {
                      const p1 = history[i];
                      const p2 = history[i + 1];
                      if (!p1 || !p2) continue;
                      const x1 = p1.x;
                      const y1 = p1.fx;
                      const x2 = p2.x;
                      const y2 = p2.fx;
                      if (!Number.isFinite(x1) || !Number.isFinite(x2)) continue;
                      const slope = (y2 - y1) / (x2 - x1 || 1e-12);
                      const left = Math.min(x1, x2) - 1;
                      const right = Math.max(x1, x2) + 1;
                      secants.push({
                        x: [left, right],
                        y: [y1 + slope * (left - x1), y1 + slope * (right - x1)],
                        mode: "lines",
                        name: `secant ${i}`,
                        line: { dash: "dash" },
                      });
                    }
                    return secants;
                  })(),
                  history.length > 0 ? { x: history.map((h) => h.x), y: history.map((h) => h.fx), mode: "markers+lines", name: "iterates", marker: { size: 6 } } : null,
                ].filter(Boolean)}
                layout={{
                  title: `Secant iterates`,
                  paper_bgcolor: theme.panelBg,
                  plot_bgcolor: theme.panelBg,
                  font: { color: theme.text },
                  height: 380,
                }}
                useResizeHandler={true}
                style={{ width: "100%" }}
              />
            </div>
          </div>

          <div className="grid grid-cols-1 sm:grid-cols-3 gap-3">
            <SummaryCard title="Iterations" value={`${history.length}`} />
            <SummaryCard title="Last x" value={`${history.length ? history[history.length - 1].x : "—"}`} />
            <SummaryCard title="Last |f(x)|" value={`${history.length ? Math.abs(history[history.length - 1].fx).toExponential(3) : "—"}`} />
          </div>

          <div className="mt-3 text-zinc-400 text-sm">
            Secant is derivative-free and often faster than bisection but lacks guaranteed global convergence. Use hybrid methods if robustness is required.
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

/* 6.4 Brent's Method (robust hybrid) */
function BrentsPanel() {
  const [expr, setExpr] = useState("Math.sin(x) - 0.5");
  const [aInput, setAInput] = useState(0.0);
  const [bInput, setBInput] = useState(2.0);
  const [tol, setTol] = useState(1e-10);
  const [maxIter, setMaxIter] = useState(80);

  const f = useMemo(() => safeEvalFn(expr, ["x"]), [expr]);

  // Implement guarded Brent-like algorithm (clean, avoid const reassignment)
  const history = useMemo(() => {
    if (!f) return [];
    // validate inputs
    let a = Number(aInput);
    let b = Number(bInput);
    if (!Number.isFinite(a) || !Number.isFinite(b) || a === b) return [];
    try {
      let fa = f(a);
      let fb = f(b);
      if (!Number.isFinite(fa) || !Number.isFinite(fb) || fa * fb > 0) return [];
      // initialize
      let c = a;
      let fc = fa;
      let s = b;
      let fs = fb;
      let d = b - a;
      let e = d;
      const hist = [];
      hist.push({ step: 0, a, b, c, fa, fb, fc, method: "init" });

      for (let iter = 1; iter <= maxIter; iter++) {
        // ensure |fb| <= |fa|
        if (Math.abs(fc) < Math.abs(fb)) {
          // rotate variables: a<-b, b<-c, c<-a (use let, reassign)
          const ta = a;
          a = b;
          b = c;
          c = ta;
          const tfa = fa;
          fa = fb;
          fb = fc;
          fc = tfa;
        }

        // convergence test
        const tol1 = 2 * Number.EPSILON * Math.abs(b) + tol / 2;
        const m = (c - b) / 2;
        hist.push({ step: iter, a, b, c, fa, fb, fc, method: "check" });

        if (Math.abs(m) <= tol1 || fb === 0) {
          break;
        }

        // attempt interpolation if useful
        let useBisection = false;
        let p = 0;
        let q = 1;

        if (Math.abs(e) >= tol1 && Math.abs(fa) > Math.abs(fb)) {
          // try inverse quadratic interpolation or secant
          if (a === c) {
            // secant
            p = (b - a) * fb;
            q = fb - fa;
          } else {
            // inverse quadratic interpolation (IQI)
            const r = fb / fc;
            const s1 = fb / fa;
            const s2 = fa / fc;
            // formula components (robust variant)
            p = s1 * ( (b - a) * ( (s2 - r) ) ) - (b - a) * (s1 - 1);
            // fallback secant if formula weird - instead use the safer typical IQI combination
            // We'll use a pragmatic secant-style expression to avoid instability
            p = ( (b * fa * fc) / ((fa - fb) * (fa - fc)) ) +
                ( (a * fb * fc) / ((fb - fa) * (fb - fc)) ) +
                ( (c * fa * fb) / ((fc - fa) * (fc - fb)) );
            // When in doubt, fallback below
          }
        } else {
          useBisection = true;
        }

        // compute s candidate; we choose to fall back to bisection for safety if candidate invalid
        let candidate = b - (fb * (b - a)) / (fb - fa || (fb - fa === 0 ? 1e-12 : 1));
        // bounds check for candidate
        if (Number.isFinite(candidate) === false || candidate <= Math.min(a, b) || candidate >= Math.max(a, b)) {
          candidate = (a + b) / 2;
          useBisection = true;
        }

        // accept candidate if interpolation plausible; else bisection
        let sNext = candidate;
        let sMethod = useBisection ? "bisection" : "interpolate";

        // ensure progress: if s is too close to previous, use bisection
        if (!useBisection) {
          const cond = Math.abs(sNext - b) >= Math.abs(b - c) / 2;
          if (cond) {
            sNext = (a + b) / 2;
            sMethod = "bisection";
          }
        }

        // evaluate sNext
        try {
          fs = f(sNext);
        } catch {
          fs = NaN;
        }

        hist.push({ step: iter, s: sNext, fs, method: sMethod });

        // update brackets
        if (Number.isFinite(fs) && fa * fs < 0) {
          b = sNext;
          fb = fs;
        } else {
          a = sNext;
          fa = fs;
        }

        // keep c as previous b to maintain bracketing geometry
        if (Math.abs(fa) < Math.abs(fb)) {
          // swap to ensure |fb| <= |fa|
          const ta = a;
          a = b;
          b = ta;
          const tfa = fa;
          fa = fb;
          fb = tfa;
        }

        // check termination
        if (Math.abs(b - a) < tol) break;
      }

      return hist;
    } catch {
      return [];
    }
  }, [f, aInput, bInput, tol, maxIter]);

  // produce plot points
  const xsPlot = useMemo(() => {
    const xs = [];
    const ys = [];
    const a = Number(aInput);
    const b = Number(bInput);
    if (!Number.isFinite(a) || !Number.isFinite(b)) return { xs, ys };
    const left = Math.min(a, b) - Math.abs(b - a) * 0.5 - 0.5;
    const right = Math.max(a, b) + Math.abs(b - a) * 0.5 + 0.5;
    const steps = 400;
    const dx = (right - left) / steps;
    for (let i = 0; i <= steps; i++) {
      const x = left + i * dx;
      xs.push(x);
      try {
        ys.push(f ? f(x) : NaN);
      } catch {
        ys.push(NaN);
      }
    }
    return { xs, ys };
  }, [f, aInput, bInput]);

  const lastStep = history.length ? history[history.length - 1] : null;

  return (
    <motion.div variants={cardMotion} initial="initial" animate="enter" whileHover="hover">
      <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700 overflow-hidden`}>
        <CardHeader>
          <CardTitle className="flex items-center gap-2 text-indigo-300">
            <RotateCw className="w-5 h-5" /> 6.4 Brent-style Hybrid (demo)
          </CardTitle>
        </CardHeader>

        <CardContent>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-3 mb-3">
            <div className="md:col-span-2">
              <label className="text-zinc-200 text-xs">f(x)</label>
              <Input value={expr} onChange={(e) => setExpr(e.target.value)} className="bg-zinc-700 text-white" />
            </div>

            <div>
              <label className="text-zinc-200 text-xs">Bracket [a, b]</label>
              <div className="flex gap-2">
                <Input value={String(aInput)} onChange={(e) => setAInput(Number(e.target.value))} className="bg-zinc-700 text-white" />
                <Input value={String(bInput)} onChange={(e) => setBInput(Number(e.target.value))} className="bg-zinc-700 text-white" />
              </div>
            </div>
          </div>

          <div className="mb-3">
            <div className="w-full overflow-hidden rounded" style={{ minHeight: 320 }}>
              <Plot
                data={[
                  { x: xsPlot.xs, y: xsPlot.ys, mode: "lines", name: "f(x)", line: { color: theme.accent1 } },
                  // show s markers:
                  ...history
                    .filter((h) => h.s !== undefined)
                    .map((h, i) => ({ x: [h.s], y: [h.fs], mode: "markers", name: `s${i}`, marker: { size: 8 } })),
                ].filter(Boolean)}
                layout={{
                  title: `Brent-style updates (demo)`,
                  paper_bgcolor: theme.panelBg,
                  plot_bgcolor: theme.panelBg,
                  font: { color: theme.text },
                  height: 380,
                }}
                useResizeHandler={true}
                style={{ width: "100%" }}
              />
            </div>
          </div>

          <div className="grid grid-cols-1 sm:grid-cols-3 gap-3">
            <SummaryCard title="Steps" value={`${history.length}`} />
            <SummaryCard title="Last s" value={`${lastStep ? (lastStep.s ?? "—") : "—"}`} />
            <SummaryCard
              title="Method mix"
              value={(() => {
                const counts = history.reduce((acc, cur) => {
                  const m = cur.method || "unknown";
                  acc[m] = (acc[m] || 0) + 1;
                  return acc;
                }, {});
                return Object.entries(counts)
                  .map(([k, v]) => `${k}:${v}`)
                  .join(" ");
              })()}
            />
          </div>

          <div className="mt-3 text-zinc-400 text-sm">
            This is a demonstration of Brent-like logic: interpolation when safe, bisection fallback for safety. Production Brent implementations include careful edge-case bookkeeping and performance optimizations.
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

/* 6.5 Multiple Roots */
function MultipleRootsPanel() {
  const [expr, setExpr] = useState("(x - 1) * (x - 1) * (x + 2)"); // (x-1)^2 (x+2)
  const [x0, setX0] = useState(0.9);
  const [m, setM] = useState(2);
  const [maxIter, setMaxIter] = useState(40);
  const [tol, setTol] = useState(1e-8);

  const f = useMemo(() => safeEvalFn(expr, ["x"]), [expr]);

  const newtonHist = useMemo(() => {
    if (!f) return [];
    const arr = [];
    let x = x0;
    for (let i = 0; i < maxIter; i++) {
      try {
        const fx = f(x);
        const dfx = approxDerivative(f, x);
        if (!Number.isFinite(dfx) || Math.abs(dfx) < 1e-14) break;
        const nx = x - fx / dfx;
        arr.push({ x: nx, fx: f(nx) });
        if (Math.abs(nx - x) < tol) break;
        x = nx;
      } catch {
        break;
      }
    }
    return arr;
  }, [f, x0, maxIter, tol]);

  const modNewtonHist = useMemo(() => {
    if (!f) return [];
    const arr = [];
    let x = x0;
    for (let i = 0; i < maxIter; i++) {
      try {
        const fx = f(x);
        const dfx = approxDerivative(f, x);
        if (!Number.isFinite(dfx) || Math.abs(dfx) < 1e-14) break;
        const nx = x - (m * fx) / dfx;
        arr.push({ x: nx, fx: f(nx) });
        if (Math.abs(nx - x) < tol) break;
        x = nx;
      } catch {
        break;
      }
    }
    return arr;
  }, [f, x0, m, maxIter, tol]);

  return (
    <motion.div variants={cardMotion} initial="initial" animate="enter" whileHover="hover">
      <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700 overflow-hidden`}>
        <CardHeader>
          <CardTitle className="flex items-center gap-2 text-indigo-300">
            <Zap className="w-5 h-5" /> 6.5 Multiple Roots — Newton vs Modified
          </CardTitle>
        </CardHeader>

        <CardContent>
          <div className="grid grid-cols-1 md:grid-cols-4 gap-3 mb-3">
            <div className="md:col-span-2">
              <label className="text-zinc-200 text-xs">f(x)</label>
              <Input value={expr} onChange={(e) => setExpr(e.target.value)} className="bg-zinc-700 text-white" />
            </div>
            <div>
              <label className="text-zinc-200 text-xs">x₀</label>
              <Input value={String(x0)} onChange={(e) => setX0(Number(e.target.value))} className="bg-zinc-700 text-white" />
            </div>
            <div>
              <label className="text-zinc-200 text-xs">multiplicity m</label>
              <Input value={String(m)} onChange={(e) => setM(Number(e.target.value))} className="bg-zinc-700 text-white" />
            </div>
          </div>

          <div className="mb-3">
            <div className="w-full overflow-hidden rounded" style={{ minHeight: 320 }}>
              <Plot
                data={[
                  {
                    x: Array.from({ length: 400 }, (_, i) => -3 + (i / 399) * 6),
                    y: Array.from({ length: 400 }, (_, i) => {
                      const x = -3 + (i / 399) * 6;
                      try {
                        return f ? f(x) : NaN;
                      } catch {
                        return NaN;
                      }
                    }),
                    mode: "lines",
                    name: "f(x)",
                    line: { color: theme.accent1 },
                  },
                  newtonHist.length > 0 ? { x: newtonHist.map((h) => h.x), y: newtonHist.map((h) => h.fx), mode: "markers+lines", name: "Newton", marker: { size: 6 } } : null,
                  modNewtonHist.length > 0 ? { x: modNewtonHist.map((h) => h.x), y: modNewtonHist.map((h) => h.fx), mode: "markers+lines", name: "Modified Newton", marker: { size: 6 } } : null,
                ].filter(Boolean)}
                layout={{
                  title: "Multiple roots — Newton vs Modified",
                  paper_bgcolor: theme.panelBg,
                  plot_bgcolor: theme.panelBg,
                  font: { color: theme.text },
                  height: 380,
                }}
                useResizeHandler={true}
                style={{ width: "100%" }}
              />
            </div>
          </div>

          <div className="grid grid-cols-1 sm:grid-cols-3 gap-3">
            <SummaryCard title="Newton iters" value={`${newtonHist.length}`} />
            <SummaryCard title="Modified iters" value={`${modNewtonHist.length}`} />
            <SummaryCard title="Last Newton x" value={`${newtonHist.length ? newtonHist[newtonHist.length - 1].x : "—"}`} />
          </div>

          <div className="mt-3 text-zinc-400 text-sm">
            For roots with multiplicity &gt;1, standard Newton can degrade to linear convergence. Modified Newton multiplies the correction by multiplicity m to recover faster behaviour (if m known).
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

/* 6.6 Systems of Nonlinear Equations (2D demo) */
function SystemsPanel() {
  const [fxExpr, setFxExpr] = useState("x*x + y*y - 4");
  const [gyExpr, setGyExpr] = useState("x - y");
  const [x0, setX0] = useState(1.0);
  const [y0, setY0] = useState(1.0);
  const [maxIter, setMaxIter] = useState(12);
  const [tol, setTol] = useState(1e-8);

  const f = useMemo(() => safeEvalFn(fxExpr, ["x", "y"]), [fxExpr]);
  const g = useMemo(() => safeEvalFn(gyExpr, ["x", "y"]), [gyExpr]);

  function jacobian(x, y) {
    const h = 1e-6;
    try {
      const fx = f(x, y);
      const gyv = g(x, y);
      const dfdx = (f(x + h, y) - f(x - h, y)) / (2 * h);
      const dfdy = (f(x, y + h) - f(x, y - h)) / (2 * h);
      const dgdx = (g(x + h, y) - g(x - h, y)) / (2 * h);
      const dgdy = (g(x, y + h) - g(x, y - h)) / (2 * h);
      return { J: [[dfdx, dfdy], [dgdx, dgdy]], F: [fx, gyv] };
    } catch {
      return null;
    }
  }

  const history = useMemo(() => {
    if (!f || !g) return [];
    const arr = [];
    let x = x0;
    let y = y0;
    arr.push({ x, y });
    for (let i = 0; i < maxIter; i++) {
      const jac = jacobian(x, y);
      if (!jac) break;
      const { J, F } = jac;
      const [dfdx, dfdy] = J[0];
      const [dgdx, dgdy] = J[1];
      const det = dfdx * dgdy - dfdy * dgdx;
      if (!Number.isFinite(det) || Math.abs(det) < 1e-14) break;
      // inverse J times F (2x2)
      const inv = [
        [dgdy / det, -dfdy / det],
        [-dgdx / det, dfdx / det],
      ];
      const dx = inv[0][0] * F[0] + inv[0][1] * F[1];
      const dy = inv[1][0] * F[0] + inv[1][1] * F[1];
      const nx = x - dx;
      const ny = y - dy;
      arr.push({ x: nx, y: ny, F });
      if (Math.hypot(nx - x, ny - y) < tol) break;
      x = nx;
      y = ny;
    }
    return arr;
  }, [f, g, x0, y0, maxIter, tol]);

  return (
    <motion.div variants={cardMotion} initial="initial" animate="enter" whileHover="hover">
      <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700 overflow-hidden`}>
        <CardHeader>
          <CardTitle className="flex items-center gap-2 text-cyan-300">
            <Layers className="w-5 h-5" /> 6.6 Systems of Nonlinear Equations (2D demo)
          </CardTitle>
        </CardHeader>

        <CardContent>
          <div className="grid grid-cols-1 md:grid-cols-4 gap-3 mb-3">
            <div className="md:col-span-2">
              <label className="text-zinc-200 text-xs">f(x,y)</label>
              <Input value={fxExpr} onChange={(e) => setFxExpr(e.target.value)} className="bg-zinc-700 text-white" />
              <div className="text-zinc-400 text-xs mt-1">Example: x*x + y*y - 4</div>
            </div>

            <div>
              <label className="text-zinc-200 text-xs">g(x,y)</label>
              <Input value={gyExpr} onChange={(e) => setGyExpr(e.target.value)} className="bg-zinc-700 text-white" />
            </div>

            <div>
              <label className="text-zinc-200 text-xs">Initial (x₀,y₀)</label>
              <div className="flex gap-2">
                <Input value={String(x0)} onChange={(e) => setX0(Number(e.target.value))} className="bg-zinc-700 text-white" />
                <Input value={String(y0)} onChange={(e) => setY0(Number(e.target.value))} className="bg-zinc-700 text-white" />
              </div>
            </div>
          </div>

          <div className="mb-3">
            <div className="w-full overflow-hidden rounded" style={{ minHeight: 360 }}>
              <Plot
                data={[
                  history.length > 0 ? { x: history.map((h) => h.x), y: history.map((h) => h.y), mode: "markers+lines", name: "Newton path", marker: { size: 6 } } : null,
                ].filter(Boolean)}
                layout={{
                  title: `Newton path (x,y)`,
                  paper_bgcolor: theme.panelBg,
                  plot_bgcolor: theme.panelBg,
                  font: { color: theme.text },
                  xaxis: { title: "x" },
                  yaxis: { title: "y" },
                  height: 360,
                }}
                useResizeHandler={true}
                style={{ width: "100%" }}
              />
            </div>
          </div>

          <div className="grid grid-cols-1 sm:grid-cols-3 gap-3">
            <SummaryCard title="Iterations" value={`${history.length}`} />
            <SummaryCard
              title="Last (x,y)"
              value={`${history.length ? `${history[history.length - 1].x.toFixed(6)}, ${history[history.length - 1].y.toFixed(6)}` : "—"}`}
            />
            <SummaryCard
              title="Residual"
              value={`${history.length && history[history.length - 1].F ? history[history.length - 1].F.map((v) => Math.abs(v).toExponential(2)).join(", ") : "—"}`}
            />
          </div>

          <div className="mt-3 text-zinc-400 text-sm">
            Newton for systems requires solving a linear system at every step. For larger systems, use robust linear algebra and consider quasi-Newton or Jacobian-free methods.
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ---------------- Documentation (LaTeX-ready strings) ----------------
const docsText = {
  "6.1 Simple Fixed-Point Iteration": [
    `**Definition:** A fixed-point iteration solves \\(f(x)=0\\) by rewriting into \\(x=g(x)\\) and applying successive substitution.`,
    `**Algorithm:**\n\n$$x_{k+1} = g(x_k).$$`,
    `**Convergence Condition:** If \\(g\\) is continuous on \\([a,b]\\) and \\(|g'(r)|<1\\) at the root \\(r\\), then the sequence converges.`,
    `**Error Analysis:**\n\n$$e_{k+1} = g'(r) e_k + O(e_k^2).$$\n\nThus convergence is *linear* if \\(|g'(r)|<1\\).`,
    `**Acceleration:**\n- Aitken's \\(\\Delta^2\\) process.\n- Steffensen’s method for quadratic convergence without explicit derivatives.`,
    `**Practical Tip:** Choose rearrangement carefully; avoid forms where \\(|g'(x)|>1\\).`,
  ],

  "6.2 The Newton-Raphson Method": [
    `**Definition:** Iterative method using tangent line approximation.`,
    `**Update Rule:**\n\n$$x_{k+1} = x_k - \\frac{f(x_k)}{f'(x_k)}.$$`,
    `**Geometric Meaning:** Root is found by intersecting tangent line at \\(x_k\\) with x-axis.`,
    `**Convergence:** Quadratic if root is simple. Error relation:\n\n$$e_{k+1} \\approx C e_k^2.$$`,
    `**Drawbacks:** Requires derivative. May fail if \\(f'(x_k)\\approx 0\\). Sensitive to initial guess.`,
    `**Variants:** Damped Newton:\n\n$$x_{k+1}=x_k-\\lambda\\frac{f(x_k)}{f'(x_k)}, \\quad 0<\\lambda<1.$$`,
    `**Example:** Solving \\(x^2-2=0\\). Newton converges to \\(\\sqrt{2}\\).`,
  ],

  "6.3 The Secant Method": [
    `**Definition:** Approximates Newton’s derivative by finite difference between two points.`,
    `**Update Formula:**\n\n$$
x_{k+1} = x_k - f(x_k) \\cdot \\frac{x_k-x_{k-1}}{f(x_k)-f(x_{k-1})}.
$$`,
    `**Convergence:** Superlinear with order \\(\\varphi=(1+\\sqrt{5})/2 \\approx 1.618\\).`,
    `**Advantages:** Does not require derivative. Efficient when f’ expensive.`,
    `**Disadvantages:** Needs two starting points. Can fail if denominator is small.`,
    `**Error Relation:**\n\n$$
e_{k+1} \\approx C e_k^{\\varphi}.
$$`,
    `**Practical Tip:** Combine with bisection for safety.`,
  ],

  "6.4 Brent's Method": [
    `**Definition:** A hybrid root-finding algorithm combining:\n- Bisection (guaranteed convergence)\n- Secant (faster)\n- Inverse quadratic interpolation (very fast)`,
    `**Interpolation Step:**\n\n$$
x_{k+1} =
  \\frac{x_{k-2} f_{k-1} f_k}{(f_{k-2}-f_{k-1})(f_{k-2}-f_k)} +
  \\frac{x_{k-1} f_{k-2} f_k}{(f_{k-1}-f_{k-2})(f_{k-1}-f_k)} +
  \\frac{x_k f_{k-2} f_{k-1}}{(f_k-f_{k-2})(f_k-f_{k-1})}.
$$`,
    `**Advantages:** Robust like bisection, faster like Newton/secant.`,
    `**Convergence:** Superlinear. Widely implemented (e.g., SciPy, MATLAB).`,
  ],

  "6.5 Multiple Roots": [
    `**Definition:** A root \\(r\\) with multiplicity \\(m>1\\).`,
    `**Condition:**\n\n$$
f(r)=0, \\quad f'(r)=0, \\quad …, \\quad f^{(m-1)}(r)=0, \\quad f^{(m)}(r)\\neq 0.
$$`,
    `**Problem:** Standard Newton converges only linearly for multiple roots.`,
    `**Modified Newton:**\n\n$$
x_{k+1}=x_k - m\\cdot \\frac{f(x_k)}{f'(x_k)}.
$$`,
    `**Alternative:** Deflation (factor out \\((x-r)\\)) or higher derivative methods.`,
    `**Example:** For \\(f(x)=(x-1)^3\\), standard Newton is slow, modified Newton restores quadratic convergence.`,
  ],

  "6.6 Systems of Nonlinear Equations": [
    `**Definition:** Solve \\(F(x)=0\\), where \\(F: \\mathbb{R}^n \\to \\mathbb{R}^n\\).`,
    `**Newton’s Method for Systems:**\n\n$$
\\mathbf{x}_{k+1} = \\mathbf{x}_k - J(\\mathbf{x}_k)^{-1} F(\\mathbf{x}_k).
$$`,
    `**Linear Solve per Iteration:**\n\n$$
J(\\mathbf{x}_k) \\Delta x = -F(\\mathbf{x}_k).
$$`,
    `**Convergence:** Quadratic near a simple root if Jacobian is nonsingular.`,
    `**Challenges:** Computing \\(J\\) may be expensive. Sensitive to initial guess.`,
    `**Variants:**\n- Approximate Jacobian by finite differences.\n- Quasi-Newton (Broyden).\n- Damped Newton with line search.`,
    `**Example:**\n\n$$
f_1(x,y)=x^2+y^2-1=0, \\\\
f_2(x,y)=x^2-y=0.
$$`,
  ],
};


 function ChapterDocs() {
  const [open, setOpen] = useState("6.1 Simple Fixed-Point Iteration");

  return (
    <div className="w-full mt-6 bg-zinc-900 border border-zinc-700 rounded-xl shadow-lg overflow-hidden">
      {/* Header */}
      <div className="p-3 border-b border-zinc-700 flex items-center justify-between">
        <div className="flex items-center gap-3">
          <BookOpen className="w-5 h-5 text-cyan-400" />
          <div>
            <div className="text-zinc-100 font-semibold">Chapter 6 — Documentation</div>
            <div className="text-zinc-400 text-xs">
              Theory, practical notes, pseudocode, and usage tips for open methods
            </div>
          </div>
        </div>
        <div className="text-zinc-400 text-sm">
          Interactive examples above — experiment and read the notes
        </div>
      </div>

      {/* Grid layout */}
      <div className="grid grid-cols-1 md:grid-cols-4">
        {/* Sidebar */}
        <div className="col-span-1 border-r border-zinc-800">
          {Object.keys(docsText).map((k) => (
            <button
              key={k}
              onClick={() => setOpen(k)}
              className={`w-full p-3 text-left hover:bg-zinc-800/40 flex items-center justify-between ${
                open === k ? "bg-zinc-800/20" : ""
              }`}
            >
              <div className="flex items-center gap-2">
                <List className="w-4 h-4 text-zinc-300" />
                <div className="text-zinc-100">{k}</div>
              </div>
              {open === k ? (
                <ChevronDown className="w-4 h-4 text-zinc-400" />
              ) : (
                <ChevronRight className="w-4 h-4 text-zinc-400" />
              )}
            </button>
          ))}
        </div>

        {/* Content */}
        <div className="col-span-3 p-4 overflow-auto" style={{ maxHeight: 420 }}>
          <motion.h3
            initial={{ opacity: 0, y: -6 }}
            animate={{ opacity: 1, y: 0 }}
            className="text-xl font-semibold text-zinc-100 mb-3"
          >
            {open}
          </motion.h3>
          <div className="text-zinc-300 space-y-3 prose prose-invert max-w-none">
            {docsText[open].map((para, i) => (
              <div key={i} className="text-sm leading-relaxed">
                <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>
                  {para}
                </ReactMarkdown>
              </div>
            ))}
            <div className="mt-3 text-zinc-400 text-sm">
              References: Burden & Faires, Atkinson, Press et al. *Numerical Recipes*, and modern numerical libraries.
            </div>
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
        <CardTitle className="text-emerald-300 flex items-center gap-2">
          <List className="w-5 h-5" /> Problems & Exercises
        </CardTitle>
      </CardHeader>
      <CardContent>
        <ul className="list-decimal list-inside text-zinc-300 text-sm space-y-2">
          <li>Use fixed-point iteration to solve \\(x=\\cos x\\). Find a suitable g(x) and show iterations.</li>
          <li>Apply Newton–Raphson to find a root of \\(x^3 - 2x - 5\\). Verify quadratic convergence numerically.</li>
          <li>Compare secant vs Newton on \\(\\arctan(x)-0.5\\) to reach \\(|f(x)| &lt; 10^{-8}\\).</li>
          <li>Use Brent (or a library brentq) to solve \\(\\sin(x)-0.5\\) in \\([0,\\pi]\\) and compare speed.</li>
          <li>Demonstrate modified Newton on \\((x-1)^3\\) and compare with standard Newton.</li>
          <li>Solve a 2x2 nonlinear system with Newton and discuss Jacobian conditioning.</li>
        </ul>
      </CardContent>
    </Card>
  );
}

// ---------------- Main page ----------------
export default function Chapter6() {
  return (
    <div className={`min-h-screen px-4 md:px-6 py-6 bg-zinc-950 text-zinc-100`}>
      <motion.header initial={{ opacity: 0, y: -8 }} animate={{ opacity: 1, y: 0 }} className="max-w-7xl mx-auto mb-6">
        <div className="flex flex-col md:flex-row items-start justify-between gap-4">
          <div className="max-w-3xl">
            <h1 className="text-2xl md:text-3xl text-emerald-400 font-bold">Open Methods — Fixed-Point, Newton, Secant, Brent & Systems</h1>
            <p className="text-zinc-300 mt-2">
              Interactive visualizations exploring convergence behaviour of open root-finding methods. Use initial guesses, tolerances and compare
              methods side-by-side.
            </p>
          </div>
        </div>
      </motion.header>

      <main className="max-w-7xl mx-auto grid grid-cols-1 gap-6">
        {/* Top overview cards */}
        <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
          <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700`}>
            <CardHeader>
              <CardTitle className="text-cyan-300 flex items-center gap-2">
                <BookOpen className="w-5 h-5" /> Overview
              </CardTitle>
            </CardHeader>
            <CardContent>
              <div className="text-zinc-300 text-sm">Open methods use local information (derivatives or secants) to converge rapidly — they require good initial guesses and safety fallbacks.</div>
            </CardContent>
          </Card>

          <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700`}>
            <CardHeader>
              <CardTitle className="text-emerald-300 flex items-center gap-2">
                <Layers className="w-5 h-5" /> Visual intuition
              </CardTitle>
            </CardHeader>
            <CardContent>
              <div className="text-zinc-300 text-sm">Tangent and secant visuals clarify Newton and Secant mechanics. Compare iteration residuals and paths for insight.</div>
            </CardContent>
          </Card>

          <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700`}>
            <CardHeader>
              <CardTitle className="text-amber-300 flex items-center gap-2">
                <Zap className="w-5 h-5" /> Practical notes
              </CardTitle>
            </CardHeader>
            <CardContent>
              <ul className="text-zinc-300 text-sm list-disc list-inside">
                <li>Use robust fallbacks (bisection / Brent) when derivatives or guesses are poor.</li>
                <li>Guard against small derivatives and ill-conditioned Jacobians (use pivoting solvers).</li>
                <li>For systems, prefer quasi-Newton for large problems or use sparse solvers.</li>
              </ul>
            </CardContent>
          </Card>
        </div>

        {/* Main panels */}
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          <div className="space-y-6">
            <FixedPointPanel />
            <NewtonPanel />
            <SecantPanel />
          </div>

          <div className="space-y-6">
            <BrentsPanel />
            <MultipleRootsPanel />
            <SystemsPanel />
          </div>
        </div>

        <ChapterDocs />
        <ProblemsCard />

        <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700`}>
          <CardHeader>
            <CardTitle className="text-cyan-300 flex items-center gap-2">
              <Layers className="w-5 h-5" /> Visualization — particle hint
            </CardTitle>
          </CardHeader>
          <CardContent>
            <div className="w-full h-48">
              <Canvas camera={{ position: [0, 0, 4] }}>
                <ambientLight intensity={0.6} />
                <pointLight position={[10, 10, 10]} />
                <ParticleCloud points={64} size={0.05} />
                <OrbitControls />
              </Canvas>
            </div>
            <div className="text-zinc-400 text-sm mt-3">Subtle 3D widget to make the page feel alive — purely decorative.</div>
          </CardContent>
        </Card>
        <BottomBar/>
      </main>
    </div>
  );
}
