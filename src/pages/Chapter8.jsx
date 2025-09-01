// src/pages/Chapter8.jsx
// ======================================================================
// CHAPTER 8 — Case Studies: Roots of Equations
// Four engineering case studies (8.1 → 8.4). Structure and style follow Chapter7.jsx:
// - Dark theme, Plotly visualizations (2D/3D where appropriate)
// - KaTeX via react-markdown + remark-math + rehype-katex
// - Framer-motion micro-animations
// - Exercises & Problems at the end
// - Decorative particle cloud (three-fiber) for continuity with Chapter 7
// ======================================================================

import React, { useMemo, useState, useRef } from "react";
import Plot from "react-plotly.js";
import { motion } from "framer-motion";
import ReactMarkdown from "react-markdown";
import remarkMath from "remark-math";
import rehypeKatex from "rehype-katex";
import "katex/dist/katex.min.css";

import { Card, CardHeader, CardContent, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Button } from "@/components/ui/button";

import {
  BookOpen,
  List,
  ChevronRight,
  ChevronDown,
  FlaskConical,
  CloudSun,
  CircuitBoard,
  Gauge,
  Beaker,
  LineChart,
  Wand2,
  Atom,
  Target,
  Calculator,
  ThermometerSun,
} from "lucide-react";
import BottomBar from "../components/BottomBar";
import { Canvas, useFrame } from "@react-three/fiber";
import { OrbitControls } from "@react-three/drei";

// ======================================================================
// Theme & motion presets (keeps parity with Chapter7's look)
// ======================================================================
const theme = {
  bg: "bg-zinc-950",
  panel: "bg-zinc-900/60 border border-zinc-700",
  text: "text-zinc-200",
  subtext: "text-zinc-400",
  accent: "#60a5fa", // cyan-ish for chapter 8
  accent2: "#34d399", // emerald
  warn: "#f59e0b",
};

const fadeUp = {
  initial: { opacity: 0, y: 12 },
  animate: { opacity: 1, y: 0 },
  transition: { duration: 0.32 },
};

// ======================================================================
// Utility functions used across the chapter
// ======================================================================

/** linspace helper */
function linspace(a, b, n) {
  if (n <= 1) return [a];
  const dx = (b - a) / (n - 1);
  const arr = new Array(n);
  for (let i = 0; i < n; i++) arr[i] = a + i * dx;
  return arr;
}

/** safe number conversion */
function safeN(n, fallback = 0) {
  const v = Number(n);
  return Number.isFinite(v) ? v : fallback;
}

/** sample function for plotting */
function sampleFunction(f, lo, hi, steps = 400) {
  const xs = [];
  const ys = [];
  if (!Number.isFinite(lo) || !Number.isFinite(hi) || lo >= hi) return { xs, ys };
  const dx = (hi - lo) / steps;
  for (let i = 0; i <= steps; i++) {
    const x = lo + i * dx;
    xs.push(x);
    try {
      ys.push(f ? f(x) : NaN);
    } catch {
      ys.push(NaN);
    }
  }
  return { xs, ys };
}

/** Numeric central derivative */
function dfdx(f, x, h = 1e-6) {
  try {
    return (f(x + h) - f(x - h)) / (2 * h);
  } catch {
    return NaN;
  }
}

/** Bisection method returning history for plotting */
function bisection(f, a, b, tol = 1e-8, maxIter = 80) {
  const hist = [];
  try {
    let fa = f(a);
    let fb = f(b);
    if (!Number.isFinite(fa) || !Number.isFinite(fb) || fa * fb > 0) return hist;
    let left = a;
    let right = b;
    for (let i = 0; i < maxIter; i++) {
      const mid = 0.5 * (left + right);
      const fm = f(mid);
      hist.push({ x: mid, fx: fm, a: left, b: right });
      if (!Number.isFinite(fm)) break;
      if (Math.abs(fm) < tol || Math.abs(right - left) < tol) break;
      if (fa * fm < 0) {
        right = mid;
        fb = fm;
      } else {
        left = mid;
        fa = fm;
      }
    }
  } catch {
    // ignore
  }
  return hist;
}

/** Newton's method returning history for plotting */
function newton(f, x0, tol = 1e-10, maxIter = 60) {
  const hist = [];
  try {
    let x = x0;
    for (let i = 0; i < maxIter; i++) {
      const fx = f(x);
      const dfx = dfdx(f, x);
      if (!Number.isFinite(fx) || !Number.isFinite(dfx) || Math.abs(dfx) < 1e-14) break;
      const x1 = x - fx / dfx;
      hist.push({ x: x1, fx: f(x1), df: dfx });
      if (Math.abs(x1 - x) < tol) break;
      x = x1;
    }
  } catch {
    // ignore
  }
  return hist;
}

/** Secant method returning history */
function secant(f, x0, x1, tol = 1e-10, maxIter = 80) {
  const hist = [];
  try {
    let a = x0;
    let b = x1;
    for (let i = 0; i < maxIter; i++) {
      const fa = f(a);
      const fb = f(b);
      const denom = fb - fa;
      if (!Number.isFinite(fa) || !Number.isFinite(fb) || Math.abs(denom) < 1e-14) break;
      const c = b - fb * (b - a) / denom;
      const fc = f(c);
      hist.push({ x: c, fx: fc });
      if (Math.abs(c - b) < tol) break;
      a = b;
      b = c;
    }
  } catch {
    // ignore
  }
  return hist;
}

// ======================================================================
// Small UI helpers (Equation, Note, Labeled, Summary, Exercise)
// Similar to Chapter7's helpers for consistent experience
// ======================================================================
function Equation({ children }) {
  return (
    <div className="my-2">
      <div className="px-3 py-2 rounded-lg overflow-auto border text-white border-zinc-700 bg-zinc-900/80">
        <div className="text-xs uppercase text-white tracking-widest mb-1" style={{ color: theme.accent }}>
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
        <div key={i} className="rounded-lg bg-zinc-900/70 border border-zinc-700 px-3 py-2 flex flex-col">
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
        <Button variant="ghost" className="bg-white cursor-pointer" size="sm" onClick={() => setOpen((o) => !o)}>
          {open ? <ChevronDown className="w-4 h-4" /> : <ChevronRight className="w-4 h-4" />}
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
// Decorative particle cloud component — reused across chapters
// ======================================================================
function ParticleHalo({ points = 72 }) {
  const ref = useRef();
  useFrame(({ clock }) => {
    if (!ref.current) return;
    const t = clock.getElapsedTime();
    ref.current.rotation.y = 0.07 * t;
    ref.current.rotation.x = 0.03 * Math.sin(t * 0.6);
  });

  const verts = useMemo(() => {
    const arr = [];
    for (let i = 0; i < points; i++) {
      const phi = (i / points) * Math.PI * 2;
      const r = 1.05 + 0.25 * Math.sin(3 * phi);
      arr.push(r * Math.cos(phi), r * Math.sin(phi), 0.14 * Math.cos(2 * phi));
    }
    return new Float32Array(arr);
  }, [points]);

  return (
    <points ref={ref}>
      <bufferGeometry>
        <bufferAttribute attach="attributes-position" array={verts} count={verts.length / 3} itemSize={3} />
      </bufferGeometry>
      <pointsMaterial size={0.04} sizeAttenuation color={theme.accent2} />
    </points>
  );
}

// ======================================================================
// Docs text for Chapter 8 (rendered in docs panel)
// ======================================================================
const docsText = {
  "8.1 Ideal & Nonideal Gas Laws (Chemical/Bio Eng.)": [
    `**Goal:** Solve for molar volume $V_m$ from Van der Waals equation or compare with ideal gas.`,
    `**Van der Waals (per mole):**  
$$
\\left(P + \\frac{a}{V_m^2}\\right) (V_m - b) = R T
$$  
Set $f(V_m) = (P + a/V_m^2)(V_m - b) - RT$ and solve $f(V_m)=0$.`,
    `**Comparison with Ideal Gas:**  
$$
PV_m = RT
$$  
For large $V_m$, corrections vanish and the Van der Waals law reduces to the ideal gas law.`,
    `**Critical Constants:**  
From Van der Waals EOS, the critical temperature, pressure, and volume are:  
$$
T_c = \\frac{8a}{27Rb}, \\quad P_c = \\frac{a}{27b^2}, \\quad V_c = 3b
$$  
Useful for scaling real gas data.`,
    `**Practical:** At high pressure the non-ideal corrections (a,b) matter; multiple real roots can appear near phase equilibrium. Use bracketing if unsure.`
  ],

  "8.2 Greenhouse Gases & Rainwater (Civil/Environmental Eng.)": [
    `**Goal:** Solve for planetary equilibrium temperature from radiative balance (0D model).`,
    `**Balance Equation:**  
$$
\\sigma T^4 = \\frac{S(1-\\alpha)}{4} + \\Delta F
$$  
where $S$ is solar constant, $\\alpha$ albedo, and $\\Delta F$ additional forcing (W/m²).`,
    `**Residual Function:**  
$$
f(T) = \\sigma T^4 - Q_{in} = 0
$$  
with $Q_{in} = \\frac{S(1-\\alpha)}{4} + \\Delta F$.`,
    `**Climate Sensitivity:**  
Approximate change in temperature from forcing:  
$$
\\Delta T \\approx \\lambda \\Delta F
$$  
with $\\lambda$ climate sensitivity parameter (K per W/m²).`,
    `**Practical:** Use Newton for rapid convergence; bracket for safety if Newton misbehaves. Results help estimate warming from greenhouse gases or cooling from aerosols.`
  ],

  "8.3 Design of an Electric Circuit (Electrical Eng.)": [
    `**Goal:** Solve Kirchhoff + Shockley diode: nonlinear algebraic equation in current $I$.`,
    `**Model:** For series resistor $R$ and diode with saturation current $I_s$, ideality $n$, and thermal voltage $V_t$:  
$$
I = I_s\\left( e^{\\frac{V_s - IR}{nV_t}} - 1 \\right)
$$  
Residual:  
$$
f(I) = I_s\\Big(e^{(V_s - IR)/(nV_t)} - 1\\Big) - I
$$`,
    `**Small-Signal Resistance:**  
Dynamic diode resistance is:  
$$
r_d = \\frac{nV_t}{I + I_s}
$$  
Important for AC analysis.`,
    `**Power Dissipation:**  
$$
P = IV
$$  
Useful for checking diode thermal limits.`,
    `**Practical:** Exponential can overflow; clamp exponent or use bracketing to ensure numerical stability. For low-voltage circuits, Shockley model may require series resistance corrections.`
  ],

  "8.4 Pipe Friction (Mechanical/Aerospace Eng.)": [
    `**Goal:** Solve Colebrook equation for Darcy friction factor $f$ given Re and relative roughness $\\varepsilon/D$.`,
    `**Colebrook implicit form:**  
$$
\\frac{1}{\\sqrt{f}} = -2\\log_{10}\\left( \\frac{\\varepsilon}{3.7D} + \\frac{2.51}{Re\\sqrt{f}} \\right)
$$  
Residual:  
$$
g(f) = \\frac{1}{\\sqrt{f}} + 2\\log_{10}\\Big(\\tfrac{\\varepsilon}{3.7D} + \\tfrac{2.51}{Re\\sqrt{f}}\\Big)
$$`,
    `**Approximate Explicit Formula (Swamee–Jain):**  
$$
f \\approx \\frac{0.25}{\\Big[\\log_{10}\\big( \\tfrac{\\varepsilon}{3.7D} + \\tfrac{5.74}{Re^{0.9}} \\big)\\Big]^2}
$$  
Good initial guess for iteration.`,
    `**Moody Chart:** Graphical representation of friction factor vs. Reynolds number and roughness, widely used in engineering design.`,
    `**Practical:** Swamee–Jain gives explicit approximation; use it as an initial guess for faster convergence when solving Colebrook iteratively.`
  ],
};


// ======================================================================
// 8.1 Panel — Van der Waals (interactive, plots, comparisons)
// ======================================================================
function Section81() {
  const [P, setP] = useState(20); // bar
  const [T, setT] = useState(300); // K
  const [a, setA] = useState(1.352); // (bar·L^2)/mol^2 (example)
  const [b, setB] = useState(0.0387); // L/mol
  const [Rgas, setRgas] = useState(0.08314); // L·bar/(mol·K)
  const [method, setMethod] = useState("newton");
  const [x0, setX0] = useState(1.2);
  const [x1, setX1] = useState(5.0);
  const [aBis, setABis] = useState(0.2);
  const [bBis, setBBis] = useState(10.0);

  // residual f(V)
  const f = useMemo(() => {
    return (V) => {
      if (!Number.isFinite(V) || V <= b) return NaN;
      return (P + a / (V * V)) * (V - b) - Rgas * T;
    };
  }, [P, a, b, Rgas, T]);

  const Videal = useMemo(() => {
    if (P <= 0) return NaN;
    return (Rgas * T) / P;
  }, [Rgas, T, P]);

  const history = useMemo(() => {
    if (method === "newton") return newton(f, Number(x0), 1e-10, 80);
    if (method === "secant") return secant(f, Number(x0), Number(x1), 1e-10, 80);
    return bisection(f, Number(aBis), Number(bBis), 1e-10, 200);
  }, [f, method, x0, x1, aBis, bBis]);

  const last = history.length ? history[history.length - 1] : null;

  const range = useMemo(() => {
    const left = Math.max(b + 1e-6, Math.min(aBis || b + 0.1, bBis || (Videal || 1) * 0.2));
    const right = Math.max(b + 0.5, Math.max(bBis || 2, (Videal || 5) * 3));
    return [left, right];
  }, [b, aBis, bBis, Videal]);

  const sampled = useMemo(() => sampleFunction(f, range[0], range[1], 600), [f, range]);

  return (
    <motion.div {...fadeUp}>
      <Card className={theme.panel}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <FlaskConical className="w-5 h-5" />
            8.1 Ideal & Nonideal Gas Laws — Van der Waals (Case Study)
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-3">
          <div className="grid md:grid-cols-3 gap-3">
            <div>
              <Labeled label="Pressure P (bar)">
                <Input value={String(P)} onChange={(e) => setP(safeN(e.target.value, P))} className="bg-zinc-800 text-white" />
              </Labeled>
              <Labeled label="Temperature T (K)">
                <Input value={String(T)} onChange={(e) => setT(safeN(e.target.value, T))} className="bg-zinc-800 text-white" />
              </Labeled>
              <div className="flex gap-2 mt-2">
                <Labeled label="a (vdW)">
                  <Input value={String(a)} onChange={(e) => setA(safeN(e.target.value, a))} className="bg-zinc-800 text-white" />
                </Labeled>
                <Labeled label="b (vdW)">
                  <Input value={String(b)} onChange={(e) => setB(safeN(e.target.value, b))} className="bg-zinc-800 text-white" />
                </Labeled>
              </div>
              <div className="text-zinc-400 text-xs mt-2">V<sub>ideal</sub> = {Number.isFinite(Videal) ? Videal.toFixed(4) : "—"} L/mol</div>
            </div>

            <div>
              <Labeled label="Method">
                <div className="flex gap-2">
                  <Button className="bg-emerald-400 cursor-pointer hover:bg-emerald-500" size="sm" variant={method === "newton" ? "default" : "ghost"} onClick={() => setMethod("newton")}>
                    Newton
                  </Button>
                  <Button className="bg-emerald-400 cursor-pointer hover:bg-emerald-500" size="sm" variant={method === "secant" ? "default" : "ghost"} onClick={() => setMethod("secant")}>
                    Secant
                  </Button>
                  <Button className="bg-emerald-400 cursor-pointer hover:bg-emerald-500" size="sm" variant={method === "bisection" ? "default" : "ghost"} onClick={() => setMethod("bisection")}>
                    Bisection
                  </Button>
                </div>
              </Labeled>

              {method === "newton" && (
                <Labeled label="Initial guess V₀ (L/mol)">
                  <Input value={String(x0)} onChange={(e) => setX0(safeN(e.target.value, x0))} className="bg-zinc-800 text-white" />
                </Labeled>
              )}

              {method === "secant" && (
                <div className="space-y-2">
                  <Labeled label="Initial guesses V₀, V₁">
                    <div className="flex gap-2">
                      <Input value={String(x0)} onChange={(e) => setX0(safeN(e.target.value, x0))} className="bg-zinc-800 text-white" />
                      <Input value={String(x1)} onChange={(e) => setX1(safeN(e.target.value, x1))} className="bg-zinc-800 text-white" />
                    </div>
                  </Labeled>
                </div>
              )}

              {method === "bisection" && (
                <div className="space-y-2">
                  <Labeled label="Bracket [a, b] (L/mol)">
                    <div className="flex gap-2">
                      <Input value={String(aBis)} onChange={(e) => setABis(safeN(e.target.value, aBis))} className="bg-zinc-800 text-white" />
                      <Input value={String(bBis)} onChange={(e) => setBBis(safeN(e.target.value, bBis))} className="bg-zinc-800 text-white" />
                    </div>
                  </Labeled>
                </div>
              )}
            </div>

            <div>
              <Summary
                items={[
                  { label: "Iterations", value: history.length },
                  { label: "Root V (L/mol)", value: last ? last.x.toPrecision(6) : "—" },
                  { label: "|f(V)|", value: last ? Math.abs(last.fx).toExponential(3) : "—" },
                ]}
              />
              <div className="text-zinc-400 text-xs mt-2">Tip: Avoid brackets that include V ≤ b (singularity).</div>
            </div>
          </div>

          <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2">
            <Plot
              data={[
                { x: sampled.xs, y: sampled.ys, type: "scatter", mode: "lines", line: { color: theme.accent }, name: "f(V)" },
                history.length > 0 ? { x: history.map((h) => h.x), y: history.map((h) => h.fx), mode: "markers+lines", name: "iterates", marker: { size: 7 } } : null,
              ].filter(Boolean)}
              layout={{
                title: "Van der Waals residual vs V",
                paper_bgcolor: "rgba(0,0,0,0)",
                plot_bgcolor: "rgba(0,0,0,0)",
                font: { color: "#e5e7eb" },
                xaxis: { title: "V (L/mol)" },
                yaxis: { title: "f(V)" },
                margin: { t: 40, l: 50, r: 20, b: 50 },
                height: 360,
              }}
              style={{ width: "100%" }}
              useResizeHandler
            />
          </div>

          <div className="grid md:grid-cols-2 gap-3">
            <div>
              <Equation>
               {`Van der Waals: $\\left(P + \\dfrac{a}{V^2}\\right)(V - b) - RT = 0$`}
              </Equation>
            </div>
            <div>
              <Note>
                Compare the computed root with the ideal gas estimation. At high P the non-ideal root deviates significantly from ideal.
              </Note>
            </div>
          </div>

          <div className="space-y-2">
            <Exercise
              title="Exercise 8.1.1 — Compare with ideal gas"
              body={`Take P = 50 bar, T = 320 K with a = 2.283, b = 0.0428 (methane-like). Solve for molar volume and compare with \\(V_{ideal}=RT/P\\).`}
              solution={`Compute \\(V_{ideal} = RT/P\\), then solve Van der Waals numerically. Discuss differences; non-ideal V is often lower at high P due to attraction term a/V².`}
            />
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// 8.2 Panel — Radiative balance (GHG forcing) (interactive T root)
// ======================================================================
function Section82() {
  const [S, setS] = useState(1361); // solar constant
  const [albedo, setAlbedo] = useState(0.30);
  const [deltaF, setDeltaF] = useState(2.0);
  const sigma = 5.670374419e-8;
  const [method, setMethod] = useState("newton");
  const [T0, setT0] = useState(255);
  const [T1, setT1] = useState(300);
  const [aBis, setABis] = useState(150);
  const [bBis, setBBis] = useState(400);

  const Qin = useMemo(() => (S * (1 - albedo)) / 4 + deltaF, [S, albedo, deltaF]);

  const f = useMemo(() => {
    return (T) => {
      if (!Number.isFinite(T) || T <= 0) return NaN;
      return sigma * T ** 4 - Qin;
    };
  }, [Qin]);

  const history = useMemo(() => {
    if (method === "newton") return newton(f, Number(T0), 1e-10, 80);
    if (method === "secant") return secant(f, Number(T0), Number(T1), 1e-10, 80);
    return bisection(f, Number(aBis), Number(bBis), 1e-10, 200);
  }, [f, method, T0, T1, aBis, bBis]);

  const last = history.length ? history[history.length - 1] : null;
  const sampled = useMemo(() => sampleFunction(f, 120, 420, 600), [f]);

  return (
    <motion.div {...fadeUp}>
      <Card className={theme.panel}>
        <CardHeader>
          <CardTitle className="text-emerald-300 flex items-center gap-2">
            <CloudSun className="w-5 h-5" />
            8.2 Greenhouse Gases & Planetary Temperature — Radiative Balance
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-3">
          <div className="grid md:grid-cols-3 gap-3">
            <div>
              <Labeled label="Solar constant S (W/m²)">
                <Input value={String(S)} onChange={(e) => setS(safeN(e.target.value, S))} className="bg-zinc-800 text-white" />
              </Labeled>
              <Labeled label="Albedo α">
                <Input value={String(albedo)} onChange={(e) => setAlbedo(safeN(e.target.value, albedo))} className="bg-zinc-800 text-white" />
              </Labeled>
              <Labeled label="Forcing ΔF (W/m²)">
                <Input value={String(deltaF)} onChange={(e) => setDeltaF(safeN(e.target.value, deltaF))} className="bg-zinc-800 text-white" />
              </Labeled>
              <div className="text-zinc-400 text-xs mt-2">Qin = {Qin.toFixed(3)} W/m²</div>
            </div>

            <div>
              <Labeled label="Method">
                <div className="flex gap-2">
                  <Button className="bg-emerald-400 cursor-pointer hover:bg-emerald-500" size="sm" variant={method === "newton" ? "default" : "ghost"} onClick={() => setMethod("newton")}>
                    Newton
                  </Button>
                  <Button className="bg-emerald-400 cursor-pointer hover:bg-emerald-500" size="sm" variant={method === "secant" ? "default" : "ghost"} onClick={() => setMethod("secant")}>
                    Secant
                  </Button>
                  <Button className="bg-emerald-400 cursor-pointer hover:bg-emerald-500" size="sm" variant={method === "bisection" ? "default" : "ghost"} onClick={() => setMethod("bisection")}>
                    Bisection
                  </Button>
                </div>
              </Labeled>

              {method === "newton" && (
                <Labeled label="Initial guess T₀ (K)">
                  <Input value={String(T0)} onChange={(e) => setT0(safeN(e.target.value, T0))} className="bg-zinc-800 text-white" />
                </Labeled>
              )}

              {method === "secant" && (
                <Labeled label="Initial T₀, T₁ (K)">
                  <div className="flex gap-2">
                    <Input value={String(T0)} onChange={(e) => setT0(safeN(e.target.value, T0))} className="bg-zinc-800 text-white" />
                    <Input value={String(T1)} onChange={(e) => setT1(safeN(e.target.value, T1))} className="bg-zinc-800 text-white" />
                  </div>
                </Labeled>
              )}

              {method === "bisection" && (
                <Labeled label="Bracket [a, b] (K)">
                  <div className="flex gap-2">
                    <Input value={String(aBis)} onChange={(e) => setABis(safeN(e.target.value, aBis))} className="bg-zinc-800 text-white" />
                    <Input value={String(bBis)} onChange={(e) => setBBis(safeN(e.target.value, bBis))} className="bg-zinc-800 text-white" />
                  </div>
                </Labeled>
              )}
            </div>

            <div>
              <Summary
                items={[
                  { label: "Iterations", value: history.length },
                  { label: "T (K)", value: last ? last.x.toFixed(3) : "—" },
                  { label: "|f(T)|", value: last ? Math.abs(last.fx).toExponential(3) : "—" },
                ]}
              />
              <div className="text-zinc-400 text-xs mt-2">Newton converges quadratically near root for well-behaved derivatives.</div>
            </div>
          </div>

          <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2">
            <Plot
              data={[
                { x: sampled.xs, y: sampled.ys, type: "scatter", mode: "lines", line: { color: theme.accent2 }, name: "f(T)" },
                history.length > 0 ? { x: history.map((h) => h.x), y: history.map((h) => h.fx), mode: "markers+lines", name: "iterates", marker: { size: 7 } } : null,
              ].filter(Boolean)}
              layout={{
                title: "Radiative balance residual vs T",
                paper_bgcolor: "rgba(0,0,0,0)",
                plot_bgcolor: "rgba(0,0,0,0)",
                font: { color: "#e5e7eb" },
                xaxis: { title: "T (K)" },
                yaxis: { title: "f(T) (W/m²)" },
                margin: { t: 40, l: 50, r: 20, b: 50 },
                height: 360,
              }}
              style={{ width: "100%" }}
              useResizeHandler
            />
          </div>

          <div className="space-y-2">
            <Exercise
              title="Exercise 8.2.1 — Sensitivity to ΔF"
              body={`Vary ΔF from 0 to 4 W/m² and plot equilibrium T. How much does T change per 1 W/m² of forcing?`}
              solution={`Compute roots for a range of ΔF; slope dT/d(ΔF) near current equilibrium is approximately 1/(4σT^3), small but significant over century timescales.`}
            />
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// 8.3 Panel — Electric Circuit (Diode + Resistor) (solve for current I)
// ======================================================================
function Section83() {
  const [Vs, setVs] = useState(5.0);
  const [Rval, setRval] = useState(1000);
  const [Is, setIs] = useState(2e-12);
  const [n, setN] = useState(1.8);
  const [Vt, setVt] = useState(0.02585);
  const [method, setMethod] = useState("newton");
  const [I0, setI0] = useState(0.001);
  const [I1, setI1] = useState(0.01);
  const [aBis, setABis] = useState(0.0);
  const [bBis, setBBis] = useState(0.02);

  const f = useMemo(() => {
    return (I) => {
      if (!Number.isFinite(I)) return NaN;
      const arg = (Vs - I * Rval) / (n * Vt);
      if (!Number.isFinite(arg)) return NaN;
      // prevent overflow
      if (arg > 700) return Number.POSITIVE_INFINITY;
      const expTerm = Math.exp(arg);
      return Is * (expTerm - 1) - I;
    };
  }, [Vs, Rval, Is, n, Vt]);

  const history = useMemo(() => {
    if (method === "newton") return newton(f, Number(I0), 1e-12, 120);
    if (method === "secant") return secant(f, Number(I0), Number(I1), 1e-12, 120);
    return bisection(f, Number(aBis), Number(bBis), 1e-12, 200);
  }, [f, method, I0, I1, aBis, bBis]);

  const last = history.length ? history[history.length - 1] : null;

  const sampled = useMemo(() => sampleFunction(f, 0, Math.max(0.02, Number(bBis) || 0.02), 600), [f, bBis]);

  const Vd = useMemo(() => {
    if (!last || !Number.isFinite(last.x)) return NaN;
    return Vs - last.x * Rval;
  }, [last, Vs, Rval]);

  return (
    <motion.div {...fadeUp}>
      <Card className={theme.panel}>
        <CardHeader>
          <CardTitle className="text-amber-300 flex items-center gap-2">
            <CircuitBoard className="w-5 h-5" />
            8.3 Electric Circuit — Diode + Resistor (Solve for I)
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-3">
          <div className="grid md:grid-cols-3 gap-3">
            <div>
              <Labeled label="Supply Vs (V)">
                <Input value={String(Vs)} onChange={(e) => setVs(safeN(e.target.value, Vs))} className="bg-zinc-800 text-white" />
              </Labeled>
              <Labeled label="Resistor R (Ω)">
                <Input value={String(Rval)} onChange={(e) => setRval(safeN(e.target.value, Rval))} className="bg-zinc-800 text-white" />
              </Labeled>
              <div className="flex gap-2 mt-2">
                <Labeled label="Is (A)">
                  <Input value={String(Is)} onChange={(e) => setIs(safeN(e.target.value, Is))} className="bg-zinc-800 text-white" />
                </Labeled>
                <Labeled label="n">
                  <Input value={String(n)} onChange={(e) => setN(safeN(e.target.value, n))} className="bg-zinc-800 text-white" />
                </Labeled>
                <Labeled label="Vt (V)">
                  <Input value={String(Vt)} onChange={(e) => setVt(safeN(e.target.value, Vt))} className="bg-zinc-800 text-white" />
                </Labeled>
              </div>
            </div>

            <div>
              <Labeled label="Method">
                <div className="flex gap-2">
                  <Button className="bg-emerald-400 cursor-pointer hover:bg-emerald-500" size="sm" variant={method === "newton" ? "default" : "ghost"} onClick={() => setMethod("newton")}>
                    Newton
                  </Button>
                  <Button className="bg-emerald-400 cursor-pointer hover:bg-emerald-500" size="sm" variant={method === "secant" ? "default" : "ghost"} onClick={() => setMethod("secant")}>
                    Secant
                  </Button>
                  <Button className="bg-emerald-400 cursor-pointer hover:bg-emerald-500" size="sm" variant={method === "bisection" ? "default" : "ghost"} onClick={() => setMethod("bisection")}>
                    Bisection
                  </Button>
                </div>
              </Labeled>

              {method === "newton" && (
                <Labeled label="Initial guess I₀ (A)">
                  <Input value={String(I0)} onChange={(e) => setI0(safeN(e.target.value, I0))} className="bg-zinc-800 text-white" />
                </Labeled>
              )}

              {method === "secant" && (
                <Labeled label="Initial I₀, I₁ (A)">
                  <div className="flex gap-2">
                    <Input value={String(I0)} onChange={(e) => setI0(safeN(e.target.value, I0))} className="bg-zinc-800 text-white" />
                    <Input value={String(I1)} onChange={(e) => setI1(safeN(e.target.value, I1))} className="bg-zinc-800 text-white" />
                  </div>
                </Labeled>
              )}

              {method === "bisection" && (
                <Labeled label="Bracket [a, b] (A)">
                  <div className="flex gap-2">
                    <Input value={String(aBis)} onChange={(e) => setABis(safeN(e.target.value, aBis))} className="bg-zinc-800 text-white" />
                    <Input value={String(bBis)} onChange={(e) => setBBis(safeN(e.target.value, bBis))} className="bg-zinc-800 text-white" />
                  </div>
                </Labeled>
              )}
            </div>

            <div>
              <Summary
                items={[
                  { label: "Iterations", value: history.length },
                  { label: "I (A)", value: last ? last.x.toExponential(6) : "—" },
                  { label: "Vd (V)", value: Number.isFinite(Vd) ? Vd.toFixed(4) : "—" },
                ]}
              />
              <div className="text-zinc-400 text-xs mt-2">Tip: Use Swamee–Jain-style explicit estimates for hydraulics; for diode clamp exponent to avoid overflow.</div>
            </div>
          </div>

          <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2">
            <Plot
              data={[
                { x: sampled.xs, y: sampled.ys, type: "scatter", mode: "lines", line: { color: theme.warn }, name: "f(I)" },
                history.length > 0 ? { x: history.map((h) => h.x), y: history.map((h) => h.fx), mode: "markers+lines", name: "iterates", marker: { size: 7 } } : null,
              ].filter(Boolean)}
              layout={{
                title: "Circuit residual vs Current",
                paper_bgcolor: "rgba(0,0,0,0)",
                plot_bgcolor: "rgba(0,0,0,0)",
                font: { color: "#e5e7eb" },
                xaxis: { title: "I (A)" },
                yaxis: { title: "f(I)" },
                margin: { t: 40, l: 50, r: 20, b: 50 },
                height: 360,
              }}
              style={{ width: "100%" }}
              useResizeHandler
            />
          </div>

          <div className="space-y-2">
            <Exercise
              title="Exercise 8.3.1 — Compare methods"
              body={`Take Vs=3.3V, R=220Ω, Is=1e-12 A, n=2. Solve for I using Newton and bisection (start bracket [0,0.02]). Compare iteration counts and robustness.`}
              solution={`Newton is faster if initial guess near the solution; bisection is robust but slower. Exponential stiffness can mislead Newton if derivative near zero.`}
            />
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// 8.4 Panel — Pipe friction (Colebrook) (implicit eqn for Darcy f)
// ======================================================================
function Section84() {
  const [Re, setRe] = useState(250000);
  const [eps, setEps] = useState(0.00015);
  const [D, setD] = useState(0.15);
  const [method, setMethod] = useState("newton");
  const [f0, setF0] = useState(0.02);
  const [f1, setF1] = useState(0.03);
  const [aBis, setABis] = useState(0.005);
  const [bBis, setBBis] = useState(0.1);

  const colebrook = useMemo(() => {
    return (f) => {
      if (!Number.isFinite(f) || f <= 0) return Number.POSITIVE_INFINITY;
      const invSqrt = 1 / Math.sqrt(f);
      const arg = eps / (3.7 * D) + 2.51 / (Re * Math.sqrt(f));
      if (!Number.isFinite(arg) || arg <= 0) return Number.POSITIVE_INFINITY;
      return invSqrt + 2 * Math.log10(arg);
    };
  }, [Re, eps, D]);

  const history = useMemo(() => {
    if (method === "newton") return newton(colebrook, Number(f0), 1e-12, 120);
    if (method === "secant") return secant(colebrook, Number(f0), Number(f1), 1e-12, 120);
    return bisection(colebrook, Number(aBis), Number(bBis), 1e-12, 200);
  }, [colebrook, method, f0, f1, aBis, bBis]);

  const last = history.length ? history[history.length - 1] : null;
  const sampled = useMemo(() => sampleFunction(colebrook, 0.004, 0.12, 700), [colebrook]);

  // Swamee–Jain as reference
  const swamee = useMemo(() => {
    if (Re <= 0) return NaN;
    const A = eps / (3.7 * D);
    const B = 5.74 / Math.pow(Re, 0.9);
    const fSJ = 0.25 / Math.pow(Math.log10(A + B), 2);
    return fSJ;
  }, [Re, eps, D]);

  return (
    <motion.div {...fadeUp}>
      <Card className={theme.panel}>
        <CardHeader>
          <CardTitle className="text-violet-300 flex items-center gap-2">
            <Gauge className="w-5 h-5" />
            8.4 Pipe Friction — Colebrook (Darcy friction factor)
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-3">
          <div className="grid md:grid-cols-3 gap-3">
            <div>
              <Labeled label="Reynolds number Re">
                <Input value={String(Re)} onChange={(e) => setRe(safeN(e.target.value, Re))} className="bg-zinc-800 text-white" />
              </Labeled>
              <Labeled label="Roughness ε (m)">
                <Input value={String(eps)} onChange={(e) => setEps(safeN(e.target.value, eps))} className="bg-zinc-800 text-white" />
              </Labeled>
              <Labeled label="Pipe diameter D (m)">
                <Input value={String(D)} onChange={(e) => setD(safeN(e.target.value, D))} className="bg-zinc-800 text-white" />
              </Labeled>
            </div>

            <div>
              <Labeled label="Method">
                <div className="flex gap-2">
                  <Button className="bg-emerald-400 cursor-pointer hover:bg-emerald-500" size="sm" variant={method === "newton" ? "default" : "ghost"} onClick={() => setMethod("newton")}>
                    Newton
                  </Button>
                  <Button className="bg-emerald-400 cursor-pointer hover:bg-emerald-500" size="sm" variant={method === "secant" ? "default" : "ghost"} onClick={() => setMethod("secant")}>
                    Secant
                  </Button>
                  <Button className="bg-emerald-400 cursor-pointer hover:bg-emerald-500" size="sm" variant={method === "bisection" ? "default" : "ghost"} onClick={() => setMethod("bisection")}>
                    Bisection
                  </Button>
                </div>
              </Labeled>

              {method === "newton" && (
                <Labeled label="Initial guess f₀">
                  <Input value={String(f0)} onChange={(e) => setF0(safeN(e.target.value, f0))} className="bg-zinc-800 text-white" />
                </Labeled>
              )}

              {method === "secant" && (
                <Labeled label="Initial f₀, f₁">
                  <div className="flex gap-2">
                    <Input value={String(f0)} onChange={(e) => setF0(safeN(e.target.value, f0))} className="bg-zinc-800 text-white" />
                    <Input value={String(f1)} onChange={(e) => setF1(safeN(e.target.value, f1))} className="bg-zinc-800 text-white" />
                  </div>
                </Labeled>
              )}

              {method === "bisection" && (
                <Labeled label="Bracket [a, b]">
                  <div className="flex gap-2">
                    <Input value={String(aBis)} onChange={(e) => setABis(safeN(e.target.value, aBis))} className="bg-zinc-800 text-white" />
                    <Input value={String(bBis)} onChange={(e) => setBBis(safeN(e.target.value, bBis))} className="bg-zinc-800 text-white" />
                  </div>
                </Labeled>
              )}
            </div>

            <div>
              <Summary
                items={[
                  { label: "Iterations", value: history.length },
                  { label: "f (Darcy)", value: last ? last.x.toFixed(6) : "—" },
                  { label: "Swamee–Jain", value: Number.isFinite(swamee) ? swamee.toFixed(6) : "—" },
                ]}
              />
              <div className="text-zinc-400 text-xs mt-2">Swamee–Jain provides a close explicit approximation useful as initial guess.</div>
            </div>
          </div>

          <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2">
            <Plot
              data={[
                { x: sampled.xs, y: sampled.ys, type: "scatter", mode: "lines", line: { color: theme.accent2 }, name: "g(f)" },
                history.length > 0 ? { x: history.map((h) => h.x), y: history.map((h) => h.fx), mode: "markers+lines", name: "iterates", marker: { size: 7 } } : null,
              ].filter(Boolean)}
              layout={{
                title: "Colebrook residual vs f (Darcy)",
                paper_bgcolor: "rgba(0,0,0,0)",
                plot_bgcolor: "rgba(0,0,0,0)",
                font: { color: "#e5e7eb" },
                xaxis: { title: "f (dimensionless)" },
                yaxis: { title: "g(f)" },
                margin: { t: 40, l: 50, r: 20, b: 50 },
                height: 360,
              }}
              style={{ width: "100%" }}
              useResizeHandler
            />
          </div>

          <div className="space-y-2">
            <Exercise
              title="Exercise 8.4.1 — Compare with Swamee–Jain"
              body={`For Re=1.2e5, ε=0.0002 m, D=0.1 m, compute f by Colebrook and compare with Swamee–Jain approximation. Report relative error.`}
              solution={`Swamee–Jain gives an explicit initial guess; compute both values and error = |f_colebrook - f_SJ| / f_colebrook.`}
            />
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// Overview cards (top of page) — match Chapter7 style
// ======================================================================
function OverviewCards() {
  return (
    <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
      <Card className={theme.panel}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Beaker className="w-5 h-5" />
            Overview
          </CardTitle>
        </CardHeader>
        <CardContent>
          <div className="text-zinc-300 text-sm">
            Four engineering case studies where root-finding is central. Each study provides interactive parameters, method selection, and plots of residuals with iterates.
          </div>
        </CardContent>
      </Card>

      <Card className={theme.panel}>
        <CardHeader>
          <CardTitle className="text-emerald-300 flex items-center gap-2">
            <LineChart className="w-5 h-5" />
            Visual Intuition
          </CardTitle>
        </CardHeader>
        <CardContent>
          <div className="text-zinc-300 text-sm">
            Visualizing residuals and iterate trajectories helps understand convergence, robustness, and when to prefer bracketing methods.
          </div>
        </CardContent>
      </Card>

      <Card className={theme.panel}>
        <CardHeader>
          <CardTitle className="text-amber-300 flex items-center gap-2">
            <Wand2 className="w-5 h-5" />
            Practical Notes
          </CardTitle>
        </CardHeader>
        <CardContent>
          <div className="text-zinc-300 text-sm">
            Check units, guard against singularities (V→b), clamp exponential arguments, and prefer bracketing for safety.
          </div>
        </CardContent>
      </Card>
    </div>
  );
}

// ======================================================================
// Documentation Panel (left-nav + content) for Chapter 8
// ======================================================================
function ChapterDocs() {
  const [open, setOpen] = useState(Object.keys(docsText)[0]);

  return (
    <div className="w-full mt-6 bg-zinc-900 border border-zinc-700 rounded-xl shadow-lg overflow-hidden">
      <div className="p-3 border-b border-zinc-700 flex items-center justify-between">
        <div className="flex items-center gap-3">
          <BookOpen className="w-5 h-5 text-cyan-400" />
          <div>
            <div className="text-zinc-100 font-semibold">Chapter 8 — Documentation</div>
            <div className="text-zinc-400 text-xs">Theory, derivations, and practical guidance for the case studies</div>
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

        <div className="col-span-3 p-4 overflow-auto" style={{ maxHeight: 420 }}>
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
            <div className="mt-3 text-zinc-400 text-xs">References: standard engineering thermofluids and numerical methods texts.</div>
          </div>
        </div>
      </div>
    </div>
  );
}

// ======================================================================
// Problems section — interactive checkboxes + hints (similar to earlier ask)
// ======================================================================
const defaultProblems = [
  {
    id: "p1",
    text: "For methane at T=320 K, P=50 bar, solve Van der Waals for V and compare with ideal gas. (a=2.283, b=0.0428)",
    hint: "Compute V_ideal = RT/P. Use bracketing around V_ideal and ensure V>b.",
  },
  {
    id: "p2",
    text: "Using radiative balance, compute equilibrium T for S=1361, α=0.28, ΔF ∈ [0,4]. Plot T vs ΔF.",
    hint: "For each ΔF, find root in [200,350] K using bisection or Newton (with good initial guess).",
  },
  {
    id: "p3",
    text: "For Vs=3.3 V, R=220 Ω, Is=1e-12 A, n=2, solve for I and compare Newton vs bisection iterations.",
    hint: "Try bracket [0,0.02] and Newton initial guess near expected current.",
  },
  {
    id: "p4",
    text: "For Re=1.2e5, ε=0.0002, D=0.1 m, compute f with Colebrook and compare with Swamee–Jain.",
    hint: "Use Swamee–Jain as initial guess; compute relative error.",
  },
  {
    id: "p5",
    text: "Advanced: Replace Van der Waals with Peng–Robinson cubic EOS and discuss root multiplicity near phase boundaries.",
    hint: "PR is cubic in V — plot residual over V to inspect multiple root basins.",
  },
];

function Problems() {
  const [state, setState] = useState(defaultProblems.map((p) => ({ ...p, done: false, showHint: false })));

  function toggleDone(id) {
    setState((s) => s.map((p) => (p.id === id ? { ...p, done: !p.done } : p)));
  }
  function toggleHint(id) {
    setState((s) => s.map((p) => (p.id === id ? { ...p, showHint: !p.showHint } : p)));
  }
  function resetAll() {
    setState(defaultProblems.map((p) => ({ ...p, done: false, showHint: false })));
  }

  return (
    <Card className={theme.panel}>
      <CardHeader>
        <CardTitle className="text-emerald-300 flex items-center gap-2">
          <List className="w-5 h-5" />
          Problems & Exercises
        </CardTitle>
      </CardHeader>
      <CardContent>
        <div className="flex items-center justify-between mb-3">
          <div className="text-zinc-300 text-sm">Work through these problems; mark when complete.</div>
          <div className="flex gap-2">
            <Button className="bg-white cursor-pointer hover:bg-gray-200" variant="ghost" size="sm" onClick={resetAll}>
              Reset
            </Button>
          </div>
        </div>

        <div className="space-y-3">
          {state.map((p) => (
            <div key={p.id} className={`p-3 rounded-md ${p.done ? "bg-zinc-800/30" : "bg-zinc-800/10"}`}>
              <div className="flex items-start gap-3">
                <input type="checkbox"  checked={p.done} onChange={() => toggleDone(p.id)} className="mt-1 cursor-pointer" />
                <div className="flex-1">
                  <div className={`text-sm ${p.done ? "line-through text-zinc-400" : "text-zinc-100"}`}>{p.text}</div>
                  <div className="mt-2 flex items-center gap-2">
                    <Button variant="ghost" size="sm" className="bg-white cursor-pointer hover:bg-gray-300" onClick={() => toggleHint(p.id)}>
                      {p.showHint ? "Hide hint" : "Show hint"}
                    </Button>
                  </div>
                  {p.showHint && <div className="mt-2 text-zinc-400 text-sm">{p.hint}</div>}
                </div>
              </div>
            </div>
          ))}
        </div>
      </CardContent>
    </Card>
  );
}

// ======================================================================
// Page assembly — bring all sections together to match Chapter7 layout
// ======================================================================
export default function Chapter8() {
  return (
    <div className={`${theme.bg} p-6 space-y-6 min-h-screen`}>
      <motion.div initial={{ opacity: 0, y: -10 }} animate={{ opacity: 1, y: 0 }}>
        <h1 className="text-3xl font-bold" style={{ color: theme.accent2 }}>
        Case Studies: Roots of Equations
        </h1>
        <div className="text-zinc-400">Practical engineering examples transformed to root-finding problems with interactive visuals.</div>
      </motion.div>

      <OverviewCards />

      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        <div className="space-y-6">
          <Section81 />
          <Section82 />
        </div>

        <div className="space-y-6">
          <Section83 />
          <Section84 />
        </div>
      </div>

      <ChapterDocs />

      <Problems />

      {/* Decorative particle cloud kept for consistency */}
      <Card className={theme.panel}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Atom className="w-5 h-5" />
            Visualization — Particle Halo
          </CardTitle>
        </CardHeader>
        <CardContent>
          <div className="w-full h-48">
            <Canvas camera={{ position: [0, 0, 4] }}>
              <ambientLight intensity={0.6} />
              <pointLight position={[10, 10, 10]} />
              <ParticleHalo points={84} />
              <OrbitControls />
            </Canvas>
          </div>
          <div className="text-zinc-400 text-sm mt-3">A subtle 3D decoration to keep the visual language consistent across chapters.</div>
        </CardContent>
      </Card>
      <BottomBar/>
    </div>
  );
}
