"use client";

/**
 * src/pages/Chapter28.jsx
 *
 * CHAPTER 28 — Case Studies: Ordinary Differential Equations
 * (Responsive, robust, avoids react-markdown className error; responsive Plotly)
 *
 * - Dark theme, Tailwind responsive layouts
 * - Framer Motion micro-animations
 * - React-Markdown with components prop (no className on ReactMarkdown)
 * - Plotly plots wrapped in responsive containers to avoid overflow on small screens
 * - Many sections, extra plots, and improved layout responsiveness
 *
 * NOTE: keep dependencies installed:
 * react, react-dom, react-plotly.js, plotly.js, framer-motion, react-markdown,
 * remark-math, rehype-katex, katex, lucide-react, and your UI components (Card, Input, Button, Textarea).
 *
 * This file is intentionally long and detailed to match your request for an extensive, responsive chapter file.
 */

import React, { useState, useMemo, useEffect } from "react";
import Plot from "react-plotly.js";
import { motion } from "framer-motion";
import ReactMarkdown from "react-markdown";
import remarkMath from "remark-math";
import rehypeKatex from "rehype-katex";
import "katex/dist/katex.min.css";

import {
  BookOpen,
  FlaskConical,
  Activity,
  Zap,
  CircleDot,
  ListChecks,
  ChartBar,
  RotateCw,
  MapPin,
  Archive,
  Settings,
  ChevronDown,
  ChevronRight,
} from "lucide-react";
import BottomBar from "../components/BottomBar";
import { Card, CardContent } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Button } from "@/components/ui/button";
import { Textarea } from "@/components/ui/textarea";

/* ---------------------------------------------
   Theme constants
   --------------------------------------------- */
const theme = {
  bgClass: "bg-zinc-950",
  panelClass: "bg-zinc-900/60 border border-zinc-700",
  primaryColor: "#22d3ee",
  secondaryColor: "#34d399",
  warnColor: "#f59e0b",
  dangerColor: "#ef4444",
};

/* ---------------------------------------------
   Small presentational components
   --------------------------------------------- */
function SectionHeader({ icon: Icon, number, title, accent = theme.primaryColor }) {
  return (
    <div className="flex items-center gap-3">
      <div className="p-2 rounded-lg" style={{ background: "rgba(255,255,255,0.02)" }}>
        <Icon className="w-6 h-6" style={{ color: accent }} />
      </div>
      <div>
        <div className="text-xs text-zinc-400 uppercase tracking-wider">{number}</div>
        <div className="text-lg font-semibold text-zinc-100">{title}</div>
      </div>
    </div>
  );
}

function MD({ children, className = "", inline = false }) {
  // Use components mapping to avoid passing className to ReactMarkdown (v8+)
  const components = {
    p: ({ node, ...props }) => <p className={`text-sm text-zinc-300 leading-relaxed ${className}`} {...props} />,
    strong: ({ node, ...props }) => <strong className="text-emerald-300" {...props} />,
    em: ({ node, ...props }) => <em className="text-zinc-200" {...props} />,
    code: ({ node, inline: isInline, ...props }) =>
      isInline ? (
        <code className="bg-zinc-800/70 px-1 rounded text-xs text-cyan-200" {...props} />
      ) : (
        <pre className="bg-zinc-800/60 rounded p-2 text-xs overflow-auto" {...props} />
      ),
    li: ({ node, ...props }) => <li className="ml-4 text-sm text-zinc-300" {...props} />,
    h1: ({ node, ...props }) => <h1 className="text-xl text-zinc-100 font-semibold" {...props} />,
    h2: ({ node, ...props }) => <h2 className="text-lg text-zinc-100 font-semibold" {...props} />,
  };

  return (
    <div className="w-full">
      <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]} components={components}>
        {children}
      </ReactMarkdown>
    </div>
  );
}

/* ---------------------------------------------
   Utility numerical routines: RK4, linspace
   --------------------------------------------- */
function linspace(a, b, n) {
  if (n <= 1) return [a];
  const out = [];
  const dx = (b - a) / (n - 1);
  for (let i = 0; i < n; i++) out.push(a + i * dx);
  return out;
}

function rk4StepVec(f, t, y, h) {
  const k1 = f(t, y);
  const k2 = f(t + 0.5 * h, y.map((yi, i) => yi + 0.5 * h * k1[i]));
  const k3 = f(t + 0.5 * h, y.map((yi, i) => yi + 0.5 * h * k2[i]));
  const k4 = f(t + h, y.map((yi, i) => yi + h * k3[i]));
  return y.map((yi, i) => yi + (h / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]));
}

function integrateRK4(f, y0, t0, tf, h) {
  const nSteps = Math.max(1, Math.ceil((tf - t0) / h));
  const T = new Array(nSteps + 1);
  const Y = new Array(nSteps + 1);
  let t = t0;
  let y = y0.slice();
  T[0] = t;
  Y[0] = y.slice();
  for (let k = 1; k <= nSteps; k++) {
    const hh = Math.min(h, tf - t);
    y = rk4StepVec(f, t, y, hh);
    t += hh;
    T[k] = t;
    Y[k] = y.slice();
  }
  return { T, Y };
}

/* ---------------------------------------------
   Responsive Plot wrapper
   - ensures plots don't overflow on small screens
   - useResizeHandler true and container with max-w-full/overflow-hidden
   --------------------------------------------- */
function ResponsivePlot({ children, height = 360, title = "" }) {
  return (
    <div className="w-full max-w-full overflow-hidden rounded-lg">
      <div className="w-full" style={{ minHeight: 120 }}>
        {/* The Plot itself will set width/height; here we ensure the container is clipped */}
        {children}
      </div>
    </div>
  );
}

/* ---------------------------------------------------------------------------
   Section 28.1 — Reactor transient (extended, responsive)
   - first & second order, numeric vs analytic, rate/time surfaces
   --------------------------------------------------------------------------- */
function ReactorSection() {
  const [k1, setK1] = useState(0.8);
  const [C0, setC0] = useState(1.0);
  const [tmax, setTmax] = useState(8);
  const [h, setH] = useState(0.02);
  const [k2, setK2] = useState(0.5);
  const [useSecondOrder, setUseSecondOrder] = useState(false);
  const [showRateSurface, setShowRateSurface] = useState(false);

  // Analytic (first-order)
  const analyticFirst = useMemo(() => {
    const N = Math.max(2, Math.ceil(tmax / h));
    const T = linspace(0, tmax, N + 1);
    const C = T.map((t) => C0 * Math.exp(-k1 * t));
    return { T, C };
  }, [k1, C0, tmax, h]);

  // Numeric RK4 (general)
  const numeric = useMemo(() => {
    if (!useSecondOrder) {
      const f = (t, y) => [-k1 * y[0]];
      const res = integrateRK4(f, [C0], 0, tmax, h);
      return { T: res.T, C: res.Y.map((yy) => yy[0]) };
    } else {
      const f = (t, y) => [-k2 * y[0] * y[0]];
      const res = integrateRK4(f, [C0], 0, tmax, h);
      const analytic = res.T.map((t) => 1 / (1 / C0 + k2 * t));
      return { T: res.T, C: res.Y.map((yy) => yy[0]), analyticSecond: analytic };
    }
  }, [k1, k2, C0, tmax, h, useSecondOrder]);

  // Rate surface data (k vs C0 slice at t=0) simple contour to visualize rate magnitude
  const rateSurface = useMemo(() => {
    const kVals = linspace(Math.max(0.01, k1 - 1), k1 + 1, 36);
    const Cvals = linspace(Math.max(0.01, C0 - 1), C0 + 1, 36);
    const Z = Cvals.map((Cv) => kVals.map((kv) => -kv * Cv));
    return { x: kVals, y: Cvals, z: Z };
  }, [k1, C0]);

  const times = numeric.T;
  const numC = numeric.C;

  return (
    <motion.div
      className={`${theme.panelClass} p-4 rounded-2xl`}
      initial={{ opacity: 0, y: 8 }}
      whileInView={{ opacity: 1, y: 0 }}
      viewport={{ once: true }}
    >
      <div className="flex items-start justify-between">
        <SectionHeader icon={FlaskConical} number="28.1" title="Reactor Transient Response" />
        <div className="space-x-2 flex">
          <button
            className="text-xs px-2 py-1 bg-emerald-500 cursor-pointer rounded border border-zinc-700 text-zinc-100"
            onClick={() => setUseSecondOrder((s) => !s)}
          >
            {useSecondOrder ? "Switch to 1st-order" : "Switch to 2nd-order"}
          </button>
          <button
            className="text-xs px-2 py-1 rounded  bg-emerald-500 cursor-pointer border border-zinc-700 text-zinc-100"
            onClick={() => setShowRateSurface((s) => !s)}
          >
            {showRateSurface ? "Hide rate surface" : "Show rate surface"}
          </button>
        </div>
      </div>

      <div className="grid gap-4 mt-4 md:grid-cols-3">
        <div className="space-y-3 md:col-span-1">
          <MD>
  {useSecondOrder
    ? `**Second-order kinetics** — $\\frac{dC}{dt} = -k_2 C^2$. Analytic solution: $C(t)=\\frac{1}{1/C_0 + k_2 t}$.`
    : `**First-order kinetics** — $\\frac{dC}{dt} = -k_1 C$. Analytic solution: $C(t)=C_0 e^{-k_1 t}$.`}
</MD>


          <div className="grid grid-cols-2 gap-2">
            <div>
              <label className="text-xs text-zinc-400">k₁</label>
              <Input className="text-white" type="number" value={k1} onChange={(e) => setK1(parseFloat(e.target.value || 0.0))} />
            </div>
            <div>
              <label className="text-xs text-zinc-400">k₂</label>
              <Input className="text-white"  type="number" value={k2} onChange={(e) => setK2(parseFloat(e.target.value || 0.0))} />
            </div>
            <div>
              <label className="text-xs text-zinc-400">C₀</label>
              <Input className="text-white"  type="number" value={C0} onChange={(e) => setC0(parseFloat(e.target.value || 0.0))} />
            </div>
            <div>
              <label className="text-xs text-zinc-400">tₘₐₓ</label>
              <Input className="text-white"  type="number" value={tmax} onChange={(e) => setTmax(parseFloat(e.target.value || 0.0))} />
            </div>
            <div>
              <label className="text-xs text-zinc-400">step h</label>
              <Input className="text-white"  type="number" value={h} onChange={(e) => setH(parseFloat(e.target.value || 0.01))} />
            </div>
            <div />
          </div>

          <MD>
            {"**Note:** Compare numerical RK4 to analytic solution where available. For first-order kinetics analytic is `C(t)=C₀ e^{-k₁ t}`. For second-order with constant k₂, analytic is `C(t)=1/(1/C₀ + k₂ t)`."}
          </MD>
        </div>

        <div className="md:col-span-2 space-y-3">
          <div className="rounded-lg border border-zinc-700 p-2 bg-zinc-900/50">
            <ResponsivePlot height={360} title="Concentration vs Time (RK4 vs Analytic)">
              <Plot
                data={[
                  {
                    x: times,
                    y: numC,
                    mode: "lines",
                    name: "RK4 (numeric)",
                    line: { width: 2, color: theme.primaryColor },
                  },
                  !useSecondOrder && {
                    x: analyticFirst.T,
                    y: analyticFirst.C,
                    mode: "lines",
                    name: "Analytic 1st-order",
                    line: { dash: "dash", width: 2, color: theme.secondaryColor },
                  },
                  useSecondOrder && numeric.analyticSecond && {
                    x: numeric.T,
                    y: numeric.analyticSecond,
                    mode: "lines",
                    name: "Analytic 2nd-order",
                    line: { dash: "dash", width: 2, color: theme.secondaryColor },
                  },
                ].filter(Boolean)}
                layout={{
                  title: "",
                  autosize: true,
                  paper_bgcolor: "rgba(0,0,0,0)",
                  plot_bgcolor: "rgba(0,0,0,0)",
                  xaxis: { title: "t", color: "#9ca3af" },
                  yaxis: { title: "C(t)", color: "#9ca3af" },
                  margin: { l: 50, r: 20, t: 30, b: 40 },
                  legend: { orientation: "h", y: -0.2 },
                }}
                style={{ width: "100%", height: 360 }}
                useResizeHandler
                config={{ responsive: true }}
              />
            </ResponsivePlot>

            <div className="grid md:grid-cols-2 gap-3 mt-3">
              <div className="rounded-lg border border-zinc-700 p-2 bg-zinc-900/50">
                <h4 className="text-sm text-zinc-100 mb-2">Reaction Rate vs Time</h4>
                <ResponsivePlot height={260}>
                  <Plot
                    data={[
                      {
                        x: times,
                        y: times.map((t, i) => {
                          const c = numC[i];
                          return useSecondOrder ? -k2 * c * c : -k1 * c;
                        }),
                        mode: "lines",
                        name: "rate",
                        line: { width: 2, color: theme.warnColor },
                      },
                    ]}
                    layout={{
                      autosize: true,
                      paper_bgcolor: "rgba(0,0,0,0)",
                      plot_bgcolor: "rgba(0,0,0,0)",
                      xaxis: { title: "t", color: "#9ca3af" },
                      yaxis: { title: "r(t)", color: "#9ca3af" },
                      margin: { l: 40, r: 10, t: 20, b: 30 },
                    }}
                    style={{ width: "100%", height: 260 }}
                    useResizeHandler
                    config={{ responsive: true }}
                  />
                </ResponsivePlot>
              </div>

              <div className="rounded-lg border border-zinc-700 p-2 bg-zinc-900/50">
                <h4 className="text-sm text-zinc-100 mb-2">Rate Surface (k vs C₀)</h4>
                {showRateSurface ? (
                  <ResponsivePlot height={260}>
                    <Plot
                      data={[
                        {
                          x: rateSurface.x,
                          y: rateSurface.y,
                          z: rateSurface.z,
                          type: "contour",
                          colorscale: "Viridis",
                        },
                      ]}
                      layout={{
                        autosize: true,
                        paper_bgcolor: "rgba(0,0,0,0)",
                        plot_bgcolor: "rgba(0,0,0,0)",
                        margin: { l: 30, r: 10, t: 20, b: 30 },
                      }}
                      style={{ width: "100%", height: 260 }}
                      useResizeHandler
                      config={{ responsive: true }}
                    />
                  </ResponsivePlot>
                ) : (
                  <div className="text-sm text-zinc-400">Toggle the rate surface to view parameter dependence.</div>
                )}
              </div>
            </div>
          </div>
        </div>
      </div>
    </motion.div>
  );
}

/* ---------------------------------------------------------------------------
   Section 28.2 — Predator–Prey & Chaos (responsive)
   - Lotka–Volterra, logistic prey, phase plane, energy-like invariant,
   - optional Lorenz toy 3D plot
   --------------------------------------------------------------------------- */
function PredatorPreySection() {
  const [alpha, setAlpha] = useState(1.1);
  const [beta, setBeta] = useState(0.4);
  const [delta, setDelta] = useState(0.1);
  const [gamma, setGamma] = useState(0.4);
  const [tmax, setTmax] = useState(120);
  const [h, setH] = useState(0.02);
  const [initial, setInitial] = useState("40,9");
  const [useLogistic, setUseLogistic] = useState(false);
  const [K, setK] = useState(100);
  const [showEnergy, setShowEnergy] = useState(false);
  const [showLorenzToy, setShowLorenzToy] = useState(false);

  const [x0, z0] = useMemo(() => {
    const parts = initial.split(",").map((s) => parseFloat(s.trim()));
    return [isFinite(parts[0]) ? parts[0] : 40, isFinite(parts[1]) ? parts[1] : 9];
  }, [initial]);

  const result = useMemo(() => {
    const f = (t, y) => {
      const x = y[0];
      const z = y[1];
      if (!useLogistic) return [alpha * x - beta * x * z, delta * x * z - gamma * z];
      return [alpha * x * (1 - x / K) - beta * x * z, delta * x * z - gamma * z];
    };
    const res = integrateRK4(f, [x0, z0], 0, tmax, h);
    return { T: res.T, X: res.Y.map((v) => v[0]), Z: res.Y.map((v) => v[1]) };
  }, [alpha, beta, delta, gamma, x0, z0, tmax, h, useLogistic, K]);

  const energy = useMemo(() => {
    if (useLogistic) return null;
    const eps = 1e-8;
    return result.T.map((_, i) => {
      const x = Math.max(eps, result.X[i]);
      const z = Math.max(eps, result.Z[i]);
      return delta * x - gamma * Math.log(x) + beta * z - alpha * Math.log(z);
    });
  }, [result, useLogistic, alpha, beta, delta, gamma]);

  const lorenzToy = useMemo(() => {
    if (!showLorenzToy) return null;
    const sigma = 10,
      rho = 28,
      betaL = 8 / 3;
    const f = (t, y) => {
      const [u, v, w] = y;
      return [sigma * (v - u), u * (rho - w) - v, u * v - betaL * w];
    };
    const res = integrateRK4(f, [1.0, 1.0, 1.0], 0, 40, 0.01);
    return res;
  }, [showLorenzToy]);

  return (
    <motion.div className={`${theme.panelClass} p-4 rounded-2xl`} initial={{ opacity: 0, y: 8 }} whileInView={{ opacity: 1, y: 0 }} viewport={{ once: true }}>
      <div className="flex items-start justify-between">
        <SectionHeader icon={Activity} number="28.2" title="Predator–Prey Models & Chaos" accent={theme.secondaryColor} />
        <div className="flex gap-2">
          <button className="text-xs px-2 py-1 bg-emerald-500 cursor-pointer rounded border border-zinc-700" onClick={() => setUseLogistic((s) => !s)}>
            {useLogistic ? "Classic LV" : "Logistic prey"}
          </button>
          <button className="text-xs px-2 py-1 bg-emerald-500 cursor-pointer rounded border border-zinc-700" onClick={() => setShowEnergy((s) => !s)}>
            {showEnergy ? "Hide energy" : "Show energy"}
          </button>
          <button className="text-xs px-2 py-1 bg-emerald-500 cursor-pointer rounded border border-zinc-700" onClick={() => setShowLorenzToy((s) => !s)}>
            {showLorenzToy ? "Hide Lorenz" : "Show Lorenz toy"}
          </button>
        </div>
      </div>

      <div className="grid gap-4 mt-4 md:grid-cols-3">
        <div className="space-y-3 md:col-span-1">
          <MD>{"**Model**: Lotka–Volterra or logistic-modified prey. Adjust parameters & initial state."}</MD>

          <div className="grid grid-cols-2 gap-2">
            <div>
              <label className="text-xs text-zinc-400">α</label>
              <Input className="text-white" type="number" value={alpha} onChange={(e) => setAlpha(parseFloat(e.target.value || 0.0))} />
            </div>
            <div>
              <label className="text-xs text-zinc-400">β</label>
              <Input className="text-white" type="number" value={beta} onChange={(e) => setBeta(parseFloat(e.target.value || 0.0))} />
            </div>
            <div>
              <label className="text-xs text-zinc-400">δ</label>
              <Input className="text-white" type="number" value={delta} onChange={(e) => setDelta(parseFloat(e.target.value || 0.0))} />
            </div>
            <div>
              <label className="text-xs text-zinc-400">γ</label>
              <Input className="text-white" type="number" value={gamma} onChange={(e) => setGamma(parseFloat(e.target.value || 0.0))} />
            </div>
            <div>
              <label className="text-xs text-zinc-400">Initial (x,z)</label>
              <Input className="text-white" value={initial} onChange={(e) => setInitial(e.target.value)} />
            </div>
            <div>
              <label className="text-xs text-zinc-400">tₘₐₓ</label>
              <Input className="text-white" type="number" value={tmax} onChange={(e) => setTmax(parseFloat(e.target.value || 100))} />
            </div>
            <div>
              <label className="text-xs text-zinc-400">h</label>
              <Input className="text-white" type="number" value={h} onChange={(e) => setH(parseFloat(e.target.value || 0.02))} />
            </div>
            <div>
              <label className="text-xs text-zinc-400">K (if logistic)</label>
              <Input className="text-white" type="number" value={K} onChange={(e) => setK(parseFloat(e.target.value || 100))} />
            </div>
          </div>

          <MD>
            {"**Note:** The classical Lotka–Volterra system has neutrally stable closed orbits. Adding logistic growth or harvesting terms changes stability and can introduce limit cycles or damped behavior."}
          </MD>
        </div>

        <div className="md:col-span-2 space-y-3">
          <div className="rounded-lg border border-zinc-700 p-2 bg-zinc-900/50">
            <ResponsivePlot>
              <Plot
                data={[
                  { x: result.T, y: result.X, mode: "lines", name: "Prey x(t)", line: { color: theme.primaryColor } },
                  { x: result.T, y: result.Z, mode: "lines", name: "Predator z(t)", line: { color: theme.warnColor } },
                ]}
                layout={{
                  title: "Populations vs Time",
                  autosize: true,
                  paper_bgcolor: "rgba(0,0,0,0)",
                  plot_bgcolor: "rgba(0,0,0,0)",
                  xaxis: { title: "t", color: "#9ca3af" },
                  yaxis: { title: "population", color: "#9ca3af" },
                  margin: { t: 30, b: 30, l: 50, r: 10 },
                  legend: { orientation: "h", y: -0.2 },
                }}
                style={{ width: "100%", height: 300 }}
                useResizeHandler
                config={{ responsive: true }}
              />
            </ResponsivePlot>

            <div className="grid md:grid-cols-2 gap-3 mt-3">
              <div className="rounded-lg border border-zinc-700 p-2 bg-zinc-900/50">
                <ResponsivePlot height={320}>
                  <Plot
                    data={[{ x: result.X, y: result.Z, mode: "lines", line: { color: theme.secondaryColor } }]}
                    layout={{
                      title: "Phase plane (x vs z)",
                      xaxis: { title: "Prey x" },
                      yaxis: { title: "Predator z" },
                      paper_bgcolor: "rgba(0,0,0,0)",
                      plot_bgcolor: "rgba(0,0,0,0)",
                      margin: { l: 50, r: 10, t: 30, b: 40 },
                    }}
                    style={{ width: "100%", height: 320 }}
                    useResizeHandler
                    config={{ responsive: true }}
                  />
                </ResponsivePlot>
              </div>

              <div className="rounded-lg border border-zinc-700 p-2 bg-zinc-900/50">
                {showEnergy && energy ? (
                  <ResponsivePlot height={320}>
                    <Plot
                      data={[{ x: result.T, y: energy, mode: "lines", line: { color: "#f59e0b" } }]}
                      layout={{
                        title: "Invariant (energy-like) for Lotka–Volterra",
                        paper_bgcolor: "rgba(0,0,0,0)",
                        plot_bgcolor: "rgba(0,0,0,0)",
                        margin: { l: 50, r: 10, t: 30, b: 40 },
                      }}
                      style={{ width: "100%", height: 320 }}
                      useResizeHandler
                      config={{ responsive: true }}
                    />
                  </ResponsivePlot>
                ) : (
                  <div className="text-sm text-zinc-400 p-2">Toggle energy to show conserved quantity (LV only).</div>
                )}
              </div>
            </div>

            {showLorenzToy && lorenzToy && (
              <div className="rounded-lg border border-zinc-700 p-2 bg-zinc-900/50 mt-3">
                <h4 className="text-sm text-zinc-100 mb-2">Lorenz-like Toy (3D chaotic attractor)</h4>
                <ResponsivePlot>
                  <Plot
                    data={[
                      {
                        x: lorenzToy.T.map((_, i) => lorenzToy.Y[i][0]),
                        y: lorenzToy.T.map((_, i) => lorenzToy.Y[i][1]),
                        z: lorenzToy.T.map((_, i) => lorenzToy.Y[i][2]),
                        type: "scatter3d",
                        mode: "lines",
                        line: { width: 2, color: "#f97316" },
                      },
                    ]}
                    layout={{
                      title: "",
                      autosize: true,
                      paper_bgcolor: "rgba(0,0,0,0)",
                      margin: { l: 0, r: 0, t: 20, b: 20 },
                    }}
                    style={{ width: "100%", height: 420 }}
                    useResizeHandler
                    config={{ responsive: true }}
                  />
                </ResponsivePlot>
              </div>
            )}
          </div>
        </div>
      </div>
    </motion.div>
  );
}

/* ---------------------------------------------------------------------------
   Section 28.3 — RLC Transient (responsive)
   - integrate q and i, analytic underdamped comparison, sweep R
   --------------------------------------------------------------------------- */
function RLCCircuitSection() {
  const [L, setL] = useState(1.0);
  const [R, setR] = useState(0.5);
  const [C, setC] = useState(0.2);
  const [q0, setQ0] = useState(1.0);
  const [i0, setI0] = useState(0.0);
  const [tmax, setTmax] = useState(40);
  const [h, setH] = useState(0.02);
  const [showSweep, setShowSweep] = useState(false);

  const omega0 = useMemo(() => 1 / Math.sqrt(L * C), [L, C]);
  const zeta = useMemo(() => R / (2 * Math.sqrt(L / C)), [R, L, C]);
  const dampingType = useMemo(() => (zeta < 1 ? "underdamped" : zeta === 1 ? "critical" : "overdamped"), [zeta]);

  const result = useMemo(() => {
    const f = (t, y) => {
      const q = y[0],
        i = y[1];
      return [i, -(R / L) * i - (1 / (L * C)) * q];
    };
    const res = integrateRK4(f, [q0, i0], 0, tmax, h);
    return { T: res.T, Q: res.Y.map((v) => v[0]), I: res.Y.map((v) => v[1]) };
  }, [L, R, C, q0, i0, tmax, h]);

  const analytic = useMemo(() => {
    if (zeta >= 1) return null;
    const wn = omega0;
    const wd = wn * Math.sqrt(1 - zeta * zeta);
    const A = q0;
    const B = (i0 + zeta * wn * q0) / wd;
    const Qana = result.T.map((t) => Math.exp(-zeta * wn * t) * (A * Math.cos(wd * t) + B * Math.sin(wd * t)));
    const Iana = result.T.map((t) => {
      // compute derivative of Qana (careful with algebra)
      const qval = Math.exp(-zeta * wn * t) * (A * Math.cos(wd * t) + B * Math.sin(wd * t));
      const dqdt =
        -zeta * wn * qval +
        Math.exp(-zeta * wn * t) * (-A * wd * Math.sin(wd * t) + B * wd * Math.cos(wd * t));
      return dqdt;
    });
    return { Qana, Iana };
  }, [zeta, omega0, q0, i0, result.T]);

  const sweepData = useMemo(() => {
    if (!showSweep) return null;
    const Rs = linspace(Math.max(0.01, R - 3), R + 3, 7);
    return Rs.map((Rval) => {
      const f = (t, y) => {
        const q = y[0],
          i = y[1];
        return [i, -(Rval / L) * i - (1 / (L * C)) * q];
      };
      const res = integrateRK4(f, [q0, i0], 0, Math.min(20, tmax), Math.max(0.01, h));
      return { x: res.T, y: res.Y.map((v) => v[0]), name: `R=${Rval.toFixed(2)}` };
    });
  }, [showSweep, R, L, C, q0, i0, tmax, h]);

  return (
    <motion.div className={`${theme.panelClass} p-4 rounded-2xl`} initial={{ opacity: 0, y: 8 }} whileInView={{ opacity: 1, y: 0 }} viewport={{ once: true }}>
      <div className="flex items-start justify-between">
        <SectionHeader icon={Zap} number="28.3" title="Transient Current: R-L-C Series Circuit" accent={theme.warnColor} />
        <div className="text-sm text-zinc-300">Regime: <span className="text-zinc-100">{dampingType}</span></div>
      </div>

      <div className="grid gap-4 mt-4 md:grid-cols-3">
        <div className="space-y-3 md:col-span-1">
          <MD>
  {"**Equation:** $L q'' + R q' + \\frac{1}{C} q = 0$. Integrate for $q(t)$ and $i(t) = q'(t)$."}
</MD>


          <div className="grid grid-cols-2 gap-2">
            <div>
              <label className="text-xs text-zinc-400">L</label>
              <Input className="text-white" type="number" value={L} onChange={(e) => setL(parseFloat(e.target.value || 1))} />
            </div>
            <div>
              <label className="text-xs text-zinc-400">R</label>
              <Input className="text-white" type="number" value={R} onChange={(e) => setR(parseFloat(e.target.value || 0.5))} />
            </div>
            <div>
              <label className="text-xs text-zinc-400">C</label>
              <Input className="text-white" type="number" value={C} onChange={(e) => setC(parseFloat(e.target.value || 0.2))} />
            </div>
            <div>
              <label className="text-xs text-zinc-400">q₀</label>
              <Input className="text-white" type="number" value={q0} onChange={(e) => setQ0(parseFloat(e.target.value || 1))} />
            </div>
            <div>
              <label className="text-xs text-zinc-400">i₀</label>
              <Input className="text-white" type="number" value={i0} onChange={(e) => setI0(parseFloat(e.target.value || 0))} />
            </div>
            <div>
              <label className="text-xs text-zinc-400">tₘₐₓ</label>
              <Input className="text-white" type="number" value={tmax} onChange={(e) => setTmax(parseFloat(e.target.value || 40))} />
            </div>
            <div>
              <label className="text-xs text-zinc-400">h</label>
              <Input className="text-white" type="number" value={h} onChange={(e) => setH(parseFloat(e.target.value || 0.02))} />
            </div>
            <div />
          </div>

          <div className="flex gap-2">
            <button className="text-xs px-2 py-1 bg-emerald-500 cursor-pointer rounded border" onClick={() => setShowSweep((s) => !s)}>
              {showSweep ? "Hide sweep" : "Show R sweep"}
            </button>
          </div>

          <MD>
            {`Natural frequency: \\(\\omega_0 = 1/\\sqrt{LC}\\). Damping ratio (textbook form) relates to R, L, C and controls transient type.`}
          </MD>
        </div>

        <div className="md:col-span-2 space-y-3">
          <div className="rounded-lg border border-zinc-700 p-2 bg-zinc-900/50">
            <ResponsivePlot>
              <Plot
                data={[
                  { x: result.T, y: result.Q, mode: "lines", name: "q(t) numeric", line: { color: theme.primaryColor } },
                  analytic && { x: result.T, y: analytic.Qana, mode: "lines", name: "q(t) analytic", line: { dash: "dash", color: theme.secondaryColor } },
                ].filter(Boolean)}
                layout={{
                  title: "Charge q(t)",
                  autosize: true,
                  paper_bgcolor: "rgba(0,0,0,0)",
                  plot_bgcolor: "rgba(0,0,0,0)",
                  margin: { t: 30, b: 30, l: 50, r: 10 },
                }}
                style={{ width: "100%", height: 300 }}
                useResizeHandler
                config={{ responsive: true }}
              />
            </ResponsivePlot>

            <div className="grid md:grid-cols-2 gap-3 mt-3">
              <div className="rounded-lg border border-zinc-700 p-2 bg-zinc-900/50">
                <ResponsivePlot>
                  <Plot
                    data={[
                      { x: result.T, y: result.I, mode: "lines", name: "i(t) numeric", line: { color: theme.warnColor } },
                      analytic && { x: result.T, y: analytic.Iana, mode: "lines", name: "i(t) analytic", line: { dash: "dash", color: "#f97316" } },
                    ].filter(Boolean)}
                    layout={{
                      title: "Current i(t)",
                      autosize: true,
                      paper_bgcolor: "rgba(0,0,0,0)",
                      plot_bgcolor: "rgba(0,0,0,0)",
                      margin: { t: 30, b: 30, l: 40, r: 10 },
                    }}
                    style={{ width: "100%", height: 280 }}
                    useResizeHandler
                    config={{ responsive: true }}
                  />
                </ResponsivePlot>
              </div>

              <div className="rounded-lg border border-zinc-700 p-2 bg-zinc-900/50">
                <ResponsivePlot>
                  <Plot
                    data={[{ x: result.Q, y: result.I, mode: "lines", name: "phase", line: { color: "#a78bfa" } }]}
                    layout={{
                      title: "Phase portrait (q vs i)",
                      xaxis: { title: "q" },
                      yaxis: { title: "i" },
                      paper_bgcolor: "rgba(0,0,0,0)",
                      plot_bgcolor: "rgba(0,0,0,0)",
                      margin: { t: 30, b: 30, l: 40, r: 10 },
                    }}
                    style={{ width: "100%", height: 280 }}
                    useResizeHandler
                    config={{ responsive: true }}
                  />
                </ResponsivePlot>
              </div>
            </div>

            {showSweep && sweepData && (
              <div className="rounded-lg border border-zinc-700 p-2 bg-zinc-900/50 mt-3">
                <ResponsivePlot>
                  <Plot
                    data={sweepData.map((s) => ({ x: s.x, y: s.y, mode: "lines", name: s.name }))}
                    layout={{
                      title: "Sweep: q(t) for different R",
                      autosize: true,
                      paper_bgcolor: "rgba(0,0,0,0)",
                      plot_bgcolor: "rgba(0,0,0,0)",
                      margin: { t: 30, b: 30, l: 40, r: 10 },
                    }}
                    style={{ width: "100%", height: 320 }}
                    useResizeHandler
                    config={{ responsive: true }}
                  />
                </ResponsivePlot>
              </div>
            )}
          </div>
        </div>
      </div>
    </motion.div>
  );
}

/* ---------------------------------------------------------------------------
   Section 28.4 — Pendulum (responsive)
   - nonlinear vs small-angle, energy, phase portrait
   --------------------------------------------------------------------------- */
function PendulumSection() {
  const [theta0, setTheta0] = useState(0.5);
  const [omega0, setOmega0] = useState(0.0);
  const [tmax, setTmax] = useState(20);
  const [h, setH] = useState(0.01);
  const [g, setG] = useState(9.81);
  const [showEnergy, setShowEnergy] = useState(true);
  const [compareSmallAngle, setCompareSmallAngle] = useState(true);

  const L = 1.0;

  const nonlinear = useMemo(() => {
    const f = (t, y) => {
      const theta = y[0],
        omega = y[1];
      return [omega, -(g / L) * Math.sin(theta)];
    };
    const res = integrateRK4(f, [theta0, omega0], 0, tmax, h);
    return { T: res.T, Theta: res.Y.map((v) => v[0]), Omega: res.Y.map((v) => v[1]) };
  }, [theta0, omega0, tmax, h, g]);

  const smallAngle = useMemo(() => {
    if (!compareSmallAngle) return null;
    const wn = Math.sqrt(g / L);
    const A = theta0;
    const B = omega0 / wn;
    const T = nonlinear.T;
    const ThetaLin = T.map((t) => A * Math.cos(wn * t) + B * Math.sin(wn * t));
    const OmegaLin = T.map((t) => -A * wn * Math.sin(wn * t) + B * wn * Math.cos(wn * t));
    return { T, ThetaLin, OmegaLin, wn };
  }, [compareSmallAngle, nonlinear, theta0, omega0, g, L]);

  const energy = useMemo(() => {
    if (!showEnergy) return null;
    return nonlinear.T.map((_, i) => 0.5 * nonlinear.Omega[i] * nonlinear.Omega[i] + g * (1 - Math.cos(nonlinear.Theta[i])));
  }, [nonlinear, showEnergy, g]);

  return (
    <motion.div className={`${theme.panelClass} p-4 rounded-2xl`} initial={{ opacity: 0, y: 8 }} whileInView={{ opacity: 1, y: 0 }} viewport={{ once: true }}>
      <div className="flex items-start justify-between">
        <SectionHeader icon={CircleDot} number="28.4" title="The Swinging Pendulum" accent={theme.dangerColor} />
        <div className="flex gap-2">
          <button className="text-xs px-2 py-1 rounded bg-emerald-500 cursor-pointer border border-zinc-700" onClick={() => setShowEnergy((s) => !s)}>
            {showEnergy ? "Hide energy" : "Show energy"}
          </button>
          <button className="text-xs px-2 py-1 bg-emerald-500 cursor-pointer rounded border border-zinc-700" onClick={() => setCompareSmallAngle((s) => !s)}>
            {compareSmallAngle ? "Hide small-angle" : "Show small-angle"}
          </button>
        </div>
      </div>

      <div className="grid gap-4 mt-4 md:grid-cols-3">
        <div className="space-y-3 md:col-span-1">
          <MD>
  {`**Nonlinear pendulum**: $\\ddot{\\theta} + \\frac{g}{L} \\sin(\\theta) = 0$. 
Small-angle linearization: $\\ddot{\\theta} + \\frac{g}{L} \\theta = 0$.`}
</MD>


          <div className="grid grid-cols-2 gap-2">
            <div>
              <label className="text-xs text-zinc-400">θ₀ (rad)</label>
              <Input className="text-white" type="number" value={theta0} onChange={(e) => setTheta0(parseFloat(e.target.value || 0))} />
            </div>
            <div>
              <label className="text-xs text-zinc-400">ω₀ (rad/s)</label>
              <Input className="text-white" type="number" value={omega0} onChange={(e) => setOmega0(parseFloat(e.target.value || 0))} />
            </div>
            <div>
              <label className="text-xs text-zinc-400">g</label>
              <Input className="text-white" type="number" value={g} onChange={(e) => setG(parseFloat(e.target.value || 9.81))} />
            </div>
            <div>
              <label className="text-xs text-zinc-400">tₘₐₓ</label>
              <Input className="text-white" type="number" value={tmax} onChange={(e) => setTmax(parseFloat(e.target.value || 20))} />
            </div>
            <div>
              <label className="text-xs text-zinc-400">h</label>
              <Input className="text-white" type="number" value={h} onChange={(e) => setH(parseFloat(e.target.value || 0.01))} />
            </div>
            <div />
          </div>

          <MD>{"**Note:** Energy conservation (no damping) indicates integrator fidelity. Phase portraits visualize stability and separatrix behavior."}</MD>
        </div>

        <div className="md:col-span-2 space-y-3">
          <div className="rounded-lg border border-zinc-700 p-2 bg-zinc-900/50">
            <ResponsivePlot>
              <Plot
                data={[
                  { x: nonlinear.T, y: nonlinear.Theta, mode: "lines", name: "θ(t) nonlinear", line: { color: theme.primaryColor } },
                  compareSmallAngle && smallAngle && { x: smallAngle.T, y: smallAngle.ThetaLin, mode: "lines", name: "θ(t) linear", line: { dash: "dash", color: theme.secondaryColor } },
                ].filter(Boolean)}
                layout={{
                  title: "Pendulum angle vs time",
                  autosize: true,
                  paper_bgcolor: "rgba(0,0,0,0)",
                  plot_bgcolor: "rgba(0,0,0,0)",
                  xaxis: { title: "t", color: "#9ca3af" },
                  yaxis: { title: "θ (rad)", color: "#9ca3af" },
                  margin: { t: 30, b: 30, l: 40, r: 10 },
                }}
                style={{ width: "100%", height: 320 }}
                useResizeHandler
                config={{ responsive: true }}
              />
            </ResponsivePlot>

            <div className="grid md:grid-cols-2 gap-3 mt-3">
              <div className="rounded-lg border border-zinc-700 p-2 bg-zinc-900/50">
                <ResponsivePlot>
                  <Plot
                    data={[
                      { x: nonlinear.Theta, y: nonlinear.Omega, mode: "lines", name: "phase", line: { color: "#f97316" } },
                      compareSmallAngle && smallAngle && { x: smallAngle.ThetaLin, y: smallAngle.OmegaLin, mode: "lines", name: "phase linear", line: { dash: "dash", color: "#f59e0b" } },
                    ].filter(Boolean)}
                    layout={{
                      title: "Phase portrait (θ vs ω)",
                      xaxis: { title: "θ" },
                      yaxis: { title: "ω" },
                      paper_bgcolor: "rgba(0,0,0,0)",
                      plot_bgcolor: "rgba(0,0,0,0)",
                      margin: { t: 30, b: 30, l: 40, r: 10 },
                    }}
                    style={{ width: "100%", height: 300 }}
                    useResizeHandler
                    config={{ responsive: true }}
                  />
                </ResponsivePlot>
              </div>

              <div className="rounded-lg border border-zinc-700 p-2 bg-zinc-900/50">
                {showEnergy && energy ? (
                  <ResponsivePlot>
                    <Plot
                      data={[{ x: nonlinear.T, y: energy, mode: "lines", name: "Energy", line: { color: "#a78bfa" } }]}
                      layout={{
                        title: "Total energy (K + V)",
                        autosize: true,
                        paper_bgcolor: "rgba(0,0,0,0)",
                        plot_bgcolor: "rgba(0,0,0,0)",
                        margin: { t: 30, b: 30, l: 40, r: 10 },
                      }}
                      style={{ width: "100%", height: 300 }}
                      useResizeHandler
                      config={{ responsive: true }}
                    />
                  </ResponsivePlot>
                ) : (
                  <div className="text-sm text-zinc-400 p-2">Toggle energy to visualise total mechanical energy.</div>
                )}
              </div>
            </div>
          </div>
        </div>
      </div>
    </motion.div>
  );
}

/* ---------------------------------------------------------------------------
   Problems section (interactive, responsive)
   --------------------------------------------------------------------------- */
function ProblemsSection() {
  const [openIndex, setOpenIndex] = useState(-1);

const problems = [
  {
    id: 1,
    title: "Parallel reactions (A→B and A→C)",
    prompt: "Write ODEs for $[A]$, $[B]$, $[C]$ and simulate for $k_1 = 0.5$, $k_2 = 0.2$. Plot concentrations and verify mass balance.",
    solution:
      "$\\frac{d[A]}{dt} = -(k_1 + k_2) [A]; \\frac{d[B]}{dt} = k_1 [A]; \\frac{d[C]}{dt} = k_2 [A]$. " +
      "Analytic: $[A](t) = A_0 e^{-(k_1+k_2)t}$, $[B]$ and $[C]$ given by integrals; numeric RK4 should match.",
  },
  {
    id: 2,
    title: "Predator–Prey linear stability",
    prompt: "Linearize around equilibrium $(\\gamma/\\delta, \\alpha/\\beta)$, compute Jacobian and eigenvalues to determine stability.",
    solution:
      "Jacobian $J = \\begin{bmatrix} \\alpha - \\beta z & -\\beta x \\\\ \\delta z & \\delta x - \\gamma \\end{bmatrix}$; " +
      "evaluated at $(\\gamma/\\delta, \\alpha/\\beta)$ yields purely imaginary eigenvalues for classic LV → center (neutrally stable).",
  },
  {
    id: 3,
    title: "RLC Q-factor",
    prompt: "Compute $\\omega_0 = 1/\\sqrt{LC}$ and $Q$. Simulate $q(t)$, extract envelope & compare numeric $Q$ to analytic.",
    solution:
      "$Q = \\frac{1}{R} \\sqrt{\\frac{L}{C}}$ (for series forms); numerically extract maxima and fit $e^{-\\omega_0 t / (2Q)}$ envelope to estimate $Q$.",
  },
  {
    id: 4,
    title: "Pendulum period amplitude dependence",
    prompt: "Compute period numerically for $\\theta_0 = 0.1, 1.0, 2.0$ rad and compare to small-angle $T \\approx 2\\pi \\sqrt{L/g}$.",
    solution:
      "Nonlinear period increases with amplitude; integrate and find time between successive crossings to estimate period. " +
      "Use energy integral or direct simulation.",
  },
  {
    id: 5,
    title: "Shooting vs finite-difference for BVP",
    prompt: "Implement shooting and finite-difference for $y'' = -\\pi^2 y + \\sin(\\pi x)$, compare solutions & errors.",
    solution:
      "Shooting converts to IVP and uses root finding on slope; finite-difference yields tridiagonal system solved by Thomas algorithm. " +
      "Compare $L_\\infty$ errors on refined meshes.",
  },
];


  return (
    <motion.div className={`${theme.panelClass} p-4 rounded-2xl`} initial={{ opacity: 0, y: 8 }} whileInView={{ opacity: 1, y: 0 }} viewport={{ once: true }}>
      <div className="flex items-center justify-between">
        <SectionHeader icon={ListChecks} number="Problems" title="Exercises & Solutions" />
        <div className="text-xs text-zinc-400">Click to toggle hints/solutions</div>
      </div>

      <div className="mt-4 space-y-3">
        {problems.map((p, i) => (
          <div key={p.id} className="rounded-lg border border-zinc-700 p-3 bg-zinc-900/50">
            <div className="flex items-center justify-between">
              <div>
                <div className="text-sm font-medium text-zinc-100">{p.title}</div>
                <MD className="text-xs text-zinc-400">{p.prompt}</MD>
              </div>
              <div>
                <button
                  className="text-xs px-2 py-1 bg-emerald-500 cursor-pointer rounded border border-zinc-700"
                  onClick={() => setOpenIndex(openIndex === i ? -1 : i)}
                >
                  {openIndex === i ? "Hide solution" : "Show solution"}
                </button>
              </div>
            </div>
            {openIndex === i && (
              <div className="mt-3 bg-zinc-950/30 border border-zinc-800 rounded p-3 text-sm text-zinc-200">
                <div className="font-semibold mb-1">Solution sketch</div>
                <MD>{p.solution}</MD>
              </div>
            )}
          </div>
        ))}
      </div>
    </motion.div>
  );
}

/* ---------------------------------------------------------------------------
   Documentation panel (notes & references) — responsive
   --------------------------------------------------------------------------- */
function ChapterDocs() {
  const [active, setActive] = useState(0);
  const sections = [
    {
      title: "28.1 Reactor Transients",
      notes: [
        "Use analytic solutions (when available) to validate integrators.",
        "Be careful with stiff kinetics — implicit solvers or adaptive step-size needed.",
        "Plot reaction rates and timescales to understand dominant dynamics.",
      ],
    },
    {
      title: "28.2 Predator–Prey",
      notes: [
        "Lotka–Volterra has neutral cycles; adding logistic or harvesting terms changes stability.",
        "Phase-plane and Poincaré sections help visualize periodic orbits and chaos.",
        "Chaos emerges in higher-dimensional systems (Lorenz), showing sensitive dependence on initial conditions.",
      ],
    },
    {
      title: "28.3 RLC Circuits",
      notes: [
        "Parameter sweeps (R, L, C) show underdamped/overdamped regimes and resonance peaks.",
        "Use analytic formulas where possible; they make excellent checks for numerical solvers.",
        "Estimate Q factor by fitting exponential envelope to maxima.",
      ],
    },
    {
      title: "28.4 Pendulum",
      notes: [
        "Nonlinear pendulum period > small-angle period for large amplitude.",
        "Energy conservation indicates integrator accuracy; poor integrators introduce numerical damping.",
        "Phase space separatrix separates oscillatory and rotational motion.",
      ],
    },
  ];

  return (
    <Card className={`${theme.panelClass} rounded-2xl`}>
      <CardContent className="p-4">
        <div className="flex items-center justify-between mb-3">
          <div className="flex items-center gap-2">
            <BookOpen className="w-5 h-5 text-cyan-400" />
            <div>
              <div className="text-sm text-zinc-100 font-semibold">Chapter 28 — Notes & References</div>
              <div className="text-xs text-zinc-400">Modelling tips and recommended readings</div>
            </div>
          </div>
          <div className="text-xs text-zinc-400">Strogatz, Ogata, Numerical Recipes</div>
        </div>

        <div className="grid md:grid-cols-4 gap-3">
          <div className="space-y-2">
            {sections.map((s, i) => (
              <button
                key={i}
                onClick={() => setActive(i)}
                className={`w-full text-left px-3 py-2 text-gray-200 cursor-pointer rounded border ${i === active ? "bg-cyan-900/30 border-cyan-600" : "bg-zinc-900/40 border-zinc-700"}`}
              >
                {s.title}
              </button>
            ))}
          </div>

          <div className="md:col-span-3 bg-zinc-950/30 border border-zinc-800 rounded p-4">
            <div className="text-sm font-semibold text-zinc-100 mb-2">{sections[active].title}</div>
            <ul className="text-sm text-zinc-300 list-disc ml-5 space-y-1">
              {sections[active].notes.map((n, idx) => (
                <li key={idx}>{n}</li>
              ))}
            </ul>
            <div className="text-xs text-zinc-500 mt-3">References: Strogatz (Nonlinear Dynamics), Ogata (Control Systems), Press et al. (Numerical Recipes)</div>
          </div>
        </div>
      </CardContent>
    </Card>
  );
}

/* ---------------------------------------------------------------------------
   Chapter assembly — bring all sections together in a responsive page
   --------------------------------------------------------------------------- */
export default function Chapter28() {
  useEffect(() => {
    // keep UX nice: scroll top when opening chapter
    if (typeof window !== "undefined") window.scrollTo({ top: 0, behavior: "smooth" });
  }, []);

  return (
    <div className={`p-4 sm:p-6 lg:p-8 space-y-6 ${theme.bgClass} min-h-screen`}>
      <motion.div initial={{ opacity: 0, y: -8 }} animate={{ opacity: 1, y: 0 }}>
        <div className="flex flex-col sm:flex-row sm:items-center sm:justify-between gap-4">
          <div>
        
            <h1 className="text-2xl sm:text-3xl font-bold text-emerald-400">Case Studies: Ordinary Differential Equations</h1>
            <div className="text-zinc-400 mt-1">28.1 Reactor • 28.2 Predator–Prey • 28.3 RLC • 28.4 Pendulum • Problems</div>
          </div>
        </div>
      </motion.div>

      <div className="space-y-6">
        <ReactorSection />
        <PredatorPreySection />
        <RLCCircuitSection />
        <PendulumSection />
        <ProblemsSection />
        <ChapterDocs />
        <BottomBar/>
      </div>

      <footer className="text-xs text-zinc-500 mt-6">
        Generated chapter UI — responsive, avoids react-markdown `className` issue, and uses responsive Plotly containers.
      </footer>
    </div>
  );
}
