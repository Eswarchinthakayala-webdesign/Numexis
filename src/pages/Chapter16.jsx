// Part 1 of 2
// src/pages/Chapter16.jsx
// Expanded Chapter 16 — Case Studies: Optimization
// Sections 16.1 → 16.4. Includes long docs, 2D & 3D plots, multiple exercises.
// Paste Part 1, then Part 2 immediately after into src/pages/Chapter16.jsx

import React, { useMemo, useState, useRef } from "react";
import Plot from "react-plotly.js";
import { motion } from "framer-motion";
import ReactMarkdown from "react-markdown";
import remarkMath from "remark-math";
import rehypeKatex from "rehype-katex";
import "katex/dist/katex.min.css";

import { FlaskConical, Droplets, Zap, Atom, BookOpen, Download } from "lucide-react";
import BottomBar from "../components/BottomBar";

// Attempt to reuse existing UI components; fallback to simple elements
let Input, Button, Card;
try {
  // eslint-disable-next-line
  Input = require("@/components/ui/input").Input;
  // eslint-disable-next-line
  Button = require("@/components/ui/button").Button;
  // eslint-disable-next-line
  Card = require("@/components/ui/card").Card;
} catch (e) {
  Input = (props) => <input {...props} />;
  Button = (props) => <button {...props} />;
  Card = ({ children, className, style }) => <div className={className} style={style}>{children}</div>;
}

// Theme
const THEME = {
  bg: "bg-zinc-900",
  panel: "#1b1f22",
  text: "#e6eef6",
  muted: "#9ca3af",
  accent: "#15803d",
  accent2: "#b45309",
};

// small animation preset
const fadeUp = { initial: { opacity: 0, y: 8 }, animate: { opacity: 1, y: 0 }, transition: { duration: 0.28 } };

// Utility: safe JS function builder for single/2-variable expressions
function makeFun(expr, vars = ["x"]) {
  try {
    const varList = vars.join(",");
    // Provide common Math functions
    const wrapper = `
      const { sin, cos, tan, asin, acos, atan, exp, log, sqrt, pow, abs, PI, E } = Math;
      return (${varList}) => { return (${expr}); };
    `;
    // eslint-disable-next-line no-new-func
    return new Function(wrapper)();
  } catch (e) {
    return null;
  }
}

// Small helpers
function linspace(a, b, n) {
  const out = [];
  for (let i = 0; i < n; i++) out.push(a + (b - a) * (i / (n - 1)));
  return out;
}
function formatNum(x) {
  if (Array.isArray(x)) return `[${x.map((v) => Number(v).toPrecision(6)).join(", ")}]`;
  return Number(x).toPrecision(6);
}
function toFixed(x, p = 4) {
  return Number(x).toFixed(p);
}

// Exercise component (built-in)
function Exercise({ title, body, solution, hint }) {
  const [open, setOpen] = useState(false);
  const [showSol, setShowSol] = useState(false);
  const [showHint, setShowHint] = useState(false);
  return (
    <div className="rounded-lg border border-zinc-700 bg-zinc-800 p-4">
      <div className="flex items-start justify-between gap-3">
        <div>
          <div style={{ color: THEME.text, fontWeight: 700 }}>{title}</div>
          <div className="text-sm" style={{ color: THEME.muted, marginTop: 6 }}>
            <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>{body}</ReactMarkdown>
          </div>
        </div>
        <div className="flex flex-col gap-2">
          <Button onClick={() => setShowHint(s => !s)} className="rounded  bg-white text-black cursor-pointer hover:bg-gray-200">Hint</Button>
          <Button onClick={() => setShowSol(s => !s)} className="px-2 rounded  bg-white text-black cursor-pointer hover:bg-gray-200">Solution</Button>
        </div>
      </div>
      {showHint && hint && <div className="mt-3 text-sm text-amber-300"><strong>Hint:</strong> <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>{hint}</ReactMarkdown></div>}
      {showSol && solution && <div className="mt-3 text-sm text-zinc-200 border-t border-zinc-700 pt-3"><ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>{solution}</ReactMarkdown></div>}
    </div>
  );
}

// Section header
function SectionHeader({ Icon, title, subtitle }) {
  return (
    <div className="flex items-center gap-3 mb-4">
      <div style={{ width: 44, height: 44, borderRadius: 8, background: "#071018", display: "flex", alignItems: "center", justifyContent: "center", border: "1px solid rgba(255,255,255,0.03)" }}>
        <Icon className="w-5 h-5" />
      </div>
      <div>
        <div style={{ color: THEME.text, fontSize: 18, fontWeight: 700 }}>{title}</div>
        <div style={{ color: THEME.muted, fontSize: 13 }}>{subtitle}</div>
      </div>
    </div>
  );
}

// Long docs used across sections (KaTeX strings)
const longDocs = {
  "16.1": [
    `**Problem statement.** Minimize the material cost of a closed cylindrical tank (radius $r$, height $h$) subject to a fixed interior volume $V$.`,
    `Volume constraint: $$V = \pi r^2 h.$$`,
    `Surface area (material) to minimize: $$A = 2\pi r h + 2\pi r^2.$$`,
    `Using Lagrange multipliers: set up $\mathcal{L}(r,h,\lambda)=2\pi r h + 2\pi r^2 - \lambda(\pi r^2 h - V)$.`,
    `Stationarity yields (after algebra) the well-known result for cost-minimizing cylinder: $$h = 2r.$$`,
    `This can be seen numerically by scanning $r$ and computing $h=V/(\pi r^2)$ and the resulting area.`,
    `The 3D surface of area vs $(r,h)$ shows the valley along feasible manifold $\pi r^2 h=V$.`
  ],
  "16.2": [
    `**Case idea.** Treatment fraction $t\in[0,1]$ influences treatment cost and discharge penalties. Define a sample cost model: $$C(t)=c_1 t Q + c_2(1-t)Q + \alpha t^2,$$ where $Q$ is flow.`,
    `Optimize $t$ for minimal cost — trivial 1D problem but useful to examine sensitivity and robustness to parameters.`
  ],
  "16.3": [
    `**Maximum power transfer theorem.** For a Thevenin equivalent source of voltage $V_s$ with internal resistance $R_s$ feeding a load $R_L$, the load power is $$P_L = \dfrac{V_s^2 R_L}{(R_s+R_L)^2}.$$`,
    `Differentiating with respect to $R_L$ and setting derivative to zero yields $R_L = R_s$ for maximum power.`
  ],
  "16.4": [
    `**Minimum potential energy principle.** For conservative systems, equilibrium configurations correspond to stationary points of potential energy. For a simple single-degree-of-freedom model with potential $U(x)$, equilibrium satisfies $U'(x)=0$ and stable equilibrium requires $U''(x)>0$.`,
    `We illustrate with $U(x)=\\ k x^2 + \\frac14 \\beta x^4$.`
  ]
};


function DocsPanel() {
  const [open, setOpen] = useState(Object.keys(docs10)[0]);
  return (
    <div className="w-full mt-6 bg-zinc-900 border border-zinc-700 rounded-xl overflow-hidden">
      <div className="p-3 border-b border-zinc-700 flex items-center justify-between">
        <div className="flex items-center gap-3">
          <BookOpen className="w-5 h-5 text-cyan-400" />
          <div>
            <div className="text-zinc-100 font-semibold">Chapter 10 — Documentation</div>
            <div className="text-zinc-400 text-xs">Theory, algorithms, and numerical tips</div>
          </div>
        </div>
       
      </div>

      <div className="grid grid-cols-1 md:grid-cols-4">
        <div className="col-span-1 border-r border-zinc-800">
          {Object.keys(docs10).map((k) => (
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
          <div className="space-y-3 text-sm text-zinc-300">
            {docs10[open].map((p, i) => (
              <div key={i}>
                <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>
                  {p}
                </ReactMarkdown>
              </div>
            ))}
          </div>
        </div>
      </div>
    </div>
  );
}


// ----------------------------- Section 16.1: Tank Design -----------------------------
function Section161() {
  const [V, setV] = useState(200); // default volume
  const [rMin, setRMin] = useState(0.1);
  const [rMax, setRMax] = useState(10);
  const [rSteps, setRSteps] = useState(600); // dense sampling for smooth curve

  // sample r -> compute h from volume -> area
  const samples = useMemo(() => {
    const xs = linspace(Number(rMin), Number(rMax), Number(rSteps));
    const data = xs.map((r) => {
      const h = (Number(V) / (Math.PI * r * r));
      const A = 2 * Math.PI * r * h + 2 * Math.PI * r * r;
      return { r, h, A };
    });
    return data.filter(d => isFinite(d.A) && d.r > 0 && d.h > 0);
  }, [V, rMin, rMax, rSteps]);

  const best = useMemo(() => {
    if (!samples.length) return null;
    return samples.reduce((p, c) => (c.A < p.A ? c : p), samples[0]);
  }, [samples]);

  // 3D surface of A(r,h) over grid (for visualization)
  const surface = useMemo(() => {
    const rGrid = linspace(Number(rMin), Number(rMax), 80);
    const hGrid = linspace(0.1, Math.max(1, best ? best.h * 1.8 : 10), 80);
    const Z = [];
    for (let i = 0; i < hGrid.length; i++) {
      const row = [];
      for (let j = 0; j < rGrid.length; j++) {
        const r = rGrid[j];
        const h = hGrid[i];
        const A = 2 * Math.PI * r * h + 2 * Math.PI * r * r;
        row.push(A);
      }
      Z.push(row);
    }
    return { rGrid, hGrid, Z };
  }, [rMin, rMax, best]);

  return (
    <section className="rounded-lg p-4" style={{ background: THEME.panel }}>
      <motion.div {...fadeUp}>
        <SectionHeader Icon={FlaskConical} title="16.1 Least-Cost Design of a Tank" subtitle="Chemical/Biological Engineering case study" />

        <div className="text-sm text-zinc-200 mb-3">
          <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>{longDocs["16.1"].join("\n\n")}</ReactMarkdown>
        </div>

        <div className="rounded-lg bg-zinc-800 border border-zinc-700 p-3 mb-4">
          <div className="grid grid-cols-1 md:grid-cols-3 gap-3">
            <div>
              <label className="text-xs text-zinc-400">Volume V</label>
              <Input value={V} onChange={(e) => setV(e.target.value)} className="w-full bg-zinc-900 p-2 rounded text-zinc-100" />
              <div className="text-xs text-zinc-400 mt-2">Set the required internal volume (same units used for radius/height calculations).</div>
            </div>

            <div>
              <label className="text-xs text-zinc-400">r range</label>
              <div className="flex gap-2">
                <Input value={rMin} onChange={(e) => setRMin(e.target.value)} className="bg-zinc-900 p-2 rounded text-zinc-100" />
                <Input value={rMax} onChange={(e) => setRMax(e.target.value)} className="bg-zinc-900 p-2 rounded text-zinc-100" />
              </div>
            </div>

            <div>
              <label className="text-xs text-zinc-400">sampling points</label>
              <Input value={rSteps} onChange={(e) => setRSteps(Number(e.target.value) || 100)} className="bg-zinc-900 p-2 rounded text-zinc-100" />
              <div className="text-xs text-zinc-400 mt-1">Use denser sampling for smoother minima estimation.</div>
            </div>
          </div>
        </div>

        {/* 2D plot: area vs r */}
        <div className="rounded-lg bg-zinc-900 border border-zinc-700 p-3 mb-4">
          <div className="text-xs text-zinc-400 mb-2">Surface area vs radius (with h determined by volume)</div>
          <Plot
            data={[
              { x: samples.map(s => s.r), y: samples.map(s => s.A), mode: "lines", name: "A(r)" },
              ...(best ? [{ x: [best.r], y: [best.A], mode: "markers", name: "optimal", marker: { size: 10, color: THEME.accent2 } }] : [])
            ]}
            layout={{
              autosize: true,
              height: 360,
              margin: { l: 50, r: 10, t: 30, b: 40 },
              paper_bgcolor: "rgba(0,0,0,0)",
              plot_bgcolor: "rgba(0,0,0,0)"
            }}
            useResizeHandler
            style={{ width: "100%" }}
          />
        </div>

        {/* 3D surface of A(r,h) */}
        <div className="rounded-lg bg-zinc-900 border border-zinc-700 p-3 mb-4">
          <div className="text-xs text-zinc-400 mb-2">3D Surface: A(r,h)</div>
          <Plot
            data={[
              { x: surface.rGrid, y: surface.hGrid, z: surface.Z, type: "surface", contours: { z: { show: true } }, showscale: false }
            ]}
            layout={{ autosize: true, height: 420, scene: { xaxis: { title: "r" }, yaxis: { title: "h" }, zaxis: { title: "A" } }, paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)" }}
            useResizeHandler
            style={{ width: "100%" }}
          />
        </div>

        {/* results */}
        <div className="rounded-lg bg-zinc-800 border border-zinc-700 p-3 mb-4 grid grid-cols-1 md:grid-cols-3 gap-3">
          <div className="p-2 bg-zinc-900 rounded border border-zinc-700">
            <div className="text-xs text-zinc-400">Optimal radius r*</div>
            <div className="text-sm text-zinc-100">{best ? toFixed(best.r, 4) : "—"}</div>
          </div>
          <div className="p-2 bg-zinc-900 rounded border border-zinc-700">
            <div className="text-xs text-zinc-400">Corresponding height h*</div>
            <div className="text-sm text-zinc-100">{best ? toFixed(best.h, 4) : "—"}</div>
          </div>
          <div className="p-2 bg-zinc-900 rounded border border-zinc-700">
            <div className="text-xs text-zinc-400">Minimal surface area</div>
            <div className="text-sm text-zinc-100">{best ? toFixed(best.A, 6) : "—"}</div>
          </div>
        </div>

        {/* multiple exercises */}
        <div className="space-y-3">
          <Exercise
            title="Exercise 16.1.1 — Derive the analytic result"
            body={`Use Lagrange multipliers to show that the minimal surface area occurs when $$h = 2r$$ for the cylindrical closed tank with volume \\(V\\).`}
            hint={`Write the Lagrangian \\(\\mathcal{L}(r,h,\\lambda) = 2\\pi r h + 2\\pi r^2 - \\lambda(\\pi r^2 h - V)\\) and take partial derivatives.`}
            solution={`From \\(\\partial_r \\mathcal{L}=0\\) and \\(\\partial_h \\mathcal{L}=0\\) we obtain equations that reduce to \\(h=2r\\) when eliminating the multiplier. Substitute back into the volume constraint to find r in terms of V.`}
          />

          <Exercise
            title="Exercise 16.1.2 — Sensitivity analysis"
            body={`For the optimal design, how sensitive is the minimal area to a 10% change in required volume \\(V\\)? Compute numerically.`}
            hint={`Compute optimal r for V and 1.1V and compare A(V) vs A(1.1V).`}
            solution={`Numerically, surface area scales roughly as V^{2/3}; a 10% increase in V results in approximately a 6.7% increase in minimal area (since (1.1)^{2/3} ≈ 1.067).`}
          />

          <Exercise
            title="Exercise 16.1.3 — Cost vs material thickness"
            body={`If material cost scales with surface area and thickness, and thickness increases linearly with radius due to manufacturing constraints t= t0 + k r, how does the optimal r change qualitatively?`}
            hint={`Total cost ~ A(r) * (t0 + k r). Consider derivative or numerically plot.`}
            solution={`Including thickness that increases with r pushes the optimum toward smaller radii (tends to flatter, taller shapes) — analyze by numerically minimizing over r.`}
          />
        </div>
      </motion.div>
    </section>
  );
}

// ----------------------------- Section 16.2: Wastewater Treatment -----------------------------
function Section162() {
  // parameters
  const [Q, setQ] = useState(100); // flow
  const [c1, setC1] = useState(25); // per-unit treatment cost
  const [c2, setC2] = useState(8); // per-unit discharge cost
  const [alpha, setAlpha] = useState(120); // penalty coefficient (nonlinear)
  const [tSteps, setTSteps] = useState(501);

  // cost model C(t) = c1 t Q + c2 (1-t) Q + alpha * t^2
  const ts = useMemo(() => linspace(0, 1, Number(tSteps)), [tSteps]);
  const costs = useMemo(() => ts.map(t => (Number(c1) * t * Number(Q) + Number(c2) * (1 - t) * Number(Q) + Number(alpha) * t * t)), [ts, Q, c1, c2, alpha]);
  const bestIdx = useMemo(() => costs.indexOf(Math.min(...costs)), [costs]);
  const bestT = ts[bestIdx];
  const bestC = costs[bestIdx];

  // 3D surface: vary Q and t
  const qGrid = useMemo(() => linspace(Math.max(1, Number(Q) * 0.2), Number(Q) * 2, 60), [Q]);
  const tGrid = useMemo(() => linspace(0, 1, 60), []);
  const costSurface = useMemo(() => {
    const Z = [];
    for (let i = 0; i < tGrid.length; i++) {
      const row = [];
      for (let j = 0; j < qGrid.length; j++) {
        const qv = qGrid[j];
        const tv = tGrid[i];
        const cv = Number(c1) * tv * qv + Number(c2) * (1 - tv) * qv + Number(alpha) * tv * tv;
        row.push(cv);
      }
      Z.push(row);
    }
    return { qGrid, tGrid, Z };
  }, [qGrid, tGrid, c1, c2, alpha]);

  return (
    <section className="rounded-lg p-4" style={{ background: THEME.panel }}>
      <motion.div {...fadeUp}>
        <SectionHeader Icon={Droplets} title="16.2 Least-Cost Treatment of Wastewater" subtitle="Civil / Environmental Engineering case study" />

        <div className="text-sm text-zinc-200 mb-3">
          <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>{longDocs["16.2"].join("\n\n")}</ReactMarkdown>
        </div>

        <div className="rounded-lg bg-zinc-800 border border-zinc-700 p-3 mb-4">
          <div className="grid grid-cols-1 md:grid-cols-4 gap-3">
            <div>
              <label className="text-xs text-zinc-400">Flow Q</label>
              <Input value={Q} onChange={(e) => setQ(e.target.value)} className="bg-zinc-900 p-2 rounded text-zinc-100" />
            </div>
            <div>
              <label className="text-xs text-zinc-400">Treatment unit cost c1</label>
              <Input value={c1} onChange={(e) => setC1(e.target.value)} className="bg-zinc-900 p-2 rounded text-zinc-100" />
            </div>
            <div>
              <label className="text-xs text-zinc-400">Discharge cost c2</label>
              <Input value={c2} onChange={(e) => setC2(e.target.value)} className="bg-zinc-900 p-2 rounded text-zinc-100" />
            </div>
            <div>
              <label className="text-xs text-zinc-400">Penalty α</label>
              <Input value={alpha} onChange={(e) => setAlpha(e.target.value)} className="bg-zinc-900 p-2 rounded text-zinc-100" />
            </div>
          </div>
        </div>

        {/* 2D plot */}
        <div className="rounded-lg bg-zinc-900 border border-zinc-700 p-3 mb-4">
          <div className="text-xs text-zinc-400 mb-2">Cost(t) for fixed Q</div>
          <Plot
            data={[
              { x: ts, y: costs, mode: "lines", name: "C(t)" },
              { x: [bestT], y: [bestC], mode: "markers", marker: { size: 10, color: THEME.accent2 }, name: "optimal" }
            ]}
            layout={{ autosize: true, height: 360, margin: { l: 50, r: 10, t: 30, b: 40 }, paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)" }}
            useResizeHandler
            style={{ width: "100%" }}
          />
        </div>

        {/* 3D surface */}
        <div className="rounded-lg bg-zinc-900 border border-zinc-700 p-3 mb-4">
          <div className="text-xs text-zinc-400 mb-2">Cost surface over (Q, t)</div>
          <Plot
            data={[{ x: costSurface.qGrid, y: costSurface.tGrid, z: costSurface.Z, type: "surface", contours: { z: { show: true } } }]}
            layout={{ autosize: true, height: 420, scene: { xaxis: { title: "Q" }, yaxis: { title: "t" }, zaxis: { title: "Cost" } }, paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)" }}
            useResizeHandler
            style={{ width: "100%" }}
          />
        </div>

        <div className="rounded-lg bg-zinc-800 border border-zinc-700 p-3 mb-4 grid grid-cols-1 md:grid-cols-3 gap-3">
          <div className="p-2 bg-zinc-900 rounded border border-zinc-700">
            <div className="text-xs text-zinc-400">Optimal treatment fraction t*</div>
            <div className="text-sm text-zinc-100">{toFixed(bestT, 4)}</div>
          </div>
          <div className="p-2 bg-zinc-900 rounded border border-zinc-700">
            <div className="text-xs text-zinc-400">Cost at optimum</div>
            <div className="text-sm text-zinc-100">{toFixed(bestC, 4)}</div>
          </div>
          <div className="p-2 bg-zinc-900 rounded border border-zinc-700">
            <div className="text-xs text-zinc-400">Sensitivity</div>
            <div className="text-sm text-zinc-100">Try changing α, c1, c2 to observe shifts</div>
          </div>
        </div>

        <div className="space-y-3">
          <Exercise
            title="Exercise 16.2.1 — Parameter study"
            body={`Perform a parameter sweep over \\(\\alpha\\) in [0,500] and plot t*(\\alpha). What trend do you observe?`}
            hint={`For small α the optimal t tends to be small (cheap to discharge); for large α optimal t increases.`}
            solution={`t*(α) increases with α; beyond some threshold it approaches 1 (full treatment).`}
          />

          <Exercise
            title="Exercise 16.2.2 — Flow scaling"
            body={`For fixed α and cost coefficients, how does optimal t scale with Q? Use the 3D surface to reason.`}
            hint={`Since linear terms in Q scale cost linearly, but penalty α t^2 does not depend on Q, increasing Q tends to shift balance towards lower or higher t depending on coefficients.`}
            solution={`If treatment unit cost is high relative to discharge cost, larger Q amplifies treatment cost making t* smaller; visualize on the surface for clarity.`}
          />

          <Exercise
            title="Exercise 16.2.3 — Discrete investments"
            body={`Suppose treatment comes in discrete increments (0, 0.5, 1.0). Which is cost-optimal? Model and compare.`}
            hint={`Evaluate C(t) at allowed discrete t and choose minimal cost.`}
            solution={`Choose the discrete t that yields minimal total cost; often discrete constraints make intermediate values infeasible so selection is by evaluating each.`}
          />
        </div>
      </motion.div>
    </section>
  );
}

// ----------------------------- Section 16.3: Max Power Transfer -----------------------------
function Section163() {
  // parameters and sampling
  const [Vs, setVs] = useState(10);
  const [RLfixed, setRLfixed] = useState(5); // default load
  const RsSamples = useMemo(() => linspace(0.1, 20, 400), []);
  const PL = useMemo(() => RsSamples.map(Rs => (Vs * Vs * RLfixed) / ((Rs + RLfixed) * (Rs + RLfixed))), [RsSamples, Vs, RLfixed]);
  const RsOptIdx = useMemo(() => PL.indexOf(Math.max(...PL)), [PL]);
  const RsOpt = RsSamples[RsOptIdx];

  // 3D surface: vary Rs and RL -> power delivered to RL
  const RLgrid = useMemo(() => linspace(0.1, 20, 60), []);
  const Rsgrid = useMemo(() => linspace(0.1, 20, 60), []);
  const powerSurface = useMemo(() => {
    const Z = [];
    for (let i = 0; i < RLgrid.length; i++) {
      const row = [];
      for (let j = 0; j < Rsgrid.length; j++) {
        const RL = RLgrid[i];
        const Rs = Rsgrid[j];
        const P = (Vs * Vs * RL) / ((Rs + RL) * (Rs + RL));
        row.push(P);
      }
      Z.push(row);
    }
    return { Rsgrid, RLgrid, Z };
  }, [Rsgrid, RLgrid, Vs]);

  return (
    <section className="rounded-lg p-4" style={{ background: THEME.panel }}>
      <motion.div {...fadeUp}>
        <SectionHeader Icon={Zap} title="16.3 Maximum Power Transfer for a Circuit" subtitle="Electrical Engineering case study" />

        <div className="text-sm text-zinc-200 mb-3">
          <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>{longDocs["16.3"].join("\n\n")}</ReactMarkdown>
        </div>

        <div className="rounded-lg bg-zinc-800 border border-zinc-700 p-3 mb-4 grid grid-cols-1 md:grid-cols-3 gap-3">
          <div>
            <label className="text-xs text-zinc-400">Source voltage Vs</label>
            <Input value={Vs} onChange={(e) => setVs(Number(e.target.value))} className="bg-zinc-900 p-2 rounded text-zinc-100" />
          </div>
          <div>
            <label className="text-xs text-zinc-400">Load RL (for 2D slice)</label>
            <Input value={RLfixed} onChange={(e) => setRLfixed(Number(e.target.value))} className="bg-zinc-900 p-2 rounded text-zinc-100" />
          </div>
          <div>
            <div className="text-xs text-zinc-400">Optimal Rs for RL={RLfixed}</div>
            <div className="text-sm text-zinc-100">{toFixed(RsOpt, 4)}</div>
          </div>
        </div>

        {/* 2D plot: P vs Rs */}
        <div className="rounded-lg bg-zinc-900 border border-zinc-700 p-3 mb-4">
          <div className="text-xs text-zinc-400 mb-2">Power delivered to RL vs source resistance Rs (RL fixed)</div>
          <Plot
            data={[
              { x: RsSamples, y: PL, mode: "lines", name: "P(Rs)" },
              { x: [RsOpt], y: [PL[RsOptIdx]], mode: "markers", name: "optimal Rs", marker: { size: 10, color: THEME.accent2 } }
            ]}
            layout={{ autosize: true, height: 360, margin: { l: 50, r: 10, t: 30, b: 40 }, paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)" }}
            useResizeHandler
            style={{ width: "100%" }}
          />
        </div>

        {/* 3D surface vs Rs, RL */}
        <div className="rounded-lg bg-zinc-900 border border-zinc-700 p-3 mb-4">
          <div className="text-xs text-zinc-400 mb-2">Power surface over (Rs, RL)</div>
          <Plot
            data={[{ x: powerSurface.Rsgrid, y: powerSurface.RLgrid, z: powerSurface.Z, type: "surface", contours: { z: { show: true } } }]}
            layout={{ autosize: true, height: 420, scene: { xaxis: { title: "Rs" }, yaxis: { title: "RL" }, zaxis: { title: "P" } }, paper_bgcolor: "rgba(0,0,0,0)" }}
            useResizeHandler
            style={{ width: "100%" }}
          />
        </div>

        <div className="space-y-3">
          <Exercise
            title="Exercise 16.3.1 — Analytical derivation"
            body={`Differentiate $$P(R_L) = \\dfrac{V_s^2 R_L}{(R_s+R_L)^2}$$ with respect to \\(R_L\\) and show that the maximizer is \\(R_L = R_s\\).`}
            hint="Use quotient rule or rewrite as $$P(R_L)=V_s^2 \\dfrac{R_L}{(R_s+R_L)^2}$$ and set derivative to zero."
            solution={`Differentiating and simplifying gives (R_s - R_L) in numerator, hence R_L = R_s at optimum. Confirm second derivative negative. `}
          />

          <Exercise
            title="Exercise 16.3.2 — Efficiency at max power"
            body={`Compute efficiency η = P_L / P_generated at R_L = R_s. What fraction of generated power is delivered to the load vs dissipated internally?`}
            hint="At R_L=R_s, check voltage division and compute P_generated = V_s^2 / (4 R_s)."
            solution={`At R_L=R_s delivered power P_L = V_s^2/(4 R_s). The total generated power equals the same value, and system dissipates equal amount? In fact, at max power half the available power is dissipated in R_s; efficiency is 50%.`}
          />

          <Exercise
            title="Exercise 16.3.3 — Nonlinear loads"
            body={`If the load depends on frequency (complex impedance), how does the optimum concept extend? Consider matching magnitude of impedance for maximum transfer.`}
            hint={`For complex loads, maximum transfer (in amplitude sense) requires conjugate matching: Z_L = Z_s^*.`}
            solution={`Maximum average power transfer requires the load impedance to be the complex conjugate of source impedance; this ensures resistive part matches and reactive parts cancel.`}
          />
        </div>
      </motion.div>
    </section>
  );
}

// Part 2 of 2
// Continue Chapter16.jsx (paste immediately after Part 1)

// ----------------------------- Section 16.4: Minimum Potential Energy -----------------------------
function Section164() {
  // potential U(x) = 0.5 k x^2 + 0.25 beta x^4 (quartic stiffening)
  const [k, setK] = useState(1.0);
  const [beta, setBeta] = useState(0.2);
  const xs = useMemo(() => linspace(-3, 3, 601), []);
  const U = useMemo(() => xs.map(x => 0.5 * k * x * x + 0.25 * beta * Math.pow(x, 4)), [xs, k, beta]);
  const minIdx = useMemo(() => U.indexOf(Math.min(...U)), [U]);
  const minX = xs[minIdx];
  const minU = U[minIdx];

  // 3D surface for a two-degree-of-freedom mockup: U(x,y) = 0.5 k1 x^2 + 0.5 k2 y^2 + coupling c xy + beta1 x^4 + beta2 y^4
  const [k2, setK2] = useState(1.5);
  const [c, setC] = useState(0.4);
  const [beta1, setBeta1] = useState(0.05);
  const [beta2, setBeta2] = useState(0.05);
  const grid = useMemo(() => {
    const X = linspace(-2, 2, 80);
    const Y = linspace(-2, 2, 80);
    const Z = [];
    for (let i = 0; i < Y.length; i++) {
      const row = [];
      for (let j = 0; j < X.length; j++) {
        const xv = X[j], yv = Y[i];
        const val = 0.5 * k * xv * xv + 0.5 * k2 * yv * yv + c * xv * yv + 0.25 * beta1 * Math.pow(xv, 4) + 0.25 * beta2 * Math.pow(yv, 4);
        row.push(val);
      }
      Z.push(row);
    }
    return { X, Y, Z };
  }, [k, k2, c, beta1, beta2]);

  // find stationary points of 1D potential via numerical derivative
  function numericDerivative(arr, dx) {
    const n = arr.length;
    const d = new Array(n).fill(0);
    for (let i = 1; i < n - 1; i++) d[i] = (arr[i + 1] - arr[i - 1]) / (2 * dx);
    d[0] = (arr[1] - arr[0]) / dx;
    d[n - 1] = (arr[n - 1] - arr[n - 2]) / dx;
    return d;
  }
  const dU = useMemo(() => numericDerivative(U, xs[1] - xs[0]), [U, xs]);

  // find roots where sign change of derivative -> stationary points
  const stationaryIdx = useMemo(() => {
    const idxs = [];
    for (let i = 1; i < dU.length; i++) {
      if (dU[i - 1] === 0 || dU[i] === 0) idxs.push(i);
      else if (dU[i - 1] * dU[i] < 0) idxs.push(i);
    }
    return idxs;
  }, [dU]);

  return (
    <section className="rounded-lg p-4" style={{ background: THEME.panel }}>
      <motion.div {...fadeUp}>
        <SectionHeader Icon={Atom} title="16.4 Equilibrium & Minimum Potential Energy" subtitle="Mechanical / Aerospace Engineering case study" />

        <div className="text-sm text-zinc-200 mb-3">
          <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>
            {longDocs["16.4"].join("\n\n")}
          </ReactMarkdown>
        </div>

        <div className="rounded bg-zinc-800 border border-zinc-700 p-3 mb-4">
          <div className="grid grid-cols-1 md:grid-cols-4 gap-3">
            <div>
              <label className="text-xs text-zinc-400">k (linear stiffness)</label>
              <Input value={k} onChange={(e) => setK(Number(e.target.value))} className="bg-zinc-900 p-2 rounded text-zinc-100" />
            </div>
            <div>
              <label className="text-xs text-zinc-400">beta (nonlinear)</label>
              <Input value={beta} onChange={(e) => setBeta(Number(e.target.value))} className="bg-zinc-900 p-2 rounded text-zinc-100" />
            </div>
            <div>
              <label className="text-xs text-zinc-400">coupling c (2DOF)</label>
              <Input value={c} onChange={(e) => setC(Number(e.target.value))} className="bg-zinc-900 p-2 rounded text-zinc-100" />
            </div>
            <div>
              <label className="text-xs text-zinc-400">k2 (stiffness y)</label>
              <Input value={k2} onChange={(e) => setK2(Number(e.target.value))} className="bg-zinc-900 p-2 rounded text-zinc-100" />
            </div>
          </div>
        </div>

        {/* 2D potential plot */}
        <div className="rounded-lg bg-zinc-900 border border-zinc-700 p-3 mb-4">
          <div className="text-xs text-zinc-400 mb-2">Potential U(x) vs x</div>
          <Plot
            data={[
              { x: xs, y: U, mode: "lines", name: "U(x)" },
              { x: [minX], y: [minU], mode: "markers", marker: { size: 10, color: THEME.accent2 }, name: "minimum" },
              ...(stationaryIdx.map(i => ({ x: [xs[i]], y: [U[i]], mode: "markers", marker: { size: 6, color: "#f97316" }, name: `stat ${i}` })))
            ]}
            layout={{ autosize: true, height: 360, margin: { l: 50, r: 10, t: 30, b: 40 }, paper_bgcolor: "rgba(0,0,0,0)",plot_bgcolor: "rgba(0,0,0,0)" }}
            useResizeHandler
            style={{ width: "100%" }}
          />
        </div>

        {/* 3D potential surface for 2DOF */}
        <div className="rounded-lg bg-zinc-900 border border-zinc-700 p-3 mb-4">
          <div className="text-xs text-zinc-400 mb-2">2D potential surface U(x,y)</div>
          <Plot
            data={[
              { x: grid.X, y: grid.Y, z: grid.Z, type: "surface", contours: { z: { show: true } }, showscale: false }
            ]}
            layout={{ autosize: true, height: 420, scene: { xaxis: { title: "x" }, yaxis: { title: "y" }, zaxis: { title: "U" } }, paper_bgcolor: "rgba(0,0,0,0)" }}
            useResizeHandler
            style={{ width: "100%" }}
          />
        </div>

        <div className="grid grid-cols-1 md:grid-cols-3 gap-3 mb-4">
          <div className="p-2 bg-zinc-800 rounded border border-zinc-700">
            <div className="text-xs text-zinc-400">Global minimum x*</div>
            <div className="text-sm text-zinc-100">{toFixed(minX, 6)}</div>
          </div>
          <div className="p-2 bg-zinc-800 rounded border border-zinc-700">
            <div className="text-xs text-zinc-400">U(x*)</div>
            <div className="text-sm text-zinc-100">{toFixed(minU, 6)}</div>
          </div>
          <div className="p-2 bg-zinc-800 rounded border border-zinc-700">
            <div className="text-xs text-zinc-400">Stationary points found</div>
            <div className="text-sm text-zinc-100">{stationaryIdx.length}</div>
          </div>
        </div>

        <div className="space-y-3">
          <Exercise
            title="Exercise 16.4.1 — Stability from second derivative"
            body={`Show that for 1D potential U(x), stability at stationary point x* requires U''(x*)>0. For the quartic model compute U''(x)=k+3β x^2 and evaluate at min.`}
            hint={`Differentiate U' and evaluate sign at stationary point.`}
            solution={`U''(x)=k+3β x^2. At x*=0, U''=k>0 (stable). For other stationary points compute x* and plug in to check positivity.`}
          />

          <Exercise
            title="Exercise 16.4.2 — Energy wells and bifurcation"
            body={`Plot U(x) as β varies and observe when extra minima appear (bifurcation). Describe implications for mechanical stability.`}
            hint={`Extra minima may appear when nonlinear term dominates; examine parameter ranges where U changes shape.`}
            solution={`As β increases positive nonlinearity deepens side wells, possibly creating local minima away from zero leading to multiple equilibria.`}
          />

          <Exercise
            title="Exercise 16.4.3 — 2DOF coupling effects"
            body={`Explore how coupling c influences the location of minima in the 2DOF model. Use the 3D plot and consider slices for y=0 or x=0.`}
            hint={`Positive coupling tilts the energy landscape causing minima to shift along diagonal directions.`}
            solution={`Coupling transfers energy between modes; large c may produce nontrivial coupled equilibria; visualize via contour slices to reason visually.`}
          />
        </div>
      </motion.div>
    </section>
  );
}

// ----------------------------- Problems section (combined) -----------------------------
function SectionProblems() {
  const problems = [
    {
      title: "Problem 16.1 — Tank design (advanced)",
      body: `Consider a cylindrical tank mounted with a reinforced ring of width w that increases cost linearly with circumference: add cost term C_ring = k_r * 2\\pi r. Recompute optimal r,h for V fixed.`,
      solution: `Add k_r term to A(r) and minimize numerically; analytic solution possible but algebra more involved. Use numeric scan to find minimizer.`
    },
    {
      title: "Problem 16.2 — Multi-objective wastewater",
      body: `Formulate a bi-objective optimization: minimize cost and pollutant load simultaneously. Sketch Pareto front for sample coefficients.`,
      solution: `Construct parametric optimization varying weight between objectives and plot Pareto curve; see trade-offs between cost and pollutant load.`
    },
    {
      title: "Problem 16.3 — Circuit networks",
      body: `Extend the two-node circuit to a network: given internal resistances on multiple sources, find load arrangement maximizing total delivered power under constraints.`,
      solution: `Formulate as constrained optimization; potentially use Lagrange multipliers or numeric search for optimal allocations.`
    },
    {
      title: "Problem 16.4 — Structural energy",
      body: `Model small finite-element beam with potential energy functional U(q) quadratic in nodal DOFs plus nonlinear terms. Use gradient-based minimization to find equilibrium and discuss convergence.`,
      solution: `Assemble small stiffness matrix and use Newton or trust-region to find minimizer; observe conditioning effects and choose damping as needed.`
    }
  ];

  return (
    <section className="rounded-lg p-4" style={{ background: THEME.panel }}>
      <motion.div {...fadeUp}>
        <SectionHeader Icon={BookOpen} title="Problems — Case Studies" subtitle="Practice & extensions" />
        <div className="space-y-3">
          {problems.map((p, i) => <Exercise key={i} title={p.title} body={p.body} solution={p.solution} />)}
        </div>
      </motion.div>
    </section>
  );
}

// ----------------------------- Page assembly -----------------------------
export default function Chapter16() {
  return (
    <div className={`min-h-screen p-6 bg-zinc-950 text-zinc-100`}>
      <header className="mb-6">
        <div className="flex items-center justify-between">
          <div>
            <h1 className="text-3xl font-bold text-emerald-400" >Case Studies: Optimization</h1>
            <div className="text-sm text-zinc-400 mt-1">Applied problems across engineering disciplines with visual and numeric demos</div>
          </div>
        </div>
      </header>

      <main className="space-y-8">
        <Section161 />
        <Section162 />
        <Section163 />
        <Section164 />
        <SectionProblems />
        <BottomBar/>
      </main>

      <footer className="mt-8 text-xs text-zinc-500">Numerix Lab — Educational visualizations. For production solve with professional solvers and validated models.</footer>
    </div>
  );
}
