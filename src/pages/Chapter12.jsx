// src/pages/Chapter12.jsx
// ======================================================================
// Chapter 12 — Case Studies: Linear Algebraic Equations
// Focus: Building real systems that reduce to Ax = b and solving them
// UI: dark emerald + neon cyan theme, framer-motion micro-animations
// Plots: 2D diagrams (truss, circuits, masses) using react-plotly
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

import { Card, CardHeader, CardContent, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Button } from "@/components/ui/button";
import { Textarea } from "@/components/ui/textarea";

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
  Atom,
  Beaker,
  CircuitBoard,
  Link2,
  Grid2X2,
  WineIcon,
} from "lucide-react";
import BottomBar from "../components/BottomBar";

// ======================================================================
// Theme helpers (same spirit as Chapter 7 for visual consistency)
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

const fadeUp = {
  initial: { opacity: 0, y: 16 },
  animate: { opacity: 1, y: 0 },
  transition: { duration: 0.3 },
};

// ======================================================================
// Utilities — parsing, linear algebra, helpers
// ======================================================================

/** parseMatrix: string -> number[][] 
 * Accepts rows separated by newline or ';'
 * Columns separated by space or comma.
 * e.g. "3  -1  0\n-1  3 -1\n0 -1 3"
 */
function parseMatrix(s) {
  const rows = s
    .trim()
    .split(/[\n;]+/g)
    .filter(Boolean);
  return rows.map((row) =>
    row
      .trim()
      .split(/[,\s]+/g)
      .filter(Boolean)
      .map((v) => Number(v))
  );
}

/** parseVector: string -> number[] (space/comma/newline separated) */
function parseVector(s) {
  return s
    .trim()
    .split(/[,\s\n]+/g)
    .filter(Boolean)
    .map((v) => Number(v));
}

/** zeros */
function zeros(n) {
  return Array.from({ length: n }, () => 0);
}

/** identity */
function identity(n) {
  const I = Array.from({ length: n }, () => Array.from({ length: n }, () => 0));
  for (let i = 0; i < n; i++) I[i][i] = 1;
  return I;
}

/** matVec */
function matVec(A, x) {
  const n = A.length;
  const y = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    let s = 0;
    const Ai = A[i];
    for (let j = 0; j < Ai.length; j++) s += Ai[j] * x[j];
    y[i] = s;
  }
  return y;
}

/** vecSub */
function vecSub(a, b) {
  return a.map((v, i) => v - b[i]);
}

/** vecAdd */
function vecAdd(a, b) {
  return a.map((v, i) => v + b[i]);
}

/** vecNorm2 */
function vecNorm2(x) {
  let s = 0;
  for (let i = 0; i < x.length; i++) s += x[i] * x[i];
  return Math.sqrt(s);
}

/** deep copy */
function deepCopy(A) {
  return A.map((row) => row.slice());
}

/** format number */
function fmt(v, p = 6) {
  if (!Number.isFinite(v)) return "NaN";
  const s = v.toFixed(p);
  // trim trailing zeros
  return s.replace(/(?:\.0+|(\.\d+?)0+)$/, "$1");
}

/** Gaussian elimination with partial pivoting (in-place) 
 * returns x, plus info: {x, Pivots, swaps, ok}
 */
function solveGEPP(Ain, bin) {
  const A = deepCopy(Ain);
  const b = bin.slice();
  const n = A.length;
  for (let k = 0; k < n - 1; k++) {
    // pivot
    let piv = k;
    let maxa = Math.abs(A[k][k]);
    for (let i = k + 1; i < n; i++) {
      const val = Math.abs(A[i][k]);
      if (val > maxa) {
        maxa = val;
        piv = i;
      }
    }
    if (maxa === 0 || !Number.isFinite(maxa)) {
      return { x: Array(n).fill(NaN), ok: false, msg: "Singular or ill-conditioned pivot" };
    }
    if (piv !== k) {
      const tmp = A[k];
      A[k] = A[piv];
      A[piv] = tmp;
      const tb = b[k];
      b[k] = b[piv];
      b[piv] = tb;
    }
    // elimination
    for (let i = k + 1; i < n; i++) {
      const factor = A[i][k] / A[k][k];
      A[i][k] = 0;
      for (let j = k + 1; j < n; j++) {
        A[i][j] -= factor * A[k][j];
      }
      b[i] -= factor * b[k];
    }
  }
  if (Math.abs(A[n - 1][n - 1]) < 1e-18) {
    return { x: Array(n).fill(NaN), ok: false, msg: "Singular matrix at last pivot" };
  }
  // back-substitution
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    let s = b[i];
    for (let j = i + 1; j < n; j++) s -= A[i][j] * x[j];
    x[i] = s / A[i][i];
  }
  return { x, ok: true };
}

/** residual & relative residual */
function residualStats(A, x, b) {
  const r = vecSub(matVec(A, x), b);
  const nr = vecNorm2(r);
  const nb = Math.max(vecNorm2(b), 1e-16);
  return { r, nr, rel: nr / nb };
}

/** simple condition estimate via ||A||_inf * ||A^{-1}||_inf (rough)
 * We compute ||A^{-1}||_inf by solving A y = e_i for each i (columns of inverse).
 */
function condInfEstimate(Ain) {
  const n = Ain.length;
  // ||A||inf
  let normA = 0;
  for (let i = 0; i < n; i++) {
    let s = 0;
    for (let j = 0; j < n; j++) s += Math.abs(Ain[i][j]);
    normA = Math.max(normA, s);
  }
  if (normA === 0) return Infinity;
  let normInv = 0;
  for (let k = 0; k < n; k++) {
    const e = zeros(n);
    e[k] = 1;
    const sol = solveGEPP(Ain, e);
    if (!sol.ok) return Infinity;
    const y = sol.x;
    let s = 0;
    for (let i = 0; i < n; i++) s += Math.abs(y[i]);
    normInv = Math.max(normInv, s);
  }
  return normA * normInv;
}

/** helper to build LaTeX matrix from A */
function latexMatrix(A) {
  const rows = A
    .map((row) => row.map((v) => fmt(v, 3)).join(" & "))
    .join(" \\\\ ");
  return `\\begin{bmatrix} ${rows} \\end{bmatrix}`;
}

/** helper to build LaTeX vector from x */
function latexVector(x) {
  const rows = x.map((v) => fmt(v, 3)).join(" \\\\ ");
  return `\\begin{bmatrix} ${rows} \\end{bmatrix}`;
}

/** clamp and safe numeric input */
function safeNumber(v, d = 0) {
  const n = Number(v);
  return Number.isFinite(n) ? n : d;
}

// ======================================================================
// Shared small components (style-consistent with Chapter 7)
// ======================================================================

function Equation({ children, color = theme.accent }) {
  return (
    <div className="my-2">
      <div className="px-3 py-2 rounded-lg text-gray-200 border border-zinc-700 bg-zinc-900/80">
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
        <Button variant="ghost" className="bg-white cursor-pointer hover:bg-gray-200" size="sm" onClick={() => setOpen((o) => !o)}>
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
// 12.1 — Steady-State Analysis of a System of Reactors
// Model: N perfectly-mixed CSTRs with first-order decay, feed & interconnections
// At steady-state: A C = b (linear in concentrations)
// ======================================================================
function buildReactorSystem({ N, F, V, k, Cin, coupling }) {
  // N reactors, each: 0 = (F/V)(in - Ci) - k_i Ci + sum_j alpha_ij (Cj - Ci)
  // Rearranged for A C = b
  // Diagonal: (F/V + k_i + sum_j alpha_ij)
  // Offdiag: -alpha_ij
  // RHS: (F/V) * Cin_i
  const A = Array.from({ length: N }, () => Array.from({ length: N }, () => 0));
  const b = Array.from({ length: N }, () => 0);
  for (let i = 0; i < N; i++) {
    let sumAlpha = 0;
    for (let j = 0; j < N; j++) if (i !== j) sumAlpha += coupling[i][j];
    const diag = F / V + k[i] + sumAlpha;
    A[i][i] = diag;
    for (let j = 0; j < N; j++) {
      if (i === j) continue;
      if (coupling[i][j] !== 0) {
        A[i][j] = -coupling[i][j];
      }
    }
    b[i] = (F / V) * Cin[i];
  }
  return { A, b };
}

function Section121() {
  const [N, setN] = useState(3);
  const [F, setF] = useState(2); // volumetric flow rate
  const [V, setV] = useState(5); // volume (same for all)
  const [CinStr, setCinStr] = useState("1 0.5 0"); // inlet concentrations (per reactor feed)
  const [kStr, setKStr] = useState("0.2 0.1 0.05"); // first-order decay
  const [alphaStr, setAlphaStr] = useState(
    // symmetric weak coupling (mass transfer)
    "0 0.2 0.0\n0.2 0 0.1\n0.0 0.1 0"
  );

  const Cin = useMemo(() => {
    let c = parseVector(CinStr);
    if (c.length !== N) c = Array.from({ length: N }, (_, i) => c[i] ?? 0);
    return c;
  }, [CinStr, N]);

  const k = useMemo(() => {
    let v = parseVector(kStr);
    if (v.length !== N) v = Array.from({ length: N }, (_, i) => v[i] ?? 0);
    return v;
  }, [kStr, N]);

  const coupling = useMemo(() => {
    let M = parseMatrix(alphaStr);
    // force NxN
    if (M.length !== N) M = Array.from({ length: N }, (_, i) => Array.from({ length: N }, (__, j) => (M[i]?.[j] ?? 0)));
    for (let i = 0; i < N; i++) {
      if (M[i].length !== N) M[i] = Array.from({ length: N }, (__, j) => (M[i]?.[j] ?? 0));
      for (let j = 0; j < N; j++) M[i][j] = Number(M[i][j]);
    }
    return M;
  }, [alphaStr, N]);

  const sys = useMemo(() => buildReactorSystem({ N, F: Number(F), V: Number(V), k, Cin, coupling }), [N, F, V, k, Cin, coupling]);

  const sol = useMemo(() => solveGEPP(sys.A, sys.b), [sys]);
  const stats = useMemo(() => (sol.ok ? residualStats(sys.A, sol.x, sys.b) : { r: [], nr: NaN, rel: NaN }), [sys, sol]);
  const condEst = useMemo(() => condInfEstimate(sys.A), [sys]);

  // Simple bar chart of concentrations
  const concentrations = useMemo(() => (sol.ok ? sol.x : Array.from({ length: N }, () => NaN)), [sol, N]);

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Beaker className="w-5 h-5" />
            12.1 Steady-State Analysis of a System of Reactors
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="grid md:grid-cols-5 gap-3">
            <Labeled label="Number of reactors N">
              <Input value={N} onChange={(e) => setN(Math.max(2, Math.min(8, safeNumber(e.target.value, 3))))} className="bg-zinc-800 text-white" />
            </Labeled>
            <Labeled label="Flow rate F">
              <Input value={F} onChange={(e) => setF(e.target.value)} className="bg-zinc-800 text-white" />
            </Labeled>
            <Labeled label="Volume V">
              <Input value={V} onChange={(e) => setV(e.target.value)} className="bg-zinc-800 text-white" />
            </Labeled>
            <Labeled label="Inlet concentrations Cin (length N)">
              <Input value={CinStr} onChange={(e) => setCinStr(e.target.value)} className="bg-zinc-800 text-white" />
            </Labeled>
            <Labeled label="Decay rates k (length N)">
              <Input value={kStr} onChange={(e) => setKStr(e.target.value)} className="bg-zinc-800 text-white" />
            </Labeled>
          </div>

          <div className="grid md:grid-cols-2 gap-3">
            <Labeled label="Coupling α (N×N; α_ij is transfer rate from j→i)">
              <Textarea value={alphaStr} onChange={(e) => setAlphaStr(e.target.value)} className="bg-zinc-800 text-white h-32" />
            </Labeled>
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-3">
              <div className="text-sm text-zinc-300 mb-2">System matrix and solution</div>
             <Equation>
  {`Set up $A\\mathbf{C} = \\mathbf{b}$ with \n\n
  $$${latexMatrix(sys.A)} \\; \\mathbf{C} = ${latexVector(sys.b)}.$$`}
</Equation>

              {sol.ok ? (
                <Equation color={theme.accent2}>
              {`Solution $\\mathbf{C} \\approx ${latexVector(sol.x)}$.`}
            </Equation>

              ) : (
                <div className="text-red-400 text-sm">Solver failed: {sol.msg ?? "unknown error"}</div>
              )}
            </div>
          </div>

          <div className="grid md:grid-cols-2 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2">
              <Plot
                data={[
                  {
                    x: Array.from({ length: N }, (_, i) => `R${i + 1}`),
                    y: concentrations,
                    type: "bar",
                    name: "C steady",
                  },
                ]}
                layout={{
                  title: "Reactor steady-state concentrations",
                  paper_bgcolor: "rgba(0,0,0,0)",
                  plot_bgcolor: "rgba(0,0,0,0)",
                  font: { color: "#e5e7eb" },
                  margin: { l: 40, r: 10, b: 40, t: 40 },
                }}
                style={{ width: "100%", height: 300 }}
                useResizeHandler
              />
            </div>
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-3">
              <Summary
                items={[
                  { label: "Residual ‖A C − b‖₂", value: fmt(stats.nr) },
                  { label: "Relative residual", value: fmt(stats.rel) },
                  { label: "cond₍∞₎(A) (est.)", value: Number.isFinite(condEst) ? fmt(condEst, 3) : "∞" },
                ]}
              />
              <Note className="mt-2">
                Larger coupling or decay increases diagonal dominance and usually improves conditioning. Extreme parameter scales can harm conditioning.
              </Note>
            </div>
          </div>

          <div className="grid md:grid-cols-2 gap-3">
            <Equation>{`**Balance per reactor**\n\nAt steady-state, inflow = outflow + decay ± transfers.\n\n$$\\left(\\tfrac{F}{V} + k_i + \\sum_{j\\ne i}\\alpha_{ij}\\right)C_i - \\sum_{j\\ne i} \\alpha_{ij} C_j = \\tfrac{F}{V} C^{\\text{in}}_i.$$`}</Equation>
           
          </div>

          <div className="space-y-3">
           <Exercise
            title="Exercise 12.1.1 — Equal reactors with ring coupling"
            body={`Set $N=4$, $F=2$, $V=4$, $k_i=0.1$, $C^{in}=[1,1,0,0]$.  
          Use a ring: $\\alpha_{12}=\\alpha_{23}=\\alpha_{34}=\\alpha_{41}=0.3$ and reverse edges zero. Solve for $C$.`}
            solution={`Construct $A$ as per the formula; solve to get a flattened gradient from reactors 1→2→3→4. Exact numbers depend on parameters; verify with the app.`}
          />

          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// 12.2 — Statically Determinate Truss
// Solve member axial forces via joint equilibrium: A f = b
// We provide a simple planar truss with user-configurable loads/supports.
// ======================================================================
function buildSimpleTruss() {
  // Nodes (x,y)
  const nodes = [
    [0, 0],     // 0
    [2, 0],     // 1
    [1, Math.sqrt(3)], // 2 (top equilateral-ish)
  ];
  // Members (node_i, node_j)
  const members = [
    [0, 1], // bottom left-right
    [0, 2], // left to top
    [1, 2], // right to top
  ];
  // Supports: node 0 pinned (ux=uy=0), node 1 roller (uy=0)
  const supports = { 0: { ux: true, uy: true }, 1: { ux: false, uy: true } };
  // External loads: [Fx,Fy] at nodes
  const loads = [
    [0, 0],    // node 0
    [0, 0],    // node 1
    [0, -10],  // node 2 (downward load)
  ];
  return { nodes, members, supports, loads };
}

/** Build equilibrium equations for axial forces in members (tension positive)
 * Unknowns: member forces + support reactions (as needed for DOFs restrained)
 */
function buildTrussSystem({ nodes, members, supports, loads }) {
  const nNodes = nodes.length;
  const nMembers = members.length;

  // determine reactions: for each support, add unknowns for restrained DOFs
  const reactionMap = []; // array of {node, dir:'x'|'y'}
  for (let i = 0; i < nNodes; i++) {
    const s = supports[i] || {};
    if (s.ux) reactionMap.push({ node: i, dir: "x" });
    if (s.uy) reactionMap.push({ node: i, dir: "y" });
  }
  const nReactions = reactionMap.length;

  // total unknowns
  const nUnknowns = nMembers + nReactions;

  // equations: 2 per node (ΣFx=0, ΣFy=0) → 2nNodes
  const mEq = 2 * nNodes;

  const A = Array.from({ length: mEq }, () => Array.from({ length: nUnknowns }, () => 0));
  const b = Array.from({ length: mEq }, () => 0);

  // helper: direction cosines for each member
  const dir = members.map(([i, j]) => {
    const [xi, yi] = nodes[i];
    const [xj, yj] = nodes[j];
    const dx = xj - xi;
    const dy = yj - yi;
    const L = Math.hypot(dx, dy);
    return { cx: dx / L, cy: dy / L, i, j };
  });

  // assemble ΣFx, ΣFy at each node
  for (let n = 0; n < nodes.length; n++) {
    const rowFx = 2 * n;
    const rowFy = 2 * n + 1;

    // member contributions
    dir.forEach((d, idx) => {
      const sign = n === d.i ? +1 : n === d.j ? -1 : 0;
      if (sign !== 0) {
        A[rowFx][idx] += sign * d.cx;
        A[rowFy][idx] += sign * d.cy;
      }
    });

    // reactions for restrained DOFs at this node
    reactionMap.forEach((r, k) => {
      if (r.node === n) {
        const col = nMembers + k;
        if (r.dir === "x") A[rowFx][col] = 1;
        if (r.dir === "y") A[rowFy][col] = 1;
      }
    });

    // RHS is negative of applied loads (move to other side)
    b[rowFx] = -loads[n][0];
    b[rowFy] = -loads[n][1];
  }

  return { A, b, dir, nMembers, reactionMap };
}

function Section122() {
  const base = useMemo(() => buildSimpleTruss(), []);
  const [loadTop, setLoadTop] = useState(-10);
  const [scale, setScale] = useState(1.0);

  const truss = useMemo(() => {
    const nodes = base.nodes.map(([x, y]) => [x * scale, y * scale]);
    const loads = [[0, 0], [0, 0], [0, Number(loadTop)]];
    return { ...base, nodes, loads };
  }, [base, loadTop, scale]);

  const sys = useMemo(() => buildTrussSystem(truss), [truss]);
  const sol = useMemo(() => solveGEPP(sys.A, sys.b), [sys]);
  const stats = useMemo(() => (sol.ok ? residualStats(sys.A, sol.x, sys.b) : { r: [], nr: NaN, rel: NaN }), [sys, sol]);

  // member forces and reactions
  const memberForces = useMemo(() => {
    if (!sol.ok) return [];
    return sol.x.slice(0, sys.nMembers);
  }, [sol, sys]);

  const reactions = useMemo(() => {
    if (!sol.ok) return [];
    return sol.x.slice(sys.nMembers);
  }, [sol, sys]);

  // Visualization: color members by tension (+) vs compression (−)
  const memberColors = memberForces.map((f) => (f >= 0 ? theme.accent : theme.warn));
  const memberWidths = memberForces.map((f) => 2 + 2 * Math.min(4, Math.abs(f) / 10));

  // Build traces
  const nodeXs = truss.nodes.map((p) => p[0]);
  const nodeYs = truss.nodes.map((p) => p[1]);

  const memberTraces = sys.dir.map((d, idx) => {
    const xi = truss.nodes[d.i][0];
    const yi = truss.nodes[d.i][1];
    const xj = truss.nodes[d.j][0];
    const yj = truss.nodes[d.j][1];
    return {
      x: [xi, xj],
      y: [yi, yj],
      type: "scatter",
      mode: "lines",
      line: { width: memberWidths[idx], color: memberColors[idx] },
      name: `m${idx + 1} (${memberForces[idx] ? fmt(memberForces[idx], 3) : "?"})`,
      hoverinfo: "name",
    };
  });

  const nodeTrace = {
    x: nodeXs,
    y: nodeYs,
    type: "scatter",
    mode: "markers+text",
    marker: { size: 10 },
    text: nodeXs.map((_, i) => `N${i}`),
    textposition: "top center",
    name: "nodes",
  };

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Link2 className="w-5 h-5" />
            12.2 Analysis of a Statically Determinate Truss
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="grid md:grid-cols-4 gap-3">
            <Labeled label="Top node vertical load (− down)">
              <Input value={loadTop} onChange={(e) => setLoadTop(e.target.value)} className="bg-zinc-800 text-white" />
            </Labeled>
            <Labeled label="Geometry scale">
              <Input value={scale} onChange={(e) => setScale(e.target.value)} className="bg-zinc-800 text-white" />
            </Labeled>
            <Summary
              items={[
                { label: "Residual ‖A x − b‖₂", value: fmt(stats.nr) },
                { label: "Relative residual", value: fmt(stats.rel) },
                { label: "Unknowns", value: sys.A[0].length },
              ]}
            />
            <div className="text-sm text-zinc-300">
              <div className="mb-1">Member forces (tension +, compression −):</div>
              <div className="text-zinc-100 text-xs">
                {memberForces.map((f, i) => (
                  <div key={i}>
                    m{i + 1}: <span style={{ color: f >= 0 ? theme.accent : theme.warn }}>{fmt(f, 3)}</span>
                  </div>
                ))}
              </div>
            </div>
          </div>

          <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2">
            <Plot
              data={[...memberTraces, nodeTrace]}
              layout={{
                title: "Planar truss (member color: + tension / − compression)",
                xaxis: { title: "x", color: "#e5e7eb", scaleanchor: "y", scaleratio: 1, zeroline: false, showgrid: false },
                yaxis: { title: "y", color: "#e5e7eb", zeroline: false, showgrid: false },
                paper_bgcolor: "rgba(0,0,0,0)",
                plot_bgcolor: "rgba(0,0,0,0)",
                font: { color: "#e5e7eb" },
                margin: { l: 40, r: 10, b: 40, t: 40 },
              }}
              style={{ width: "100%", height: 330 }}
              useResizeHandler
            />
          </div>

          <div className="grid md:grid-cols-2 gap-3">
           <Equation>
            {`**Joint equilibrium**  
                For each node, $$\\sum F_x = 0, \\quad \\sum F_y = 0.$$  
                Member axial force contributes components with direction cosines $(c_x, c_y)$.  
                Reactions enter as unknowns at restrained DOFs.`}</Equation>

            <Note>For determinate planar trusses, the number of unknowns equals the number of equilibrium equations. Tension/compression sign follows the assumed positive axial direction.</Note>
          </div>

          <div className="space-y-3">
            <Exercise
              title="Exercise 12.2.1 — Symmetry check"
              body={`With symmetric geometry and load at the top node, verify equal magnitudes in the two inclined members and zero horizontal reaction at the roller support.`}
              solution={`By symmetry, inclined members carry equal force; the roller at node 1 carries only vertical reaction. The pinned support at node 0 carries vertical and horizontal reactions.`}
            />
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// 12.3 — Currents and Voltages in Resistor Circuits (Nodal Analysis)
// Build conductance matrix G and solve G V = I
// ======================================================================
function buildNodalSystem({ nNodes, resistors, currents, ground = 0 }) {
  // Nodes: 0..nNodes-1; ground held at 0V by removing its equation/unknown.
  // resistors: [{n1, n2, R}] (n2 can equal -1 to denote resistor to ground)
  // currents: injection at nodes (A), positive into node
  const G = Array.from({ length: nNodes }, () => Array.from({ length: nNodes }, () => 0));
  const I = Array.from({ length: nNodes }, (_, i) => currents[i] ?? 0);

  resistors.forEach((r) => {
    const { n1, n2, R } = r;
    const g = 1 / R;
    if (n2 >= 0) {
      G[n1][n1] += g;
      G[n2][n2] += g;
      G[n1][n2] -= g;
      G[n2][n1] -= g;
    } else {
      // to ground
      G[n1][n1] += g;
    }
  });

  // Apply ground by reducing system
  const keep = [];
  for (let i = 0; i < nNodes; i++) if (i !== ground) keep.push(i);

  const A = keep.map((ri) => keep.map((ci) => G[ri][ci]));
  const b = keep.map((ri) => I[ri]);

  return { A, b, map: keep, ground };
}

function Section123() {
  const [nNodes, setNNodes] = useState(4);
  const [ground, setGround] = useState(0);
  const [resStr, setResStr] = useState(
    // n1 n2 R ; use -1 for ground
    "0 1 100\n1 2 200\n2 3 100\n3 -1 50\n0 -1 100"
  );
  const [currStr, setCurrStr] = useState("0 0 0 0"); // current injections (A)

  const resistors = useMemo(() => {
    const rows = resStr
      .trim()
      .split(/[\n;]+/g)
      .filter(Boolean)
      .map((line) => line.trim().split(/[,\s]+/g).map(Number));
    return rows.map(([n1, n2, R]) => ({ n1, n2, R }));
  }, [resStr]);

  const currents = useMemo(() => {
    let v = parseVector(currStr);
    if (v.length !== nNodes)
      v = Array.from({ length: nNodes }, (_, i) => v[i] ?? 0);
    return v;
  }, [currStr, nNodes]);

  const sys = useMemo(
    () =>
      buildNodalSystem({
        nNodes: Number(nNodes),
        resistors,
        currents,
        ground: Number(ground),
      }),
    [nNodes, resistors, currents, ground]
  );

  const sol = useMemo(() => solveGEPP(sys.A, sys.b), [sys]);

  const stats = useMemo(
    () =>
      sol.ok
        ? residualStats(sys.A, sol.x, sys.b)
        : { r: [], nr: NaN, rel: NaN },
    [sys, sol]
  );

  // reconstruct full node voltages (insert 0V at ground)
  const voltages = useMemo(() => {
    const V = Array.from({ length: Number(nNodes) }, () => 0);
    if (sol.ok) {
      sys.map.forEach((nodeIdx, i) => (V[nodeIdx] = sol.x[i]));
      V[sys.ground] = 0;
    } else {
      sys.map.forEach((nodeIdx) => (V[nodeIdx] = NaN));
    }
    return V;
  }, [sol, sys, nNodes]);

  // simple ladder plot: nodes along a line, voltages as y
  const ladderX = Array.from({ length: Number(nNodes) }, (_, i) => i);
  const ladderY = voltages;

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <CircuitBoard className="w-5 h-5" />
            12.3 Currents and Voltages in Resistor Circuits — Nodal Analysis
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          {/* Inputs and summary */}
          <div className="grid md:grid-cols-5 gap-3">
            <Labeled label="Nodes">
              <Input
                value={nNodes}
                onChange={(e) =>
                  setNNodes(
                    Math.max(2, Math.min(10, safeNumber(e.target.value, 4)))
                  )
                }
                className="bg-zinc-800 text-white"
              />
            </Labeled>
            <Labeled label="Ground node index">
              <Input
                value={ground}
                onChange={(e) =>
                  setGround(
                    Math.max(
                      0,
                      Math.min(Number(nNodes) - 1, safeNumber(e.target.value, 0))
                    )
                  )
                }
                className="bg-zinc-800 text-white"
              />
            </Labeled>
            <Labeled label="Resistors (n1 n2 R per line; use -1 for ground)">
              <Textarea
                value={resStr}
                onChange={(e) => setResStr(e.target.value)}
                className="bg-zinc-800 text-white h-32"
              />
            </Labeled>
            <Labeled label="Current injections I (length n; + into node, A)">
              <Input
                value={currStr}
                onChange={(e) => setCurrStr(e.target.value)}
                className="bg-zinc-800 text-white"
              />
            </Labeled>
            <Summary
              items={[
                { label: "Residual ‖G V − I‖₂", value: fmt(stats.nr) },
                { label: "Relative residual", value: fmt(stats.rel) },
                { label: "Unknown voltages", value: sys.A[0].length },
              ]}
            />
          </div>

          {/* Plots and voltages */}
          <div className="grid md:grid-cols-2 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2">
              <Plot
                data={[
                  {
                    x: ladderX,
                    y: ladderY,
                    type: "scatter",
                    mode: "markers+lines",
                    name: "Node voltages",
                  },
                ]}
                layout={{
                  title: "Node voltages (ground at 0 V)",
                  xaxis: { title: "node index", color: "#e5e7eb" },
                  yaxis: { title: "V (volts)", color: "#e5e7eb" },
                  paper_bgcolor: "rgba(0,0,0,0)",
                  plot_bgcolor: "rgba(0,0,0,0)",
                  font: { color: "#e5e7eb" },
                  margin: { l: 40, r: 10, b: 40, t: 40 },
                }}
                style={{ width: "100%", height: 300 }}
                useResizeHandler
              />
            </div>
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-3">
              <div className="text-sm text-zinc-300 mb-2">Computed voltages</div>
              <Equation color={theme.accent2}>
  {`$\\mathbf{V} \\approx ${latexVector(voltages)}$`}
</Equation>

             
            </div>
          </div>

          {/* Equations and notes */}
          <div className="grid md:grid-cols-2 gap-3">
           <Equation>
  {`**Nodal analysis**\n\nConstruct the conductance matrix $G$ by superposition: for each resistor, add conductance on diagonals and subtract on off-diagonals. Then solve $$G\\mathbf{V} = \\mathbf{I}.$$`}
</Equation>

            <Note>
              Choosing a ground reference removes the singularity (rank
              deficiency). Poor conditioning can occur with extreme resistor
              ratios.
            </Note>
          </div>

          {/* Exercises */}
          <div className="space-y-3">
            <Exercise
              title="Exercise 12.3.1 — Divider"
              body="Build a series divider: 3 nodes (0,1,2), ground at 0, resistors: (1–0) 1kΩ, (2–1) 1kΩ, (2–ground) 1kΩ. Inject +1 mA at node 2. Compute node voltages."
              solution="Symmetry suggests node 1 voltage is roughly mid of node 2 and ground. Exact values follow from G V = I; verify numerically."
            />
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}


// ======================================================================
// 12.4 — Spring-Mass Systems (static)
// Linear system: K u = f
// We model a 1D chain of masses connected by springs to neighbors and walls.
// ======================================================================
function buildSpringChain({ nMasses, kLeft, kRight, kMid, fLoads, leftFixed = true, rightFixed = true }) {
  // DOFs = displacements of masses: u1..un
  // Stiffness K: tridiagonal with kMid plus wall springs
  const n = nMasses;
  const K = Array.from({ length: n }, () => Array.from({ length: n }, () => 0));
  for (let i = 0; i < n; i++) {
    // contributions from middle springs between i and i+1
    if (i > 0) {
      K[i][i] += kMid;
      K[i][i - 1] -= kMid;
      K[i - 1][i] -= kMid;
      K[i - 1][i - 1] += kMid;
    }
  }
  // wall springs
  if (leftFixed) K[0][0] += kLeft;
  if (rightFixed) K[n - 1][n - 1] += kRight;

  const f = Array.from({ length: n }, (_, i) => fLoads[i] ?? 0);
  return { K, f };
}

function Section124() {
  const [n, setN] = useState(5);
  const [kMid, setKMid] = useState(100);
  const [kLeft, setKLeft] = useState(200);
  const [kRight, setKRight] = useState(200);
  const [fStr, setFStr] = useState("0 0 -10 0 0");
  const [fixL, setFixL] = useState(true);
  const [fixR, setFixR] = useState(true);
  const [scale, setScale] = useState(1.0);

  const loads = useMemo(() => {
    let v = parseVector(fStr);
    if (v.length !== n) v = Array.from({ length: n }, (_, i) => v[i] ?? 0);
    return v;
  }, [fStr, n]);

  const sys = useMemo(
    () =>
      buildSpringChain({
        nMasses: Number(n),
        kLeft: Number(kLeft),
        kRight: Number(kRight),
        kMid: Number(kMid),
        fLoads: loads,
        leftFixed: fixL,
        rightFixed: fixR,
      }),
    [n, kLeft, kRight, kMid, loads, fixL, fixR]
  );

  const sol = useMemo(() => solveGEPP(sys.K, sys.f), [sys]);
  const stats = useMemo(() => (sol.ok ? residualStats(sys.K, sol.x, sys.f) : { r: [], nr: NaN, rel: NaN }), [sys, sol]);

  const u = sol.ok ? sol.x : Array.from({ length: n }, () => NaN);

  // Visualization: masses along x, displaced in y
  const xPos = Array.from({ length: n }, (_, i) => i);
  const yPos = u.map((v) => (Number.isFinite(v) ? v * scale : NaN));

  // springs lines
  const springTraces = [];
  for (let i = 0; i < n - 1; i++) {
    springTraces.push({
      x: [xPos[i], xPos[i + 1]],
      y: [yPos[i], yPos[i + 1]],
      type: "scatter",
      mode: "lines",
      line: { width: 3 },
      name: `spring ${i + 1}`,
      hoverinfo: "none",
    });
  }

  const massTrace = {
    x: xPos,
    y: yPos,
    type: "scatter",
    mode: "markers+text",
    marker: { size: 12 },
    text: xPos.map((i) => `m${i + 1}`),
    textposition: "top center",
    name: "masses",
  };

  return (
    <motion.div {...fadeUp}>
      <Card className={`p-4 rounded-2xl shadow-xl ${theme.panel}`}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <WineIcon className="w-5 h-5" />
            12.4 Spring–Mass Systems (Static) — Solve K u = f
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="grid md:grid-cols-6 gap-3">
            <Labeled label="Masses (n)">
              <Input value={n} onChange={(e) => setN(Math.max(2, Math.min(12, safeNumber(e.target.value, 5))))} className="bg-zinc-800 text-white" />
            </Labeled>
            <Labeled label="k_mid (neighbor springs)">
              <Input value={kMid} onChange={(e) => setKMid(e.target.value)} className="bg-zinc-800 text-white" />
            </Labeled>
            <Labeled label="k_left (wall)">
              <Input value={kLeft} onChange={(e) => setKLeft(e.target.value)} className="bg-zinc-800 text-white" />
            </Labeled>
            <Labeled label="k_right (wall)">
              <Input value={kRight} onChange={(e) => setKRight(e.target.value)} className="bg-zinc-800 text-white" />
            </Labeled>
            <Labeled label="Loads f (length n)">
              <Input value={fStr} onChange={(e) => setFStr(e.target.value)} className="bg-zinc-800 text-white" />
            </Labeled>
            <Labeled label="Viz scale (y)">
              <Input value={scale} onChange={(e) => setScale(e.target.value)} className="bg-zinc-800 text-white" />
            </Labeled>
          </div>

          <div className="flex gap-3 items-center">
            <label className="text-sm text-zinc-300 flex items-center gap-2">
              <input type="checkbox" checked={fixL} onChange={(e) => setFixL(e.target.checked)} />
              Left fixed
            </label>
            <label className="text-sm text-zinc-300 flex items-center gap-2">
              <input type="checkbox" checked={fixR} onChange={(e) => setFixR(e.target.checked)} />
              Right fixed
            </label>
          </div>

          <div className="grid md:grid-cols-2 gap-3">
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-2">
              <Plot
                data={[...springTraces, massTrace]}
                layout={{
                  title: "Static displacements (scaled)",
                  xaxis: { title: "index", color: "#e5e7eb" },
                  yaxis: { title: "u (scaled)", color: "#e5e7eb" },
                  paper_bgcolor: "rgba(0,0,0,0)",
                  plot_bgcolor: "rgba(0,0,0,0)",
                  font: { color: "#e5e7eb" },
                  margin: { l: 40, r: 10, b: 40, t: 40 },
                }}
                style={{ width: "100%", height: 300 }}
                useResizeHandler
              />
            </div>
            <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-3">
              <Summary
                items={[
                  { label: "Residual ‖K u − f‖₂", value: fmt(stats.nr) },
                  { label: "Relative residual", value: fmt(stats.rel) },
                  { label: "Size (n×n)", value: `${n}×${n}` },
                ]}
              />
             <Equation>
  {`$K\\mathbf{u} = \\mathbf{f}$, $K$ is tridiagonal with neighbor springs and wall stiffness. Displacements: $\\mathbf{u} \\approx ${latexVector(u)}$.`}
</Equation>

            </div>
          </div>

          <div className="grid md:grid-cols-2 gap-3">
           <Equation>
  {`**Stiffness assembly**\n\nEach spring contributes $k\\begin{bmatrix}1 & -1 \\\\ -1 & 1\\end{bmatrix}$ to the connected DOFs. Wall springs add $k$ to the corresponding diagonal.`}
</Equation>

            <Note>Stiffer chains produce smaller displacements for the same load. Conditioning grows with system size and stiffness contrast.</Note>
          </div>

          <div className="space-y-3">
           <Exercise
  title="Exercise 12.4.1 — Center push"
  body={`For $n=5$ with $k_{mid}=100$, wall springs $k_L=k_R=200$, apply $f=[0,0,-10,0,0]$. Verify symmetric displacement about center.`}
  solution={`By symmetry with equal wall stiffness, the displacement field is symmetric, peaking near the loaded mass.`}
/>

          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// Documentation Panel — left nav + content for Chapter 12
// ======================================================================
const docsText12 = {
  "12.1 CSTR Steady-State": [
    `**Mass balance:** At steady state, set \\(\\tfrac{dC}{dt}=0\\). This leads to algebraic equations relating concentrations in all reactors.`,
    `For reactor \\(i\\):  
    $$\\left(\\tfrac{F}{V}+k_i+\\sum_{j\\ne i}\\alpha_{ij}\\right) C_i \;-\; \\sum_{j\\ne i} \\alpha_{ij} C_j \;=\; \\tfrac{F}{V} C_i^{in}.$$`,
    `The diagonal term combines **through-flow (\\(F/V\\))**, **decay (\\(k_i\\))**, and **exchange coefficients (\\(\\alpha_{ij}\\))**.`,
    `Because diagonal entries dominate, the coefficient matrix is often **diagonally dominant**, giving numerical stability.`,
    `When interconnections are sparse, the resulting matrix is also **sparse**, making it efficient for iterative solvers.`,
  ],

  "12.2 Truss Equilibrium": [
    `At each joint, enforce static equilibrium:  
    $$\\sum F_x = 0, \\quad \\sum F_y = 0.$$`,
    `Unknowns are the **axial member forces** and **support reactions**.`,
    `Each member force is projected into x/y components using **direction cosines**.`,
    `The global equilibrium system is:  
    $$K \\mathbf{u} = \\mathbf{f},$$  
    where \\(K\\) is assembled from element stiffness matrices.`,
    `Sign convention: **Tension = positive**, **Compression = negative**.`,
    `Supports (fixed, roller, hinge) add boundary conditions that reduce the DOFs.`,
    `Redundant trusses may lead to **indeterminate systems**, requiring advanced methods (compatibility or matrix stiffness approach).`,
  ],

  "12.3 Nodal Circuits": [
    `Use **Kirchhoff’s Current Law (KCL):** sum of currents leaving a node = 0.`,
    `Each resistor contributes to the **conductance matrix (G)** as:  
    $$G_{ii} += 1/R, \\quad G_{ij} -= 1/R.$$`,
    `After assembly:  
    $$G \\mathbf{V} = \\mathbf{I},$$  
    where \\(\\mathbf{V}\\) = node voltages, \\(\\mathbf{I}\\) = current injections.`,
    `Choose a **ground reference** node to prevent singularity.`,
    `The conductance matrix \\(G\\) is **symmetric positive definite** for purely resistive networks.`,
    `Voltage sources or supernodes require constraint handling (modified nodal analysis).`,
    `For large circuits, iterative solvers (CG, BiCGSTAB) are efficient due to sparsity.`,
  ],

  "12.4 Spring–Mass (Static)": [
    `Each spring between DOFs \\(i\\) and \\(j\\) contributes:  
    $$k \\begin{bmatrix} 1 & -1 \\\\ -1 & 1 \\end{bmatrix}$$  
    to the global stiffness matrix.`,
    `Wall springs (connected to ground) contribute \\(k\\) to the corresponding diagonal entry.`,
    `Final system of equations:  
    $$K \\mathbf{u} = \\mathbf{f},$$  
    where \\(K\\) is symmetric and positive definite (SPD).`,
    `With fixed supports, displacements at constrained DOFs are set to zero.`,
    `**Numerical issues:** Large stiffness contrasts or very long chains can degrade conditioning.`,
    `Scaling stiffness values or applying preconditioning helps improve solver performance.`,
    `The solution vector \\(\\mathbf{u}\\) gives displacements; spring forces follow from Hooke’s law.`,
  ],

  "General Tips": [
    `Use **partial pivoting** in Gaussian elimination for numerical robustness.`,
    `Iterative methods (Jacobi, Gauss–Seidel, Conjugate Gradient) are effective for large sparse systems.`,
    `Monitor **residual norms** to check convergence of iterative solvers.`,
    `Estimate **condition numbers** to detect ill-conditioning.`,
    `Scaling variables and parameters often improves stability and reduces round-off errors.`,
    `Exploit **symmetry and sparsity** in matrices for computational efficiency.`,
  ],
};

function ChapterDocs12() {
  const [open, setOpen] = useState("12.1 CSTR Steady-State");
  return (
    <div className="w-full mt-6 bg-zinc-900 border border-zinc-700 rounded-xl shadow-lg overflow-hidden">
      <div className="p-3 border-b border-zinc-700 flex items-center justify-between">
        <div className="flex items-center gap-3">
          <BookOpen className="w-5 h-5 text-cyan-400" />
          <div>
            <div className="text-zinc-100 font-semibold">Chapter 12 — Documentation</div>
            <div className="text-zinc-400 text-xs">Theory, notes, and references</div>
          </div>
        </div>
        <div className="text-zinc-400 text-xs">Linear systems in engineering practice</div>
      </div>

      <div className="grid grid-cols-1 md:grid-cols-4">
        <div className="col-span-1 border-r border-zinc-800">
          {Object.keys(docsText12).map((k) => (
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
            {docsText12[open].map((para, i) => (
              <div key={i}>
                <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>
                  {para}
                </ReactMarkdown>
              </div>
            ))}
            <div className="mt-3 text-zinc-400 text-xs">
              References: Ogata (Control), Hibbeler (Mechanics of Materials), Dorf & Svoboda (Circuits), Ogden (Applied Math), Strang (Linear Algebra).
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}

// ======================================================================
// Problems List
// ======================================================================
function Problems12() {
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
          <Exercise
            title="P12.1 — Reactor Network Calibration"
            body={`Given measurements of outlet concentrations \\(C\\) for a 4-reactor network, adjust \\(k_i\\) to match data by solving a sequence of linear systems. Start with \\(k_i=0.1\\) and use least-squares updates.`}
            solution={`Linearize the map from \\(k\\) to \\(C\\) around current estimates, solve normal equations for updates, iterate until residual drops.`}
          />
          <Exercise
            title="P12.2 — Truss with Distributed Load Approximation"
            body={`Approximate a uniformly distributed load along the top chord of a longer truss with equivalent nodal loads. Build equilibrium equations and compute member forces.`}
            solution={`Distribute total load into joints proportionally to tributary lengths. Assemble \\(A\\) accordingly and solve for axial forces and reactions.`}
          />
          <Exercise
            title="P12.3 — Resistive Mesh"
            body={`Create a 3×3 node grid resistive mesh (each edge 100Ω) with ground at the corner and inject +1mA at the opposite corner. Compute node voltages.`}
            solution={`Build \\(G\\) with interior nodes having 4 neighbors (degree 4). Solve reduced system with ground removed.`}
          />
          <Exercise
            title="P12.4 — Spring Chain with Mixed Supports"
            body={`For a chain of 6 masses with only left wall fixed and right free (no wall spring), apply loads \\([0, -5, 0, -5, 0, 0]\\). Compute displacements.`}
            solution={`Remove right wall stiffness; ensure \\(K\\) remains positive definite by geometry and boundary. Solve for \\(u\\).`}
          />
        </CardContent>
      </Card>
    </motion.div>
  );
}

// ======================================================================
// Page assembly
// ======================================================================
export default function Chapter12() {
  return (
    <div className={`p-6 space-y-6 ${theme.bg}`}>
      <motion.div initial={{ opacity: 0, y: -10 }} animate={{ opacity: 1, y: 0 }}>
        <h1 className="text-3xl font-bold" style={{ color: theme.accent2 }}>
          Case Studies: Linear Algebraic Equations
        </h1>
        <div className="text-zinc-400">Real engineering systems reduced to Ax = b, solved and visualized.</div>
      </motion.div>

      {/* 12.1 */}
      <Section121 />

      {/* 12.2 */}
      <Section122 />

      {/* 12.3 */}
      <Section123 />

      {/* 12.4 */}
      <Section124 />

      {/* Docs & Problems */}
      <ChapterDocs12 />
      <Problems12 />
      <BottomBar/>
    </div>
  );
}
