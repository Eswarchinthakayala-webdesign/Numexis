// src/pages/Chapter15.jsx
// ======================================================================
// Chapter 15 — Constrained Optimization
// Sections: 15.1 Linear Programming, 15.2 Nonlinear Constrained Optimization,
// 15.3 Optimization with Software Packages, Problems
// Dark-gray theme, emerald/amber accents, stacked layout
// ======================================================================

import React, { useMemo, useState } from "react";
import Plot from "react-plotly.js";
import { motion } from "framer-motion";
import ReactMarkdown from "react-markdown";
import remarkMath from "remark-math";
import rehypeKatex from "rehype-katex";
import "katex/dist/katex.min.css";

import { Columns, Triangle, Zap, BookOpen, Code, ChevronDown, ChevronRight } from "lucide-react";
import BottomBar from "../components/BottomBar";

// Attempt to reuse your project's UI inputs/buttons; fallback to simple ones
let Input, Button;
try {
  // eslint-disable-next-line
  Input = require("@/components/ui/input").Input;
  // eslint-disable-next-line
  Button = require("@/components/ui/button").Button;
} catch (e) {
  Input = (props) => <input {...props} />;
  Button = (props) => <button {...props} />;
}

// ---------------- Theme & animation ----------------
const THEME = {
  bgClass: "bg-zinc-900",
panelBg: "#09101b",
  text: "#e6eef6",
  muted: "#9ca3af",
  accent: "#15803d", // emerald
  accent2: "#b45309", // amber
};
const fadeUp = { initial: { opacity: 0, y: 8 }, animate: { opacity: 1, y: 0 }, transition: { duration: 0.28 } };

// ---------------- Utility helpers ----------------
// Simple safe function builder for single-variable or 2-variable JS expressions
function makeFun(expr, vars = ["x"]) {
  try {
    const varList = vars.join(",");
    const wrapper = `
      const { sin, cos, tan, asin, acos, atan, exp, log, sqrt, pow, abs, PI, E } = Math;
      return (${varList}) => { return ${expr}; };
    `;
    // eslint-disable-next-line no-new-func
    return new Function(wrapper)();
  } catch (e) {
    return null;
  }
}

// small vector helpers
function addVec(a, b) { return a.map((v, i) => v + b[i]); }
function subVec(a, b) { return a.map((v, i) => v - b[i]); }
function scaleVec(a, s) { return a.map((v) => v * s); }
function dot(a, b) { return a.reduce((s, v, i) => s + v * b[i], 0); }
function norm2(v) { return Math.sqrt(v.reduce((s, x) => s + x * x, 0)); }

// ---------------- Small components ----------------
function SectionHeader({ Icon, title, subtitle }) {
  return (
    <div className="flex items-center gap-3 mb-3">
      <div style={{ width: 44, height: 44, background: "#071018", display: "flex", alignItems: "center", justifyContent: "center", borderRadius: 8 }}>
        <Icon className="w-5 h-5" />
      </div>
      <div>
        <div style={{ color: THEME.text, fontSize: 18, fontWeight: 700 }}>{title}</div>
        <div style={{ color: THEME.muted, fontSize: 13 }}>{subtitle}</div>
      </div>
    </div>
  );
}

function Exercise({ title, body, solution }) {
  const [open, setOpen] = useState(false);
  return (
    <div className="rounded-xl border border-zinc-700 bg-zinc-900/60 p-3">
      <div className="flex items-center justify-between gap-2">
        <div className="text-zinc-100 font-medium">{title}</div>
        <Button variant="ghost" className="bg-white cursor-pointer py-2 px-2 rounded text-black flex items-center hover:bg-gray-200" size="sm" onClick={() => setOpen((o) => !o)}>
          {open ? <ChevronDown className="w-4 h-4" /> : <ChevronRight className="w-4 h-4" /> } 
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

// ---------------- Documentation (KaTeX-ready) ----------------
const docs = {
  "15.1 Linear Programming": [
    `Linear programming (LP) solves problems of the form\n\n$$\\min\\; c^T x \\quad\\text{s.t. }\\; Ax \\le b,\\; x\\ge0.$$`,
    `The simplex method moves along vertices of the feasible polytope; interior-point methods traverse the interior using barrier or primal-dual systems.`,
  ],
  "15.2 Nonlinear Constrained Optimization": [
    `Nonlinear constrained optimization requires handling equality and/or inequality constraints. KKT (Karush–Kuhn–Tucker) conditions generalize Lagrange multipliers to inequality-constrained problems.`,
    `Penalty and augmented Lagrangian methods convert constrained problems into a sequence of unconstrained subproblems.`,
  ],
  "15.3 Optimization with Software Packages": [
    `Use mature packages: for LP use GLPK, CPLEX, Gurobi, or SciPy's linprog; for nonlinear constrained use IPOPT, NLopt, or SciPy's ` + "`minimize(..., method='SLSQP')`.",
  ],
};

// ---------------- 15.1 Linear Programming: Simplex (educational) ----------------
//
// This is an educational simplex implementation for small problems in standard form:
// maximize c^T x subject to Ax <= b, x>=0 (we will convert to standard maximize with slack variables).
// NOTE: This is not intended as a production solver — use libraries for large/robust solves.
//
function simplexSolve(c, A, b) {
  // Convert to maximize with slack variables: A x + I s = b, x,s >= 0
  const m = A.length;
  const n = c.length;
  const tableau = []; // rows: m constraints + 1 objective row
  // build tableau: [A | I | b]
  for (let i = 0; i < m; i++) {
    const row = A[i].slice();
    for (let j = 0; j < m; j++) row.push(i === j ? 1 : 0); // slack columns
    row.push(b[i]);
    tableau.push(row);
  }
  // objective row (note: simplex tableau uses -c for maximization)
  const objRow = c.map((v) => -v).concat(Array(m).fill(0)).concat([0]);
  tableau.push(objRow);
  const numCols = tableau[0].length;
  const basis = Array.from({ length: m }, (_, i) => n + i); // slack variables initially

  function pivot(col, row) {
    const piv = tableau[row][col];
    // normalize pivot row
    for (let j = 0; j < numCols; j++) tableau[row][j] /= piv;
    // eliminate other rows
    for (let i = 0; i < tableau.length; i++) {
      if (i === row) continue;
      const factor = tableau[i][col];
      for (let j = 0; j < numCols; j++) tableau[i][j] -= factor * tableau[row][j];
    }
    basis[row] = col;
  }

  const trace = [];
  for (let iter = 0; iter < 200; iter++) {
    // find entering column (most negative in objective row)
    const obj = tableau[tableau.length - 1];
    let enter = -1;
    let minVal = 0;
    for (let j = 0; j < numCols - 1; j++) {
      if (obj[j] < minVal) { minVal = obj[j]; enter = j; }
    }
    if (enter === -1) break; // optimal

    // ratio test: min positive b_i / a_i_enter
    let bestRatio = Infinity;
    let leaveRow = -1;
    for (let i = 0; i < m; i++) {
      const a = tableau[i][enter];
      const bi = tableau[i][numCols - 1];
      if (a > 1e-12) {
        const ratio = bi / a;
        if (ratio < bestRatio) {
          bestRatio = ratio;
          leaveRow = i;
        }
      }
    }
    if (leaveRow === -1) {
      return { status: "unbounded", trace };
    }
    pivot(enter, leaveRow);
    // snapshot
    const x = Array(c.length + m).fill(0);
    for (let i = 0; i < m; i++) {
      if (basis[i] < x.length) x[basis[i]] = tableau[i][numCols - 1];
    }
    const val = tableau[tableau.length - 1][numCols - 1];
    trace.push({ iter, x: x.slice(), obj: val });
  }

  // assemble solution
  const x = Array(c.length + A.length).fill(0);
  const numColsFinal = tableau[0].length;
  for (let i = 0; i < A.length; i++) {
    if (basis[i] < x.length) x[basis[i]] = tableau[i][numColsFinal - 1];
  }
  const objVal = tableau[tableau.length - 1][numColsFinal - 1];
  return { status: "optimal", x: x.slice(0, c.length), obj: objVal, trace };
}

// Plot feasible region (2D only) for linear inequalities: A x <= b, x >= 0
function feasibleRegion2D(A, b) {
  // We'll sample grid and test constraints
  const xs = [];
  const ys = [];
  // find bounding box: crude
  let xmax = 5, ymax = 5;
  // attempt to find rough bounds from intersections of lines
  // sample grid
  const N = 80;
  const grid = [];
  for (let i = 0; i <= N; i++) {
    const x = (i / N) * xmax;
    for (let j = 0; j <= N; j++) {
      const y = (j / N) * ymax;
      let ok = true;
      for (let k = 0; k < A.length; k++) {
        const lhs = A[k][0] * x + A[k][1] * y;
        if (lhs > b[k] + 1e-6) { ok = false; break; }
      }
      if (ok) { xs.push(x); ys.push(y); grid.push({ x, y }); }
    }
  }
  return { xs, ys, grid, xmax, ymax };
}

// ---------------- Section 15.1 UI ----------------
function Section151() {
  // default LP: maximize 3x1 + 2x2 s.t. x1 + x2 <= 4; x1 <= 3; x2 <= 2; x>=0
  const [cText, setCText] = useState("3,2");
  const [Atext, setAtext] = useState("1 1\n1 0\n0 1");
  const [bText, setBText] = useState("4\n3\n2");

  const c = useMemo(() => cText.split(/[\s,]+/).map(Number).filter(n => !Number.isNaN(n)), [cText]);
  const A = useMemo(() => Atext.trim().split(/[\r\n]+/).map(r => r.trim().split(/[\s,]+/).map(Number)), [Atext]);
  const b = useMemo(() => bText.trim().split(/[\r\n]+/).map(r => Number(r)), [bText]);

  const feasible = useMemo(() => {
    if (A.length === 0 || c.length < 2 || A[0].length < 2) return null;
    // use only first 2 variables for visualization
    const A2 = A.map(row => [row[0] || 0, row[1] || 0]);
    const b2 = b.slice();
    return feasibleRegion2D(A2, b2);
  }, [A, b, c]);

  const simplexResult = useMemo(() => {
    try {
      if (c.length === 0 || A.length === 0) return null;
      return simplexSolve(c, A, b);
    } catch (e) {
      return null;
    }
  }, [c, A, b]);

  return (
    <section className="rounded-lg p-4" style={{ background: THEME.panelBg }}>
      <motion.div {...fadeUp}>
        <SectionHeader Icon={Columns} title="15.1 Linear Programming" subtitle="Simplex demo and 2D feasible region visualization" />

        <div className="rounded bg-zinc-800 border border-zinc-700 p-3 mb-3">
          <div className="grid grid-cols-1 md:grid-cols-3 gap-3">
            <div className="md:col-span-2">
              <label className="text-xs text-zinc-400">Objective coefficients c (comma separated)</label>
              <Input value={cText} onChange={(e) => setCText(e.target.value)} className="w-full bg-zinc-900 p-2 rounded text-zinc-100" />
              <div className="text-xs text-zinc-400 mt-2">Example: <code>3,2</code> for maximize 3x1 + 2x2</div>
            </div>
            <div>
              <div className="text-xs text-zinc-400">Notes</div>
              <div className="text-sm text-zinc-200 mt-1">This demo visualizes only the first two variables. For robust LP solving use a proper library.</div>
            </div>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-3 gap-3 mt-3">
            <div className="md:col-span-2">
              <label className="text-xs text-zinc-400">Constraints A (rows, space-separated)</label>
              <textarea rows={4} value={Atext} onChange={(e) => setAtext(e.target.value)} className="w-full bg-zinc-900 p-2 rounded text-zinc-100" />
            </div>
            <div>
              <label className="text-xs text-zinc-400">Right-hand side b (rows)</label>
              <textarea rows={4} value={bText} onChange={(e) => setBText(e.target.value)} className="w-full bg-zinc-900 p-2 rounded text-zinc-100" />
            </div>
          </div>
        </div>

        <div className="rounded border border-zinc-700 bg-zinc-900 p-3 mb-3">
          <div className="text-xs text-zinc-400 mb-2">Feasible region (2D projection)</div>
          {feasible ? (
            <Plot
              data={[
                {
                  x: feasible.xs,
                  y: feasible.ys,
                  mode: "markers",
                  marker: { size: 4 },
                  name: "feasible points"
                },
                ...(simplexResult && simplexResult.status === "optimal" ? [{
                  x: [simplexResult.x[0]],
                  y: [simplexResult.x[1]],
                  mode: "markers",
                  marker: { size: 12, symbol: "diamond" },
                  name: "optimal x"
                }] : [])
              ]}
              layout={{ autosize: true, height: 420, paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)", xaxis: { title: "x1" }, yaxis: { title: "x2" } }}
              useResizeHandler
              style={{ width: "100%" }}
            />
          ) : <div className="text-sm text-zinc-400">Provide at least 2 variables to visualize feasible region</div>}
        </div>

        <div className="rounded bg-zinc-800 border border-zinc-700 p-3">
          <div className="text-sm text-zinc-200 mb-2">Simplex result (educational)</div>
          {simplexResult ? (
            simplexResult.status === "optimal" ? (
              <div className="grid grid-cols-1 md:grid-cols-3 gap-3">
                <div className="rounded p-2 bg-zinc-900 border border-zinc-700">
                  <div className="text-xs text-zinc-400">Status</div>
                  <div className="text-sm text-zinc-100">Optimal</div>
                </div>
                <div className="rounded p-2 bg-zinc-900 border border-zinc-700">
                  <div className="text-xs text-zinc-400">Objective</div>
                  <div className="text-sm text-zinc-100">{Number(simplexResult.obj).toPrecision(6)}</div>
                </div>
                <div className="rounded p-2 bg-zinc-900 border border-zinc-700">
                  <div className="text-xs text-zinc-400">x</div>
                  <div className="text-sm text-zinc-100">[{simplexResult.x.map(v => Number(v).toFixed(4)).join(", ")}]</div>
                </div>
              </div>
            ) : (
              <div className="text-sm text-zinc-200">Status: {simplexResult.status}</div>
            )
          ) : <div className="text-sm text-zinc-400">Enter an LP and run simplex</div>}
        </div>

        <div className="mt-3">
          <Exercise
            title="Exercise 15.1.1 — Geography of optimum"
            body={`Change the objective coefficients and observe how the optimal vertex moves along the feasible polygon. Explain why optimal solutions lie at vertices for LP.`}
            solution={`Linear objective over a convex polytope attains maximum at an extreme point (vertex) because the level sets are hyperplanes and the feasible region is convex.`}
          />
        </div>
      </motion.div>
    </section>
  );
}

// ---------------- 15.2 Nonlinear Constrained Optimization ----------------
// We'll show: Lagrange multiplier demo for equality constraint (2D), simple quadratic penalty method for inequality conversion
function Section152() {
  const [fExpr, setFExpr] = useState(" (x0-1)*(x0-1) + 2*(x1-0.5)*(x1-0.5) ");
  const [gExpr, setGExpr] = useState(" x0 + x1 - 1 "); // constraint g(x)=0
  const fFun = useMemo(() => makeFun(fExpr, ["x0", "x1"]), [fExpr]);
  const gFun = useMemo(() => makeFun(gExpr, ["x0", "x1"]), [gExpr]);

  // Lagrange multiplier via solving gradient(f) = lambda * gradient(g) and g=0 (2 eqs + lambda)
  // We'll use numerical finite differences and small Newton solve for demonstration
  function gradApprox2(f, x, h = 1e-6) {
    const x0 = x.slice();
    const g = [];
    for (let i = 0; i < x.length; i++) {
      const xp = x0.slice(); xp[i] += h;
      const xm = x0.slice(); xm[i] -= h;
      g.push((f(xp) - f(xm)) / (2 * h));
    }
    return g;
  }

  function lagrangeSolve(f, g, x0 = [0.5, 0.5]) {
    // solve for (x0,x1,lambda) using simple Newton on system: grad f(x) - lambda grad g(x) = 0; g(x)=0
    let x = x0.slice();
    let lambda = 1;
    const trace = [];
    for (let iter = 0; iter < 60; iter++) {
      const gf = gradApprox2(f, x);
      const gg = gradApprox2(g, x);
      // system size 3: [H_f - lambda H_g  | -grad_g] [dx, dl] = -[grad_f - lambda grad_g, g]
      // Using finite-diff Jacobians is messy; do a simple damped gradient approach:
      const resid = subVec(gf.map((v, i) => v - lambda * gg[i]), [0, 0]);
      const gval = g(x);
      // update x by moving opposite to resid projected
      const step = 0.1;
      x = x.map((xi, i) => xi - step * resid[i]);
      // update lambda with simple rule to satisfy g(x)=0
      lambda = lambda + 0.5 * gval;
      trace.push({ iter, x: x.slice(), lambda, resid, g: gval });
      if (norm2(resid) < 1e-6 && Math.abs(gval) < 1e-6) break;
    }
    return { x, lambda, trace };
  }

  // Quadratic penalty method for inequality constraints h_i(x) <= 0
  function penaltyMethod(f, hList, x0 = [0, 0], mu0 = 1.0, beta = 10, tol = 1e-6, maxIter = 50) {
    let x = x0.slice();
    let mu = mu0;
    const trace = [];
    for (let k = 0; k < maxIter; k++) {
      // minimize F(x) = f(x) + 0.5*mu * sum(max(0, h_i(x))^2)
      // We'll do a few gradient descent steps (educational)
      for (let t = 0; t < 20; t++) {
        const grad = gradApprox2((xx) => {
          let val = f(xx);
          for (let hi of hList) {
            const hv = hi(xx);
            if (hv > 0) val += 0.5 * mu * hv * hv;
          }
          return val;
        }, x, 1e-6);
        const step = 0.05;
        x = x.map((xi, i) => xi - step * grad[i]);
      }
      // check constraints
      const violations = hList.map(hf => hf(x)).filter(v => v > 1e-6);
      trace.push({ k, x: x.slice(), violations: violations.slice(), mu });
      if (violations.length === 0) break;
      mu *= beta;
    }
    return { x, trace };
  }

  const lagResult = useMemo(() => {
    if (!fFun || !gFun) return null;
    try {
      return lagrangeSolve(fFun, gFun, [0.2, 0.8]);
    } catch (e) {
      return null;
    }
  }, [fFun, gFun]);

  const penaltyResult = useMemo(() => {
    if (!fFun) return null;
    try {
      // example: inequality constraint x0>=0, x1>=0 -> h1 = -x0, h2 = -x1 (so h <= 0 means x>=0)
      const hList = [
        (xx) => Math.max(0, -xx[0]), // violation positive if x0<0
        (xx) => Math.max(0, -xx[1]),
      ];
      return penaltyMethod((xx) => fFun(xx[0], xx[1]), hList, [ -0.5, 0.5 ], 1.0, 10, 1e-6, 20);
    } catch (e) {
      return null;
    }
  }, [fFun]);

  // prepare contour for plotting
  const surface = useMemo(() => {
    if (!fFun) return null;
    const lo = -2, hi = 2, N = 80;
    const xs = Array.from({ length: N }, (_, i) => lo + (hi - lo) * (i / (N - 1)));
    const ys = xs.slice();
    const z = Array.from({ length: N }, () => Array(N).fill(0));
    for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) z[i][j] = fFun(xs[j], ys[i]);
    return { xs, ys, z };
  }, [fFun]);

  return (
    <section className="rounded-lg p-4" style={{ background: THEME.panelBg }}>
      <motion.div {...fadeUp}>
        <SectionHeader Icon={Triangle} title="15.2 Nonlinear Constrained Optimization" subtitle="Lagrange multipliers demo and penalty methods" />

        <div className="rounded bg-zinc-800 border border-zinc-700 p-3 mb-3">
          <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
            <div>
              <label className="text-xs text-zinc-400">Objective f(x0,x1)</label>
              <Input value={fExpr} onChange={(e) => setFExpr(e.target.value)} className="w-full bg-zinc-900 p-2 rounded text-zinc-100" />
              <div className="text-xs text-zinc-400 mt-1">Example: <code>(x0-1)*(x0-1)+2*(x1-0.5)*(x1-0.5)</code></div>
            </div>
            <div>
              <label className="text-xs text-zinc-400">Equality constraint g(x0,x1) = 0</label>
              <Input value={gExpr} onChange={(e) => setGExpr(e.target.value)} className="w-full bg-zinc-900 p-2 rounded text-zinc-100" />
              <div className="text-xs text-zinc-400 mt-1">Example: <code>x0 + x1 - 1</code></div>
            </div>
          </div>
        </div>

        <div className="grid grid-cols-1 lg:grid-cols-3 gap-3 mb-3">
          <div className="lg:col-span-2 rounded border border-zinc-700 bg-zinc-900 p-3">
            <div className="text-xs text-zinc-400 mb-2">Contour & methods</div>
            {surface ? (
              <Plot
                data={[
                  { x: surface.xs, y: surface.ys, z: surface.z, type: "contour", contours: { coloring: "heatmap" } },
                  ...(lagResult ? [{ x: [lagResult.x[0]], y: [lagResult.x[1]], mode: "markers", name: "Lagrange solution", marker: { size: 10, symbol: "diamond" } }] : []),
                  ...(penaltyResult ? [{ x: [penaltyResult.x[0]], y: [penaltyResult.x[1]], mode: "markers", name: "Penalty solution", marker: { size: 10, symbol: "circle" } }] : []),
                ]}
                layout={{ autosize: true, height: 420, paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)" }}
                useResizeHandler
                style={{ width: "100%" }}
              />
            ) : <div className="text-sm text-zinc-400">Invalid objective</div>}
          </div>

          <div className="rounded border border-zinc-700 bg-zinc-800 p-3">
            <div className="text-xs text-zinc-400 mb-2">Method outputs</div>
            <div className="text-sm text-zinc-100">Lagrange approx: x = {lagResult ? `(${lagResult.x.map(v => Number(v).toFixed(4)).join(", ")})` : "—"}</div>
            <div className="text-sm text-zinc-100 mt-1">lambda ≈ {lagResult ? Number(lagResult.lambda).toPrecision(6) : "—"}</div>
            <div className="text-sm text-zinc-100 mt-2">Penalty method result: {penaltyResult ? `(${penaltyResult.x.map(v => Number(v).toFixed(4)).join(", ")})` : "—"}</div>
          </div>
        </div>

        <div>
         <Exercise
  title="Exercise 15.2.1 — KKT & multipliers"
  body={`For the equality-constrained problem minimize $f(x)$ subject to $g(x)=0$, write the Lagrangian $L(x,\\lambda)=f(x)+\\lambda g(x)$ and derive the KKT conditions. Use the demo to verify numerically.`}
  solution={`KKT conditions (for equality only): $\\nabla_x L(x,\\lambda)=0$ and $g(x)=0$. Solve the system of equations for $(x,\\lambda)$.`}
/>

        </div>
      </motion.div>
    </section>
  );
}

// ---------------- 15.3 Optimization with Software Packages ----------------
function Section153() {
  const snippet = `# Python (SciPy)
from scipy.optimize import linprog, minimize

# LP
res = linprog(c=-c, A_ub=A, b_ub=b, bounds=(0, None))  # maximize c^T x by minimizing -c^T x

# Nonlinear constrained (SLSQP)
cons = ({'type': 'eq', 'fun': lambda x: g(x)})
res2 = minimize(f, x0, method='SLSQP', constraints=cons)`;

  return (
    <section className="rounded-lg p-4" style={{ background: THEME.panelBg }}>
      <motion.div {...fadeUp}>
        <SectionHeader Icon={BookOpen} title="15.3 Optimization with Software Packages" subtitle="Practical recommendations and pseudocode" />

        <div className="rounded bg-zinc-800 border border-zinc-700 p-3 mb-3">
          <div className="text-sm text-zinc-200 mb-3">
            <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>{docs["15.3 Optimization with Software Packages"].join("\n\n")}</ReactMarkdown>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
            <div>
              <div className="text-xs text-zinc-400">When to use which</div>
              <ul className="list-disc pl-4 text-sm text-zinc-300">
                <li>LP: use industrial solvers (Gurobi, CPLEX) or open-source GLPK for large problems.</li>
                <li>Nonlinear constraints: IPOPT is a high-quality interior-point solver; SLSQP is suitable for smaller problems.</li>
                <li>Use automatic differentiation where possible for reliable derivatives (CasADi, JAX, autograd).</li>
              </ul>
            </div>
            <div>
              <div className="text-xs text-zinc-400">SciPy pseudocode</div>
              <pre className="bg-zinc-900 p-3 rounded text-zinc-200 text-sm" style={{ overflowX: "auto" }}>{snippet}</pre>
            </div>
          </div>
        </div>

        <Exercise
          title="Exercise 15.3.1 — Try a package"
          body={`Try solving a small LP using SciPy's linprog or GLPK. Compare runtime & solution accuracy. Report how many constraints and variables you tested.`}
          solution={`Expected: small LPs solve in milliseconds; accuracy depends on solver tolerances. For large LPs, industrial solvers are orders of magnitude faster.`}
        />
      </motion.div>
    </section>
  );
}

// ---------------- Problems ----------------
function SectionProblems() {
  const problems = [
    {
      title: "Problem 15.1 — LP geometry",
      body: `Given an LP with several inequalities, sketch the feasible polygon and identify corner points. Show which corner optimizes a given objective.`,
      solution: `Compute all intersections of constraint lines within bounds, evaluate objective at vertices, and pick max/min.`,
    },
    {
      title: "Problem 15.2 — KKT verification",
      body: `Solve a small nonlinear constrained problem using KKT conditions (analytically) and verify numerically using the demo page.`,
      solution: `Derive Lagrangian and solve the stationarity + primal feasibility + complementary slackness conditions; check numerically.`,
    },
    {
      title: "Problem 15.3 — Compare solvers",
      body: `Solve a constrained nonlinear problem with SciPy (SLSQP) and IPOPT (if available). Compare iterations, function evaluations and final objective.`,
      solution: `Document solver statistics; IPOPT often needs fewer iterations for large-scale problems due to second-order information.`,
    },
  ];

  return (
    <section className="rounded-lg p-4" style={{ background: THEME.panelBg }}>
      <motion.div {...fadeUp}>
        <SectionHeader Icon={Code} title="Problems (Chapter 15)" subtitle="Practice and comparison tasks" />
        <div className="grid grid-cols-1 gap-3">
          {problems.map((p, i) => <Exercise key={i} title={p.title} body={p.body} solution={p.solution} />)}
        </div>
      </motion.div>
    </section>
  );
}

// ---------------- Page assembly ----------------
export default function Chapter15() {
  return (
    <div className={`min-h-screen p-6 bg-zinc-950 text-zinc-100`}>
      <header className="mb-6">
        <div className="flex items-center justify-between">
          <div>
            <h1 className="text-3xl font-bold text-emerald-400" >Constrained Optimization</h1>
            <div className="text-sm text-zinc-400 mt-1">Linear programming, nonlinear constraints, and software recommendations</div>
          </div>
        </div>
      </header>

      <main className="space-y-6">
        <Section151 />
        <Section152 />
        <Section153 />
        <SectionProblems />
        <BottomBar/>
      </main>

      <footer className="mt-8 text-xs text-zinc-500">Tip: For production problems always rely on well-tested solvers and pay attention to scaling, conditioning and accurate derivatives.</footer>
    </div>
  );
}
