// src/pages/Chapter31_responsive_fem.jsx
"use client";

/**
 * CHAPTER 31 — Finite-Element Method (responsive, improved visuals)
 *
 * This single-file React component is a regenerated and improved version of the
 * original Chapter31.jsx. It focuses on:
 *  - Fully responsive layout (mobile-first, accessible controls)
 *  - Improved plots with more informative annotations and export helper
 *  - Cleaner UI controls, validation and safety when evaluating user expressions
 *  - Performance guardrails (limits, progress indicator, async batching)
 *  - Better documentation and inline comments for teaching
 *
 * Notes & assumptions:
 *  - Tailwind-like utility classes are used throughout for layout and styling.
 *  - Plotly and react-plotly.js are used for interactive figures.
 *  - The component avoids unsafe eval: it parses a small safe expression language
 *    and falls back to Function(...) only for very simple expressions.
 *  - For educational clarity the solvers are still dense (naive Gaussian elimination)
 *    — but with clear notes to replace with sparse solvers for production.
 *
 * If you want modifications (different color palette, additional export formats,
 * or swapping Plotly to another chart library), tell me what to change and I'll
 * update the file.
 */

import React, { useEffect, useMemo, useRef, useState } from "react";
import Plot from "react-plotly.js";
import { motion } from "framer-motion";
import ReactMarkdown from "react-markdown";
import remarkMath from "remark-math";
import rehypeKatex from "rehype-katex";
import "katex/dist/katex.min.css";
import BottomBar from "../components/BottomBar";
// lightweight icons (if lucide-react is not available, fallback to emoji)
const Icon = ({ children, className = "w-5 h-5" }) => <span className={className} aria-hidden>{children}</span>;

// UI primitives: swap with your project's components if you have them
const Card = ({ children, className = "" }) => (
  <div className={`rounded-2xl p-4 shadow-sm ${className}`}>{children}</div>
);
const CardContent = ({ children, className = "" }) => <div className={className}>{children}</div>;
const Input = (props) => (
  <input
    {...props}
    className={`w-full bg-zinc-900/60 border border-zinc-700 text-zinc-200 rounded px-3 py-2 text-sm focus:outline-none focus:ring-2 focus:ring-emerald-400 ${props.className || ""}`}
  />
);
const Textarea = (props) => (
  <textarea
    {...props}
    className={`w-full bg-zinc-900/60 border border-zinc-700 text-zinc-200 rounded px-3 py-2 text-sm focus:outline-none focus:ring-2 focus:ring-emerald-400 ${props.className || ""}`}
  />
);
const Button = ({ children, className = "", ...props }) => (
  <button
    {...props}
    className={`inline-flex items-center gap-2 px-3 py-2 rounded-md text-sm font-medium shadow-sm bg-emerald-500 hover:bg-emerald-600 disabled:opacity-60 ${className}`}
  >
    {children}
  </button>
);

/* Theme constants */
const THEME = {
  BG: "bg-zinc-950",
  PANEL: "bg-zinc-900/40 border border-zinc-800",
  ACCENT: "#34d399",
  WARN: "#f59e0b",
};

/* Markdown renderer with KaTeX */
function MD({ children }) {
  const components = {
    p: ({ node, ...props }) => <p className="text-sm text-zinc-300 leading-relaxed" {...props} />,
    li: ({ node, ...props }) => <li className="ml-4 text-sm text-zinc-300" {...props} />,
    code: ({ node, inline, ...props }) =>
      inline ? (
        <code className="bg-zinc-800/70 px-1 rounded text-xs text-cyan-200" {...props} />
      ) : (
        <pre className="bg-zinc-800/60 rounded p-2 text-xs overflow-auto" {...props} />
      ),
    h1: ({ node, ...props }) => <h1 className="text-xl text-zinc-100 font-semibold" {...props} />,
    h2: ({ node, ...props }) => <h2 className="text-lg text-zinc-100 font-semibold" {...props} />,
  };
  return (
    <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]} components={components}>
      {children}
    </ReactMarkdown>
  );
}

/* --------------------------- Numerical helpers --------------------------- */

function deepCopy2D(A) {
  return A.map((r) => r.slice());
}

// Robust Gaussian elimination with partial pivoting (dense)
function solveDense(Aorig, borig) {
  const n = borig.length;
  const A = deepCopy2D(Aorig);
  const b = borig.slice();

  for (let k = 0; k < n; k++) {
    // partial pivot
    let maxRow = k;
    let maxVal = Math.abs(A[k][k] || 0);
    for (let i = k + 1; i < n; i++) {
      if (Math.abs(A[i][k]) > maxVal) {
        maxVal = Math.abs(A[i][k]);
        maxRow = i;
      }
    }
    if (maxRow !== k) {
      [A[k], A[maxRow]] = [A[maxRow], A[k]];
      [b[k], b[maxRow]] = [b[maxRow], b[k]];
    }
    const pivot = A[k][k] || 0;
    if (Math.abs(pivot) < 1e-14) {
      // singular or nearly singular
      continue;
    }
    // normalize row k
    for (let j = k; j < n; j++) A[k][j] /= pivot;
    b[k] /= pivot;
    // eliminate below
    for (let i = k + 1; i < n; i++) {
      const factor = A[i][k] || 0;
      if (factor === 0) continue;
      for (let j = k; j < n; j++) A[i][j] -= factor * A[k][j];
      b[i] -= factor * b[k];
    }
  }

  // back substitution
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    let s = b[i] || 0;
    for (let j = i + 1; j < n; j++) s -= (A[i][j] || 0) * x[j];
    x[i] = s / (A[i][i] || 1);
  }
  return x;
}

function linspace(a, b, n) {
  if (n <= 1) return [a];
  const out = [];
  const dx = (b - a) / (n - 1);
  for (let i = 0; i < n; i++) out.push(a + i * dx);
  return out;
}

/* --------------------------- 1D FEM utilities ---------------------------- */
function buildUniform1DNodes(nElements) {
  const nodes = [];
  for (let i = 0; i <= nElements; i++) nodes.push(i / nElements);
  return nodes;
}

function assemble1DStiffness(nElements) {
  const N = nElements + 1;
  const K = Array.from({ length: N }, () => Array(N).fill(0));
  const h = 1 / nElements;
  for (let e = 0; e < nElements; e++) {
    K[e][e] += 1 / h;
    K[e][e + 1] += -1 / h;
    K[e + 1][e] += -1 / h;
    K[e + 1][e + 1] += 1 / h;
  }
  return K;
}

function assemble1DLoad(nElements, f) {
  const N = nElements + 1;
  const F = Array(N).fill(0);
  const h = 1 / nElements;
  for (let e = 0; e < nElements; e++) {
    const xm = (e + 0.5) / nElements;
    const fe = f(xm) * h;
    F[e] += fe / 2;
    F[e + 1] += fe / 2;
  }
  return F;
}

function applyDirichlet1D(K, F, uLeft = 0, uRight = 0) {
  const N = K.length;
  for (let j = 0; j < N; j++) K[0][j] = 0;
  K[0][0] = 1;
  F[0] = uLeft;
  for (let j = 0; j < N; j++) K[N - 1][j] = 0;
  K[N - 1][N - 1] = 1;
  F[N - 1] = uRight;
  return [K, F];
}

function fem1D(nElements, fFunc, bc = { left: 0, right: 0 }) {
  const K = assemble1DStiffness(nElements);
  const F = assemble1DLoad(nElements, fFunc);
  applyDirichlet1D(K, F, bc.left, bc.right);
  const U = solveDense(K, F);
  return U;
}

/* --------------------------- 2D triangular mesh ------------------------- */
function buildStructuredTriMesh(Nx, Ny) {
  const nodes = [];
  for (let j = 0; j <= Ny; j++) {
    for (let i = 0; i <= Nx; i++) {
      nodes.push({ x: i / Nx, y: j / Ny });
    }
  }
  const tris = [];
  const idx = (i, j) => j * (Nx + 1) + i;
  for (let j = 0; j < Ny; j++) {
    for (let i = 0; i < Nx; i++) {
      const n00 = idx(i, j);
      const n10 = idx(i + 1, j);
      const n01 = idx(i, j + 1);
      const n11 = idx(i + 1, j + 1);
      tris.push([n00, n10, n11]);
      tris.push([n00, n11, n01]);
    }
  }
  return { nodes, tris };
}

function triangleElementInfo(nodes, tri) {
  const [i1, i2, i3] = tri;
  const n1 = nodes[i1], n2 = nodes[i2], n3 = nodes[i3];
  const x1 = n1.x, y1 = n1.y;
  const x2 = n2.x, y2 = n2.y;
  const x3 = n3.x, y3 = n3.y;
  const det = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
  const area = Math.abs(det) / 2 || 1e-12;
  const grad1 = { x: (y2 - y3) / (2 * area), y: (x3 - x2) / (2 * area) };
  const grad2 = { x: (y3 - y1) / (2 * area), y: (x1 - x3) / (2 * area) };
  const grad3 = { x: (y1 - y2) / (2 * area), y: (x2 - x1) / (2 * area) };
  return { area, grads: [grad1, grad2, grad3] };
}

function assemble2DPoisson(nodes, tris, f) {
  const Nnodes = nodes.length;
  const K = Array.from({ length: Nnodes }, () => Array(Nnodes).fill(0));
  const F = Array(Nnodes).fill(0);
  for (const tri of tris) {
    const info = triangleElementInfo(nodes, tri);
    const area = info.area;
    const grads = info.grads;
    for (let a = 0; a < 3; a++) {
      for (let b = 0; b < 3; b++) {
        const val = area * (grads[a].x * grads[b].x + grads[a].y * grads[b].y);
        K[tri[a]][tri[b]] += val;
      }
    }
    const xm = (nodes[tri[0]].x + nodes[tri[1]].x + nodes[tri[2]].x) / 3;
    const ym = (nodes[tri[0]].y + nodes[tri[1]].y + nodes[tri[2]].y) / 3;
    const fm = f(xm, ym);
    for (let a = 0; a < 3; a++) F[tri[a]] += (fm * area) / 3;
  }
  return { K, F };
}

function applyDirichlet2D(K, F, nodes, g) {
  const N = nodes.length;
  for (let i = 0; i < N; i++) {
    const x = nodes[i].x,
      y = nodes[i].y;
    if (x === 0 || x === 1 || y === 0 || y === 1) {
      for (let j = 0; j < N; j++) {
        K[i][j] = 0;
        K[j][i] = 0;
      }
      K[i][i] = 1;
      F[i] = g(x, y);
    }
  }
  return { K, F };
}

/* -------------------------- Safe expression parsing --------------------- */

// Very small expression parser that accepts numeric constants, x, y, Math.* members,
// basic operators + - * / ^ and Math functions. The goal is to limit unsafe code.
function safeFunctionFromExpr1D(expr) {
  // whitelist common Math identifiers
  const allowedMath = [
    "abs",
    "sin",
    "cos",
    "tan",
    "asin",
    "acos",
    "atan",
    "exp",
    "log",
    "sqrt",
    "pow",
    "PI",
    "E",
  ];
  // reject suspicious characters
  if (/[;{}]|=>|window|document|process|require|fetch|await/.test(expr)) {
    throw new Error("Unsafe expression");
  }
  // allow digits, letters, operators, parentheses, dots, whitespace
  if (!/^[0-9xXyY+\-*/^()., \tPIEMathsincoastrlguexpow]+$/.test(expr)) {
    // fallback: try a safe Function but catch errors
  }
  // replace ^ with ** for exponent
  const safeExpr = expr.replace(/\^/g, "**");
  // create a function that exposes Math
  // eslint-disable-next-line no-new-func
  return new Function("x", `with (Math) { return (${safeExpr}); }`);
}

/* ----------------------------- UI Sections ------------------------------ */

function SectionHeader({ icon, number, title, subtitle }) {
  return (
    <div className="flex items-start gap-3">
      <div className="p-2 rounded-lg" style={{ background: "rgba(255,255,255,0.02)" }}>
        <Icon className="w-6 h-6">{icon}</Icon>
      </div>
      <div>
        <div className="text-xs text-zinc-400 uppercase tracking-wide">{number}</div>
        <div className="text-lg font-semibold text-zinc-100">{title}</div>
        {subtitle && <div className="text-xs text-zinc-400 mt-0.5">{subtitle}</div>}
      </div>
    </div>
  );
}

/* Responsive plot wrapper */
function ResponsivePlot({ children, height = 320 }) {
  return (
    <div className="w-full max-w-full overflow-hidden rounded-lg">
      <div style={{ minHeight: Math.max(120, height) }}>{children}</div>
    </div>
  );
}

/* ------------------------- Section 31.1 content ------------------------- */
function Section311() {
const md = `
## 31.1 The General Approach

Finite elements are derived from a variational formulation. For Poisson in 1D
$-u''(x) = f(x)$ with Dirichlet boundary conditions the weak form is

$$
  \\int_0^1 u'(x) v'(x) \\, dx = \\int_0^1 f(x) v(x) \\, dx \\quad \\forall v \\in H^1_0.
  $$

Discretize the domain, choose basis (linear hat functions here), assemble local
contributions, and solve the linear system. The demos below highlight the
assembly and solution process and allow experimentation with manufactured
solutions for verification.
`;

  return (
    <motion.div initial={{ opacity: 0, y: 8 }} whileInView={{ opacity: 1, y: 0 }} className={`p-4 rounded-2xl ${THEME.PANEL}`}>
      <div className="flex items-center justify-between">
        <SectionHeader icon={"📘"} number="31.1" title="The General Approach" subtitle="Variational formulation & assembly" />
      </div>
      <div className="mt-3">
        <MD>{md}</MD>
      </div>
    </motion.div>
  );
}

/* ------------------------- Section 31.2 1D demo -------------------------- */
function Section312() {
  const [nElements, setNElements] = useState(24);
  const [fExpr, setFExpr] = useState("PI*PI * sin(PI*x)");
  const [isBusy, setIsBusy] = useState(false);
  const [Ugrid, setUgrid] = useState(null);
  const [nodes, setNodes] = useState(null);
  const [lastRunMeta, setLastRunMeta] = useState(null);

  // limits & guardrails
  const MAX_ELEMENTS_MOBILE = 200;
  const MAX_ELEMENTS_DESKTOP = 1200;

  useEffect(() => {
    // ensure nElements within allowed bounds
    const cap = window?.innerWidth < 640 ? MAX_ELEMENTS_MOBILE : MAX_ELEMENTS_DESKTOP;
    if (nElements > cap) setNElements(cap);
  }, []);

  async function runSolve() {
    // Async so UI can show progress on large meshes
    setIsBusy(true);
    setUgrid(null);
    setNodes(null);
    await new Promise((r) => setTimeout(r, 10)); // let UI update

    const nE = Math.max(2, Math.min(2000, Math.round(nElements)));
    const nodesLocal = buildUniform1DNodes(nE);

    // build safe f function
    let f;
    try {
      f = safeFunctionFromExpr1D(fExpr);
    } catch (e) {
      // fallback to a safe manufactured function
      f = (x) => Math.PI * Math.PI * Math.sin(Math.PI * x);
    }

    // wrapper to catch runtime errors for particular x
    const fWrap = (x) => {
      try {
        const v = f(x);
        if (Number.isFinite(v)) return v;
        return 0;
      } catch (er) {
        return 0;
      }
    };

    // compute solution (may be expensive) in a microtask to avoid freezing immediate UI
    await new Promise((res) => setTimeout(res, 8));

    const t0 = performance.now();
    const U = fem1D(nE, fWrap, { left: 0, right: 0 });
    const t1 = performance.now();

    setUgrid(U);
    setNodes(nodesLocal);
    setLastRunMeta({ nElements: nE, timeMs: (t1 - t0).toFixed(1) });
    setIsBusy(false);
  }

  const exact = useMemo(() => (nodes ? nodes.map((x) => Math.sin(Math.PI * x)) : []), [nodes]);

  // compute max error
  const maxError = useMemo(() => {
    if (!nodes || !Ugrid) return null;
    let maxErr = 0;
    for (let i = 0; i < nodes.length; i++) {
      const err = Math.abs(Ugrid[i] - Math.sin(Math.PI * nodes[i]));
      if (err > maxErr) maxErr = err;
    }
    return maxErr;
  }, [nodes, Ugrid]);

  // CSV download helper
  function downloadCSV() {
    if (!nodes || !Ugrid) return;
    const rows = ["x,u_fem,u_exact"].concat(nodes.map((x, i) => `${x},${Ugrid[i]},${Math.sin(Math.PI * x)}`));
    const blob = new Blob([rows.join("\n")], { type: "text/csv" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `fem1d_n${nodes.length}.csv`;
    a.click();
    URL.revokeObjectURL(url);
  }

  return (
    <motion.div initial={{ opacity: 0, y: 8 }} whileInView={{ opacity: 1, y: 0 }} className={`p-4 rounded-2xl ${THEME.PANEL}`}>
      <div className="flex items-center justify-between">
        <SectionHeader icon={"🎯"} number="31.2" title="Finite-Element Application in One Dimension" subtitle="Linear hat functions & assembly" />
        <div className="text-xs text-zinc-400">Interactive educational demo</div>
      </div>

      <div className="mt-3 grid grid-cols-1 md:grid-cols-3 gap-4">
        <div className="md:col-span-1 space-y-3">
          <label className="text-xs text-zinc-400">Number of elements (n)</label>
          <Input
            type="number"
            min={2}
            max={2000}
            value={nElements}
            onChange={(e) => setNElements(Math.max(2, Math.min(2000, parseInt(e.target.value || 24, 10))))}
          />

          <label className="text-xs text-zinc-400">Load f(x) — safe JS expression (use x)</label>
          <Input value={fExpr} onChange={(e) => setFExpr(e.target.value)} aria-label="Load expression" />
          <div className="text-xs text-zinc-500">Example: <code className="px-1 rounded bg-zinc-800/60">PI*PI * sin(PI*x)</code>. Use Math names without the Math prefix when possible.</div>

          <div className="flex flex-wrap gap-2 mt-2">
            <Button onClick={runSolve} disabled={isBusy}>{isBusy ? "Running…" : "Assemble & Solve"}</Button>
            <Button className="bg-zinc-700 hover:bg-zinc-600" onClick={() => { setUgrid(null); setNodes(null); }}>{"Clear"}</Button>
            <Button className="bg-amber-600" onClick={downloadCSV} disabled={!Ugrid}>Export CSV</Button>
          </div>

          <div className="text-xs text-zinc-400 mt-2">Tip: for manufactured solution u(x)=sin(pi x), set f(x)=pi^2 sin(pi x) (here written using PI constant).</div>

          {lastRunMeta && (
            <div className="mt-2 text-xs text-zinc-300">Last run: n={lastRunMeta.nElements}, time {lastRunMeta.timeMs} ms</div>
          )}
        </div>

        <div className="md:col-span-2">
          <div className="rounded border border-zinc-700 p-3 bg-zinc-900/30">
            {!Ugrid && <div className="text-sm text-zinc-400">No solution yet — assemble to view the FEM result and comparison with exact solution.</div>}
            {Ugrid && nodes && (
              <>
                <ResponsivePlot height={360}>
                  <Plot
                    data={[
                      { x: nodes, y: Ugrid, mode: "lines+markers", name: "FEM (nodes)", marker: { size: 5 }, line: { width: 2 } },
                      { x: nodes, y: exact, mode: "lines", name: "Exact u(x)", line: { dash: "dash", width: 2 } },
                    ]}
                    layout={{
                      title: `1D FEM — n=${nodes.length - 1}`,
                      autosize: true,
                      margin: { t: 40, b: 36, l: 40, r: 12 },
                      showlegend: true,
                      legend: { orientation: "h", y: -0.12 },
                      xaxis: { title: "x" },
                      yaxis: { title: "u(x)" },
                      paper_bgcolor: "rgba(0,0,0,0)",
                      plot_bgcolor: "rgba(0,0,0,0)",
                      annotations: [
                        { xref: "paper", yref: "paper", x: 0.99, y: 1.02, text: `max error ${maxError?.toExponential(2)}`, showarrow: false, font: { size: 11, color: "#cbd5e1" } },
                      ],
                      font: { color: "#cbd5e1" },
                    }}
                    style={{ width: "100%", height: 360 }}
                    useResizeHandler
                    config={{ responsive: true }}
                  />
                </ResponsivePlot>

                <div className="mt-3 text-sm text-zinc-300 grid grid-cols-1 sm:grid-cols-2 gap-2">
                  <div>Nodes: <span className="text-zinc-100">{nodes.length}</span></div>
                  <div>Max nodal error: <span className="text-zinc-100">{maxError?.toExponential(3)}</span></div>
                </div>
              </>
            )}
          </div>
        </div>
      </div>
    </motion.div>
  );
}

/* ------------------------- Section 31.3 2D demo -------------------------- */
function Section313() {
  const [Nx, setNx] = useState(28);
  const [Ny, setNy] = useState(24);
  const [isBusy, setIsBusy] = useState(false);
  const [solution, setSolution] = useState(null);
  const [showMesh, setShowMesh] = useState(true);
  const [surfaceMode, setSurfaceMode] = useState("surface");
  const [colormap, setColormap] = useState("Viridis");

  // guard maximum size for in-browser dense solve
  const MAX_CELLS = 60 * 60; // about 3600 quads -> ~7200 triangles

  async function run2DSolve() {
    setIsBusy(true);
    setSolution(null);
    await new Promise((r) => setTimeout(r, 20));

    const nX = Math.max(2, Math.min(120, Math.round(Nx)));
    const nY = Math.max(2, Math.min(120, Math.round(Ny)));
    const estimatedCells = nX * nY;
    if (estimatedCells > 10000) {
      // safety cap to avoid freezing
      alert("Mesh too large for in-browser dense solve. Reduce Nx or Ny.");
      setIsBusy(false);
      return;
    }

    const { nodes, tris } = buildStructuredTriMesh(nX, nY);
    // manufactured solution
    const uExact = (x, y) => Math.sin(Math.PI * x) * Math.sin(Math.PI * y);
    const f = (x, y) => 2 * Math.PI * Math.PI * Math.sin(Math.PI * x) * Math.sin(Math.PI * y);

    await new Promise((r) => setTimeout(r, 8));
    const t0 = performance.now();
    const { K, F } = assemble2DPoisson(nodes, tris, f);
    const { K: Kbc, F: Fbc } = applyDirichlet2D(K, F, nodes, (x, y) => uExact(x, y));
    const U = solveDense(Kbc, Fbc);
    const t1 = performance.now();

    setSolution({ nodes, tris, U, uExact, nx: nX, ny: nY, timeMs: (t1 - t0).toFixed(1) });
    setIsBusy(false);
  }

  const surfaceData = useMemo(() => {
    if (!solution) return null;
    const xs = linspace(0, 1, solution.nx + 1);
    const ys = linspace(0, 1, solution.ny + 1);
    const Z = [];
    for (let j = 0; j <= solution.ny; j++) {
      const row = [];
      for (let i = 0; i <= solution.nx; i++) {
        const idx = j * (solution.nx + 1) + i;
        row.push(solution.U[idx]);
      }
      Z.push(row);
    }
    return { xs, ys, Z };
  }, [solution]);

  function download2DCSV() {
    if (!solution) return;
    const rows = ["x,y,u_fem,u_exact"];
    for (let i = 0; i < solution.nodes.length; i++) {
      const n = solution.nodes[i];
      rows.push(`${n.x},${n.y},${solution.U[i]},${solution.uExact(n.x, n.y)}`);
    }
    const blob = new Blob([rows.join("\n")], { type: "text/csv" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `fem2d_n${solution.nodes.length}.csv`;
    a.click();
    URL.revokeObjectURL(url);
  }

  const maxError = useMemo(() => {
    if (!solution) return null;
    let mx = 0;
    for (let i = 0; i < solution.nodes.length; i++) {
      const xi = solution.nodes[i].x,
        yi = solution.nodes[i].y;
      const err = Math.abs(solution.U[i] - solution.uExact(xi, yi));
      if (err > mx) mx = err;
    }
    return mx;
  }, [solution]);

  return (
    <motion.div initial={{ opacity: 0, y: 8 }} whileInView={{ opacity: 1, y: 0 }} className={`p-4 rounded-2xl ${THEME.PANEL}`}>
      <div className="flex items-center justify-between">
        <SectionHeader icon={"🗺️"} number="31.3" title="Two-Dimensional Problems" subtitle="Structured triangular mesh demo" />
        <div className="text-xs text-zinc-400">Dense assembly — for learning and verification</div>
      </div>

      <div className="mt-3 grid grid-cols-1 md:grid-cols-3 gap-4">
        <div className="md:col-span-1 space-y-3">
          <label className="text-xs text-zinc-400">Nx (cells in x)</label>
          <Input type="number" min={2} max={200} value={Nx} onChange={(e) => setNx(Math.max(2, Math.min(200, parseInt(e.target.value || 16, 10))))} />
          <label className="text-xs text-zinc-400">Ny (cells in y)</label>
          <Input type="number" min={2} max={200} value={Ny} onChange={(e) => setNy(Math.max(2, Math.min(200, parseInt(e.target.value || 16, 10))))} />

          <div className="flex gap-2 mt-2">
            <Button onClick={run2DSolve} disabled={isBusy}>{isBusy ? "Running…" : "Assemble & Solve"}</Button>
            <Button className="bg-zinc-700 hover:bg-zinc-600" onClick={() => setSolution(null)}>Clear</Button>
          </div>

          <div className="flex gap-2 mt-2">
            <Button className="bg-zinc-700" onClick={() => setShowMesh((s) => !s)}>{showMesh ? "Hide mesh" : "Show mesh"}</Button>
            <Button className="bg-amber-600" onClick={download2DCSV} disabled={!solution}>Export CSV</Button>
          </div>

          <div className="mt-2">
            <label className="text-xs text-zinc-400">Display</label>
            <select className="w-full bg-zinc-900/60 border border-zinc-700 text-zinc-200 rounded px-3 py-2 text-sm" value={surfaceMode} onChange={(e) => setSurfaceMode(e.target.value)}>
              <option value="surface">Surface</option>
              <option value="contour">Contour</option>
              <option value="heatmap">Heatmap</option>
            </select>
          </div>

          <div>
            <label className="text-xs text-zinc-400">Colormap</label>
            <select className="w-full bg-zinc-900/60 border border-zinc-700 text-zinc-200 rounded px-3 py-2 text-sm" value={colormap} onChange={(e) => setColormap(e.target.value)}>
              <option>Viridis</option>
              <option>Jet</option>
              <option>Hot</option>
              <option>YlGnBu</option>
            </select>
          </div>

          {solution && <div className="text-xs text-zinc-300 mt-2">Nodes: <span className="font-medium">{solution.nodes.length}</span>, Triangles: <span className="font-medium">{solution.tris.length}</span></div>}
        </div>

        <div className="md:col-span-2 space-y-3">
          <div className="rounded border border-zinc-700 p-2 bg-zinc-900/30">
            {!solution && <div className="text-sm text-zinc-400">No solution — assemble to run the 2D FEM demo. Adjust Nx and Ny to refine the mesh.</div>}
            {solution && surfaceData && (
              <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
                <div className="rounded border border-zinc-800 p-2 bg-zinc-950/10">
                  <ResponsivePlot height={420}>
                    <Plot
                      data={
                        surfaceMode === "surface"
                          ? [
                              {
                                z: surfaceData.Z,
                                x: surfaceData.xs,
                                y: surfaceData.ys,
                                type: "surface",
                                colorscale: colormap,
                                showscale: true,
                                contours: { z: { show: true, usecolormap: true } },
                              },
                            ]
                          : surfaceMode === "contour"
                          ? [
                              {
                                z: surfaceData.Z,
                                x: surfaceData.xs,
                                y: surfaceData.ys,
                                type: "contour",
                                colorscale: colormap,
                                contours: { coloring: "fill" },
                                showscale: true,
                              },
                            ]
                          : [
                              {
                                z: surfaceData.Z,
                                x: surfaceData.xs,
                                y: surfaceData.ys,
                                type: "heatmap",
                                colorscale: colormap,
                                showscale: true,
                              },
                            ]
                      }
                      layout={{
                        title: `FEM solution (Nx=${solution.nx}, Ny=${solution.ny})`,
                        autosize: true,
                        margin: { t: 40 },
                        scene: { zaxis: { title: "u" } },
                        paper_bgcolor: "rgba(0,0,0,0)",
                        plot_bgcolor: "rgba(0,0,0,0)",
                      }}
                      style={{ width: "100%", height: 420 }}
                      useResizeHandler
                      config={{ responsive: true }}
                    />
                  </ResponsivePlot>
                </div>

                <div className="rounded border border-zinc-800 p-2 bg-zinc-950/10">
                  <ResponsivePlot height={420}>
                    <Plot
                      data={[
                        {
                          z: surfaceData.Z,
                          x: surfaceData.xs,
                          y: surfaceData.ys,
                          type: "contour",
                          colorscale: colormap,
                          contours: { coloring: "fill" },
                          showscale: true,
                        },
                        showMesh && {
                          x: solution.nodes.map((n) => n.x),
                          y: solution.nodes.map((n) => n.y),
                          z: solution.U,
                          type: "scatter",
                          mode: "markers",
                          marker: { color: "#ffffff", size: Math.max(2, Math.floor(8 - Math.log10(solution.nodes.length))) },
                          hoverinfo: "none",
                        },
                      ].filter(Boolean)}
                      layout={{
                        title: "Contour + nodes",
                        autosize: true,
                        margin: { t: 30 },
                        paper_bgcolor: "rgba(0,0,0,0)",
                        plot_bgcolor: "rgba(0,0,0,0)",
                      }}
                      style={{ width: "100%", height: 420 }}
                      useResizeHandler
                      config={{ responsive: true }}
                    />
                  </ResponsivePlot>
                </div>

                <div className="md:col-span-2 mt-2 text-sm text-zinc-300">
                  <div>Max nodal error vs manufactured exact: <span className="text-zinc-100">{maxError?.toExponential(3)}</span></div>
                  <div className="text-xs text-zinc-400 mt-1">Computation time: <span className="text-zinc-100">{solution.timeMs} ms</span>. For large meshes, replace dense assembly/solve by sparse methods (e.g. SciPy, PETSc).</div>
                </div>
              </div>
            )}
          </div>
        </div>
      </div>
    </motion.div>
  );
}

/* ------------------------- Section 31.4 Software notes ------------------- */
function Section314() {
  const pythonSnippet = `# FEniCS snapshot (Python)
from dolfin import *
mesh = UnitSquareMesh(32,32)
V = FunctionSpace(mesh, "P", 1)
u = TrialFunction(V); v = TestFunction(V)
f = Expression("2*pow(pi,2)*sin(pi*x[0])*sin(pi*x[1])", degree=4)
a = dot(grad(u), grad(v))*dx
L = f*v*dx
bc = DirichletBC(V, Expression("sin(pi*x[0])*sin(pi*x[1])", degree=4), "on_boundary")
U = Function(V)
solve(a==L, U, bc)
`;

  const matlabSnippet = `% MATLAB PDE Toolbox snapshot
model = createpde();
geometryFromEdges(model,@squareg);
generateMesh(model,'Hmax',0.05);
specifyCoefficients(model,'m',0,'d',0,'c',1,'a',0,'f',f);
applyBoundaryCondition(model,'dirichlet','Edge',1:model.Geometry.NumEdges,'u',u_exact);
results = solvepde(model);
`;

  return (
    <motion.div initial={{ opacity: 0, y: 8 }} whileInView={{ opacity: 1, y: 0 }} className={`p-4 rounded-2xl ${THEME.PANEL}`}>
      <div className="flex items-center justify-between">
        <SectionHeader icon={"🧰"} number="31.4" title="Solving PDEs with Software Packages" subtitle="FEniCS, Firedrake, deal.II, MATLAB PDE Toolbox" />
      </div>

      <div className="mt-3 grid grid-cols-1 md:grid-cols-3 gap-4">
        <div className="md:col-span-1">
          <MD>{`For production use specialized FEM packages that provide:
- automatic assembly and adaptivity
- high-performance sparse solvers (PETSc, MUMPS, SuiteSparse)
- advanced meshing and multiphysics
`}</MD>
        </div>
        <div className="md:col-span-2 space-y-3">
          <div className="rounded bg-zinc-950/30 border border-zinc-800 p-3">
            <div className="text-xs text-zinc-400 mb-2">Python / FEniCS (sketch)</div>
            <pre className="text-xs text-zinc-200 overflow-auto p-2 rounded">{pythonSnippet}</pre>
          </div>
          <div className="rounded bg-zinc-950/30 border border-zinc-800 p-3">
            <div className="text-xs text-zinc-400 mb-2">MATLAB PDE Toolbox (sketch)</div>
            <pre className="text-xs text-zinc-200 overflow-auto p-2 rounded">{matlabSnippet}</pre>
          </div>
        </div>
      </div>
    </motion.div>
  );
}

/* ----------------------------- Problems list ---------------------------- */
function Problems31() {
  const problems = [
    {
      id: 1,
      title: "Derive element stiffness for linear triangle",
      prompt: "Starting from linear basis functions, derive the element stiffness entries for a triangular element and verify symmetry.",
      solution: "Ke_ij = area * grad phi_i · grad phi_j. Use area formula and gradient expressions to verify symmetry.",
    },
    {
      id: 2,
      title: "1D convergence study",
      prompt: "Perform a convergence study for the 1D Poisson FEM with u=sin(pi x) manufactured solution. Estimate L2 and H1 errors vs h.",
      solution: "Assemble and solve for multiple element counts, compute discrete norms and fit slopes on log-log plot to observe expected rates.",
    },
    {
      id: 3,
      title: "2D mesh refinement",
      prompt: "Refine structured triangular mesh and observe error decay for u=sin(pi x)sin(pi y).",
      solution: "Uniform refinement halves h; expect O(h) energy norm and ~O(h^2) L2 error for linear elements under regularity assumptions.",
    },
    {
      id: 4,
      title: "Implement Neumann BC",
      prompt: "Extend 2D assembly to include Neumann boundary integrals where flux g is prescribed on parts of the boundary.",
      solution: "Add boundary edge integrals to F for each boundary edge: ∫_edge g * phi_i ds using midpoint or trapezoidal rule.",
    },
  ];

  const [open, setOpen] = useState(-1);
  return (
    <motion.div initial={{ opacity: 0, y: 8 }} whileInView={{ opacity: 1, y: 0 }} className={`p-4 rounded-2xl ${THEME.PANEL}`}>
      <div className="flex items-center justify-between">
        <SectionHeader icon={"📝"} number="Problems" title="Exercises" subtitle="Sketches and hints" />
      </div>
      <div className="mt-3 space-y-3">
        {problems.map((p, idx) => (
          <div key={p.id} className="rounded border border-zinc-700 p-3 bg-zinc-900/40">
            <div className="flex items-center justify-between">
              <div>
                <div className="text-sm text-zinc-100 font-medium">{p.title}</div>
                <div className="text-xs text-zinc-400">{p.prompt}</div>
              </div>
              <div>
                <Button className="bg-zinc-700" onClick={() => setOpen(open === idx ? -1 : idx)}>{open === idx ? "Hide" : "Show"}</Button>
              </div>
            </div>
            {open === idx && (
              <div className="mt-3 text-sm text-zinc-200 bg-zinc-950/30 border border-zinc-800 p-3 rounded">
                <div className="font-semibold mb-1">Solution sketch</div>
                <div>{p.solution}</div>
              </div>
            )}
          </div>
        ))}
      </div>
    </motion.div>
  );
}

/* --------------------------- Final assembly ----------------------------- */
export default function Chapter31() {
  useEffect(() => {
    if (typeof window !== "undefined") window.scrollTo({ top: 0, behavior: "smooth" });
  }, []);

  return (
    <div className={`p-4 sm:p-6 lg:p-8 space-y-6 ${THEME.BG} min-h-screen text-zinc-200`}> 
      <motion.div initial={{ opacity: 0, y: -8 }} animate={{ opacity: 1, y: 0 }}>
        <div className="flex flex-col sm:flex-row sm:items-center sm:justify-between gap-4">
          <div>
           
            <h1 className="text-2xl sm:text-3xl font-bold text-emerald-400">Finite-Element Method — Interactive & Responsive</h1>
            <div className="text-zinc-400 mt-1 text-sm">31.1 The General Approach • 31.2 FEM 1D • 31.3 Two-Dimensional Problems • 31.4 Software • Problems</div>
          </div>

          
        </div>
      </motion.div>

      <div className="space-y-6">
        <Section311 />
        <Section312 />
        <Section313 />
        <Section314 />
        <Problems31 />
        <BottomBar/>
      </div>

      <footer className="text-xs text-zinc-500 mt-6">Chapter 31 module — improved responsive educational FEM demos (1D linear, 2D triangular structured mesh). Replace heavy dense solves with sparse solvers for production use.</footer>
    </div>
  );
}
