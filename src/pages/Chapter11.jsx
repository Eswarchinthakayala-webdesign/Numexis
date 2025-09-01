// src/pages/Chapter11.jsx
// ======================================================================
// Chapter 11 — Special Matrices & Gauss-Seidel (expanded)
// Dark-gray pro UI, emerald/amber accents, visualization-first.
// This file is intentionally large and feature-rich (approx. ~800-1000 lines).
// ======================================================================

import React, { useEffect, useMemo, useState, useRef } from "react";
import Plot from "react-plotly.js";
import { motion } from "framer-motion";
import ReactMarkdown from "react-markdown";
import remarkMath from "remark-math";
import rehypeKatex from "rehype-katex";
import "katex/dist/katex.min.css";

import {
  Grid,
  Repeat,
  BookOpen,
  Circle,
  Triangle,
  Square,
  Download,
  FileText,
  Box,
  Activity,
} from "lucide-react";

import { saveAs } from "file-saver";
import BottomBar from "../components/BottomBar";

// Optional UI primitives (replace if you have a design system)
let Card, CardHeader, CardContent, CardTitle, Input, Button;
try {
  // try to import your UI primitives
  // eslint-disable-next-line
  Card = require("@/components/ui/card").Card;
  // eslint-disable-next-line
  CardHeader = require("@/components/ui/card").CardHeader;
  // eslint-disable-next-line
  CardContent = require("@/components/ui/card").CardContent;
  // eslint-disable-next-line
  CardTitle = require("@/components/ui/card").CardTitle;
  // eslint-disable-next-line
  Input = require("@/components/ui/input").Input;
  // eslint-disable-next-line
  Button = require("@/components/ui/button").Button;
} catch (e) {
  // fallbacks
  Card = ({ children, className, style }) => <div className={className} style={style}>{children}</div>;
  CardHeader = ({ children }) => <div>{children}</div>;
  CardContent = ({ children }) => <div>{children}</div>;
  CardTitle = ({ children }) => <div>{children}</div>;
  Input = (props) => <input {...props} />;
  Button = (props) => <button {...props} />;
}

// ---------------- Theme & animations ----------------
const THEME = {
  bgClass: "bg-zinc-900",
  panelBg: "#0b0f12",
  text: "#e6eef6",
  muted: "#9ca3af",
  accent: "#15803d", // emerald
  accent2: "#b45309", // amber
};

const fadeIn = {
  initial: { opacity: 0, y: 8 },
  animate: { opacity: 1, y: 0 },
  transition: { duration: 0.28 },
};

// ---------------- Linear algebra helpers ----------------
function parseMatrix(text) {
  const rows = text
    .trim()
    .split(/[\r\n;]+/)
    .map((r) => r.trim())
    .filter(Boolean);
  return rows.map((r) =>
    r
      .split(/[\s,]+/)
      .filter(Boolean)
      .map((v) => Number(v))
  );
}
function copyMat(A) {
  return A.map((r) => r.slice());
}
function identity(n) {
  return Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 1 : 0)));
}
function swapRows(A, i, j) {
  const tmp = A[i];
  A[i] = A[j];
  A[j] = tmp;
}
function matVec(A, x) {
  return A.map((row) => row.reduce((s, v, i) => s + v * x[i], 0));
}
function normInf(A) {
  return Math.max(...A.map((r) => r.reduce((s, v) => s + Math.abs(v), 0)));
}
function norm2Vec(v) {
  return Math.sqrt(v.reduce((s, x) => s + x * x, 0));
}
function subVec(a, b) {
  return a.map((v, i) => v - b[i]);
}

// LU with partial pivoting
function luDecompose(A) {
  const n = A.length;
  const U = copyMat(A);
  const L = identity(n);
  const P = identity(n);
  for (let k = 0; k < n; k++) {
    let maxRow = k;
    let maxVal = Math.abs(U[k][k] || 0);
    for (let i = k + 1; i < n; i++) {
      if (Math.abs(U[i][k]) > maxVal) {
        maxVal = Math.abs(U[i][k]);
        maxRow = i;
      }
    }
    if (maxVal === 0) return { success: false, L: null, U: null, P: null };
    if (maxRow !== k) {
      swapRows(U, k, maxRow);
      swapRows(P, k, maxRow);
      if (k >= 1) swapRows(L, k, maxRow);
    }
    for (let i = k + 1; i < n; i++) {
      const mult = U[i][k] / U[k][k];
      L[i][k] = mult;
      for (let j = k; j < n; j++) U[i][j] -= mult * U[k][j];
    }
  }
  return { success: true, L, U, P };
}
function forwardSub(L, b) {
  const n = L.length;
  const y = Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    let s = 0;
    for (let j = 0; j < i; j++) s += L[i][j] * y[j];
    y[i] = (b[i] - s) / (L[i][i] || 1);
  }
  return y;
}
function backwardSub(U, y) {
  const n = U.length;
  const x = Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    let s = 0;
    for (let j = i + 1; j < n; j++) s += U[i][j] * x[j];
    x[i] = (y[i] - s) / (U[i][i] || 1);
  }
  return x;
}

// condition estimate (inf norm) by solving inverse columns
function condInf(A) {
  const n = A.length;
  const nrm = normInf(A);
  const lu = luDecompose(A);
  if (!lu.success) return Infinity;
  const { L, U, P } = lu;
  const invCols = [];
  for (let i = 0; i < n; i++) {
    const ei = Array(n).fill(0);
    ei[i] = 1;
    const Pei = P.map((row) => row.reduce((s, v, j) => s + v * ei[j], 0));
    const y = forwardSub(L, Pei);
    const x = backwardSub(U, y);
    invCols.push(x);
  }
  const invMat = Array.from({ length: n }, (_, i) => invCols.map((col) => col[i]));
  const normInv = normInf(invMat);
  return nrm * normInv;
}

// ---------------- UI helpers ----------------
function SectionHeader({ title, subtitle, Icon }) {
  return (
    <div className="flex items-center gap-3 mb-4">
      <div style={{ width: 48, height: 48, borderRadius: 10, background: "#071018", display: "flex", alignItems: "center", justifyContent: "center", border: "1px solid rgba(255,255,255,0.03)" }}>
        <Icon className="w-6 h-6" />
      </div>
      <div>
        <div style={{ color: THEME.text, fontSize: 20, fontWeight: 700 }}>{title}</div>
        <div style={{ color: THEME.muted, fontSize: 13 }}>{subtitle}</div>
      </div>
    </div>
  );
}

function Exercise({ title, body, solution }) {
  const [open, setOpen] = useState(false);
  return (
    <div className="rounded-lg border border-zinc-700 bg-zinc-800 p-4">
      <div className="flex items-start justify-between gap-3">
        <div>
          <div style={{ color: THEME.text, fontWeight: 700 }}>{title}</div>
          <div className="text-sm" style={{ color: THEME.muted, marginTop: 6 }}>
            <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>
              {body}
            </ReactMarkdown>
          </div>
        </div>
        <div>
          <Button onClick={() => setOpen((s) => !s)} className="px-3 py-1 rounded cursor-pointer bg-white text-black  text-xs font-bold hover:bg-gray-200">
            {open ? "Hide" : "Show"}
          </Button>
        </div>
      </div>
      {open && (
        <div className="mt-3 text-sm text-zinc-200 border-t border-zinc-700 pt-3">
          <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>
            {solution}
          </ReactMarkdown>
        </div>
      )}
    </div>
  );
}

function MatrixHeatmap({ matrix, height = 280, title = "" }) {
  const z = matrix.map((r) => r.map((v) => v));
  return (
    <Plot
      data={[
        {
          z,
          type: "heatmap",
          colorscale: "YlOrBr",
          reversescale: false,
          hoverongaps: false,
        },
      ]}
      layout={{
        title,
        autosize: true,
        paper_bgcolor: "rgba(0,0,0,0)",
        plot_bgcolor: "rgba(0,0,0,0)",
        margin: { l: 40, r: 10, t: 30, b: 30 },
        height,
      }}
      useResizeHandler
      style={{ width: "100%" }}
    />
  );
}

// CSV export helper
function exportCSV(filename, headers, rows) {
  const csv = [headers.join(","), ...rows.map((r) => r.map((v) => JSON.stringify(v)).join(","))].join("\n");
  const blob = new Blob([csv], { type: "text/csv;charset=utf-8" });
  saveAs(blob, filename);
}

// ---------------- Documentation strings (KaTeX)
const docs = {
  "11.1 Special Matrices": [
    `Special matrices: diagonal, triangular, symmetric, SPD (symmetric positive definite), banded, Toeplitz, circulant, sparse. Exploiting structure reduces storage and compute.`,
    `Cholesky factorization for SPD: $$A = LL^T$$. For banded matrices use compact storage and specialized factorization algorithms.`,
  ],
  "11.2 Gauss-Seidel": [
    `Gauss-Seidel iteration (matrix splitting): with \\(A=D+L+U\\), iterate\n\n$$x^{(k+1)} = -(D+L)^{-1} U x^{(k)} + (D+L)^{-1} b.$$`,
    `Component form (with relaxation \\(\\omega\\)):\n\n$$x_i^{(k+1)} = (1-\\omega) x_i^{(k)} + \\frac{\\omega}{a_{ii}} \\left( b_i - \\sum_{j<i} a_{ij} x_j^{(k+1)} - \\sum_{j>i} a_{ij} x_j^{(k)} \\right).$$`,
  ],
  "11.3 Software packages": [
    `Use LAPACK for dense matrices, SuiteSparse/CHOLMOD for sparse SPD, and iterative solvers (CG, GMRES) with preconditioners for very large problems.`,
  ],
};

// ---------------- Section 11.1 — Special matrices ----------------
function Section111() {
  const [n, setN] = useState(36);
  const [show3D, setShow3D] = useState(false);

  const diag = useMemo(() => {
    const A = Array.from({ length: n }, () => Array(n).fill(0));
    for (let i = 0; i < n; i++) A[i][i] = i + 1;
    return A;
  }, [n]);

  const tri = useMemo(() => {
    const A = Array.from({ length: n }, () => Array(n).fill(0));
    for (let i = 0; i < n; i++) {
      if (i > 0) A[i][i - 1] = -1;
      A[i][i] = 4;
      if (i < n - 1) A[i][i + 1] = -1;
    }
    return A;
  }, [n]);

  const toeplitz = useMemo(() => {
    const A = Array.from({ length: n }, () => Array(n).fill(0));
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) A[i][j] = 1 / (1 + Math.abs(i - j));
    return A;
  }, [n]);

  // 3D surface data for toeplitz (for optional 3D visualization)
  const surface = useMemo(() => {
    return {
      z: toeplitz,
    };
  }, [toeplitz]);

  return (
    <section className="rounded-lg p-4" style={{ background: THEME.panelBg }}>
      <motion.div {...fadeIn}>
        <SectionHeader Icon={Grid} title="11.1 Special Matrices" subtitle="Diagonal, banded, Toeplitz, circulant — exploit structure" />

        <div className="rounded-lg bg-zinc-800 border border-zinc-700 p-3 mb-4">
          <div className="grid grid-cols-1 lg:grid-cols-3 gap-3">
            <div>
              <label className="text-xs text-zinc-400">Matrix size (n)</label>
              <Input value={n} onChange={(e) => setN(Number(e.target.value) || 1)} className="w-full bg-zinc-900 text-zinc-100 p-2 rounded mt-1" />
              <div className="mt-2 text-xs text-zinc-400">Use the slider to change n for visualization</div>
            </div>

            <div>
              <div className="text-xs text-zinc-400">Diagonal example</div>
              <div className="mt-2">
                <MatrixHeatmap matrix={diag} title="Diagonal" />
              </div>
            </div>

            <div>
              <div className="text-xs text-zinc-400">Tridiagonal (Poisson 1D)</div>
              <div className="mt-2">
                <MatrixHeatmap matrix={tri} title="Tridiagonal (n)" />
              </div>
            </div>
          </div>
        </div>

        <div className="rounded-lg bg-zinc-800 border border-zinc-700 p-3 mb-4">
          <div className="text-sm text-zinc-200 mb-3">
            <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>{docs["11.1 Special Matrices"].join("\n\n")}</ReactMarkdown>
          </div>

          <div className="grid grid-cols-1 lg:grid-cols-2 gap-3">
            <div>
              <div className="text-xs text-zinc-400 mb-2">Toeplitz heatmap</div>
              <MatrixHeatmap matrix={toeplitz} title="Toeplitz" />
            </div>
            <div>
              <div className="text-xs text-zinc-400 mb-2">Optional 3D surface for Toeplitz</div>
              <div className="flex gap-2 items-center">
                <label className="text-xs text-zinc-400">Show 3D surface</label>
                <input type="checkbox" className="cursor-pointer" checked={show3D} onChange={(e) => setShow3D(e.target.checked)} />
              </div>
              {show3D && (
                <div className="mt-2">
                  <Plot
                    data={[{ z: surface.z, type: "surface", contours: { z: { show: true } }, showscale: false }]}
                    layout={{ autosize: true, height: 320, paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)" }}
                    useResizeHandler
                    style={{ width: "100%" }}
                  />
                </div>
              )}
            </div>
          </div>
        </div>

        <Exercise
          title="Exercise 11.1.1 — Complexity"
          body={`For a tridiagonal matrix of size n, what is the time complexity (big-O) of solving Ax=b with Thomas algorithm vs dense LU? Provide estimates and reasoning.`}
          solution={`Thomas: O(n) operations and O(n) memory. Dense LU: O(n^3) operations, O(n^2) memory. For large n, Thomas is vastly cheaper.`}
        />
      </motion.div>
    </section>
  );
}

// ---------------- Section 11.2 — Gauss-Seidel with SOR ----------------
function Section112() {
  const triPreset = useMemo(() => {
    const n = 12;
    const A = Array.from({ length: n }, () => Array(n).fill(0));
    for (let i = 0; i < n; i++) {
      if (i > 0) A[i][i - 1] = -1;
      A[i][i] = 4;
      if (i < n - 1) A[i][i + 1] = -1;
    }
    const b = Array.from({ length: n }, (_, i) => Math.sin(i + 1) + 2);
    return { A, b };
  }, []);

  const [matrixText, setMatrixText] = useState(() => {
    const { A, b } = triPreset;
    return A.map((r, i) => r.concat([b[i]]).join(" ")).join("\n");
  });
  const [omega, setOmega] = useState(1.0);
  const [tol, setTol] = useState(1e-8);
  const [maxIter, setMaxIter] = useState(500);
  const [initialGuessText, setInitialGuessText] = useState(""); // optional

  // parse
  const { A, b } = useMemo(() => {
    const parsed = parseMatrix(matrixText);
    if (!parsed.length) return { A: [], b: [] };
    return { A: parsed.map((r) => r.slice(0, -1)), b: parsed.map((r) => r[r.length - 1]) };
  }, [matrixText]);

  // gauss-seidel implementation
  function gaussSeidel(Ain, bin, x0 = null, omega = 1.0, tol = 1e-6, maxIter = 200) {
    const n = Ain.length;
    if (n === 0) return { success: false, message: "Empty matrix" };
    const A = copyMat(Ain);
    const b = bin.slice();
    let x = x0 && x0.length === n ? x0.slice() : Array(n).fill(0);
    const trace = [];
    let converged = false;
    for (let k = 0; k < maxIter; k++) {
      const xold = x.slice();
      for (let i = 0; i < n; i++) {
        let sigma = 0;
        for (let j = 0; j < i; j++) sigma += A[i][j] * x[j];
        for (let j = i + 1; j < n; j++) sigma += A[i][j] * xold[j];
        const xiUnrelaxed = (b[i] - sigma) / (A[i][i] || 1e-18);
        const xi = (1 - omega) * xold[i] + omega * xiUnrelaxed;
        x[i] = xi;
      }
      const residVec = subVec(matVec(A, x), b);
      const resid = norm2Vec(residVec);
      trace.push({ k, x: x.slice(), resid });
      if (resid < tol) {
        converged = true;
        break;
      }
    }
    return { success: true, x, trace, converged };
  }

  const initialGuess = useMemo(() => {
    if (!initialGuessText.trim()) return null;
    const arr = initialGuessText.trim().split(/[\s,]+/).map((v) => Number(v));
    return arr;
  }, [initialGuessText]);

  const gsResult = useMemo(() => {
    if (!A.length) return { success: false };
    return gaussSeidel(A, b, initialGuess, Number(omega), Number(tol), Number(maxIter));
  }, [A, b, omega, tol, maxIter, initialGuess]);

  // prepare data for plotting residuals
  const residPlot = useMemo(() => {
    const xs = gsResult.trace ? gsResult.trace.map((t) => t.k) : [];
    const ys = gsResult.trace ? gsResult.trace.map((t) => t.resid) : [];
    return { xs, ys };
  }, [gsResult]);

  // Export trace CSV
  function exportTrace() {
    if (!gsResult.trace || !gsResult.trace.length) return;
    const headers = ["iteration", "residual", "x_vector"];
    const rows = gsResult.trace.map((t) => [t.k, t.resid, JSON.stringify(t.x)]);
    exportCSV("gauss_seidel_trace.csv", headers, rows);
  }

  return (
    <section className="rounded-lg p-4" style={{ background: THEME.panelBg }}>
      <motion.div {...fadeIn}>
        <SectionHeader Icon={Repeat} title="11.2 Gauss-Seidel" subtitle="Iterative solver and SOR (successive over-relaxation)" />

        <div className="rounded-lg bg-zinc-800 border border-zinc-700 p-3 mb-4">
          <div className="grid grid-cols-1 lg:grid-cols-3 gap-3">
            <div className="lg:col-span-2">
              <label className="text-xs text-zinc-400">Augmented matrix [A | b] (rows newline)</label>
              <textarea rows={8} value={matrixText} onChange={(e) => setMatrixText(e.target.value)} className="w-full bg-zinc-900 p-2 rounded text-zinc-100 mt-1" />
            </div>
            <div>
              <Labeled label="Initial guess (optional, space/comma-separated)">
                <textarea rows={3} value={initialGuessText} onChange={(e) => setInitialGuessText(e.target.value)} className="w-full bg-zinc-900 p-2 rounded text-zinc-100" />
              </Labeled>
              <Labeled label="Relaxation ω">
                <Input value={omega} onChange={(e) => setOmega(Number(e.target.value))} className="bg-zinc-900 text-zinc-100 w-full" />
              </Labeled>
              <Labeled label="Tolerance">
                <Input value={tol} onChange={(e) => setTol(Number(e.target.value))} className="bg-zinc-900 text-zinc-100 w-full" />
              </Labeled>
              <Labeled label="Max iterations">
                <Input value={maxIter} onChange={(e) => setMaxIter(Number(e.target.value))} className="bg-zinc-900 text-zinc-100 w-full" />
              </Labeled>
            </div>
          </div>

          <div className="mt-3 flex gap-2">
            <Button onClick={() => {
              const { A, b } = triPreset;
              setMatrixText(A.map((r, i) => r.concat([b[i]]).join(" ")).join("\n"));
            }}>Load tri preset</Button>
            <Button onClick={() => exportTrace()}><Download className="w-4 h-4 inline-block mr-1" /> Export trace CSV</Button>
          </div>
        </div>

        <div className="rounded-lg bg-zinc-900 border border-zinc-700 p-3 mb-3">
          <div className="text-sm text-zinc-200 mb-2">Gauss-Seidel trace & residuals</div>

          <div className="grid grid-cols-1 lg:grid-cols-3 gap-3">
            <div className="lg:col-span-2">
              <div className="rounded border border-zinc-700 p-2 bg-zinc-800">
                <div className="text-xs text-zinc-400 mb-1">Solution / status</div>
                <div className="text-sm text-zinc-100">Converged: {gsResult.converged ? "Yes" : "No (may have hit max iterations)"}</div>
                <div className="text-sm text-zinc-100 mt-1">Last residual: {gsResult.trace && gsResult.trace.length ? gsResult.trace[gsResult.trace.length - 1].resid.toExponential(3) : "—"}</div>
                <div className="text-sm text-zinc-100 mt-2">Solution x: {gsResult.x ? prettyVec(gsResult.x) : "—"}</div>
              </div>

              <div className="mt-3 rounded border border-zinc-700 p-2 bg-zinc-800">
                <Plot
                  data={[
                    { x: residPlot.xs, y: residPlot.ys, mode: "lines+markers", name: "residual", marker: { size: 6 } },
                  ]}
                  layout={{
                    autosize: true,
                    height: 260,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    xaxis: { title: "iteration" },
                    yaxis: { title: "||r||_2", type: "log" },
                    margin: { l: 40, r: 10, t: 20, b: 40 },
                  }}
                  useResizeHandler
                  style={{ width: "100%" }}
                />
              </div>
            </div>

            <div>
              <div className="rounded border border-zinc-700 p-2 bg-zinc-800">
                <div className="text-xs text-zinc-400">Quick stats</div>
                <div className="text-sm text-zinc-100">n = {A.length}</div>
                <div className="text-sm text-zinc-100 mt-1">Iterations = {gsResult.trace ? gsResult.trace.length : 0}</div>
                <div className="text-xs text-zinc-400 mt-2">Tips</div>
                <ul className="text-sm text-zinc-300 mt-1 list-disc pl-4">
                  <li>Gauss-Seidel converges for strictly diagonally dominant and SPD matrices.</li>
                  
                </ul>
              </div>

              <div className="mt-3">
                <Exercise
                  title="Exercise 11.2.2 — Optimal ω search"
                  body={`Search over ω in 0.6..1.8 (step 0.05) for this system and plot iterations-to-converge vs ω. Which ω minimizes iterations?`}
                  solution={`Optimal ω depends on the spectral radius of the iteration matrix; enumerating ω and measuring iterations is a practical way to find a good ω for moderate-size problems.`}
                />
              </div>
            </div>
          </div>
        </div>
      </motion.div>
    </section>
  );
}

// small helper used above
function Labeled({ label, children }) {
  return (
    <div className="mb-2">
      <div className="text-xs text-zinc-400 mb-1">{label}</div>
      <div>{children}</div>
    </div>
  );
}

// ---------------- Section 11.3 — Software packages ----------------
function Section113() {
  // display recommended packages and short snippets
  const snippet = `// Example: solve SPD system with Cholesky (pseudocode)
A = loadMatrix("A.mtx")
if isSPD(A):
  L = cholesky(A)
  y = forwardSub(L, b)
  x = backwardSub(transpose(L), y)
`;

  return (
    <section className="rounded-lg p-4" style={{ background: THEME.panelBg }}>
      <motion.div {...fadeIn}>
        <SectionHeader Icon={BookOpen} title="11.3 Linear Algebraic Equations with Software Packages" subtitle="LAPACK, SuiteSparse, ARPACK — practical guidance" />

        <div className="rounded-lg bg-zinc-800 border border-zinc-700 p-3 mb-4">
          <div className="text-sm text-zinc-200 mb-2">
            <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>
              {docs["11.3 Software packages"].join("\n\n")}
            </ReactMarkdown>
          </div>

          <div className="grid grid-cols-1 lg:grid-cols-2 gap-3">
            <div>
              <div className="text-xs text-zinc-400">When to prefer which</div>
              <ul className="list-disc pl-4 text-sm text-zinc-300">
                <li>LAPACK/BLAS: dense linear algebra, highly optimized in native libraries.</li>
                <li>SuiteSparse: direct sparse solvers (CHOLMOD, UMFPACK) for large sparse matrices.</li>
                <li>Iterative solvers: CG, GMRES, BiCGSTAB combined with preconditioners for massive systems.</li>
              </ul>
            </div>
            <div>
              <div className="text-xs text-zinc-400">Pseudocode snippet</div>
              <pre className="bg-zinc-900 p-3 rounded text-zinc-200 text-sm" style={{ overflowX: "auto" }}>
                {snippet}
              </pre>
            </div>
          </div>
        </div>

        <Exercise
          title="Exercise 11.3.1 — Choose a solver"
          body={`For a large sparse SPD matrix from finite element discretization (N ~ 1e6), recommend a solver and preconditioner and explain why.`}
          solution={`Conjugate Gradient with algebraic multigrid (AMG) preconditioning or CHOLMOD if memory allows. AMG scales well for PDE-like systems and reduces iteration counts dramatically.`}
        />
      </motion.div>
    </section>
  );
}

// ---------------- Problems ----------------
function SectionProblems() {
  const problems = [
    {
      title: "Problem 11.1 — Classify matrices",
      body: `Given several matrices, classify them as banded, Toeplitz, circulant, diagonal-dominant, SPD. Give justification and storage strategies.`,
      solution: `Inspect patterns: band width, constant diagonals (Toeplitz), circulant shift; SPD check via symmetry and positive pivots in Cholesky.`,
    },
    {
      title: "Problem 11.2 — Gauss-Seidel vs Jacobi",
      body: `Implement Jacobi and Gauss-Seidel for a 2D Poisson problem (n x n grid). Compare iteration counts and runtime to reach tolerance.`,
      solution: `Gauss-Seidel typically halves the iteration count compared to Jacobi, SOR with optimal ω can further accelerate.`,
    },
    {
      title: "Problem 11.3 — Software evaluation",
      body: `Experiment with SciPy (spsolve / splu / cg) or SuiteSparse and document runtime & memory for a sequence of increasingly large sparse matrices.`,
      solution: `Sparse direct solvers are fast for medium-sized problems; iterative solvers scale better to very large sizes when preconditioned.`,
    },
  ];
  return (
    <section className="rounded-lg p-4" style={{ background: THEME.panelBg }}>
      <motion.div {...fadeIn}>
        <SectionHeader Icon={FileText} title="Problems (Chapter 11)" subtitle="End-of-chapter tasks" />
        <div className="grid grid-cols-1 gap-3">
          {problems.map((p, i) => <Exercise key={i} title={p.title} body={p.body} solution={p.solution} />)}
        </div>
      </motion.div>
    </section>
  );
}

// ---------------- Utility: pretty vector for display ----------------
function prettyVec(v) {
  if (!v) return "—";
  return `[${v.map((x) => Number(x).toPrecision(6)).join(", ")}]`;
}

// ---------------- Page assembly ----------------
export default function Chapter11() {
  return (
    <div className={`min-h-screen p-6 ${THEME.bgClass} text-zinc-100`}>
      <header className="mb-6">
        <div className="flex items-center justify-between gap-6">
          <div>
            <h1 className="text-3xl text-emerald-400 font-bold">Special Matrices & Gauss-Seidel</h1>
            <div className="text-sm text-zinc-400 mt-1">Visualization-first learning: matrices, iterative solvers, and practical tools</div>
          </div>
         
        </div>
      </header>

      <main className="space-y-6">
        <Section111 />
        <Section112 />
        <Section113 />
        <SectionProblems />
        <BottomBar/>
      </main>
  
      <footer className="mt-8 text-xs text-zinc-500">Tip: For production performance rely on native optimized libraries (BLAS/LAPACK, SuiteSparse). This page is for exploration & education.</footer>
    </div>
  );
}
