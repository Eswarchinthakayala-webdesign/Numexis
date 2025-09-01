// src/pages/Chapter9.jsx
// ======================================================================
// Chapter 9 — Gauss Elimination
// Stacked section layout, dark-gray pro UI, KaTeX docs, Plotly visualizations
// ======================================================================

import React, { useMemo, useState, useEffect } from "react";
import Plot from "react-plotly.js";
import { motion } from "framer-motion";
import ReactMarkdown from "react-markdown";
import remarkMath from "remark-math";
import rehypeKatex from "rehype-katex";
import "katex/dist/katex.min.css";
import BottomBar from "../components/BottomBar";
import {
  BookOpen,
  List,
  ChevronRight,
  FilePlus,
  Zap,
  Sliders,
  CornerUpLeft,
  RotateCcw,
  Grid,
  BoxSelect,
  Code,
  Menu,
} from "lucide-react";

// If you don't have these UI components available, replace with your own.
import { Card, CardHeader, CardContent, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Button } from "@/components/ui/button";

// ----------------------------------------------------------------------
// Theme and small animations
// ----------------------------------------------------------------------
const theme = {
  bg: "bg-zinc-900",
panelBg: "#18181b", // zinc-900

  text: "#e6eef6",
  accent: "#22d3ee",
  accent2: "#34d399",
  muted: "#94a3b8",
};

const fadeIn = {
  initial: { opacity: 0, y: 8 },
  animate: { opacity: 1, y: 0 },
  transition: { duration: 0.28 },
};

// ----------------------------------------------------------------------
// Linear algebra helpers (self-contained)
// ----------------------------------------------------------------------
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

// ----------------------------------------------------------------------
// Small UI primitives used across sections
// ----------------------------------------------------------------------
function Exercise({ title, body, solution }) {
  const [open, setOpen] = useState(false);
  return (
    <div className="rounded-lg border border-zinc-700 bg-zinc-800 p-4">
      <div className="flex items-start justify-between gap-3">
        <div>
          <div className="text-sm font-semibold text-zinc-100">{title}</div>
          <div className="text-xs text-zinc-300 mt-1">
            <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]}>
              {body}
            </ReactMarkdown>
          </div>
        </div>
        <div>
          <Button onClick={() => setOpen((s) => !s)} className="cursor-pointer rounded bg-white text-black hover:bg-gray-300 text-xs border border-zinc-700">
            {open ? "Hide solution" : "Show solution"}
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

function Labeled({ label, children }) {
  return (
    <div className="mb-2">
      <div className="text-xs text-zinc-400 mb-1">{label}</div>
      {children}
    </div>
  );
}

function Summary({ items = [] }) {
  return (
    <div className="grid grid-cols-1 sm:grid-cols-3 gap-3 mt-2">
      {items.map((it, i) => (
        <div key={i} className="rounded-lg border border-zinc-700 p-3 bg-zinc-800">
          <div className="text-xs text-zinc-400">{it.label}</div>
          <div className="text-sm text-zinc-100 mt-1">{it.value}</div>
        </div>
      ))}
    </div>
  );
}

function MatrixHeatmap({ matrix, title = "" }) {
  const z = matrix.map((r) => r.map((v) => v));
  return (
    <Plot
      data={[{ z, type: "heatmap", colorscale: "RdBu", reversescale: true }]}
      layout={{
        title,
        autosize: true,
        paper_bgcolor: "rgba(0,0,0,0)",
        plot_bgcolor: "rgba(0,0,0,0)",
        margin: { l: 40, r: 10, t: 30, b: 30 },
        height: 300,
      }}
      useResizeHandler
      style={{ width: "100%" }}
    />
  );
}

// ----------------------------------------------------------------------
// Documentation strings (KaTeX-ready)
// ----------------------------------------------------------------------
const docsText = {
  "9.1 Solving Small Numbers of Equations": [
    `Direct solution for small systems often uses substitution or elimination. Example system:\n\n$$\\begin{cases}2x+y-z=8\\\\ -3x-y+2z=-11\\\\ -2x+y+2z=-3\\end{cases}$$`,
    `Form: \\(A\\mathbf{x}=\\mathbf{b}\\). Use elimination to form \\(U\\mathbf{x}=\\mathbf{c}\\) and back-substitute.`,
  ],
  "9.2 Naive Gauss Elimination": [
    `Elimination uses multipliers $$m_{i,k}=\\frac{a_{i,k}}{a_{k,k}}$$ and row operation \\(R_i\\leftarrow R_i - m_{i,k}R_k\\).`,
  ],
  "9.3 Pitfalls of Elimination Methods": [
    `Pivoting needed to avoid division by small pivots. Ill-conditioning amplifies relative errors: $$\\frac{\\|\\delta x\\|}{\\|x\\|}\\le\\kappa(A)\\frac{\\|\\delta A\\|}{\\|A\\|}.$$`,
  ],
  "9.4 Techniques for Improving Solutions": [
    `Use partial pivoting, scaled pivoting, LU with pivoting, and iterative refinement (compute residual and correct).`,
  ],
  "9.5 Complex Systems": [
    `Sparse/banded systems (tridiagonal) allow specialized algorithms (Thomas algorithm) with O(n) cost.`,
  ],
  "9.6 Nonlinear Systems of Equations": [
    `Newton for systems: \\(x_{k+1}=x_k - J(x_k)^{-1}F(x_k)\\). A linear solve is required at each iteration.`,
  ],
  "9.7 Gauss-Jordan": [
    `Gauss-Jordan reduces to RREF; useful for demonstration, pseudoinverse, but expensive for repeated solves.`,
  ],
  "9.8 Summary": [
    `Direct elimination with pivoting and iterative refinement is standard for dense systems; use sparse/iterative solvers for very large problems.`,
  ],
};

// ----------------------------------------------------------------------
// SECTION 9.1: Solving small systems
// ----------------------------------------------------------------------
function Section91() {
  const [matrixText, setMatrixText] = useState("2 1 -1 8\n-3 -1 2 -11\n-2 1 2 -3");
  const Aaug = useMemo(() => parseMatrix(matrixText), [matrixText]);
  const A = useMemo(() => Aaug.map((r) => r.slice(0, -1)), [Aaug]);

  const trace = useMemo(() => {
    const M = copyMat(Aaug);
    const steps = [{ mat: copyMat(M), desc: "Initial augmented matrix" }];
    const n = M.length;
    for (let k = 0; k < n - 1; k++) {
      const pivot = M[k][k];
      if (Math.abs(pivot) < 1e-15) {
        steps.push({ mat: copyMat(M), desc: `Zero pivot at column ${k}` });
        break;
      }
      for (let i = k + 1; i < n; i++) {
        const mult = M[i][k] / pivot;
        for (let j = k; j < M[0].length; j++) M[i][j] -= mult * M[k][j];
      }
      steps.push({ mat: copyMat(M), desc: `Eliminated below pivot at col ${k}` });
    }
    // back-sub
    const U = M.map((r) => r.slice(0, -1));
    const bb = M.map((r) => r[r.length - 1]);
    const x = Array(U.length).fill(0);
    for (let i = U.length - 1; i >= 0; i--) {
      let s = 0;
      for (let j = i + 1; j < U.length; j++) s += U[i][j] * x[j];
      x[i] = (bb[i] - s) / (U[i][i] || 1);
    }
    steps.push({ mat: copyMat(M), desc: "Back substitution (naive)", sol: x });
    return steps;
  }, [Aaug]);

  return (
    <section className="rounded-lg p-4" style={{ background: theme.panelBg }}>
      <motion.div {...fadeIn}>
        <div className="flex items-center gap-3 mb-3">
          <FilePlus className="w-6 h-6 text-cyan-300" />
          <div>
            <h2 className="text-lg font-semibold text-zinc-100">9.1 Solving Small Numbers of Equations</h2>
            <div className="text-xs text-zinc-400">Direct elimination demo with an augmented matrix example</div>
          </div>
        </div>

        <div className="mb-3">
          <div className="rounded-lg border border-zinc-700 p-3 bg-zinc-800">
            <div className="grid grid-cols-1 md:grid-cols-3 gap-3">
              <div className="md:col-span-2">
                <Labeled label="Augmented matrix [A | b] (rows newline)">
                  <textarea rows={4} value={matrixText} onChange={(e) => setMatrixText(e.target.value)} className="w-full bg-zinc-900 text-zinc-100 p-2 rounded" />
                </Labeled>
              </div>
              <div>
                <Summary items={[{ label: "n", value: A.length }, { label: "cond_inf", value: Number.isFinite(condInf(A)) ? condInf(A).toFixed(4) : "∞" }]} />
              </div>
            </div>
          </div>
        </div>

        <div className="space-y-3">
          {trace.map((s, i) => (
            <div key={i} className="rounded border border-zinc-700 p-3 bg-zinc-900">
              <div className="text-xs text-zinc-400 mb-2">{s.desc}</div>
              <div className="overflow-auto mb-2">
                <table className="table-auto text-sm">
                  <tbody>
                    {s.mat.map((r, ri) => (
                      <tr key={ri}>
                        {r.map((c, ci) => (
                          <td key={ci} className="px-2 py-1" style={{ minWidth: 56, color: ci === r.length - 1 ? "#93c5fd" : "#e5e7eb" }}>
                            {Number.isFinite(c) ? Number(c).toFixed(6) : c}
                          </td>
                        ))}
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
              {s.sol && <div className="text-sm">Solution: <span className="font-mono">[{s.sol.map((v) => v.toFixed(6)).join(", ")}]</span></div>}
            </div>
          ))}
        </div>

        <div className="mt-4">
          <Exercise
            title="Exercise 9.1.1 — Manual elimination"
            body={`Perform elimination manually on the provided matrix and verify the computed solution.`}
            solution={`Follow the elimination steps shown above. Final solution is the vector displayed in the last step.`}
          />
        </div>
      </motion.div>
    </section>
  );
}

// ----------------------------------------------------------------------
// SECTION 9.2: Naive Gauss Elimination interactive stepper
// ----------------------------------------------------------------------
function Section92() {
  const [matrixText, setMatrixText] = useState("1e-8 1 1 3\n1 1 1 6\n2 1 -1 3");
  const Aaug = useMemo(() => parseMatrix(matrixText), [matrixText]);
  const [trace, setTrace] = useState([{ mat: copyMat(Aaug), step: "Initial matrix" }]);
  const [usePivot, setUsePivot] = useState(true);

  useEffect(() => {
    setTrace([{ mat: copyMat(Aaug), step: "Initial matrix" }]);
  }, [matrixText]);

  function nextStep() {
    const curr = copyMat(trace[trace.length - 1].mat);
    const n = curr.length;
    // find pivot column k
    let k;
    for (k = 0; k < n - 1; k++) {
      let hasNonZero = false;
      for (let r = k + 1; r < n; r++) if (Math.abs(curr[r][k]) > 1e-12) { hasNonZero = true; break; }
      if (hasNonZero) break;
    }
    if (k >= n - 1) {
      setTrace((t) => [...t, { mat: curr, step: "Finished elimination" }]);
      return;
    }
    if (usePivot) {
      let maxRow = k;
      let maxVal = Math.abs(curr[k][k] || 0);
      for (let r = k + 1; r < n; r++) {
        if (Math.abs(curr[r][k]) > maxVal) {
          maxVal = Math.abs(curr[r][k]);
          maxRow = r;
        }
      }
      if (maxRow !== k) {
        swapRows(curr, k, maxRow);
        setTrace((t) => [...t, { mat: copyMat(curr), step: `Pivot swap rows ${k} <-> ${maxRow}` }]);
        return;
      }
    }
    const pivot = curr[k][k];
    if (Math.abs(pivot) < 1e-14) {
      setTrace((t) => [...t, { mat: copyMat(curr), step: `Tiny pivot at col ${k}` }]);
      return;
    }
    for (let i = k + 1; i < n; i++) {
      const mult = curr[i][k] / pivot;
      for (let j = k; j < curr[0].length; j++) curr[i][j] -= mult * curr[k][j];
    }
    setTrace((t) => [...t, { mat: copyMat(curr), step: `Eliminated column ${k}` }]);
  }

  function backSub() {
    const M = copyMat(trace[trace.length - 1].mat);
    const U = M.map((r) => r.slice(0, -1));
    const b = M.map((r) => r[r.length - 1]);
    const n = U.length;
    const x = Array(n).fill(0);
    for (let i = n - 1; i >= 0; i--) {
      let s = 0;
      for (let j = i + 1; j < n; j++) s += U[i][j] * x[j];
      x[i] = (b[i] - s) / (U[i][i] || 1e-18);
    }
    setTrace((t) => [...t, { mat: M, step: "Back substitution", sol: x }]);
  }

  return (
    <section className="rounded-lg p-4" style={{ background: theme.panelBg }}>
      <motion.div {...fadeIn}>
        <div className="flex items-center gap-3 mb-3">
          <Zap className="w-6 h-6 text-cyan-300" />
          <div>
            <h2 className="text-lg font-semibold text-zinc-100">9.2 Naive Gauss Elimination — Interactive Stepper</h2>
            <div className="text-xs text-zinc-400">Step through elimination; toggle partial pivoting</div>
          </div>
        </div>

        <div className="rounded-lg border border-zinc-700 p-3 bg-zinc-800 mb-3">
          <div className="grid grid-cols-1 md:grid-cols-3 gap-3">
            <div className="md:col-span-2">
              <Labeled label="Augmented matrix [A|b]">
                <textarea rows={4} value={matrixText} onChange={(e) => setMatrixText(e.target.value)} className="w-full bg-zinc-900 p-2 rounded text-zinc-100" />
              </Labeled>
            </div>
            <div className="space-y-2">
              <div className="flex gap-2">
                <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={() => setTrace([{ mat: copyMat(Aaug), step: "Initial matrix" }])}>Reset</Button>
                <Button className="cursor-pointer bg-emerald-300/50 hover:bg-emerald-400 text-black" onClick={nextStep}>Next Step</Button>
                <Button className="cursor-pointer bg-emerald-300/50 hover:bg-emerald-400 text-black" onClick={backSub}>Back-sub</Button>
              </div>
              <div className="flex items-center gap-2 mt-2">
                <input id="pivot" type="checkbox" className="cursor-pointer" checked={usePivot} onChange={(e) => setUsePivot(e.target.checked)} />
                <label htmlFor="pivot" className="text-xs text-zinc-400">Use partial pivoting</label>
              </div>
            </div>
          </div>
        </div>

        <div className="space-y-3">
          {trace.map((s, i) => (
            <div key={i} className="rounded border border-zinc-700 p-3 bg-zinc-900">
              <div className="text-xs text-zinc-400 mb-2">{s.step}</div>
              <div className="overflow-auto">
                <table className="table-auto text-sm">
                  <tbody>
                    {s.mat.map((r, ri) => (
                      <tr key={ri}>
                        {r.map((c, ci) => (
                          <td key={ci} className="px-2 py-1" style={{ minWidth: 52, color: ci === r.length - 1 ? "#93c5fd" : "#e5e7eb" }}>
                            {Number.isFinite(c) ? Number(c).toExponential(4) : c}
                          </td>
                        ))}
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
              {s.sol && <div className="mt-2 text-sm">Solution: <span className="font-mono">[{s.sol.map((v) => v.toExponential(6)).join(", ")}]</span></div>}
            </div>
          ))}
        </div>

        <div className="mt-4">
          <Exercise
            title="Exercise 9.2.1 — Pivot experiment"
            body={`Toggle partial pivoting off and on. Use the sample matrix to observe how pivoting changes the elimination steps and solution stability.`}
            solution={`Partial pivoting swaps rows to place the largest pivot; this avoids large multipliers and increases numerical stability.`}
          />
        </div>
      </motion.div>
    </section>
  );
}

// ----------------------------------------------------------------------
// SECTION 9.3: Pitfalls — Hilbert matrix & conditioning
// ----------------------------------------------------------------------
function Section93() {
  const n = 6;
  const H = useMemo(() => Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => 1 / (i + j + 1))), [n]);
  const b = useMemo(() => Array.from({ length: n }, () => 1), [n]);
  const cond = useMemo(() => condInf(H), [H]);
  const lu = useMemo(() => luDecompose(H), [H]);
  let sol = [];
  if (lu.success) {
    const y = forwardSub(lu.L, matVec(lu.P, b));
    sol = backwardSub(lu.U, y);
  }

  return (
    <section className="rounded-lg p-4" style={{ background: theme.panelBg }}>
      <motion.div {...fadeIn}>
        <div className="flex items-center gap-3 mb-3">
          <Sliders className="w-6 h-6 text-cyan-300" />
          <div>
            <h2 className="text-lg font-semibold text-zinc-100">9.3 Pitfalls of Elimination Methods</h2>
            <div className="text-xs text-zinc-400">Ill-conditioning, rounding, and pivoting cautions</div>
          </div>
        </div>

        <div className="grid grid-cols-1 md:grid-cols-3 gap-3 mb-3">
          <div className="md:col-span-1">
            <div className="text-xs text-zinc-400 mb-1">Hilbert matrix heatmap</div>
            <MatrixHeatmap matrix={H} title="Hilbert (n=6)" />
          </div>
          <div className="md:col-span-2">
            <div className="rounded-lg border border-zinc-700 p-3 bg-zinc-800">
              <div className="text-xs text-zinc-400">Condition estimate (inf)</div>
              <div className="text-sm text-zinc-100 mt-1">{Number.isFinite(cond) ? cond.toFixed(4) : "∞"}</div>
              <div className="text-xs text-zinc-400 mt-2">LU solution (partial pivot)</div>
              <div className="text-sm font-mono mt-1">{lu.success ? `[${sol.map((v) => v.toPrecision(6)).join(", ")}]` : "Decomposition failed"}</div>
              <div className="text-xs text-zinc-400 mt-3">Remarks</div>
              <div className="text-sm text-zinc-300 mt-1">Hilbert matrices are classic examples of extreme ill-conditioning; small rounding errors can cause large changes in the solution.</div>
            </div>
          </div>
        </div>

        <div className="mt-3">
      <Exercise
        title="Exercise 9.3.1 — Perturbation experiment"
        body={`Perturb the RHS by $10^{-8}$ in the first component. Solve again and compute relative change $\\frac{\\lVert \\delta x \\rVert}{\\lVert x \\rVert}$. Compare to estimated condition number.`}
        solution={`If $\\kappa(A)$ is large, the relative solution perturbation will be roughly $\\kappa(A)$ times the relative data perturbation (order-of-magnitude).`}
      />

        </div>
      </motion.div>
    </section>
  );
}

// ----------------------------------------------------------------------
// SECTION 9.4: Techniques — pivoting and iterative refinement
// ----------------------------------------------------------------------
function Section94() {
  const [matrixText, setMatrixText] = useState("4 2 3 1\n2 1 1 0\n2 2 3 1\n1 0 0 1");
  const Aaug = useMemo(() => parseMatrix(matrixText), [matrixText]);
  const A = useMemo(() => Aaug.map((r) => r.slice(0, -1)), [Aaug]);
  const b = useMemo(() => Aaug.map((r) => r[r.length - 1]), [Aaug]);

  const lu = useMemo(() => luDecompose(A), [A]);
  let x = [];
  if (lu.success) {
    const y = forwardSub(lu.L, matVec(lu.P, b));
    x = backwardSub(lu.U, y);
  }

  function refine(xk) {
    const Ax = matVec(A, xk);
    const r = b.map((bi, i) => bi - Ax[i]);
    const y = forwardSub(lu.L, matVec(lu.P, r));
    const dx = backwardSub(lu.U, y);
    const xnew = xk.map((xi, i) => xi + dx[i]);
    return { dx, xnew, r };
  }

  const firstRef = lu.success ? refine(x) : null;

  return (
    <section className="rounded-lg p-4" style={{ background: theme.panelBg }}>
      <motion.div {...fadeIn}>
        <div className="flex items-center gap-3 mb-3">
          <CornerUpLeft className="w-6 h-6 text-cyan-300" />
          <div>
            <h2 className="text-lg font-semibold text-zinc-100">9.4 Techniques for Improving Solutions</h2>
            <div className="text-xs text-zinc-400">Pivoting strategies and iterative refinement</div>
          </div>
        </div>

        <div className="rounded-lg border border-zinc-700 p-3 bg-zinc-800 mb-3">
          <Labeled label="Augmented matrix [A|b]">
            <textarea rows={4} value={matrixText} onChange={(e) => setMatrixText(e.target.value)} className="w-full bg-zinc-900 p-2 rounded text-zinc-100" />
          </Labeled>
          <div className="mt-2">
            <Summary items={[{ label: "LU status", value: lu.success ? "OK" : "Failed" }, { label: "Initial solution", value: lu.success ? `[${x.map((v) => v.toPrecision(6)).join(", ")}]` : "—" }]} />
          </div>
        </div>

        <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
          <div className="rounded border border-zinc-700 p-3 bg-zinc-900">
            <div className="text-xs text-zinc-400">Iterative refinement (one step)</div>
            {firstRef ? (
              <div className="text-sm text-zinc-200 mt-2">
                Residual r = [{firstRef.r.map((v) => v.toExponential(3)).join(", ")}] <br />
                Correction Δx = [{firstRef.dx.map((v) => v.toExponential(3)).join(", ")}] <br />
                Updated x = [{firstRef.xnew.map((v) => v.toPrecision(6)).join(", ")}]
              </div>
            ) : (
              <div className="text-xs text-zinc-400">LU failed — cannot refine</div>
            )}
          </div>
          <div className="rounded border border-zinc-700 p-3 bg-zinc-900">
            <div className="text-xs text-zinc-400">Pivoting strategies</div>
            <div className="text-sm text-zinc-200 mt-2">
              Partial pivoting swaps rows for largest pivot in column. Scaled pivoting divides row entries by row scale. Complete pivoting swaps rows and columns.
            </div>
          </div>
        </div>

        <div className="mt-3">
          <Exercise
            title="Exercise 9.4.1 — Refinement steps"
            body={`Perform two iterative refinement steps on the system and compare residual norms before and after.`}
            solution={`Refinement typically reduces residual for well-conditioned matrices; benefits depend on LU factor quality and machine precision.`}
          />
        </div>
      </motion.div>
    </section>
  );
}

// ----------------------------------------------------------------------
// SECTION 9.5: Complex systems (sparse / banded)
// ----------------------------------------------------------------------
function Section95() {
  const n = 60;
  const A = useMemo(() => {
    const m = Array.from({ length: n }, () => Array.from({ length: n }, () => 0));
    for (let i = 0; i < n; i++) {
      if (i > 0) m[i][i - 1] = -1;
      m[i][i] = 2;
      if (i < n - 1) m[i][i + 1] = -1;
    }
    return m;
  }, [n]);

  const bandwidth = useMemo(() => {
    let bw = 0;
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) if (A[i][j] !== 0) bw = Math.max(bw, Math.abs(i - j));
    return bw;
  }, [A, n]);

  return (
    <section className="rounded-lg p-4" style={{ background: theme.panelBg }}>
      <motion.div {...fadeIn}>
        <div className="flex items-center gap-3 mb-3">
          <Grid className="w-6 h-6 text-cyan-300" />
          <div>
            <h2 className="text-lg font-semibold text-zinc-100">9.5 Complex Systems — Sparse & Banded</h2>
            <div className="text-xs text-zinc-400">Tridiagonal example and bandwidth visualization</div>
          </div>
        </div>

        <div className="rounded border border-zinc-700 p-3 bg-zinc-800 mb-3">
          <MatrixHeatmap matrix={A} title="Tridiagonal (n=60)" />
          <div className="mt-2">
            <Summary items={[{ label: "Bandwidth", value: bandwidth }, { label: "Recommendation", value: "Use Thomas algorithm or sparse solvers" }]} />
          </div>
        </div>

        <div>
          <Exercise
            title="Exercise 9.5.1 — Thomas algorithm"
            body={`Implement Thomas algorithm for tridiagonal systems and compare runtime to dense LU for n=500.`}
            solution={`Thomas algorithm runs in O(n) time and O(n) memory; dense LU is O(n^3) and uses much more memory.`}
          />
        </div>
      </motion.div>
    </section>
  );
}

// ----------------------------------------------------------------------
// SECTION 9.6: Nonlinear systems (Newton's method) + 2D iteration scatter
// ----------------------------------------------------------------------
function Section96() {
  const [x0, setX0] = useState(1);
  const [y0, setY0] = useState(0);

  const F = (v) => {
    const [x, y] = v;
    return [x * x + y * y - 4, Math.exp(x) + y - 1];
  };
  const J = (v) => {
    const [x, y] = v;
    return [
      [2 * x, 2 * y],
      [Math.exp(x), 1],
    ];
  };

  const iterTrace = useMemo(() => {
    const steps = [];
    let xk = [Number(x0), Number(y0)];
    for (let k = 0; k < 12; k++) {
      const Fx = F(xk);
      const Jx = J(xk);
      const lu = luDecompose(Jx);
      let dx = [0, 0];
      if (lu.success) {
        const rhs = Fx.map((v) => -v);
        const y = forwardSub(lu.L, matVec(lu.P, rhs));
        dx = backwardSub(lu.U, y);
      } else break;
      const xnew = [xk[0] + dx[0], xk[1] + dx[1]];
      steps.push({ k, x: [xk[0], xk[1]], Fx, Jx, dx, xnew });
      xk = xnew;
      if (Math.hypot(dx[0], dx[1]) < 1e-12) break;
    }
    return steps;
  }, [x0, y0]);

  const scatter = useMemo(() => {
    const xs = iterTrace.map((s) => s.x[0]);
    const ys = iterTrace.map((s) => s.x[1]);
    const labels = iterTrace.map((s) => `k=${s.k}`);
    return { xs, ys, labels };
  }, [iterTrace]);

  return (
    <section className="rounded-lg p-4" style={{ background: theme.panelBg }}>
      <motion.div {...fadeIn}>
        <div className="flex items-center gap-3 mb-3">
          <BoxSelect className="w-6 h-6 text-cyan-300" />
          <div>
            <h2 className="text-lg font-semibold text-zinc-100">9.6 Nonlinear Systems of Equations</h2>
            <div className="text-xs text-zinc-400">Newton's method for systems — each step solves a linear system</div>
          </div>
        </div>

        <div className="rounded border border-zinc-700 p-3 bg-zinc-800 mb-3">
          <div className="grid grid-cols-1 md:grid-cols-3 gap-3">
            <div>
              <Labeled label="Initial guess x0, y0">
                <div className="flex gap-2">
                  <Input value={x0} onChange={(e) => setX0(e.target.value)} className="w-24 bg-zinc-900 text-zinc-100" />
                  <Input value={y0} onChange={(e) => setY0(e.target.value)} className="w-24 bg-zinc-900 text-zinc-100" />
                </div>
              </Labeled>
            </div>
            <div className="md:col-span-2 text-xs text-zinc-400">The iteration trace below shows Δx computed by solving J(x_k) Δx = -F(x_k) via LU decomposition.</div>
          </div>
        </div>

        <div className="grid grid-cols-1 lg:grid-cols-3 gap-3 mb-3">
          <div className="lg:col-span-2 space-y-2">
            {iterTrace.map((s) => (
              <div key={s.k} className="rounded border border-zinc-700 p-2 bg-zinc-900">
                <div className="text-xs text-zinc-400">k = {s.k}</div>
                <div className="text-sm text-zinc-200">x_k = [{s.x.map((v) => Number(v).toPrecision(6)).join(", ")}]</div>
                <div className="text-sm text-zinc-200">F(x_k) = [{s.Fx.map((v) => Number(v).toPrecision(6)).join(", ")}]</div>
                <div className="text-sm text-zinc-200">Δx = [{s.dx.map((v) => Number(v).toPrecision(6)).join(", ")}]</div>
                <div className="text-sm text-zinc-200">x<sub>{s.k + 1}</sub> = [{s.xnew.map((v) => Number(v).toPrecision(6)).join(", ")}]</div>
              </div>
            ))}
          </div>

          <div className="rounded border border-zinc-700 p-2 bg-zinc-900">
            <div className="text-xs text-zinc-400 mb-2">Iteration Path (2D)</div>
            <Plot
              data={[
                { x: scatter.xs, y: scatter.ys, mode: "markers+lines+text", text: scatter.labels, textposition: "top center", marker: { size: 8 } },
              ]}
              layout={{
                autosize: true,
                height: 300,
                margin: { l: 40, r: 10, t: 30, b: 30 },
                paper_bgcolor: "rgba(0,0,0,0)",
                plot_bgcolor: "rgba(0,0,0,0)",
                xaxis: { title: "x" },
                yaxis: { title: "y" },
              }}
              useResizeHandler
              style={{ width: "100%" }}
            />
          </div>
        </div>

        <div>
          <Exercise
            title="Exercise 9.6.1 — Newton for systems"
            body={`Try initial guesses (1,0), (0,1), (-1,0). For each guess, run iterations and report if they converge and to what root.`}
            solution={`Newton is locally quadratically convergent near a nonsingular root. Poor initial guesses may diverge or converge to a different root.`}
          />
        </div>
      </motion.div>
    </section>
  );
}

// ----------------------------------------------------------------------
// SECTION 9.7: Gauss-Jordan (RREF stepper)
// ----------------------------------------------------------------------
function Section97() {
  const [matrixText, setMatrixText] = useState("2 1 -1 8\n-3 -1 2 -11\n-2 1 2 -3");
  const Aaug = useMemo(() => parseMatrix(matrixText), [matrixText]);
  const [states, setStates] = useState([copyMat(Aaug)]);
  const [index, setIndex] = useState(0);

  useEffect(() => {
    setStates([copyMat(Aaug)]);
    setIndex(0);
  }, [matrixText]);

  function runGaussJordan() {
    const M = copyMat(Aaug);
    const m = M.length;
    const cols = M[0].length;
    let r = 0;
    for (let c = 0; c < cols - 1 && r < m; c++) {
      let piv = r;
      while (piv < m && Math.abs(M[piv][c]) < 1e-14) piv++;
      if (piv === m) continue;
      if (piv !== r) swapRows(M, piv, r);
      const pivot = M[r][c];
      for (let j = c; j < cols; j++) M[r][j] /= pivot;
      setStates((s) => [...s, copyMat(M)]);
      for (let i = 0; i < m; i++) {
        if (i === r) continue;
        const factor = M[i][c];
        for (let j = c; j < cols; j++) M[i][j] -= factor * M[r][j];
      }
      setStates((s) => [...s, copyMat(M)]);
      r++;
    }
    setStates((s) => [...s, copyMat(M)]);
    setIndex((i) => i + 1);
  }

  return (
    <section className="rounded-lg p-4" style={{ background: theme.panelBg }}>
      <motion.div {...fadeIn}>
        <div className="flex items-center gap-3 mb-3">
          <RotateCcw className="w-6 h-6 text-cyan-300" />
          <div>
            <h2 className="text-lg font-semibold text-zinc-100">9.7 Gauss-Jordan — RREF</h2>
            <div className="text-xs text-zinc-400">Step-by-step reduction to reduced row echelon form</div>
          </div>
        </div>

        <div className="rounded-lg border border-zinc-700 p-3 bg-zinc-800 mb-3">
          <Labeled label="Augmented matrix [A|b]">
            <textarea rows={4} value={matrixText} onChange={(e) => setMatrixText(e.target.value)} className="w-full bg-zinc-900 p-2 rounded text-zinc-100" />
          </Labeled>

          <div className="flex gap-2 mt-2">
            <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={() => { setStates([copyMat(Aaug)]); setIndex(0); }}>Reset</Button>
            <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={runGaussJordan}>Run Gauss-Jordan</Button>
            <Button className="cursor-pointer bg-emerald-300/50 hover:bg-emerald-400 text-black" onClick={() => setIndex((i) => Math.max(0, i - 1))}>Prev</Button>
            <Button className="cursor-pointer bg-emerald-300/50 hover:bg-emerald-400 text-black" onClick={() => setIndex((i) => Math.min(states.length - 1, i + 1))}>Next</Button>
          </div>

          <div className="text-xs text-zinc-400 mt-2">States: {states.length} • Current: {index}</div>
        </div>

        <div className="rounded border border-zinc-700 p-3 bg-zinc-900">
          <div className="text-xs text-zinc-400 mb-2">RREF state</div>
          <div className="overflow-auto">
            <table className="table-auto text-sm">
              <tbody>
                {states[index]?.map((r, ri) => (
                  <tr key={ri}>
                    {r.map((c, ci) => (
                      <td key={ci} className="px-2 py-1" style={{ minWidth: 56, color: ci === r.length - 1 ? "#93c5fd" : "#e5e7eb" }}>
                        {Number.isFinite(c) ? Number(c).toFixed(6) : c}
                      </td>
                    ))}
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
          <div className="text-xs text-zinc-400 mt-2">Gauss-Jordan shows free variables and inconsistency directly in RREF.</div>
        </div>

        <div className="mt-3">
          <Exercise
            title="Exercise 9.7.1 — RREF for consistency"
            body={`Use Gauss-Jordan to determine if a given augmented system is consistent. Identify free variables if any.`}
            solution={`If RREF yields a row like [0 ... 0 | c] with c ≠ 0, system is inconsistent. Otherwise, free variables correspond to non-pivot columns.`}
          />
        </div>
      </motion.div>
    </section>
  );
}

// ----------------------------------------------------------------------
// SECTION 9.8: Summary & problems
// ----------------------------------------------------------------------
function Section98() {
  const problems = [
    {
      title: "Problem 9.1 — Solve the small system",
      body: `Solve the system:
$$
\\begin{cases}
2x + y - z = 8 \\\\
-3x - y + 2z = -11 \\\\
-2x + y + 2z = -3
\\end{cases}
$$ 
Show elimination steps and final solution.`,
      solution: `Follow elimination: final solution is $x = 2$, $y = 3$, $z = 1$.`,
    },
    {
      title: "Problem 9.2 — Pivoting necessity",
      body: `Use:
$$
\\begin{bmatrix}
10^{-8} & 1 & 1 \\\\
1 & 1 & 1 \\\\
2 & 1 & -1
\\end{bmatrix}
$$
Show naive elimination fails and partial pivoting fixes it.`,
      solution: `Partial pivoting places a larger pivot to avoid division by a tiny number, stabilizing elimination.`,
    },
    {
      title: "Problem 9.3 — LU reuse",
      body: `Factor $A$ once via LU (with pivoting) and solve for 3 different RHS vectors. 
Explain the computational savings versus recomputing elimination.`,
      solution: `LU factorization is $O(n^3)$; each solve is $O(n^2)$. 
Reusing LU greatly reduces cost for multiple RHS problems.`,
    },
    {
      title: "Problem 9.4 — Conditioning",
      body: `Compute condition estimate for $\\text{Hilbert}(8)$ and show that perturbing 
the RHS by $10^{-8}$ causes large relative solution error.`,
      solution: `Hilbert matrices have very large $\\kappa(A)$. 
Even tiny RHS perturbations cause large relative solution changes.`,
    },
  ];


  return (
    <section className="rounded-lg p-4" style={{ background: theme.panelBg }}>
      <motion.div {...fadeIn}>
        <div className="flex items-center gap-3 mb-3">
          <Code className="w-6 h-6 text-cyan-300" />
          <div>
            <h2 className="text-lg font-semibold text-zinc-100">9.8 Summary & Problems</h2>
            <div className="text-xs text-zinc-400">Exercises and end-of-chapter problems</div>
          </div>
        </div>

        <div className="grid grid-cols-1 gap-3">
          {problems.map((p, i) => (
            <Exercise key={i} title={p.title} body={p.body} solution={p.solution} />
          ))}
        </div>
      </motion.div>
    </section>
  );
}

// ----------------------------------------------------------------------
// Page assembly: all sections stacked
// ----------------------------------------------------------------------
export default function Chapter9Page() {
  return (
    <div className={`min-h-screen p-6 bg-zinc-950 text-zinc-100`}>
      <header className="mb-6">
        <div className="flex items-center justify-between gap-4">
          <div>
            <h1 className="text-3xl font-bold" style={{ color: theme.accent2 }}>Gauss Elimination</h1>
            <div className="text-sm text-zinc-400 mt-1">Interactive, visualization-first exploration of direct linear solvers</div>
          </div>
          
        </div>
      </header>

      <main className="space-y-6">
        <Section91 />
        <Section92 />
        <Section93 />
        <Section94 />
        <Section95 />
        <Section96 />
        <Section97 />
        <Section98 />
        <BottomBar/>
      </main>

      <footer className="mt-8 text-xs text-zinc-500">Note: This educational demo is for visualization and teaching. For production, use optimized numerical libraries (LAPACK, Eigen, SuiteSparse).</footer>
    </div>
  );
}
