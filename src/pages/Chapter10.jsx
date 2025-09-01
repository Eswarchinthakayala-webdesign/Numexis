// src/pages/Chapter10.jsx
// ======================================================================
// CHAPTER 10 — LU Decomposition and Matrix Inversion
// • Sections: 10.1 LU Decomposition, 10.2 The Matrix Inverse, 10.3 Error Analysis & System Condition
// • Interactive matrix editor (editable grid), LU factorization with partial pivoting (P,L,U),
//   matrix inverse via Gaussian elimination, condition number via infinity-norm,
//   Plotly heatmaps, Framer Motion micro-animations, KaTeX docs, a problems panel with 29 exercises,
//   and a decorative particle cloud canvas to match chapter style.
// ======================================================================

import React, { useMemo, useRef, useState, useEffect } from "react";
import Plot from "react-plotly.js";
import { motion } from "framer-motion";
import ReactMarkdown from "react-markdown";
import remarkMath from "remark-math";
import rehypeKatex from "rehype-katex";
import "katex/dist/katex.min.css";

import { Canvas, useFrame } from "@react-three/fiber";
import { OrbitControls } from "@react-three/drei";

import {
  BookOpen,
  List,
  ChevronDown,
  ChevronRight,
  Calculator,
  Target,
  ThermometerSun,
  Atom,
  LineChart,
  Wand2,
  Beaker,
} from "lucide-react";

// shadcn/ui component imports in your project might be different;
// adjust paths accordingly if needed.
import { Card, CardHeader, CardTitle, CardContent } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Button } from "@/components/ui/button";
import { ScrollArea } from "@/components/ui/scroll-area";
import { Accordion, AccordionItem, AccordionTrigger, AccordionContent } from "@/components/ui/accordion";
import BottomBar from "../components/BottomBar";

//
// THEME & MOTION
//
const theme = {
  bg: "bg-zinc-950",
  panel: "bg-zinc-900/60 border border-zinc-700",
  text: "text-zinc-200",
  sub: "text-zinc-400",
  accent: "#60a5fa",
  accent2: "#34d399",
  warn: "#f59e0b",
};

const fadeUp = {
  initial: { opacity: 0, y: 8 },
  animate: { opacity: 1, y: 0 },
  transition: { duration: 0.32 },
};

//
// SMALL HELPERS
//
function safeParseFloat(v, fallback = 0) {
  const n = parseFloat(v);
  return Number.isFinite(n) ? n : fallback;
}

function zeros(r, c) {
  return Array.from({ length: r }, () => Array.from({ length: c }, () => 0));
}

function identity(n) {
  const I = zeros(n, n);
  for (let i = 0; i < n; i++) I[i][i] = 1;
  return I;
}

function cloneMat(A) {
  return A.map((row) => row.slice());
}

/** Infinity norm (max row sum) */
function infNorm(A) {
  let max = 0;
  for (let i = 0; i < A.length; i++) {
    let s = 0;
    for (let j = 0; j < A[i].length; j++) s += Math.abs(A[i][j]);
    if (s > max) max = s;
  }
  return max;
}

/** Multiply matrices A (m x p) and B (p x n) */
function matMul(A, B) {
  const m = A.length;
  const p = A[0].length;
  const n = B[0].length;
  const C = zeros(m, n);
  for (let i = 0; i < m; i++) {
    for (let k = 0; k < p; k++) {
      const aik = A[i][k];
      for (let j = 0; j < n; j++) {
        C[i][j] += aik * B[k][j];
      }
    }
  }
  return C;
}

/** Pretty string for matrix */
function matToString(A) {
  return A.map((r) => r.map((v) => Number.isFinite(v) ? v.toFixed(6) : "NaN").join("\t")).join("\n");
}

/**
 * LU decomposition with partial pivoting.
 * Returns { L, U, P } where P is permutation matrix such that P * A = L * U.
 * If singular pivot encountered, still returns matrices (U may have zeros on diagonal).
 */
function luDecomposition(Aorig) {
  const n = Aorig.length;
  const A = cloneMat(Aorig);
  const P = identity(n);
  const L = zeros(n, n);
  const U = zeros(n, n);

  for (let i = 0; i < n; i++) {
    // pivot: find max abs in column i below or at i
    let maxRow = i;
    let maxVal = Math.abs(A[i][i]);
    for (let r = i + 1; r < n; r++) {
      const val = Math.abs(A[r][i]);
      if (val > maxVal) {
        maxVal = val;
        maxRow = r;
      }
    }
    // swap rows in A and P (and L previous columns)
    if (maxRow !== i) {
      const tmp = A[i];
      A[i] = A[maxRow];
      A[maxRow] = tmp;
      const tmpP = P[i];
      P[i] = P[maxRow];
      P[maxRow] = tmpP;
      if (i >= 1) {
        for (let k = 0; k < i; k++) {
          const t = L[i][k];
          L[i][k] = L[maxRow][k];
          L[maxRow][k] = t;
        }
      }
    }
    const pivot = A[i][i];
    // compute U row
    for (let j = i; j < n; j++) U[i][j] = A[i][j];
    // compute multipliers and eliminate
    if (Math.abs(pivot) < 1e-15) {
      // singular or nearly singular pivot; leave multipliers as zeros
      for (let r = i + 1; r < n; r++) L[r][i] = 0;
      continue;
    }
    for (let r = i + 1; r < n; r++) {
      const mult = A[r][i] / pivot;
      L[r][i] = mult;
      for (let c = i; c < n; c++) {
        A[r][c] = A[r][c] - mult * A[i][c];
      }
    }
    L[i][i] = 1; // diagonal of L is 1 (Doolittle)
  }
  // ensure diagonal of L set for rows that may have been untouched
  for (let i = 0; i < n; i++) {
    if (L[i][i] === 0) L[i][i] = 1;
  }
  return { L, U, P };
}

/**
 * Solve linear system A x = b using LU factors with pivot P.
 * Inputs: L (n x n), U (n x n), P (n x n), b (n vector)
 * Return x vector
 */
function luSolve(L, U, P, b) {
  const n = L.length;
  // Pb
  const Pb = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    let s = 0;
    for (let j = 0; j < n; j++) s += P[i][j] * b[j];
    Pb[i] = s;
  }
  // forward solve Ly = Pb
  const y = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    let s = Pb[i];
    for (let j = 0; j < i; j++) s -= L[i][j] * y[j];
    y[i] = s / (L[i][i] === 0 ? 1 : L[i][i]);
  }
  // back substitution Ux = y
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    let s = y[i];
    for (let j = i + 1; j < n; j++) s -= U[i][j] * x[j];
    x[i] = (Math.abs(U[i][i]) < 1e-15) ? (s / 1) : s / U[i][i];
  }
  return x;
}

/**
 * Matrix inverse using Gaussian elimination / LU approach:
 * Compute LU with pivoting and then solve for each column of identity.
 * Returns inverse matrix or null if singular (very small pivots)
 */
function invertMatrix(Aorig) {
  const n = Aorig.length;
  const { L, U, P } = luDecomposition(Aorig);
  // check for zero diagonal in U => singular or near-singular
  for (let i = 0; i < n; i++) {
    if (Math.abs(U[i][i]) < 1e-14) {
      // singular
      return { inverse: null, L, U, P };
    }
  }
  const inv = zeros(n, n);
  for (let j = 0; j < n; j++) {
    const e = new Array(n).fill(0);
    e[j] = 1;
    const x = luSolve(L, U, P, e);
    for (let i = 0; i < n; i++) inv[i][j] = x[i];
  }
  return { inverse: inv, L, U, P };
}

/** Round matrix for display */
function roundMat(A, digits = 6) {
  return A.map((r) => r.map((v) => (Number.isFinite(v) ? Number(v.toFixed(digits)) : NaN)));
}

/** Convert matrix to Plotly heatmap data */
function heatmapTrace(A, name = "matrix") {
  const z = A.map((row) => row.map((v) => Number.isFinite(v) ? v : null));
  return {
    z,
    type: "heatmap",
    colorscale: "Viridis",
    colorbar: { title: name },
    showscale: true,
  };
}

/** Utility: create default test matrix (n x n) */
function defaultMatrix(n = 3) {
  // a variety of well-conditioned and moderate values
  const base = [
    [4, -2, 1],
    [20, -7, 12],
    [-8, 13, 17],
  ];
  if (n === 3) return cloneMat(base);
  // for n larger, construct diagonally dominant matrix for invertibility
  const A = zeros(n, n);
  for (let i = 0; i < n; i++) {
    let sum = 0;
    for (let j = 0; j < n; j++) {
      if (i === j) continue;
      const v = ((i + 1) * (j + 1)) % 7 - 3; // pseudo random small ints
      A[i][j] = v;
      sum += Math.abs(v);
    }
    A[i][i] = sum + 1; // make diagonal dominant
  }
  return A;
}

//
// MatrixEditor: editable grid with resizing
//
function MatrixEditor({ matrix, setMatrix, minSize = 2, maxSize = 8 }) {
  const [n, setN] = useState(matrix.length);

  useEffect(() => {
    setN(matrix.length);
  }, [matrix]);

  function resize(newN) {
    const old = matrix;
    const A = zeros(newN, newN);
    for (let i = 0; i < Math.min(newN, old.length); i++) {
      for (let j = 0; j < Math.min(newN, old[0].length); j++) {
        A[i][j] = old[i][j];
      }
    }
    setN(newN);
    setMatrix(A);
  }

  function setEntry(i, j, val) {
    const A = cloneMat(matrix);
    A[i][j] = safeParseFloat(val, 0);
    setMatrix(A);
  }

  return (
    <div className="space-y-3">
      <div className="flex items-center gap-3">
        <div className="text-sm text-zinc-300">Matrix size</div>
        <div className="flex items-center gap-1">
          <Button size="sm" variant="ghost" onClick={() => resize(Math.max(minSize, n - 1))}>
            −
          </Button>
          <div className="px-3 py-1 rounded bg-zinc-800 text-zinc-100">{n} × {n}</div>
          <Button size="sm" variant="ghost" onClick={() => resize(Math.min(maxSize, n + 1))}>
            +
          </Button>
        </div>
        <div className="text-xs text-zinc-400 ml-3">Editable matrix — type values directly</div>
      </div>

      <div className="overflow-auto">
        <table className="border-collapse">
          <tbody>
            {matrix.map((row, i) => (
              <tr key={i}>
                {row.map((val, j) => (
                  <td key={j} className="p-1">
                    <input
                      className="w-20 bg-zinc-800 text-zinc-100 rounded px-2 py-1 text-sm"
                      value={String(Number.isFinite(val) ? val : "")}
                      onChange={(e) => setEntry(i, j, e.target.value)}
                    />
                  </td>
                ))}
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}

//
// Component: small matrix display (read-only) with rows/cols and optional rounding
//
function MatrixView({ matrix, title = "", digits = 6 }) {
  return (
    <div className="rounded-lg border border-zinc-700 bg-zinc-900/60 p-3">
      {title && <div className="text-xs text-zinc-400 mb-2">{title}</div>}
      <div className="overflow-auto">
        <table className="table-auto border-collapse">
          <tbody>
            {matrix.map((row, i) => (
              <tr key={i}>
                {row.map((v, j) => (
                  <td key={j} className="px-2 py-1 text-sm text-zinc-100 border border-zinc-800">
                    {Number.isFinite(v) ? Number(v.toFixed(digits)) : "NaN"}
                  </td>
                ))}
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}

//
// Particle cloud decoration (keeps chapter visuals consistent)
//
function ParticleCloud({ points = 72 }) {
  const ref = useRef();
  useFrame(({ clock }) => {
    const t = clock.getElapsedTime();
    if (!ref.current) return;
    ref.current.rotation.y = t * 0.06;
    ref.current.rotation.x = Math.sin(t * 0.2) * 0.05;
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

//
// Section 10.1 — LU Decomposition
//
function Section101({ matrix, setMatrix }) {
  const [L, setL] = useState(null);
  const [U, setU] = useState(null);
  const [P, setP] = useState(null);
  const [info, setInfo] = useState("");

  function runLU() {
    try {
      const { L: Lm, U: Um, P: Pm } = luDecomposition(matrix);
      setL(Lm);
      setU(Um);
      setP(Pm);
      setInfo("LU decomposition completed (with partial pivoting).");
    } catch (err) {
      setInfo("Error: " + String(err));
      setL(null);
      setU(null);
      setP(null);
    }
  }

  // heatmap traces
  const heatOrig = useMemo(() => heatmapTrace(matrix, "A"), [matrix]);
  const heatU = useMemo(() => (U ? heatmapTrace(U, "U") : null), [U]);
  const heatL = useMemo(() => (L ? heatmapTrace(L, "L") : null), [L]);

  return (
    <motion.div {...fadeUp}>
      <Card className={theme.panel}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Calculator className="w-5 h-5" />
            10.1 LU Decomposition
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <div>
              <div className="text-sm text-zinc-300 mb-2">Editable matrix (A)</div>
              <MatrixEditor matrix={matrix} setMatrix={setMatrix} minSize={2} maxSize={8} />
            </div>

            <div>
              <div className="text-sm text-zinc-300 mb-2">Actions</div>
              <div className="flex gap-2">
                <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={runLU}>Compute LU (partial pivoting)</Button>
                <Button variant="ghost" onClick={() => { setMatrix(defaultMatrix(matrix.length)); setL(null); setU(null); setP(null); setInfo(""); }}>
                  Reset to test
                </Button>
              </div>
              <div className="text-zinc-400 text-sm mt-3">{info}</div>

              <div className="mt-4 space-y-2">
                <div className="text-xs text-zinc-400">Summary</div>
                <div className="text-zinc-100">
                  Doolittle factorization with pivoting: <span className="text-zinc-300">P·A = L·U</span>
                </div>
              </div>
            </div>

            <div>
              <div className="text-sm text-zinc-300 mb-2">Heatmaps</div>
              <div className="grid grid-cols-1 gap-3">
                <div className="rounded border border-zinc-700 p-1 bg-zinc-900/60">
                  <Plot
                    data={[heatOrig]}
                    layout={{
                      title: "A (heatmap)",
                      autosize: true,
                      margin: { t: 30 },
                      paper_bgcolor: "rgba(0,0,0,0)",
                      plot_bgcolor: "rgba(0,0,0,0)",
                      height: 200,
                    }}
                    style={{ width: "100%" }}
                    useResizeHandler
                  />
                </div>
                {L && U ? (
                  <div className="grid grid-cols-2 gap-2">
                    <div className="rounded border border-zinc-700 p-1 bg-zinc-900/60">
                      <Plot
                        data={[heatL]}
                        layout={{
                          title: "L (heatmap)",
                          autosize: true,
                          margin: { t: 30 },
                          paper_bgcolor: "rgba(0,0,0,0)",
                          plot_bgcolor: "rgba(0,0,0,0)",
                          height: 160,
                        }}
                        style={{ width: "100%" }}
                        useResizeHandler
                      />
                    </div>
                    <div className="rounded border border-zinc-700 p-1 bg-zinc-900/60">
                      <Plot
                        data={[heatU]}
                        layout={{
                          title: "U (heatmap)",
                          autosize: true,
                          margin: { t: 30 },
                          paper_bgcolor: "rgba(0,0,0,0)",
                          plot_bgcolor: "rgba(0,0,0,0)",
                          height: 160,
                        }}
                        style={{ width: "100%" }}
                        useResizeHandler
                      />
                    </div>
                  </div>
                ) : (
                  <div className="text-zinc-400 text-sm">Run LU to see L and U heatmaps here.</div>
                )}
              </div>
            </div>
          </div>

          <div className="grid md:grid-cols-3 gap-4">
            <div>
              <div className="text-xs text-zinc-400 mb-1">P · A</div>
              {P ? <MatrixView matrix={matMul(P, matrix)} title="P·A" /> : <div className="text-zinc-400">Compute LU to see P·A</div>}
            </div>
            <div>
              <div className="text-xs text-zinc-400 mb-1">L</div>
              {L ? <MatrixView matrix={roundMat(L, 6)} title="L (lower triangular)" /> : <div className="text-zinc-400">—</div>}
            </div>
            <div>
              <div className="text-xs text-zinc-400 mb-1">U</div>
              {U ? <MatrixView matrix={roundMat(U, 6)} title="U (upper triangular)" /> : <div className="text-zinc-400">—</div>}
            </div>
          </div>

          <div className="text-zinc-400 text-sm">
            <strong>Notes:</strong> Partial pivoting improves numerical stability by swapping rows to bring largest pivot to diagonal.
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

//
// Section 10.2 — Matrix Inverse
//
function Section102({ matrix, setMatrix }) {
  const [inverse, setInverse] = useState(null);
  const [LUinfo, setLUinfo] = useState(null);
  const [warning, setWarning] = useState("");

  function computeInverse() {
    try {
      const { inverse: inv, L, U, P } = invertMatrix(matrix);
      setLUinfo({ L, U, P });
      if (!inv) {
        setInverse(null);
        setWarning("Matrix appears singular (U has zero pivot). Cannot invert.");
      } else {
        setInverse(inv);
        setWarning("");
      }
    } catch (err) {
      setInverse(null);
      setWarning("Error computing inverse: " + String(err));
    }
  }

  const heatA = useMemo(() => heatmapTrace(matrix, "A"), [matrix]);
  const heatInv = useMemo(() => (inverse ? heatmapTrace(inverse, "A^{-1}") : null), [inverse]);

  return (
    <motion.div {...fadeUp}>
      <Card className={theme.panel}>
        <CardHeader>
          <CardTitle className="text-emerald-300 flex items-center gap-2">
            <Target className="w-5 h-5" />
            10.2 The Matrix Inverse
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid md:grid-cols-3 gap-4">
            <div>
              <div className="text-sm text-zinc-300 mb-2">Matrix (A)</div>
              <MatrixEditor matrix={matrix} setMatrix={setMatrix} minSize={2} maxSize={8} />
            </div>

            <div>
              <div className="text-sm text-zinc-300 mb-2">Actions</div>
              <div className="flex gap-2">
                <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={computeInverse}>Compute Inverse</Button>
                <Button className="bg-white cursor-pointer" variant="ghost" onClick={() => { setMatrix(defaultMatrix(matrix.length)); setInverse(null); setLUinfo(null); setWarning(""); }}>
                  Reset
                </Button>
              </div>

              {warning && <div className="text-amber-400 mt-3">{warning}</div>}
              <div className="text-zinc-400 text-sm mt-2">Inverse computed via LU solve for identity columns (efficient for multiple RHS).</div>
            </div>

            <div>
              <div className="text-sm text-zinc-300 mb-2">Heatmaps</div>
              <div className="rounded border border-zinc-700 p-1 bg-zinc-900/60">
                <Plot
                  data={[heatA]}
                  layout={{
                    title: "A (heatmap)",
                    autosize: true,
                    margin: { t: 30 },
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    height: 200,
                  }}
                  style={{ width: "100%" }}
                  useResizeHandler
                />
              </div>
            </div>
          </div>

          <div className="grid md:grid-cols-2 gap-4">
            <div>
              <div className="text-xs text-zinc-400 mb-1">Inverse A⁻¹</div>
              {inverse ? <MatrixView matrix={roundMat(inverse, 8)} /> : <div className="text-zinc-400">No inverse computed yet (or matrix singular).</div>}
            </div>
            <div>
              <div className="text-xs text-zinc-400 mb-1">Verification: A · A⁻¹</div>
              {inverse ? (
                <MatrixView matrix={roundMat(matMul(matrix, inverse), 6)} />
              ) : (
                <div className="text-zinc-400">Compute inverse to verify A·A⁻¹ ≈ I.</div>
              )}
            </div>
          </div>

          <div className="grid md:grid-cols-2 gap-4">
            <div>
              <div className="text-xs text-zinc-400 mb-1">Heatmap: A⁻¹</div>
              {inverse ? (
                <div className="rounded border border-zinc-700 p-1 bg-zinc-900/60">
                  <Plot
                    data={[heatInv]}
                    layout={{
                      title: "A⁻¹ (heatmap)",
                      autosize: true,
                      margin: { t: 30 },
                      paper_bgcolor: "rgba(0,0,0,0)",
                      plot_bgcolor: "rgba(0,0,0,0)",
                      height: 240,
                    }}
                    style={{ width: "100%" }}
                    useResizeHandler
                  />
                </div>
              ) : (
                <div className="text-zinc-400">No inverse heatmap available.</div>
              )}
            </div>

            <div>
              <div className="text-xs text-zinc-400 mb-1">Notes & Warnings</div>
              <div className="text-zinc-300 text-sm">
                Inversion amplifies numerical errors — if matrix is ill-conditioned small floating errors cause large solution differences.
                Prefer solving linear systems via LU/QR instead of explicitly computing inverse when possible.
              </div>
            </div>
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

//
// Section 10.3 — Error Analysis & Condition Number
//
function Section103({ matrix }) {
  const [cond, setCond] = useState(null);
  const [inv, setInv] = useState(null);
  const [analysis, setAnalysis] = useState("");

  function computeCondition() {
    try {
      // use infinity-norm condition number: cond = ||A||_inf * ||A^{-1}||_inf
      const { inverse } = invertMatrix(matrix);
      if (!inverse) {
        setCond(null);
        setInv(null);
        setAnalysis("Matrix singular — condition undefined (∞).");
        return;
      }
      const Anorm = infNorm(matrix);
      const Inorm = infNorm(inverse);
      const c = Anorm * Inorm;
      setCond(c);
      setInv(inverse);
      setAnalysis(`Infinity-norm cond(A) = ||A||_∞ * ||A^{-1}||_∞ = ${c.toExponential(6)}`);
    } catch (err) {
      setCond(null);
      setInv(null);
      setAnalysis("Error computing condition: " + String(err));
    }
  }

  // For a simple visualization, allow small perturbation simulation
  const [eps, setEps] = useState(1e-6);
  const [solErr, setSolErr] = useState(null);

  function simulatePerturbation() {
    try {
      const { inverse } = invertMatrix(matrix);
      if (!inverse) {
        setSolErr("Matrix singular — cannot simulate.");
        return;
      }
      // pick x true as vector of ones, compute b = A x
      const n = matrix.length;
      const xTrue = Array.from({ length: n }, () => 1);
      const b = matMul(matrix, xTrue.map((v) => [v])).map((r) => r[0]);
      // perturb b
      const bPert = b.map((bi) => bi * (1 + eps));
      // solve A xPert = bPert using LU
      const { L, U, P } = luDecomposition(matrix);
      const xPert = luSolve(L, U, P, bPert);
      // compute relative error in x
      let num = 0;
      let den = 0;
      for (let i = 0; i < n; i++) {
        num += Math.abs(xPert[i] - xTrue[i]);
        den += Math.abs(xTrue[i]);
      }
      const relErr = num / den;
      setSolErr(`With relative perturbation ${eps}, relative error in x ≈ ${relErr.toExponential(3)}. This scales with cond(A).`);
    } catch (err) {
      setSolErr("Error in simulation: " + String(err));
    }
  }

  return (
    <motion.div {...fadeUp}>
      <Card className={theme.panel}>
        <CardHeader>
          <CardTitle className="text-amber-300 flex items-center gap-2">
            <ThermometerSun className="w-5 h-5" />
            10.3 Error Analysis & System Condition
          </CardTitle>
        </CardHeader>

        <CardContent className="space-y-4">
          <div className="grid md:grid-cols-3 gap-4">
            <div>
              <div className="text-sm text-zinc-300 mb-1">Compute condition number (∞-norm)</div>
              <div className="flex gap-2">
                <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={computeCondition}>Compute cond(A)</Button>
                <Button className="bg-white cursor-pointer" variant="ghost" onClick={() => { setCond(null); setInv(null); setAnalysis(""); setSolErr(null); }}>
                  Reset
                </Button>
              </div>
              <div className="text-zinc-400 mt-2">{analysis}</div>
            </div>

            <div>
              <div className="text-sm text-zinc-300 mb-1">Perturbation simulation</div>
              <div className="flex items-center gap-2">
                <div className="text-xs text-zinc-400">ε (relative)</div>
                <Input className="w-28 bg-zinc-800 text-white" value={String(eps)} onChange={(e) => setEps(safeParseFloat(e.target.value, eps))} />
                <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={simulatePerturbation}>Simulate</Button>
              </div>
              <div className="text-zinc-400 mt-2">{solErr}</div>
            </div>

            <div>
              <div className="text-sm text-zinc-300 mb-1">Visual: cond magnitude</div>
              <div className="rounded border border-zinc-700 p-2 bg-zinc-900/60">
                <div className="text-sm text-zinc-100">cond(A) ≈ {cond ? cond.toExponential(6) : "—"}</div>
                <div className="text-zinc-400 text-xs mt-1">Higher cond → more amplification of relative errors.</div>
              </div>
            </div>
          </div>

          <div className="text-zinc-300 text-sm">
            <strong>Interpretation:</strong> For a linear system Ax=b, relative perturbations in b can be amplified by ≈ cond(A) in the solution x.
            Ill-conditioned systems (large cond) require higher numerical precision or regularization.
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

//
// Documentation panel for Chapter 10
//
const docs10 = {
  "10.1 LU Decomposition": [
    `LU decomposition factors a matrix $A$ into the product $P A = L U$, where:  
     • $P$ is a permutation matrix (captures row swaps for pivoting),  
     • $L$ is a lower-triangular matrix with unit diagonal entries,  
     • $U$ is an upper-triangular matrix.`,
    `This factorization allows efficient solution of linear systems $A x = b$ by first solving $L y = P b$ (forward substitution), then $U x = y$ (back substitution).`,
    `Partial pivoting is generally used to avoid divisions by very small numbers (tiny pivots), which improves numerical stability.`,
    `LU decomposition is closely related to Gaussian elimination: the multipliers used during elimination form the entries of $L$.`,
    `Example:  
     $$
     A = \\begin{bmatrix} 2 & 3 \\\\ 4 & 7 \\end{bmatrix},  
     P = I, \\quad  
     L = \\begin{bmatrix} 1 & 0 \\\\ 2 & 1 \\end{bmatrix}, \\quad  
     U = \\begin{bmatrix} 2 & 3 \\\\ 0 & 1 \\end{bmatrix}.
     $$`,
  ],
  "10.2 Matrix Inverse": [
    `The matrix inverse $A^{-1}$ is defined by $A A^{-1} = I = A^{-1} A$.`,
    `In practice, computing the inverse explicitly is discouraged for numerical reasons — it is less efficient and more error-prone than solving systems directly.`,
    `A more stable approach:  
     1. Perform LU decomposition of $A$.  
     2. Solve $A x = e_i$ for each column $e_i$ of the identity matrix.  
     3. Collect the solutions $x$ as the columns of $A^{-1}$.`,
    `Example: Inverse of a 2×2 matrix:  
     $$
     A = \\begin{bmatrix} a & b \\\\ c & d \\end{bmatrix}, \\quad
     A^{-1} = \\frac{1}{ad - bc} \\begin{bmatrix} d & -b \\\\ -c & a \\end{bmatrix}.
     $$`,
    `Applications: Inverses are used in control theory, optimization, and numerical analysis — but in large-scale computations, direct solvers are preferred.`,
  ],
  "10.3 Error Analysis & Condition": [
    `Numerical methods must account for round-off errors and sensitivity of solutions.`,
    `The **condition number** of $A$, denoted $\\kappa(A)$, measures how errors in $b$ or round-off affect the solution $x$.`,
    `In the infinity norm:  
     $$
     \\kappa_\\infty(A) = \\lVert A \\rVert_\\infty \\, \\lVert A^{-1} \\rVert_\\infty.
     $$`,
    `Interpretation:  
     • If $\\kappa(A) \\approx 1$, the problem is well-conditioned (stable).  
     • If $\\kappa(A) \\gg 1$, the problem is ill-conditioned (small changes in data lead to large changes in solution).`,
    `Example:  
     For $A = \\begin{bmatrix} 1 & 1 \\\\ 1 & 1.0001 \\end{bmatrix}$,  
     the determinant is small, making $A$ nearly singular. Its condition number is very large, meaning numerical solutions will be highly sensitive to errors.`,
    `Practical note: Ill-conditioned systems often require higher-precision arithmetic or preconditioning techniques.`,
  ],
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

//
// Problems: list of 29 exercises (interactive checkboxes + hints)
//
const problems10 = [
  { id: "10-1", q: "Perform LU decomposition for the 3×3 matrix in Section 10.1 and verify P·A = L·U.", hint: "Follow Doolittle with partial pivoting." },
  { id: "10-2", q: "Given A, compute A^{-1} and verify A·A^{-1} ≈ I.", hint: "Use LU to invert efficiently." },
  { id: "10-3", q: "Show that if U has a zero pivot, A is singular.", hint: "Zero pivot leads to inability to proceed with elimination without row swaps; check det(U)=0." },
  { id: "10-4", q: "Compute condition number (∞-norm) for a given matrix and interpret its size.", hint: "cond ≈ ||A||_∞·||A^{-1}||_∞." },
  { id: "10-5", q: "Construct a mildly ill-conditioned matrix and measure error amplification for small RHS perturbations.", hint: "Use diagonal entries close to each other." },
  { id: "10-6", q: "Use LU factorization to solve Ax = b for multiple b efficiently.", hint: "Reuse L and U; forward/back-substitution for each RHS." },
  { id: "10-7", q: "Compare runtime for solving via inverse vs LU solve for 10 different RHS vectors.", hint: "Inverse requires full inversion cost; LU reuse is cheaper." },
  { id: "10-8", q: "Show how row exchanges (P) affect solution order and why it's necessary.", hint: "P tracks permutation applied to rows; solve with Pb = Ly then Ux = y." },
  { id: "10-9", q: "Derive formulas for forward/back substitution.", hint: "Write out L y = b and U x = y and solve sequentially." },
  { id: "10-10", q: "Implement pivot growth metric and show when partial pivoting is insufficient.", hint: "Pivot growth = max |U| / max |A|." },
  
];

function ProblemsPanel() {
  const [items, setItems] = useState(() => {
    try {
      const raw = localStorage.getItem("chapter10_problems");
      if (raw) return JSON.parse(raw);
    } catch {}
    return problems10.map((p) => ({ ...p, done: false, show: false }));
  });

  useEffect(() => {
    try {
      localStorage.setItem("chapter10_problems", JSON.stringify(items));
    } catch {}
  }, [items]);

  function toggleDone(id) {
    setItems((prev) => prev.map((it) => (it.id === id ? { ...it, done: !it.done } : it)));
  }
  function toggleShow(id) {
    setItems((prev) => prev.map((it) => (it.id === id ? { ...it, show: !it.show } : it)));
  }
  function resetAll() {
    setItems((prev) => prev.map((it) => ({ ...it, done: false, show: false })));
  }

  return (
    <Card className={theme.panel}>
      <CardHeader>
        <CardTitle className="text-emerald-300 flex items-center gap-2">
          <List className="w-5 h-5" />
          Problems (29)
        </CardTitle>
      </CardHeader>
      <CardContent className="space-y-3">
        <div className="flex items-center justify-between">
          <div className="text-zinc-300 text-sm">Work through problems; state persists in localStorage.</div>
          <div className="flex gap-2">
            <Button variant="ghost" size="sm" onClick={resetAll}>
              Reset
            </Button>
          </div>
        </div>

        <div className="space-y-2">
          {items.map((p) => (
            <div key={p.id} className={`p-3 rounded border border-zinc-800 ${p.done ? "bg-zinc-800/30" : "bg-zinc-800/10"}`}>
              <div className="flex items-start gap-3">
                <input type="checkbox" checked={p.done} onChange={() => toggleDone(p.id)} className="mt-1 cursor-pointer" />
                <div className="flex-1">
                  <div className={`text-sm ${p.done ? "line-through text-zinc-400" : "text-zinc-100"}`}>{p.q}</div>
                  <div className="mt-2 flex items-center gap-2">
                    <Button variant="ghost" className="bg-white cursor-pointer" size="sm" onClick={() => toggleShow(p.id)}>{p.show ? "Hide hint" : "Show hint"}</Button>
                  </div>
                  {p.show && <div className="mt-2 text-zinc-400 text-sm">{p.hint}</div>}
                </div>
              </div>
            </div>
          ))}
        </div>
      </CardContent>
    </Card>
  );
}

//
// Assembled Chapter10 Page
//
export default function Chapter10() {
  const [matrix, setMatrix] = useState(() => defaultMatrix(3));

  // persist matrix to localStorage
  useEffect(() => {
    try {
      const raw = localStorage.getItem("chapter10_matrix");
      if (raw) setMatrix(JSON.parse(raw));
    } catch {}
  }, []);

  useEffect(() => {
    try {
      localStorage.setItem("chapter10_matrix", JSON.stringify(matrix));
    } catch {}
  }, [matrix]);

  return (
    <div className={`${theme.bg} p-6 min-h-screen space-y-6`}>
      <motion.header initial={{ opacity: 0, y: -8 }} animate={{ opacity: 1, y: 0 }}>
        <h1 className="text-3xl font-bold text-emerald-400">LU Decomposition and Matrix Inversion</h1>
        <p className="text-zinc-400 mt-2">covers LU decomposition, the matrix inverse, and an introduction to error analysis and system condition.</p>
      </motion.header>

      {/* Overview */}
      <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
        <Card className={theme.panel}>
          <CardHeader>
            <CardTitle className="text-cyan-300 flex items-center gap-2">
              <Beaker className="w-5 h-5" /> Overview
            </CardTitle>
          </CardHeader>
          <CardContent>
            <div className="text-zinc-300 text-sm">
              Practical algorithms for solving linear systems, computing inverses, and understanding numerical stability.
            </div>
          </CardContent>
        </Card>

        <Card className={theme.panel}>
          <CardHeader>
            <CardTitle className="text-emerald-300 flex items-center gap-2">
              <LineChart className="w-5 h-5" /> Applications
            </CardTitle>
          </CardHeader>
          <CardContent>
            <div className="text-zinc-300 text-sm">Use LU in engineering (circuit analysis, structural systems, discretized PDEs) for efficient repeated solving.</div>
          </CardContent>
        </Card>

        <Card className={theme.panel}>
          <CardHeader>
            <CardTitle className="text-amber-300 flex items-center gap-2">
              <Wand2 className="w-5 h-5" /> Numerical Notes
            </CardTitle>
          </CardHeader>
          <CardContent>
            <div className="text-zinc-300 text-sm">
              Pivoting improves stability. Condition numbers predict amplification of errors. Prefer factorization over explicit inversion for solving systems.
            </div>
          </CardContent>
        </Card>
      </div>

      {/* Sections */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        <div className="space-y-6">
          <Section101 matrix={matrix} setMatrix={setMatrix} />
          <Section103 matrix={matrix} />
        </div>
        <div className="space-y-6">
          <Section102 matrix={matrix} setMatrix={setMatrix} />
          <DocsPanel />
        </div>
      </div>

      <ProblemsPanel />

      {/* Decorative particle cloud */}
      <Card className={theme.panel}>
        <CardHeader>
          <CardTitle className="text-cyan-300 flex items-center gap-2">
            <Atom className="w-5 h-5" /> Visualization — Particle Cloud
          </CardTitle>
        </CardHeader>
        <CardContent>
          <div className="w-full h-48">
            <Canvas camera={{ position: [0, 0, 4] }}>
              <ambientLight intensity={0.6} />
              <pointLight position={[10, 10, 10]} />
              <ParticleCloud points={84} />
              <OrbitControls />
            </Canvas>
          </div>
          <div className="text-zinc-400 text-sm mt-3">Decorative 3D element for visual continuity across chapters.</div>
        </CardContent>
      </Card>
      <BottomBar/>
    </div>
  );
}
