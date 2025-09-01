"use client";

/**
 * src/pages/Chapter29.jsx
 *
 * Regenerated: responsive improvements, consistent spacing, and plot responsiveness.
 * Maintains original solver logic (Jacobi, Gauss-Seidel, SOR), BCs, control-volume demo, and problems.
 *
 * Requirements (same as before):
 *  - react, react-dom
 *  - react-plotly.js, plotly.js
 *  - framer-motion
 *  - react-markdown, remark-math, rehype-katex, katex
 *  - lucide-react
 *  - "@/components/ui/*" (Card, CardContent, Input, Button, Textarea) — adapt imports if different
 */

import React, { useEffect, useMemo, useRef, useState } from "react";
import Plot from "react-plotly.js";
import { motion } from "framer-motion";
import ReactMarkdown from "react-markdown";
import remarkMath from "remark-math";
import rehypeKatex from "rehype-katex";
import "katex/dist/katex.min.css";
import BottomBar from "../components/BottomBar";

import {
  BookOpen,
  Grid,
  Zap,
  Settings,
  Terminal,
  BarChart2,
  ListChecks,
  SunMedium,
} from "lucide-react";

// Replace or adjust these imports if your project layout differs
import { Card, CardContent } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Button } from "@/components/ui/button";
import { Textarea } from "@/components/ui/textarea";

/* ---------------------------------------------------------------------------
   Theme and small helpers
   --------------------------------------------------------------------------- */
const THEME = {
  BG: "bg-zinc-950",
  PANEL: "bg-zinc-900/60 border border-zinc-700",
  ACCENT: "#22d3ee",
  ACCENT2: "#34d399",
  WARN: "#f59e0b",
  DANGER: "#ef4444",
};

/* Markdown wrapper with KaTeX support */
function MD({ children, className = "" }) {
  const components = {
    p: ({ node, ...props }) => <p className={`text-sm text-zinc-300 leading-relaxed ${className}`} {...props} />,
    strong: ({ node, ...props }) => <strong className="text-emerald-300" {...props} />,
    em: ({ node, ...props }) => <em className="text-zinc-200" {...props} />,
    code: ({ node, inline, ...props }) =>
      inline ? (
        <code className="bg-zinc-800/70 px-1 rounded text-xs text-cyan-200" {...props} />
      ) : (
        <pre className="bg-zinc-800/60 rounded p-2 text-xs overflow-auto" {...props} />
      ),
    li: ({ node, ...props }) => <li className="ml-4 text-sm text-zinc-300" {...props} />,
    h1: ({ node, ...props }) => <h1 className="text-xl text-zinc-100 font-semibold" {...props} />,
    h2: ({ node, ...props }) => <h2 className="text-lg text-zinc-100 font-semibold" {...props} />,
  };

  return (
    <ReactMarkdown remarkPlugins={[remarkMath]} rehypePlugins={[rehypeKatex]} components={components}>
      {children}
    </ReactMarkdown>
  );
}

function SectionHeader({ icon: Icon, number, title, accent = THEME.ACCENT }) {
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

/* ---------------------------------------------------------------------------
   Numerical helpers / finite-difference solver
   --------------------------------------------------------------------------- */

/**
 * Build grid coordinates for N subdivisions -> (N+1) grid points [0..N]
 */
function buildGrid(N, a = 0, b = 1) {
  const h = (b - a) / N;
  const xs = [];
  for (let i = 0; i <= N; i++) xs.push(a + i * h);
  return { N, h, xs };
}

function createGrid(N, initVal = 0) {
  const arr = new Array(N + 1);
  for (let i = 0; i <= N; i++) {
    arr[i] = new Array(N + 1).fill(initVal);
  }
  return arr;
}

function copyGrid(grid) {
  return grid.map((row) => row.slice());
}

/**
 * Apply Dirichlet BCs: bc.left(y), bc.right(y), bc.bottom(x), bc.top(x)
 * Mutates provided grid.
 */
function applyDirichletBCs(grid, xs, N, bc) {
  if (!bc) return;
  for (let j = 0; j <= N; j++) {
    const y = xs[j];
    if (typeof bc.left === "function") grid[0][j] = bc.left(y);
    if (typeof bc.right === "function") grid[N][j] = bc.right(y);
    if (typeof bc.bottom === "function") grid[j][0] = bc.bottom(xs[j]);
    if (typeof bc.top === "function") grid[j][N] = bc.top(xs[j]);
  }
}

/**
 * Solve Laplace (∇² u = 0) with simple Neumann support (ghost points).
 *
 * options: { N, bc, method, tol, maxIter, omega, verbose }
 */
function solveLaplaceFD({ N = 40, bc = {}, method = "jacobi", tol = 1e-6, maxIter = 5000, omega = 1.5, verbose = false }) {
  const { h, xs } = buildGrid(N);
  const grid = createGrid(N, 0);
  applyDirichletBCs(grid, xs, N, bc);

  const residualHistory = [];
  let iter = 0;

  function laplaceUpdate(u, i, j) {
    return 0.25 * (u[i + 1][j] + u[i - 1][j] + u[i][j + 1] + u[i][j - 1]);
  }

  function applyNeumannGhosts(u) {
    if (!bc.neumann) return;
    const gLeft = bc.neumann.left;
    const gRight = bc.neumann.right;
    const gBottom = bc.neumann.bottom;
    const gTop = bc.neumann.top;

    for (let j = 0; j <= N; j++) {
      const y = xs[j];
      if (typeof gLeft === "function") u[0][j] = u[1][j] - h * gLeft(y);
      if (typeof gRight === "function") u[N][j] = u[N - 1][j] + h * gRight(y);
    }
    for (let i = 0; i <= N; i++) {
      const x = xs[i];
      if (typeof gBottom === "function") u[i][0] = u[i][1] - h * gBottom(x);
      if (typeof gTop === "function") u[i][N] = u[i][N - 1] + h * gTop(x);
    }
  }

  function computeResidual(u) {
    let maxRes = 0;
    for (let i = 1; i < N; i++) {
      for (let j = 1; j < N; j++) {
        const r = Math.abs(u[i + 1][j] + u[i - 1][j] + u[i][j + 1] + u[i][j - 1] - 4 * u[i][j]);
        if (r > maxRes) maxRes = r;
      }
    }
    return maxRes;
  }

  let uOld = copyGrid(grid);

  while (iter < maxIter) {
    iter += 1;

    if (method === "jacobi") {
      const uNew = copyGrid(grid);
      applyNeumannGhosts(uNew);

      for (let i = 1; i < N; i++) {
        for (let j = 1; j < N; j++) {
          uNew[i][j] = laplaceUpdate(uOld, i, j);
        }
      }

      applyDirichletBCs(uNew, xs, N, bc);
      uOld = uNew;
    } else if (method === "gs" || method === "sor") {
      applyNeumannGhosts(grid);

      for (let i = 1; i < N; i++) {
        for (let j = 1; j < N; j++) {
          const sigma = grid[i + 1][j] + grid[i - 1][j] + grid[i][j + 1] + grid[i][j - 1];
          const newVal = 0.25 * sigma;
          if (method === "gs") {
            grid[i][j] = newVal;
          } else {
            grid[i][j] = grid[i][j] + omega * (newVal - grid[i][j]);
          }
        }
      }
      applyDirichletBCs(grid, xs, N, bc);
      uOld = copyGrid(grid);
    } else {
      throw new Error("Unknown method: " + method);
    }

    const res = computeResidual(uOld);
    residualHistory.push(res);

    if (res < tol) break;
    // safety: prevent infinite loops
    if (iter >= maxIter) break;
  }

  const outputGrid = method === "jacobi" ? uOld : grid;
  return { u: outputGrid, iterations: iter, residualHistory };
}

/* ---------------------------------------------------------------------------
   Analytic tests and error metrics
   --------------------------------------------------------------------------- */
function analyticTestSolution(name, x, y) {
  if (name === "sinhsin") {
    return Math.sin(Math.PI * x) * Math.sinh(Math.PI * y) / Math.sinh(Math.PI);
  }
  if (name === "zero") return 0;
  return 0;
}

function computeLinftyError(uGrid, N, analyticName) {
  const { xs } = buildGrid(N);
  let maxErr = 0;
  for (let i = 0; i <= N; i++) {
    for (let j = 0; j <= N; j++) {
      const x = xs[i];
      const y = xs[j];
      const ua = analyticTestSolution(analyticName, x, y);
      const err = Math.abs(uGrid[i][j] - ua);
      if (err > maxErr) maxErr = err;
    }
  }
  return maxErr;
}

/* ---------------------------------------------------------------------------
   Convert grid to Plotly arrays
   --------------------------------------------------------------------------- */
function gridToPlotly(uGrid, N, a = 0, b = 1) {
  const xs = [];
  const ys = [];
  const h = (b - a) / N;
  for (let i = 0; i <= N; i++) {
    xs.push(a + i * h);
    ys.push(a + i * h);
  }
  const Z = new Array(N + 1);
  for (let j = 0; j <= N; j++) {
    Z[j] = new Array(N + 1);
    for (let i = 0; i <= N; i++) {
      Z[j][i] = uGrid[i][j];
    }
  }
  return { xs, ys, Z };
}

/* ---------------------------------------------------------------------------
   ResponsivePlot wrapper (keeps plot inside responsive card)
   --------------------------------------------------------------------------- */
function ResponsivePlot({ children, height = 360 }) {
  return (
    <div className="w-full max-w-full overflow-hidden rounded-lg">
      <div style={{ minHeight: Math.max(120, height) }}>{children}</div>
    </div>
  );
}

/* ---------------------------------------------------------------------------
   Section 29.1: Intro
   --------------------------------------------------------------------------- */
function Section291() {
  return (
    <motion.div initial={{ opacity: 0, y: 8 }} whileInView={{ opacity: 1, y: 0 }} className={`${THEME.PANEL} p-4 rounded-2xl`}>
      <div className="flex items-center justify-between">
        <SectionHeader icon={Grid} number="29.1" title="The Laplace Equation" />
        <div className="text-xs text-zinc-400">Equation: ∇²u = 0 (elliptic)</div>
      </div>

      <div className="grid gap-4 grid-cols-1 md:grid-cols-3 mt-4">
        <div className="md:col-span-1 space-y-3">
          <MD>
  {`The Laplace equation
$$\\nabla^2 u(x,y) = 0$$
models steady-state potentials (temperature, electrostatics, streamfunction) where no internal sources exist.`}
</MD>

<MD>
  {`**Boundary conditions** determine the solution uniquely for elliptic equations. Typical BCs:

- Dirichlet: specify $u$ on boundary.
- Neumann: specify $\\frac{\\partial u}{\\partial n}$ (flux) on boundary.
- Mixed (Robin): linear combination.`}
</MD>

        </div>

        <div className="md:col-span-2 space-y-3">
          <div className="rounded-lg border border-zinc-700 p-3 bg-zinc-900/50">
           <MD>
  {`**Analytic test**: choose a manufactured solution such as 
$u(x,y) = \\frac{\\sinh(\\pi y) \\sin(\\pi x)}{\\sinh(\\pi)}$ 
which satisfies Laplace on the unit square with Dirichlet BCs.`}
</MD>

          </div>
        </div>
      </div>
    </motion.div>
  );
}

/* ---------------------------------------------------------------------------
   Section 29.2: Solver UI + Plots
   --------------------------------------------------------------------------- */
function Section292() {
  const [N, setN] = useState(40);
  const [method, setMethod] = useState("sor");
  const [tol, setTol] = useState(1e-6);
  const [maxIter, setMaxIter] = useState(5000);
  const [omega, setOmega] = useState(1.5);
  const [initGuess, setInitGuess] = useState(0.0);
  const [runNow, setRunNow] = useState(false);

  const [bcType] = useState("dirichlet");
  const [analyticName, setAnalyticName] = useState("sinhsin");

  const bcDirichlet = {
    left: (y) => 0,
    right: (y) => 0,
    bottom: (x) => 0,
    top: (x) => Math.sin(Math.PI * x),
  };

  const [result, setResult] = useState(null);
  const [running, setRunning] = useState(false);

  useEffect(() => {
    if (!runNow) return;
    (function runSolver() {
      setRunning(true);
      const Ncl = Math.max(8, Math.min(200, Math.floor(N)));
      const options = { N: Ncl, bc: bcDirichlet, method, tol: parseFloat(tol), maxIter: parseInt(maxIter, 10), omega: parseFloat(omega) };
      const gridInit = createGrid(Ncl, parseFloat(initGuess));
      const { xs } = buildGrid(Ncl);
      applyDirichletBCs(gridInit, xs, Ncl, bcDirichlet);

      // run solver synchronously (JS single thread) but we show spinner state
      const res = solveLaplaceFD(options);
      const linf = computeLinftyError(res.u, Ncl, analyticName);
      setResult({ res, N: Ncl, linf });
      setRunning(false);
      setRunNow(false);
    })();
  }, [runNow]); // eslint-disable-line

  const plotData = useMemo(() => {
    if (!result) return null;
    const { u } = result.res;
    const Nlocal = result.N;
    return gridToPlotly(u, Nlocal);
  }, [result]);

  return (
    <motion.div initial={{ opacity: 0, y: 8 }} whileInView={{ opacity: 1, y: 0 }} className={`${THEME.PANEL} p-4 rounded-2xl`}>
      <div className="flex items-center justify-between">
        <SectionHeader icon={Settings} number="29.2" title="Solution Technique (Jacobi, GS, SOR)" />
        <div className="text-xs text-zinc-400">Iterative solvers & convergence</div>
      </div>

      <div className="grid gap-4 grid-cols-1 md:grid-cols-4 mt-4">
        <div className="md:col-span-1 space-y-3">
          <MD>{`Choose grid resolution and method. SOR uses relaxation factor ω (1<ω<2 typically).`}</MD>

          <div className="grid grid-cols-1 gap-2">
            <div>
              <label className="text-xs text-zinc-400">N (subdivisions)</label>
              <Input
                type="number"
                value={N}
                className="text-white"
                onChange={(e) => {
                  const v = parseInt(e.target.value || 8, 10);
                  setN(isNaN(v) ? 8 : Math.max(8, v));
                  
                }}
              />
            </div>

            <div>
              <label className="text-xs text-zinc-400">Method</label>
              <select value={method} onChange={(e) => setMethod(e.target.value)} className="bg-zinc-800/90 text-zinc-100 rounded border border-zinc-700 p-2 w-full">
                <option value="jacobi">Jacobi</option>
                <option value="gs">Gauss-Seidel</option>
                <option value="sor">SOR</option>
              </select>
            </div>

            <div className="grid grid-cols-2 gap-2">
              <div>
                <label className="text-xs text-zinc-400">tol</label>
                <Input
                  type="number"
                  value={tol}
                   className="text-white"
                  onChange={(e) => setTol(parseFloat(e.target.value || 1e-6))}
                />
              </div>

              <div>
                <label className="text-xs text-zinc-400">maxIter</label>
                <Input
                  type="number"
                  value={maxIter}
                   className="text-white"
                  onChange={(e) => setMaxIter(parseInt(e.target.value || 5000, 10))}
                />
              </div>
            </div>

            <div className="grid grid-cols-2 gap-2">
              <div>
                <label className="text-xs text-zinc-400">ω (for SOR)</label>
                <Input  className="text-white" type="number" value={omega} onChange={(e) => setOmega(parseFloat(e.target.value || 1.5))} />
              </div>

              <div>
                <label className="text-xs text-zinc-400">Init guess</label>
                <Input  className="text-white" type="number" value={initGuess} onChange={(e) => setInitGuess(parseFloat(e.target.value || 0))} />
              </div>
            </div>
          </div>

          <div className="flex gap-2 mt-2">
            <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={() => setRunNow(true)}>{running ? "Running..." : "Run solver"}</Button>
            <Button className="bg-white cursor-pointer" variant="ghost" onClick={() => { setResult(null); setRunNow(false); }}>Reset</Button>
          </div>

          <div className="text-xs text-zinc-400 mt-2">Analytic test:</div>
          <select value={analyticName} onChange={(e) => setAnalyticName(e.target.value)} className="bg-zinc-800/90 text-zinc-100 rounded border border-zinc-700 p-2 w-full mt-1">
            <option value="sinhsin">sin(pi x) sinh(pi y)/sinh(pi)</option>
            <option value="zero">Zero solution</option>
          </select>
        </div>

        <div className="md:col-span-3 space-y-3">
          <div className="rounded-lg border border-zinc-700 p-3 bg-zinc-900/50">
            {!plotData && (
              <div className="text-sm text-zinc-400">No result yet — adjust parameters and click "Run solver".</div>
            )}

            {plotData && result && (
              <>
                <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
                  <div className="rounded-lg border border-zinc-700 p-2 bg-zinc-950/20">
                    <ResponsivePlot height={360}>
                      <Plot
                        data={[
                          {
                            z: plotData.Z,
                            x: plotData.xs,
                            y: plotData.ys,
                            type: "surface",
                            colorscale: "Viridis",
                            showscale: true,
                          },
                        ]}
                        layout={{
                          title: "Numerical solution (surface)",
                          autosize: true,
                          paper_bgcolor: "rgba(0,0,0,0)",
                          plot_bgcolor: "rgba(0,0,0,0)",
                          margin: { t: 36, b: 40, l: 40, r: 10 },
                        }}
                        style={{ width: "100%", height: 360 }}
                        useResizeHandler
                        config={{ responsive: true }}
                      />
                    </ResponsivePlot>
                  </div>

                  <div className="rounded-lg border border-zinc-700 p-2 bg-zinc-950/20">
                    <ResponsivePlot height={360}>
                      <Plot
                        data={[
                          {
                            z: plotData.Z,
                            x: plotData.xs,
                            y: plotData.ys,
                            type: "contour",
                            colorscale: "Viridis",
                            contours: { coloring: "fill" },
                          },
                        ]}
                        layout={{
                          title: "Numerical solution (contour)",
                          autosize: true,
                          paper_bgcolor: "rgba(0,0,0,0)",
                          plot_bgcolor: "rgba(0,0,0,0)",
                          margin: { t: 36, b: 40, l: 40, r: 10 },
                        }}
                        style={{ width: "100%", height: 360 }}
                        useResizeHandler
                        config={{ responsive: true }}
                      />
                    </ResponsivePlot>
                  </div>
                </div>

                <div className="mt-3 text-sm text-zinc-300 grid grid-cols-1 md:grid-cols-2 gap-3">
                  <div>Iterations: <span className="text-zinc-100 font-medium">{result.res.iterations}</span></div>
                  <div>Linfty error vs analytic: <span className="text-zinc-100 font-medium">{result.linf.toExponential(3)}</span></div>
                </div>

                <div className="mt-3">
                  <h4 className="text-sm text-zinc-100 mb-2">Residual history</h4>
                  <ResponsivePlot height={200}>
                    <Plot
                      data={[{ x: result.res.residualHistory.map((_, i) => i + 1), y: result.res.residualHistory, mode: "lines", line: { color: THEME.ACCENT } }]}
                      layout={{
                        title: "",
                        autosize: true,
                        paper_bgcolor: "rgba(0,0,0,0)",
                        plot_bgcolor: "rgba(0,0,0,0)",
                        xaxis: { title: "iteration", color: "#9ca3af" },
                        yaxis: { title: "residual (∞-norm)", type: "log", color: "#9ca3af" },
                        margin: { t: 6, b: 30, l: 50, r: 10 },
                      }}
                      style={{ width: "100%", height: 200 }}
                      useResizeHandler
                      config={{ responsive: true }}
                    />
                  </ResponsivePlot>
                </div>
              </>
            )}
          </div>
        </div>
      </div>
    </motion.div>
  );
}

/* ---------------------------------------------------------------------------
   Section 29.3: Boundary conditions interactive
   --------------------------------------------------------------------------- */
function Section293() {
  const [N, setN] = useState(50);
  const [method, setMethod] = useState("sor");
  const [omega, setOmega] = useState(1.6);
  const [mode, setMode] = useState("dirichletSquare");
  const [run, setRun] = useState(false);
  const [result, setResult] = useState(null);
  const [running, setRunning] = useState(false);

  useEffect(() => {
    if (!run) return;
    (function runScenario() {
      setRunning(true);
      const Ncl = Math.max(8, Math.min(160, Math.floor(N)));
      let bc = {};
      if (mode === "dirichletSquare") {
        bc = {
          left: (y) => 0,
          right: (y) => 0,
          bottom: (x) => 0,
          top: (x) => Math.sin(Math.PI * x),
        };
      } else if (mode === "neumannFlux") {
        bc = {
          left: (y) => Math.sin(Math.PI * y),
          right: (y) => 0,
          bottom: null,
          top: null,
          neumann: {
            bottom: (x) => 0,
            top: (x) => 0,
          },
        };
      } else {
        bc = {
          left: (y) => 0,
          bottom: (x) => 0,
          top: (x) => 1,
          right: null,
          neumann: { right: (y) => 0 },
        };
      }
      const res = solveLaplaceFD({ N: Ncl, bc, method, tol: 1e-7, maxIter: 5000, omega });
      setResult({ res, N: Ncl, bcMode: mode });
      setRunning(false);
      setRun(false);
    })();
  }, [run]); // eslint-disable-line

  const plotData = useMemo(() => {
    if (!result) return null;
    return gridToPlotly(result.res.u, result.N);
  }, [result]);

  return (
    <motion.div initial={{ opacity: 0, y: 8 }} whileInView={{ opacity: 1, y: 0 }} className={`${THEME.PANEL} p-4 rounded-2xl`}>
      <div className="flex items-center justify-between">
        <SectionHeader icon={Zap} number="29.3" title="Boundary Conditions: Dirichlet, Neumann, Mixed" accent={THEME.WARN} />
        <div className="text-xs text-zinc-400">See how BCs change solution structure</div>
      </div>

      <div className="grid gap-4 grid-cols-1 md:grid-cols-4 mt-4">
        <div className="md:col-span-1 space-y-3">
          <MD>{`Pick a boundary scenario and run the solver. Neumann conditions are approximated with ghost points (first-order).`}</MD>

          <div className="grid gap-2">
            <label className="text-xs text-zinc-400">N</label>
            <Input  className="text-white" type="number" value={N} onChange={(e) => { const v = parseInt(e.target.value || 8, 10); setN(isNaN(v) ? 8 : Math.max(8, v)); }} />
            <label className="text-xs text-zinc-400">Method</label>
            <select value={method} onChange={(e) => setMethod(e.target.value)} className="bg-zinc-800/90 text-zinc-100 rounded border border-zinc-700 p-2 w-full">
              <option value="gs">Gauss-Seidel</option>
              <option value="jacobi">Jacobi</option>
              <option value="sor">SOR</option>
            </select>
            <label className="text-xs text-zinc-400">ω (for SOR)</label>
            <Input  className="text-white" type="number" value={omega} onChange={(e) => setOmega(parseFloat(e.target.value || 1.5))} />
            <label className="text-xs text-zinc-400">Scenario</label>
            <select value={mode} onChange={(e) => setMode(e.target.value)} className="bg-zinc-800/90 text-zinc-100 rounded border border-zinc-700 p-2 w-full">
              <option value="dirichletSquare">Dirichlet: sin top</option>
              <option value="neumannFlux">Neumann: insulated top/bottom</option>
              <option value="mixedCorner">Mixed: top=1, others mixed</option>
            </select>

            <div className="flex gap-2">
              <Button className="cursor-pointer bg-emerald-300 hover:bg-emerald-400 text-black" onClick={() => setRun(true)}>{running ? "Running..." : "Run scenario"}</Button>
              <Button className="bg-white cursor-pointer" variant="ghost" onClick={() => setResult(null)}>Clear</Button>
            </div>
          </div>
        </div>

        <div className="md:col-span-3 space-y-3">
          <div className="rounded-lg border border-zinc-700 p-3 bg-zinc-900/50">
            {!plotData && <div className="text-sm text-zinc-400">No scenario solved yet.</div>}
            {plotData && result && (
              <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
                <div className="rounded-lg border border-zinc-700 p-2 bg-zinc-950/20">
                  <ResponsivePlot height={360}>
                    <Plot
                      data={[{ x: plotData.xs, y: plotData.ys, z: plotData.Z, type: "surface", colorscale: "Viridis", showscale: true }]}
                      layout={{
                        title: `Solution surface (${result.bcMode}, N=${result.N})`,
                        autosize: true,
                        paper_bgcolor: "rgba(0,0,0,0)",
                        plot_bgcolor: "rgba(0,0,0,0)",
                        margin: { t: 36, b: 40, l: 40, r: 10 },
                      }}
                      style={{ width: "100%", height: 360 }}
                      useResizeHandler
                      config={{ responsive: true }}
                    />
                  </ResponsivePlot>
                </div>

                <div className="rounded-lg border border-zinc-700 p-2 bg-zinc-950/20">
                  <ResponsivePlot height={360}>
                    <Plot
                      data={[{ x: plotData.xs, y: plotData.ys, z: plotData.Z, type: "contour", colorscale: "Viridis", contours: { coloring: "fill" } }]}
                      layout={{
                        title: "Contour",
                        autosize: true,
                        paper_bgcolor: "rgba(0,0,0,0)",
                        plot_bgcolor: "rgba(0,0,0,0)",
                        margin: { t: 36, b: 40, l: 40, r: 10 },
                      }}
                      style={{ width: "100%", height: 360 }}
                      useResizeHandler
                      config={{ responsive: true }}
                    />
                  </ResponsivePlot>
                </div>
              </div>
            )}
          </div>
        </div>
      </div>
    </motion.div>
  );
}

/* ---------------------------------------------------------------------------
   Section 29.4: Control-Volume Approach
   --------------------------------------------------------------------------- */
function Section294() {
  const [N, setN] = useState(32);
  const [showCells, setShowCells] = useState(true);

  const bc = { left: () => 0, right: () => 0, bottom: () => 0, top: () => 1 };
  const uResult = useMemo(() => {
    const Ncl = Math.max(8, Math.min(200, N));
    return solveLaplaceFD({ N: Ncl, bc, method: "sor", tol: 1e-6, maxIter: 5000, omega: 1.6 });
  }, [N]);

  const plotData = useMemo(() => {
    if (!uResult) return null;
    return gridToPlotly(uResult.u, N);
  }, [uResult, N]);

  return (
    <motion.div className={`${THEME.PANEL} p-4 rounded-2xl`} initial={{ opacity: 0, y: 8 }} whileInView={{ opacity: 1, y: 0 }}>
      <div className="flex items-center justify-between">
        <SectionHeader icon={BarChart2} number="29.4" title="Control-Volume Approach (cell-centered)" accent="#f97316" />
        <div className="text-xs text-zinc-400">Cell-centered vs node-centered discretizations</div>
      </div>

      <div className="grid gap-4 grid-cols-1 md:grid-cols-3 mt-4">
        <div className="md:col-span-1 space-y-3">
          <MD>
            {"Finite-volume places unknowns at cell centers and enforces flux balance over small control volumes. This yields similar algebraic equations to finite-difference for structured grids."}
          </MD>

          <div className="grid gap-2">
            <label className="text-xs text-zinc-400">N</label>
            <Input className="text-white"  type="number" value={N} onChange={(e) => { const v = parseInt(e.target.value || 8, 10); setN(isNaN(v) ? 8 : Math.max(8, Math.min(160, v))); }} />
            <div className="flex gap-2">
              <button className="text-xs px-2 py-1 bg-emerald-400 cursor-pointer rounded border" onClick={() => setShowCells((s) => !s)}>
                {showCells ? "Hide cells" : "Show cells"}
              </button>
            </div>
          </div>
        </div>

        <div className="md:col-span-2 space-y-3">
          <div className="rounded-lg border border-zinc-700 p-3 bg-zinc-900/50">
            {plotData && (
              <ResponsivePlot height={420}>
                <Plot
                  data={[
                    { z: plotData.Z, x: plotData.xs, y: plotData.ys, type: "surface", colorscale: "Cividis", showscale: true },
                    {
                      x: plotData.xs.flatMap((x) => plotData.xs),
                      y: plotData.ys.flatMap((y) => plotData.ys.map(() => y)).slice(0, (N + 1) * (N + 1)),
                      z: plotData.Z.flat(),
                      type: "scatter3d",
                      mode: "markers",
                      marker: { size: 2, color: "rgba(255,255,255,0.2)" },
                      visible: showCells ? true : "legendonly",
                    },
                  ]}
                  layout={{
                    title: "Solution surface with optional cell-centers overlay",
                    autosize: true,
                    paper_bgcolor: "rgba(0,0,0,0)",
                    plot_bgcolor: "rgba(0,0,0,0)",
                    margin: { t: 30 },
                  }}
                  style={{ width: "100%", height: 420 }}
                  useResizeHandler
                  config={{ responsive: true }}
                />
              </ResponsivePlot>
            )}
            {!plotData && <div className="text-sm text-zinc-400">Run solver in Section 29.2 to produce solution.</div>}
          </div>
        </div>
      </div>
    </motion.div>
  );
}

/* ---------------------------------------------------------------------------
   Section 29.5: Software notes
   --------------------------------------------------------------------------- */
function Section295() {
  const [showSnippet] = useState(true);
  const pythonSnippet = `# Solve Laplace in 2D using scipy (solve_bvp or sparse linear system)
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

N = 50
h = 1.0 / N
# build 5-point Laplacian (interior unknowns M=(N-1)^2) ...
# Use Dirichlet BCs to form RHS, then solve sparse system A u = b
`;

  const matlabSnippet = `% MATLAB: use built-in PDE toolbox or build sparse Laplacian
N = 50;
h = 1/N;
% build sparse matrix using spdiags and solve
`;

  return (
    <motion.div className={`${THEME.PANEL} p-4 rounded-2xl`} initial={{ opacity: 0, y: 8 }} whileInView={{ opacity: 1, y: 0 }}>
      <div className="flex items-center justify-between">
        <SectionHeader icon={Terminal} number="29.5" title="Software to Solve Elliptic Equations" accent={THEME.ACCENT2} />
        <div className="text-xs text-zinc-400">Python / MATLAB / Specialized libs</div>
      </div>

      <div className="grid gap-4 grid-cols-1 md:grid-cols-3 mt-4">
        <div className="md:col-span-1 space-y-3">
          <MD>{`Common approaches:\n\n- Direct sparse solves (sparse LU) for moderate size.\n- Iterative solvers (CG, GMRES) with preconditioners for large systems.\n- Multigrid for optimal complexity.`}</MD>
        </div>

        <div className="md:col-span-2 space-y-3">
          <div className="rounded-lg border border-zinc-700 p-2 bg-zinc-900/50">
            <div className="flex items-center gap-3 mb-2">
              <SunMedium className="w-5 h-5 text-cyan-300" />
              <div className="text-sm text-zinc-100 font-semibold">Code snippets</div>
            </div>

            <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
              <div className="rounded bg-zinc-950/30 overflow-auto p-3 border border-zinc-800">
                <div className="text-xs text-zinc-400 mb-2">Python (SciPy sparse)</div>
                <pre className="text-xs text-zinc-200 overflow-auto bg-transparent p-2 rounded">{pythonSnippet}</pre>
              </div>

              <div className="rounded bg-zinc-950/30 p-3 border border-zinc-800">
                <div className="text-xs text-zinc-400 mb-2">MATLAB</div>
                <pre className="text-xs text-zinc-200 overflow-auto bg-transparent p-2 rounded">{matlabSnippet}</pre>
              </div>
            </div>

            <MD>{`**Recommendation:** start with a sparse direct solve for moderate N to validate discretization, then profile and switch to iterative/preconditioned solvers or multigrid for large-scale problems.`}</MD>
          </div>
        </div>
      </div>
    </motion.div>
  );
}

/* ---------------------------------------------------------------------------
   Problems section
   --------------------------------------------------------------------------- */
function Problems29() {
  const [open, setOpen] = useState(-1);
  const problems = [
    {
      id: 1,
      title: "Laplace on square with sin top BC",
      prompt: "Solve ∇²u=0 on (0,1)×(0,1), u(x,0)=0, u(x,1)=sin(π x), u(0,y)=u(1,y)=0. Compare numeric vs analytic manufactured solution.",
      solution: "Use the 5-point stencil and SOR; compute Linfty error and show convergence with grid refinement.",
    },
    {
      id: 2,
      title: "Neumann BC (insulated sides)",
      prompt: "Simulate Laplace with insulated (∂u/∂n=0) on top/bottom and Dirichlet left/right. How to implement ghost points?",
      solution: "Approximate ghost using forward/backward differences: u_ghost = u_interior ± h g, derived from du/dn ≈ (u_ghost - u_interior)/h.",
    },
    {
      id: 3,
      title: "Control-Volume discretization",
      prompt: "Derive finite-volume discretization on cell-centered grid and show equivalence to 5-point stencil on uniform mesh.",
      solution: "Balance fluxes across faces; central difference approximations recover the 5-point scheme.",
    },
    {
      id: 4,
      title: "Performance: Jacobi vs SOR vs GS",
      prompt: "Time the convergence for N=100 and compare number of iterations to reach tol for each method. Discuss relaxation parameter choice.",
      solution: "Jacobi converges slowly; GS faster; SOR with optimal ω significantly accelerates convergence (optimal ω depends on spectral radius / mode).",
    },
  ];

  return (
    <motion.div className={`${THEME.PANEL} p-4 rounded-2xl`} initial={{ opacity: 0, y: 8 }} whileInView={{ opacity: 1, y: 0 }}>
      <div className="flex items-center justify-between">
        <SectionHeader icon={ListChecks} number="Problems" title="Exercises" />
        <div className="text-xs text-zinc-400">Solutions sketched below</div>
      </div>

      <div className="mt-4 space-y-3">
        {problems.map((p, i) => (
          <div key={p.id} className="rounded-lg border border-zinc-700 p-3 bg-zinc-900/50">
            <div className="flex items-center justify-between">
              <div>
                <div className="text-sm text-zinc-100 font-medium">{p.title}</div>
                <div className="text-xs text-zinc-400">{p.prompt}</div>
              </div>
              <div>
                <Button variant="ghost" className="bg-white cursor-pointer" onClick={() => setOpen(open === i ? -1 : i)}>
                  {open === i ? "Hide" : "Show"}
                </Button>
              </div>
            </div>
            {open === i && (
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

/* ---------------------------------------------------------------------------
   Main chapter assembly
   --------------------------------------------------------------------------- */
export default function Chapter29() {
  useEffect(() => {
    if (typeof window !== "undefined") window.scrollTo({ top: 0, behavior: "smooth" });
  }, []);

  return (
    <div className={`p-4 sm:p-6 lg:p-8 space-y-6 ${THEME.BG} min-h-screen`}>
      <motion.div initial={{ opacity: 0, y: -8 }} animate={{ opacity: 1, y: 0 }}>
        <div className="flex flex-col sm:flex-row sm:items-center sm:justify-between gap-4">
          <div>
            
            <h1 className="text-2xl sm:text-3xl font-bold text-emerald-400">Finite Difference: Elliptic Equations</h1>
            <div className="text-zinc-400 mt-1">29.1 The Laplace Equation • 29.2 Solution Technique • 29.3 Boundary Conditions • 29.4 Control-Volume Approach • 29.5 Software • Problems</div>
          </div>


        </div>
      </motion.div>

      <div className="space-y-6">
        <Section291 />
        <Section292 />
        <Section293 />
        <Section294 />
        <Section295 />
        <Problems29 />
        <BottomBar/>
      </div>

      <footer className="text-xs text-zinc-500 mt-6">
        Regenerated Chapter 29 module — interactive finite-difference solver and responsive visualizations.
      </footer>
    </div>
  );
}
