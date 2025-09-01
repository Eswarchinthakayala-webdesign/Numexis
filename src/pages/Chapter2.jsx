// src/pages/Chapter2_NumericalMethods.jsx
import React, { useState, useMemo, useRef, useEffect, Suspense } from 'react';
import { motion } from 'framer-motion';
import Plot from 'react-plotly.js';
import { Canvas, useFrame } from '@react-three/fiber';
import { OrbitControls, Html } from '@react-three/drei';
import * as THREE from 'three';
import { ChevronDown, ChevronRight, BookOpen } from 'lucide-react';
import { Card, CardHeader, CardTitle, CardContent } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import BottomBar from '../components/BottomBar';

// ---------------- Theme: zinc-dark inspired ----------------
const theme = {
  background: 'bg-zinc-900',
  panelBg: '#0f1720', // slightly darker than zinc-900
  text: '#e6eef6',
  accent1: '#60a5fa', // blue-400
  accent2: '#34d399', // emerald-400
  accent3: '#f59e0b', // amber-500
};

// ---------------- Framer motion helpers ----------------
const cardMotion = {
  initial: { opacity: 0, y: 14, scale: 0.995 },
  enter: { opacity: 1, y: 0, scale: 1, transition: { duration: 0.45, ease: [0.2, 0.8, 0.2, 1] } },
  hover: { scale: 1.02, y: -6, transition: { duration: 0.25 } },
};

// ---------------- Small utility: parse user function safely-ish ----------------
function parseFn(fnStr) {
  try {
    // create a function of x that returns the expression
    // NOTE: this uses the Function constructor; only use for trusted/sandboxed inputs
    // we wrap Math to provide common math functions
    // allow 'x' and 't' variables
    return new Function('x', 't', `with (Math) { return (${fnStr}); }`);
  } catch (e) {
    return null;
  }
}

// ---------------- Rotating responsive 3D object ----------------
function ResponsiveCube({ color = theme.accent1 }) {
  const ref = useRef();
  useFrame((state, delta) => {
    if (!ref.current) return;
    // make motion respond to cursor position for interactivity
    const { mouse } = state;
    ref.current.rotation.x = THREE.MathUtils.lerp(ref.current.rotation.x, mouse.y * 0.5, 0.05);
    ref.current.rotation.y += delta * 0.6;
  });
  return (
    <mesh ref={ref} castShadow>
      <boxGeometry args={[1.2, 1.2, 1.2]} />
      <meshStandardMaterial color={color} metalness={0.4} roughness={0.25} />
    </mesh>
  );
}

// ---------------- Numerical methods implementations ----------------
function bisection(f, a, b, tol = 1e-6, maxIter = 100) {
  const res = [];
  let fa = f(a, 0), fb = f(b, 0);
  if (fa === undefined || fb === undefined) return res;
  if (fa * fb > 0) return res;
  let left = a, right = b;
  for (let i = 0; i < maxIter; i++) {
    const mid = (left + right) / 2;
    const fm = f(mid, 0);
    res.push({ iter: i + 1, left, right, mid, fm });
    if (Math.abs(fm) < tol || (right - left) / 2 < tol) break;
    if (fa * fm <= 0) right = mid; else left = mid, fa = fm;
  }
  return res;
}

function newtonRaphson(f, df, x0, tol = 1e-6, maxIter = 50) {
  const res = [];
  let x = x0;
  for (let i = 0; i < maxIter; i++) {
    const y = f(x, 0);
    const dy = df ? df(x, 0) : null;
    res.push({ iter: i + 1, x, y });
    if (Math.abs(y) < tol) break;
    if (!dy) break;
    x = x - y / dy;
  }
  return res;
}

function trapezoidal(f, a, b, n) {
  const h = (b - a) / n;
  let s = 0.5 * (f(a, 0) + f(b, 0));
  for (let i = 1; i < n; i++) s += f(a + i * h, 0);
  return h * s;
}

function simpson(f, a, b, n) {
  if (n % 2 === 1) n++; // n must be even
  const h = (b - a) / n;
  let s = f(a, 0) + f(b, 0);
  for (let i = 1; i < n; i++) {
    s += f(a + i * h, 0) * (i % 2 === 0 ? 2 : 4);
  }
  return (h / 3) * s;
}

function finiteDiff(f, x, h = 1e-4) {
  // central difference
  return (f(x + h, 0) - f(x - h, 0)) / (2 * h);
}

function explicitEuler(dt, tMax, f, y0 = 0) {
  const result = [];
  let y = y0;
  for (let t = 0; t <= tMax + 1e-9; t += dt) {
    result.push({ t, y });
    y += dt * f(y, t);
  }
  return result;
}

function rungeKutta4(dt, tMax, f, y0 = 0) {
  const result = [];
  let y = y0;
  for (let t = 0; t <= tMax + 1e-9; t += dt) {
    result.push({ t, y });
    const k1 = f(y, t);
    const k2 = f(y + dt * k1 / 2, t + dt / 2);
    const k3 = f(y + dt * k2 / 2, t + dt / 2);
    const k4 = f(y + dt * k3, t + dt);
    y += dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
  }
  return result;
}

function gaussElimination(A, b) {
  const n = A.length;
  const M = A.map((r) => r.slice());
  const B = b.slice();
  for (let k = 0; k < n; k++) {
    // partial pivot
    let max = k;
    for (let i = k + 1; i < n; i++) if (Math.abs(M[i][k]) > Math.abs(M[max][k])) max = i;
    [M[k], M[max]] = [M[max], M[k]];
    [B[k], B[max]] = [B[max], B[k]];
    for (let i = k + 1; i < n; i++) {
      const factor = M[i][k] / M[k][k];
      for (let j = k; j < n; j++) M[i][j] -= factor * M[k][j];
      B[i] -= factor * B[k];
    }
  }
  // back substitution
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    let sum = B[i];
    for (let j = i + 1; j < n; j++) sum -= M[i][j] * x[j];
    x[i] = sum / M[i][i];
  }
  return x;
}

// ---------------- UI Sections ----------------
function SectionRootFinding() {
  const [fnText, setFnText] = useState('Math.sin(x) - 0.5');
  const [a, setA] = useState(0);
  const [b, setB] = useState(2);
  const [guess, setGuess] = useState(1);
  const [bisectionRes, setBisectionRes] = useState([]);
  const [nrRes, setNrRes] = useState([]);

  const f = useMemo(() => parseFn(fnText), [fnText]);
  const df = useMemo(() => parseFn('(x) => ((' + fnText + ') - ((' + fnText + ')))'), [fnText]);
  // NOTE: user can provide explicit derivative if they want by editing the code below

  useEffect(() => {
    if (!f) return;
    try {
      const bis = bisection((x) => f(x), Number(a), Number(b), 1e-6, 80);
      setBisectionRes(bis);
    } catch (e) { setBisectionRes([]); }
    try {
      // try to compute derivative numerically
      const nr = newtonRaphson((x) => f(x), (x) => finiteDiff((xx) => f(xx), x), Number(guess), 1e-7, 30);
      setNrRes(nr);
    } catch (e) { setNrRes([]); }
  }, [f, a, b, guess]);

  const xPlot = useMemo(() => {
    const arr = [];
    if (!f) return arr;
    for (let x = Number(a) - 1; x <= Number(b) + 1; x += 0.02) arr.push(x);
    return arr;
  }, [f, a, b]);

  const yPlot = xPlot.map((x) => {
    try { return f(x); } catch { return null; }
  });

  return (
    <motion.div variants={cardMotion} initial="initial" animate="enter" whileHover="hover">
      <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700`}> 
        <CardHeader>
          <CardTitle className="text-blue-300">Root Finding — Bisection & Newton-Raphson</CardTitle>
        </CardHeader>
        <CardContent>
          <p className="text-zinc-200 mb-2">Enter a function of <code className="font-mono">x</code>. Example: <code className="font-mono">Math.sin(x) - 0.5</code></p>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-3 mb-3">
            <Input value={fnText} onChange={(e) => setFnText(e.target.value)} className="col-span-2 bg-zinc-600 text-white" />
            <div className="flex gap-2">
              <Input  className="bg-zinc-600 text-white" value={a} onChange={(e) => setA(e.target.value)} placeholder="a" />
              <Input  className="bg-zinc-600 text-white" value={b} onChange={(e) => setB(e.target.value)} placeholder="b" />
              <Input  className="bg-zinc-600 text-white" value={guess} onChange={(e) => setGuess(e.target.value)} placeholder="x0" />
            </div>
          </div>

          <div className="mb-4">
            <Plot data={[
              { x: xPlot, y: yPlot, type: 'scatter', mode: 'lines', name: 'f(x)' },
              { x: bisectionRes.map(r => r.mid), y: bisectionRes.map(() => 0), type: 'scatter', mode: 'markers', marker: { size: 8 }, name: 'Bisection mid' },
              { x: nrRes.map(r => r.x), y: nrRes.map(r => r.y), type: 'scatter', mode: 'markers+lines', name: 'Newton iter' }
            ]}
              layout={{
                title: 'Function & Root search',
                paper_bgcolor: theme.panelBg,
                plot_bgcolor: theme.panelBg,
                font: { color: theme.text },
                xaxis: { zeroline: true },
                yaxis: { zeroline: true },
                autosize: true,
                margin: { t: 40, r: 20, l: 40, b: 40 }
              }} style={{ width: '100%', height: 320 }} useResizeHandler={true} />
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
            <div className="text-zinc-200">
              <h4 className="font-semibold">Bisection Steps</h4>
              <div className="text-sm max-h-40 overflow-auto mt-2 bg-zinc-800 p-2 rounded">{bisectionRes.length ? (
                bisectionRes.map((r) => <div key={r.iter}>#{r.iter}: mid={r.mid.toFixed(6)}, f(mid)={r.fm.toExponential(2)}</div>)
              ) : <div className="opacity-70">No valid bisection result (check sign change)</div>}</div>
            </div>
            <div className="text-zinc-200">
              <h4 className="font-semibold">Newton Steps</h4>
              <div className="text-sm max-h-40 overflow-auto mt-2 bg-zinc-800 p-2 rounded">{nrRes.length ? (
                nrRes.map((r) => <div key={r.iter}>#{r.iter}: x={r.x.toFixed(6)}, f={r.y.toExponential(2)}</div>)
              ) : <div className="opacity-70">Newton did not converge / derivative missing</div>}</div>
            </div>
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

function SectionInterpolation() {
  const [points, setPoints] = useState([[ -1, -1 ], [0, 0], [1, 0.8]]);
  const [newX, setNewX] = useState(2);
  const [newY, setNewY] = useState(0);

  const addPoint = () => { setPoints([...points, [Number(newX), Number(newY)]]); setNewX(0); setNewY(0); };

  // Lagrange interpolation
  function lagrange(x) {
    const n = points.length;
    let s = 0;
    for (let i = 0; i < n; i++) {
      let xi = points[i][0], yi = points[i][1];
      let Li = 1;
      for (let j = 0; j < n; j++) if (j !== i) Li *= (x - points[j][0]) / (xi - points[j][0]);
      s += yi * Li;
    }
    return s;
  }

  const xs = useMemo(() => {
    const min = Math.min(...points.map(p => p[0])) - 1;
    const max = Math.max(...points.map(p => p[0])) + 1;
    const arr = [];
    for (let x = min; x <= max; x += (max - min) / 200) arr.push(x);
    return arr;
  }, [points]);

  const ys = xs.map(x => lagrange(x));

  return (
    <motion.div variants={cardMotion} initial="initial" animate="enter" whileHover="hover">
      <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700`}>
        <CardHeader>
          <CardTitle className="text-emerald-300">Interpolation — Lagrange</CardTitle>
        </CardHeader>
        <CardContent>
          <p className="text-zinc-200 mb-2">Add data points and see the interpolating polynomial. Works best with small number of points.</p>
          <div className="flex gap-2 mb-3">
            <Input  className="bg-zinc-600 text-white" value={newX} onChange={(e) => setNewX(e.target.value)} placeholder="x" />
            <Input  className="bg-zinc-600 text-white" value={newY} onChange={(e) => setNewY(e.target.value)} placeholder="y" />
            <Button onClick={addPoint} className="bg-emerald-500 cursor-pointer hover:bg-emerald-600">Add Point</Button>
          </div>

          <div className="mb-3">
            <Plot data={[
              { x: xs, y: ys, type: 'scatter', mode: 'lines', name: 'interp' },
              { x: points.map(p => p[0]), y: points.map(p => p[1]), type: 'scatter', mode: 'markers', name: 'data', marker: { size: 8 } }
            ]}
              layout={{ paper_bgcolor: theme.panelBg, plot_bgcolor: theme.panelBg, font: { color: theme.text }, title: 'Lagrange Interpolation', autosize: true, margin: { t: 40, b: 40 } }} useResizeHandler={true} style={{ width: '100%', height: 340 }} />
          </div>

          <div className="text-zinc-200 text-sm">
            <strong>Points:</strong>
            <div className="mt-2 grid grid-cols-2 sm:grid-cols-4 gap-2">
              {points.map((p, i) => (
                <div key={i} className="bg-zinc-800 p-2 rounded">({p[0]}, {p[1]})</div>
              ))}
            </div>
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

function SectionIntegrationDifferentiation() {
  const [fnText, setFnText] = useState('Math.exp(-x*x)');
  const [a, setA] = useState(0);
  const [b, setB] = useState(1);
  const [n, setN] = useState(100);
  const [xEval, setXEval] = useState(0.5);

  const f = useMemo(() => parseFn(fnText), [fnText]);

  const trap = useMemo(() => { if (!f) return null; return trapezoidal((x)=>f(x), Number(a), Number(b), Number(n)); }, [f, a, b, n]);
  const simp = useMemo(() => { if (!f) return null; return simpson((x)=>f(x), Number(a), Number(b), Number(n)); }, [f, a, b, n]);
  const derivative = useMemo(() => { if (!f) return null; return finiteDiff((x)=>f(x), Number(xEval)); }, [f, xEval]);

  const xs = useMemo(() => {
    const arr = [];
    for (let x = Number(a); x <= Number(b); x += (Number(b)-Number(a))/200) arr.push(x);
    return arr;
  }, [a, b]);
  const ys = xs.map(x => f ? f(x) : null);

  return (
    <motion.div variants={cardMotion} initial="initial" animate="enter" whileHover="hover">
      <Card className={`p-4 bg-zinc-900/50 ${theme.panelBg} border border-zinc-700`}> 
        <CardHeader>
          <CardTitle className="text-amber-300">Numerical Integration & Differentiation</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-3 mb-3">
            <Input value={fnText} onChange={(e) => setFnText(e.target.value)} className="col-span-2 bg-zinc-600 text-white" />
            <div className="flex gap-2">
              <Input  className="bg-zinc-600 text-white" value={a} onChange={(e) => setA(e.target.value)} placeholder="a" />
              <Input  className="bg-zinc-600 text-white" value={b} onChange={(e) => setB(e.target.value)} placeholder="b" />
              <Input  className="bg-zinc-600 text-white" value={n} onChange={(e) => setN(e.target.value)} placeholder="n" />
            </div>
          </div>

          <div className="mb-3">
            <Plot data={[{ x: xs, y: ys, type: 'scatter', mode: 'lines', name: 'f(x)' }]} layout={{ title: 'f(x)', paper_bgcolor: theme.panelBg, plot_bgcolor: theme.panelBg, font: { color: theme.text }, autosize: true }} style={{ width: '100%', height: 300 }} useResizeHandler={true} />
          </div>

          <div className="grid grid-cols-1 md:grid-cols-3 gap-3 text-zinc-200">
            <div className="bg-zinc-800 p-3 rounded">Trapezoidal ≈ <strong>{trap !== null ? trap.toFixed(6) : '—'}</strong></div>
            <div className="bg-zinc-800 p-3 rounded">Simpson ≈ <strong>{simp !== null ? simp.toFixed(6) : '—'}</strong></div>
            <div className="bg-zinc-800 p-3 rounded">f'({xEval}) ≈ <strong>{derivative !== null ? derivative.toFixed(6) : '—'}</strong></div>
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

function SectionODE() {
  const [dt, setDt] = useState(0.1);
  const [tMax, setTMax] = useState(6);
  const [y0, setY0] = useState(1);
  const [fnText, setFnText] = useState('y - t*t + 1');

  // f has signature f(y,t) so we create a parser wrapper
  const fWrapped = useMemo(() => {
    const inner = parseFn(fnText);
    if (!inner) return null;
    return (y, t) => {
      try { return inner(y, t); } catch { return 0; }
    };
  }, [fnText]);

  const euler = useMemo(() => fWrapped ? explicitEuler(Number(dt), Number(tMax), fWrapped, Number(y0)) : [], [fWrapped, dt, tMax, y0]);
  const rk4 = useMemo(() => fWrapped ? rungeKutta4(Number(dt), Number(tMax), fWrapped, Number(y0)) : [], [fWrapped, dt, tMax, y0]);

  return (
    <motion.div variants={cardMotion} initial="initial" animate="enter" whileHover="hover">
      <Card className={`p-4 bg-zinc-900/50 ${theme.panelBg} border border-zinc-700`}> 
        <CardHeader>
          <CardTitle className="text-blue-300">ODE Solvers — Euler vs RK4</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="grid grid-cols-1 md:grid-cols-4 gap-2 mb-3">
            <Input value={fnText} onChange={(e) => setFnText(e.target.value)} className="col-span-2 bg-zinc-600 text-white" />
            <Input className="bg-zinc-600 text-white" value={dt} onChange={(e) => setDt(e.target.value)} placeholder="dt" />
            <Input  className="bg-zinc-600 text-white" value={tMax} onChange={(e) => setTMax(e.target.value)} placeholder="tMax" />
            <Input  className="bg-zinc-600 text-white" value={y0} onChange={(e) => setY0(e.target.value)} placeholder="y0" />
          </div>

          <Plot data={[
            { x: euler.map(p => p.t), y: euler.map(p => p.y), type: 'scatter', mode: 'lines+markers', name: 'Euler' },
            { x: rk4.map(p => p.t), y: rk4.map(p => p.y), type: 'scatter', mode: 'lines+markers', name: 'RK4' }
          ]}
            layout={{ title: 'ODE Solution', paper_bgcolor: theme.panelBg, plot_bgcolor: theme.panelBg, font: { color: theme.text }, autosize: true }} useResizeHandler={true} style={{ width: '100%', height: 340 }} />

        </CardContent>
      </Card>
    </motion.div>
  );
}

function SectionLinearSystems() {
  const [Atext, setAtext] = useState('[[2,1],[1,3]]');
  const [btext, setBtext] = useState('[3,7]');
  const [solution, setSolution] = useState(null);

  useEffect(() => {
    try {
      const A = JSON.parse(Atext);
      const b = JSON.parse(btext);
      const sol = gaussElimination(A, b);
      setSolution(sol);
    } catch (e) { setSolution(null); }
  }, [Atext, btext]);

  return (
    <motion.div variants={cardMotion} initial="initial" animate="enter" whileHover="hover">
      <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50  border border-zinc-700`}> 
        <CardHeader>
          <CardTitle className="text-emerald-300">Linear Systems — Gaussian Elimination</CardTitle>
        </CardHeader>
        <CardContent>
          <p className="text-zinc-200 mb-2">Enter matrix <code className="font-mono">A</code> and vector <code className="font-mono">b</code> as JSON arrays.</p>
          <div className="grid grid-cols-1 md:grid-cols-2 gap-2 mb-3">
            <Input className="bg-zinc-500 border border-zinc-700" value={Atext} onChange={(e) => setAtext(e.target.value)} />
            <Input className="bg-zinc-500 border border-zinc-700" value={btext} onChange={(e) => setBtext(e.target.value)} />
          </div>

          <div className="bg-zinc-800 p-3 rounded text-zinc-200">Solution: {solution ? '[' + solution.map(v => v.toFixed(6)).join(', ') + ']' : 'Invalid input'}</div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

function Section3DPlotInteractive() {
  // create a 3D surface (z = sin(sqrt(x^2+y^2))/r) to demonstrate responsive 3D plotting
  const [range, setRange] = useState(4);
  const grid = useMemo(() => {
    const n = 40;
    const xs = [], ys = [], zs = [];
    for (let i = 0; i < n; i++) {
      xs.push(-range + (2 * range * i) / (n - 1));
      ys.push(-range + (2 * range * i) / (n - 1));
    }
    for (let i = 0; i < n; i++) {
      const row = [];
      for (let j = 0; j < n; j++) {
        const x = xs[j], y = ys[i];
        const r = Math.sqrt(x * x + y * y) + 1e-6;
        row.push(Math.sin(r) / r);
      }
      zs.push(row);
    }
    return { xs, ys, zs };
  }, [range]);

  return (
    <motion.div variants={cardMotion} initial="initial" animate="enter" whileHover="hover">
      <Card className={`p-4 ${theme.panelBg}  bg-zinc-900/50 border border-zinc-700`}> 
        <CardHeader>
          <CardTitle className="text-blue-300">Interactive 3D Plot — Responsive</CardTitle>
        </CardHeader>
        <CardContent>
          <p className="text-zinc-200 mb-2">This 3D surface is responsive across devices. Use the range slider to zoom.</p>
          <div className="mb-3">
            <Plot data={[{
              z: grid.zs,
              x: grid.xs,
              y: grid.ys,
              type: 'surface',
              contours: { z: { show: true } },
              colorscale: 'Viridis'
            }]}
              layout={{
                title: '3D Surface (responsive)',
                autosize: true,
                scene: { camera: { eye: { x: 1.2, y: 1.2, z: 0.8 } }, xaxis: { title: 'x' }, yaxis: { title: 'y' }, zaxis: { title: 'z' } },
                paper_bgcolor: theme.panelBg,
                font: { color: theme.text },
                margin: { t: 40, b: 40 }
              }} useResizeHandler={true} style={{ width: '100%', height: 420 }} />
          </div>

          <div className="flex gap-2 items-center">
            <Input type="range" min={2} max={8} value={range} onChange={(e) => setRange(Number(e.target.value))} />
            <div className="text-zinc-200">range: {range}</div>
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
}

const notes = {
  "Root Finding": `
  • Bisection method halves the interval each step until root found.
  • Requires f(a) and f(b) with opposite signs.
  • Newton-Raphson uses tangent slope, faster but depends on derivative and good initial guess.`,

  "Interpolation": `
  • Lagrange interpolation constructs polynomial through all given points.
  • Avoids solving linear systems explicitly.
  • Works well for small datasets but oscillates for large points (Runge's phenomenon).`,

  "Integration & Differentiation": `
  • Trapezoidal rule approximates area under curve with trapezoids. Error ~ O(1/n^2).
  • Simpson’s rule uses parabolas, more accurate. Error ~ O(1/n^4).
  • Numerical differentiation often uses central difference for better accuracy.`,

  "ODE Solvers": `
  • Euler method: simple, but error accumulates quickly.
  • RK4 (Runge-Kutta 4th order): much more accurate for same step size.
  • Both solve y' = f(y,t) with initial condition y(0) = y0.`,

  "Linear Systems": `
  • Gaussian elimination reduces system Ax = b into upper triangular form.
  • Then solve by back substitution.
  • Pivoting improves numerical stability.`,

  "3D Plot": `
  • Interactive 3D plots help visualize surfaces like z = sin(r)/r.
  • Use rotation, zoom, and scaling to understand behavior across dimensions.`
};

function DocumentationPanel() {
  const [openTopic, setOpenTopic] = useState(null);

  return (
    <div className="w-full mt-10 bg-zinc-900 border border-zinc-700 rounded-xl shadow-lg">
      <div className="p-3 border-b border-zinc-700 text-zinc-100 font-semibold text-lg flex items-center gap-2">
        <BookOpen size={18} /> Documentation
      </div>
      <div className="max-h-[400px] overflow-y-auto">
        {Object.keys(notes).map((topic) => (
          <div key={topic} className="border-b border-zinc-800">
            <button
              onClick={() => setOpenTopic(openTopic === topic ? null : topic)}
              className="w-full text-left px-4 py-3 hover:bg-zinc-800/60 flex justify-between items-center"
            >
              <span className="text-zinc-100">{topic}</span>
              {openTopic === topic ? (
                <ChevronDown className="text-zinc-400" size={18} />
              ) : (
                <ChevronRight className="text-zinc-400" size={18} />
              )}
            </button>
            {openTopic === topic && (
              <div className="px-5 py-3 text-zinc-300 whitespace-pre-line">
                {notes[topic]}
              </div>
            )}
          </div>
        ))}
      </div>
    </div>
  );
}
// ---------------- Main Component ----------------
export default function Chapter2() {
  return (
    <div className={`min-h-screen px-6 py-8 bg-zinc-950 text-zinc-100`}> 
      <motion.header initial={{ opacity: 0, y: -10 }} animate={{ opacity: 1, y: 0 }} transition={{ duration: 0.5 }} className="max-w-7xl mx-auto mb-6">
        <h1 className="text-3xl font-bold text-emerald-400">Numerical Methods</h1>
        <p className="text-zinc-300 mt-1">Root finding, interpolation, integration, ODE solvers, linear systems — interactive and responsive for all devices.</p>
      </motion.header>

      <main className="max-w-7xl mx-auto grid grid-cols-1 gap-6">
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          <SectionRootFinding />
          <SectionInterpolation />
        </div>

        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          <SectionIntegrationDifferentiation />
          <SectionODE />
        </div>

        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          
          <div className='grid grid-cols-1 gap-2'>
          <SectionLinearSystems />
           <DocumentationPanel />
           </div>
          <Section3DPlotInteractive />
        </div>

       
    <BottomBar/>
      </main>
    </div>
  );
}
