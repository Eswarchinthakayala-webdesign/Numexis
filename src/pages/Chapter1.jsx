// src/pages/Chapter1.jsx
import React, { useMemo, useState, useRef, useEffect, Suspense } from "react";
import { motion } from "framer-motion";
import Plot from "react-plotly.js";
import {
  LineChart,
  Line,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  Legend,
  ResponsiveContainer,
} from "recharts";
import { Canvas, useFrame } from "@react-three/fiber";
import { OrbitControls, Html, Stats } from "@react-three/drei";
import * as THREE from "three";

import { Card, CardHeader, CardTitle, CardContent } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Slider } from "@/components/ui/slider";
import { Tooltip as UITooltip } from "@/components/ui/tooltip";

import "katex/dist/katex.min.css";
import { BlockMath, InlineMath } from "react-katex";
import BottomBar from "../components/BottomBar";

/* ---------- Physics Helpers ---------- */
const G = 9.81;

function analyticVelocity(t, m, c) {
  if (c === 0) return G * t;
  return (G * m) / c * (1 - Math.exp(-(c / m) * t));
}

function explicitEulerSeries(tMax, dt, m, c) {
  const steps = Math.max(1, Math.ceil(tMax / dt));
  const out = [];
  let y = 0;
  let t = 0;
  for (let i = 0; i <= steps; i++) {
    out.push({ t: +t.toFixed(6), y });
    const dydt = G - (c / m) * y;
    y = y + dt * dydt;
    t += dt;
  }
  return out;
}

/* ---------- Three.js Animated Parachutist ---------- */
function AnimatedParachutist({ analyticData, loopSeconds = 8 }) {
  const ref = useRef();
  const start = useRef(performance.now());

  useFrame(() => {
    if (!analyticData || analyticData.length === 0) return;
    const elapsed = ((performance.now() - start.current) / 1000) % loopSeconds;
    const tmax = analyticData[analyticData.length - 1].t || loopSeconds;
    const t = (elapsed / loopSeconds) * tmax;
    let i = analyticData.findIndex((d) => d.t >= t);
    if (i === -1) i = analyticData.length - 1;
    const vel = analyticData[Math.max(0, i)].analytic || 0;
    if (ref.current) {
      ref.current.position.y = THREE.MathUtils.lerp(ref.current.position.y, 0.2 + vel * 0.05, 0.12);
      ref.current.position.x = THREE.MathUtils.lerp(ref.current.position.x, -1.5 + (t / tmax) * 3, 0.08);
      ref.current.rotation.y += 0.01;
    }
  });

  return (
    <group>
      <mesh ref={ref} castShadow>
        <sphereGeometry args={[0.12, 28, 28]} />
        <meshStandardMaterial color="#00ff00" metalness={0.3} roughness={0.4} />
      </mesh>
      <mesh position={[0, -0.9, 0]} rotation-x={-Math.PI / 2}>
        <cylinderGeometry args={[3, 3, 0.02, 64]} />
        <meshStandardMaterial color="#585858" />
      </mesh>
    </group>
  );
}

/* ---------- Framer Motion Config ---------- */
const cardMotion = {
  initial: { opacity: 0, y: 12, scale: 0.995 },
  enter: { opacity: 1, y: 0, scale: 1, transition: { duration: 0.45, ease: [0.2, 0.8, 0.2, 1] } },
};

/* ---------- Newton-Raphson Helper ---------- */
function newtonSteps(x0 = 1.5, tol = 1e-8, maxIter = 20) {
  let x = x0;
  const steps = [];
  for (let i = 0; i < maxIter; i++) {
    const fx = x * x * x - x - 2;
    const dfx = 3 * x * x - 1;
    const x1 = x - fx / dfx;
    steps.push({ iter: i + 1, x, fx, x1 });
    if (Math.abs(x1 - x) < tol) break;
    x = x1;
  }
  return steps;
}

/* ---------- Main Component ---------- */
export default function Chapter1() {
  const [mass, setMass] = useState(68.1);
  const [drag, setDrag] = useState(12.5);
  const [tMax, setTMax] = useState(12);
  const [dt, setDt] = useState(1);
  const [showThreeHelpers, setShowThreeHelpers] = useState(false);

  const [nrSteps, setNrSteps] = useState([]);

  const analytic = useMemo(() => {
    const N = 200;
    const arr = [];
    for (let i = 0; i <= N; i++) {
      const t = (tMax * i) / N;
      arr.push({ t: +t.toFixed(6), analytic: analyticVelocity(t, mass, drag) });
    }
    return arr;
  }, [mass, drag, tMax]);

  const euler = useMemo(() => explicitEulerSeries(tMax, dt, mass, drag), [mass, drag, tMax, dt]);

  const merged = useMemo(() => {
    if (!analytic || analytic.length === 0) return [];
    return analytic.map((a, idx) => {
      const pos = Math.round((a.t / tMax) * (euler.length - 1));
      const e = euler[Math.min(euler.length - 1, Math.max(0, pos))];
      const error = Math.abs(a.analytic - (e?.y ?? 0));
      return { t: a.t, analytic: +a.analytic.toFixed(6), euler: e?.y ? +e.y.toFixed(6) : null, error: +error.toFixed(6) };
    });
  }, [analytic, euler, tMax]);

  const maxError = useMemo(() => (merged.length ? merged.reduce((acc, r) => Math.max(acc, r.error), 0) : 0), [merged]);

  const sampleTimes = [0, 2, 4, 6, 8, 10, 12];
  const analyticSamples = useMemo(() => sampleTimes.map((t) => ({ t, v: analyticVelocity(t, mass, drag) })), [mass, drag]);
  const eulerSamples = useMemo(() => sampleTimes.map((t) => {
    const idx = Math.min(euler.length - 1, Math.round(t / dt));
    return { t, v: +(euler[idx]?.y ?? 0) };
  }), [euler, dt]);

  const surfaceData = useMemo(() => {
    const N = 60;
    const xs = new Array(N).fill(0).map((_, i) => -6 + (12 * i) / (N - 1));
    const ys = xs;
    const z = xs.map((x) => ys.map((y) => {
      const r = Math.hypot(x, y) || 1e-6;
      return Math.sin(r) / r;
    }));
    return { xs, ys, z };
  }, []);

  function handleRunNewton() { setNrSteps(newtonSteps(1.5, 1e-10, 15)); }

  const prefersReducedMotion = typeof window !== "undefined" && window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  return (
    <div className="min-h-screen bg-zinc-950 text-[#DDE6ED] px-6 py-8">
      <motion.header
        initial={{ opacity: 0, y: -8 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.5 }}
        className="max-w-7xl mx-auto mb-6 flex flex-col lg:flex-row items-start lg:items-center justify-between gap-4"
      >
        <div>
          <h1 className="text-3xl font-semibold text-emerald-400">Mathematical Modeling & Engineering Problem Solving</h1>
          <p className="text-sm text-[#A0AEC0] mt-1">Interactive lesson: parachutist example — analytic vs explicit Euler, plus visual intuition.</p>
        </div>

        <div className="flex items-center gap-3 flex-wrap">
          <Button variant="ghost" className="hover:bg-white text-black cursor-pointer bg-gray-400" onClick={() => { setDt(Math.max(0.02, Math.round((dt / 2) * 100) / 100)); }}>Halve dt</Button>
          <Button className="hover:bg-white text-black cursor-pointer bg-gray-400" onClick={() => { setDt(Math.min(2, Math.round((dt * 2) * 100) / 100)); }}>Double dt</Button>
          <Button variant="outline" className="bg-emerald-400 text-black border border-zinc-500 cursor-pointer hover:bg-emerald-300" onClick={() => { setMass(68.1); setDrag(12.5); setDt(1); setTMax(12); }}>Reset</Button>
        </div>
      </motion.header>

      <main className="max-w-7xl mx-auto grid grid-cols-1 lg:grid-cols-3 gap-6">
        {/* Left Column: theory + params */}
        <motion.div {...cardMotion} className="col-span-1" initial="initial" animate="enter">
          <Card className="bg-zinc-900/50 border border-zinc-500">
            <CardHeader>
              <CardTitle className="text-white">Model & Theory</CardTitle>
            </CardHeader>
            <CardContent>
              <p className="text-sm text-[#DDE6ED] mb-3">
                We model the falling object (parachutist) with linear drag:
              </p>

              <div className="mb-3 text-[#DDE6ED]">
                <BlockMath math={`\\frac{dy}{dt} = g - \\frac{c}{m} y`} />
                <div className="mt-2 text-xs text-[#A0AEC0]">With <InlineMath math="y(0)=0" /> the analytic solution is:</div>
                <BlockMath math={`y(t)=\\frac{gm}{c}\\left(1-e^{-\\frac{c}{m}t}\\right)`} />
                <div className="mt-2 text-xs text-[#A0AEC0]">Compare analytic vs explicit Euler numerical integration.</div>
              </div>

              <div className="space-y-3">
                <label className="block text-xs text-[#A0AEC0]">Mass (kg)</label>
                <Input type="number" value={mass} onChange={(e) => setMass(Number(e.target.value))} className="bg-[#0D1B1E] text-[#DDE6ED] border-zinc-500" />

                <label className="block text-xs text-[#A0AEC0] mt-2">Drag coefficient c (kg/s)</label>
                <Input type="number" value={drag} onChange={(e) => setDrag(Number(e.target.value))} className="bg-[#0D1B1E] text-[#DDE6ED] border-zinc-500" />

                <div className="mt-2">
                  <label className="block text-xs text-[#A0AEC0]">Simulation time t<sub>max</sub> (s)</label>
                  <Input type="number" value={tMax} onChange={(e) => setTMax(Number(e.target.value))} className="bg-[#0D1B1E] text-[#DDE6ED] border-zinc-500" />
                </div>

                <div className="mt-2">
                  <label className="block text-xs text-[#A0AEC0]">Euler step dt (s)</label>
                  <div className="flex items-center gap-3">
                    <input
                      type="range"
                      min="0.02"
                      max="2"
                      step="0.02"
                      value={dt}
                      onChange={(e) => setDt(Number(e.target.value))}
                      className="w-full accent-[#717473]"
                    />
                    <div className="w-16 text-right text-xs text-[#DDE6ED]">{dt}s</div>
                  </div>
                </div>

                <div className="mt-3 flex gap-2 flex-wrap">
                  <Button className="bg-zinc-600 cursor-pointer hover:bg-zinc-400" onClick={() => { setDt(1); setTMax(12); }}>Example (dt=1)</Button>
                  <Button className="bg-white text-black cursor-pointer hover:bg-gray-300" onClick={() => setShowThreeHelpers((v) => !v)}>{showThreeHelpers ? "Hide Helpers" : "Show Helpers"}</Button>
                </div>
              </div>
            </CardContent>
          </Card>

          {/* Sample Tables */}
          <motion.div className="mt-4" initial={{ opacity: 0 }} animate={{ opacity: 1 }} transition={{ delay: 0.12 }}>
            <Card className="bg-[#0B1B20]/50 border border-zinc-500">
              <CardHeader>
                <CardTitle className="text-white">Example Samples</CardTitle>
              </CardHeader>
              <CardContent>
                <div className="text-xs text-[#A0AEC0] mb-2">Analytic (t = 0,2,...,12 s)</div>
                <div className="space-y-1 text-sm text-[#DDE6ED]">
                  {analyticSamples.map((s) => (
                    <div key={s.t} className="flex justify-between">
                      <div>t = {s.t}s</div>
                      <div>{s.v.toFixed(4)} m/s</div>
                    </div>
                  ))}
                </div>

                <div className="text-xs text-[#A0AEC0] mt-3 mb-2">Euler (dt = {dt}s)</div>
                <div className="space-y-1 text-sm text-[#DDE6ED]">
                  {eulerSamples.map((s) => (
                    <div key={s.t} className="flex justify-between">
                      <div>t = {s.t}s</div>
                      <div>{s.v.toFixed(4)} m/s</div>
                    </div>
                  ))}
                </div>

                <div className="mt-3 text-xs text-[#A0AEC0]">Max abs error: <span className="text-[#E67E22] ml-2 font-medium">{maxError.toFixed(4)} m/s</span></div>
              </CardContent>
            </Card>
          </motion.div>
        </motion.div>

        {/* Middle/Right Visuals */}
        <motion.div className="col-span-1 lg:col-span-2 grid grid-cols-1 gap-6" initial="initial" animate="enter">

          {/* Recharts: Velocity vs Time */}
          <Card className="bg-[#0B1B20]/60 border border-zinc-500">
            <CardHeader>
              <CardTitle className="text-white">Velocity vs Time — Analytic vs Euler</CardTitle>
            </CardHeader>
            <CardContent className="h-[340px]">
              <ResponsiveContainer width="100%" height="100%">
                <LineChart data={merged} margin={{ top: 12, right: 18, left: 6, bottom: 6 }}>
                  <CartesianGrid stroke="#1f2937" />
                  <XAxis dataKey="t" tick={{ fill: "#DDE6ED" }} label={{ value: "t (s)", position: "insideBottom", offset: -6, fill: "#DDE6ED" }} />
                  <YAxis tick={{ fill: "#DDE6ED" }} label={{ value: "v (m/s)", angle: -90, position: "insideLeft", fill: "#DDE6ED" }} />
                  <Tooltip contentStyle={{ background: "#0B1B20", border: "1px solid #16A085" }} />
                  <Legend wrapperStyle={{ color: "#DDE6ED" }} />
                  <Line name="Analytic" dataKey="analytic" stroke="#1ABC9C" strokeWidth={2} dot={false} />
                  <Line name="Euler" dataKey="euler" stroke="#E67E22" strokeWidth={2} dot={false} strokeDasharray="6 4" />
                </LineChart>
              </ResponsiveContainer>
            </CardContent>
          </Card>

          {/* Grid: Error + Plotly 3D */}
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">

            {/* Absolute Error */}
            <Card className="bg-[#0B1B20]/60 border border-zinc-500">
              <CardHeader>
                <CardTitle className="text-white">Absolute Error over Time</CardTitle>
              </CardHeader>
              <CardContent className="h-[260px]">
                <ResponsiveContainer width="100%" height="100%">
                  <LineChart data={merged}>
                    <CartesianGrid stroke="#1f2937" />
                    <XAxis dataKey="t" tick={{ fill: "#DDE6ED" }} />
                    <YAxis tick={{ fill: "#DDE6ED" }} />
                    <Tooltip contentStyle={{ background: "#0B1B20", border: "1px solid #16A085" }} />
                    <Line type="monotone" dataKey="error" stroke="#EF4444" strokeWidth={2} dot={false} />
                  </LineChart>
                </ResponsiveContainer>
              </CardContent>
            </Card>

            {/* Plotly 3D Surface */}
            <Card className="bg-[#0B1B20]/60 border border-zinc-500">
              <CardHeader>
                <CardTitle className="text-white">3D Surface — Illustration</CardTitle>
              </CardHeader>
              <CardContent className="h-[260px]">
                <Plot
                  data={[{ z: surfaceData.z, x: surfaceData.xs, y: surfaceData.ys, type: "surface", colorscale: "Viridis", showscale: false }]}
                  layout={{
                    autosize: true,
                    margin: { t: 10, b: 30, l: 30, r: 10 },
                    paper_bgcolor: "#0D1B1E",
                    plot_bgcolor: "#0D1B1E",
                    scene: { xaxis: { title: "x", color: "#DDE6ED" }, yaxis: { title: "y", color: "#DDE6ED" }, zaxis: { title: "z", color: "#DDE6ED" } },
                  }}
                  style={{ width: "100%", height: "100%" }}
                  config={{ responsive: true, displayModeBar: false }}
                />
              </CardContent>
            </Card>

          </div>

          {/* Full-width Three.js */}
          <Card className="bg-[#0B1B20]/60 border border-zinc-500">
            <CardHeader>
              <CardTitle className="text-white">3D Scene — Animated Parachutist</CardTitle>
            </CardHeader>
            <CardContent className="h-[420px] p-0 overflow-hidden">
              <div className="h-[420px] w-full bg-gray-600">
                <Canvas camera={{ position: [0, 2.2, 6], fov: 50 }}>
                  <ambientLight intensity={0.6} />
                  <directionalLight position={[5, 10, 5]} intensity={1.0} castShadow />
                  <Suspense fallback={<Html center><div className="text-[#DDE6ED]">Loading 3D...</div></Html>}>
                    <AnimatedParachutist analyticData={analytic} loopSeconds={Math.max(6, tMax)} />
                  </Suspense>
                  <OrbitControls enablePan enableZoom enableRotate />
                  {showThreeHelpers && (
                    <>
                      <gridHelper args={[10, 10, "#1ABC9C", "#0D1B1E"]} />
                      <axesHelper args={[3]} />
                      <Stats />
                    </>
                  )}
                </Canvas>
              </div>
            </CardContent>
          </Card>

          {/* Newton + Code Snippets */}
          <Card className="bg-[#0B1B20]/60 border border-zinc-500">
            <CardHeader>
              <CardTitle className="text-white">Implementation & Exercises</CardTitle>
            </CardHeader>
            <CardContent>
              <div className="grid grid-cols-1 md:grid-cols-2 gap-6">

                <div>
                  <h4 className="text-sm font-semibold text-[#DDE6ED]">Newton–Raphson (tiny demo)</h4>
                  <div className="text-xs text-[#A0AEC0] mb-2">Solve f(x) = x^3 - x - 2</div>
                  <Button className="bg-gray-500 cursor-pointer hover:bg-white text-black" onClick={handleRunNewton}>Run Newton–Raphson</Button>

                  {nrSteps.length > 0 && (
                    <div className="mt-3 text-xs space-y-1">
                      {nrSteps.map((s) => ( 
                        <div key={s.iter} className="flex justify-between text-white bg-[#0B1B20]/30 p-2 rounded">
                          <div>Iter {s.iter}</div>
                          <div>x ≈ {s.x1.toFixed(8)}</div>
                        </div>
                      ))}
                    </div>
                  )}
                </div>

                <div>
                  <h4 className="text-sm font-semibold text-[#DDE6ED]">Code snippets</h4>
                  <div className="mt-2">
                    <div className="text-xs font-medium text-[#A0AEC0]">Explicit Euler (JS)</div>
                    <pre className="mt-2 bg-[#0D1B1E]/60 p-3 rounded text-xs overflow-auto text-[#DDE6ED]">{`// y_{n+1} = y_n + dt * [g - (c/m) * y_n]
function explicitEuler(dt, tMax, m, c) {
  const steps = [];
  let y=0, t=0;
  while (t<=tMax+1e-9) {
    steps.push({t, y});
    y = y + dt * (g - (c/m)*y);
    t += dt;
  }
  return steps;
}`}</pre>
                  </div>
                </div>

              </div>

              <div className="mt-4 text-xs text-[#A0AEC0]">
                <strong>Exercises</strong>
                <ol className="list-decimal list-inside mt-2 space-y-1">
                  <li>Implement explicit Euler for the parachutist and reproduce Example 1.2 tables with dt=2s and dt=1s.</li>
                  <li>Plot analytic vs numeric for dt = 2, 1, 0.5 and discuss how dt affects accuracy and computational cost.</li>
                  <li>Use Newton–Raphson to solve a sample nonlinear equation and plot convergence vs iteration.</li>
                </ol>
              </div>
            </CardContent>
          </Card>

        </motion.div>
        <BottomBar/>
      </main>
    </div>
  );
}
