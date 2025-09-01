

import React, {
  useEffect,
  useMemo,
  useRef,
  useState,
  Suspense,
  useCallback,
} from "react";
import { Canvas, useFrame } from "@react-three/fiber";
import { OrbitControls, Html, Preload } from "@react-three/drei";
import { motion, AnimatePresence } from "framer-motion";
import {
  Rocket,
  Zap,
  Menu,
  Globe,
  Boxes,
  Users,
  GitBranch,
  ArrowRight,
  ChevronDown,
  Github,
  Mail,
  Sparkles,
  Play,
  Code,
  BookOpen,
  HelpCircle,
  Calculator,
  X,
  Home,
  BarChart3,
  FlaskConical,
  CreditCard,
} from "lucide-react";

/* adjust shadcn/ui imports to match your project paths */
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Input } from "@/components/ui/input";
import { useNavigate } from "react-router-dom";

/* ============================
   THEME COLORS & HELPERS
   ============================ */
const THEME = {
  bg: "#0D1B1E",
  text: "#DDE6ED",
  accent: "#1ABC9C",
  accent2: "#16A085",
  warn: "#E67E22",
  glass: "rgba(13,27,30,0.6)",
};

const CONTAINER = "mx-auto max-w-[1380px] px-6 lg:px-12";
const DARK_GLASS =
  "backdrop-blur-md bg-zinc-900/40 border border-zinc-800 shadow-[0_0_40px_rgba(16,185,129,0.08)]";

/* tiny safe evaluator for user-supplied functions (used in visualizers) */
function makeSafeFn(expr) {
  try {
    const fn = new Function("Math", "x", `return (${expr});`);
    return (x) => {
      try {
        const v = fn(Math, x);
        return typeof v === "number" && isFinite(v) ? v : NaN;
      } catch {
        return NaN;
      }
    };
  } catch {
    return () => NaN;
  }
}

/* small loader component used in Suspense fallbacks */
function Loader({ text = "Rendering…" }) {
  return (
    <div className="flex items-center justify-center h-40 bg-zinc-900/30 rounded">
      <div className="text-zinc-400">{text}</div>
    </div>
  );
}

/* ============================
   3D SCENE PIECES (many)
   - ContourField
   - GradientSwarm
   - RelaxationSurface
   - EigenArrows
   - RootBeacons
   - MatrixGrid (visual for linear systems)
   - ParticleField (background filler)
   - WaveField (animated wave for PDE)
   ============================ */

/* ---------- ContourField (rotating level sets) ---------- */
function ContourField({ rings = 18, a = 1.0, b = 0.35, c = 0.15, scale = 1 }) {
  const group = useRef();
  useFrame(({ clock }) => {
    if (!group.current) return;
    const t = clock.elapsedTime;
    group.current.rotation.z = Math.sin(t * 0.2) * 0.12;
    group.current.rotation.y = Math.cos(t * 0.12) * 0.06;
    const s = 1 + Math.sin(t * 0.18) * 0.02;
    group.current.scale.set(s * scale, s * scale, s * scale);
  });

  const ringData = useMemo(() => {
    const arr = [];
    for (let i = 0; i < rings; i++) arr.push(0.3 + i * 0.35);
    return arr;
  }, [rings]);

  return (
    <group ref={group} position={[0, -1.05, -8]}>
      {ringData.map((v, i) => {
        const points = [];
        const N = 320;
        for (let k = 0; k <= N; k++) {
          const t = (k / N) * Math.PI * 2;
          const u = Math.cos(t),
            w = Math.sin(t);
          const denomMajor = Math.max(0.12, a - c);
          const denomMinor = Math.max(0.12, b - c);
          const rx = Math.sqrt((2 * v) / denomMajor);
          const ry = Math.sqrt((2 * v) / denomMinor);
          const x = rx * u;
          const y = ry * w + (i % 2 === 0 ? 0 : 0.06);
          points.push(x, y, -i * 0.01);
        }
        const array = new Float32Array(points);
        const color = i % 2 ? THEME.accent : THEME.accent2;
        return (
          <line key={i} frustumCulled={false}>
            <bufferGeometry>
              <bufferAttribute
                attach="attributes-position"
                array={array}
                itemSize={3}
                count={array.length / 3}
              />
            </bufferGeometry>
            <lineBasicMaterial transparent opacity={0.72} color={color} linewidth={1} />
          </line>
        );
      })}
    </group>
  );
}

/* ---------- GradientSwarm (particles following -grad f) ---------- */
function GradientSwarm({ n = 220, step = 0.007, radius = 5, seed = 0 }) {
  const group = useRef();
  const state = useMemo(() => {
    const pts = [];
    for (let i = 0; i < n; i++) {
      const r = 1 + Math.random() * radius;
      const t = Math.random() * Math.PI * 2;
      pts.push({ x: Math.cos(t) * r, y: Math.sin(t) * r, z: -3 - Math.random() * 2, v: 0.3 + Math.random() });
    }
    return pts;
  }, [n, radius, seed]);

  const u = 1.05,
    v = 0.6,
    k = 0.22;

  useFrame(({ clock }) => {
    const t = clock.elapsedTime;
    if (!group.current) return;
    group.current.rotation.y = Math.sin(t * 0.05) * 0.08;
    for (const p of state) {
      const gx = u * p.x + k * p.y;
      const gy = v * p.y + k * p.x;
      p.x -= step * gx;
      p.y -= step * gy;
      p.z += Math.sin((p.x + p.y + t) * 0.6) * 0.002;
      if (Math.hypot(p.x, p.y) < 0.06) {
        const r = 3 + Math.random() * 5;
        const th = Math.random() * Math.PI * 2;
        p.x = Math.cos(th) * r;
        p.y = Math.sin(th) * r;
        p.z = -3 - Math.random() * 2;
      }
    }
  });

  return (
    <group ref={group} position={[0, -0.15, -10]}>
      {state.map((p, i) => (
        <mesh key={i} position={[p.x, p.y * 0.7, p.z]}>
          <sphereGeometry args={[0.06 + (i % 4) * 0.01, 12, 12]} />
          <meshStandardMaterial
            color={i % 2 ? THEME.accent : THEME.accent2}
            metalness={0.5}
            roughness={0.3}
            emissive={"#021617"}
            emissiveIntensity={0.25}
          />
        </mesh>
      ))}
    </group>
  );
}

/* ---------- RelaxationSurface (Jacobi-like iterative smoothing) ---------- */
function RelaxationSurface({ W = 80, H = 56, amplitude = 1.0 }) {
  const meshRef = useRef();
  const fieldRef = useRef(null);

  useEffect(() => {
    const f = new Float32Array((W + 1) * (H + 1)).fill(0);
    const poke = (ix, iy, amp) => {
      const idx = iy * (W + 1) + ix;
      f[idx] += amp;
    };
    poke(Math.floor(W * 0.18), Math.floor(H * 0.4), 1.0 * amplitude);
    poke(Math.floor(W * 0.72), Math.floor(H * 0.62), -0.85 * amplitude);
    poke(Math.floor(W * 0.52), Math.floor(H * 0.24), 0.6 * amplitude);
    fieldRef.current = f;
  }, [W, H, amplitude]);

  useFrame(({ clock }) => {
    if (!meshRef.current) return;
    const geom = meshRef.current.geometry;
    const pos = geom.attributes.position;
    if (!pos) return;
    const f = fieldRef.current;
    if (!f) return;

    const next = new Float32Array(f.length);
    const alpha = 0.94;
    for (let y = 1; y < H; y++) {
      for (let x = 1; x < W; x++) {
        const i = y * (W + 1) + x;
        const avg = (f[i - 1] + f[i + 1] + f[i - (W + 1)] + f[i + (W + 1)]) * 0.25;
        next[i] = alpha * avg;
      }
    }

    for (let x = 0; x <= W; x++) {
      next[x] = 0;
      next[H * (W + 1) + x] = 0;
    }
    for (let y = 0; y <= H; y++) {
      next[y * (W + 1)] = 0;
      next[y * (W + 1) + W] = 0;
    }
    fieldRef.current = next;

    let idx = 0;
    for (let y = 0; y <= H; y++) {
      for (let x = 0; x <= W; x++) {
        const z = next[idx] * 0.9;
        pos.setZ(idx, z);
        idx++;
      }
    }
    pos.needsUpdate = true;

    meshRef.current.rotation.x = -0.92;
    meshRef.current.rotation.z = Math.sin(clock.elapsedTime * 0.12) * 0.06;
  });

  return (
    <mesh ref={meshRef} position={[0, -1.28, -7]}>
      <planeGeometry args={[18, 14, W, H]} />
      <meshStandardMaterial
        color={THEME.accent}
        metalness={0.08}
        roughness={0.9}
        wireframe
        transparent
        opacity={0.28}
      />
    </mesh>
  );
}

/* ---------- EigenArrows (power-method visualization) ---------- */
function EigenArrows() {
  const group = useRef();
  const A = [
    [1.0, 0.28],
    [0.12, 0.65],
  ];
  const vRef = useRef([0.3, 1.0]);
  useFrame(({ clock }) => {
    const t = clock.elapsedTime;
    if (!group.current) return;
    let [x, y] = vRef.current;
    const nx = A[0][0] * x + A[0][1] * y;
    const ny = A[1][0] * x + A[1][1] * y;
    const norm = Math.hypot(nx, ny) || 1;
    vRef.current = [nx / norm, ny / norm];
    const angle = Math.atan2(vRef.current[1], vRef.current[0]);
    const len = 2.6 + Math.sin(t * 0.9) * 0.2;
    group.current.rotation.z = angle;
    group.current.scale.set(len, 1, 1);
  });

  return (
    <group position={[0, 0.1, -6.5]}>
      <group ref={group}>
        <mesh position={[1.3, 0, 0]}>
          <boxGeometry args={[2.6, 0.06, 0.06]} />
          <meshStandardMaterial
            color="#93c5fd"
            metalness={0.6}
            roughness={0.2}
            emissive="#06202a"
            emissiveIntensity={0.45}
          />
        </mesh>
        <mesh position={[2.7, 0, 0]} rotation={[0, 0, Math.PI / 2]}>
          <coneGeometry args={[0.18, 0.36, 12]} />
          <meshStandardMaterial color="#6ee7b7" metalness={0.6} roughness={0.2} />
        </mesh>
      </group>
    </group>
  );
}

/* ---------- RootBeacons (pulsing root markers) ---------- */
function Beacon({ position = [0, 0, -7.5], color = THEME.accent }) {
  const ringRef = useRef();
  useFrame(({ clock }) => {
    const t = clock.elapsedTime;
    if (!ringRef.current) return;
    const s = 0.7 + (Math.sin(t * 2.0) + 1) * 0.35;
    ringRef.current.scale.set(s, s, s);
    if (ringRef.current.material) ringRef.current.material.opacity = 0.35 + (Math.sin(t * 2.0) + 1) * 0.2;
  });
  return (
    <group position={position}>
      <mesh>
        <sphereGeometry args={[0.14, 16, 16]} />
        <meshStandardMaterial color={color} metalness={0.82} roughness={0.15} />
      </mesh>
      <mesh ref={ringRef} rotation={[Math.PI / 2, 0, 0]} position={[0, -0.01, 0]}>
        <torusGeometry args={[0.48, 0.016, 8, 80]} />
        <meshBasicMaterial color={color} transparent opacity={0.45} />
      </mesh>
    </group>
  );
}

function RootBeacons({ roots = [[-3.2, -0.4, -7.5], [0.0, 0.2, -7.5], [3.1, -0.1, -7.5]] }) {
  return (
    <group>
      {roots.map((p, i) => (
        <Beacon key={i} position={p} color={i % 2 ? THEME.accent : THEME.accent2} />
      ))}
    </group>
  );
}

/* ---------- ParticleField (background star-like) ---------- */
function ParticleField({ count = 1500 }) {
  const ref = useRef();
  const positions = useMemo(() => {
    const pos = new Float32Array(count * 3);
    for (let i = 0; i < count; i++) {
      pos[i * 3 + 0] = (Math.random() - 0.5) * 40;
      pos[i * 3 + 1] = (Math.random() - 0.5) * 20;
      pos[i * 3 + 2] = -5 - Math.random() * 40;
    }
    return pos;
  }, [count]);

  useFrame(({ clock }) => {
    if (!ref.current) return;
    const t = clock.elapsedTime;
    ref.current.rotation.y = t * 0.02;
  });

  return (
    <group>
      <points ref={ref}>
        <bufferGeometry>
          <bufferAttribute attach="attributes-position" array={positions} itemSize={3} count={positions.length / 3} />
        </bufferGeometry>
        <pointsMaterial size={0.04} sizeAttenuation depthWrite={false} transparent opacity={0.6} color={THEME.accent} />
      </points>
    </group>
  );
}

/* ---------- WaveField (animated sine waves for PDE demos) ---------- */
function WaveField({ cols = 80, rows = 28, amp = 0.6 }) {
  const ref = useRef();
  useFrame(({ clock }) => {
    const t = clock.elapsedTime;
    if (!ref.current) return;
    const pos = ref.current.geometry.attributes.position;
    let idx = 0;
    for (let j = 0; j <= rows; j++) {
      for (let i = 0; i <= cols; i++) {
        const x = (i / cols - 0.5) * 14;
        const y = (j / rows - 0.5) * 8;
        const z = Math.sin((i + t * 6) * 0.4) * Math.cos((j + t * 5) * 0.35) * 0.25 * amp;
        pos.setZ(idx, z);
        idx++;
      }
    }
    pos.needsUpdate = true;
    ref.current.rotation.x = -0.95;
    ref.current.rotation.z = Math.sin(t * 0.13) * 0.04;
  });

  return (
    <mesh position={[0, -1.2, -7]} ref={ref}>
      <planeGeometry args={[16, 10, cols, rows]} />
      <meshStandardMaterial color={THEME.accent2} wireframe transparent opacity={0.32} />
    </mesh>
  );
}

/* ============================
   SVG Plot helper for visualizer previews
   ============================ */
function SVGPlot({ fn, domain = [-6, 6], w = 520, h = 260, points = 420, markers = [] }) {
  const [path, setPath] = useState("");
  const [ymin, ymax] = useMemo(() => {
    let min = Infinity,
      max = -Infinity;
    for (let i = 0; i < points; i++) {
      const x = domain[0] + (i / (points - 1)) * (domain[1] - domain[0]);
      const y = fn(x);
      if (Number.isFinite(y)) {
        min = Math.min(min, y);
        max = Math.max(max, y);
      }
    }
    if (!isFinite(min) || !isFinite(max)) return [-1, 1];
    if (min === max) {
      min -= 1;
      max += 1;
    }
    return [min, max];
  }, [fn, domain, points]);

  useEffect(() => {
    const pts = [];
    for (let i = 0; i < points; i++) {
      const x = domain[0] + (i / (points - 1)) * (domain[1] - domain[0]);
      const y = fn(x);
      const sx = ((x - domain[0]) / (domain[1] - domain[0])) * w;
      const sy = h - ((y - ymin) / (ymax - ymin)) * h;
      if (!Number.isFinite(sy)) continue;
      pts.push(`${sx},${sy}`);
    }
    setPath(pts.join(" "));
  }, [fn, domain, points, w, h, ymin, ymax]);

  return (
    <svg width={w} height={h} viewBox={`0 0 ${w} ${h}`} className="rounded">
      <rect width={w} height={h} fill="transparent" />
      <line x1="0" x2={w} y1={h / 2} y2={h / 2} stroke="#071417" strokeWidth="1" />
      <polyline fill="none" stroke="#60a5fa" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" points={path} />
      {markers.map((m, i) => {
        const mx = ((m.x - domain[0]) / (domain[1] - domain[0])) * w;
        const my = h - ((m.y - ymin) / (ymax - ymin)) * h;
        if (!Number.isFinite(mx) || !Number.isFinite(my)) return null;
        return <circle key={i} cx={mx} cy={my} r={4} fill={m.color || THEME.warn} stroke="#071417" />;
      })}
    </svg>
  );
}

/* ============================
   VISUALIZERS (Newton & Bisection) — Enhanced
   ============================ */
function NewtonVisualizer() {
  const [expr, setExpr] = useState("Math.pow(x,3)-2*x-5");
  const [x0, setX0] = useState(2);
  const [tol, setTol] = useState(1e-6);
  const [maxIter, setMaxIter] = useState(30);
  const [steps, setSteps] = useState([]);
  const fnRef = useRef(() => NaN);

  useEffect(() => {
    fnRef.current = makeSafeFn(expr);
    setSteps([]);
  }, [expr]);

  const deriv = (f, x) => {
    const h = 1e-6;
    return (f(x + h) - f(x - h)) / (2 * h);
  };

  const compute = useCallback(() => {
    const f = fnRef.current;
    const seq = [];
    let x = Number(x0);
    for (let i = 0; i < maxIter; i++) {
      const fx = f(x);
      if (!Number.isFinite(fx)) break;
      seq.push({ x, fx });
      const d = deriv(f, x);
      if (!Number.isFinite(d) || Math.abs(d) < 1e-12) break;
      const nx = x - fx / d;
      if (!Number.isFinite(nx)) break;
      if (Math.abs(nx - x) < tol) {
        seq.push({ x: nx, fx: f(nx) });
        break;
      }
      x = nx;
    }
    setSteps(seq);
  }, [x0, tol, maxIter]);

  const animate = async () => {
    const f = fnRef.current;
    let x = Number(x0);
    const seq = [];
    for (let i = 0; i < maxIter; i++) {
      const fx = f(x);
      if (!Number.isFinite(fx)) break;
      seq.push({ x, fx });
      const d = deriv(f, x);
      if (!Number.isFinite(d) || Math.abs(d) < 1e-12) break;
      const nx = x - fx / d;
      if (!Number.isFinite(nx)) break;
      if (Math.abs(nx - x) < tol) {
        seq.push({ x: nx, fx: f(nx) });
        break;
      }
      x = nx;
    }
    setSteps([]);
    for (let i = 0; i < seq.length; i++) {
      setSteps((s) => [...s, seq[i]]);
      // small delay to animate
      // eslint-disable-next-line no-await-in-loop
      await new Promise((res) => setTimeout(res, 280));
    }
  };

  const markers = steps.map((s, i) => ({
    x: s.x,
    y: s.fx,
    color: i === steps.length - 1 ? THEME.accent2 : THEME.warn,
  }));

  return (
    <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
      <div className="space-y-3">
        <label className="text-sm text-zinc-400">f(x)</label>
        <input
          className="w-full mt-2 px-3 py-2 rounded bg-zinc-900 border border-zinc-800 text-sm"
          value={expr}
          onChange={(e) => setExpr(e.target.value)}
        />
        <div className="flex gap-3 items-end overflow-auto">
          <div>
            <label className="text-xs text-zinc-400">x₀</label>
            <input
              className="mt-1 px-3 py-2 rounded bg-zinc-900 text-sm w-28 border border-zinc-800"
              value={x0}
              onChange={(e) => setX0(Number(e.target.value))}
            />
          </div>
          <div>
            <label className="text-xs text-zinc-400">tol</label>
            <input
              className="mt-1 px-3 py-2 rounded bg-zinc-900 text-sm w-36 border border-zinc-800"
              value={tol}
              onChange={(e) => setTol(Number(e.target.value))}
            />
          </div>
          <div>
            <label className="text-xs text-zinc-400">maxIter</label>
            <input
              className="mt-1 px-3 py-2 rounded bg-zinc-900 text-sm w-24 border border-zinc-800"
              value={maxIter}
              onChange={(e) => setMaxIter(Number(e.target.value))}
            />
          </div>
        </div>

        <div className="flex gap-2">
          <Button onClick={compute}>Compute</Button>
          <Button onClick={animate}>Animate</Button>
          <Button variant="outline" onClick={() => setSteps([])}>
            Reset
          </Button>
        </div>

        <div className="text-sm text-zinc-400">
          Steps: {steps.length}{" "}
          {steps.length > 0 && (
            <span>
              • Last x: {Number(steps[steps.length - 1].x).toFixed(6)}
            </span>
          )}
        </div>
      </div>

      <div style={{ maxHeight: 340, overflow: "auto" }} className="rounded border border-zinc-800 p-2 bg-zinc-900/30">
        <SVGPlot fn={fnRef.current} domain={[-6, 6]} markers={markers} w={520} h={280} />
      </div>
    </div>
  );
}

function BisectionVisualizer() {
  const [expr, setExpr] = useState("Math.cos(x)-x");
  const [a, setA] = useState(0);
  const [b, setB] = useState(1);
  const [tol, setTol] = useState(1e-6);
  const [maxIter, setMaxIter] = useState(50);
  const [steps, setSteps] = useState([]);
  const fnRef = useRef(() => NaN);

  useEffect(() => {
    fnRef.current = makeSafeFn(expr);
    setSteps([]);
  }, [expr]);

  const run = useCallback(() => {
    const f = fnRef.current;
    let fa = f(a),
      fb = f(b);
    if (!isFinite(fa) || !isFinite(fb) || fa * fb > 0) {
      setSteps([]);
      return;
    }
    let lo = a,
      hi = b;
    const seq = [];
    for (let i = 0; i < maxIter; i++) {
      const mid = (lo + hi) / 2;
      const fm = f(mid);
      seq.push({ lo, hi, mid, fm });
      if (!isFinite(fm)) break;
      if (Math.abs(fm) < tol || (hi - lo) / 2 < tol) break;
      if (fa * fm <= 0) {
        hi = mid;
        fb = fm;
      } else {
        lo = mid;
        fa = fm;
      }
    }
    setSteps(seq);
  }, [a, b, tol, maxIter]);

  const animate = async () => {
    const f = fnRef.current;
    let fa = f(a),
      fb = f(b);
    if (!isFinite(fa) || !isFinite(fb) || fa * fb > 0) {
      setSteps([]);
      return;
    }
    let lo = a,
      hi = b;
    const seq = [];
    for (let i = 0; i < maxIter; i++) {
      const mid = (lo + hi) / 2;
      const fm = f(mid);
      seq.push({ lo, hi, mid, fm });
      setSteps([...seq]);
      // eslint-disable-next-line no-await-in-loop
      await new Promise((res) => setTimeout(res, 340));
      if (!isFinite(fm)) break;
      if (Math.abs(fm) < tol || (hi - lo) / 2 < tol) break;
      if (fa * fm <= 0) {
        hi = mid;
        fb = fm;
      } else {
        lo = mid;
        fa = fm;
      }
    }
  };

  const markers = steps.map((s, i) => ({ x: s.mid, y: s.fm, color: i === steps.length - 1 ? THEME.accent2 : THEME.warn }));

  return (
    <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
      <div className="space-y-3">
        <label className="text-sm text-zinc-400">f(x)</label>
        <input
          className="w-full mt-2 px-3 py-2 rounded bg-zinc-900 border border-zinc-800 text-sm"
          value={expr}
          onChange={(e) => setExpr(e.target.value)}
        />
        <div className="flex gap-3 overflow-auto items-end">
          <div>
            <label className="text-xs text-zinc-400">a</label>
            <input
              className="mt-1 px-3 py-2 rounded bg-zinc-900 text-sm w-24 border border-zinc-800"
              value={a}
              onChange={(e) => setA(Number(e.target.value))}
            />
          </div>
          <div>
            <label className="text-xs text-zinc-400">b</label>
            <input
              className="mt-1 px-3 py-2 rounded bg-zinc-900 text-sm w-24 border border-zinc-800"
              value={b}
              onChange={(e) => setB(Number(e.target.value))}
            />
          </div>
          <div>
            <label className="text-xs text-zinc-400">tol</label>
            <input
              className="mt-1 px-3 py-2 rounded bg-zinc-900 text-sm w-36 border border-zinc-800"
              value={tol}
              onChange={(e) => setTol(Number(e.target.value))}
            />
          </div>
        </div>
        <div className="flex gap-2">
          <Button onClick={run}>Compute</Button>
          <Button onClick={animate}>Animate</Button>
          <Button variant="outline" onClick={() => setSteps([])}>
            Reset
          </Button>
        </div>
      </div>

      <div style={{ maxHeight: 340, overflow: "auto" }} className="rounded border border-zinc-800 p-2 bg-zinc-900/30">
        <SVGPlot fn={fnRef.current} domain={[a - 2, b + 2]} markers={markers} w={520} h={280} />
      </div>
    </div>
  );
}

/* ============================
   UI Sections (Hero, Features, Visualizers, Implementation, Calculations, Demo, FAQ, Developer, Footer)
   Each section has its own 3D canvas background bound to the logic above.
   ============================ */

/* ---------- small decorative cursor ---------- */
function CustomCursor() {
  const dotRef = useRef();
  const ringRef = useRef();
  useEffect(() => {
    const onMove = (e) => {
      const x = e.clientX;
      const y = e.clientY;
      if (dotRef.current) dotRef.current.style.transform = `translate3d(${x}px, ${y}px, 0)`;
      if (ringRef.current) ringRef.current.style.transform = `translate3d(${x - 14}px, ${y - 14}px, 0)`;
    };
    window.addEventListener("mousemove", onMove, { passive: true });
    return () => window.removeEventListener("mousemove", onMove);
  }, []);
  return (
    <>
      <div
        ref={ringRef}
        aria-hidden
        className="pointer-events-none fixed z-[9999] h-8 w-8 rounded-full border"
        style={{ borderColor: "rgba(26,188,156,0.18)", transform: "translate3d(-50px,-50px,0)", transition: "transform 80ms linear" }}
      />
      <div
        ref={dotRef}
        aria-hidden
        className="pointer-events-none fixed z-[9999] h-2 w-2 rounded-full"
        style={{ background: "#86efac", transform: "translate3d(-50px,-50px,0)", transition: "transform 25ms linear" }}
      />
    </>
  );
}

/* ---------- NAVBAR ---------- */
function Navbar({ onJump }) {
  const [isOpen, setIsOpen] = useState(false);

  const navItems = [
    { label: "Home", icon: <Home className="h-4 w-4 mr-2" /> },
    { label: "Why", icon: <Sparkles className="h-4 w-4 mr-2" /> },
    { label: "Visualizers", icon: <BarChart3 className="h-4 w-4 mr-2" /> },
    { label: "Labs", icon: <FlaskConical className="h-4 w-4 mr-2" /> },
    { label: "Calculations", icon: <CreditCard className="h-4 w-4 mr-2" /> },
   
  ];
  const navigate=useNavigate()

  return (
    <div className="fixed top-0 inset-x-0 z-50 bg-zinc-900 backdrop-blur-md border-b border-emerald-800/40">
      <div className={`${CONTAINER} flex items-center justify-between py-3`}>
        {/* Logo */}
        <div className="flex items-center gap-3">
          <div className="h-10 w-10 rounded-xl bg-zinc-900/60 flex items-center justify-center border border-zinc-800">
            <Rocket className="h-5 w-5 text-emerald-300" />
          </div>
          <div>
            <div className="text-sm font-semibold text-zinc-100">Numexis</div>
            <div className="text-xs text-zinc-400">
              Numerical Methods • Futuristic
            </div>
          </div>
        </div>

        {/* Desktop Nav */}
        <div className="hidden md:flex items-center gap-4">
          {navItems.map((item, idx) => (
            <button
              key={item.label}
              onClick={() => onJump(idx)}
              className="flex items-center text-sm px-3 py-2 rounded-md text-zinc-300 hover:text-white transition"
            >
              {item.icon}
              {item.label}
            </button>
          ))}
        </div>

        {/* Right side */}
        <div className="flex items-center gap-3">
         
          <Button
            onClick={() =>navigate("/chapters")}
            className="bg-emerald-300 cursor-pointer hover:bg-emerald-400 text-black hidden md:inline-flex"
          >
            Chapters
          </Button>

          {/* Mobile Hamburger */}
          <button
            className="md:hidden p-2 rounded-md text-zinc-300 hover:text-white focus:outline-none"
            onClick={() => setIsOpen(!isOpen)}
          >
            {isOpen ? <X className="h-6 w-6" /> : <Menu className="h-6 w-6" />}
          </button>
        </div>
      </div>

      {/* Mobile Menu Overlay */}
      <AnimatePresence>
        {isOpen && (
          <motion.div
            initial={{ x: "100%" }}
            animate={{ x: 0 }}
            exit={{ x: "100%" }}
            transition={{ type: "spring", stiffness: 300, damping: 30 }}
            className="fixed top-0 right-0 h-screen w-64 bg-zinc-900 backdrop-blur-md border-l border-emerald-800/40 shadow-xl md:hidden z-50"
          >
            <div className="flex flex-col p-6 space-y-6">
              {navItems.map((item, idx) => (
                <button
                  key={item.label}
                  onClick={() => {
                    onJump(idx);
                    setIsOpen(false);
                  }}
                  className="flex items-center text-lg font-medium text-zinc-200 hover:text-emerald-300 transition"
                >
                  {item.icon}
                  {item.label}
                </button>
              ))}

              <Button
                onClick={() => {
                  navigate("/chapters")
                  setIsOpen(false);
                }}
                className="bg-emerald-300 hover:bg-emerald-400 cursor-pointer text-black mt-4"
              >
                Chapters
              </Button>
            </div>
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  );
}
/* ---------- HERO ---------- */
function Hero({ onJump }) {
  return (
    <section className="relative min-h-screen flex items-center snap-start">
      <div className="absolute inset-0 -z-20">
        <Canvas camera={{ position: [0, 1.0, 9], fov: 52 }}>
          <ambientLight intensity={0.9} />
          <directionalLight intensity={0.7} position={[6, 8, 2]} />
          <Suspense fallback={<Loader text="Rendering…" />}>
            <ParticleField count={800} />
            <ContourField rings={22} scale={1.05} />
            <GradientSwarm n={260} radius={6} />
            <RootBeacons />
            <OrbitControls enablePan={false} enableZoom={false} autoRotate autoRotateSpeed={0.04} />
            <Preload all />
          </Suspense>
        </Canvas>
      </div>

      <div className={`${CONTAINER} relative z-20 py-28`}>
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-10 items-start">
          <div className="space-y-6">
            <div className={`${DARK_GLASS} inline-flex items-center gap-2 rounded-full px-3 py-1 text-xs text-emerald-300`}>
              <span className="h-2 w-2 rounded-full bg-emerald-400" />
              Numerical Methods • 3D • Interactive
            </div>

            <motion.h1 initial={{ y: 10, opacity: 0 }} animate={{ y: 0, opacity: 1 }} transition={{ duration: 0.9 }} className="text-6xl md:text-[76px] font-extrabold leading-tight text-zinc-100">
              Futuristic Visuals for <span className="text-emerald-300">Numerical Methods</span>
            </motion.h1>

            <motion.p initial={{ y: 8, opacity: 0 }} animate={{ y: 0, opacity: 1 }} transition={{ delay: 0.12 }} className="text-zinc-300 max-w-2xl text-lg">
              Explore convergence, optimization, PDE relaxation, eigenvectors, and root finding — brought to life with custom Three.js animations.
            </motion.p>

            <div className="flex flex-wrap gap-4 mt-4">
              <Button  onClick={() => onJump(2)} className="px-6 bg-zinc-700 cursor-pointer hover:bg-zinc-600 py-3">
                Open Visualizers <ArrowRight />
              </Button>
              <Button className="bg-white text-black hover:bg-gray-300 cursor-pointer" variant="outline" onClick={() => onJump(1)}>
                Why Numexis
              </Button>
            </div>

            <div className="mt-8 grid grid-cols-3 gap-3 max-w-md">
              <div className={`${DARK_GLASS} rounded-lg p-3`}>
                <div className="text-xs text-zinc-400">Algorithms</div>
                <div className="font-semibold text-zinc-100">30+</div>
              </div>
              <div className={`${DARK_GLASS} rounded-lg p-3`}>
                <div className="text-xs text-zinc-400">Live Demos</div>
                <div className="font-semibold text-zinc-100">20+</div>
              </div>
              <div className={`${DARK_GLASS} rounded-lg p-3`}>
                <div className="text-xs text-zinc-400">Exportable Labs</div>
                <div className="font-semibold text-zinc-100">Yes</div>
              </div>
            </div>
          </div>

          <div>
            <Card className="rounded-2xl overflow-hidden bg-zinc-900/50 border border-zinc-800">
              <div className="relative aspect-[16/9] bg-zinc-900/20">
                <Canvas camera={{ position: [1.2, 3.0, 3.8], fov: 55 }}>
                  <ambientLight intensity={1.0} />
                  <directionalLight position={[4, 4, 4]} intensity={0.6} />
                  <Suspense fallback={<Loader text="Preview…" />}>
                    <ContourField rings={10} scale={0.9} />
                    <OrbitControls enableZoom={false} enablePan={false} />
                  </Suspense>
                </Canvas>

                <div className="absolute left-4 bottom-4 text-sm text-zinc-300">
                  <div className="text-xs text-zinc-400">Featured demo</div>
                  <div className="font-semibold text-zinc-100">Gradient Descent</div>
                </div>
              </div>
              <CardContent className="bg-zinc-900/20 border-t border-zinc-800">
                <div className="flex items-center justify-between">
                  <div className="text-sm text-zinc-300">Visual convergence & landscapes</div>
                  <div>
                    <Button onClick={() => onJump(2)}>Open</Button>
                  </div>
                </div>
              </CardContent>
            </Card>
          </div>
        </div>
      </div>

      <div className="absolute bottom-10 left-1/2 -translate-x-1/2 z-30">
        <motion.div animate={{ y: [0, 8, 0] }} transition={{ repeat: Infinity, duration: 2 }} className="flex items-center gap-2 text-sm text-zinc-400">
          <ChevronDown />
          <span>Scroll to explore</span>
        </motion.div>
      </div>
    </section>
  );
}

/* ---------- FEATURES SECTION (each tile contains its own mini canvas) ---------- */
function FeaturesSection() {
  const features = [
    {
      title: "Optimization Landscapes",
      desc: "Level sets, gradient descent, and saddle behavior.",
      icon: Zap,
      scene: (idx) => <ContourField rings={12 + (idx % 3) * 2} scale={0.9} />,
    },
    {
      title: "PDE Relaxation",
      desc: "Jacobi & Gauss-Seidel-like smoothing in real time.",
      icon: Globe,
      scene: () => <RelaxationSurface W={48} H={34} amplitude={0.9} />,
    },
    {
      title: "Root Finding",
      desc: "Visual root markers and step-by-step iterations.",
      icon: Boxes,
      scene: () => <RootBeacons />,
    },
    {
      title: "Eigen Methods",
      desc: "Power method convergence toward dominant directions.",
      icon: Users,
      scene: () => <EigenArrows />,
    },
    {
      title: "High Performance",
      desc: "GPU-friendly rendering and optimized loops.",
      icon: Rocket,
      scene: () => <GradientSwarm n={100} radius={4} />,
    },
    {
      title: "Reproducible",
      desc: "Code-first demos and exportable labs.",
      icon: GitBranch,
      scene: () => <ContourField rings={8} scale={0.8} />,
    },
  ];

  return (
    <section className="relative min-h-screen flex items-center snap-start">
      <div className="absolute inset-0 -z-20">
        <Canvas camera={{ position: [0, 1.2, 9], fov: 58 }}>
          <ambientLight intensity={0.6} />
          <directionalLight position={[6, 8, 2]} intensity={0.5} />
          <Suspense fallback={null}>
            <ParticleField count={500} />
            <ContourField rings={10} scale={1.0} />
            <OrbitControls enablePan={false} enableZoom={false} autoRotate autoRotateSpeed={0.02} />
            <Preload all />
          </Suspense>
        </Canvas>
      </div>

      <div className={`${CONTAINER} py-20 relative z-10`}>
        <div className="text-center">
          <div className="text-xs uppercase tracking-wider text-zinc-400">Why Numexis</div>
          <h2 className="text-3xl font-semibold text-zinc-100 mt-2">A pro toolkit for teaching numerical methods</h2>
          <p className="mt-2 text-zinc-400 max-w-2xl mx-auto">From visual intuition to reproducible labs — fast, modern, and extensible.</p>
        </div>

        <div className="mt-8 grid grid-cols-1 md:grid-cols-3 gap-6">
          {features.map((f, i) => (
            <div key={i} className={`${DARK_GLASS} p-6 rounded-xl`}>
              <div className="flex items-center gap-3">
                <div className="h-12 w-12 rounded-lg bg-zinc-900/50 flex items-center justify-center border border-zinc-800">
                  {React.createElement(f.icon, { className: "h-5 w-5 text-emerald-300" })}
                </div>
                <div>
                  <div className="font-semibold text-zinc-100">{f.title}</div>
                  <div className="text-sm text-zinc-400 mt-1">{f.desc}</div>
                </div>
              </div>

              <div className="mt-4 h-28 relative rounded-md overflow-hidden border border-zinc-800">
                <Canvas camera={{ position: [0, 0, 6], fov: 60 }}>
                  <ambientLight intensity={0.8} />
                  <directionalLight position={[4, 4, 4]} intensity={0.5} />
                  <Suspense fallback={<Loader text="Mini scene…" />}>
                    {f.scene(i)}
                    <OrbitControls enablePan={false} enableZoom={false} />
                  </Suspense>
                </Canvas>
              </div>
            </div>
          ))}
        </div>
      </div>
    </section>
  );
}

/* ---------- VISUALIZERS ---------- */
function VisualizersSection() {
  return (
    <section className="relative min-h-screen flex items-start snap-start">
      <div className="absolute inset-0 -z-20">
        <Canvas camera={{ position: [0, 1.0, 9], fov: 56 }}>
          <ambientLight intensity={0.6} />
          <directionalLight position={[6, 8, 2]} intensity={0.5} />
          <Suspense fallback={null}>
            <RootBeacons />
            <GradientSwarm n={120} radius={4} />
            <OrbitControls enablePan={false} enableZoom={false} autoRotate autoRotateSpeed={0.015} />
            <Preload all />
          </Suspense>
        </Canvas>
      </div>

      <div className={`${CONTAINER} py-20 relative z-10`}>
        <div className="text-center">
          <div className="text-xs uppercase tracking-wider text-zinc-400">Demos</div>
          <h2 className="text-3xl font-semibold text-zinc-100 mt-2">Interactive visualizers in your browser</h2>
          <p className="mt-2 text-zinc-400 max-w-2xl mx-auto">Tune parameters, animate iterations, and export results.</p>
        </div>

        <div className="mt-8 grid grid-cols-1 lg:grid-cols-2 gap-8">
          <div className={`${DARK_GLASS} rounded-xl p-6`}>
            <div className="text-sm text-zinc-400 flex items-center gap-2">
              <Play /> Newton-Raphson
            </div>
            <div className="mt-4">
              <NewtonVisualizer />
            </div>
          </div>

          <div className={`${DARK_GLASS} rounded-xl p-6`}>
            <div className="text-sm text-zinc-400 flex items-center gap-2">
              <Play /> Bisection
            </div>
            <div className="mt-4">
              <BisectionVisualizer />
            </div>
          </div>
        </div>

        <div className="mt-10 grid grid-cols-1 md:grid-cols-3 gap-6">
          <div className={`${DARK_GLASS} p-6 rounded-xl`}>
            <div className="font-semibold text-zinc-100">Integrator</div>
            <div className="text-sm text-zinc-400 mt-2">Trapezoid & Simpson demos (coming soon)</div>
          </div>
          <div className={`${DARK_GLASS} p-6 rounded-xl`}>
            <div className="font-semibold text-zinc-100">Linear Solvers</div>
            <div className="text-sm text-zinc-400 mt-2">Gauss, LU, Jacobi, and comparisons</div>
          </div>
          <div className={`${DARK_GLASS} p-6 rounded-xl`}>
            <div className="font-semibold text-zinc-100">PDE Suite</div>
            <div className="text-sm text-zinc-400 mt-2">Relaxation, diffusion, and wave equation explorers</div>
          </div>
        </div>
      </div>
    </section>
  );
}

/* ---------- Implementation / Code Snippets ---------- */
function ImplementationSection() {
  return (
    <section className="relative min-h-screen flex items-center snap-start">
      <div className="absolute inset-0 -z-20">
        <Canvas camera={{ position: [0, 1.0, 8], fov: 52 }}>
          <ambientLight intensity={0.6} />
          <directionalLight position={[6, 8, 2]} intensity={0.5} />
          <Suspense>
            <EigenArrows />
            <ParticleField count={300} />
            <OrbitControls enablePan={false} enableZoom={false} autoRotate autoRotateSpeed={0.01} />
            <Preload all />
          </Suspense>
        </Canvas>
      </div>

      <div className={`${CONTAINER} py-20 relative  z-10`}>
        <div className="text-center">
          <div className="text-xs uppercase tracking-wider text-zinc-400">Implementation</div>
          <h2 className="text-3xl font-semibold text-zinc-100 mt-2">Pseudocode & snippets</h2>
          <p className="mt-2 text-zinc-400 max-w-2xl mx-auto">Copy-ready pseudocode, JS and Python references for each algorithm.</p>
        </div>

        <div className="mt-8 grid grid-cols-1 md:grid-cols-2 gap-6">
          <Card className="p-4 bg-zinc-900/50 border border-zinc-600">
            <CardHeader>
              <CardTitle className="text-white">Newton-Raphson (Pseudocode)</CardTitle>
            </CardHeader>
            <CardContent>
              <pre className="text-xs overflow-auto bg-gray-300 p-3 rounded">
{`Given f(x), f'(x), initial guess x0, tolerance eps
for k = 1..maxIter:
  x1 = x0 - f(x0)/f'(x0)
  if |x1 - x0| < eps: return x1
  x0 = x1
return x0`}
              </pre>
            </CardContent>
          </Card>

          <Card className="p-4 bg-zinc-900/50 border border-zinc-600">
            <CardHeader>
              <CardTitle className="text-white">Trapezoidal Rule (JavaScript)</CardTitle>
            </CardHeader>
            <CardContent>
              <pre className="text-xs overflow-auto bg-gray-300 p-3 rounded">
{`function trapezoid(f, a, b, n) {
  const h = (b - a) / n;
  let s = 0.5 * (f(a) + f(b));
  for (let i = 1; i < n; i++) s += f(a + i*h);
  return h * s;
}`}
              </pre>
            </CardContent>
          </Card>
        </div>
      </div>
    </section>
  );
}

/* ---------- Calculations / Mini calculators ---------- */
function CalculationsSection() {
  const [fx, setFx] = useState("x*x - 2");
  const [a, setA] = useState(1);
  const [b, setB] = useState(2);
  const [result, setResult] = useState(null);

  const computeBisection = () => {
    const f = makeSafeFn(fx);
    let lo = Number(a), hi = Number(b);
    if (!isFinite(f(lo)) || !isFinite(f(hi))) {
      setResult({ error: "Invalid function at endpoints" });
      return;
    }
    let mid = lo;
    for (let i = 0; i < 50; i++) {
      mid = (lo + hi) / 2;
      const fm = f(mid);
      if (Math.abs(fm) < 1e-6) break;
      if (f(lo) * fm <= 0) hi = mid;
      else lo = mid;
    }
    setResult({ root: mid });
  };

  return (
    <section className="relative min-h-screen flex items-center snap-start">
      <div className="absolute inset-0 -z-20">
        <Canvas camera={{ position: [0, 1.0, 8], fov: 52 }}>
          <ambientLight intensity={0.6} />
          <directionalLight position={[6, 8, 2]} intensity={0.5} />
          <Suspense>
            <RootBeacons />
            <WaveField cols={80} rows={26} amp={0.8} />
            <OrbitControls enablePan={false} enableZoom={false} autoRotate autoRotateSpeed={0.01} />
            <Preload all />
          </Suspense>
        </Canvas>
      </div>

      <div className={`${CONTAINER} py-20 relative z-10`}>
        <div className="text-center">
          <div className="text-xs uppercase tracking-wider text-zinc-400">Calculations</div>
          <h2 className="text-3xl font-semibold text-zinc-100 mt-2">Built-in calculators</h2>
          <p className="mt-2 text-zinc-400 max-w-2xl mx-auto">Quick calculators for roots, interpolation, and integration.</p>
        </div>

        <div className="mt-8 max-w-3xl  mx-auto">
          <Card className="p-6 bg-zinc-900/50 border border-zinc-500">
            <div className="grid grid-cols-1 md:grid-cols-3 gap-3">
              <Input value={fx} onChange={(e) => setFx(e.target.value)} className="col-span-2 bg-gray-400" />
              <div className="flex gap-2">
                <Input className="bg-gray-400" value={a} onChange={(e) => setA(e.target.value)} />
                <Input className="bg-gray-400" value={b} onChange={(e) => setB(e.target.value)} />
              </div>
            </div>
            <div className="mt-4 flex gap-3">
              <Button className="bg-white text-black hover:bg-gray-200 cursor-pointer" onClick={computeBisection}>Compute Bisection</Button>
              <Button className="cursor-pointer" variant="outline" onClick={() => { setFx("x*x-2"); setA(1); setB(2); setResult(null); }}>Reset</Button>
            </div>

            <div className="mt-4 text-sm text-zinc-300">
              {result ? (result.error ? <span className="text-amber-400">{result.error}</span> : <span>Root ≈ <strong style={{ color: THEME.warn }}>{Number(result.root).toFixed(6)}</strong></span>) : <span>Enter a function and range, then compute.</span>}
            </div>
          </Card>
        </div>
      </div>
    </section>
  );
}

/* ---------- Demo / Walkthrough ---------- */
function DemoSection() {
  return (
    <section className="relative min-h-screen flex items-center snap-start">
      <div className="absolute inset-0 -z-20">
        <Canvas camera={{ position: [0, 0.9, 8], fov: 52 }}>
          <ambientLight intensity={0.6} />
          <directionalLight position={[6, 8, 2]} intensity={0.6} />
          <Suspense>
            <ParticleField count={600} />
            <ContourField rings={14} scale={0.95} />
            <OrbitControls enablePan={false} enableZoom={false} autoRotate autoRotateSpeed={0.02} />
            <Preload all />
          </Suspense>
        </Canvas>
      </div>

      <div className={`${CONTAINER} py-20 relative z-10`}>
        <div className="text-center">
          <div className="text-xs uppercase tracking-wider text-zinc-400">Demo</div>
          <h2 className="text-3xl font-semibold text-zinc-100 mt-2">Guided walkthrough</h2>
          <p className="mt-2 text-zinc-400 max-w-2xl mx-auto">Short videos and guided tours to quickly onboard instructors and students.</p>
        </div>

        <div className="mt-8 grid grid-cols-1 md:grid-cols-3 gap-6">
          <Card className="p-4 bg-zinc-900/50 border-zinc-500">
            <div className="text-sm text-zinc-300">Intro</div>
            <div className="mt-2 text-zinc-100 font-semibold">What is Numexis</div>
            <div className="mt-3 text-sm text-zinc-400">A short 3-minute overview of features and mission.</div>
            <div className="mt-4"><Button className="bg-white cursor-pointer text-black hover:bg-gray-300">Play</Button></div>
          </Card>

          <Card className="p-4 bg-zinc-900/50 border-zinc-500">
            <div className="text-sm text-zinc-300">Deep Dive</div>
            <div className="mt-2 text-zinc-100 font-semibold">Visualizer walkthrough</div>
            <div className="mt-3 text-sm text-zinc-400">How to use the visualizer, edit equations, and export results.</div>
            <div className="mt-4"><Button className="bg-white cursor-pointer text-black hover:bg-gray-300">Play</Button></div>
          </Card>

          <Card className="p-4 bg-zinc-900/50 border-zinc-500">
            <div className="text-sm text-zinc-300">Instructor Pack</div>
            <div className="mt-2 text-zinc-100 font-semibold">Syllabus alignment</div>
            <div className="mt-3 text-sm text-zinc-400">Exercises, lab notes, and slides mapped to curricula.</div>
            <div className="mt-4"><Button className="bg-white cursor-pointer text-black hover:bg-gray-300" >Get Pack</Button></div>
          </Card>
        </div>
      </div>
    </section>
  );
}

/* ---------- FAQ ---------- */
function FAQSection() {
  const faqs = [
    { q: "Can I run these demos in class?", a: "Yes — everything runs in modern browsers. You can toggle heavy scenes for low-end devices." },
    { q: "Can I export data?", a: "Yes — demos provide CSV/JSON export and downloadable notebooks." },
    { q: "Is it open-source?", a: "Core visualizers are open source; see licensing for educational & commercial use." },
    { q: "Does it need WebGL 2?", a: "Standard WebGL is fine; performance scales with GPU capability." },
  ];

  return (
    <section className="relative min-h-screen flex items-center snap-start">
      <div className="absolute inset-0 -z-20">
        <Canvas camera={{ position: [0, 1.0, 8], fov: 56 }}>
          <ambientLight intensity={0.6} />
          <directionalLight position={[6, 8, 2]} intensity={0.45} />
          <Suspense>
            <GradientSwarm n={80} radius={3} />
            <OrbitControls enablePan={false} enableZoom={false} autoRotate autoRotateSpeed={0.015} />
            <Preload all />
          </Suspense>
        </Canvas>
      </div>

      <div className={`${CONTAINER} py-20 relative z-10`}>
        <div className="text-center">
          <div className="text-xs uppercase tracking-wider text-zinc-400">FAQ</div>
          <h2 className="text-3xl font-semibold text-zinc-100 mt-2">Common questions</h2>
          <p className="mt-2 text-zinc-400 max-w-2xl mx-auto">The essentials for instructors and integrators.</p>
        </div>

        <div className="mt-8 grid grid-cols-1 md:grid-cols-2 gap-6">
          {faqs.map((f, i) => (
            <Card key={i} className="p-6 bg-zinc-900/50 border-zinc-500">
              <div className="font-semibold text-zinc-100">{f.q}</div>
              <div className="mt-2 text-sm text-zinc-400">{f.a}</div>
            </Card>
          ))}
        </div>
      </div>
    </section>
  );
}

/* ---------- Developer / Contact ---------- */
function DeveloperSection() {
  const [msg, setMsg] = useState("");
  const submit = () => {
    // stub: integrate with email or backend later
    alert("Message sent (stub): " + msg.slice(0, 140));
    setMsg("");
  };

  return (
    <section className="relative min-h-screen flex items-center snap-start">
      <div className="absolute inset-0 -z-20">
        <Canvas camera={{ position: [0, 0.9, 8], fov: 52 }}>
          <ambientLight intensity={0.5} />
          <directionalLight position={[6, 8, 2]} intensity={0.5} />
          <Suspense>
            <ParticleField count={400} />
            <ContourField rings={8} scale={0.9} />
            <OrbitControls enablePan={false} enableZoom={false} autoRotate autoRotateSpeed={0.008} />
            <Preload all />
          </Suspense>
        </Canvas>
      </div>

      <div className={`${CONTAINER} py-20 relative z-10`}>
        <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
          <Card className="md:col-span-2 p-6 flex-1 bg-zinc-900/50 border-zinc-600">
            <div className="text-xs text-zinc-400 uppercase">Developer</div>
            <div className="text-2xl font-semibold text-zinc-100 mt-2">Eswar Chinthakayala</div>
            <div className="mt-2 text-sm text-zinc-400">Full-Stack / Frontend / Backend / Software Developer</div>

            <div className="mt-6">
              <div className="text-xs text-zinc-400">Contact</div>
              <div className="mt-2 flex items-center gap-3 text-sm text-zinc-300">
                <Mail className="h-4 w-4" /> eswarch2004y@gmail.com
              </div>

              <div className="mt-4">
                <div className="text-xs text-zinc-400">Message</div>
                <textarea value={msg} onChange={(e) => setMsg(e.target.value)} rows={4} className="w-full mt-2 p-3 rounded bg-gray-800 border border-zinc-800 text-sm" />
                <div className="mt-3 flex gap-2"><Button className="bg-zinc-600 cursor-pointer hover:bg-zinc-500" onClick={submit}>Say Hello</Button><Button variant="outline" className="cursor-pointer">View GitHub</Button></div>
              </div>
            </div>
          </Card>

          <Card className="p-6  bg-zinc-900/50 border-zinc-500">
            <div className="text-xs text-zinc-400">Resources</div>
            <div className="mt-3 flex flex-col gap-2">
              <a className="text-sm text-emerald-200">Docs</a>
              <a className="text-sm text-emerald-200">GitHub</a>
              <a className="text-sm text-emerald-200">Instructor Pack</a>
            </div>
          </Card>
        </div>
      </div>
    </section>
  );
}

/* ---------- Footer ---------- */
function Footer() {
  return (
    <footer className="bg-zinc-950/50 border-t border-zinc-800 text-zinc-400">
      <div className={`${CONTAINER} py-10`}>
        <div className="flex flex-col md:flex-row items-start md:items-center justify-between gap-4">
          <div>
            <div className="font-semibold text-zinc-100">Numexis</div>
            <div className="text-sm text-zinc-400">Numerical Methods • Futuristic 3D</div>
          </div>
          <div className="flex gap-3">
            <a className="text-sm hover:text-white">Docs</a>
            <a className="text-sm hover:text-white">GitHub</a>
            <a className="text-sm hover:text-white">Contact</a>
          </div>
        </div>
        <div className="mt-6 text-sm text-zinc-500">© {new Date().getFullYear()} Numexis — Built for teaching & research</div>
      </div>
    </footer>
  );
}

/* ============================
   MAIN LandingPage Component (one giant file)
   - uses IntersectionObserver for section nav
   - assembles all sections above
   ============================ */

export default function LandingPage() {
  const sections = [
    { id: "hero", ref: useRef(null) },
    { id: "why", ref: useRef(null) },
    { id: "visualizers", ref: useRef(null) },
    { id: "implementation", ref: useRef(null) },
    { id: "calculations", ref: useRef(null) },
    { id: "demo", ref: useRef(null) },
    { id: "faq", ref: useRef(null) },
    { id: "dev", ref: useRef(null) },
  ];

  const [active, setActive] = useState(0);
  useEffect(() => {
    const obs = new IntersectionObserver(
      (entries) => {
        entries.forEach((entry) => {
          const idx = sections.findIndex((s) => s.ref.current === entry.target);
          if (entry.isIntersecting && idx >= 0) setActive(idx);
        });
      },
      { threshold: 0.55 }
    );

    sections.forEach((s) => {
      if (s.ref.current) obs.observe(s.ref.current);
    });
    return () => obs.disconnect();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  const jumpTo = (idx) => {
    const s = sections[idx];
    if (s && s.ref.current) s.ref.current.scrollIntoView({ behavior: "smooth", block: "start" });
  };

  return (
    <div style={{ color: THEME.text }} className="min-h-screen bg-zinc-950 text-zinc-100">
      <CustomCursor />
      <Navbar onJump={(i) => jumpTo(i)} />

      {/* floating side dot nav */}
      <div className="fixed right-6 top-1/2 -translate-y-1/2 z-50 hidden md:flex flex-col gap-3">
        {sections.map((s, i) => (
          <button
            key={s.id}
            onClick={() => jumpTo(i)}
            className={`h-3 w-3 rounded-full transition-all ${active === i ? "h-6 w-1.5 rounded-md bg-zinc-100" : "bg-zinc-600 hover:bg-zinc-400"}`}
            aria-label={s.id}
          />
        ))}
      </div>

      <main className="snap-y snap-mandatory h-screen overflow-y-auto scroll-smooth">
        <div ref={sections[0].ref} data-section="hero">
          <Hero onJump={(i) => jumpTo(i)} />
        </div>

        <div ref={sections[1].ref} data-section="why">
          <FeaturesSection />
        </div>

        <div ref={sections[2].ref} data-section="visualizers">
          <VisualizersSection />
        </div>

        <div ref={sections[3].ref} data-section="implementation">
          <ImplementationSection />
        </div>

        <div ref={sections[4].ref} data-section="calculations">
          <CalculationsSection />
        </div>

        <div ref={sections[5].ref} data-section="demo">
          <DemoSection />
        </div>

        <div ref={sections[6].ref} data-section="faq">
          <FAQSection />
        </div>

        <div ref={sections[7].ref} data-section="dev">
          <DeveloperSection />
        </div>

        <div>
          <Footer />
        </div>
      </main>
    </div>
  );
}
