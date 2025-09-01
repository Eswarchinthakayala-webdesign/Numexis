// src/pages/Chapter3_Approximations.jsx
import React, { useState, useMemo, useRef, useEffect } from "react";
import { motion } from "framer-motion";
import Plot from "react-plotly.js";
import { Canvas, useFrame } from "@react-three/fiber";
import { OrbitControls } from "@react-three/drei";
import * as THREE from "three";
import {
  BookOpen,
  ChevronRight,
  ChevronDown,
  Calculator,
  Eye,
  AlertTriangle,
  CircleDot,
  Zap,
  Layers,
  List,
  Terminal,
  GitMerge,
  Sparkles,
} from "lucide-react";

import { Card, CardHeader, CardTitle, CardContent } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import BottomBar from "../components/BottomBar";

/**
 * Chapter 3 — Approximations & Round-Off Errors
 *
 * Design: follows Chapter 2 zinc-dark theme and layout.
 * Features:
 * - Extensive documentation (like a real textbook chapter)
 * - Interactive calculators:
 *   - Significant figures formatter
 *   - Rounding simulator (round vs truncate)
 *   - Error calculator (absolute/relative/percent)
 *   - Catastrophic cancellation demo
 *   - Summation precision demo (pairwise vs naive)
 * - Visualizations: Plotly charts and a simple 3D "error bars" canvas
 * - Uses lucide-react icons for UI affordances
 *
 * NOTE: This file is intentionally verbose and content-rich to behave like a chapter.
 */

// ---------------- Theme (keep consistent with Chapter 2) ----------------
const theme = {
  background: "bg-zinc-900",
  panelBg: "#0f1720",
  text: "#e6eef6",
  accent1: "#60a5fa",
  accent2: "#34d399",
  accent3: "#f59e0b",
};

// ---------------- motion helpers ----------------
const cardMotion = {
  initial: { opacity: 0, y: 12, scale: 0.995 },
  enter: { opacity: 1, y: 0, scale: 1, transition: { duration: 0.45, ease: [0.2, 0.8, 0.2, 1] } },
  hover: { scale: 1.02, y: -6, transition: { duration: 0.25 } },
};

// ---------------- small numeric helpers ----------------
function roundToSignificantFigures(value, sig) {
  if (!isFinite(value) || value === 0) return value;
  const sign = Math.sign(value) || 1;
  const abs = Math.abs(value);
  const power = Math.floor(Math.log10(abs));
  const factor = Math.pow(10, power - sig + 1);
  return sign * Math.round(abs / factor) * factor;
}

function truncateToSignificantFigures(value, sig) {
  if (!isFinite(value) || value === 0) return value;
  const sign = Math.sign(value) || 1;
  const abs = Math.abs(value);
  const power = Math.floor(Math.log10(abs));
  const factor = Math.pow(10, power - sig + 1);
  return sign * Math.floor(abs / factor) * factor;
}

function toFixedSig(value, sig) {
  if (!isFinite(value)) return String(value);
  if (value === 0) return "0".padEnd(sig, "0");
  const sign = value < 0 ? "-" : "";
  const a = Math.abs(value);
  const power = Math.floor(Math.log10(a));
  const shift = sig - 1 - power;
  const scaled = Math.round(a * Math.pow(10, shift));
  const s = String(scaled);
  if (power >= sig - 1) {
    // no decimal point
    return sign + s + "0".repeat(power - (sig - 1));
  } else {
    const intPartLen = power + 1;
    const intPart = s.slice(0, intPartLen) || "0";
    const fracPart = s.slice(intPartLen).padEnd(sig - intPartLen, "0");
    return sign + intPart + "." + fracPart;
  }
}

// naive summation vs pairwise (Kahan could be added)
function naiveSum(arr) {
  let s = 0;
  for (let i = 0; i < arr.length; i++) s += arr[i];
  return s;
}

function pairwiseSum(arr) {
  // simple divide and conquer pairwise summation
  function recursive(a, l, r) {
    if (r - l === 1) return a[l];
    const m = Math.floor((l + r) / 2);
    return recursive(a, l, m) + recursive(a, m, r);
  }
  if (arr.length === 0) return 0;
  return recursive(arr, 0, arr.length);
}

function kahanSum(arr) {
  let sum = 0.0;
  let c = 0.0;
  for (let i = 0; i < arr.length; i++) {
    const y = arr[i] - c;
    const t = sum + y;
    c = (t - sum) - y;
    sum = t;
  }
  return sum;
}

// ---------------- 3D error bars for visualization ----------------
function ErrorBars3D({ seed = 42 }) {
  const ref = useRef();
  useFrame(({ clock }) => {
    if (!ref.current) return;
    ref.current.rotation.y = Math.sin(clock.elapsedTime * 0.2) * 0.15;
    ref.current.rotation.x = Math.sin(clock.elapsedTime * 0.12) * 0.08;
  });

  // deterministic pseudo-random generator so visualization is stable
  const prng = useMemo(() => {
    let s = seed;
    return () => {
      s = (s * 9301 + 49297) % 233280;
      return s / 233280;
    };
  }, [seed]);

  const bars = useMemo(() => {
    const arr = [];
    for (let i = 0; i < 12; i++) {
      arr.push({
        x: (i - 5.5) * 0.45,
        h: 0.2 + prng() * 1.6,
        color: `rgba(${120 + Math.floor(prng() * 130)}, ${160 + Math.floor(prng() * 60)}, ${80 + Math.floor(prng() * 120)}, 1)`,
      });
    }
    return arr;
  }, [prng]);

  return (
    <group ref={ref}>
      {bars.map((b, i) => (
        <mesh key={i} position={[b.x, b.h / 2, 0]}>
          <boxGeometry args={[0.35, b.h, 0.35]} />
          <meshStandardMaterial color={b.color} />
        </mesh>
      ))}
      <ambientLight intensity={0.8} />
      <directionalLight position={[2, 4, 2]} intensity={0.6} />
    </group>
  );
}

// ---------------- Expanded Documentation (textbook-style) ----------------
const docText = {
  "chapterTitle": "Chapter 3 — Approximations and Round-Off Errors",
  "overview": `This chapter explores sources of numerical error that appear in computing: rounding, truncation, propagation, and cancellation. We develop an intuition for significant figures, accuracy vs precision, formal error definitions, and practical techniques to mitigate loss of precision. Visualizations and experiments show how computations behave in finite precision. Practical exercises and calculators are included to reinforce concepts.`,
  "sections": {
    "3.1 Significant Figures": `
Definition
• Significant figures (sig figs) count the digits that carry meaning contributing to precision.
• They are used in engineering and experimental reporting to indicate the trustworthiness of a measurement.

Rules:
1. Non-zero digits are significant: 123 → 3 sig figs
2. Zeros between non-zero digits are significant: 1002 → 4 sig figs
3. Leading zeros are NOT significant: 0.00456 → 3 sig figs
4. Trailing zeros after a decimal are significant: 45.600 → 5 sig figs
5. Trailing zeros without a decimal are ambiguous: 1500 → could be 2,3, or 4 sig figs. Use scientific notation to clarify (1.5 × 10^3 → 2 sig figs).

Practice:
• Convert and round numbers to N significant figures
• Understand how arithmetic (especially subtraction) can change effective sig figs.`,
    "3.2 Accuracy and Precision": `
Accuracy
• How close a measured or computed value is to the true value.
Precision
• The repeatability of measurements — how close repeated measurements are to each other.
Interpretation
• High accuracy + high precision = desirable
• High precision + low accuracy = systematic bias
• Low precision + high accuracy = inconsistent measurement but unbiased on average

Visual metaphors:
• Dartboard: cluster around bullseye (accurate & precise), cluster away from center (precise but biased), scattered (imprecise).
Applications:
• When repeatedly computing results (Monte Carlo or iterative solvers), estimate both mean (bias) and variance (precision).`,
    "3.3 Error Definitions": `
True value (T): theoretical exact number (often unknown)
Approximate value (A): measured or computed number

Absolute Error = |T - A|
Relative Error = |T - A| / |T|
Percentage Error = Relative Error × 100%

Stability, Consistency, Convergence (short)
• Stability: algorithm's sensitivity to small perturbations (rounding or input noise)
• Consistency: discretization error goes to zero as step sizes go to zero (for numerical PDEs/ODEs)
• Convergence: computed solution approaches true solution as resolution increases

Propagation:
• Errors propagate through arithmetic; sometimes they amplify (unstable) or dampen.
• Use condition numbers (for linear systems) to estimate amplification.`,
    "3.4 Round-Off Errors": `
Finite precision:
• Computers represent floating-point numbers with finite mantissa length (e.g., IEEE 754 double ~ 53 bits mantissa).
• Not all decimal fractions can be represented exactly (e.g., 0.1)

Round-off:
• When a real number cannot be represented exactly, it's rounded to nearest representable number.
• Truncation is another source: using finite series terms, finite step sizes, or discretization.

Catastrophic cancellation:
• Subtraction of nearly equal numbers removes leading significant digits and amplifies relative error.
• Example: computing f(x) = sqrt(x+1) - sqrt(x) for large x; instead compute via algebraic manipulation to avoid cancellation.

Practical mitigation:
• Reorder summations: sum small numbers first or use pairwise/Kahan summation
• Use higher precision in intermediate steps
• Reformulate expressions to avoid subtracting nearly equal quantities`,
    "problems": `
A. Compute the number of significant figures for these:
  1) 0.000450
  2) 3.200 × 10^4
  3) 1200 (interpretation)

B. Calculate errors:
  True = π ≈ 3.141592653589793, Approx = 22/7
  Compute absolute, relative, percent.

C. Summation experiment:
  Sum 1e6 times 1e-6 using naive and Kahan summation. Compare exact and computed totals.

D. Catastrophic cancellation example:
  Evaluate f(x) = sqrt(x^2 + 1) - x for large x using direct subtraction and algebraic manipulation:
    f = 1 / (sqrt(x^2 + 1) + x)
  Compare results for x = 1e6.`,
  },
};

// ---------------- Documentation Panel (expanded & interactive) ----------------
function DocumentationPanelLarge({ startOpen = "3.1 Significant Figures" }) {
  const [openTopic, setOpenTopic] = useState(startOpen);

  // detailed content prepared as long strings (makes this panel large)
  const bigNotes = useMemo(() => {
    return {
      "3.1 Significant Figures": `
Significant figures (sig figs) express the precision of a measurement or computed quantity.
They tell the reader which digits are trustworthy.

Why they matter:
• When reporting results, sig figs indicate the uncertainty.
• When chaining computations, number of sig figs can decrease (especially after subtraction).

Example rules and exercises:
• Round 0.004567 to 2 sig figs → 0.0046
• Round 12345 to 3 sig figs → 12300 (or 1.23 × 10^4 to avoid ambiguity)

Rounded vs truncated:
• Rounding: round half up (or banker’s rounding in some contexts)
• Truncation: simply cut off digits — this introduces bias

Best practice:
• Use scientific notation for clarity when trailing zeros are significant.
• Propagate uncertainty using absolute/relative error rather than sig figs when doing analysis.`,
      "3.2 Accuracy and Precision": `
Accuracy and precision are distinct:
• Accuracy: closeness to true value (low bias)
• Precision: repeatability (low variance)

In numerical computation:
• Accuracy: is our algorithm consistent and convergent? Is there bias?
• Precision: how many bits of mantissa are we using? Are our iterative steps consistent?

Recommendations:
• Always report both point estimates and uncertainty (stdev, CI)
• Use multiple independent methods to validate results (different algorithms, different step sizes)`,
      "3.3 Error Definitions": `
Formal definitions recap:
• Absolute error: |T - A|
• Relative error: |T - A| / |T|
• Percentage error: * 100%

Condition number (example for linear systems):
• For Ax = b, the condition number κ(A) approximates amplification of input error to solution error.
• A large κ(A) => problem is ill-conditioned.

Error propagation example:
• Suppose A and B are measured with small errors δA, δB. For C = A+B, absolute errors add; for multiplication relative errors add approximately.`,
      "3.4 Round-Off Errors": `
Round-off emerges from limited mantissa length and base conversion:
• Base-10 decimals may not map exactly to binary floating point.

Important behavior:
• Rounding error is bounded by machine epsilon (ε): the difference between 1 and the next representable float.
• For IEEE double, ε ≈ 2.22e-16.

Strategies:
1) Use stable algorithms (algebraic reformulation).
2) Accumulate sums in increasing magnitude or use compensated summation.
3) Choose an appropriate precision for problem size.

Examples and exercises follow in the interactive widgets below.
`,
      "Problems & Exercises": `
1) Sig figs exercises: identify and round given numbers to N sig figs
2) Error calculation tasks for familiar constants (π, e, sqrt(2))
3) Summation experiment — see interactive widget
4) Catastrophic cancellation — see interactive widget`,
    };
  }, []);

  return (
    <div className="w-full mt-8 bg-zinc-900 border border-zinc-700 rounded-xl shadow-lg">
      <div className="p-4 border-b border-zinc-700 flex items-center justify-between">
        <div className="flex items-center gap-3">
          <BookOpen className="w-5 h-5 text-cyan-400" />
          <div>
            <div className="text-zinc-100 font-semibold">Chapter 3 Documentation</div>
            <div className="text-zinc-400 text-xs">Approximations, round-off, and practical mitigation</div>
          </div>
        </div>
        <div className="text-zinc-400 text-sm">Prefer visualization & calculation → practical learning</div>
      </div>

      <div className="grid grid-cols-1 md:grid-cols-3">
        {/* Left: Topics list */}
        <div className="col-span-1 border-r border-zinc-800">
          {Object.keys(bigNotes).map((topic) => (
            <div key={topic} className="border-b border-zinc-800">
              <button
                onClick={() => setOpenTopic(openTopic === topic ? null : topic)}
                className={`w-full text-left px-4 py-3 hover:bg-zinc-800/50 flex items-center justify-between ${openTopic === topic ? "bg-zinc-800/20" : ""}`}
              >
                <div className="flex items-center gap-3">
                  <List className="w-4 h-4 text-zinc-300" />
                  <div className="text-zinc-100">{topic}</div>
                </div>
                {openTopic === topic ? <ChevronDown className="w-4 h-4 text-zinc-400" /> : <ChevronRight className="w-4 h-4 text-zinc-400" />}
              </button>
            </div>
          ))}
        </div>

        {/* Middle: Content */}
        <div className="col-span-2 md:col-span-2 p-4">
          {openTopic ? (
            <>
              <motion.h3 initial={{ opacity: 0, y: -6 }} animate={{ opacity: 1, y: 0 }} className="text-xl font-semibold text-zinc-100 mb-3">
                {openTopic}
              </motion.h3>
              <motion.div initial={{ opacity: 0 }} animate={{ opacity: 1 }} className="text-zinc-300 whitespace-pre-line leading-relaxed">
                {bigNotes[openTopic]}
              </motion.div>

              {/* contextual small examples for each topic */}
              {openTopic === "3.1 Significant Figures" && (
                <div className="mt-4 grid grid-cols-1 md:grid-cols-2 gap-3">
                  <div className="bg-zinc-800 p-3 rounded">
                    <div className="text-zinc-200 font-medium mb-2">Formatter — Round to N sig figs</div>
                    <SigFigFormatter />
                  </div>
                  <div className="bg-zinc-800 p-3 rounded">
                    <div className="text-zinc-200 font-medium mb-2">Truncate vs Round (visual)</div>
                    <TruncateVsRoundDemo />
                  </div>
                </div>
              )}

              {openTopic === "3.3 Error Definitions" && (
                <div className="mt-4 grid grid-cols-1 md:grid-cols-2 gap-3">
                  <div className="bg-zinc-800 p-3 rounded">
                    <div className="text-zinc-200 font-medium mb-2">Error calculator</div>
                    <ErrorCalculator />
                  </div>
                  <div className="bg-zinc-800 p-3 rounded">
                    <div className="text-zinc-200 font-medium mb-2">Condition number demo (simple)</div>
                    <ConditionNumberDemo />
                  </div>
                </div>
              )}

              {openTopic === "3.4 Round-Off Errors" && (
                <div className="mt-4 grid grid-cols-1 md:grid-cols-3 gap-3">
                  <div className="bg-zinc-800 p-3 rounded col-span-1">
                    <div className="text-zinc-200 font-medium mb-2">Catastrophic Cancel Demo</div>
                    <CatastrophicCancellationDemo />
                  </div>
                  <div className="bg-zinc-800 p-3 rounded col-span-2">
                    <div className="text-zinc-200 font-medium mb-2">Summation Precision Demo (Naive, Pairwise, Kahan)</div>
                    <SummationPrecisionDemo />
                  </div>
                </div>
              )}
            </>
          ) : (
            <div className="text-zinc-400">Select a topic from the left to read details, examples and interactive exercises.</div>
          )}
        </div>
      </div>
    </div>
  );
}

// ---------------- UI Widgets ----------------

// 1) Significant Figures formatter UI
function SigFigFormatter() {
  const [val, setVal] = useState("12345.6789");
  const [sig, setSig] = useState(4);
  const [rounded, setRounded] = useState("");
  const [truncated, setTruncated] = useState("");
  const [toFixed, setToFixed] = useState("");

  useEffect(() => {
    let v = Number(val);
    if (!isFinite(v)) {
      setRounded("Invalid");
      setTruncated("Invalid");
      setToFixed("Invalid");
      return;
    }
    setRounded(String(roundToSignificantFigures(v, sig)));
    setTruncated(String(truncateToSignificantFigures(v, sig)));
    setToFixed(String(toFixedSig(v, sig)));
  }, [val, sig]);

  return (
    <div>
      <div className="flex gap-2">
        <Input value={val} onChange={(e) => setVal(e.target.value)} className="bg-zinc-700 text-white" />
        <Input type="number" value={sig} onChange={(e) => setSig(Number(e.target.value))} className="w-24 bg-zinc-700 text-white" />
      </div>

      <div className="mt-3 text-sm text-zinc-300 space-y-2">
        <div>Round → <span className="text-zinc-100 font-medium">{rounded}</span></div>
        <div>Truncate → <span className="text-zinc-100 font-medium">{truncated}</span></div>
        <div>Formatted → <span className="text-zinc-100 font-medium">{toFixed}</span></div>
      </div>
    </div>
  );
}

// 2) Truncate vs round demo with a small plot
function TruncateVsRoundDemo() {
  const [base, setBase] = useState(0.00456789);
  const [maxSig, setMaxSig] = useState(8);
  const xs = useMemo(() => {
    const arr = [];
    for (let s = 1; s <= maxSig; s++) {
      arr.push(s);
    }
    return arr;
  }, [maxSig]);

  const rounded = xs.map((s) => roundToSignificantFigures(base, s));
  const truncated = xs.map((s) => truncateToSignificantFigures(base, s));

  return (
    <div>
      <div className="flex gap-2">
        <Input value={String(base)} onChange={(e) => setBase(Number(e.target.value))} className="bg-zinc-700 text-white" />
        <Input type="number" value={maxSig} onChange={(e) => setMaxSig(Number(e.target.value))} className="w-24 bg-zinc-700 text-white" />
      </div>

      <div className="mt-3">
        <Plot
          data={[
            { x: xs, y: rounded, mode: "lines+markers", name: "rounded" },
            { x: xs, y: truncated, mode: "lines+markers", name: "truncated" },
          ]}
          layout={{ paper_bgcolor: theme.panelBg, plot_bgcolor: theme.panelBg, font: { color: theme.text }, title: "Round vs Truncate (value vs sig figs)", height: 260 }}
          useResizeHandler={true}
          style={{ width: "100%" }}
        />
      </div>
    </div>
  );
}

// 3) Error calculator
function ErrorCalculator() {
  const [trueVal, setTrueVal] = useState("3.141592653589793");
  const [approxVal, setApproxVal] = useState("22/7");
  const compute = useMemo(() => {
    let T = Number(trueVal);
    let A;
    try {
      // allow expressions like 22/7
      A = eval(approxVal);
    } catch (e) {
      A = Number(approxVal);
    }
    if (!isFinite(T) || !isFinite(A)) return null;
    const abs = Math.abs(T - A);
    const rel = Math.abs(T - A) / Math.abs(T || 1);
    const pct = rel * 100;
    return { T, A, abs, rel, pct };
  }, [trueVal, approxVal]);

  return (
    <div>
      <div className="flex gap-2">
        <Input value={trueVal} onChange={(e) => setTrueVal(e.target.value)} className="bg-zinc-700 text-white" />
        <Input value={approxVal} onChange={(e) => setApproxVal(e.target.value)} className="bg-zinc-700 text-white" />
      </div>

      {compute ? (
        <div className="mt-3 text-sm text-zinc-300 space-y-1">
          <div>True (T): <span className="text-zinc-100">{compute.T}</span></div>
          <div>Approx (A): <span className="text-zinc-100">{compute.A}</span></div>
          <div>Absolute Error: <span className="text-zinc-100">{compute.abs}</span></div>
          <div>Relative Error: <span className="text-zinc-100">{compute.rel}</span></div>
          <div>Percentage Error: <span className="text-zinc-100">{compute.pct.toFixed(6)}%</span></div>
        </div>
      ) : (
        <div className="mt-3 text-sm text-zinc-400">Invalid inputs</div>
      )}
    </div>
  );
}

// 4) Condition number demo (very simple for 2x2)
function ConditionNumberDemo() {
  const [Atext, setAtext] = useState("[[2,1],[1,3]]");
  const [btext, setBtext] = useState("[3,7]");
  const [cond, setCond] = useState(null);
  const [sol, setSol] = useState(null);

  useEffect(() => {
    try {
      const A = JSON.parse(Atext);
      const b = JSON.parse(btext);
      // compute solution by simple elimination (provided earlier function)
      const x = gaussElimination(A, b);
      setSol(x);
      // approximate condition number via ratio of norms: ||A|| * ||A^{-1}||
      // compute inverse for 2x2 or general via gauss (crude)
      const inv = invertMatrix(A);
      if (inv) {
        const normA = matrixNormInf(A);
        const normInv = matrixNormInf(inv);
        setCond(normA * normInv);
      } else {
        setCond(null);
      }
    } catch (e) {
      setCond(null);
      setSol(null);
    }
  }, [Atext, btext]);

  return (
    <div>
      <div className="flex gap-2">
        <Input value={Atext} onChange={(e) => setAtext(e.target.value)} className="bg-zinc-700 text-white" />
        <Input value={btext} onChange={(e) => setBtext(e.target.value)} className="bg-zinc-700 text-white" />
      </div>

      <div className="mt-3 text-sm text-zinc-300">
        <div>Solution x: <span className="text-zinc-100">{sol ? "[" + sol.map(v => v.toFixed(6)).join(", ") + "]" : "—"}</span></div>
        <div>Estimated condition number (∞-norm): <span className="text-zinc-100">{cond ? cond.toFixed(6) : "—"}</span></div>
        <div className="text-zinc-400 text-xs mt-2">Note: condition number estimates sensitivity of Ax=b to perturbations.</div>
      </div>
    </div>
  );
}

// helper matrix utilities used by above
function matrixNormInf(A) {
  let max = 0;
  for (let i = 0; i < A.length; i++) {
    let s = 0;
    for (let j = 0; j < A[i].length; j++) s += Math.abs(A[i][j]);
    if (s > max) max = s;
  }
  return max;
}

function invertMatrix(A) {
  const n = A.length;
  // only handle square matrices
  if (!A.every(row => row.length === n)) return null;
  // build augmented matrix
  let M = A.map((r, i) => [...r, ...identityRow(n, i)]);
  // gaussian elimination
  for (let k = 0; k < n; k++) {
    // pivot
    let pivot = k;
    for (let i = k + 1; i < n; i++) if (Math.abs(M[i][k]) > Math.abs(M[pivot][k])) pivot = i;
    if (Math.abs(M[pivot][k]) < 1e-16) return null;
    [M[k], M[pivot]] = [M[pivot], M[k]];
    const div = M[k][k];
    for (let j = 0; j < 2 * n; j++) M[k][j] /= div;
    for (let i = 0; i < n; i++) if (i !== k) {
      const mult = M[i][k];
      for (let j = 0; j < 2 * n; j++) M[i][j] -= mult * M[k][j];
    }
  }
  // extract inverse
  return M.map(r => r.slice(n));
}

function identityRow(n, idx) {
  const row = new Array(n).fill(0); row[idx] = 1; return row;
}

// 5) Catastrophic cancellation demo
function CatastrophicCancellationDemo() {
  const [x, setX] = useState(1e6);
  const [direct, setDirect] = useState(null);
  const [stable, setStable] = useState(null);

  // Compute single values for the chosen x
  useEffect(() => {
    const dx = Math.sqrt(x * x + 1) - x; // unstable form
    const alt = 1 / (Math.sqrt(x * x + 1) + x); // stable form
    setDirect(dx);
    setStable(alt);
  }, [x]);

  // Generate dynamic range of xs around user input
  const { xs, directVals, stableVals } = useMemo(() => {
    const xs = [];
    const directVals = [];
    const stableVals = [];
    // create ~20 points log-spaced from x/100 to x*100
    for (let i = -2; i <= 2; i++) {
      const val = x * Math.pow(10, i / 5); // finer log spacing
      xs.push(val);
      directVals.push(Math.sqrt(val * val + 1) - val);
      stableVals.push(1 / (Math.sqrt(val * val + 1) + val));
    }
    return { xs, directVals, stableVals };
  }, [x]);

  return (
    <div className="bg-zinc-900 p-5 rounded-xl shadow-lg border border-zinc-700">
      <h3 className="text-xl font-semibold text-zinc-100 mb-3">
        Catastrophic Cancellation Demo
      </h3>

      {/* Input + Button */}
      <div className="flex items-center gap-2">
        <Input
          type="number"
          value={x}
          onChange={(e) => setX(Number(e.target.value))}
          className="bg-zinc-800 text-white border-zinc-600"
        />
        <Button
          onClick={() => setX(prev => prev * 10)}
          className="bg-zinc-700 hover:bg-zinc-600"
        >
          Multiply X by 10
        </Button>
      </div>

      {/* Results */}
      <div className="mt-4 text-sm text-zinc-300 space-y-1">
        <div>
          Direct: <span className="text-emerald-400">{direct}</span>
        </div>
        <div>
          Stable (rearranged): <span className="text-sky-400">{stable}</span>
        </div>
        <div className="text-zinc-500 text-xs">
          Note: For large x the direct subtraction loses significant digits;
          algebraic rearrangement is numerically stable.
        </div>
      </div>

      {/* Plot */}
      <div className="mt-5">
        <Plot
          data={[
            {
              x: xs,
              y: directVals,
              mode: "lines+markers",
              name: "Direct",
              line: { color: "#22c55e" },
              marker: { color: "#22c55e" }
            },
            {
              x: xs,
              y: stableVals,
              mode: "lines+markers",
              name: "Stable",
              line: { color: "#0ea5e9" },
              marker: { color: "#0ea5e9" }
            }
          ]}
          layout={{
            title: `Catastrophic Cancellation (around x=${x})`,
            paper_bgcolor: "#18181b",
            plot_bgcolor: "#18181b",
            font: { color: "#e4e4e7" },
            xaxis: { type: "log", title: "x" },
            yaxis: { type: "log", title: "f(x)" },
            height: 320,
            margin: { t: 50, l: 60, r: 10, b: 50 }
          }}
          useResizeHandler={true}
          style={{ width: "100%" }}
        />
      </div>
    </div>
  );
}

// 6) Summation precision demo
function SummationPrecisionDemo() {
  // generate a list with many small numbers and one large
  const [n, setN] = useState(100000);
  const [small, setSmall] = useState(1e-6);
  const [large, setLarge] = useState(1.0);
  const [results, setResults] = useState(null);

  useEffect(() => {
    // build array: many small numbers + large
    // to keep compute reasonable in browser, clamp n
    const N = Math.max(1000, Math.min(200000, n));
    const arr = new Array(N).fill(small);
    // append large
    arr.push(large);
    // run naive, pairwise, Kahan
    const naive = naiveSum(arr);
    const pair = pairwiseSum(arr);
    const kahan = kahanSum(arr);
    setResults({ naive, pair, kahan, exact: small * N + large, N });
  }, [n, small, large]);

  return (
    <div>
      <div className="flex gap-2 items-center">
        <Input type="number" value={n} onChange={(e) => setN(Number(e.target.value))} className="bg-zinc-700 text-white" />
        <Input value={String(small)} onChange={(e) => setSmall(Number(e.target.value))} className="bg-zinc-700 text-white" />
        <Input value={String(large)} onChange={(e) => setLarge(Number(e.target.value))} className="bg-zinc-700 text-white" />
      </div>

      {results && (
        <div className="mt-3 text-sm text-zinc-300 space-y-2">
          <div>Elements N: <span className="text-zinc-100">{results.N}</span></div>
          <div>Exact sum (ideal): <span className="text-zinc-100">{results.exact}</span></div>
          <div>Naive sum: <span className="text-zinc-100">{results.naive}</span></div>
          <div>Pairwise sum: <span className="text-zinc-100">{results.pair}</span></div>
          <div>Kahan sum: <span className="text-zinc-100">{results.kahan}</span></div>
          <div className="text-zinc-400 text-xs">Observe how naive summation may lose small contributions when a large value is present. Use compensated or pairwise summation when needed.</div>
        </div>
      )}
    </div>
  );
}

// ---------------- Reuse Gauss Elimination already present in your app (replicated here) ----------------
function gaussElimination(A, b) {
  const n = A.length;
  const M = A.map((r) => r.slice());
  const B = b.slice();
  for (let k = 0; k < n; k++) {
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
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    let sum = B[i];
    for (let j = i + 1; j < n; j++) sum -= M[i][j] * x[j];
    x[i] = sum / M[i][i];
  }
  return x;
}

// ---------------- Main Chapter 3 Page ----------------
export default function Chapter3() {
  return (
    <div className={`min-h-screen px-6 py-8 bg-zinc-950 text-zinc-100`}>
      <motion.header initial={{ opacity: 0, y: -10 }} animate={{ opacity: 1, y: 0 }} transition={{ duration: 0.5 }} className="max-w-7xl mx-auto mb-6">
        <div className="flex items-center justify-between">
          <div>
            <h1 className="text-3xl font-bold text-emerald-400">Approximations & Round-Off Errors</h1>
            <p className="text-zinc-300 mt-1 max-w-2xl">Significant figures, accuracy vs precision, error definitions, round-off errors — prioritized documentation, rich visualizations, and hands-on calculations.</p>
          </div>
         
        </div>
      </motion.header>

      <main className="max-w-7xl mx-auto grid grid-cols-1 gap-6">
        {/* top cards summarizing chapter */}
        <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
          <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700`}>
            <CardHeader>
              <CardTitle className="text-cyan-300 flex items-center gap-2"><BookOpen className="w-5 h-5" /> Chapter overview</CardTitle>
            </CardHeader>
            <CardContent>
              <div className="text-zinc-300 text-sm whitespace-pre-line">{docText.overview}</div>
            </CardContent>
          </Card>

          <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700`}>
            <CardHeader>
              <CardTitle className="text-emerald-300 flex items-center gap-2"><Eye className="w-5 h-5" /> Visualization & Practice</CardTitle>
            </CardHeader>
            <CardContent>
              <div className="text-zinc-300 text-sm">This chapter emphasizes practical visualization, hands-on calculators, and experiments to build intuition about numerical accuracy and precision.</div>
              <div className="mt-3 text-zinc-400 text-sm">Use the interactive documentation below to compute examples, run experiments, and visualize error propagation.</div>
            </CardContent>
          </Card>

          <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700`}>
            <CardHeader>
              <CardTitle className="text-amber-300 flex items-center gap-2"><AlertTriangle className="w-5 h-5" /> Key takeaways</CardTitle>
            </CardHeader>
            <CardContent>
              <ul className="list-disc list-inside text-zinc-300 text-sm space-y-1">
                <li>Report uncertainties; don’t overstate precision.</li>
                <li>Reformulate expressions to avoid catastrophic cancellation.</li>
                <li>Use compensated summation for large accumulations.</li>
                <li>Estimate condition numbers for sensitive linear problems.</li>
              </ul>
            </CardContent>
          </Card>
        </div>

        {/* Visualization canvas + interactive content */}
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          {/* left: 3D error bars */}
          <div className="lg:col-span-1">
            <motion.div variants={cardMotion} initial="initial" animate="enter" whileHover="hover">
              <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700`}>
                <CardHeader>
                  <CardTitle className="text-blue-300 flex items-center gap-2"><Layers className="w-5 h-5" /> Error Bars — 3D</CardTitle>
                </CardHeader>
                <CardContent>
                  <p className="text-zinc-200 text-sm mb-2">A small 3D visualization to show distributions of error magnitudes and their visual "weight". Rotate and inspect.</p>
                  <div className="h-56 bg-zinc-950 border border-zinc-800 rounded">
                    <Canvas camera={{ position: [3, 2, 3], fov: 50 }}>
                      <OrbitControls />
                      <ambientLight intensity={0.6} />
                      <ErrorBars3D />
                    </Canvas>
                  </div>
                </CardContent>
              </Card>
            </motion.div>
          </div>

          {/* center: catastrophic cancellation + summation */}
          <div className="lg:col-span-2">
            <motion.div variants={cardMotion} initial="initial" animate="enter" whileHover="hover">
              <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700`}>
                <CardHeader>
                  <CardTitle className="text-emerald-300 flex items-center gap-2"><Sparkles className="w-5 h-5" /> Demos: Cancellation & Summation</CardTitle>
                </CardHeader>
                <CardContent>
                  <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                    <div>
                      <div className="text-zinc-200 font-medium mb-2">Catastrophic Cancellation</div>
                      <CatastrophicCancellationDemo />
                    </div>
                    <div>
                      <div className="text-zinc-200 font-medium mb-2">Summation Precision</div>
                      <SummationPrecisionDemo />
                    </div>
                  </div>
                </CardContent>
              </Card>
            </motion.div>
          </div>
        </div>

        {/* mid: deeper explanations and code-like blocks */}
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          <motion.div variants={cardMotion} initial="initial" animate="enter" whileHover="hover">
            <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700`}>
              <CardHeader>
                <CardTitle className="text-amber-300 flex items-center gap-2"><Terminal className="w-5 h-5" /> Implementation Notes</CardTitle>
              </CardHeader>
              <CardContent>
                <div className="text-zinc-300 text-sm whitespace-pre-line">
                  Practical coding tips:
                  1. Avoid subtracting nearly equal floats — algebraically transform expressions.
                  2. Use compensated summation (Kahan) for long accumulations.
                  3. Use pairwise summation or higher precision in intermediate computations.
                  4. When reporting results from simulations, include standard deviation and sample size.
                </div>
              </CardContent>
            </Card>
          </motion.div>

          <motion.div variants={cardMotion} initial="initial" animate="enter" whileHover="hover">
            <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700`}>
              <CardHeader>
                <CardTitle className="text-blue-300 flex items-center gap-2"><GitMerge className="w-5 h-5" /> Use cases & examples</CardTitle>
              </CardHeader>
              <CardContent>
                <div className="text-zinc-300 text-sm">
                  Examples of where round-off matters:
                  <ul className="list-disc list-inside mt-2 text-zinc-400">
                    <li>Large-scale linear algebra (ill-conditioned matrices)</li>
                    <li>Iterative solvers accumulating error over millions of iterations</li>
                    <li>Signal processing when subtracting baseline from near-zero amplitudes</li>
                    <li>Geometry and computational physics where differences of large coordinates appear</li>
                  </ul>
                </div>
              </CardContent>
            </Card>
          </motion.div>
        </div>

        {/* documentation large integrated into page (not sticky) */}
        <DocumentationPanelLarge startOpen={"3.1 Significant Figures"} />

        {/* concluding problems and "try yourself" section */}
        <motion.div variants={cardMotion} initial="initial" animate="enter" whileHover="hover">
          <Card className={`p-4 ${theme.panelBg} bg-zinc-900/50 border border-zinc-700`}>
            <CardHeader>
              <CardTitle className="text-emerald-300 flex items-center gap-2"><List className="w-5 h-5" /> Problems & Exercises</CardTitle>
            </CardHeader>
            <CardContent>
              <div className="text-zinc-300 text-sm whitespace-pre-line">
                {docText.sections.problems}
              </div>

              <div className="mt-4 grid grid-cols-1 md:grid-cols-2 gap-4">
                <div className="bg-zinc-800 p-3 rounded">
                  <div className="text-zinc-200 font-medium mb-2">Problem A — Sig figs</div>
                  <div className="text-zinc-300 text-sm">
                    Try the SigFig formatter above to answer A.1–A.3. Use scientific notation if necessary.
                  </div>
                </div>

                <div className="bg-zinc-800 p-3 rounded">
                  <div className="text-zinc-200 font-medium mb-2">Problem C — Summation experiment</div>
                  <div className="text-zinc-300 text-sm">
                    Use the Summation Precision demo: set N large and compare naive vs Kahan. Observe differences and write short notes on why they differ.
                  </div>
                </div>
              </div>
            </CardContent>
          </Card>
        </motion.div>

<BottomBar/>
      </main>
    </div>
  );
}
