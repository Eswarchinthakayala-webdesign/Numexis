// src/pages/ChaptersPage.jsx
/* ChaptersPage.jsx
   - Full replacement file
   - Left panel updates right panel (no navigation)
   - Suggestions navigate and select
   - Syncs with URL /chapter{id}
   - Back to top works (scrolls right detail)
   - Subtle 3D animated background using react-three-fiber + drei
   - Tailwind + shadcn/ui + lucide-react + framer-motion
   - Persist notes/favorites/checklists to localStorage
*/

/* ============================
   IMPORTS
   ============================ */
import React, { useEffect, useMemo, useRef, useState } from "react";
import { useNavigate, useLocation } from "react-router-dom";
import CHAPTERS from "../data/chapters"; // path: adjust if needed
import BottomBar2 from "../components/BottomBar2";
// shadcn/ui components - update import paths if your project differs
import { Input } from "@/components/ui/input";
import { Button } from "@/components/ui/button";
import { Card, CardHeader, CardContent } from "@/components/ui/card";
import { Checkbox } from "@/components/ui/checkbox";
import { Separator } from "@/components/ui/separator";
import { Avatar } from "@/components/ui/avatar";

// lucide icons
import {
  BookOpen,
  Search,
  Star,
  Menu,
  X,
  Star as StarOutline,
  CheckCircle2,
  ChevronRight,
  ArrowUp,
  Filter,
  Edit3,
  Play,
  Layers,
  Moon,
  SunMoon,
  Repeat,
  ListRestart,
  RefreshCcw,
  StarHalf,
} from "lucide-react";
import { toast } from "sonner"
// framer-motion
import { motion, AnimatePresence } from "framer-motion";

// three.js react renderer
import { Canvas, useFrame } from "@react-three/fiber";
import { OrbitControls, Points, PointMaterial } from "@react-three/drei";
import * as THREE from "three";

/* ============================
   UTILITIES & SMALL HELPERS
   ============================ */

/* Simple small progress circle used in left and right panels */
function SimpleProgressCircle({ value = 0, size = 64, stroke = 8 }) {
  const r = (size - stroke) / 2;
  const c = Math.PI * 2 * r;
  const dash = (value / 100) * c;

  // Unique gradient id (important when rendering multiple circles)
  const gradientId = `emeraldGradient-${size}-${stroke}`;

  return (
    <svg width={size} height={size} viewBox={`0 0 ${size} ${size}`}>
      {/* Background circle */}
      <circle
        cx={size / 2}
        cy={size / 2}
        r={r}
        strokeWidth={stroke}
        stroke="rgba(255,255,255,0.08)"
        fill="transparent"
      />

      {/* Progress circle */}
      <circle
        cx={size / 2}
        cy={size / 2}
        r={r}
        strokeWidth={stroke}
        stroke={`url(#${gradientId})`}
        fill="transparent"
        strokeLinecap="round"
        strokeDasharray={c}
        strokeDashoffset={c - dash}
        style={{ transition: "stroke-dashoffset 0.6s ease" }}
      />

      {/* Gradient definition */}
      <defs>
        <linearGradient
          id={gradientId}
          x1="0"
          y1="0"
          x2={size}
          y2="0"
          gradientUnits="userSpaceOnUse"
        >
          <stop offset="0%" stopColor="#10b981" /> {/* Emerald */}
          <stop offset="100%" stopColor="#065f46" /> {/* Dark Emerald */}
        </linearGradient>
      </defs>

      {/* Text inside circle */}
      <text
        x="50%"
        y="50%"
        textAnchor="middle"
        dominantBaseline="central"
        fontSize={size * 0.22}
        fill="rgba(255,255,255,0.92)"
        fontWeight={700}
      >
        {value}%
      </text>
    </svg>
  );
}


/* tiny PlayIcon fallback to avoid additional imports */
function PlayIconSmall({ className = "w-4 h-4" }) {
  return (
    <svg className={className} viewBox="0 0 24 24" fill="none">
      <path d="M5 3v18l15-9L5 3z" fill="currentColor" />
    </svg>
  );
}

/* ============================
   3D BACKDROP COMPONENTS
   ============================ */

/* Floating TorusKnot mesh - subtle slow rotation */
function FloatingKnot({ position = [0, 0, 0], color = "#0ea5e9" }) {
  const ref = React.useRef();
  useFrame((state, delta) => {
    if (ref.current) {
      ref.current.rotation.x += 0.01 * delta * 60;
      ref.current.rotation.y += 0.008 * delta * 60;
    }
  });

  // torusKnotGeometry is available in three core and works fine with r3f
  return (
    <mesh ref={ref} position={position} scale={1.1}>
      <torusKnotGeometry args={[1.4, 0.25, 128, 32]} />
      <meshStandardMaterial color={color} metalness={0.9} roughness={0.2} />
    </mesh>
  );
}

/* Particle cloud using Points for a soft background ambience */
function ParticleField({ count = 1200, radius = 8 }) {
  const pointsRef = React.useRef();
  const positions = useMemo(() => {
    const arr = new Float32Array(count * 3);
    for (let i = 0; i < count; i++) {
      const phi = Math.random() * Math.PI * 2;
      const costheta = Math.random() * 2 - 1;
      const u = Math.random();
      const r = radius * Math.cbrt(u);
      const x = r * Math.cos(phi) * Math.sqrt(1 - costheta * costheta);
      const y = r * costheta;
      const z = r * Math.sin(phi) * Math.sqrt(1 - costheta * costheta);
      arr[i * 3 + 0] = x;
      arr[i * 3 + 1] = y;
      arr[i * 3 + 2] = z;
    }
    return arr;
  }, [count, radius]);

  useFrame((state) => {
    if (!pointsRef.current) return;
    pointsRef.current.rotation.y += 0.0008;
    pointsRef.current.rotation.x += 0.0003;
  });

  return (
    <points ref={pointsRef}>
      <bufferGeometry>
        <bufferAttribute attach="attributes-position" count={positions.length / 3} array={positions} itemSize={3} />
      </bufferGeometry>
      <pointsMaterial
        size={0.015}
        sizeAttenuation={true}
        color={"#9ca3af"}
        transparent={true}
        opacity={0.12}
      />
    </points>
  );
}

/* Full Canvas backdrop placed behind the page (absolute) */
function ThreeBackdrop() {
  return (
    <div className="pointer-events-none fixed inset-0 -z-50 opacity-90">
      <Canvas camera={{ position: [0, 0, 8], fov: 45 }} style={{ width: "100%", height: "100%" }}>
        <ambientLight intensity={0.6} />
        <directionalLight position={[5, 5, 5]} intensity={0.4} />
        <FloatingKnot position={[0, -0.8, -1.2]} color={"#0ea5e9"} />
        <FloatingKnot position={[3.2, 1.2, -2.5]} color={"#7c3aed"} />
        <ParticleField count={1000} radius={12} />
        <OrbitControls enableZoom={false} enablePan={false} autoRotate={false} />
      </Canvas>
    </div>
  );
}

/* ============================
   MAIN COMPONENT
   ============================ */

export default function ChaptersPage() {
  const navigate = useNavigate();
  const location = useLocation();

  // persisted state
  const [favorites, setFavorites] = useState(() => {
    try {
      return JSON.parse(localStorage.getItem("fav_chapters")) || [];
    } catch {
      return [];
    }
  });
  const [notesMap, setNotesMap] = useState(() => {
    try {
      return JSON.parse(localStorage.getItem("chapter_notes")) || {};
    } catch {
      return {};
    }
  });
  const [checklists, setChecklists] = useState(() => {
    try {
      return JSON.parse(localStorage.getItem("chapter_checklists")) || {};
    } catch {
      return {};
    }
  });

  // UI state
  const [query, setQuery] = useState("");
  const [filter, setFilter] = useState("All"); // All / Completed / Pending / Favorites
  const [drawerOpen, setDrawerOpen] = useState(false);

  // Selected chapter is controlled by selectedChapterId
  const [selectedChapterId, setSelectedChapterId] = useState(CHAPTERS[0]?.id || 1);

  // Right detail scroll ref
  const detailRef = useRef(null);

  /* Persist to localStorage whenever these change */
  useEffect(() => localStorage.setItem("fav_chapters", JSON.stringify(favorites)), [favorites]);
  useEffect(() => localStorage.setItem("chapter_notes", JSON.stringify(notesMap)), [notesMap]);
  useEffect(() => localStorage.setItem("chapter_checklists", JSON.stringify(checklists)), [checklists]);

  /* Sync selectedChapterId with URL pathname /chapter{id}:
     - When the location changes (back/forward or direct navigation),
       parse the id and set selectedChapterId so right panel updates.
     - When there is no id in the URL, do NOT navigate away — keep current selection.
  */
  useEffect(() => {
    const match = location.pathname.match(/\/chapter(\d+)/);
    if (match) {
      const id = Number(match[1]);
      if (!Number.isNaN(id) && CHAPTERS.some((c) => c.id === id)) {
        setSelectedChapterId(id);
      }
    }
    // If pathname has no chapter, keep previous selection (don't reset to 1)
  }, [location.pathname]);

  /* Ensure checklist exists when a chapter becomes selected */
  useEffect(() => {
    if (!checklists[selectedChapterId]) {
      const ch = CHAPTERS.find((c) => c.id === selectedChapterId);
      if (ch) {
        setChecklists((prev) => ({ ...prev, [selectedChapterId]: ch.subs.map(() => false) }));
      }
    }
  }, [selectedChapterId]);

  /* Suggestions: first 5 entries while typing (non-empty query) */
  const suggestions = useMemo(() => {
    const q = query.trim().toLowerCase();
    if (!q) return [];
    const list = CHAPTERS.filter((c) => {
      if (c.title.toLowerCase().includes(q)) return true;
      if (c.subs.some((s) => s.toLowerCase().includes(q))) return true;
      const notes = (notesMap[c.id] || "").toLowerCase();
      if (notes.includes(q)) return true;
      return false;
    });
    return list.slice(0, 5);
  }, [query, notesMap]);

  /* PROGRESS COMPUTATIONS */
  function chapterProgressPercent(chId) {
    const arr = checklists[chId] || [];
    if (!arr || arr.length === 0) return 0;
    const done = arr.filter(Boolean).length;
    return Math.round((done / arr.length) * 100);
  }

  function overallProgressPercent() {
    if (!CHAPTERS || CHAPTERS.length === 0) return 0;
    const sum = CHAPTERS.reduce((acc, c) => acc + chapterProgressPercent(c.id), 0);
    return Math.round(sum / CHAPTERS.length);
  }

  /* ACTIONS */
  function toggleFavorite(id) {
    setFavorites((prev) => (prev.includes(id) ? prev.filter((x) => x !== id) : [...prev, id]));
  }

  function updateNoteForChapter(id, text) {
    setNotesMap((prev) => ({ ...prev, [id]: text }));
  }

  function toggleChecklistItem(chId, idx) {
    setChecklists((prev) => {
      const arr = prev[chId] ? [...prev[chId]] : CHAPTERS.find((c) => c.id === chId).subs.map(() => false);
      arr[idx] = !arr[idx];
      return { ...prev, [chId]: arr };
    });
  }

  /* When user clicks a chapter item in the LEFT panel:
     - Only select the chapter (update right panel)
     - Don't navigate to /chapter{id}
  */
  function selectChapterOnly(id) {
    setSelectedChapterId(id);
    setDrawerOpen(false);
    // keep URL unchanged
    if (detailRef.current) {
      // scroll details to top on selecting a different chapter
      detailRef.current.scrollTo({ top: 0, behavior: "smooth" });
    }
  }

  /* When suggestions are clicked:
     - Navigate to /chapter{id} (this changes URL)
     - Also select the chapter so right panel updates immediately
  */
  function suggestionClick(id) {
    setQuery("");
    setSelectedChapterId(id);
    navigate(`/chapter${id}`);
    setDrawerOpen(false);
    // right panel will be updated by location effect too
    if (detailRef.current) {
      detailRef.current.scrollTo({ top: 0, behavior: "smooth" });
    }
  }

  /* "Open full chapter page" link click: keep SPA navigation */
  function openFullChapterPage(id) {
    // navigate and also set selection
    setSelectedChapterId(id);
    navigate(`/chapter${id}`);
  }

  /* Back to top scrolls the right detail pane (which is scrollable) */
  function scrollDetailToTop() {
    if (detailRef.current) {
      detailRef.current.scrollTo({ top: 0, behavior: "smooth" });
    } else {
      window.scrollTo({ top: 0, behavior: "smooth" });
    }
  }

  /* Filtered chapters for the sidebar, based on filter buttons */
  const filteredChapters = useMemo(() => {
    return CHAPTERS.filter((c) => {
      if (filter === "All") return true;
      if (filter === "Favorites") return favorites.includes(c.id);
      const perc = chapterProgressPercent(c.id);
      if (filter === "Completed") return perc === 100;
      if (filter === "Pending") return perc < 100;
      return true;
    });
  }, [filter, favorites, checklists]);

  /* Selected chapter object */
  const selectedChapter = useMemo(() => CHAPTERS.find((c) => c.id === selectedChapterId) || CHAPTERS[0], [selectedChapterId]);

  /* small helper to reset checklist for the selected chapter */
  function resetChecklistForSelected() {
    setChecklists((prev) => ({ ...prev, [selectedChapterId]: selectedChapter.subs.map(() => false) }));
  }

  /* UI: small responsive layout values */
  // nothing else

  /* ============================
     RENDER
     ============================ */
  
  return (
    <div className="min-h-screen w-full bg-zinc-950 text-zinc-100 antialiased relative">
      {/* 3D backdrop behind everything */}
      <ThreeBackdrop />

      {/* HEADER */}
      <header className="sticky top-0 z-40 bg-zinc-950/70 backdrop-blur-sm border-b border-zinc-800">
        <div className="max-w-7xl mx-auto px-5 py-4 flex sm:flex-row flex-col items-center gap-4">
          {/* LEFT: Mobile menu + branding */}
          <div className="flex sm:items-center items-start gap-3">
            <button onClick={() => setDrawerOpen(true)} className="md:hidden p-2 rounded-lg cursor-pointer hover:bg-zinc-900 transition" aria-label="Open chapters">
              <Menu className="w-5 h-5" />
            </button>

            <div className="flex  sm:items-center items-start gap-3">
              <div className="w-10 h-10 rounded-lg bg-gradient-to-br from-zinc-900 to-zinc-800 flex items-center justify-center">
                <BookOpen className="w-5 h-5 text-emerald-300" />
              </div>
              <div>
                <div className="text-lg font-extrabold tracking-tight">Numerical Methods — Chapters</div>
            
              </div>
            </div>
          </div>

          {/* CENTER: Search (desktop centered) */}
          <div className="sm:flex-1 w-full mx-6">
            <div className="relative max-w-3xl mx-auto">
              <Search className="absolute left-3 top-3.5 w-4 h-4 text-zinc-400" />
              <Input
                value={query}
                onChange={(e) => setQuery(e.target.value)}
                placeholder="Search chapters, subtopics, notes..."
                className="pl-10 pr-4 py-3 bg-zinc-900/60"
                aria-label="Search chapters"
              />

              {/* Suggestions dropdown (first 5) */}
              <AnimatePresence>
                {suggestions.length > 0 && (
                  <motion.div initial={{ opacity: 0, y: -6 }} animate={{ opacity: 1, y: 0 }} exit={{ opacity: 0, y: -6 }} className="absolute left-0 right-0 mt-2 bg-zinc-900 rounded-lg border border-zinc-800 shadow-lg overflow-hidden">
                    {suggestions.map((s) => (
                      <button
                        key={s.id}
                        onClick={() => suggestionClick(s.id)}
                        className="w-full text-left px-4 py-3 cursor-pointer hover:bg-zinc-800 flex items-center gap-3"
                      >
                        <Avatar>
                          <div className="text-xs flex items-center">{s.id}</div>
                        </Avatar>
                        <div className="flex-1">
                          <div className="font-medium text-sm">{s.title}</div>
                          <div className="text-xs text-zinc-400">{s.subs[0] ?? ""}</div>
                        </div>
                        <div className="text-xs text-zinc-500">{s.subs.length} sub</div>
                      </button>
                    ))}
                  </motion.div>
                )}
              </AnimatePresence>
            </div>
          </div>

          {/* RIGHT: desktop filters */}
          <div className="hidden md:flex items-center gap-3">
            <div className="flex items-center gap-2 rounded-full bg-zinc-900/40 p-1 border border-zinc-800">
              {["All", "Completed", "Pending", "Favorites"].map((f) => (
                <button
                  key={f}
                  onClick={() => setFilter(f)}
                  className={`px-3 py-1 rounded-full text-sm ${filter === f ? "bg-zinc-800 text-white ring-1 cursor-pointer ring-emerald-500" : "text-zinc-400 cursor-pointer hover:bg-zinc-900"}`}
                >
                  {f}
                </button>
              ))}
            </div>
          </div>
        </div>
      </header>

      {/* MAIN GRID */}
      <div className="max-w-7xl mx-auto px-5 py-6 grid grid-cols-12 gap-6">
        {/* LEFT SIDEBAR for desktop */}
        <aside className="hidden md:block col-span-4 lg:col-span-3">
          <Card className="rounded-2xl bg-zinc-900/40 border border-zinc-800 overflow-hidden">
            <CardHeader className="p-6">
              <div className="flex items-center justify-between">
                <div>
                  <div className="text-sm text-zinc-400">Contents</div>
                  <div className="mt-2 text-2xl text-emerald-400 font-bold">Overall Progress</div>
                  <div className="text-sm text-zinc-400">Track your course progress</div>
                </div>
                <div className="flex items-center">
                  <SimpleProgressCircle value={overallProgressPercent()} size={84} />
                </div>
              </div>
            </CardHeader>

            <Separator />

            <CardContent className="p-0">
              <div className="max-h-[56vh] overflow-auto custom-scroll px-3 py-4 space-y-3">
                {/* Chapter list - clicking here only updates right panel (no navigation) */}
                {filteredChapters.map((ch) => (
                  <motion.div
                    layout
                    key={ch.id}
                    onClick={() => selectChapterOnly(ch.id)} // IMPORTANT: no nav here
                    whileHover={{ scale: 1.01 }}
                    className={`flex items-center gap-3 p-3 rounded-lg cursor-pointer ${ch.id === selectedChapterId ? "bg-zinc-800 ring-1 ring-emerald-500" : "hover:bg-zinc-900"}`}
                  >
                    <div className="w-10 h-10 flex items-center justify-center rounded-md bg-zinc-800 text-gray-200 text-sm font-semibold">
                      {ch.id}
                    </div>
                    <div className="flex-1">
                      <div className="font-medium text-gray-200 text-sm">{ch.title}</div>
                      <div className="text-xs text-zinc-400 mt-1">{ch.subs.slice(0, 2).join(" • ")}</div>
                    </div>
                    <div className="flex flex-col items-end gap-1">
                      <div className="text-xs text-zinc-400">{chapterProgressPercent(ch.id)}%</div>
                      <div className="text-xs text-zinc-400">{ch.subs.length} sub</div>
                    </div>
                  </motion.div>
                ))}
              </div>
            </CardContent>
          </Card>
        </aside>

        {/* MOBILE drawer */}
        <AnimatePresence>
          {drawerOpen && (
            <motion.aside initial={{ x: "-100%" }} animate={{ x: 0 }} exit={{ x: "-100%" }} className="fixed inset-y-0 left-0 z-50 w-80 bg-zinc-950 border-r border-zinc-800 p-4 md:hidden">
              <div className="flex items-center justify-between mb-4">
                <div className="flex items-center gap-2">
                  <BookOpen className="w-5 h-5 text-emerald-300" />
                  <div className="font-bold">Chapters</div>
                </div>
                <button onClick={() => setDrawerOpen(false)} className="p-2 rounded-lg hover:bg-zinc-900">
                  <X className="w-5 h-5" />
                </button>
              </div>

              <div className="space-y-3 overflow-auto max-h-[80vh]">
                {filteredChapters.map((ch) => (
                  <button
                    key={ch.id}
                    onClick={() => selectChapterOnly(ch.id)}
                    className={`w-full text-left p-3 rounded-md ${ch.id === selectedChapterId ? "bg-zinc-800" : "hover:bg-zinc-900"} flex items-center gap-3`}
                  >
                    <div className="w-8 h-8 rounded bg-zinc-800 flex items-center justify-center text-sm font-semibold">{ch.id}</div>
                    <div className="flex-1">
                      <div className="font-medium text-sm">{ch.title}</div>
                      <div className="text-xs text-zinc-400">{ch.subs.length} sub</div>
                    </div>
                    <div className="text-xs text-zinc-400">{chapterProgressPercent(ch.id)}%</div>
                  </button>
                ))}
              </div>
            </motion.aside>
          )}
        </AnimatePresence>

        {/* RIGHT: main content / chapter detail */}
        <main className="col-span-12 md:col-span-8 lg:col-span-9">
          <Card className="rounded-2xl bg-zinc-900/40 border border-zinc-800 overflow-hidden">
            <CardHeader className="p-6 flex items-start justify-between gap-6">
              <div className="flex-1">
                <div className="inline-flex items-center gap-2 px-2 py-1 bg-emerald-800/50 border border-emerald-500/50 rounded-full text-xs text-emerald-300">Chapter {selectedChapter?.id}</div>
                <h2 className="mt-3 text-2xl text-emerald-400 font-extrabold leading-tight">{selectedChapter?.title}</h2>
                <p className="mt-2 text-sm text-zinc-400">{selectedChapter?.title} — concise chapter summary goes here.</p>
              </div>

              <div className="flex flex-col items-end gap-4">
                <div className="flex items-center gap-4">
                  <div className="w-24 h-24">
                    <SimpleProgressCircle value={chapterProgressPercent(selectedChapterId)} size={96} />
                  </div>

                  <div className="flex flex-col gap-3 items-end">
                    <button onClick={() => toggleFavorite(selectedChapterId)} className="p-2 rounded-full bg-zinc-900 cursor-pointer hover:bg-zinc-800 transition" aria-label="Toggle favorite">
                      {favorites.includes(selectedChapterId) ? <Star className="w-5 h-5 text-yellow-400" /> : <StarOutline className="w-5 h-5 text-white" />}
                    </button>

                    <button onClick={resetChecklistForSelected} className="px-3 py-2 rounded-md bg-emerald-800 flex items-center cursor-pointer hover:bg-emerald-700 text-sm">
                     <RefreshCcw className="w-4 h-4"/> Reset 
                    </button>

                    <button onClick={() => openFullChapterPage(selectedChapterId)} className="text-xs text-emerald-300 cursor-pointer hover:text-emerald-400 flex items-center gap-1">
                      Open full chapter page <ChevronRight className="w-4 h-4" />
                    </button>
                  </div>
                </div>
              </div>
            </CardHeader>

            <Separator />

            <CardContent className="p-6 grid grid-cols-1 md:grid-cols-3 gap-6">
              {/* LEFT (main): subtopics */}
              <section className="md:col-span-2 space-y-4" ref={detailRef} style={{ maxHeight: "56vh", overflowY: "auto" }}>
                <h3 className="text-lg text-emerald-500 font-semibold">Subtopics</h3>
                <div className="space-y-3">
                  {selectedChapter?.subs.map((s, idx) => {
                    const isChecked = (checklists[selectedChapterId] || [])[idx] || false;
                    return (
                      <motion.div key={idx} layout initial={{ opacity: 0, y: 6 }} animate={{ opacity: 1, y: 0 }} className="flex items-center justify-between gap-3 p-4 bg-zinc-900/30 border border-zinc-800 rounded-lg">
                        <div className="flex items-start gap-3">
                          <Checkbox className="data-[state=checked]:border-emerald-600 data-[state=checked]:bg-emerald-600 data-[state=checked]:text-white dark:data-[state=checked]:border-emerald-700 cursor-pointer dark:data-[state=checked]:bg-emerald-700" checked={isChecked} onCheckedChange={() => toggleChecklistItem(selectedChapterId, idx)} />
                          <div>
                            <div className={`font-medium text-gray-400 ${isChecked ? "line-through text-zinc-300" : ""}`}>{s}</div>
                            <div className="text-xs text-zinc-400">{idx % 2 === 0 ? "10m" : "15m"} • {isChecked ? "Done" : "Pending"}</div>
                          </div>
                        </div>

                        <div className="flex items-center gap-3 text-zinc-400">
                         
                          <button onClick={() => openFullChapterPage(selectedChapterId)} className="p-2 cursor-pointer bg-zinc-700/50 rounded-md hover:bg-zinc-800"><ChevronRight className="w-5 h-5" /></button>
                        </div>
                      </motion.div>
                    );
                  })}
                </div>
              </section>

              {/* RIGHT: notes & quick stats */}
              <aside className="md:col-span-1 space-y-4">
                <div>
                  <h4 className="text-sm font-semibold mb-2 flex text-gray-300 items-center gap-2"><Edit3 className="w-4 h-4 text-emerald-300" /> Notes</h4>
                  <textarea value={notesMap[selectedChapterId] || ""} onChange={(e) => updateNoteForChapter(selectedChapterId, e.target.value)} placeholder="Write private notes for this chapter..." className="w-full h-36 p-3 rounded-lg bg-zinc-900 border border-zinc-800 text-sm text-zinc-100 resize-none" />
                </div>

                <Separator />

                <div className="space-y-3">
                  <div className="flex items-center justify-between">
                    <div className="text-sm text-emerald-500">Chapter Progress</div>
                    <div className="text-sm font-bold">{chapterProgressPercent(selectedChapterId)}%</div>
                  </div>
                  <div className="w-full h-3 bg-zinc-800 rounded-full overflow-hidden">
                    <motion.div className="h-full" style={{ background: "linear-gradient(90deg,#10b921,#059669)" }} initial={{ width: 0 }} animate={{ width: `${chapterProgressPercent(selectedChapterId)}%` }} transition={{ duration: 0.7, ease: "easeInOut" }} />
                  </div>
                </div>

                <div className="space-y-2">
                 
                  <Button onClick={() => toast("Export not implemented. Notes are saved locally in your browser.")} className="w-full cursor-pointer bg-zinc-700 hover:bg-zinc-600">Export Notes</Button>
                </div>
              </aside>
            </CardContent>
          </Card>

     
         <BottomBar2/>
        </main>
      </div>
    </div>
  );
}
