// src/components/BottomBar.jsx
import { useNavigate } from "react-router-dom";
import { ArrowLeft } from "lucide-react";

export default function BottomBar() {
  const navigate = useNavigate();

  return (
    <footer className="fixed bottom-4 left-1/2 -translate-x-1/2 z-50 w-[95%] sm:w-[80%] md:w-[65%]">
      <div
        className="flex items-center justify-between px-6 py-3
        rounded-2xl shadow-xl border border-emerald-700/60
        bg-zinc-900/70 backdrop-blur-xl
        transition-all duration-300"
      >
        {/* Website Name */}
        <div className="text-sm sm:text-base font-semibold tracking-tight 
          text-emerald-400 drop-shadow-sm select-none">
          Numexis
        </div>

        {/* Back Button */}
        <button
          onClick={() => navigate("/chapters")}
          className="flex items-center gap-2 px-4 py-1 rounded-xl
            bg-zinc-800/60 backdrop-blur-md border border-zinc-700
            text-zinc-300 font-medium text-sm sm:text-base
            shadow-md hover:bg-zinc-700/70 hover:text-white
            hover:border-zinc-500 active:scale-95
            transition-all duration-200 cursor-pointer"
        >
          <ArrowLeft className="w-4 h-4 sm:w-5 sm:h-5" />
          Back
        </button>
      </div>
    </footer>
  );
}
