import { BrowserRouter as Router, Routes, Route } from "react-router-dom";
import LandingPage from "./pages/LandingPage";
import Chapter1 from "./pages/Chapter1";
import Chapter2 from "./pages/Chapter2";
import Chapter3 from "./pages/Chapter3";
import Chapter4 from "./pages/Chapter4";
import Chapter5 from "./pages/Chapter5";
import Chapter6 from "./pages/Chapter6";
import Chapter7 from "./pages/Chapter7";
import Chapter8 from "./pages/Chapter8";
import Chapter9 from "./pages/Chapter9";
import Chapter10 from "./pages/Chapter10";
import Chapter11 from "./pages/Chapter11";
import Chapter12 from "./pages/Chapter12";
import Chapter13 from "./pages/Chapter13";
import Chapter14 from "./pages/Chapter14";
import Chapter15 from "./pages/Chapter15";
import Chapter16 from "./pages/Chapter16";
import Chapter17 from "./pages/Chapter17";
import Chapter18 from "./pages/Chapter18";
import Chapter19 from "./pages/Chapter19";
import Chapter20 from "./pages/Chapter20";
import Chapter21 from "./pages/Chapter21";
import Chapter22 from "./pages/Chapter22";
import Chapter23 from "./pages/Chapter23";
import Chapter24 from "./pages/Chapter24"
import Chapter25 from "./pages/Chapter25"
import Chapter26 from "./pages/Chapter26"
import Chapter27 from "./pages/Chapter27";
import Chapter28 from "./pages/Chapter28";
import Chapter29 from "./pages/Chapter29";
import Chapter30 from "./pages/Chapter30";
import Chapter31 from "./pages/Chapter31";
import Chapter32 from "./pages/Chapter32";
import ChaptersPage from "./pages/ChaptersPage";
import { Toaster } from "@/components/ui/sonner"
function App() {
  return (
    <Router>
      <Routes>
        {/* Landing page route */}
        <Route path="/" element={<LandingPage />} />
        <Route path="/chapters" element={<ChaptersPage />} />
        <Route path="/chapter1" element={<Chapter1/>}/>
         <Route path="/chapter2" element={<Chapter2/>}/>
         <Route path="/chapter3" element={<Chapter3/>}/>
          <Route path="/chapter4" element={<Chapter4/>}/>
           <Route path="/chapter5" element={<Chapter5/>}/>
           <Route path="/chapter6" element={<Chapter6/>}/>
           <Route path="/chapter7" element={<Chapter7/>}/>
            <Route path="/chapter8" element={<Chapter8/>}/>
            <Route path="/chapter9" element={<Chapter9/>}/>
            ̥<Route path="/chapter10" element={<Chapter10/>}/>
            <Route path="/chapter11" element={<Chapter11/>}/>
             <Route path="/chapter12" element={<Chapter12/>}/>
            <Route path="/chapter13" element={<Chapter13/>}/>  
             <Route path="/chapter14" element={<Chapter14/>}/>  
        <Route path="/chapter15" element={<Chapter15/>}/> 
        <Route path="/chapter16" element={<Chapter16/>}/>
        <Route path="/chapter17" element={<Chapter17/>}/>
         <Route path="/chapter18" element={<Chapter18/>}/>
         <Route path="/chapter19" element={<Chapter19/>}/>
         <Route path="/chapter20" element={<Chapter20/>}/>
          <Route path="/chapter21" element={<Chapter21/>}/>
        <Route path="/chapter22" element={<Chapter22/>}/>
        <Route path="/chapter23" element={<Chapter23/>}/>
        <Route path="/chapter24" element={<Chapter24/>}/>
        <Route path="/chapter25" element={<Chapter25/>}/>
        <Route path="/chapter26" element={<Chapter26/>}/>
        <Route path="/chapter27" element={<Chapter27/>}/>
        <Route path="/chapter28" element={<Chapter28/>}/>
         <Route path="/chapter29" element={<Chapter29/>}/>
         <Route path="/chapter30" element={<Chapter30/>}/>
         <Route path="/chapter31" element={<Chapter31/>}/>
          <Route path="/chapter32" element={<Chapter32/>}/>
         

      </Routes>
<Toaster
  position="top-right"
  toastOptions={{
    classNames: {
      toast: "bg-gradient-to-r from-cyan-500 to-violet-600 text-white shadow-lg rounded-lg",
      title: "font-semibold",
      description: "text-zinc-100",
    },
  }}
/>


    </Router>
  );
}

export default App;
