import { useTheme } from "next-themes"
import { Toaster as Sonner } from "sonner";

const Toaster = ({ ...props }) => {
  const { theme = "system" } = useTheme()

  return (
    <Sonner
      theme={theme}
      position="top-right"
      className="toaster group"
      style={{
        // 🌿 Emerald gradient for normal toast
        "--normal-bg": "linear-gradient(90deg, #10b981, #065f46)", 
        "--normal-text": "white",
        "--normal-border": "transparent",

        // ✅ Success (brighter emerald)
        "--success-bg": "linear-gradient(90deg, #34d399, #059669)", 
        "--success-text": "white",

        // ❌ Error (red gradient)
        "--error-bg": "linear-gradient(90deg, #ef4444, #b91c1c)", 
        "--error-text": "white",

        // ℹ️ Info (emerald → teal mix)
        "--info-bg": "linear-gradient(90deg, #14b8a6, #0d9488)",
        "--info-text": "white",
      }}
      {...props}
    />
  );
}

export { Toaster }
