/** @type {import('tailwindcss').Config} */
module.exports = {
  content: ["./src/**/*.{html,ts}"],
  theme: {
    extend: {
      boxShadow: {
        academic: "0 14px 40px rgba(15, 23, 42, 0.08)",
        "academic-lg": "0 24px 64px rgba(15, 23, 42, 0.14)"
      }
    }
  },
  plugins: []
};
