// Electron-based browsers may expose CommonJS globals, causing the UMD bundle
// to export to module.exports instead of window. Normal browsers skip this.
if (
  typeof window.html2canvas !== "function" &&
  typeof module !== "undefined" &&
  typeof module.exports === "function"
) {
  window.html2canvas = module.exports;
}

