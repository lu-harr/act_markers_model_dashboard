(function () {
  function notifyError(message) {
    if (window.Shiny) {
      Shiny.setInputValue("screenshot_error", message, { priority: "event" });
    }
  }

  document.addEventListener("click", function (event) {
    var button = event.target.closest("#download_map");
    if (!button) return;

    if (typeof window.html2canvas !== "function") {
      notifyError("The screenshot library has not loaded. Check the network connection and try again.");
      return;
    }

    var map = document.getElementById("map-shell");
    if (!map) {
      notifyError("The map panel could not be found.");
      return;
    }

    button.disabled = true;
    button.classList.add("disabled");

    window.html2canvas(map, {
      useCORS: true,
      allowTaint: false,
      backgroundColor: "#e3d8d6",
      logging: false,
      scale: Math.min(window.devicePixelRatio || 1, 2)
    }).then(function (canvas) {
      var now = new Date();
      var stamp = now.toISOString().replace(/[:.]/g, "-");
      var link = document.createElement("a");
      link.download = "marcse-prediction-map-" + stamp + ".png";
      link.href = canvas.toDataURL("image/png");
      link.click();
    }).catch(function (error) {
      notifyError("Could not capture the map: " + error.message);
    }).finally(function () {
      button.disabled = false;
      button.classList.remove("disabled");
    });
  });
})();
