(function () {
  "use strict";

  // Display the numeric sentinel as a semantic range endpoint for pre-2000 data.
  function configureObservationYearSlider() {
    var sliderElement = document.getElementById("observation_years");
    if (!sliderElement || !window.jQuery) return;

    var slider = window.jQuery(sliderElement).data("ionRangeSlider");
    if (!slider) {
      window.setTimeout(configureObservationYearSlider, 50);
      return;
    }

    slider.update({
      prettify_enabled: true,
      prettify: function (value) {
        return Number(value) === 1999 ? "<2000" : String(value);
      }
    });
  }

  document.addEventListener("DOMContentLoaded", configureObservationYearSlider);
  if (window.jQuery) {
    window.jQuery(document).on("shiny:connected", configureObservationYearSlider);
  }
  window.setTimeout(configureObservationYearSlider, 250);
})();
