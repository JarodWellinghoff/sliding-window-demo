let mode = "stepSize";

function optimizeWindowSize(totalImages, numWindows, overlapPercent) {
  let bestWindowSize = 1;
  let bestCoverage = 0;

  for (let windowSize = 1; windowSize <= totalImages; windowSize++) {
    const overlapSize = Math.floor((windowSize * overlapPercent) / 100);
    const stepSize = windowSize - overlapSize;

    if (stepSize <= 0) continue;

    const lastWindowStart = (numWindows - 1) * stepSize;
    const lastWindowEnd = lastWindowStart + windowSize - 1;
    const coverage = Math.min(lastWindowEnd + 1, totalImages);
    const coveragePercent = coverage / totalImages;

    if (coveragePercent > bestCoverage && lastWindowEnd < totalImages) {
      bestCoverage = coveragePercent;
      bestWindowSize = windowSize;
    }
  }
  console.debug("Best Window Size:", bestWindowSize);
  return bestWindowSize;
}

// Nelder-Mead optimization function for target matching
function optimizeWithNelderMead(totalImages, targets, sliceInterval, weights) {
  // --- normalize weights to sum=1 (keeping zeros) ---
  const keys = [
    "window_length",
    "window_number",
    "step_or_overlap",
    "percentage",
  ];
  let wsum = 0;
  for (const k of keys) wsum += weights[k] || 0;
  if (wsum === 0) {
    for (const k of keys) weights[k] = 0.25;
    wsum = 1;
  }
  for (const k of keys) weights[k] = (weights[k] || 0) / wsum;

  // Reasonable scales for relative error (prevents one unit dominating)
  const scale = {
    window_length: Math.max(
      targets.window_length || 5 * sliceInterval,
      sliceInterval
    ),
    window_number: Math.max(targets.window_number || 4, 2),
    step_or_overlap: Math.max(
      targets.step_or_overlap || 2 * sliceInterval,
      sliceInterval
    ),
    percentage: Math.max(targets.percentage || 80, 1),
  };

  function clampInt(x, lo, hi) {
    return Math.max(lo, Math.min(hi, Math.floor(x)));
  }

  function objectiveFunction(params) {
    let windowCm = Math.max(params[0], sliceInterval);
    let stepCm = Math.max(params[2], sliceInterval);
    let pctTgt = Math.max(params[3], 0);

    // to indices
    const windowImages = clampInt(windowCm / sliceInterval, 1, totalImages);
    const stepImages = clampInt(
      stepCm / sliceInterval,
      1,
      Math.max(1, windowImages - 1)
    );

    // derived num windows
    const numWindows =
      Math.floor((totalImages - windowImages) / stepImages) + 1;
    if (numWindows < 2) return 1e3; // infeasible

    // coverage
    const lastStart = (numWindows - 1) * stepImages;
    const lastEnd = lastStart + windowImages - 1;
    if (lastEnd >= totalImages) return 1e3;

    const coverage = Math.min(lastEnd + 1, totalImages);
    const coveragePct = (coverage / totalImages) * 100;

    // errors (scaled, mostly relative)
    const v = {
      window_length: windowCm,
      window_number: numWindows,
      step_or_overlap: stepCm,
      percentage: coveragePct,
    };

    let cost = 0;
    for (const k of keys) {
      const t = targets[k] ?? 0;
      const s = scale[k] || 1;
      // smooth relative error (~Huber-like) so big moves are rewarded
      const e = Math.abs(v[k] - t) / s;
      const h = e < 1 ? 0.5 * e * e : e - 0.5; // quadratic near 0, linear far away
      cost += (weights[k] || 0) * h;
    }

    // gentle bound penalties (smooth, avoid cliffs)
    if (stepImages >= windowImages) cost += 5;
    return cost * 100;
  }

  // --- initial guesses (seeded + random) ---
  const seed = [
    targets.window_length,
    targets.window_number,
    targets.step_or_overlap,
    targets.percentage,
  ];

  function randomGuess() {
    const wImg = clampInt(
      (Math.random() * 0.8 + 0.2) * Math.min(totalImages, 200),
      2,
      totalImages
    );
    const sImg = clampInt(Math.random() * (wImg - 1) + 1, 1, wImg - 1);
    const pct = Math.random() * 100;
    return [wImg * sliceInterval, 0, sImg * sliceInterval, pct];
  }

  function jitterGuess(g, factor = 0.6) {
    // multiply by (0.4 .. 1.6) to allow big exploration
    return g.map((x, i) => {
      if (i === 1) return 0; // window_number placeholder, derived
      const m = 1 - factor + Math.random() * (2 * factor);
      return Math.max(sliceInterval, x * m);
    });
  }

  const guesses = [seed, ...Array.from({ length: 8 }, randomGuess)].map((g) =>
    jitterGuess(g)
  );

  let best = null;
  for (const g of guesses) {
    try {
      const res = fmin.nelderMead(objectiveFunction, g);
      if (!best || res.fx < best.fx) best = res;
    } catch (e) {
      // ignore failed starts
    }
  }
  if (!best) throw new Error("All optimization starts failed.");

  const x = best.x;
  const windowCm = Math.max(x[0], sliceInterval);
  const stepCm = Math.max(x[2], sliceInterval);

  const windowImages = Math.max(1, Math.floor(windowCm / sliceInterval));
  const stepImages = Math.max(
    1,
    Math.min(Math.floor(stepCm / sliceInterval), windowImages - 1)
  );
  const numWindows = Math.floor((totalImages - windowImages) / stepImages) + 1;

  const lastStart = (numWindows - 1) * stepImages;
  const lastEnd = lastStart + windowImages - 1;
  const coverage = Math.min(lastEnd + 1, totalImages) / totalImages;

  return {
    windowSize: windowImages,
    stepOrOverlap: stepImages,
    numWindows,
    coverage,
    windowCm,
    stepCm,
    optimizationScore: best.fx,
    converged: best.fx < 0.05,
  };
}

function calculateWindows() {
  const totalImages = parseInt(
    document.getElementById("totalImages").value,
    10
  );
  const totalLength = parseFloat(document.getElementById("totalLength").value);
  const sliceInterval = totalLength / totalImages;

  let windowSize,
    stepOrOverlap,
    numWindows,
    optimizationResult = null;

  // Targets (in cm / units as labeled)
  const targetWindowCm = parseFloat(
    document.getElementById("targetWindowLength").value
  );
  const rawStepOrOverlap = parseFloat(
    document.getElementById("targetStepOrOverlap").value
  );
  const targetPercentage = parseFloat(
    document.getElementById("targetPercentage").value
  );
  const targetWindowNumber = parseFloat(
    document.getElementById("targetWindowNumber").value
  );

  const targetStepCm =
    mode === "overlap"
      ? targetWindowCm * (1 - rawStepOrOverlap / 100) // % → cm
      : rawStepOrOverlap;
  const targets = {
    window_length: targetWindowCm,
    window_number: targetWindowNumber,
    step_or_overlap: targetStepCm,
    percentage: targetPercentage,
  };

  const weights = {
    window_length: parseFloat(
      document.getElementById("targetWindowLengthWeight").value
    ),
    window_number: parseFloat(
      document.getElementById("targetWindowNumberWeight").value
    ),
    step_or_overlap: parseFloat(
      document.getElementById("targetStepOrOverlapWeight").value
    ),
    percentage: parseFloat(
      document.getElementById("targetPercentageWeight").value
    ),
  };

  // Nelder-Mead target optimization
  optimizationResult = optimizeWithNelderMead(
    totalImages,
    targets,
    sliceInterval,
    { ...weights }
  );

  windowSize = optimizationResult.windowSize;
  stepOrOverlap = optimizationResult.stepOrOverlap;
  numWindows = optimizationResult.numWindows;

  const windows = [];
  const centerPoints = [];

  for (let i = 0; i < numWindows; i++) {
    const start = i * stepOrOverlap;
    const end = Math.min(start + windowSize - 1, totalImages - 1);

    if (start >= totalImages) break;

    const centerIndex =
      start + (Math.min(windowSize, totalImages - start) - 1) / 2;
    const centercm = centerIndex * sliceInterval;

    windows.push({ start, end, center: centerIndex, centercm });
    centerPoints.push(centercm);
  }

  return {
    windows,
    centerPoints,
    stepOrOverlap,
    totalImages,
    windowSize,
    sliceInterval,
    actualNumWindows: windows.length,
    optimizationResult,
  };
}

function updateVisualization() {
  const data = calculateWindows();
  const {
    windows,
    centerPoints,
    stepOrOverlap,
    totalImages,
    windowSize,
    sliceInterval,
    optimizationResult,
  } = data;

  console.debug("Visualization Data:", data);
  if (!windows.length) return;

  // Calculate stats
  const lastWindow = windows[windows.length - 1];
  const firstWindow = windows[0];
  const windowSizeCm = windowSize * sliceInterval;
  const stepOrOverlapCm = stepOrOverlap * sliceInterval;
  const seriesLength = totalImages;
  const seriesLengthCm = totalImages * sliceInterval;
  const coverageLength = lastWindow.end - firstWindow.start + 1;
  const coverageLengthCm = coverageLength * sliceInterval; // FIXED
  const coverage = (coverageLength / seriesLength) * 100;

  let actualOverlap = 0;
  for (let i = 0; i < windows.length - 1; i++) {
    const currentEnd = windows[i].end;
    const nextStart = windows[i + 1].start;
    actualOverlap += currentEnd + 1 - nextStart;
  }
  actualOverlap = windows.length > 1 ? actualOverlap / (windows.length - 1) : 0;
  const actualOverlapCm = actualOverlap * sliceInterval;
  // Update stats
  document.getElementById(
    "windowSizeImagesInfo"
  ).textContent = `${windowSize} images`;
  document.getElementById(
    "windowSizeCmInfo"
  ).textContent = `${windowSizeCm.toFixed(2)} cm`;
  document.getElementById("coverageInfo").textContent = `${coverage.toFixed(
    2
  )}%`;
  document.getElementById(
    "stepSizeInfo"
  ).textContent = `${stepOrOverlap} images`;
  document.getElementById(
    "stepSizeCmInfo"
  ).textContent = `${stepOrOverlapCm.toFixed(2)} cm`;
  document.getElementById(
    "seriesLengthInfo"
  ).textContent = `${seriesLength} images`;
  document.getElementById(
    "seriesLengthCmInfo"
  ).textContent = `${seriesLengthCm.toFixed(2)} cm`;
  document.getElementById(
    "coverageLengthInfo"
  ).textContent = `${coverageLength} images`;
  document.getElementById(
    "coverageLengthCmInfo"
  ).textContent = `${coverageLengthCm.toFixed(2)} cm`;
  document.getElementById(
    "overlapActualInfo"
  ).textContent = `${actualOverlap} images`;
  document.getElementById("coverageFill").style.width = `${coverage}%`;
  document.getElementById(
    "overlapSetInfo"
  ).innerHTML = `${actualOverlapCm.toFixed(1)} cm<br>`;

  // Update target comparisons and show optimization details
  if (!optimizationResult) {
    document.getElementById("windowSizeComparison").textContent =
      "Optimization failed";
    document.getElementById("stepSizeComparison").textContent =
      "Using fallback values";
  }

  // Create windows visualization
  const windowsContainer = document.getElementById("windowsContainer");
  windowsContainer.innerHTML = "";

  windows.forEach((w, i) => {
    const row = document.createElement("div");
    row.className = "window-row";

    const hue = ((i * 360) / windows.length) % 360;
    const color = `hsl(${hue}, 70%, 55%)`;

    row.innerHTML = `
      <div class="window-label">Window ${i + 1}</div>
      <div class="window-timeline">
        <div class="window-bar" style="
          left: ${(w.start / totalImages) * 100}%;
          width: ${((w.end - w.start + 1) / totalImages) * 100}%;
          background: ${color};
        ">
          ${w.start + 1}-${w.end + 1}
        </div>
      </div>
      <div class="window-info">
        Images: ${w.start + 1} to ${w.end + 1}<br>
        Center: ${w.centercm.toFixed(2)} cm
      </div>
    `;
    windowsContainer.appendChild(row);
  });

  // Enhanced center timeline
  const centerTimeline = document.getElementById("centerTimeline");
  centerTimeline.innerHTML = "";

  // Gradients
  windows.forEach((w, i) => {
    const gradient = document.createElement("div");
    gradient.className = "window-gradient";
    gradient.dataset.windowIndex = i;

    const hue = ((i * 360) / windows.length) % 360;
    const color1 = `hsla(${hue}, 70%, 65%)`;
    const color2 = `hsla(${hue}, 70%, 45%)`;

    gradient.style.left = `${(w.start / totalImages) * 100}%`;
    gradient.style.width = `${((w.end - w.start + 1) / totalImages) * 100}%`;
    gradient.style.background = `linear-gradient(to right, ${color1}, ${color2})`;
    gradient.title = `Window ${i + 1}: ${w.start}-${
      w.end
    } (${w.centercm.toFixed(2)} cm center)`;

    centerTimeline.appendChild(gradient);
  });

  // Overlap indicators
  for (let i = 0; i < windows.length - 1; i++) {
    const A = windows[i];
    const B = windows[i + 1];

    if (A.end >= B.start) {
      const overlapStart = B.start;
      const overlapEnd = Math.min(A.end, B.end);

      if (overlapStart <= overlapEnd) {
        const overlap = document.createElement("div");
        overlap.className = "overlap-indicator";
        overlap.style.left = `${(overlapStart / totalImages) * 100}%`;
        overlap.style.width = `${
          ((overlapEnd - overlapStart + 1) / totalImages) * 100
        }%`;
        overlap.title = `Overlap between Window ${i + 1} and Window ${
          i + 2
        }: ${overlapStart}-${overlapEnd}`;
        centerTimeline.appendChild(overlap);
      }
    }
  }

  // Center points
  windows.forEach((w, i) => {
    const point = document.createElement("div");
    point.className = "center-point";
    point.style.left = `${(w.center / totalImages) * 100}%`;
    point.style.background = `hsla(${(i * 360) / windows.length}, 70%, 50%)`;
    point.title = `Window ${i + 1} center: ${w.centercm.toFixed(2)} cm`;

    point.addEventListener("mouseenter", () => {
      const g = centerTimeline.querySelector(`[data-window-index="${i}"]`);
      if (g) g.classList.add("highlighted");
    });
    point.addEventListener("mouseleave", () => {
      const g = centerTimeline.querySelector(`[data-window-index="${i}"]`);
      if (g) g.classList.remove("highlighted");
    });

    const label = document.createElement("div");
    label.className = "center-label";
    label.textContent = `${w.centercm.toFixed(2)}cm`;
    label.style.left = `${(w.center / totalImages) * 100}%`;

    centerTimeline.appendChild(point);
    centerTimeline.appendChild(label);
  });

  // Optimization info
  let optimizationInfo = "";
  if (!optimizationResult) {
    optimizationInfo = `<br><strong>Nelder-Mead:</strong> <span style="color: #ff9800;">No optimization result available</span>`;
  } else if (optimizationResult.error) {
    optimizationInfo = `<br><strong>Nelder-Mead:</strong> <span style="color: #ff9800;">${optimizationResult.error}</span>`;
  } else {
    optimizationInfo = `<br><strong>Nelder-Mead:</strong> Score ${optimizationResult.optimizationScore.toFixed(
      4
    )}, ${optimizationResult.converged ? "✓" : "⚠"}`;
  }

  document.getElementById("centerDetails").innerHTML = `
    <p><strong>Center Points:</strong> ${centerPoints
      .map((p) => p.toFixed(2) + "cm")
      .join(", ")}</p>
    <p><strong>Average Distance:</strong> ${
      windows.length > 1
        ? (
            (windows[windows.length - 1].centercm - windows[0].centercm) /
            (windows.length - 1)
          ).toFixed(2) + " cm"
        : "N/A"
    }${optimizationInfo}</p>
  `;
}

// Event listeners
document
  .getElementById("totalImages")
  .addEventListener("input", updateVisualization);
document
  .getElementById("totalLength")
  .addEventListener("input", updateVisualization);
document
  .getElementById("targetWindowNumber")
  .addEventListener("input", updateVisualization);
document
  .getElementById("targetWindowLength")
  .addEventListener("input", updateVisualization);
document
  .getElementById("targetStepOrOverlap")
  .addEventListener("input", updateVisualization);
document
  .getElementById("targetPercentage")
  .addEventListener("input", updateVisualization);
document
  .getElementById("targetWindowNumberWeight")
  .addEventListener("input", updateVisualization);
document
  .getElementById("targetWindowLengthWeight")
  .addEventListener("input", updateVisualization);
document
  .getElementById("targetStepOrOverlapWeight")
  .addEventListener("input", updateVisualization);
document
  .getElementById("targetPercentageWeight")
  .addEventListener("input", updateVisualization);

document.getElementById("targetType").addEventListener("change", (e) => {
  const selectedValue = e.target.value; // "stepSize" | "overlap"
  mode = selectedValue;

  const windowLength = parseFloat(
    document.getElementById("targetWindowLength").value
  );
  const field = document.getElementById("targetStepOrOverlap");

  if (selectedValue === "stepSize") {
    // value is currently % -> convert to cm
    const oldOverlapPercent = parseFloat(field.value);
    const stepSizeCm = windowLength - (oldOverlapPercent / 100) * windowLength;
    field.min = "0.01";
    field.max = "50";
    field.step = "0.01";
    field.value = stepSizeCm.toFixed(3);
  } else {
    // value is currently cm -> convert to %
    const oldStepCm = parseFloat(field.value);
    const overlapPercent = (1 - oldStepCm / windowLength) * 100;
    field.min = "0";
    field.max = "100";
    field.step = "0.1";
    field.value = overlapPercent.toFixed(1);
  }

  updateVisualization();
});

// Initialize
function initializeApplication() {
  updateVisualization();
}
initializeApplication();
