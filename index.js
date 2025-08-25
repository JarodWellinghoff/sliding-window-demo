function getSlicePositions(totalImages, totalLength, sliceThickness) {
  const spaceBetweenSlices = totalLength / (totalImages - 1);
  const positions = {};
  for (let i = totalImages; i > 0; i--) {
    positions[i] = {
      rel_center_mm: (totalImages - i) * spaceBetweenSlices,
      thickness_mm: sliceThickness,
    };
  }
  return positions;
}

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
  const keys = ["window_coverage", "window_number", "step_size", "percentage"];
  let wsum = 0;
  for (const k of keys) wsum += weights[k] || 0;
  if (wsum === 0) {
    for (const k of keys) weights[k] = 0.25;
    wsum = 1;
  }
  for (const k of keys) weights[k] = (weights[k] || 0) / wsum;

  // Reasonable scales for relative error (prevents one unit dominating)
  const scale = {
    window_coverage: Math.max(
      targets.window_coverage || 5 * sliceInterval,
      sliceInterval
    ),
    window_number: Math.max(targets.window_number || 4, 2),
    step_size: Math.max(targets.step_size || 2 * sliceInterval, sliceInterval),
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
      window_coverage: windowCm,
      window_number: numWindows,
      step_size: stepCm,
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
    targets.window_coverage,
    targets.window_number,
    targets.step_size,
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
    stepSize: stepImages,
    numWindows,
    coverage,
    windowCm,
    stepCm,
    optimizationScore: best.fx,
    converged: best.fx < 0.05,
  };
}

function coverage_mm(start, end = null) {
  if (end === null) return start["thickness_mm"];
  const distance = distance_mm(start, end);
  const thickness = thickness_mm(start, end);
  return distance + thickness;
}

function distance_mm(start, end = null) {
  if (end === null) return 0;
  const start_pos = start["rel_center_mm"];
  const end_pos = end["rel_center_mm"];
  return Math.abs(end_pos - start_pos);
}

function thickness_mm(start, end = null) {
  if (end === null) return start["thickness_mm"];
  return (start["thickness_mm"] + end["thickness_mm"]) / 2;
}

function get_future_positions(positions, start_index) {
  const future_positions_mm = {};
  for (const [key, pos] of Object.entries(positions)) {
    if (key > start_index) {
      future_positions_mm[key] = pos;
    }
  }
  return future_positions_mm;
}

function get_asdf(positions, start_index, target_mm) {
  const future_positions_mm = get_future_positions(positions, start_index);
  const coverages_mm = {};
  for (const [key, pos] of Object.entries(future_positions_mm)) {
    coverages_mm[key] = coverage_mm(positions[start_index], pos);
  }
  const errors = {};
  for (const [key, c] of Object.entries(coverages_mm)) {
    errors[key] = (Math.abs(c - target_mm) / target_mm) * 100;
  }
  const min_error_key = parseInt(
    Object.entries(errors).reduce((min, curr) =>
      curr[1] < min[1] ? curr : min
    )[0]
  );
  const min_error = errors[min_error_key];
  const min_error_coverage = coverages_mm[min_error_key];
  return [min_error_key, min_error, min_error_coverage];
}

function get_asdf1(positions, start_index, target_mm) {
  const future_positions_mm = get_future_positions(positions, start_index);
  const distances_mm = {};
  for (const [key, pos] of Object.entries(future_positions_mm)) {
    distances_mm[key] = distance_mm(positions[start_index], pos);
  }
  const errors = {};
  for (const [key, c] of Object.entries(distances_mm)) {
    errors[key] = (Math.abs(c - target_mm) / target_mm) * 100;
  }
  const min_error_key = parseInt(
    Object.entries(errors).reduce((min, curr) =>
      curr[1] < min[1] ? curr : min
    )[0]
  );
  const min_error = errors[min_error_key];
  const min_error_distance = distances_mm[min_error_key];
  return [min_error_key, min_error, min_error_distance];
}

function calculateWindows() {
  const totalImages = parseInt(
    document.getElementById("totalImages").value,
    10
  );
  const totalLength =
    parseFloat(document.getElementById("totalLength").value) * 10;
  const sliceThickness =
    parseFloat(document.getElementById("sliceThickness").value) * 10;

  let windowSize,
    stepSize,
    numWindows,
    optimizationResult = null;

  // Targets (in cm / units as labeled)

  const slicePositions = getSlicePositions(
    totalImages,
    totalLength,
    sliceThickness
  );

  const smallestKey = Math.min(...Object.keys(slicePositions));
  const largestKey = Math.max(...Object.keys(slicePositions));

  const targets = {
    window_coverage:
      parseFloat(document.getElementById("targetWindowLength").value) * 10,
    step_size: parseFloat(document.getElementById("targetStepSize").value) * 10,
    total_coverage:
      (parseFloat(document.getElementById("targetPercentage").value) / 100) *
      coverage_mm(slicePositions[smallestKey], slicePositions[largestKey]),
  };

  const raw_weights = {
    window_coverage: parseFloat(
      document.getElementById("targetWindowLengthWeight").value
    ),
    step_size: parseFloat(
      document.getElementById("targetStepSizeWeight").value
    ),
    total_coverage: parseFloat(
      document.getElementById("targetPercentageWeight").value
    ),
  };
  let weights;
  const totalWeight = Object.values(raw_weights).reduce(
    (acc, val) => acc + val,
    0
  );

  if (totalWeight === 0) {
    const weight_count = Object.keys(raw_weights).length;
    weights = {};
    for (const [key, val] of Object.entries(raw_weights)) {
      weights[key] = 1 / weight_count;
    }
  } else {
    weights = {};
    for (const [key, val] of Object.entries(raw_weights)) {
      weights[key] = val / totalWeight;
    }
  }
  const pt = [targets["window_coverage"], targets["step_size"]];

  function get_windows(x) {
    let window_start_img = smallestKey;
    let [window_end_img, e1, c] = get_asdf(
      slicePositions,
      window_start_img,
      x[0]
    );
    let [next_window_start_img, e2, d] = get_asdf1(
      slicePositions,
      window_start_img,
      x[1]
    );
    let center_idx = Math.floor((window_start_img + window_end_img) / 2);
    const windows = [
      {
        window: {
          start_img: window_start_img,
          start_idx: window_start_img - 1,
          start_mm: slicePositions[window_start_img].rel_center_mm,
          end_img: window_end_img,
          end_idx: window_end_img - 1,
          end_mm: slicePositions[window_end_img].rel_center_mm,
          error: e1,
          coverage_mm: c,
          coverage_imgs: window_end_img - window_start_img + 1,
          center_idx: center_idx,
          center_mm: slicePositions[center_idx].rel_center_mm,
        },
        step: {
          img: next_window_start_img,
          idx: next_window_start_img - 1,
          mm: slicePositions[next_window_start_img].rel_start_mm,
          error: e2,
          distance_mm: d,
          distance_imgs: next_window_start_img - window_start_img + 1,
        },
      },
    ];
    while (true) {
      window_start_img = next_window_start_img;
      [window_end_img, e1, c] = get_asdf(
        slicePositions,
        window_start_img,
        x[0]
      );
      if (e1 > 1) {
        break;
      }
      [next_window_start_img, e2] = get_asdf1(
        slicePositions,
        window_start_img,
        x[1]
      );
      center_idx = Math.floor((window_start_img + window_end_img) / 2);
      const window = {
        window: {
          start_img: window_start_img,
          start_idx: window_start_img - 1,
          start_mm: slicePositions[window_start_img].rel_center_mm,
          end_img: window_end_img,
          end_idx: window_end_img - 1,
          end_mm: slicePositions[window_end_img].rel_center_mm,
          error: e1,
          coverage_mm: c,
          coverage_imgs: window_end_img - window_start_img + 1,
          center_idx: center_idx,
          center_mm: slicePositions[center_idx].rel_center_mm,
        },
        step: {
          img: next_window_start_img,
          idx: next_window_start_img - 1,
          mm: slicePositions[next_window_start_img].rel_start_mm,
          error: e2,
          distance_mm: d,
          distance_imgs: next_window_start_img - window_start_img + 1,
        },
      };
      windows.push(window);
    }
    return windows;
  }

  function get_stats(windows) {
    const first_window = windows[0];
    const last_window = windows[windows.length - 1];
    const coverage = coverage_mm(
      slicePositions[first_window["window"]["start_img"]],
      slicePositions[last_window["window"]["end_img"]]
    );
    const window_coverage =
      windows.reduce(
        (acc, window) => acc + window["window"]["coverage_mm"],
        0
      ) / windows.length;
    const step_size =
      windows.reduce((acc, window) => acc + window["step"]["distance_mm"], 0) /
      windows.length;
    return {
      total_coverage: coverage,
      window_coverage: window_coverage,
      step_size: step_size,
    };
  }

  function get_error(stats, key) {
    return (Math.abs(stats[key] - targets[key]) / targets[key]) * 100;
  }

  function get_errors(stats) {
    return {
      total_coverage: get_error(stats, "total_coverage"),
      window_coverage: get_error(stats, "window_coverage"),
      step_size: get_error(stats, "step_size"),
    };
  }

  function get_cost(errors, key) {
    return weights[key] * errors[key];
  }

  function get_cost_sum(errors) {
    return Object.keys(errors).reduce(
      (acc, key) => acc + get_cost(errors, key),
      0
    );
  }

  function get_costs(errors) {
    const c = {
      total_coverage: get_cost(errors, "total_coverage"),
      window_coverage: get_cost(errors, "window_coverage"),
      step_size: get_cost(errors, "step_size"),
    };
    c["sum"] = get_cost_sum(errors);
    return c;
  }

  function get_results(x) {
    const _windows = get_windows(x);
    const stats = get_stats(_windows);
    const errors = get_errors(stats);
    const cost = get_costs(errors);
    return {
      windows: _windows,
      stats,
      errors,
      cost,
    };
  }

  function objective(x) {
    const results = get_results(x);
    return results["cost"]["sum"];
  }

  const result = fmin.nelderMead(objective, pt);
  const results = get_results(result.x);
  results["targets"] = targets;
  results["weights"] = weights;
  results["targets"]["totalImages"] = totalImages;

  // Nelder-Mead target optimization
  //   optimizationResult = optimizeWithNelderMead(
  //     totalImages,
  //     targets,
  //     sliceInterval,
  //     { ...weights }
  //   );

  //   const windows = results.windows;
  //   const centerPoints = [];
  //   windowSize = optimizationResult.windowSize;
  //   stepSize = optimizationResult.stepSize;
  //   numWindows = optimizationResult.numWindows;

  //   for (let i = 0; i < numWindows; i++) {
  //     const start = i * stepSize;
  //     const end = Math.min(start + windowSize - 1, totalImages - 1);

  //     if (start >= totalImages) break;

  //     const centerIndex =
  //       start + (Math.min(windowSize, totalImages - start) - 1) / 2;
  //     const centercm = centerIndex * sliceInterval;

  //     _windows.push({ start, end, center: centerIndex, centercm });
  //     centerPoints.push(centercm);
  //   }
  return results;
  //   return {
  //     windows: _windows,
  //     centerPoints,
  //     stepSize,
  //     totalImages,
  //     windowSize,
  //     sliceInterval,
  //     actualNumWindows: _windows.length,
  //     optimizationResult,
  //   };
}

function updateVisualization() {
  const data = calculateWindows();
  const { cost, errors, stats, targets, weights, windows } = data;
  //   const {
  //     windows,
  //     centerPoints,
  //     stepSize,
  //     totalImages,
  //     windowSize,
  //     sliceInterval,
  //     optimizationResult,
  //   } = data;

  console.debug("Visualization Data:", data);
  if (!windows.length) return;

  // Calculate stats
  const lastWindow = windows[windows.length - 1];
  const firstWindow = windows[0];
  const windowSizeCm = stats.window_coverage / 10;
  const stepSizeCm = stats.step_size / 10;
  const seriesLength = targets.totalImages;
  const seriesLengthCm = targets.total_coverage / 10;
  const coverageLength =
    lastWindow.window.end_img - firstWindow.window.start_img + 1;
  const coverageLengthCm = stats.total_coverage / 10;
  const coverage = (coverageLength / seriesLength) * 100;

  let actualOverlap = 0;
  let actualOverlapCm = 0;
  let windowSize = 0;
  let stepSize = 0;
  for (let i = 0; i < windows.length - 1; i++) {
    const currentEndImg = windows[i].window.end_img;
    const currentEndCm = windows[i].window.end_mm / 10;
    const currentStartImg = windows[i].window.start_img;
    const nextStartImg = windows[i + 1].window.start_img;
    const nextStartCm = windows[i + 1].window.start_mm / 10;
    const currentStepSize = windows[i].step.distance_imgs;

    actualOverlap += currentEndImg - nextStartImg;
    actualOverlapCm += currentEndCm - nextStartCm;
    windowSize += currentEndImg - currentStartImg;
    stepSize += currentStepSize;
  }
  actualOverlap = windows.length > 1 ? actualOverlap / (windows.length - 1) : 0;
  actualOverlapCm =
    windows.length > 1 ? actualOverlapCm / (windows.length - 1) : 0;
  windowSize = windows.length > 1 ? windowSize / (windows.length - 1) : 0;
  stepSize = windows.length > 1 ? stepSize / (windows.length - 1) : 0;

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
  document.getElementById("stepSizeInfo").textContent = `${stepSize} images`;
  document.getElementById("stepSizeCmInfo").textContent = `${stepSizeCm.toFixed(
    2
  )} cm`;
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
  ).textContent = `${actualOverlap.toFixed(1)} images`;
  document.getElementById("coverageFill").style.width = `${coverage}%`;
  document.getElementById(
    "overlapSetInfo"
  ).innerHTML = `${actualOverlapCm.toFixed(1)} cm`;

  // Create windows visualization
  const windowsContainer = document.getElementById("windowsContainer");
  windowsContainer.innerHTML = "";

  windows.forEach((w, i) => {
    const row = document.createElement("div");
    const hue = ((i * 360) / windows.length) % 360;
    const color = `hsl(${hue}, 70%, 55%)`;
    const left = (w.window.start_idx / targets.totalImages) * 100;
    const width =
      ((w.window.end_img - w.window.start_img) / targets.totalImages) * 100;
    const center_cm = w.window.center_mm / 10;
    row.className = "card white z-depth-1";
    row.style.padding = "0.75rem";
    row.style.marginBottom = "0.5rem";
    row.innerHTML = `
            <div style="display: flex; align-items: center; gap: 1rem;">
              <div class="window-timeline" style="flex: 1;">
                <div class="window-bar" style="
                  left: ${left}%;
                  width: ${width}%;
                  background: ${color};
                ">${width > 8 ? `W${i + 1}` : ""}
                </div>
              </div>
              <div class="window-info" style="min-width: 140px; text-align: right;">
                <div style="font-weight: 500;">Images: ${
                  w.window.start_img
                } to ${w.window.end_img}</div>
                <div style="color: #666; font-size: 0.8rem;">Center: ${center_cm.toFixed(
                  2
                )} cm</div>
              </div>
            </div>
          `;
    windowsContainer.appendChild(row);
  });

  // Show center timeline section
  document.getElementById("centerTimelineSection").style.display = "block";

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

    gradient.style.left = `${
      (w.window.start_idx / targets.totalImages) * 100
    }%`;
    gradient.style.width = `${
      ((w.window.end_img - w.window.start_img + 1) / targets.totalImages) * 100
    }%`;
    gradient.style.background = `linear-gradient(to right, ${color1}, ${color2})`;
    gradient.title = `Window ${i + 1}: ${w.window.start_img}-${
      w.window.end_img
    } (${(w.window.center_mm * 10).toFixed(2)} cm center)`;

    centerTimeline.appendChild(gradient);
  });

  // Overlap indicators
  for (let i = 0; i < windows.length - 1; i++) {
    const A = windows[i];
    const B = windows[i + 1];

    if (A.window.end_idx >= B.window.start_idx) {
      const overlapStart = B.window.start_idx;
      const overlapEnd = Math.min(A.window.end_idx, B.window.end_idx);

      if (overlapStart <= overlapEnd) {
        const overlap = document.createElement("div");
        overlap.className = "overlap-indicator";
        overlap.style.left = `${(overlapStart / targets.totalImages) * 100}%`;
        overlap.style.width = `${
          ((overlapEnd - overlapStart + 1) / targets.totalImages) * 100
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
    point.style.left = `${(w.window.center_idx / targets.totalImages) * 100}%`;
    point.style.background = `hsla(${(i * 360) / windows.length}, 70%, 50%)`;
    point.title = `Window ${i + 1} center: ${(w.window.center_mm * 10).toFixed(
      2
    )} cm`;

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
    label.textContent = `${(w.window.center_mm * 10).toFixed(2)}cm`;
    label.style.left = `${(w.window.center_idx / targets.totalImages) * 100}%`;

    centerTimeline.appendChild(point);
    centerTimeline.appendChild(label);
  });

  // Optimization info
  //   let optimizationInfo = "";
  //   if (!optimizationResult) {
  //     optimizationInfo = `<p><strong>Nelder-Mead:</strong> <span class="orange-text">No optimization result available</span></p>`;
  //   } else if (optimizationResult.error) {
  //     optimizationInfo = `<p><strong>Nelder-Mead:</strong> <span class="orange-text">${optimizationResult.error}</span></p>`;
  //   } else {
  //     const statusIcon = optimizationResult.converged ? "✓" : "⚠";
  //     const statusColor = optimizationResult.converged
  //       ? "green-text"
  //       : "orange-text";
  //     optimizationInfo = `<p><strong>Nelder-Mead:</strong> Score ${optimizationResult.optimizationScore.toFixed(
  //       4
  //     )} <span class="${statusColor}">${statusIcon}</span></p>`;
  //   }

  //   document.getElementById("centerDetails").innerHTML = `
  //           <h6 class="blue-text text-darken-2"><i class="material-icons left small">info</i>Details</h6>
  //           <p><strong>Center Points:</strong>
  //             ${centerPoints
  //               .map(
  //                 (p, i) =>
  //                   `<span class="chip" style="background: hsla(${
  //                     (i * 360) / centerPoints.length
  //                   }, 70%, 85%)">${p.toFixed(2)}cm</span>`
  //               )
  //               .join("")}
  //           </p>
  //           <p><strong>Average Distance:</strong> ${
  //             windows.length > 1
  //               ? (
  //                   (windows[windows.length - 1].centercm - windows[0].centercm) /
  //                   (windows.length - 1)
  //                 ).toFixed(2) + " cm"
  //               : "N/A"
  //           }</p>
  //           ${optimizationInfo}
  //         `;
}

// Event listeners
document
  .getElementById("totalImages")
  .addEventListener("input", updateVisualization);
document
  .getElementById("totalLength")
  .addEventListener("input", updateVisualization);
document
  .getElementById("sliceThickness")
  .addEventListener("input", updateVisualization);
document
  .getElementById("targetWindowLength")
  .addEventListener("input", updateVisualization);
document
  .getElementById("targetStepSize")
  .addEventListener("input", updateVisualization);
document
  .getElementById("targetPercentage")
  .addEventListener("input", updateVisualization);
document
  .getElementById("targetWindowLengthWeight")
  .addEventListener("input", updateVisualization);
document
  .getElementById("targetStepSizeWeight")
  .addEventListener("input", updateVisualization);
document
  .getElementById("targetPercentageWeight")
  .addEventListener("input", updateVisualization);

// Initialize MaterializeCSS components
document.addEventListener("DOMContentLoaded", function () {
  // Initialize select elements
  var elems = document.querySelectorAll("select");
  var instances = M.FormSelect.init(elems);

  // Initialize application
  updateVisualization();
});
