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

function calculateWindows() {
  const totalImages = parseInt(document.getElementById("totalImages").value);
  const numWindows = parseInt(document.getElementById("numWindows").value);
  const overlapPercent = parseInt(
    document.getElementById("overlapPercent").value
  );
  const sliceInterval = parseFloat(
    document.getElementById("sliceInterval").value
  );

  const windowSize = optimizeWindowSize(
    totalImages,
    numWindows,
    overlapPercent
  );

  const overlapSize = Math.floor((windowSize * overlapPercent) / 100);
  const stepSize = windowSize - overlapSize;

  const windows = [];
  const centerPoints = [];

  for (let i = 0; i < numWindows; i++) {
    const start = i * stepSize;
    const end = Math.min(start + windowSize - 1, totalImages - 1);

    if (start >= totalImages) break;

    const centerIndex =
      start + (Math.min(windowSize, totalImages - start) - 1) / 2;
    const centercm = centerIndex * sliceInterval;

    windows.push({
      start,
      end,
      center: centerIndex,
      centercm,
    });
    centerPoints.push(centercm);
  }

  return {
    windows,
    centerPoints,
    stepSize,
    totalImages,
    windowSize,
    sliceInterval,
    actualNumWindows: windows.length,
  };
}

function updateVisualization() {
  const data = calculateWindows();
  const {
    windows,
    centerPoints,
    stepSize,
    totalImages,
    windowSize,
    sliceInterval,
    actualNumWindows,
  } = data;

  console.debug("Visualization Data:", data);
  const overlapPercent = document.getElementById("overlapPercent").value;
  // Update value displays
  document.getElementById(
    "totalImagesValue"
  ).textContent = `${totalImages} images`;
  document.getElementById(
    "numWindowsValue"
  ).textContent = `${actualNumWindows} windows`;
  document.getElementById(
    "overlapPercentValue"
  ).textContent = `${overlapPercent}% overlap`;

  // Calculate stats
  const lastWindow = windows[windows.length - 1];
  const firstWindow = windows[0];
  const windowSizeCm = windowSize * sliceInterval;
  const stepSizeCm = stepSize * sliceInterval;
  const seriesLength = totalImages;
  const seriesLengthCm = totalImages * sliceInterval;
  const coverageLength = lastWindow.end - firstWindow.start + 1;
  const coverageLengthCm = (lastWindow.end + 1) * sliceInterval;
  const coverage = (coverageLength / seriesLength) * 100;
  let actualOverlap = 0;
  for (let i = 0; i < windows.length - 1; i++) {
    const currentEnd = windows[i].end;
    const nextStart = windows[i + 1].start;
    actualOverlap += (currentEnd + 1 - nextStart) / windowSize;
  }
  actualOverlap = (actualOverlap / (windows.length - 1)) * 100;
  console.debug("Actual Overlap:", actualOverlap);
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
  ).textContent = `${actualOverlap.toFixed(2)}% (Actual)`;
  document.getElementById(
    "overlapSetInfo"
  ).textContent = `${overlapPercent}% (Set)`;
  document.getElementById("coverageFill").style.width = `${coverage}%`;

  // Create windows visualization
  const windowsContainer = document.getElementById("windowsContainer");
  windowsContainer.innerHTML = "";

  windows.forEach((window, i) => {
    const row = document.createElement("div");
    row.className = "window-row";

    const hue = ((i * 360) / windows.length) % 360;
    const color = `hsl(${hue}, 70%, 55%)`;

    row.innerHTML = `
                        <div class="window-label">Window ${i + 1}</div>
                        <div class="window-timeline">
                            <div class="window-bar" style="
                                left: ${(window.start / totalImages) * 100}%;
                                width: ${
                                  ((window.end - window.start + 1) /
                                    totalImages) *
                                  100
                                }%;
                                background: ${color};
                            ">
                                ${window.start}-${window.end}
                            </div>
                        </div>
                        <div class="window-info">
                            Images: ${window.start} to ${window.end}<br>
                            Center: ${window.centercm.toFixed(2)} cm
                        </div>
                    `;

    windowsContainer.appendChild(row);
  });

  // Create enhanced center points visualization with gradients and overlaps
  const centerTimeline = document.getElementById("centerTimeline");
  centerTimeline.innerHTML = "";

  // Add window gradients first (bottom layer)
  windows.forEach((window, i) => {
    const gradient = document.createElement("div");
    gradient.className = "window-gradient";
    gradient.dataset.windowIndex = i; // Add data attribute for easy reference

    const hue = ((i * 360) / windows.length) % 360;
    const color1 = `hsla(${hue}, 70%, 65%)`;
    const color2 = `hsla(${hue}, 70%, 45%)`;

    gradient.style.left = `${(window.start / totalImages) * 100}%`;
    gradient.style.width = `${
      ((window.end - window.start + 1) / totalImages) * 100
    }%`;
    gradient.style.background = `linear-gradient(to right, ${color1}, ${color2})`;
    gradient.title = `Window ${i + 1}: ${window.start}-${
      window.end
    } (${window.centercm.toFixed(2)} cm center)`;

    centerTimeline.appendChild(gradient);
  });

  // Add overlap indicators (middle layer)
  for (let i = 0; i < windows.length - 1; i++) {
    const currentWindow = windows[i];
    const nextWindow = windows[i + 1];

    // Check if windows overlap
    if (currentWindow.end >= nextWindow.start) {
      const overlapStart = nextWindow.start;
      const overlapEnd = Math.min(currentWindow.end, nextWindow.end);

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

  // Add center points (top layer)
  windows.forEach((window, i) => {
    const point = document.createElement("div");
    point.className = "center-point";
    point.style.left = `${(window.center / totalImages) * 100}%`;
    point.style.background = `hsla(${(i * 360) / windows.length}, 70%, 50%)`;
    point.title = `Window ${i + 1} center: ${window.centercm.toFixed(2)} cm`;

    // Add hover events to highlight corresponding window gradient
    point.addEventListener("mouseenter", () => {
      const gradient = centerTimeline.querySelector(
        `[data-window-index="${i}"]`
      );
      if (gradient) {
        gradient.classList.add("highlighted");
      }
    });

    point.addEventListener("mouseleave", () => {
      const gradient = centerTimeline.querySelector(
        `[data-window-index="${i}"]`
      );
      if (gradient) {
        gradient.classList.remove("highlighted");
      }
    });

    const label = document.createElement("div");
    label.className = "center-label";
    label.textContent = `${window.centercm.toFixed(2)}cm`;
    label.style.left = `${(window.center / totalImages) * 100}%`;

    centerTimeline.appendChild(point);
    centerTimeline.appendChild(label);
  });

  // Center details
  document.getElementById("centerDetails").innerHTML = `
                    <p><strong>Center Points:</strong> ${centerPoints
                      .map((p) => p.toFixed(2) + "cm")
                      .join(", ")}</p>
                    <p><strong>Average Distance:</strong> ${
                      windows.length > 1
                        ? (seriesLengthCm / (windows.length - 1)).toFixed(2) +
                          " cm"
                        : "N/A"
                    }</p>
                `;
}

// Add event listeners
document
  .getElementById("totalImages")
  .addEventListener("input", updateVisualization);
document
  .getElementById("numWindows")
  .addEventListener("input", updateVisualization);
document
  .getElementById("overlapPercent")
  .addEventListener("input", updateVisualization);
document
  .getElementById("sliceInterval")
  .addEventListener("input", updateVisualization);

// Initial visualization
updateVisualization();
