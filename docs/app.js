"use strict";

const $ = (id) => document.getElementById(id);
const form = $("calculatorForm");
const analyzeBtn = $("analyzeBtn");
const statusEl = $("runtimeStatus");
const errorBox = $("errorBox");
const worker = new Worker("./worker.js");

function numberOrNull(id) {
  const value = $(id).value.trim();
  return value === "" ? null : Number(value);
}

function payloadFromForm() {
  return {
    creatinine: numberOrNull("creatinine"),
    bilirubin: numberOrNull("bilirubin"),
    inr: numberOrNull("inr"),
    sodium: numberOrNull("sodium"),
    albumin: numberOrNull("albumin"),
    weight: numberOrNull("weight"),
    ascites: $("ascites").value,
    he: Number($("he").value),
    pmn: numberOrNull("pmn"),
    baseline_cr: numberOrNull("baselineCr"),
    female: $("female").checked,
    dialysis: $("dialysis").checked,
    vasopressors: $("vasopressors").checked,
    pf_ratio: numberOrNull("pfRatio"),
    pulmonary_hypertension: $("pulmHtn").checked,
    heart_failure: $("heartFailure").checked,
    hrs_albumin_no_response: $("hrsAlbumin").checked,
    hrs_no_shock_nephrotoxins: $("hrsNoShock").checked,
    hrs_no_structural_kidney_signs: $("hrsNoStructural").checked
  };
}

function setText(id, value) {
  $(id).textContent = value;
}

function render(data) {
  const meld = data.meld;
  setText("meld3", meld.meld_3_0 ?? "—");
  setText("meldNa", meld.meld_na);
  setText("meldBand", meld.allocation_tier);
  setText("meldOriginal", `Original MELD ${meld.original_meld}`);

  setText("childClass", `Class ${data.child_pugh.ctp_class}`);
  setText("childPoints", `${data.child_pugh.total_points} points`);

  setText("aclfGrade", `Grade ${data.aclf.aclf_grade}`);
  setText("organFailures", `${data.aclf.organ_failures.total_failures_count} organ failure(s)`);

  if (data.sbp) {
    const schedule = data.sbp.albumin_dosing_schedule;
    setText(
      "sbpResult",
      data.sbp.is_sbp_confirmed
        ? `PMN ≥250/mm³. Reference albumin calculation: ${schedule.day_1_grams} g day 1; ${schedule.day_3_grams} g day 3.`
        : "PMN <250/mm³; SBP threshold not met by this value."
    );
  } else {
    setText("sbpResult", "No ascitic PMN value supplied.");
  }

  if (data.hrs) {
    setText(
      "hrsResult",
      data.hrs.is_hrs_aki_suspected
        ? `Screen positive; KDIGO stage ${data.hrs.kdigo_aki_stage}. Verify full HRS-AKI criteria clinically.`
        : `Screen not positive/incomplete; KDIGO stage ${data.hrs.kdigo_aki_stage} from supplied creatinine values.`
    );
  } else {
    setText("hrsResult", "No baseline creatinine supplied.");
  }

  setText(
    "tipsResult",
    data.tips.is_candidate
      ? `${data.tips.risk_level}. No explicit absolute contraindication supplied; this is not a candidacy determination.`
      : `${data.tips.risk_level}. Specialist assessment is required.`
  );

  const list = $("alertsList");
  list.replaceChildren();
  const alerts = data.alerts.length ? data.alerts : ["No high-acuity flag generated from supplied inputs."];
  for (const alert of alerts) {
    const item = document.createElement("li");
    item.textContent = alert;
    list.appendChild(item);
  }

  setText("resultHint", "Calculated locally with the repository Python module.");
}

function showError(message) {
  errorBox.textContent = message;
  errorBox.hidden = false;
}

worker.addEventListener("message", (event) => {
  const message = event.data;
  if (message.type === "ready") {
    statusEl.textContent = "Python ready";
    statusEl.classList.add("ready");
    analyzeBtn.disabled = false;
    return;
  }
  if (message.type === "result") {
    analyzeBtn.disabled = false;
    analyzeBtn.textContent = "Analyze";
    errorBox.hidden = true;
    render(message.data);
    return;
  }
  if (message.type === "error") {
    analyzeBtn.disabled = false;
    analyzeBtn.textContent = "Analyze";
    showError(message.message || "Calculation failed.");
  }
});

worker.addEventListener("error", () => {
  statusEl.textContent = "Runtime failed";
  showError("Python runtime failed to load. Check your network connection and reload the page.");
});

form.addEventListener("submit", (event) => {
  event.preventDefault();
  errorBox.hidden = true;
  if (!form.reportValidity()) return;
  analyzeBtn.disabled = true;
  analyzeBtn.textContent = "Calculating…";
  worker.postMessage({ type: "calculate", payload: payloadFromForm() });
});

$("exampleBtn").addEventListener("click", () => {
  $("creatinine").value = "2.4";
  $("bilirubin").value = "5.8";
  $("inr").value = "2.2";
  $("sodium").value = "127";
  $("albumin").value = "2.4";
  $("weight").value = "68";
  $("ascites").value = "moderate_severe_refractory";
  $("he").value = "2";
  $("pmn").value = "420";
  $("baselineCr").value = "1.2";
  $("female").checked = true;
});

$("resetBtn").addEventListener("click", () => {
  form.reset();
  errorBox.hidden = true;
  ["meld3", "meldNa", "childClass", "aclfGrade"].forEach((id) => setText(id, "—"));
  setText("meldBand", "Awaiting input");
  setText("meldOriginal", "Original MELD —");
  setText("childPoints", "— points");
  setText("organFailures", "— organ failures");
  setText("sbpResult", "Not calculated.");
  setText("hrsResult", "Not calculated.");
  setText("tipsResult", "Not calculated.");
  $("alertsList").replaceChildren(Object.assign(document.createElement("li"), { textContent: "No result yet." }));
  setText("resultHint", "Run the calculation to populate results.");
});

const savedTheme = localStorage.getItem("theme");
if (savedTheme === "dark") document.documentElement.dataset.theme = "dark";
$("themeToggle").addEventListener("click", () => {
  const next = document.documentElement.dataset.theme === "dark" ? "light" : "dark";
  document.documentElement.dataset.theme = next;
  localStorage.setItem("theme", next);
});
