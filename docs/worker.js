"use strict";

const PYODIDE_BASE = new URL("./pyodide/", self.location.href).href;
let pyodide;

async function initialize() {
  importScripts(new URL("pyodide.js", PYODIDE_BASE).href);
  pyodide = await loadPyodide({ indexURL: PYODIDE_BASE });

  const response = await fetch("./cirrhosis_decompensation.py", { cache: "no-store" });
  if (!response.ok) {
    throw new Error("Unable to load calculator module.");
  }

  pyodide.runPython(await response.text());
  self.postMessage({ type: "ready" });
}

const ready = initialize().catch((error) => {
  const detail = String(error?.message || error);
  self.postMessage({
    type: "error",
    message: "Bundled Python runtime failed to load. Reload the page after the latest deployment. " + detail
  });
  throw error;
});

self.addEventListener("message", async (event) => {
  if (event.data?.type !== "calculate") return;

  try {
    await ready;
    pyodide.globals.set("payload_json", JSON.stringify(event.data.payload));

    const resultJson = pyodide.runPython(`
import json
from dataclasses import asdict

p = json.loads(payload_json)
engine = CirrhosisDecompensationEngine()
dossier = engine.evaluate_patient_case(
    serum_creatinine_mg_dl=p["creatinine"],
    total_bilirubin_mg_dl=p["bilirubin"],
    inr=p["inr"],
    serum_sodium_mmol_l=p["sodium"],
    serum_albumin_g_dl=p["albumin"],
    patient_weight_kg=p["weight"],
    is_female=p["female"],
    on_dialysis=p["dialysis"],
    ascites=AscitesDegree(p["ascites"]),
    encephalopathy=EncephalopathyGrade(p["he"]),
    ascitic_pmn_count=p["pmn"],
    baseline_creatinine_mg_dl=p["baseline_cr"],
    requires_vasopressors=p["vasopressors"],
    pao2_fio2_ratio=p["pf_ratio"],
    has_severe_pulm_htn=p["pulmonary_hypertension"],
    has_severe_heart_failure=p["heart_failure"],
    hrs_no_response_to_albumin=p["hrs_albumin_no_response"],
    hrs_no_shock_or_nephrotoxins=p["hrs_no_shock_nephrotoxins"],
    hrs_no_structural_kidney_signs=p["hrs_no_structural_kidney_signs"],
)

out = {
    "meld": asdict(dossier.meld_suite),
    "child_pugh": asdict(dossier.child_pugh),
    "aclf": asdict(dossier.aclf_status),
    "sbp": asdict(dossier.sbp_protocol) if dossier.sbp_protocol else None,
    "hrs": asdict(dossier.hrs_aki_protocol) if dossier.hrs_aki_protocol else None,
    "tips": asdict(dossier.tips_eligibility),
    "alerts": dossier.clinical_alerts,
}

json.dumps(out, default=lambda x: x.value if hasattr(x, "value") else str(x))
`);

    self.postMessage({ type: "result", data: JSON.parse(resultJson) });
  } catch (error) {
    self.postMessage({ type: "error", message: String(error?.message || error) });
  }
});
