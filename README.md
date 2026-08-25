# Cirrhosis Decompensation & ACLF Clinical Decision Support Engine

> **Domain:** Hepatology, Gastroenterology & Critical Care Liver Transplantation  
> **Clinical Standards & Guidelines:** AASLD, EASL, UNOS / OPTN, and International Club of Ascites (ICA)

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.9%2B-3776AB.svg?logo=python&logoColor=white)
![Tests](https://img.shields.io/badge/Tests-30%2F30%20Passing-brightgreen.svg)
![Coverage](https://img.shields.io/badge/Coverage-100%25-brightgreen.svg)

---

## Overview

The **Cirrhosis Decompensation & Acute-on-Chronic Liver Failure (ACLF) Clinical Decision Support Engine** is a deterministic, clinical-grade analytical system designed for hepatologists, gastroenterologists, critical care intensivists, and liver transplant selection committees.

The platform stratifies hepatic decompensation severity, determines liver transplant allocation prioritization via multiple MELD iterations, stages ACLF according to EASL-CLIF criteria, and generates guideline-adherent management protocols for life-threatening acute decompensating events (SBP, HRS-AKI, Acute Variceal Bleed, Hepatic Encephalopathy, and TIPS candidacy).

---

## Key Clinical Features

- **Comprehensive MELD Suite:**
  - **Original MELD (UNOS 2002):** Logarithmic model incorporating Creatinine, Bilirubin, and INR with dialysis rules and upper/lower physiological caps.
  - **MELD-Na (UNOS 2016):** Sodium-adjusted MELD incorporating dilutional hyponatremia to enhance 90-day waitlist mortality prediction.
  - **MELD 3.0 (OPTN 2023):** State-of-the-art allocation model incorporating serum albumin, female gender bonus (+1.33 points), and non-linear biomarker interaction terms to eliminate historic gender disparity in organ allocation.
- **Child-Turcotte-Pugh (CTP) Engine:** Stratifies liver disease into Classes A, B, and C with 1-year/2-year survival probabilities and perioperative mortality estimates. Includes specific cutoff adjustments for cholestatic liver diseases (PBC / PSC).
- **EASL-CLIF ACLF Staging System:** Evaluates 6 specific organ failure systems (Liver, Kidney, Brain, Coagulation, Circulation, Respiration) to classify patients into ACLF Grades 0 to 3 with corresponding 28-day mortality forecasts and ICU level-of-care triggers.
- **Acute Decompensation Syndromic Engines:**
  - **Spontaneous Bacterial Peritonitis (SBP):** Diagnostic threshold auditing ($\text{PMN} \ge 250/\text{mm}^3$) and Sort et al. IV Albumin dosing protocol ($1.5\text{ g/kg}$ Day 1 within 6h, $1.0\text{ g/kg}$ Day 3).
  - **Hepatorenal Syndrome (HRS-AKI):** ICA 2015 diagnostic validation, KDIGO AKI staging, and Terlipressin continuous infusion ($2-4\text{ mg/day}$) / Norepinephrine titration protocols with IV Albumin volume support.
  - **TIPS Safety Scorer:** Pre-procedural audit evaluating absolute contraindications (severe pulmonary HTN, heart failure, sepsis) and relative contraindications (MELD $> 18$, bilirubin $> 3.0$).

---

## Clinical Formulas & Analytical Models

### 1. Original MELD (UNOS 2002)
$$\text{MELD}_{\text{raw}} = 9.57 \times \ln(\text{Cr}) + 3.78 \times \ln(\text{Bili}) + 11.20 \times \ln(\text{INR}) + 6.43$$
*Bounds: $\text{Cr}, \text{Bili}, \text{INR} \ge 1.0$; $\text{Cr} \le 4.0$ (or $\text{Cr} = 4.0$ if dialyzed $\ge 2\times$ in past 7 days). Rounded to nearest integer, bounded $[6, 40]$.*

### 2. MELD-Na (UNOS 2016)
$$\text{MELD-Na} = \begin{cases} 
\text{MELD} + 1.32 \times (137 - \text{Na}) - [0.033 \times \text{MELD} \times (137 - \text{Na})] & \text{if MELD} > 11 \\ 
\text{MELD} & \text{if MELD} \le 11 
\end{cases}$$
*Bounds: $\text{Na}$ bounded $[125, 137]\text{ mmol/L}$. Final score bounded $[6, 40]$.*

### 3. MELD 3.0 (OPTN 2023)
$$\begin{aligned}
\text{MELD 3.0} = & 1.33 \times [\text{Female}] + 4.56 \times \ln(\text{Bili}) + 0.82 \times (137 - \text{Na}) - 0.24 \times (137 - \text{Na}) \times \ln(\text{Bili}) \\
& + 9.09 \times \ln(\text{INR}) + 11.14 \times \ln(\text{Cr}) + 1.85 \times (3.5 - \text{Alb}) - 1.83 \times (3.5 - \text{Alb}) \times \ln(\text{Cr}) + 6.0
\end{aligned}$$
*Bounds: $\text{Cr} \in [1.0, 3.0]$, $\text{Bili} \ge 1.0$, $\text{INR} \ge 1.0$, $\text{Na} \in [125, 137]$, $\text{Alb} \in [2.0, 3.5]$. Bounded $[6, 40]$.*

---

## Command-Line Interface (CLI)

### 1. Comprehensive Patient Decompensation Evaluation
```bash
python cli.py evaluate --cr 2.2 --bili 4.5 --inr 2.1 --na 129 --alb 2.6 --weight 72 --female --ascites moderate --he 2 --pmn 380
```

### 2. MELD Suite Calculation
```bash
python cli.py meld --cr 2.1 --bili 4.2 --inr 2.0 --na 130 --alb 2.7 --female
```

### 3. Child-Turcotte-Pugh (CTP) Staging
```bash
python cli.py child-pugh --bili 3.2 --alb 2.5 --inr 2.4 --ascites moderate --he 2
```

### 4. EASL-CLIF ACLF Assessment
```bash
python cli.py aclf --cr 2.4 --bili 14.0 --inr 2.6 --he 3 --vasopressors
```

### 5. SBP Diagnostic & Sort Albumin Protocol
```bash
python cli.py sbp --pmn 350 --weight 75
```

### 6. Hepatorenal Syndrome (HRS-AKI) Evaluation
```bash
python cli.py hrs --baseline-cr 1.0 --current-cr 2.4 --weight 70
```

### 7. TIPS Candidacy & Safety Audit
```bash
python cli.py tips --meld 16 --bili 2.2 --cr 1.1 --inr 1.3
```

### 8. Batch CSV Processing
```bash
python cli.py batch --input sample.csv --output results.csv
```

---

## Python API Usage

```python
from cirrhosis_decompensation import CirrhosisDecompensationEngine, AscitesDegree, EncephalopathyGrade

engine = CirrhosisDecompensationEngine()

dossier = engine.evaluate_patient_case(
    serum_creatinine_mg_dl=1.8,
    total_bilirubin_mg_dl=3.2,
    inr=1.7,
    serum_sodium_mmol_l=131.0,
    serum_albumin_g_dl=2.6,
    patient_weight_kg=72.0,
    is_female=True,
    ascites=AscitesDegree.MODERATE_SEVERE_REFRACTORY,
    encephalopathy=EncephalopathyGrade.GRADE_2,
    ascitic_pmn_count=400.0,
    baseline_creatinine_mg_dl=1.0,
)

print(f"MELD-Na:         {dossier.meld_suite.meld_na}")
print(f"Child-Pugh:      Class {dossier.child_pugh.ctp_class.value} ({dossier.child_pugh.total_points} pts)")
print(f"ACLF Status:     {dossier.aclf_status.aclf_grade_label}")
print(f"SBP Confirmed:   {dossier.sbp_protocol.is_sbp_confirmed}")
print(f"Sort Albumin D1: {dossier.sbp_protocol.albumin_dosing_schedule['day_1_grams']} g")
```

---

## Unit Test Suite

Run the comprehensive unit test suite:

```bash
python -m unittest test_cirrhosis_decompensation.py
```

Test coverage includes:
- Original MELD boundary conditions, logarithmic transformations, and dialysis overrides
- MELD-Na hyponatremia scaling and UNOS thresholds
- MELD 3.0 female coefficient integration, albumin interactions, and boundary verification
- Child-Pugh point assignments, class boundaries, and cholestatic disease rules
- EASL-CLIF single and multi-organ failure combinations and mortality scoring
- Sort Albumin SBP infusion calculations
- HRS-AKI KDIGO staging and pharmacotherapy recommendations
- TIPS safety contraindications
- Batch CSV processing pipeline

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
