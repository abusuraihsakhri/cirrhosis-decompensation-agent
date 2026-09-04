# Cirrhosis Decompensation, MELD Suite & ACLF Clinical Decision Support Engine

> **Domain:** Hepatology, Critical Care Gastroenterology & Liver Transplantation  
> **Clinical Guidelines & Standards:** AASLD 2021 Practice Guidance on Prevention and Management of Cirrhosis Complications, EASL Clinical Practice Guidelines on Decompensated Cirrhosis (2018), OPTN MELD 3.0 Policy, EASL-CLIF Consortium ACLF Definitions, International Club of Ascites (ICA) HRS-AKI Diagnostic Criteria

---

## 📖 Clinical Overview

The **Cirrhosis Decompensation Agent** provides multi-dimensional clinical decision support for patients experiencing acute decompensation of cirrhosis. It computes the entire MELD family (Original MELD, MELD-Na, and MELD 3.0 with female sex coefficient and serum albumin integration), grades Child-Turcotte-Pugh status, diagnoses and protocols Spontaneous Bacterial Peritonitis (SBP ascitic PMN $\ge 250/\mu\text{L}$ with Sort Albumin dosing), stages Hepatorenal Syndrome (HRS-AKI with Terlipressin/Albumin protocols), audits EASL-CLIF Acute-on-Chronic Liver Failure (ACLF Grades 1–3), and screens transjugular intrahepatic portosystemic shunt (TIPS) candidacy.

### Key Clinical Protocols & Formulas

#### 1. MELD 3.0 Formula (OPTN/UNOS 2023 Implementation)
$$\begin{aligned}
\text{MELD 3.0} = & 1.33 \times (\text{Female}) + [4.56 \times \ln(\text{Bilirubin})] + [0.82 \times (137 - \text{Na})] - [0.24 \times (137 - \text{Na}) \times \ln(\text{Bilirubin})] \\
& + [9.09 \times \ln(\text{INR})] + [11.14 \times \ln(\text{Creatinine})] + [1.85 \times (3.5 - \text{Albumin})] - [1.83 \times (3.5 - \text{Albumin}) \times \ln(\text{Creatinine})] + 6.64
\end{aligned}$$
*(Upper limit capped at 40; variables bound to physiological boundaries).*

#### 2. Acute Complication Management Protocols
- **Spontaneous Bacterial Peritonitis (SBP):** Ascitic PMN $\ge 250/\mu\text{L} \implies$ 3rd-gen cephalosporin (Cefotaxime 2g IV q8h) + IV Albumin ($1.5\,\text{g/kg}$ within 6 hours, $1.0\,\text{g/kg}$ on Day 3).
- **HRS-AKI Protocol:** ICA-AKI criteria met without response to 48h diuretic cessation and albumin volume expansion $\implies$ Terlipressin (1 mg IV bolus q4-6h) + Albumin ($20 - 40\,\text{g/day}$).
- **EASL-CLIF ACLF Classification:** Evaluates failure across 6 organ systems (Liver, Kidney, Brain, Coagulation, Circulation, Respiration) grading ACLF 1 through 3 with ICU allocation alerts.

---

## 💻 CLI Quickstart & Usage

### 1. Comprehensive Decompensation Audit
```bash
python cli.py evaluate --cr 1.9 --bili 3.8 --inr 1.8 --na 131.0 --albumin 2.7 --female --ascites moderate --he 1 --pmn 280
```

### 2. Isolated MELD Suite Calculation
```bash
python cli.py meld --cr 2.4 --bili 5.8 --inr 2.2 --na 127.0 --albumin 2.4 --female
```

### 3. SBP Paracentesis Evaluation
```bash
python cli.py sbp --pmn 420 --weight 68.0 --creatinine 2.4 --bilirubin 5.8
```

### 4. Batch Process Cirrhosis Cohort CSV
```bash
python cli.py batch -i sample.csv -o out_results.csv
```

---

## 🧪 Verification & Testing

Execute comprehensive unit tests via pytest:
```bash
python -m pytest -p no:zarr
```
