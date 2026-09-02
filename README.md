# Cirrhosis Decompensation Agent

> **Domain:** Gastroenterology, Hepatology & Clinical Nutrition  
> **Reference Guidelines & Standards:** `AASLD & ACG Clinical Practice Guidelines`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

**Cirrhosis Decompensation Agent** is an advanced analytical and computational platform implementing MELD 3.0, SBP Ascitic Analysis & Hepatorenal Syndrome Staging.

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Core Algorithmic & Evaluation Engines

- **`AscitesDegree`** — dedicated module for ascites degree evaluation and state verification.
- **`EncephalopathyGrade`** — dedicated module for encephalopathy grade evaluation and state verification.
- **`ChildPughClass`** — dedicated module for child pugh class evaluation and state verification.
- **`ACLFGrade`** — dedicated module for a c l f grade evaluation and state verification.
- **`MELDResult`** — dedicated module for m e l d result evaluation and state verification.
- **`ChildPughResult`** — dedicated module for child pugh result evaluation and state verification.

---

## 📐 Mathematical Formulation & Logic

```text
  implementing standard international clinical formulas and consensus guidelines:
  meld_score = int(round(meld_raw))
  score = int(round(meld_na_raw))
  score = original_meld
  score = int(round(meld_3_raw))
```

---

## 💻 CLI Quickstart & Usage

### 1. Guided Interactive Mode
```bash
python cli.py
```

### 2. Direct Parameterized Evaluation
```bash
python cli.py --cr <value> --bili <value> --inr <value> --na <value>
```

### Parameter Reference
- `--cr`: Specifies input measurement or parameter value.
- `--bili`: Specifies input measurement or parameter value.
- `--inr`: Specifies input measurement or parameter value.
- `--na`: Specifies input measurement or parameter value.
- `--alb`: Specifies input measurement or parameter value.
- `--weight`: Specifies input measurement or parameter value.
- `--female`: Specifies input measurement or parameter value.
- `--ascites`: Specifies input measurement or parameter value.
- `--he`: Specifies input measurement or parameter value.
- `--pmn`: Specifies input measurement or parameter value.

### Input Data Schema

| Field | Description | Requirement |
|:------|:------------|:------------|
| `case_id` | Parameter / observation metric | Required |
| `patient_id` | Parameter / observation metric | Required |
| `creatinine` | Parameter / observation metric | Required |
| `bilirubin` | Parameter / observation metric | Required |
| `inr` | Parameter / observation metric | Required |
| `sodium` | Parameter / observation metric | Required |
| `albumin` | Parameter / observation metric | Required |
| `weight_kg` | Parameter / observation metric | Required |

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active AST and regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Air-Gapped LLM Reasoning Adapter:** Agnostic integration for local Ollama instances (`llama3`, `mistral`), Claude 3.5 Sonnet, GPT-4o, and deterministic test mocks.
* **Active Learning Bayesian Calibration:** Dynamic tracker updating worker reliability weights and monitoring Brier calibration drift.
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py --tasks 1000 --concurrency 8
```

---

## 🐳 Container Deployment

```bash
docker build -t cirrhosis-decompensation-agent .
docker run -p 8000:8000 cirrhosis-decompensation-agent
```
