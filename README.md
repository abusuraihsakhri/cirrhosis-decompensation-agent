# Cirrhosis Decompensation Agent

### [Open the Live Application →](https://abusuraihsakhri.github.io/cirrhosis-decompensation-agent/)

Reference calculators and screening utilities for common cirrhosis and acute-decompensation assessments. The same Python calculation module is available through a command-line interface and runs in the browser through Pyodide.

> **Clinical-use note:** This project is for reference, education, and software validation. It is not an official OPTN calculator and does not establish a diagnosis, transplant status, treatment plan, or TIPS candidacy. Use current local guidance and clinical judgment for patient care.

## Features

- Original MELD, MELD-Na, and adult MELD 3.0 reference calculations
- Child-Turcotte-Pugh score and class
- EASL-CLIF-style ACLF organ-failure staging implemented by this repository
- SBP ascitic PMN threshold check with reference albumin calculation
- HRS-AKI screening from explicitly supplied AKI/exclusion criteria
- TIPS pre-procedure risk and contraindication flags
- Single-case CLI commands and CSV batch processing
- Responsive browser interface with light/dark themes
- In-browser Python execution; no application backend is required

## Browser application

Open the live application, enter the laboratory and clinical values, and select **Analyze**. The Python module is loaded into a Web Worker through Pyodide so calculations run locally without blocking the page.

Patient values entered in the browser are not submitted to this repository or to an application server. The page loads the Pyodide runtime from jsDelivr, so normal network requests for those static runtime assets still occur.

The interface is intended for current desktop and mobile browsers with WebAssembly and Web Worker support. An internet connection is required when the Pyodide runtime is not already cached.

## CLI

Python 3.9 or later is supported.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

Examples:

```bash
cirrhosis-decompensation meld \
  --cr 2.4 --bili 5.8 --inr 2.2 --na 127 --alb 2.4 --female

cirrhosis-decompensation evaluate \
  --cr 1.9 --bili 3.8 --inr 1.8 --na 131 --alb 2.7 \
  --weight 70 --ascites moderate --he 1 --pmn 280

cirrhosis-decompensation sbp --pmn 420 --weight 68

cirrhosis-decompensation hrs \
  --baseline-cr 1.0 --current-cr 2.2 --weight 70 \
  --albumin-no-response --no-shock-nephrotoxins \
  --no-structural-kidney-signs

cirrhosis-decompensation batch -i sample.csv -o results.csv
```

Run `cirrhosis-decompensation --help` or a subcommand with `--help` for all available inputs.

### Batch CSV

Batch processing requires numeric creatinine, bilirubin, INR, sodium, albumin, and weight values (accepted column aliases are defined in `process_batch_csv`). Missing required laboratory values are rejected rather than silently replaced with defaults.

`sample.csv` provides an example input layout.

## Development and verification

Install the development tools and run the same checks used by CI:

```bash
python -m pip install -e ".[dev]"
ruff check .
pytest -q
python -m build
python -m compileall -q cirrhosis_decompensation.py cli.py
node --check docs/app.js
node --check docs/worker.js
```

GitHub Actions tests Python 3.9, 3.10, 3.11, 3.12, and 3.13. The Pages workflow packages the static browser application from `docs/` together with the Python calculation module.

## Technical notes

The core calculation code uses the Python standard library only. The browser layer is plain HTML, CSS, and JavaScript plus Pyodide. Results are rendered with DOM text nodes rather than injected HTML.

MELD-related allocation policy changes over time. For current U.S. allocation calculations and policy context, use the official HRSA/OPTN resources:

- [HRSA/OPTN MELD calculator](https://www.hrsa.gov/optn/data-calculators/allocation-calculators/meld-calculator)
- [North American practice-based TIPS recommendations](https://pmc.ncbi.nlm.nih.gov/articles/PMC8760361/)
- [Pyodide deployment documentation](https://pyodide.org/en/stable/usage/downloading-and-deploying.html)

## License

MIT. See [LICENSE](LICENSE).
