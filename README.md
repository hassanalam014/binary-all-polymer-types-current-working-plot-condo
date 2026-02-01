# binary-all-polymer-types-current-working-plot-condo

This repository contains Python scripts for **thermodynamic modeling, numerical fitting, and visualization of polymer and binary polymer–mixture systems**, with a primary focus on **glass transition behavior (Tg)**, **PVT properties**, and **Condo / CHV / DHV–based theoretical formulations**.

The codebase is research-oriented and intended for use in **polymer physics, chemical engineering, and materials science** studies.

---

## 📖 Description

The project provides tools to:

- Analyze **pure polymers** and **binary polymer mixtures**
- Fit experimental data using:
  - Condo original formulation
  - CHV (Cell–Hole–Volume) model
  - DHV model
- Solve nonlinear thermodynamic equations numerically
- Process experimental datasets (including CO₂–polymer systems)
- Generate 2D and 3D plots for direct comparison with experiments

---
---

## 🧪 Main Capabilities

### Polymer & Mixture Modeling
- Binary mixture property calculation
- Glass transition temperature (Tg) prediction
- Entropy and thermodynamic variable analysis

### Numerical Methods
- Nonlinear equation solving (`nsolve`)
- Bisection methods for stable convergence
- Discontinuity detection in thermodynamic curves

### Data Processing
- Load experimental polymer and CO₂ data
- Automatic splitting into:
  - Isotherms
  - Regions above and below Tg
- Interpolation of sparse datasets

### Visualization
- 2D comparison plots (model vs experiment)
- 3D thermodynamic surface plots
- Publication-ready figures

---

## 🧠 Core Files Explained

| File | Purpose |
|----|----|
| `All_Functions.py` | Shared mathematical and thermodynamic utilities |
| `Parameters_of_Different_Polymers.py` | Polymer-specific constants |
| `Parameters_for_Mixtures_and_Tg.py` | Mixture and Tg parameters |
| `Tait_Parameters_of_Different_Polymers.py` | PVT modeling using Tait equation |
| `calculateBinaryVariablesCHV_Condo_nsolve.py` | Binary mixture solver |
| `fit_mixtureCondo_Original_nsolve.py` | Condo model fitting |
| `fit_mixtureDHV_bisect_method.py` | DHV fitting via bisection |
| `loadExperimentalData.py` | Experimental data loader |
| `plot_mixture.py` | Mixture property visualization |

---

## ⚙️ Requirements

Python 3.x with the following libraries:

```bash
pip install numpy scipy matplotlib pandas openpyxl
```
