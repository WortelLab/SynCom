# Environmentally Mediated Interactions Predict Community Assembly and Invasion Success in a Gut Microbiota Synthetic Community

This repository contains all code and data files accompanying the paper:

> **"Environmentally mediated interactions predict community assembly and invasion success in a gut microbiota synthetic community"**

The repository provides full reproducibility of the data processing pipeline, model parameterization, numerical simulations, and statistical analyses presented in the article. It is organized into four self-contained modules, each corresponding to a distinct stage of the analysis.

---

## Repository Structure

### 1. `Growthrate_calculations/`
**Language:** R

This module processes raw optical density (OD) measurements collected from monoculture growth experiments. Starting from the raw OD time-series data (`Rawdata/`), the code fits growth curves to extract for  each bacterial species:

- **Final population densities**,
- **Growth rates**,

measured both in standard YCFA medium and in conditioned media derived from other bacterial species.

The conditioned media experiments are central to characterizing environmentally mediated interactions, as they capture how one species byproducts alter the growth of another.

---

### 2. `Parameters/`
**Language:** Python (Jupyter Notebook)

Using the output from `Growthrate_calculations/`, this module estimates all model parameters required for simulating community dynamics:

- **Monoculture parameters** (growth rate, carrying capacity) for each individual species
- **Interaction parameters** quantifying the effect of conditioned media on pairwise species interactions, as shown in **Figures 1B and 1C** of the article

These parameters are used to parametrised the GLV model used throughout the rest of the analysis.

---

### 3. `QPCR_data/`
**Language:** Python (Jupyter Notebook)

This module processes quantitative PCR (qPCR) data from co-culture experiments to determine the final population densities of each species under competitive conditions. The analysis covers:

- **Biculture experiments** (two-species communities)
- **Triculture and quadriculture experiments** (three- and four-species communities)
- **Invasion experiments**, in which a focal species is introduced into an established community

---

### 4. `Main/`
**Language:** Python (Jupyter Notebook)

This is the core analysis module of the repository. It integrates the parameters and experimental data from the previous modules to perform:

- **Analysis of interaction parameters**, including visualization and interpretation of the interaction matrices.
- **Sensitivity analysis** to assess the robustness of model predictions to parameter uncertainty.
- **Simulation of community assembly dynamics**, using a modified version of the generalised LV model.
- **Simulation of invasion dynamics**, predicting whether and under what conditions a random invading species can establish in a resident community.
- **Comparison of simulated and experimental results**.
- **Data analysis of invasion outcomes** across different community compositions and invasion scenarios.

To reduce computational load, pre-computed output files (`.pkl` format) are included alongside the notebooks, allowing users to reproduce figures and analyses without re-running all simulations from scratch.
