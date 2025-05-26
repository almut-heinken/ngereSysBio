# Flux Analysis Scripts for SIRT1-Related Metabolic Modeling

This repository contains all Python scripts used to generate context-specific models and run flux variability analyses (FVA) related to SIRT1 regulation, diet, and butyrate exposure. The models are based on **Recon3D** and use MEWpy, pyFASTCORE, and custom preprocessing to explore metabolic behavior under various conditions.

## Contents

### 1. `fastcore_GTEx.py`

This script builds **context-specific metabolic models** using the FASTCORE algorithm for multiple GTEx samples based on RNA-seq expression data.

* Loads gene expression data and core reaction definitions.
* Computes gene and reaction confidence scores per sample.
* Runs FASTCORE for each sample in parallel to construct Recon3D-based models.
* Saves each context-specific model in SBML format.

> Output: `Recon3D_<tissue/sample>.xml` files in the GTEx results folder.

---

### 2. `fva_sirt1.py`

Performs **FVA across a range of SIRT1 expression levels** (from 0.0 to 1.0) using a model regulated by a Boolean TRN.

* Uses a standard medium defined in `NewDiets.xlsx`.
* Sets flux bounds for irreversible and exchange reactions.
* Dynamically adjusts TRN based on SIRT1 input.
* Saves FVA results for each expression level.

> Output: `SIRT1_<expression_level>.csv` files under `results/SIRT1/fva_experiments/fva_955_NEW`.

---

### 3. `fva_diets_mp.py`

Runs **FVA simulations across multiple diet-specific models** (e.g., Ketogenic, Western) with variable SIRT1 expression.

* Iterates over predefined diets and SIRT1 levels.
* Reads corresponding diet models and TRN.
* Performs FVA in parallel across diets and SIRT1 conditions.

> Output: `diet_SIRT1_<expression_level>.csv` per diet and condition in `results/SIRT1/fva_experiments/fva_diets_955`.

---

### 4. `caco2_fva_mp.py`

Runs **FVA simulations in Caco-2 cells** under varying **butyrate flux** conditions.

* Reads exchange medium data from DMEM spreadsheet.
* Adjusts the expression of SIRT1 based on a linear regression with butyrate.
* Applies Boolean TRN and performs FVA with MEWpy.

> Output: `butyrate_<flux>_SIRT1_<expression>.csv` in `results/caco2/fva_simulations`.

---

### 5. `tissues_fva_mp.py`

(Similar to `caco2_fva_mp.py`, possibly older version or duplicate)

* Performs FVA in Caco-2 models under different butyrate exposures.
* Adjusts model bounds using DMEM medium data.
* Links butyrate uptake to SIRT1 activity via linear approximation.

> Output: Same as `caco2_fva_mp.py`.

---

## Requirements

* Python ≥ 3.8
* [`mewpy`](https://github.com/BioSystemsUM/mewpy)
* [`pyfastcore`](https://github.com/BioSystemsUM/pyfastcore)
* `cobra`, `pandas`, `openpyxl`, `multiprocessing`, `re`

Install dependencies via pip:

```bash
pip install cobra pandas openpyxl mewpy
```

## Notes

* Adjust paths depending on your local or server directory structure.
* Models are SBML files; regulatory networks are CSV with Boolean rules.
* Run scripts individually; parallelization is included where applicable.