# iSIRT1_regulatoryModeling

This repository contains all scripts and data needed to reproduce the simulations and analyses in:

> **Roma Pi et al.**  
> *"Genome-scale modeling reveals regulation of human metabolism by the histone deacetylase SIRT1"*

---

## Project Overview

This project explores the regulatory role of **SIRT1**, a histone deacetylase, in controlling human metabolism using genome-scale metabolic modeling. It integrates:

- **Transcriptomics data** from GTEx  
- **Dietary constraints**  
- **Microbiome-derived butyrate fluxes**  
- **Boolean transcriptional regulatory networks (TRNs)**  
- **Flux Balance Analysis (FBA)** and **SRFBA** simulations  

The models are based on **Recon3D**, and context-specific models are constructed using **FASTCORE**.

---

## Repository Structure

- `code`: Python scripts for model building, FVA, SRFBA  
- `notebook`: Jupyter notebooks for running simulations and visualizations  
- `results`: Output files, figures, and simulation results (organized by topic)  
- `data`: Input models, TRNs, gene sets, dietary and microbiome data  

---

## Setup & Dependencies

This project requires Python 3.8+ and standard scientific computing libraries.

Other dependencies (e.g., `pyfastcore`, `libsbml`, `mewpy`) may be required depending on your system setup.

---

## Reproducing the Results

1. **Clone the repository**  
   Run the following in your terminal:  
   bash  
   git clone https://github.com/your-username/iSIRT1_regulatoryModeling.git  
   cd iSIRT1_regulatoryModeling

2. **Run the notebooks** from the `notebook` folder. Key examples include:  
   - `tissues_batch_srfba.ipynb`: SRFBA simulations across GTEx tissue models  
   - `infant_microbiome_colon_srFBA.ipynb`: Simulates butyrate influence on infant colon  
   - `diet_models_robustness_FBA.ipynb`: Diet-specific flux robustness and PCA  
   - `omnipath_sirt1_trn_build.ipynb`: Builds TRN from OmniPath data

3. **Check the results**  
   Outputs will be saved in the `results` folder, organized by analysis domain.

---

## Key Features

- Context-specific model construction using **FASTCORE**  
- TRN integration via Boolean logic  
- SRFBA and standard FVA simulations  
- Robustness and statistical analysis of metabolic fluxes  
- Domain-specific modeling: tissue, diet, microbiome  

---

## Citation

If you use this repository, please cite:

> **Roma Pi et al.**  
> *"Genome-scale modeling reveals regulation of human metabolism by the histone deacetylase SIRT1"*

---

## Contact

**Almut Heinken**  
Université de Lorraine — Systems Biology, NGERE  
almut-heinken@univ-lorraine.fr