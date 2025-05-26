# 📓 Notebooks for SIRT1 Regulatory Modeling

This folder contains all the Jupyter notebooks used for simulations, model building, and data exploration in the **iSIRT1_regulatoryModeling** project. Each notebook focuses on a specific domain of analysis, including tissue-specific FBA, infant microbiome experiments, Caco-2 simulations, and SIRT1 transcriptional network construction.

---

## 📁 Contents

| Notebook | Purpose |
|----------|---------|
| **`Add_DM_reactions.ipynb`** | Prepares Recon3D models by adding or modifying exchange reactions for custom media. |
| **`butyrate_experiments_infant_data-genericR3D_prev.ipynb`** | Legacy version using generic Recon3D to simulate SIRT1 flux changes with butyrate. |
| **`butyrate_infant_colon_models.ipynb`** | Simulates infant colon-specific models under butyrate perturbations with TRN and SRFBA. |
| **`caco2_medium_parser.ipynb`** | Loads and inspects Caco-2 DMEM medium components from Excel files. |
| **`caco2_robustness_analysis.ipynb`** | Performs statistical analysis on FVA/SRFBA results in Caco-2, generating correlation plots and clustermaps. |
| **`diet_models_robustness_FBA.ipynb`** | Simulates various diet-specific models and analyzes robustness of flux profiles using PCA and clustering. |
| **`infant_microbiome_colon_srFBA.ipynb`** | Applies SRFBA to infant colon models with real butyrate/SIRT1 data, correlating fluxes with microbiome profiles. |
| **`omnipath_sirt1_trn_build.ipynb`** | Builds and filters a SIRT1 transcriptional regulatory network from OmniPath. |
| **`tissue_cs_models_fba_analysis.ipynb`** | Performs standard FBA across context-specific tissue models generated from GTEx expression data. |
| **`tissues_batch_srfba.ipynb`** | Applies TRN-constrained SRFBA in batch mode to all tissues using pre-built models. |

---

## 📂 Additional Files

- **Figures** generated from these notebooks (e.g., PCA plots, clustermaps) have been moved to `/results/` under appropriate subfolders (`/results/caco2/`, `/results/infant_microbiome/`, etc.).
- **TRN logs** from `pypath` are stored in `/notebook/pypath_log/`.

---

## 🔄 Recommendations

- To run these notebooks, install the dependencies listed in the main project `README.md`.
- Ensure access to the `/data/` and `/results/` directories for loading models and saving output.