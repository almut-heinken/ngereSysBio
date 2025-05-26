# Data Directory — iSIRT1_regulatoryModeling

This folder contains all input data used for model construction, simulations, regulatory inference, and analysis throughout the project. Data is organized by theme (models, media, constraints, etc.) and supports all workflows in `/notebook/` and `/code/`.

---

## Folder Overview

### `caco2/`
- Experimental and processed media data for Caco-2 colon cell simulations.
- Files:
  - `DMEM.xlsx`: Uptake flux values from DMEM medium.
  - `caco2_experiment.xlsx`, `caco2_full.xlsx`, `caco2_two_samples.xlsx`: Experimental flux data or subsets for simulation comparison.

---

### `diets/`
- **Diet composition data** used in simulations.
- File:
  - `NewDiets.xlsx`: Defines various dietary media (e.g., Western, Ketogenic) used for FVA/srFBA.

---

### `fva_utils/`
- Utility files supporting FVA selection or filtering.
- File:
  - `1CM_reactions.tsv`: Custom reaction list used in the 1CM analysis.

---

### `infant_microbiomes/`
- Raw and processed **infant microbiome data** from published studies.
- Files:
  - `43856_2024_715_MOESM7_ESM.xlsx`, `43856_2024_715_MOESM15_ESM.xlsx`: Supplementary data from microbiome studies.
  - `microbiome_fluxes_recon3d.xlsx`: Processed microbial butyrate fluxes mapped to Recon3D.

---

### `models/`
Metabolic models in SBML format.

#### `models/sirt1_recon3d/`
- `Recon3D_SIRT1_generic_gpr_fix_HAM_medium_sinks_closed_irreversibles.xml`: Cleaned model with constraints for SIRT1 regulation.
- `iSIRT1.xml`: SIRT1-integrated version of Recon3D.

---

### `nutrition/`
(Empty for now or used to store intermediate nutrition files)

---

### `organ_atlas/`
Tissue-specific models derived from the **Harvey/Harvetta organ atlas**.

- `female/` and `male/`:
  - `Colon.xml`, `Liver.xml`: XML versions of metabolic models for organ-specific simulations.
  - `.mat` files: MATLAB-format versions of the same models.

---

### `pypath/`
Boolean regulatory networks for SIRT1 and supporting TRNs.

- `grouped_sirt1_trn_manually_curated.csv`: Final TRN used in most simulations.
- `regulatory_rules.csv`: Earlier or alternate rule format for Boolean inference.

---

### `reaction_constraints/`
Constraint files for closing reactions and setting irreversibility.

- `closed.xlsx`: List of reactions to close (flux = 0).
- `irreversible.xlsx`: List of reactions treated as irreversible.

---

### `signor/`
- Regulatory information downloaded from the SIGNOR database.
- File:
  - `all_data_03_01_24.tsv`: Full signaling and interaction data for TRN generation.

---

### `subsystems_gene_set/`
Custom gene sets organized by metabolic subsystems.

- Used for pathway enrichment or filtering during flux analysis.
- Files: `glycolysis.xlsx`, `gluconeogenesis.xlsx`, etc.
- Also contains a `README.txt` with original gene set info.

---

### `tissue_models/`
SBML models for **GTEx-derived** tissues, used for context-specific simulations.

- Files: `Recon3D_<Tissue Name>.xml`
  - e.g., `Recon3D_Liver.xml`, `Recon3D_Colon - Sigmoid.xml`
- These models are built using expression data and the FASTCORE algorithm.

---

### `vmh/`
- Supporting data from the **VMH (Virtual Metabolic Human)** database.
- File:
  - `recon-store-genes-1.tsv`: Mapping of Recon3D gene IDs to external sources.

---

## Notes

- **File naming** is consistent across tissues and tools to support reproducibility.
- These files are used throughout `/code/` scripts and `/notebook/` workflows.