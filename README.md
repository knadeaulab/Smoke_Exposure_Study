# Single-Cell Exposomics in Human Blood Reveal Aberrant Immune Impacts of Fire Smoke Exposure
## Abstract 
Given the global impact of increasing wildfires associated with climate change, understanding the specific mechanisms of how wildfire smoke affects health is critical to developing preventative and targeted therapies. This study investigates the impacts of fire smoke exposure, including PM2.5, toxic metals, and per- and polyfluoroalkyl substances (PFAS), on the immune systems of firefighters using a single-cell exposomics approach. We collected blood samples from smoke exposed individuals (n=31) and matched non-smoke exposed individuals (n=29), performing immunophenotyping, DNA methylation profiling, and mass cytometry analyses to identify cellular and epigenetic changes associated with exposure. Key findings include increased immune cell subsets coupled to toxic metal isotopes in exposed individuals and several epigenetic modifications associated with metal isotopes across multiple chromosomes. These results highlight potential molecular targets for intervention in smoke exposure-induced immune dysregulation.
## Repository Structure
This repository provides the necessary R scripts, data, and intermediate files for reproducing the results and figures presented in the study. Raw & Processed datasets are available in Harvard Dataverse repository associated with this study. 
### Folder Structure
- **DNAmethylation/**
  - `DNAmethylation.R` - Main analysis R script (annotated). 
- **CyTOF/**
  - `LiveCell.R` - Main analysis R script for analyzing Total live cells.
  - `CD8_Tcells.R` - Main analysis R script for analyzing handgated CD8+ T cells.
  - `CD4_Tcell.R` - Main analysis R script for analyzing handgated CD4+ T cells.
  - `NK_cells.R` - Main analysis R script for analyzing handgated NK cells.
  - `CD11c.R` - Main analysis R script for analyzing handgated CD11c+ Myeloid cells.
  - `B_cell.R` - Main analysis R script for analyzing handgated CD19+ B cells.
- **Files/**
  - Metadata and additional CSV/Tab-delimited files for replication of results presented in this study.
## Instructions
The up-to-date instructions to download associated raw/processed data are available within the each R script.
## Contact
For questions, please contact Abhinav Kaushik at akaushik [at] hsph.harvard.edu

