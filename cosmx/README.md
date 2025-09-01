# CosMx

### 📝 Description
Scripts for the analysis of Nanostring CosMx data.

### 🎯 Aim
Perform QC, pre-processing, clustering and visualisation of CosMx data.

### Getting data from AtoMx
- Export Seurat object, tiledb object, flat files, and raw data from AtoMx

### 🔍 Details

#### R Scripts:
- These scripts use the tileDB output from AtoMx.
- Load tileDB data object the R environment.
- Raw counts, normalised counts and meta data can all be extracted from the tileDB object.
- `CosMx_QC.R` 
  - QC metrics are visualised and saved to files
- `CosMx_PreProcess.R`
  - UMAP and initial clusters calculated by AtoMx are extracted and visualised.
  - As the data is not filtered the initial UMAP and clustering results are not accurate.
  - For instance, an excessive number of clusters are identified. We note that the number of clusters is reduced when these cells are removed.
  - Cells are visualised in their spatial context - in 2D space on the slide.
  - Cells can be visualised by QC flag or cluster ID.
  - Cells are scored for cell type signatures - in order to identify hybrid cells
  - All metadata is saved and seurat objects saved as RDS files
- `CosMx_PostProcess.R`
  - RDS files saved in previous script are loaded
  - Add hybrid status to the seurat object meta data nd then remove hybrid cells
  - Re-process data - find variable features, scale, normalise, cluster etc...
  - Use scType for automated cell type annotation
  - Compare EGFR+ vs EGFR- tumour cells

#### squipy scripts
- See squidpy docs: https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_nanostring.html
- This script uses raw data exported from AtoMx to perform analysis in jupyter notebook
- Raw data exported from AtoMx needs to be treated according to documentation here:
  - https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/squidpy-essentials/squidpy-essentials.html#sec-noimage
  - In a single folder for each sample we need the following (copy files and folders as required, see CosMx_Data_Prep.sh for functions needed to copy over files - note this is a guide not a script.):
    - Flat files
    - CellComposite folder
      - Found in RawFiles/<flowcell>/<slide>/CellStatsDir
    - CellLabels folder
      - CellLabels files are in individual FOV folders in the CellStatsDir - need to copy them all into a folder
    - CellOverlay folder
      - Found in RawFiles/<flowcell>/<slide>/CellStatsDir
    - CompartmentLabels folder
      - CellLabels files are in individual FOV folders in the CellStatsDir - need to copy them all into a folder
- Once data is in a single folder in the correct folder format the squidpy script can be used
- Launch an OnDemand session on Apocrita using jupyter notebooks - a lot of memory is required - e..g. 2 cores, 128GB memory
- Install relevant packages needed and then run the notebook


