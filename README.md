# CosMx
### 📝 Description
Scripts for the analysis of Nanostring CosMx data.
### 🎯 Aim
Perform QC, pre-processing, clustering and visualisation of CosMx data.
### 🔍 Details
- This script uses the output from AtoMx (where CosMx data is initially deposited from the machine) for further analysis.
- A tileDB data object is loaded into the R environment.
- Raw counts, normalised counts and meta data can all be extracted from the tileDB object.
- QC metrics are visualised and saved to file.
- UMAP and initial clusters calculated by AtoMx are extracted and visualised.
  - As the data is not filtered the initial UMAP and clustering results are not accurate.
  - For instance, an excessive number of clusters are identified. We note that the number of clusters is reduced when these cells are removed.
- Cells are visualised in their spatial context - in 2D space on the slide.
  - Cells can be visualised by QC flag or cluster ID.
- An example plot is also made for a random FOV to show how data / gene expression can be visualised on an FOV level.
- Finally, at the end of the script the cells with QC flags are filtered out and a new dataset is created in a Seurat object, using normalised data.
- Seurat functions are then used to process and recluster the data:
  - Find variable features, PCA, cluster, UMAP, identify cluster markers
- After re-clustering, results are further visualised and saved to file. 
