# Proteomics Analysis Pipeline

## Overview
This repository contains an analysis pipeline for processing and analyzing proteomics data. The workflow includes data preprocessing, normalization, PCA analysis, outlier detection, differential expression analysis, and visualization of significant results. The analysis is tailored for comparing different groups, such as controls and patients at various time points.

---

## Features
1. **Data Preprocessing**:
   - Handles missing values.
   - Converts raw proteomics counts to integers for downstream analysis.
   - Ensures alignment between metadata and dataset.
2. **Normalization**:
   - Implements Total Ion Current (TIC) normalization.
   - Filters out proteins with high proportions of zero values.
   - Log2 normalization for stabilization.
3. **Principal Component Analysis (PCA)**:
   - Performs PCA for both pre- and post-outlier removal datasets.
   - Generates interactive 3D PCA plots.
4. **Outlier Detection and Removal**:
   - Identifies and removes outlier samples to improve data quality.
5. **Differential Expression Analysis**:
   - Utilizes `limma` and `voom` for linear modeling and statistical testing.
   - Defines contrasts for multiple group comparisons.
6. **Visualization**:
   - Generates PCA plots, volcano plots, boxplots, and heatmaps.
   - Creates Venn diagrams to compare significant results across groups.
7. **Output Files**:
   - Exports significant results, unique proteins, and differential expression results.

---

## Data Requirements
- **Proteomics Data**:
  - CSV file containing raw proteomics counts. Each row represents a protein, and each column represents a sample.
- **Metadata**:
  - CSV file containing sample annotations, such as group labels and time points.

---

## Workflow
### 1. Input Data Preparation
- **Input Files**:
  - `../Raw_Counts_serum_proteomics.csv`: Raw proteomics data.
  - `../Metadata_TimePoints.csv`: Metadata for the samples.
- Preprocess the raw counts by handling missing values, reshaping, and aligning metadata with the dataset.

### 2. Normalization
- Apply Total Ion Current (TIC) normalization to scale sample counts.
- Filter out proteins with more than 80% zero values.
- Apply log2 normalization for downstream analysis.

### 3. Outlier Detection and Removal
- Identify outlier samples using PCA and remove them.
- Reapply normalization steps to the filtered dataset.

### 4. Differential Expression Analysis
- Define contrasts for group comparisons:
  - Control vs Patient groups at different time points.
- Use `voom` and `limma` to identify differentially expressed proteins.
- Extract results for each contrast, including p-values, adjusted p-values, and log fold changes.

### 5. Visualization
- **PCA Plots**:
  - Generate pre- and post-normalization PCA plots.
  - Create interactive 3D PCA plots using `plotly`.
- **Volcano Plots**:
  - Highlight significant up- and down-regulated proteins.
- **Heatmaps**:
  - Visualize significant protein expression across groups.
- **Boxplots**:
  - Show expression levels of significant proteins.
- **Venn Diagrams**:
  - Compare significant results across multiple contrasts.

### 6. Output Files
- Differentially expressed proteins for each comparison:
  - `../Results_Files/Differential_Proteins_<Comparison>.csv`
- Unique proteins for each comparison:
  - `../Results_Files/Unique_Proteins_<Comparison>.csv`
- Venn diagrams and other plots:
  - Saved as PDFs in the `../Plots/` directory.

---

## Repository Structure
```
.
├── Scripts
│   └── analysis_pipeline.R  # Main analysis script
├── Data
│   ├── Raw_Counts_serum_proteomics.csv
│   └── Metadata_TimePoints.csv
├── Results_Files
│   ├── Differential_Proteins_<Comparison>.csv
│   └── Unique_Proteins_<Comparison>.csv
├── Plots
│   ├── PCA_Plots.pdf
│   ├── Volcano_Plots.pdf
│   └── Venn_TopGenes_ggvenn.pdf
└── README.md
```

---

## Dependencies
The pipeline requires the following R libraries:
- `readxl`
- `dplyr`
- `limma`
- `ggplot2`
- `tidyr`
- `pheatmap`
- `RColorBrewer`
- `bit64`
- `ggrepel`
- `plotly`
- `data.table`
- `htmlwidgets`
- `tidyverse`
- `edgeR`
- `ggvenn`
- `ComplexHeatmap`
- `circlize`

Install dependencies using:
```R
install.packages(c("readxl", "dplyr", "ggplot2", "tidyr", "pheatmap", "RColorBrewer", "bit64", "ggrepel", "plotly", "data.table", "htmlwidgets", "tidyverse", "edgeR"))
BiocManager::install(c("limma", "ComplexHeatmap", "circlize", "ggvenn"))
```

---

## Usage
1. Clone the repository:
   ```bash
   git clone <repository_url>
   ```

2. Set the working directory to the location of the scripts and data files.

3. Run the analysis pipeline in R:
   ```R
   source("RScripts/DEP_Analysis.R")
   ```

4. Check the `Results_Files` and `Plots` directories for outputs.

---

## Outputs
- **PCA Plots**: Visualizations of variance among samples.
- **Volcano Plots**: Highlighting significant proteins with fold changes and p-values.
- **Heatmaps**: Clustered visualizations of significant protein expression.
- **Differential Expression Results**: CSV files containing p-values, adjusted p-values, and log fold changes.
- **Venn Diagrams**: Overlap of significant proteins across comparisons.

---

## Citation
If you use this pipeline in your research, please cite this repository and relevant R libraries.

---

## Contact
For questions or issues, please contact:
- **Name**: [Arun Boddapati]
- **Email**: [arunbodd (at) outlook.com]

---

