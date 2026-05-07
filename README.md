
**Multi-omics Profiling of Cell-Cycle Spatial Niches in High-Grade Serous Carcinoma**

This repository contains the code utilized in the Master's thesis: "Multi-omics profiling of cell-cycle spatial niches in high-grade serous carcinoma."

📖 Overview

This project investigates the mechanisms operating within the spatially cell-cycle defined highly proliferative niches in High-Grade Serous Carcinoma (HGSC). To achieve this, the study integrates multi-omics data by leveraging spatial transcriptomics (GeoMx) alongside tissue-based cyclic immunofluorescence (t-CyCIF).

🗂️ Repository Structure

The code is organized into three main analytical stages:

1. Tribus Phenotyping and Gating
Directory: /tribus_phenotyping_main_scripts/

Purpose: Phenotypes the cells according to the details predefined in a logic table.

Input: Quantified segmented data, logic table.

Output: Assigned cell phenotypes.

Logic Table Utilized:

<img width="540" height="300" alt="Screenshot 2026-05-06 at 16 43 23" src="https://github.com/user-attachments/assets/309f1370-38b5-4d42-bd6f-3617422b9b19" />


2. Image Visualization and ROI Transfer
Directory: /Rerun_visualization_and_ROI_transfer/

Purpose: Maps and extracts cells corresponding to specific Regions of Interest (ROIs). The ROIs from adjacent slides profiled with GeoMx were manually transferred to t-CyCIF images and saved as CSV files. Using the ROI coordinates in these CSVs, the cells lying within them were extracted and annotated to their corresponding ROI. Cells outside these defined ROIs were assigned an NA value and saved.

Input: Coordinates of the ROIs (CSV files), Quantified segmented data.

Output: Annotated single-cell data with ROI assignments.

3. Downstream Analyses
This stage includes the primary statistical, computational, and visual analyses of the multi-omics data, split across three directories:

/GSEA/

Contains scripts to prepare pseudo-bulk expression data for Gene Set Enrichment Analysis. This includes formatting the expression matrix, defining phenotypic categories, and structuring other prerequisite datasets required for compatibility with Broad Institute GSEA software.

/Analysis_in_R/

Contains R scripts utilized for Differential Gene Expression (DGE) analysis and Differential Pathway Enrichment, primarily leveraging the limma package.

/downstream_analysis_thesis/

Houses the general downstream plotting and statistical analysis pipelines, including:

- Visualization of DGE and Differential Pathway Enrichment results (generated in R).

- Generation of Voronoi plots based on Tribus guessed and gated phenotypes within the defined ROIs.

- Statistical comparisons of the counts and proportions of proliferative cells within specific proliferative niches.

- Plotting the expression profiles of proliferative signature genes within these niches.

- Implementation of Machine Learning models (Random Forest) to cluster and classify transcriptional states.



✍️ Author

Nika Mikhailava

Master's Thesis Project
