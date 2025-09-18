# Wound age project

The main goal of this project is to caracterize different timmings of wound healing, and find specific biomarkers (Manuscript in preparation). 

This repostitory hosts the code for all the sections of the project, including data exploration, hierachical clustering, exploration of the DEGs, enrichement analysis, finding of the biomarkers and ploting of the final figures. 
The identification of DEGs was performed using TAC software, and its no included in this repository. 

## Structure of the repository
The repository is structured as follows:
- `data/`: Folder containing the data used in the analysis (Minus the expression files).
- `code/`: Folder containing all the scripts used for the analysis.

The code folder is structured as follows:
    - `01_data_exploration_analysis.R`: Code for the data exploration analysis. Includes Hierarchical clustering, PCA analysis and exploration of the covariates.
    - `02_DEG_analysis.R`: Code for the exploration of the DEGs. Run after TAC
    - `03_enrichment_analysis.R`: Code for the enrichment analysis of the DEGs.
    - `04_orsum_summary.sh`: Code for the summary of the ORSUM analysis.
    - `05_plot_enrichement_results.R.R`: Code for the plotting of the enrichment analysis.
    - `06_biomarker_analysis.R`: Code for the identification of the biomarkers.


Paper figures were generated using these scripts and edited in Inkescape when required. 