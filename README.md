# Simultaneous integration of gene expression and nutrient availability for studying metabolism of hepatocellular carcinoma

This repository contains all files and data to reproduce the publication "Simultaneous Integration of Gene Expression and Nutrient Availability for Studying the Metabolism of Hepatocellular Carcinoma Cell Lines" by Ewelina Węglarz-Tomczak, Thierry D.G.A. Mondeel, Diewertje G.E. Piebes and Hans V. Westerhoff. 

If you use or reference the results from this work please cite the publication https://doi.org/10.3390/biom11040490.

# Steps to reproduce the analysis in the publication

The data and code in this repository should allow full reproduction of the data and figures discussed in the publication. 

All results can be reproduced by executing the various Jupyter Python notebooks and MATLAB files (In the "Code" folder) in this order:

* Microarray_analysis.ipynb
    - Output: microarray_data_with_entrez_genes.csv
* RNAseq_analysis.ipynb
    - Uses microarray_data_with_entrez_genes.csv to fill zeros
    - output: RNAseq_data_with_entrez_genes.csv
* patch_Recon3D_301.m
    - Generates Recon3DModel_301_patched.mat
* Generate_medium_models.m
    - Generates the medium specific models (without RAS bounds, here we just set exchange reaction bounds)
    - Generates Recon3DModel_301_patched_medium_1.xlsx
* calculate_RAS_and_constrain_medium_models.m
    - Starts from RNAseq_data_with_entrez_genes.csv
    - Uses getScore.m
    - Generates RAS_scores.mat
    - Produces the biomass vs. alpha plot
    - Saves the constrained models for NNR = 0.003607
    - Performs the FVA analyses and makes the resulting plots
* Export_Recon3D_301_patched_to_json.ipynb
    - Generates Recon3DModel_301_patched.json
* Plot_TPM_and_RAS_comparisons.ipynb
    - Starts from RNAseq_data_with_entrez_genes.csv and RAS_scores.mat
    - Uses Recon3DModel_301_patched.json
    - Generates RAS_outliers.xlsx and TPM_RAS_comparison figure
* Generate_supplementary_tables.ipynb
    - Creates 3 suppl. Excel tables: Biomass_composition.xlsx, Metabolites_allowed_to_be_secreted.xlsx, Oxphos_reactions.xlsx
