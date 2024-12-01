#!/bin/bash

# for NCBI its essential to specify the user-agent
# some weird bot-filter gives a 403 error otherwise

cd data/01_raw
wget --user-agent="Mozilla" pmc.ncbi.nlm.nih.gov/articles/instance/6280667/bin/NIHMS1504008-supplement-Supplementary_Tables_S1-S22.xlsx
wget --user-agent="Mozilla" https://github.com/biodev/beataml2.0_data/raw/main/beataml_wv1to4_clinical.xlsx
wget --user-agent="Mozilla" https://github.com/biodev/beataml2.0_data/raw/main/beataml_waves1to4_sample_mapping.xlsx
wget --user-agent="Mozilla" https://github.com/biodev/beataml2.0_data/raw/main/beataml_wes_wv1to4_mutations_dbgap.txt
wget --user-agent="Mozilla" https://github.com/biodev/beataml2.0_data/raw/main/beataml_waves1to4_norm_exp_dbgap.txt
wget --user-agent="Mozilla" https://github.com/biodev/beataml2.0_data/raw/main/beataml_probit_curve_fits_v4_dbgap.txt
wget --user-agent="Mozilla" https://github.com/biodev/beataml2.0_data/raw/main/beataml_drug_families.xlsx
wget --user-agent="Mozilla" https://github.com/biodev/beataml2.0_data/raw/main/beataml_wv1to4_raw_inhibitor_v4_dbgap.txt
wget --user-agent="Mozilla" https://github.com/biodev/beataml2.0_data/raw/main/beataml2_manuscript_vg_cts.xlsx