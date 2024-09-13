#!/bin/bash

# Elegans

Rscript check_concordance_V2.R \
'data/raw/20240821_c_elegans_V3/WI_strain_snps_alt_gts.tsv' \
'data/raw/20240821_c_elegans_V3/seq_WI_fil_gts.tsv' \
'data/raw/20240821_c_elegans_V3/samples.txt' \
'data/proc/20240821_c_elegans_V3'

# Briggsae

Rscript check_concordance_V2.R \
'data/raw/20240821_c_briggsae_V3/WI_strain_snps_alt_gts.tsv' \
'data/raw/20240821_c_briggsae_V3/seq_WI_fil_gts.tsv' \
'data/raw/20240821_c_briggsae_V3/samples.txt' \
'data/proc/20240821_c_briggsae_V3'
