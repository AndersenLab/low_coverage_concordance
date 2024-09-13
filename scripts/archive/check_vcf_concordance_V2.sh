#!/bin/bash

module load bcftools

# Define variables
ref_vcf="$1"
seq_vcf="$2"  # Update with the actual path to your sequence VCF file
samples="$3"  # Path to file containing list of samples to include in analysis
out_dir_name="$4"

current_date=$(date +"%Y%m%d")

## Replace with parameter
output_dir="${current_date}_${out_dir_name}"

# Create the output directory if it doesn't exist
mkdir -p ${output_dir}


# Filter the WI vcf to just SNPs and samples of interest and then pull sites where at least one sample has the ALT allele
output_file="${output_dir}/WI_strain_snps.vcf.gz"
if [ ! -f "${output_file}" ] || [ ! -s "${output_file}" ]; then
    bcftools view -S ${samples} -i 'TYPE="snp"' -Ou ${ref_vcf} | bcftools view -i 'N_PASS(GT="alt")>=1' -Oz > ${output_file}
    if [ $? -ne 0 ]; then
        echo "Error: Failed to filter reference VCF" 
        exit 1
    fi
    echo "Filtered reference VCF for samples and saved to ${output_file}"
else
    echo "File ${output_file} already exists. Skipping this step."
fi

# Pull a list of sites and sample GTs from the filtered VCF
# Filter the WI vcf to just SNPs and samples of interest and then pull sites where at least one sample has the ALT allele
output_file="${output_dir}/WI_strain_snps_alt_gts.tsv"
if [ ! -f "${output_file}" ] || [ ! -s "${output_file}" ]; then
    bcftools query -f'%CHROM\t%POS[\t%SAMPLE=%GT]\n' "${output_dir}/WI_strain_snps.vcf.gz" > ${output_file}
    if [ $? -ne 0 ]; then
        echo "Error: Failed to filter reference VCF" 
        exit 1
    fi
    echo "Filtered reference VCF for samples and saved to ${output_file}"
else
    echo "File ${output_file} already exists. Skipping this step."
fi

## Process to just get the positons to filter the seq VCF
# Pull only the first two columns (CHROM and POS) from the SNP sites where the sample has the alt genotype
output_file="${output_dir}/WI_strain_snps_alt_pos.tsv"
if [ ! -f "${output_file}" ] || [ ! -s "${output_file}" ]; then
    cut -f1,2 "${output_dir}/WI_strain_snps_alt_gts.tsv" > ${output_file}
    if [ $? -ne 0 ]; then
        echo "Error: Failed to extract CHROM and POS columns"
        exit 1
    fi
    echo "Extracted CHROM and POS columns and saved to ${output_file}"
else
    echo "File ${output_file} already exists. Skipping this step."
fi

# Filter the Seq vcf to just the same samples and sites from the WI vcf
output_file="${output_dir}/seq_WI_fil.vcf.gz"
if [ ! -f "${output_file}" ] || [ ! -s "${output_file}" ]; then
    bcftools view -S ${samples} -T "${output_dir}/WI_strain_snps_alt_pos.tsv" -Oz ${seq_vcf} > ${output_file}
    if [ $? -ne 0 ]; then
        echo "Error: Failed to filter sequence VCF using reference SNP sites"
        exit 1
    fi
    echo "Filtered sequence VCF using reference SNP sites and saved to ${output_file}"
else
    echo "File ${output_file} already exists. Skipping this step."
fi

# Query the SNP sites where the sample has the alt genotype and report the genotype
output_file="${output_dir}/seq_WI_fil_gts.tsv"
if [ ! -f "${output_file}" ] || [ ! -s "${output_file}" ]; then
    bcftools query -f'%CHROM\t%POS[\t%SAMPLE=%GT]\n' --print-header "${output_dir}/seq_WI_fil.vcf.gz" > ${output_file}
    if [ $? -ne 0 ]; then
        echo "Error: Failed to query SNP sites for alt genotype"
        exit 1
    fi
    echo "Queried SNP sites for alt genotype and saved to ${output_file}"
else
    echo "File ${output_file} already exists. Skipping this step."
fi

echo "Concordance analysis completed for sample ${sample}. Check the stats files for details."