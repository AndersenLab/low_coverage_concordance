#!/bin/bash
# Load necessary modules
module load bcftools

# Define variables
ref_vcf="$1"
seq_vcf="$2"  # Update with the actual path to your sequence VCF file
sample="$3"  # Sample name passed as an argument
current_date=$(date +"%Y%m%d")
output_dir="data/temp/${current_date}_${sample}"

# Create the output directory if it doesn't exist
mkdir -p ${output_dir}

##### Process the reference VCF #####
# Filter the reference VCF to include only the specified sample
output_file="${output_dir}/SNP.vcf.gz"
if [ ! -f "${output_file}" ]; then
    bcftools view -s ${sample} -i 'TYPE="snp"' -Oz ${ref_vcf} > ${output_file}
    if [ $? -ne 0 ]; then
        echo "Error: Failed to filter reference VCF for sample ${sample}"
        exit 1
    fi
    echo "Filtered reference VCF for sample ${sample} and saved to ${output_file}"
else
    echo "File ${output_file} already exists. Skipping this step."
fi

# Generate stats for the filtered reference VCF
stats_file="${output_dir}/SNP.stats.txt"
if [ ! -f "${stats_file}" ]; then
    bcftools stats ${output_file} > ${stats_file}
    if [ $? -ne 0 ]; then
        echo "Error: Failed to generate stats for ${output_file}"
        exit 1
    fi
    echo "Generated stats for ${output_file} and saved to ${stats_file}"
else
    echo "File ${stats_file} already exists. Skipping this step."
fi

# Filter the SNP VCF to include only sites where the sample is `alt`
output_file="${output_dir}/SNP_alt.vcf.gz"
if [ ! -f "${output_file}" ]; then
    bcftools view -i 'GT[0]="alt"' -Oz ${output_dir}/SNP.vcf.gz > ${output_file}
    if [ $? -ne 0 ]; then
        echo "Error: Failed to filter SNP VCF for alt sites"
        exit 1
    fi
    echo "Filtered SNP VCF for alt sites and saved to ${output_file}"
else
    echo "File ${output_file} already exists. Skipping this step."
fi

# Generate stats for the filtered SNP VCF
stats_file="${output_dir}/SNP_alt.stats.txt"
if [ ! -f "${stats_file}" ]; then
    bcftools stats ${output_file} > ${stats_file}
    if [ $? -ne 0 ]; then
        echo "Error: Failed to generate stats for ${output_file}"
        exit 1
    fi
    echo "Generated stats for ${output_file} and saved to ${stats_file}"
else
    echo "File ${stats_file} already exists. Skipping this step."
fi

# Query the SNP sites where the sample has the alt genotype and report the genotype
output_file="${output_dir}/SNP_alt_gts.txt"
if [ ! -f "${output_file}" ]; then
    bcftools query -f'%CHROM\t%POS[\t%SAMPLE=%GT]\n' ${output_dir}/SNP_alt.vcf.gz > ${output_file}
    if [ $? -ne 0 ]; then
        echo "Error: Failed to query SNP sites for alt genotype"
        exit 1
    fi
    echo "Queried SNP sites for alt genotype and saved to ${output_file}"
else
    echo "File ${output_file} already exists. Skipping this step."
fi

# Pull only the first two columns (CHROM and POS) from the SNP sites where the sample has the alt genotype
output_file="${output_dir}/SNP_alt.txt"
if [ ! -f "${output_file}" ]; then
    cut -f1,2 ${output_dir}/SNP_alt_gts.txt > ${output_file}
    if [ $? -ne 0 ]; then
        echo "Error: Failed to extract CHROM and POS columns"
        exit 1
    fi
    echo "Extracted CHROM and POS columns and saved to ${output_file}"
else
    echo "File ${output_file} already exists. Skipping this step."
fi

#### Process the sequence VCF #####
# Filter the sequence VCF to include only the specified sample
output_file="${output_dir}/seq.vcf.gz"
if [ ! -f "${output_file}" ]; then
    bcftools view -s ${sample} -Oz ${seq_vcf} > ${output_file}
    if [ $? -ne 0 ]; then
        echo "Error: Failed to filter sequence VCF for sample ${sample}"
        exit 1
    fi
    echo "Filtered sequence VCF for sample ${sample} and saved to ${output_file}"
else
    echo "File ${output_file} already exists. Skipping this step."
fi

# Generate stats for the filtered sequence VCF
stats_file="${output_dir}/seq.stats.txt"
if [ ! -f "${stats_file}" ]; then
    bcftools stats ${output_file} > ${stats_file}
    if [ $? -ne 0 ]; then
        echo "Error: Failed to generate stats for ${output_file}"
        exit 1
    fi
    echo "Generated stats for ${output_file} and saved to ${stats_file}"
else
    echo "File ${stats_file} already exists. Skipping this step."
fi

# Filter the sequence VCF to include only SNPs
output_file="${output_dir}/seq_SNP.vcf.gz"
if [ ! -f "${output_file}" ]; then
    bcftools view -i 'TYPE="snp"' -Oz ${output_dir}/seq.vcf.gz > ${output_file}
    if [ $? -ne 0 ]; then
        echo "Error: Failed to filter sequence VCF for SNPs"
        exit 1
    fi
    echo "Filtered sequence VCF for SNPs and saved to ${output_file}"
else
    echo "File ${output_file} already exists. Skipping this step."
fi

# Generate stats for the filtered sequence SNP VCF
stats_file="${output_dir}/seq_SNP.stats.txt"
if [ ! -f "${stats_file}" ]; then
    bcftools stats ${output_file} > ${stats_file}
    if [ $? -ne 0 ]; then
        echo "Error: Failed to generate stats for ${output_file}"
        exit 1
    fi
    echo "Generated stats for ${output_file} and saved to ${stats_file}"
else
    echo "File ${stats_file} already exists. Skipping this step."
fi

# Filter the sequence SNP VCF to include only sites where the sample is `alt`
output_file="${output_dir}/seq_SNP_alt.vcf.gz"
if [ ! -f "${output_file}" ]; then
    bcftools view -i 'GT[0]="alt"' -Oz ${output_dir}/seq_SNP.vcf.gz > ${output_file}
    if [ $? -ne 0 ]; then
        echo "Error: Failed to filter sequence SNP VCF for alt sites"
        exit 1
    fi
    echo "Filtered sequence SNP VCF for alt sites and saved to ${output_file}"
else
    echo "File ${output_file} already exists. Skipping this step."
fi

# Generate stats for the filtered sequence SNP VCF
stats_file="${output_dir}/seq_SNP_alt.stats.txt"
if [ ! -f "${stats_file}" ]; then
    bcftools stats ${output_file} > ${stats_file}
    if [ $? -ne 0 ]; then
        echo "Error: Failed to generate stats for ${output_file}"
        exit 1
    fi
    echo "Generated stats for ${output_file} and saved to ${stats_file}"
else
    echo "File ${stats_file} already exists. Skipping this step."
fi

# Use the SNP sites from the reference VCF to filter the sequence VCF
output_file="${output_dir}/seq_WI_fil.vcf.gz"
if [ ! -f "${output_file}" ]; then
    bcftools view -T ${output_dir}/SNP_alt.txt -Oz "${output_dir}/seq.vcf.gz" > ${output_file}
    if [ $? -ne 0 ]; then
        echo "Error: Failed to filter sequence VCF using reference SNP sites"
        exit 1
    fi
    echo "Filtered sequence VCF using reference SNP sites and saved to ${output_file}"
else
    echo "File ${output_file} already exists. Skipping this step."
fi

# Generate stats for the filtered sequence VCF
stats_file="${output_dir}/seq_SNP_alt_WI_fil.stats.txt"
if [ ! -f "${stats_file}" ]; then
    bcftools stats ${output_file} > ${stats_file}
    if [ $? -ne 0 ]; then
        echo "Error: Failed to generate stats for ${output_file}"
        exit 1
    fi
    echo "Generated stats for ${output_file} and saved to ${stats_file}"
else
    echo "File ${stats_file} already exists. Skipping this step."
fi

# Query the SNP sites where the sample has the alt genotype and report the genotype
output_file="${output_dir}/seq_WI_fil_gts.txt"
if [ ! -f "${output_file}" ]; then
    bcftools query -f'%CHROM\t%POS[\t%SAMPLE=%GT]\n' "${output_dir}/seq_WI_fil.vcf.gz" > ${output_file}
    if [ $? -ne 0 ]; then
        echo "Error: Failed to query SNP sites for alt genotype"
        exit 1
    fi
    echo "Queried SNP sites for alt genotype and saved to ${output_file}"
else
    echo "File ${output_file} already exists. Skipping this step."
fi

echo "Concordance analysis completed for sample ${sample}. Check the stats files for details."