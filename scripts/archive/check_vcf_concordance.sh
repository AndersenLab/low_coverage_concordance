#!/bin/bash

module load bcftools

in_vcf="$1"
ref_vcf="$2"
pairs_file="$3"
out="$4"
sp="$5"

# Print the input arguments
echo "Input VCF file: $in_vcf"
echo "Reference VCF file: $ref_vcf"
echo "Pairs file: $pairs_file"
echo "Output file: $out"
echo "Species: $sp"

current_date=$(date +"%Y%m%d_%H%M")
#echo $current_date

# Successfully writes a snp vcf
## save an output file name
query_vcf_path="data/temp/$sp/${current_date}_query_snp.vcf.gz"
query_vcf_index_path="${query_vcf_path}.csi"

# Check if the directory for the query VCF path exists, and create it if it doesn't
query_vcf_dir=$(dirname "$query_vcf_path")
if [ ! -d "$query_vcf_dir" ]; then
    echo "Directory $query_vcf_dir does not exist. Creating it..."
    mkdir -p "$query_vcf_dir"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to create directory $query_vcf_dir"
        exit 1
    fi
fi

echo "Saving Query SNP vcf to: $query_vcf_path"

#bcftools view -i 'TYPE="snp"' --output "$query_vcf_path" -Oz "$in_vcf" 

# Extract SNPs from the input VCF and save to the query VCF path
bcftools view -i 'TYPE="snp"' --output "$query_vcf_path" -Oz "$in_vcf"
if [ $? -ne 0 ]; then
    echo "Error: Failed to extract SNPs from $in_vcf"
    exit 1
fi

# Index the query VCF
#bcftools index -o "$query_vcf_index_path" "$query_vcf_path"
bcftools index -o "$query_vcf_index_path" "$query_vcf_path"
if [ $? -ne 0 ]; then
    echo "Error: Failed to index $query_vcf_path"
    exit 1
fi
#bcftools gtcheck -H --genotype "$ref_vcf" -s 'gt:N2' -s 'qry:N2' data/temp/query_snp.vcf.gz

# Run the gtcheck with a GT-vs-GT comparison
#bcftools gtcheck -H --genotype "$ref_vcf" -s 'gt:N2' -s 'qry:N2' -u 'GT' data/temp/query_snp.vcf.gz > "$out"

bcftools gtcheck -H \
--genotype "$ref_vcf" \
-p AF16,AF16,BRC20232,BRC20232,ECA1657,ECA1657,ECA2115,ECA2115,ECA2299,ECA2299,ECA2667,ECA2667,EG6260,EG6260,EG6267,EG6267,EG6271,EG6271,JU1346,JU1346,JU1348,JU1348,JU2472,JU2472,NIC1052,NIC1052,NIC1645,NIC1645,NIC1667,NIC1667,NIC827,NIC827,NIC828,NIC828,PE887,PE887,QG2908,QG2908,QG4097,QG4097,QG4233,QG4233,QX1410,QX1410,VT847,VT847 \
-e 0 \
-u GT,GT \
"$query_vcf_path" > "$out"

if [ $? -ne 0 ]; then
    echo "Error: Failed to run gtcheck"
    exit 1
fi

echo "gtcheck completed successfully. Output saved to $out"