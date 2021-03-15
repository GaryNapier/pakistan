#! /bin/bash

# master_script_pakistan.sh

set -e
set -u
set -o pipefail

# ------------------------------------------------------------------------------

# Master script for Pakistan data analysis

# RUN
# Assume run from pakistan/ :
# master_script_pakistan.sh

# ------------------------------------------------------------------------------

# Setup - GLOBAL

cd ~/pakistan

chmod +x ~/pipeline/install.sh

~/pipeline/install.sh && echo "installed pipeline scripts to ~/bin/"
printf "\n"

# Variables
study_accession=PAKISTAN_ALL
gvcf_file_suffix=.g.vcf.gz

# Directories
main_metadata_dir=~/metadata/
local_metadata_dir=metadata/
vcf_dir=~/vcf/
ref_dir=~/refgenome/
dist_and_pca_dir=dist_and_pca/
fasta_dir=fasta/
newick_output_dir=newick/
tmp_dir=tmp/

# Files
main_metadata_file=${main_metadata_dir}tb_data_18_02_2021.csv
pakistan_unpublished_metadata=${main_metadata_dir}pakistan_data_non_mixed.csv
pakistan_metadata_file=${local_metadata_dir}pakistan_metadata.csv
sample_list_file=${tmp_dir}${study_accession}.samps.csv
ref_fasta_file=${ref_dir}MTB-h37rv_asm19595v2-eg18.fa
ex_loci_file=${ref_dir}excluded_loci_rep_regions_dr_regions.bed
val_multi_vcf_file=${vcf_dir}${study_accession}.val.gt.g.vcf.gz
filt_multi_vcf_file=${vcf_dir}${study_accession}.filt.val.gt.g.vcf.gz
dist_file=${dist_and_pca_dir}${study_accession}.dist.dist
dist_id_file=${dist_and_pca_dir}${study_accession}.dist.dist.id
unfilt_fasta_file=${fasta_dir}${study_accession}.filt.val.gt.g.snps.fa
fasta_file_base=$(basename -- ${unfilt_fasta_file})
iqtree_file=${newick_output_dir}${fasta_file_base}.iqtree
treefile=${newick_output_dir}${fasta_file_base}*treefile

# Make directories
if [ ! -d ${tmp_dir} ]; then
    mkdir ${tmp_dir}
fi

# ------------------------------------------------------------------------------

# Metadata wrangling

echo "------------------------------------------------------------------------------"

echo "Getting and cleaning Pakistan metadata"
printf "\n"
echo "Running r_scripts/clean_metadata.R - outputs ${pakistan_metadata_file}"
printf "\n"
set -x
# Rscript r_scripts/clean_metadata.R <main_metadata_file>  <other_metadata_file>            <outfile>
Rscript   r_scripts/clean_metadata.R ${main_metadata_file} ${pakistan_unpublished_metadata} ${pakistan_metadata_file}
set +x
echo "------------------------------------------------------------------------------"
printf "\n"

# ------------------------------------------------------------------------------

# Concat VCFs

if [ ! -f ${val_multi_vcf_file} ]; then


    # Get list of samples to pass to variant_calling_and_concat_gvcfs.sh
    cat ${pakistan_metadata_file} | \
    csvtk grep -f study_accession_word -p ${study_accession} | \
    csvtk cut -f wgs_id | tail -n +2 > ${sample_list_file}

    echo "------------------------------------------------------------------------------"

    echo "Variant calling and concat VCFs of samples"
    printf "\n"
    echo "Running shell_scripts/variant_calling_and_concat_gvcfs.sh - outputs file ${val_multi_vcf_file}"
    set -x
    # variant_calling_and_concat_gvcfs.sh <study_accession>   <sample_list_file>  <vcf_dir>  <gvcf_file_suffix>  <ref_file>          <threads>
    variant_calling_and_concat_gvcfs.sh   ${study_accession}  ${sample_list_file} ${vcf_dir} ${gvcf_file_suffix} ${ref_fasta_file}   10
    set +x
    echo "------------------------------------------------------------------------------"
    printf "\n"
else
    echo "------------------------------------------------------------------------------"
    echo "File ${val_multi_vcf_file} already exists, skipping shell_scripts/variant_calling_and_concat_gvcfs.sh"
    echo "------------------------------------------------------------------------------"
    printf "\n"
fi

# ------------------------------------------------------------------------------

# Variant filtering

if [ ! -f ${filt_multi_vcf_file} ]; then

    echo "------------------------------------------------------------------------------"

    echo "Variant filtering"
    printf "\n"
    echo "Running shell_scripts/variant_filtering.sh - outputs file ${filt_multi_vcf_file}"
    printf "\n"
    set -x
    # shell_scripts/variant_filtering.sh    <multi-sample vcf file name>  <output file name>     <bed file name>  <ref fasta>
    # shell_scripts/variant_filtering.sh      ${val_multi_vcf_file}         ${filt_multi_vcf_file} ${ex_loci_file}  ${ref_fasta_file}
    variant_filtering.sh      ${val_multi_vcf_file}         ${filt_multi_vcf_file} ${ex_loci_file}  ${ref_fasta_file}
    set +x
    echo "------------------------------------------------------------------------------"
    printf "\n"

else
    echo "------------------------------------------------------------------------------"
    echo "File ${filt_multi_vcf_file} already exists, skipping shell_scripts/variant_filtering.sh"
    echo "------------------------------------------------------------------------------"
    printf "\n"
fi


# ------------------------------------------------------------------------------

# Plink distances

# Define output file for test of existence - these files are used by r_scripts/transmission_clusters.R

if [ ! -f ${dist_file} ] || [ ! -f ${dist_id_file} ]; then

    echo "------------------------------------------------------------------------------"

    echo "Plink distances"
    printf "\n"

    echo "Running Plink distances command - outputs files to ${dist_and_pca_dir}"
    printf "\n"
    set -x
    # shell_scripts/plink_dist_and_pca.sh   <study_accession>   <vcf_input_file>        <output_dir>
    # shell_scripts/plink_dist_and_pca.sh     ${study_accession}  ${filt_multi_vcf_file}  ${dist_and_pca_dir}
    plink_dist_and_pca.sh     ${study_accession}  ${filt_multi_vcf_file}  ${dist_and_pca_dir}
    set +x
    echo "------------------------------------------------------------------------------"
    printf "\n"
else
    echo "------------------------------------------------------------------------------"
    echo "Files ${dist_file} and ${dist_id_file} already exist, skipping shell_scripts/plink_dist_and_pca.sh"
    echo "------------------------------------------------------------------------------"
    printf "\n"

fi

# ------------------------------------------------------------------------------

# VCF to fasta

# Define output file for test - files used by shell_scripts/iqtree.sh


if [ ! -f ${unfilt_fasta_file} ]; then

    echo "------------------------------------------------------------------------------"

    echo "VCF to fasta"
    printf "\n"

    echo "Running VCF to fasta command - outputs file ${unfilt_fasta_file}"
    set -x
    # Nb vcf2fasta.py (run by vcf2fasta.sh) needs datamash - conda install -c bioconda datamash
    # shell_scripts/vcf2fasta.sh <study_accession>  <vcf_file>              <fasta_output_dir>   <ref_fasta>
    # shell_scripts/vcf2fasta.sh   ${study_accession} ${filt_multi_vcf_file}  ${fasta_dir}         ${ref_fasta_file}
    vcf2fasta.sh   ${study_accession} ${filt_multi_vcf_file}  ${fasta_dir}         ${ref_fasta_file}
    set +x
    echo "------------------------------------------------------------------------------"
    printf "\n"
else
    echo "------------------------------------------------------------------------------"
    echo "File ${unfilt_fasta_file} exists, skipping shell_scripts/vcf2fasta.sh"
    echo "------------------------------------------------------------------------------"
    printf "\n"
fi

# IQ tree

if [ ! -f ${iqtree_file} ] || [ ! -f ${treefile} ]; then

    echo "------------------------------------------------------------------------------"

    echo "IQ tree"
    printf "\n"

    echo "Running shell_scripts/iqtree.sh - outputs files ${iqtree_file}, ${treefile}"
    set -x
    # shell_scripts/iqtree.sh <study_accession>   <fasta_dir>           <newick_output_dir>
    # shell_scripts/iqtree.sh   ${study_accession}  ${unfilt_fasta_file}  ${newick_output_dir}
    iqtree.sh   ${study_accession}  ${unfilt_fasta_file}  ${newick_output_dir}
    set +x
    echo "------------------------------------------------------------------------------"
    printf "\n"
else
    echo "------------------------------------------------------------------------------"
    echo "Files ${iqtree_file} and ${treefile} exist, skipping shell_scripts/iqtree.sh"
    echo "------------------------------------------------------------------------------"
    printf "\n"
fi


# Clean up
rm -r ${tmp_dir}
