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

cd ~/pipeline
chmod +x install.sh
./install.sh && echo "installed pipeline scripts to ~/bin/"
printf "\n"

cd ~/pakistan

# Variables
date_change=false
study_accession=PAKISTAN_ALL
gvcf_file_suffix=.g.vcf.gz
date_column=year
id_column=wgs_id
cut_year=50

# Directories

# Existing directories
main_metadata_dir=~/metadata/
vcf_dir=~/vcf/
ref_dir=~/refgenome/

# Out directories
local_metadata_dir=metadata/
dist_and_pca_dir=dist_and_pca/
fasta_dir=fasta/
newick_output_dir=newick/
itol_dir=itol/
tmp_dir=tmp/
beast_results_dir=beast_results/
beast_xml_dir=beast_xml/
transphylo_results_dir=transphylo_results/

# Make directories if not exist
if [ ! -d ${local_metadata_dir} ]; then
    mkdir ${local_metadata_dir}
fi

if [ ! -d ${dist_and_pca_dir} ]; then
    mkdir ${dist_and_pca_dir}
fi

if [ ! -d ${fasta_dir} ]; then
    mkdir ${fasta_dir}
fi

if [ ! -d ${newick_output_dir} ]; then
    mkdir ${newick_output_dir}
fi

if [ ! -d ${itol_dir} ]; then
    mkdir ${itol_dir}
fi

if [ ! -d ${tmp_dir} ]; then
    mkdir ${tmp_dir}
fi

if [ ! -d ${beast_results_dir} ]; then
    mkdir ${beast_results_dir}
fi

if [ ! -d ${beast_xml_dir} ]; then
    mkdir ${beast_xml_dir}
fi

if [ ! -d ${transphylo_results_dir} ]; then
    mkdir ${transphylo_results_dir}
fi


# Files

# Existing files
main_metadata_file=${main_metadata_dir}tb_data_18_02_2021.csv
pakistan_unpublished_metadata=${main_metadata_dir}pakistan_data_non_mixed.csv
ref_fasta_file=${ref_dir}MTB-h37rv_asm19595v2-eg18.fa
ex_loci_file=${ref_dir}excluded_loci_rep_regions_dr_regions.bed


# Out files

# metadata/samples
pakistan_metadata_file=${local_metadata_dir}pakistan_metadata.csv
pakistan_metadata_file_dated=${local_metadata_dir}pakistan_metadata_dated.csv
tbprofiler_file=${local_metadata_dir}tbprofiler_results.pakistan.txt
sample_list_file=${tmp_dir}${study_accession}.samps.csv
dated_samples_file=${tmp_dir}dated_samples.txt
beast_clusters_file=${local_metadata_dir}${study_accession}.clusters.csv
# vcf
val_multi_vcf_file=${vcf_dir}${study_accession}.val.gt.g.vcf.gz
filt_multi_vcf_file=${vcf_dir}${study_accession}.filt.val.gt.g.vcf.gz
dated_samps_vcf=${vcf_dir}${study_accession}.dated_samps.vcf.gz
dated_samps_vcf_unzip=${vcf_dir}${study_accession}.dated_samps.vcf
# dist
dist_file=${dist_and_pca_dir}${study_accession}.dist.dist
dist_id_file=${dist_and_pca_dir}${study_accession}.dist.dist.id
# fasta
unfilt_fasta_file=${fasta_dir}${study_accession}.filt.val.gt.g.snps.fa
fasta_file_base=$(basename -- ${unfilt_fasta_file})
fasta_dated_samples_file=${fasta_dir}${study_accession}.dated_samps.snps.fa
fasta_file_annotated_with_dates=${fasta_dir}${study_accession}.dated.fa
# newick
iqtree_file=${newick_output_dir}${fasta_file_base}.iqtree
treefile=${newick_output_dir}${fasta_file_base}*treefile
# xml files for BEAST
xml_file=${beast_xml_dir}${study_accession}.xml
xml_file_plus_const=${beast_xml_dir}${study_accession}_plus_const.xml
# Beast output files - n.b. saves the .log and .trees files to the current directory and takes the study_accession
# So files are ./<study_accession>.log and ./<study_accession>.trees
original_beast_log_file=*.log
original_beast_trees_file=*.trees
original_beast_state_file=*.xml.state

new_beast_log_file=${beast_results_dir}${study_accession}*.log
new_beast_trees_file=${beast_results_dir}${study_accession}*.trees
new_beast_state_file=${beast_results_dir}${study_accession}*.xml.state
# mcc
mcc_tree=${beast_results_dir}${study_accession}.mcc.tree
# TransPhylo
transphylo_es_table_file=${transphylo_results_dir}${study_accession}.es_table.csv

# ------------------------------------------------------------------------------

# Remove files that change if the dates change

#  set date_change=true in variables section

if [ ${date_change} = true ] ; then
    echo "-------------------------------------------------"
    echo "Removing files that change when the dates change"
    echo "-------------------------------------------------"
    set -x
    rm -f ${pakistan_metadata_file} \
     ${dated_samps_vcf} \
     ${dated_samps_vcf_unzip} \
     ${fasta_dated_samples_file} \
     ${fasta_file_annotated_with_dates} \
     ${xml_file} \
     ${xml_file_plus_const} \
     ${new_beast_log_file} \
     ${new_beast_trees_file} \
     ${new_beast_state_file} \
     ${mcc_tree} \
     ${transphylo_results_dir}*
    set +x
    printf "\n"
fi

# ------------------------------------------------------------------------------

# Metadata wrangling

echo "------------------------------------------------------------------------------"

if [ ! -f ${pakistan_metadata_file} ]; then

    echo "Getting and cleaning Pakistan metadata"
    printf "\n"
    echo "Running r_scripts/clean_metadata.R - outputs ${pakistan_metadata_file}"
    printf "\n"
    set -x
    # Rscript r_scripts/clean_metadata.R    <main_metadata_file>   <other_metadata_file>            <tbprofiler_results.pakistan.txt> <outfile>
    Rscript r_scripts/clean_metadata.R      ${main_metadata_file}  ${pakistan_unpublished_metadata} ${tbprofiler_file}                ${pakistan_metadata_file}
    set +x
    echo "------------------------------------------------------------------------------"
    printf "\n"
else
    echo "------------------------------------------------------------------------------"
    echo "File ${pakistan_metadata_file} already exists, skipping clean_metadata.R"
    echo "------------------------------------------------------------------------------"
    printf "\n"
fi

# ------------------------------------------------------------------------------

# NEED TO SUBSET SAMPLES TO JUST THOSE WITH DATES

# Subset metadata to dated samples only
cat ${pakistan_metadata_file} | csvtk grep -f ${date_column} -p "NA" -v > ${pakistan_metadata_file_dated}

# Get those samples for subsetting
cat ${pakistan_metadata_file_dated} | csvtk cut -f wgs_id | tail -n +2 > ${dated_samples_file}

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
    # shell_scripts/variant_filtering.sh    ${val_multi_vcf_file}         ${filt_multi_vcf_file} ${ex_loci_file}  ${ref_fasta_file}
    variant_filtering.sh                    ${val_multi_vcf_file}         ${filt_multi_vcf_file} ${ex_loci_file}  ${ref_fasta_file}
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

# Create VCF of just dated samples for Beast and everything else dated

if [ ! -f ${dated_samps_vcf} ]; then
    echo "------------------------------------------------------------------------------"
    echo "Creating VCF for dated samples only - outputs ${dated_samps_vcf}"
    set -x
    bcftools view -S ${dated_samples_file} ${filt_multi_vcf_file} | bcftools filter -e 'AC==0 || AC==AN' - -o ${dated_samps_vcf} -Oz
    set +x
    echo "------------------------------------------------------------------------------"
    printf "\n"
else
    echo "------------------------------------------------------------------------------"
    echo "File ${dated_samps_vcf} already exists, skipping creation of VCF for dated samples only"
    echo "------------------------------------------------------------------------------"

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
    # shell_scripts/plink_dist_and_pca.sh   ${study_accession}  ${filt_multi_vcf_file}  ${dist_and_pca_dir}
    plink_dist_and_pca.sh                   ${study_accession}  ${filt_multi_vcf_file}  ${dist_and_pca_dir}
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
    # vcf2fasta.sh   <study_accession>  <vcf_file>              <fasta_output_dir>   <ref_fasta>
    vcf2fasta.sh     ${study_accession} ${filt_multi_vcf_file}  ${fasta_dir}         ${ref_fasta_file}
    set +x
    echo "------------------------------------------------------------------------------"
    printf "\n"
else
    echo "------------------------------------------------------------------------------"
    echo "File ${unfilt_fasta_file} exists, skipping shell_scripts/vcf2fasta.sh"
    echo "------------------------------------------------------------------------------"
    printf "\n"
fi

# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------

# Tree annotation - itol

if [ ! -f ${itol_dir}itol.dr.txt ] || [ ! -f ${itol_dir}itol.lineages.txt ] || [ ! -f ${itol_dir}itol.major_lineages.txt ]; then
    echo "------------------------------------------------------------------------------"
    echo "Running itol_templates.R"
    printf "\n"
    # itol_templates.R  <itol_templates_location>
    itol_templates.R    ${itol_dir}
    echo "itol_templates.R done"
    printf "\n"
else
    echo "------------------------------------------------------------------------------"
    echo "${itol_dir} template files exist, skipping itol_templates.R"
    echo "------------------------------------------------------------------------------"
    printf "\n"
fi

itol_out_files=${itol_dir}$(basename ${pakistan_metadata_file%.csv})*

if compgen -G ! ${itol_out_files} > /dev/null; then
    echo "------------------------------------------------------------------------------"
    echo "Running itol_annotation.R"
    printf "\n"
    # itol_annotation.R <metadata_file>             <itol_location>
    itol_annotation.R   ${pakistan_metadata_file}   ${itol_dir}
    echo "itol_annotation.R done"
    printf "\n"
else
    echo "------------------------------------------------------------------------------"
    echo ${itol_out_files} "files exist, skipping itol_annotation.R"
    echo "------------------------------------------------------------------------------"
    printf "\n"
fi


# ------------------------------------------------------------------------------

# Create fasta file from VCF file - dated samples

if [ ! -f ${fasta_dated_samples_file} ]; then

    echo "------------------------------------------------------------------------------"
    echo "Running vcf2fasta.py - outputs ${fasta_dated_samples_file}"
    printf "\n"
    set -x
    vcf2fasta.py --vcf ${dated_samps_vcf} --ref ${ref_fasta_file} --threads 20 --snps --snps-no-filt
    mv ${vcf_dir}$(basename -- ${fasta_dated_samples_file} ) ${fasta_dir}
    set +x
    echo "------------------------------------------------------------------------------"
    printf "\n"
else
    echo "------------------------------------------------------------------------------"
    echo "${fasta_dated_samples_file} file exists, skipping vcf2fasta.sh"
    echo "------------------------------------------------------------------------------"
    printf "\n"
fi

# Filter fasta file to just dated samples
# grep -f ${dated_samples_file} -A1 ${unfilt_fasta_file} > ${fasta_dated_samples_file}
# seqtk subseq ${unfilt_fasta_file} ${dated_samples_file} > ${fasta_dated_samples_file}

# Filter out redundant SNPs from dated samples fasta (filtering on those samples only produces many monomorphic sites)
# snp-sites -m -o ${fasta_dated_samples_file}.tmp ${fasta_dated_samples_file} && mv ${fasta_dated_samples_file}.tmp ${fasta_dated_samples_file}

# ------------------------------------------------------------------------------

#  Annotate fasta file of dated samples with dates

if [ ! -f ${fasta_file_annotated_with_dates} ]; then

    echo "------------------------------------------------------------------------------"
    echo "Dated fasta file"
    printf "\n"

    echo "Running fasta_add_annotations.py - outputs ${fasta_file_annotated_with_dates}"
    printf "\n"

    # Run fasta_add_annotations.py to concatenate the metadata collection date with the sample ID in the fasta file
    # e.g. >ERR1234 -> >ERR1234_01_01_10
    # - see https://github.com/pathogenseq/pathogenseq-scripts/tree/master/scripts
    set -x
    fasta_add_annotations.py --fasta ${fasta_dated_samples_file} --csv ${pakistan_metadata_file_dated} --out ${fasta_file_annotated_with_dates} --annotations ${date_column} --id-key ${id_column}
    set +x
    echo "------------------------------------------------------------------------------"
    printf "\n"
else
    echo "------------------------------------------------------------------------------"
    echo "File ${fasta_file_annotated_with_dates} exists, skipping fasta_add_annotations.py"
    echo "------------------------------------------------------------------------------"
    printf "\n"
fi

# ------------------------------------------------------------------------------

# BEAST

# ------------------------------------------------------------------------------

# Make xml file

if [ ! -f ${xml_file} ]; then
    echo "------------------------------------------------------------------------------"
    echo "Creating XML file r_scripts/pakistan_babette.R - outputs ${xml_file}"
    set -x
    Rscript r_scripts/pakistan_babette.R -s ${study_accession} -f ${fasta_file_annotated_with_dates} -i ${local_metadata_dir} -o ${xml_file}
    set +x
    echo "------------------------------------------------------------------------------"
    printf "\n"
else
    echo "------------------------------------------------------------------------------"
    echo "File ${xml_file} exists, skipping r_scripts/pakistan_babette.R"
    echo "------------------------------------------------------------------------------"
    printf "\n"

fi

# Add "constant site weights" - https://github.com/andersgs/beast2_constsites

if [ ! -f ${xml_file_plus_const} ]; then
    echo "------------------------------------------------------------------------------"
    echo "Adding constant site weights to XML file with run_b2cs - outputs ${xml_file_plus_const}"
    # Uncompress because run_b2cs needs ordinary vcf
    gunzip ${dated_samps_vcf}
    run_b2cs -x ${xml_file} beast2 ${ref_fasta_file} ${dated_samps_vcf_unzip}
else
    echo "------------------------------------------------------------------------------"
    echo "File ${xml_file_plus_const} exists, skipping run_b2cs"
    echo "------------------------------------------------------------------------------"
    printf "\n"
fi

# Stop here if xml file does not exist

if [ ! -f ${xml_file_plus_const} ]; then
    echo "Stopping script before Beast section - xml file ${xml_file_plus_const} does not exist"
    exit 1
fi

# ------------------------------------------------------------------------------


# Run Beast

if [ ! -f ${new_beast_log_file} ] || [ ! -f ${new_beast_trees_file} ] || [ ! -f ${new_beast_state_file} ];then

    echo "------------------------------------------------------------------------------"
    echo "Run Beast"
    printf "\n"

    echo "Running BEAST - outputs ${original_beast_log_file}, ${original_beast_trees_file} and ${original_beast_state_file}"
    set -x
    # Run Beast and clean up - put in right folder
    beast -threads 20 ${xml_file_plus_const}
    mv ${original_beast_log_file} ${beast_results_dir}
    mv ${original_beast_trees_file} ${beast_results_dir}
    mv ${original_beast_state_file} ${beast_results_dir}
    set +x
    echo "------------------------------------------------------------------------------"
    printf "\n"
else
    echo "------------------------------------------------------------------------------"
    echo "Files ${new_beast_log_file}, ${new_beast_trees_file} and ${new_beast_state_file} exist, skipping beast"
    echo "------------------------------------------------------------------------------"
    printf "\n"

fi

# ------------------------------------------------------------------------------

# Run TreeAnnotator to get MCC tree (input to TransPhylo)

if [ ! -f ${mcc_tree} ]; then

    echo "------------------------------------------------------------------------------"
    echo "Running TreeAnnotator on ${new_beast_trees_file} - outputs ${mcc_tree}"
    set -x
    treeannotator -heights ca -burnin 10 ${new_beast_trees_file} ${mcc_tree}
    set +x
    echo "------------------------------------------------------------------------------"
    printf "\n"
else
    echo "------------------------------------------------------------------------------"
    echo "Files ${mcc_tree} exists, skipping TreeAnnotator"
    echo "------------------------------------------------------------------------------"
    printf "\n"

fi

# ------------------------------------------------------------------------------

# Get clusters from MCC tree - cut at year and find clusters (discards single samples)

if [ ! -f ${beast_clusters_file} ]; then
    echo "------------------------------------------------------------------------------"
    echo "Running python_scripts/cut_tree.py on ${mcc_tree} - outputs ${beast_clusters_file}"
    set -x
    python python_scripts/cut_tree.py --infile ${mcc_tree} --outfile ${beast_clusters_file} --cut ${cut_year}
    set +x
    echo "------------------------------------------------------------------------------"
    printf "\n"
else
    echo "------------------------------------------------------------------------------"
    echo "Files ${beast_clusters_file} exists, skipping cut_tree.py"
    echo "------------------------------------------------------------------------------"
    printf "\n"
fi

# ------------------------------------------------------------------------------

# Transphylo

# if [ ! -f ${transphylo_es_table_file} ]; then
#     echo "------------------------------------------------------------------------------"
#     echo "Running TransPhylo on clusters output from cut_tree.py - outputs ${transphylo_es_table_file} and pdfs in ${transphylo_results_dir}"
#     set -x
#     transphylo.R -s ${study_accession} -t ${mcc_tree} -o ${transphylo_results_dir} -c ${beast_clusters_file} -m 1000000
#     set +x
#     echo "------------------------------------------------------------------------------"
#     printf "\n"
# else
#     echo "------------------------------------------------------------------------------"
#     echo "Files ${transphylo_es_table_file} exists, skipping transphylo.R"
#     echo "------------------------------------------------------------------------------"
#     printf "\n"
# fi

# Clean up
rm -r ${tmp_dir}


# TEST
