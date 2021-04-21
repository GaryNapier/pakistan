#!/usr/bin/env Rscript

# RUN:
# Rscript r_scripts/.R -s <study_acc>  -t <template_file_name>      -f <dated_fasta_file_name> -o <output_file_name>


# Checklist for Babette
# See XML file in (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7725332/)
# Upload XML file to Beauti to get parameters:

# Check if Babette is parsing years (Tip Dates tab > Use Tip Dates = T > Auto-configure > after last "_")
#  - Tip date file and create inference model
# create_inference_model
# tipdates_filename 
#  - MRCA prior
#  - Node at the crown age

# CLOCK MODEL
# - strict clock
# - clock rate = 0.0000001 (1.0E-7)

# PRIORS
# - Tree.t - Coalescent Constant Population
#   - Pop Size = 100 

# - popSize.t
#   - Log Normal
#   - Lower = 0
#   - Upper = 200
#   - Value = 100

# - Add Prior
#   - Taxon set label = TreePrior
#   - Distribution = Laplace Distribution
#   - Monophyletic = True
#   - Mu = [DATE] - HOW IS THIS DATE OBTAINED??

# SITE MODELS
# - GTR

# MANUAL ADDITIONS:
# beast2 correction for ascertainment bias, specifying the number of invariant A, C, G and T sites as 758511 1449901 1444524 758336. 
# <taxa id="TaxonSet.snps_london_351_dna" spec="TaxonSet">
#   <alignment id="snps_london_351_dna" spec="FilteredAlignment" filter="-">
#     <data idref="snps_london_351_dna_original"/>
#     <constantSiteWeights id="IntegerParameter.0" spec="parameter.IntegerParameter" dimension="4" lower="0" upper="0">758511 1449901 1444524 758336</constantSiteWeights>
#   </alignment>
# </taxa>

# MCMC
# Chain Length = 10000000
# tracelog
# - File Name = <accession><clust>.log
# - Log Every = 10000
# screenlog
# - File Name = <accession><clust>.screen
# - Log Every = 10000
# treelog
# - <accession><clust>.trees
# - Log Every = 10000


# Setup ----

rm(list=ls())

remotes::install_github("ropensci/beautier")

library(babette)
library(seqinr)
library(optparse)

heaD <- function(x,...){
  head(x, ...)
}

len_str <- function(string){
  length(unlist(strsplit(string, split = "")))
}

hs <- function(x, ...){
  print(head(x, ...))
  print("---")
  str(x, ...)
}


remove_tail <- function(x, sep = "_", del = 1){
  sapply(strsplit(x, split = sep, fixed = TRUE),
         function(i) paste(head(i, -del), collapse = sep))
}

# Arguments ----

option_list = list(
  make_option(c("-s", "--STUDY_ACCESSION"), type="character", default=NULL,
              help="input study accession number", metavar="character"),
  
  make_option(c("-f", "--dated_fasta_file_name"), type="character", default=NULL,
              help="input dated fasta file name", metavar="character"),
  
  make_option(c("-i", "--id_date_file_loc"), type="character", default="./",
              help="file location to save dataframe of ids and dates file to", metavar="character"),
  
  make_option(c("-o", "--output_file_name"), type="character", default="./",
              help="file location to save output xml file to", metavar="character"),
  
  make_option(c("-C", "--CLOCKRATE"), type="character", default="1.0E-7",
              help="input clockrate, default = 1.0E-7", metavar="character"),
  
  make_option(c("-F", "--FREQPARAMETER"), type="character", default="0.25",
              help="input freqparameter, default = 0.25", metavar="character"), 
  
  make_option(c("-E", "--every"), type="character", default="10000",
              help="input frequency of sampling - sample MCMC every x, default = 10000", metavar="character"), 
  
  make_option(c("-H", "--chain_length"), type="character", default="1e+08",
              help="input chainlength of MCMC, default = 1e+08", metavar="character"), 
  
  make_option(c("-M", "--mrca_prior_mu"), type="character", default="0",
              help="input prior for MRCA date of all samples, default = 0", metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print("ARGUMENTS:")
print(opt)
print("---")
print(str(opt))

# DIRECTORIES

# fasta_dir <- "fasta/"
# metadata_local_dir <- "metadata/"
# xml_dir <- "beast_xml/"

id_date_file_loc <- opt$id_date_file_loc

# FILES
study_accession <- opt$STUDY_ACCESSION # "PAKISTAN_ALL"
fasta_file <- opt$dated_fasta_file_name # paste0(fasta_dir, study_accession, ".dated.fa")
fasta_id_date_df_outfile <- paste0(id_date_file_loc, remove_tail(basename(fasta_file), sep = "."), ".txt")
xml_outfile <- opt$output_file_name # paste0(xml_dir, study_accession, ".xml")


# READ IN FILES
fasta <- seqinr::read.fasta(file = fasta_file, forceDNAtolower = F)

# Get fasta sample names and parse into name and date
fasta_names <- names(fasta)

fasta_id_date_df <- do.call(rbind, lapply(strsplit(fasta_names ,"_"), function(x){
  data.frame(id = paste(x[1:(length(x))], collapse = "_"), year = x[length(x)])
}))

# Save df of ids and year as csv for model
write.table(fasta_id_date_df, file = fasta_id_date_df_outfile, quote = F, row.names = F, col.names = F, sep = "\t")


# SITE MODEL
site_model <- create_gtr_site_model()


# CLOCK MODEL
# clock_rate <- 0.0000001
clock_rate <- opt$CLOCKRATE
clock_model <- create_strict_clock_model(clock_rate_param = create_clock_rate_param(value = clock_rate),
                                         clock_rate_distr = create_log_normal_distr(value = clock_rate, m = 1, s = 1.25))

# MCMC
# every <- 10000
every <- opt$every
chain_length <- opt$chain_length
mcmc <- create_mcmc(
  chain_length = chain_length,
  tracelog = beautier::create_tracelog(log_every = every),
  screenlog = beautier::create_screenlog(log_every = every),
  treelog = beautier::create_treelog(log_every = every)
)

# TREE PRIOR
tree_prior <- create_ccp_tree_prior(
  pop_size_distr = create_log_normal_distr(
    m = 1,
    s = 1.25,
    value = 100.0,
    lower = 0.0,
    upper = 200.0
  ))

# MRCA PRIOR
mrca_prior_mu <- opt$mrca_prior_mu
mrca_prior <- create_mrca_prior(
  is_monophyletic = TRUE, 
  mrca_distr = create_laplace_distr(mu = mrca_prior_mu), 
  name = "TreePrior")


# MAKE XML
create_beast2_input_file(
  fasta_file,
  xml_outfile,
  site_model = site_model,
  clock_model = clock_model,
  tree_prior = tree_prior,
  mrca_prior = mrca_prior,
  mcmc = mcmc,
  tipdates_filename = fasta_id_date_df_outfile, 
  beauti_options = beautier::create_beauti_options(nucleotides_uppercase = T)
)
























