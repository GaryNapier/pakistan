#!/usr/bin/env Rscript

# Finds samples from global metadata that belong to same sublineages as pakistan data. 
# Removes countries with less than n samples and downsamples these.
# Have to read in a country code lookup file to clean the countries.
# Outputs a list of samples to be passed to vcf to make fasta to make a tree.
# Outputs the metadata of the global samples so that the later tree can be annotated.

# RUN:
# pull_global_samps.R -m <metadata_file> -c <country_code_lookup_file> -s <sample_list_outfile> -g <global_metadata_outfile>

# Setup ----

rm(list=ls())

library(optparse)

# Arguments ----

option_list = list(
  make_option(c("-m", "--metadata_file"), type="character", default=NULL,
              help="input metadata file name", metavar="character"),
  make_option(c("-c", "--country_code_lookup_file"), type="character", default=NULL,
              help="input metadata file name", metavar="character"),
  make_option(c("-s", "--sample_list_outfile"), type="character", default=NULL,
              help="input outfile name", metavar="character"),
  make_option(c("-g", "--global_metadata_outfile"), type="character", default=NULL,
              help="input outfile name", metavar="character"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print("ARGUMENTS:")
print(opt)
print("---")
print(str(opt))

# FILES

metadata_34k_file <- opt$metadata_file
country_code_lookup_file <- opt$country_code_lookup_file
sample_list_outfile <- opt$sample_list_outfile
global_metadata_outfile <- opt$global_metadata_outfile


# rm(list=ls())
# methods_path <- "../methods/"
# main_metadata_path <- "../../metadata/"
# metadata_path <- "../metadata/"
# newick_path <- "../newick/"
# plots_path <- "../plots/"
# db_path <- "../../pipeline/db/"
# dist_path <- "../dist_and_pca/"
# beast_results_path <- "../beast_results/"
# beast_test_path <- "../beast_test/"
# metadata_34k_file <- paste0(main_metadata_path, "tb_data_18_02_2021.csv")
# country_code_lookup_file <- paste0(main_metadata_path, "country_code_lookup.txt")
# sample_list_outfile <- "sample_list_outfile.txt"
# global_metadata_outfile <- "global_metadata_outfile.csv"


# READ IN FILES

metadata_34k <- read.csv(metadata_34k_file, header = T, stringsAsFactors = F)
country_code_lookup <- read.delim(country_code_lookup_file, na.strings = "")

# CLEAN

# Strip 'lineage' from sub_lineage col
metadata_34k$sub_lineage <- gsub("lineage", "", metadata_34k$sub_lineage)

# Drop NA country codes
metadata_34k <- metadata_34k[!(metadata_34k[, "country_code"] == "N/A"), ]

# Merge in clean country names
metadata_34k <- merge(metadata_34k, country_code_lookup, by.x = "country_code", by.y = "country_code_lower", all.x = T, sort = F)
# Replace country
metadata_34k$country <- metadata_34k$country.y 
# Drop the two old country cols
metadata_34k <- metadata_34k[, !(names(metadata_34k) %in% c("country.x", "country.y"))]

# Get unique lineages from pakistan data
uniq_lins <- unique(metadata_34k[metadata_34k["country_code"] == "pk", "sub_lineage"])

# Get the global samples that are the same sublin as the pakisatan samples
global_samps <- metadata_34k[metadata_34k$sub_lineage %in% uniq_lins, ]

# Remove the pakistan samples
# global_samps <- global_samps[!(global_samps$country_code == "pk"), ]

# DO ANALYSIS

n <- 10

# Remove countries that have less than 20 samples
country_tab <- table(global_samps$country_code)
country_keep <- names(subset(country_tab, country_tab > n)) # Keep these
# Remove NAs
country_keep <- country_keep[!(country_keep == "N/A")]

# Check
# x <- unique(global_samps[global_samps$country_code %in% country_keep, c("country_code", "country", "region")])

# Remove countries not in list
global_samps <- global_samps[global_samps$country_code %in% country_keep, ]

# Loop over the split data and downsample to ten rows per sublin
global_samps_split <- split(global_samps, global_samps$country_code)
global_samps_downsamp <- list()
for(i in seq(global_samps_split)){
  x <- global_samps_split[[i]]
  if(length(x$country_code) < n){
    next
  }else{
    set.seed(1)
    global_samps_downsamp[[i]] <- x[sample(1:nrow(x), n), ]
  }
}
global_samps_downsamp <- do.call("rbind", global_samps_downsamp)

# Check
# x <- unique(global_samps_downsamp[global_samps_downsamp$country_code %in% country_keep, 
#                                   c("country_code", "country", "subregion", "region")])

# WRITE OUT

# Write out metadata
# Put togehter with pakistan metadata
write.csv(global_samps_downsamp, global_metadata_outfile, row.names = F, quote = F)

# Take the sample names
# Put together with pakistan samples
write.table(global_samps_downsamp$wgs_id, sample_list_outfile, sep = "\t", row.names = F, col.names = F, quote = F)






















  