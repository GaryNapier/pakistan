#!/usr/bin/env Rscript

# RUN:
# Rscript r_scripts/clean_metadata.R <main_metadata_file>  <other_metadata_file>            <outfile>

# Setup ----

args <- commandArgs(trailingOnly=TRUE)

library(plyr)
library(dplyr)
library(stringr)

heaD <- function(x, ...){
  head(x, ...)
}

clear <- function(){rm(list=ls(envir = .GlobalEnv),envir = .GlobalEnv); }

options(stringsAsFactors = F)

# Paths ----

# metadata_path <- "~/Documents/metadata/"
# output_path <- "~/Documents/pakistan/metadata/"

# Files ----

# metadata_file <- paste0(metadata_path, "tb_data_18_02_2021.csv")
# pakistan_data_file <- paste0(metadata_path, "pakistan_data_non_mixed.csv")
# pakistan_data_outfile <- paste0(output_path, "pakistan_metadata.csv")

metadata_file <- args[1]
pakistan_data_file <- args[2]
pakistan_data_outfile <- args[3]

# Read in data ----

metadata <- read.csv(metadata_file)
pakistan_data <- read.csv(pakistan_data_file)

# Subset ----

metadata <- subset(metadata, country_code == "pk")


# Merge in Pakistan unpublished data ----

# Rename cols
pakistan_data <- dplyr::rename(pakistan_data, genotypic_drtype = drtype)

# Add in country code and country columns:
pakistan_data$country_code <- rep("pk", nrow(pakistan_data))
pakistan_data$country <- rep("Pakistan", nrow(pakistan_data))

# rbind together both datasets - n.b. plyr::rbind.fill function
metadata <- plyr::rbind.fill(metadata, pakistan_data)


# Study accession codes ----

metadata$study_accession_word <- rep("PAKISTAN_ALL", nrow(metadata))


# Drug resistance ----

# Rename para.aminosalicylic_acid
metadata <- dplyr::rename(metadata, para_aminosalicylic_acid = para.aminosalicylic_acid)

# Classifications

# https://www.sciencedirect.com/topics/agricultural-and-biological-sciences/prothionamide

# First line
fl <- c("rifampicin", "isoniazid", "ethambutol", "pyrazinamide", "streptomycin")
# Second Line Injectable
sli <- c("amikacin", "kanamycin", "capreomycin")
# Fluoroquinolones (second line)
flq <- c("ofloxacin", "moxifloxacin", "levofloxacin", "ciprofloxacin")
# Other second line
osl <- c("prothionamide", "ethionamide", "para_aminosalicylic_acid", "cycloserine")
# All second line
sl <- c("ofloxacin", "moxifloxacin", "levofloxacin", "amikacin", "kanamycin", "capreomycin", "ciprofloxacin", "prothionamide", 
        "ethionamide", "para_aminosalicylic_acid", "cycloserine")
# Third line
tl <- c("clarithromycin", "clofazimine", "bedaquiline", "linezolid", "rifabutin", "delamanid")

# DR tests
# dr_test <- function(x){any(x[fl] == 1, na.rm = T)}
# mdr_test <- function(x){ifelse(is.na(x["rifampicin"]) || is.na(x["isoniazid"]), 
#                                FALSE, 
#                                x["rifampicin"] == 1 && x["isoniazid"] == 1)}
# xdr_test <- function(x){mdr_test(x) && any(x[sl] == 1, na.rm = T)}

# Susceptible: Nothing
sus_test <- function(x){all(x[c("rifampicin", "isoniazid", sli, flq)] == 0, na.rm = T)}
# Pre-MDR: rif OR inh
pre_mdr_test <- function(x){ifelse(is.na(x["rifampicin"]) || is.na(x["isoniazid"]), FALSE, x["rifampicin"] == 1 || x["isoniazid"] == 1)}
# MDR: rif AND inh
mdr_test <- function(x){
  ifelse(is.na(x["rifampicin"]) || is.na(x["isoniazid"]), FALSE, x["rifampicin"] == 1 && x["isoniazid"] == 1)
}
# Pre-XDR: rif AND inh AND (SLI OR flq)
pre_xdr_test <- function(x){
  mdr_test(x) && 
    ifelse( (all(is.na(x[sli])) && all(is.na(x[flq]))), FALSE, 
            (any(x[sli] == 1, na.rm = T) || any(x[flq] == 1, na.rm = T)) )
}
# XDR: rif AND inh AND sli AND flq
xdr_test <- function(x){
  mdr_test(x) && 
    ifelse( (all(is.na(x[sli])) && all(is.na(x[flq]))), FALSE, 
            (any(x[sli] == 1, na.rm = T) && any(x[flq] == 1, na.rm = T)) )
}
# Other: Not any of the above

# Put all tests together in if/else
resistance_tests <- function(x){
  drugs <- c("rifampicin", "isoniazid", "ethambutol", "pyrazinamide", "streptomycin", 
             "ofloxacin", "moxifloxacin", "levofloxacin", "amikacin", "kanamycin", "capreomycin", "ciprofloxacin", "prothionamide", 
             "ethionamide", "para_aminosalicylic_acid", "cycloserine",
             "clarithromycin", "clofazimine", "bedaquiline", "linezolid", "rifabutin", "delamanid")
  ifelse(all(is.na(x[drugs])), NA, 
         ifelse(xdr_test(x), "XDR",  
                ifelse(pre_xdr_test(x), "Pre-XDR", 
                       ifelse(mdr_test(x), "MDR", 
                              ifelse(pre_mdr_test(x), "Pre-MDR", 
                                     ifelse(sus_test(x), "Sensitive", "Other") ) ) ) ) )
}

# Loop through and assign test. Can't get apply to work for some bloody reason. 
dr_status <- vector()
for(i in 1:nrow(metadata)){
  dr_status[i] <- resistance_tests(metadata[i, ])
}

# cbind to data
metadata$dr_status <- dr_status

# Clean DR status - bring across genotypic if na in the clinical data
metadata$dr_status <- ifelse(is.na(metadata$dr_status), metadata$genotypic_drtype, metadata$dr_status)


# Lineages ---- 

# Make SURE there are no mixed samples coming through
metadata <- metadata[grep(";", metadata$main_lineage, invert = T), ]
metadata <- metadata[grep(";", metadata$sub_lineage, invert = T), ]

# Make sure there are no blanks in the lineage columns
metadata <- metadata[!(metadata$main_lineage == ""), ]
metadata <- metadata[!(metadata$sub_lineage == ""), ]

# Make sure bov/cap/ory main lineages are correct
metadata$main_lineage <- ifelse(metadata$main_lineage == "M", metadata$sub_lineage, metadata$main_lineage)

# Remove 'lineage' from main_lineage and sub_lineage
metadata$main_lineage <- sub("lineage", "", metadata$main_lineage)
metadata$sub_lineage <- sub("lineage", "", metadata$sub_lineage)

# Make "major lineage" column - first decimal place of sublineage
metadata$major_lineage <- substr(metadata$sub_lineage, 1, 3)

# Clean animal strains
animal <- c("M.bovis", "M.caprae", "M.orygis")
metadata$major_lineage <- ifelse(metadata$sub_lineage %in% animal, metadata$sub_lineage, metadata$major_lineage)


# Dates ----

# Add year column
metadata$year <- stringr::str_extract(metadata$collection_date, pattern = "\\d{4}")

# Misc. clean ----

# Replace all commas with semi-colon
metadata <- data.frame(lapply(metadata, function(x) {gsub(",", ";", x)}))

# Write data ----

write.csv(metadata, file = pakistan_data_outfile, quote = F, row.names = F)











