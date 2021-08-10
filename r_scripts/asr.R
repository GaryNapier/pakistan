
library(colorspace)

# PATHS

rm(list=ls())
methods_path <- "../methods/"
main_metadata_path <- "../../metadata/"
metadata_path <- "../metadata/"
newick_path <- "../newick/"
plots_path <- "../plots/"

# FILES

metadata_34k_file <- paste0(main_metadata_path, "tb_data_18_02_2021.csv")
country_code_lookup_file <- paste0(main_metadata_path, "country_code_lookup.txt")
sample_list_outfile <- "sample_list_outfile.txt"
global_metadata_file <- paste0(metadata_path, "global_metadata.txt")
global_tree_file <- paste0(newick_path, "GLOBAL.filt.val.gt.g.snps.fa.treefile")
plot_file <- paste0(plots_path, "asr_tree.png")


# READ IN FILES

global_metadata <- read.delim(global_metadata_file, header = T, stringsAsFactors = F)
global_tree <- read.tree(global_tree_file)


# CLEAN

# Midpoint root the tree
global_tree <- phytools::midpoint.root(global_tree)


# COLOURS.. here we go.. 

# Get df of unique regions and subregions
country_region <- unique(global_metadata[, c("subregion", "region")])
# Sort
country_region <- country_region[order(country_region$region, country_region$region), ]
# Split on main region
country_region_split <- split(country_region, country_region$region)
# Get the unique main regions
regions <- sort(unique(country_region$region))
# Get the main colours for the main regions
region_cols <- rainbow(length(regions), alpha = 0.9)
# Define function to loop over the split df and lighten the colours by the number of subreegions
get_cols_split <- function(splits){
  
  top_cols <- rainbow(length(splits))
  top_cols <- darken(top_cols, 0.3, space = "HLS")
  show_col(top_cols)
  
  all_cols <- vector()
  for(split in seq(splits)){
    print(split)
    print(splits[[split]])
    
    # len_split <- length(splits[[split]])
    len_split <- nrow(splits[[split]])
    split_top_col <- top_cols[split]
    if(len_split == 1){
      all_cols <- append(all_cols, split_top_col)
    }else{
      lighten_to_col <- colorspace::lighten(split_top_col, 1-(1/len_split))
      fc <- colorRampPalette(c(split_top_col, lighten_to_col))
      split_cols <- fc(len_split)
      all_cols <- append(all_cols, split_cols)
    }
  }
  # names(all_cols) <- unlist(splits)
  # cols_df <- data.frame(lin = unlist(splits), col = all_cols)
  return(all_cols)
}
# Apply function
all_cols <- get_cols_split(country_region_split)
# Check
show_col(all_cols)
# Add colours to df
country_region$col <- all_cols
# Add Pakistan
country_region[nrow(country_region)+1, ] <- c("Pakistan", "Asia", "#808080")
# Make pakistan grey
# country_region[country_region[, "country"] == "Pakistan", "col"] <- "#808080"


# http://www.phytools.org/eqg2015/asr.html

# Pull countries as vector and add sample names 
# countries <- global_metadata[, "country"]
# names(countries) <- global_metadata[, "wgs_id"]
x <- global_metadata[, "subregion"] # 'x' is used as the variable of the ASR in the examples on the webpage
names(x) <- global_metadata[, "wgs_id"]

# Get pk samps and change value
pk_samps <- global_metadata[global_metadata[, "country_code"] == "pk", "wgs_id"]
x[names(x) %in% pk_samps] <- "Pakistan"

# # Same for colours
# # cols <- country_region$col
# # names(cols) <- country_region$country
# cols <- rainbow(length(unique(x)), alpha = 0.9)
# cols <- colorspace::darken(cols)
# names(cols) <- sort(unique(x))
# # Change pakistan col to grey
# cols[names(cols) == "Pakistan"] <- "#808080"

cols <- country_region[, "col"]
names(cols) <- country_region[, "subregion"]

# fitER <- ace(x, anoletree, model="ER",type="discrete")
# fitER <- ace(countries, global_tree, model="ER",type="discrete")
fitER <- ace(x, global_tree, model="ER",type="discrete")

## simulate single stochastic character map using empirical Bayes method
mtree <- make.simmap(global_tree, x, model="ER")

mtree$tip.label <- rep("", length(mtree$tip.label))

# plotTree(anoletree,type="fan",fsize=0.8,ftype="i")
# plot(global_tree, show.tip.label = F)
png(filename = plot_file, width = 2000, height = 2000, res = 300)
plot(mtree, cols, lwd = 1)
nodelabels(node=1:global_tree$Nnode+Ntip(global_tree), pie = fitER$lik.anc, piecol = cols, cex = 0.3)
add.simmap.legend(colors = cols, prompt = F, fsize = 0.8, y = 550)
dev.off()










