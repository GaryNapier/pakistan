







# ------------------------------------------------------------------------

# Chi-sq posthoc


trans_dr_chi_sq
trans_dr_chi_sq_ph

trans_loc_chi_sq
trans_loc_chi_sq_ph


# ------------------------------------------------------------------------

# Have a look at where transmission samples are on the tree


tree_all_samps <- midpoint.root(tree_all_samps)

# Tree setup
line_sz <- 0.25
n <- length(tree_all_samps$tip.label)
width <- 0.1
angle <- 45
font_sz <- 2

# Remove years from lineage data
lin_data_no_year <- lin_data
row.names(lin_data_no_year) <- unlist(lapply(strsplit(row.names(lin_data_no_year), "_"), function(x){ paste0(x[1:(length(x)-1)], collapse = "_") }))

# # Remove years from DR data 
# row.names(dr_status_data) <- unlist(lapply(strsplit(row.names(dr_status_data), "_"), function(x){ paste0(x[1:(length(x)-1)], collapse = "_") }))
# 
# # Remove years from location data
# row.names(loc_data) <- unlist(lapply(strsplit(row.names(loc_data), "_"), function(x){ paste0(x[1:(length(x)-1)], collapse = "_") }))

# Transmission wrangling
row.names(trans_non_trans_df) <- trans_non_trans_df$id
trans_non_trans_df <- drop_cols(trans_non_trans_df, "id")
trans_non_trans_df$trans_status <- as.factor(trans_non_trans_df$trans_status)


ggtree_all_samps <- ggtree(tree_all_samps, size = line_sz, layout="circular")

# Add lin data
lin_hm <- gheatmap(ggtree_all_samps, lin_data_no_year,
                   width = width,
                   # offset = offset,
                   color = NULL,
                   colnames_position = "top",
                   colnames_angle = angle, 
                   colnames_offset_y = 1,
                   hjust = 0,
                   font.size = 2) +
  scale_fill_manual(values = lin_colours, breaks = names(lin_colours) ) +
  labs(fill = "Lineage")+
  legend_spec

# Do this to add new gheatmap for some reason
# See "7.3.1 Visualize tree with multiple associated matrix" https://yulab-smu.top/treedata-book/chapter7.html
lin_hm <- lin_hm + ggnewscale::new_scale_fill()


gheatmap(lin_hm, trans_non_trans_df,
         offset = width*0.02,
         width = width,
         # color = NULL,
         low="white", high="black", color="black",
         colnames_position = "top",
         colnames_angle = angle, colnames_offset_y = 1,
         hjust = 0,
         font.size = font_sz,
         legend_title = "llw") +
  scale_fill_manual(values=c("white", "black"), labels = c("transmission", "non-transmission"))+
  labs(fill = "trans")+
  legend_spec



# ------------------------------------------------------------------------

# GWAS

gwas_dir <- "../gwas_results/"

snp_res_file <- paste0(gwas_dir, "trans_status.assoc.txt")
gene_res_file <- paste0(gwas_dir, "trans_status.genesum.assoc.txt")

snp_res <- read.delim(snp_res_file)
gene_res <- read.delim(gene_res_file) 

snp_res$log_p <- -log10(snp_res$p_wald)
gene_res$log_p <- -log10(gene_res$p_wald) 

snp_res$rs <- gsub("[a-z, A-Z, _]", "", snp_res$rs)

snp_res_sort <- snp_res[order(snp_res$log_p, decreasing = T), ]
gene_res_sort <- gene_res[order(gene_res$log_p, decreasing = T), ]

# Get p val of 10th snp/gene
snp_res_10 <- snp_res_sort[10, "log_p"]
gene_res_10 <- gene_res_sort[10, "log_p"]


plot(snp_res$log_p, pch = 19, cex = 0.2)
abline(h = snp_res_10, col = "red")
text(x = seq(snp_res$log_p), y = snp_res$log_p, labels = ifelse(snp_res$log_p >= snp_res_10, snp_res$rs, ""), offset = 0.5, cex=1, srt=-25)

plot(gene_res$log_p, pch = 19, cex = 0.2, ylab = "-log10 p-value", xlab = "Gene index")
abline(h = gene_res_10, col = "red")
text(x = seq(gene_res$log_p), y = gene_res$log_p, labels = ifelse(gene_res$log_p >= gene_res_10, gene_res$rs, ""), 
     adj = 0, cex=1, srt=-25)


# ------------------------------------------------------------------------


# Any transmissions involving >1 province?


in_trans_metadata

table(in_trans_metadata$province)
table(in_trans_metadata$location)


in_trans_metadata[, c("province", "location")]




table(metadata$sub_lineage)



# ------------------------------------------------------------------------

# TABLE 3 - join global freqs

global_file <- paste0(metadata_path, "novel_mutations_global_freqs.csv")

global <- read.csv(global_file)


names(global)[names(global) %in% "drug_resistant_."] <- "drug_resistant"


global[global[, "change"] == "c.-92T>G", ]
global[global[, "change"] == "p.Ala288Asp", ]

global_subset <- global[, c("drug", "change", 
                            "global_freq", 
                            "Sensitive_freq", 
                            "drug_resistant", 
                            "Pre.MDR_freq", 
                            "MDR_freq", 
                            "Pre.XDR_freq", 
                            "XDR_freq", 
                            "Other_freq")]

x <- merge(FN_results_table_pivot, global_subset, 
      by = c("drug", "change"), all.x = T, sort = F)

x <- x %>% mutate_if(is.numeric, round, 3)

x <- x[order(x$drug), ]

x <- x[, c("drug", 
      "gene", 
      "change",
      "N", 
      "global_freq", 
      "Sensitive_freq", 
      "drug_resistant", 
      "Pre.MDR_freq", 
      "MDR_freq", 
      "Pre.XDR_freq", 
      "XDR_freq", 
      "Other_freq")]

names(x) <- c("Drug", 
              "Gene", 
              "Change",
              "N", 
              "Global N", 
              "Sensitive freq", 
              "Drug resistant freq", 
              "Pre MDR freq", 
              "MDR freq", 
              "Pre XDR freq", 
              "XDR freq", 
              "Other freq")

x <- unique(x)

write.csv(x, paste0(metadata_path, "novel_mutations_and_global_freq_merged.csv"), quote = F, row.names = F)






























