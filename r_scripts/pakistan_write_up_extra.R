







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

gene_res_sort[1:5, c("rs", "beta", "log_p")]


plot(snp_res$log_p, pch = 19, cex = 0.2)
abline(h = snp_res_10, col = "red")
text(x = seq(snp_res$log_p), y = snp_res$log_p, labels = ifelse(snp_res$log_p >= snp_res_10, snp_res$rs, ""), offset = 0.5, cex=1, srt=-25)

plot(gene_res$log_p, pch = 19, cex = 0.2, ylab = "-log10 p-value", xlab = "Location")
abline(h = gene_res_10, col = "red")
text(x = seq(gene_res$log_p), y = gene_res$log_p, labels = ifelse(gene_res$log_p >= gene_res_10, gene_res$rs, ""), 
     adj = 0, cex=1, srt=-25)

# Put on tree

genesum_file <- "../gwas_results/PAKISTAN_ALL.filt.val.gt.g.ann.genesum"
samples_file <- "../gwas_results/samples.txt"

genesum <- read.delim(genesum_file, header = F)
samples <- read.delim(samples_file, header = F)

genesum <- genesum[, -c(2, 3)]

names(genesum) <- c("gene", samples[,1])

top_genes <- gene_res_sort[1:5, "rs"]

genesum_sub <- genesum[genesum[, "gene"] %in% top_genes, ]

#         nusG
# ERR123  0
# ERR321  0
# ERR456  1

genesum_list <- list()
for(i in seq(nrow(genesum_sub))){
  gene_df <- genesum_sub[i, ]
  df <- setNames(data.frame(as.numeric(as.vector(gene_df[2:ncol(gene_df)]))), gene_df["gene"])
  row.names(df) <- samples[, 1]
  genesum_list[[i]] <- df
}

genesum_list <- lapply(genesum_list, function(x) {
  apply(x, 2, as.factor)
})

legend_spec <- theme(legend.title = element_text(size = 10),
                     legend.text = element_text(size = 10),
                     legend.key.size = unit(0.3, "cm"))

ggtree_all_samps <- ggtree(tree_all_samps, size = line_sz, layout="fan", open.angle = 10)

trans_tree <- gheatmap(ggtree_all_samps, trans_non_trans_df,
         width = width,
         offset = 0,
         low="white", high="black", color="black",
         colnames_position = "top",
         colnames_angle = angle, colnames_offset_y = 1,
         hjust = 0,
         font.size = font_sz,
         legend_title = "Transmission status") +
  labs(fill = "Transmission status")+
  scale_fill_manual(values=c("white", "black"), labels = c("Transmission", "Non-transmission"))+
  legend_spec

trans_tree <- trans_tree + ggnewscale::new_scale_fill()

df <- genesum_list[[1]]
gene <- colnames(df)

gene_tree <- gheatmap(trans_tree, df,
                       width = width,
                       offset = 0.003,
                       low="white", high="black", color="black",
                       colnames_position = "top",
                       colnames_angle = angle, colnames_offset_y = 1,
                       hjust = 0,
                       font.size = font_sz,
                       legend_title = gene) +
  labs(fill = gene)+
  scale_fill_gradient2()+
  scale_fill_manual(values=c("white", "black"), labels = c("0", "1"))+
  legend_spec

gene_tree <- gene_tree + ggnewscale::new_scale_fill()

df <- genesum_list[[2]]
gene <- colnames(df)

gene_tree <- gheatmap(gene_tree, df,
                      width = width,
                      offset = 0.006,
                      low="white", high="black", color="black",
                      colnames_position = "top",
                      colnames_angle = angle, colnames_offset_y = 1,
                      hjust = 0,
                      font.size = font_sz,
                      legend_title = gene) +
  labs(fill = gene)+
  scale_fill_gradient2()+
  scale_fill_manual(values=c("white", "black"), labels = c("0", "1"))+
  legend_spec

gene_tree <- gene_tree + ggnewscale::new_scale_fill()

df <- genesum_list[[3]]
gene <- colnames(df)

gene_tree <- gheatmap(gene_tree, df,
                      width = width,
                      offset = 0.009,
                      low="white", high="black", color="black",
                      colnames_position = "top",
                      colnames_angle = angle, colnames_offset_y = 1,
                      hjust = 0,
                      font.size = font_sz,
                      legend_title = gene) +
  labs(fill = gene)+
  scale_fill_gradient2()+
  scale_fill_manual(values=c("white", "black"), labels = c("0", "1"))+
  legend_spec

gene_tree <- gene_tree + ggnewscale::new_scale_fill()

df <- genesum_list[[4]]
gene <- colnames(df)

gene_tree <- gheatmap(gene_tree, df,
                      width = width,
                      offset = 0.012,
                      low="white", high="black", color="black",
                      colnames_position = "top",
                      colnames_angle = angle, colnames_offset_y = 1,
                      hjust = 0,
                      font.size = font_sz,
                      legend_title = gene) +
  labs(fill = gene)+
  scale_fill_gradient2()+
  scale_fill_manual(values=c("white", "black"), labels = c("0", "1"))+
  legend_spec

gene_tree <- gene_tree + ggnewscale::new_scale_fill()

df <- genesum_list[[5]]
gene <- colnames(df)

gheatmap(gene_tree, df,
                      width = width,
                      offset = 0.015,
                      low="white", high="black", color="black",
                      colnames_position = "top",
                      colnames_angle = angle, colnames_offset_y = 1,
                      hjust = 0,
                      font.size = font_sz,
                      legend_title = gene) +
  labs(fill = gene)+
  scale_fill_gradient2()+
  scale_fill_manual(values=c("white", "black"), labels = c("0", "1"))+
  legend_spec




genesum <- as.data.frame(do.call("cbind", genesum_list))

genesum$wgs_id <- row.names(genesum)

df <- merge(metadata[, c("wgs_id", "trans_status")], genesum, by = "wgs_id", all.x = T, sort = F)

table(df$trans_status, df$nusG)
table(df$trans_status, df$Rv0914c)
table(df$trans_status, df$Rv1896c)
table(df$trans_status, df$Rv2102)
table(df$trans_status, df$Rv2184c)



# Put nusG mutations on tree

nusg_mut <- read.delim("../gwas_results/nusG_samples_mutations.txt", header = F)

names(nusg_mut) <- c("wgs_id", "mut")

nusg_mut <- merge(metadata["wgs_id"], nusg_mut, by = "wgs_id", all.x = T, sort = F)

nusg_cols <- rainbow(length(unique(nusg_mut$mut)))

nusg_cols[length(nusg_cols)] <- "#FFFFFF"

names(nusg_cols) <- unique(nusg_mut$mut)

row.names(nusg_mut) <- nusg_mut$wgs_id

nusg_mut <- nusg_mut[!(names(nusg_mut) %in% "wgs_id")]

ggtree_all_samps <- ggtree(tree_all_samps, size = line_sz, layout="fan", open.angle = 10)

trans_tree <- gheatmap(ggtree_all_samps, trans_non_trans_df,
                       width = width,
                       offset = 0,
                       low="white", high="black", color="black",
                       colnames_position = "top",
                       colnames_angle = angle, colnames_offset_y = 1,
                       hjust = 0,
                       font.size = 3,
                       legend_title = "Transmission status") +
  labs(fill = "Transmission status")+
  scale_fill_manual(values=c("white", "black"), labels = c("Transmission", "Non-transmission"))+
  legend_spec

trans_tree <- trans_tree + ggnewscale::new_scale_fill()

png(filename = "../plots/nusG.png", width = 2000, height = 2000, res = 300)
gheatmap(trans_tree, nusg_mut,
         width = width,
         offset = 0.003,
         color = NULL,
         colnames_position = "top",
         colnames_angle = angle, 
         colnames_offset_y = 1,
         hjust = 0,
         font.size = 3) +
  scale_fill_manual(values = nusg_cols, breaks = names(nusg_cols) ) +
  labs(fill = "nusG mutations")+
  legend_spec
dev.off()




# ------------------------------------------------------------------------


# Any transmissions involving >1 province?


in_trans_metadata

table(in_trans_metadata$province)
table(in_trans_metadata$location)


in_trans_metadata[, c("province", "location")]




table(metadata$sub_lineage)


# ------------------------------------------------------------------------

# Join global freqs to known mutations table (dr_variants_pivot)

# Read in global data
global_file <- paste0(metadata_path, "mutations_global_freqs.csv")
global <- read.csv(global_file)

# Clean
names(global)[names(global) %in% "drug_resistant_."] <- "drug_resistant_pc"

# Checks/have a look
global[global[, "change"] == "c.-92T>G", ]
global[global[, "change"] == "p.Ala288Asp", ]

# Re-arrange cols, take out gene col
global_subset <- global[, c("drug", "change", 
                            "global_freq", 
                            "Sensitive_freq", 
                            "drug_resistant_pc", 
                            "Pre.MDR_freq", 
                            "MDR_freq", 
                            "Pre.XDR_freq", 
                            "XDR_freq", 
                            "Other_freq")]

# Merge
known_mut_global <- merge(dr_variants_pivot, global_subset, by = c("drug", "change"), all.x = T, sort = F)

# Clean
known_mut_global <- known_mut_global %>% mutate_if(is.numeric, round, 3)
known_mut_global <- known_mut_global[order(known_mut_global$drug), ]
known_mut_global <- known_mut_global[, c("drug", 
           "gene", 
           "change",
           "N", 
           "global_freq", 
           "Sensitive_freq", 
           "drug_resistant_pc", 
           "Pre.MDR_freq", 
           "MDR_freq", 
           "Pre.XDR_freq", 
           "XDR_freq", 
           "Other_freq")]

names(known_mut_global) <- c("Drug", 
              "Gene", 
              "Change",
              "N", 
              "Global N", 
              "Sensitive freq", 
              "Drug resistant %", 
              "Pre MDR freq", 
              "MDR freq", 
              "Pre XDR freq", 
              "XDR freq", 
              "Other freq")

known_mut_global <- data.frame(lapply(known_mut_global, function(x) {gsub(",", ";", x)}))

write.csv(known_mut_global, file = paste0(metadata_path, "known_mutations_and_global_freq_merged.csv"),
          row.names = F, quote = F)


# ------------------------------------------------------------------------

# Join global freqs to novel mutations table (FN_results_table_pivot)

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
      "drug_resistant_pc", 
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



none_global <- x[x[, "Global N"] == 0, ]

length(unique(none_global$Change))


# ------------------------------------------------------------------------


















