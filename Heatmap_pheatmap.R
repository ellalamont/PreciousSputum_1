# Try a heatmap with pheatmap
# E. Lamont 
# 7/29/25

# https://rpubs.com/tgjohnst/heatmaps_testing_1

source("Import_data.R")

# Start with All_tpm, which is all the samples passing filter




###########################################################
###################### PROCESS DATA #######################

# Filter so there are fewer genes to deal with right now
# Lets filter for 100 just so it's easier for right now!
All_tpm_filtered100 <- All_tpm %>%
  filter(if_all(where(is.numeric), ~ .x >= 100))
# Now there are only 534 rows! 

# Change to a matrix?
# All_tpm_filtered100_2 <- All_tpm_filtered100 %>% column_to_rownames("X")
# All_tpm_matrix <- as.matrix(All_tpm_filtered100_2)

All_tpm2 <- All_tpm %>% column_to_rownames("X") %>%
  select(-contains("THP1"))

# Shorten the All_pipeSummary to see if it helps
pipeSummary_2 <- All_pipeSummary %>% filter(SampleID %in% c(SputumSampleList, BrothSampleList, MarmSampleList, MimicSampleList, RabbitSampleList))

my_annotation_colors <- list(
  Type2 = c("W0 sputum (cure)" = "#0072B2",
            "W0 sputum (relapse)" = "red", 
            "Caseum mimic" = "green4",
            "Marmoset" = "#6A3D9A", 
            "Rabbit" = "#E69F00", 
            "Broth" = "#999999")
)

###########################################################
######################## PHEATMAP #########################

pheatmap(All_tpm_matrix[1:10,], scale = "row")

# Need to shorten the names first
names(All_tpm2) <- gsub(x = names(All_tpm2), pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)
pipeSummary_3 <- pipeSummary_2 %>%
  mutate(SampleID = gsub("_S.*", "", SampleID))

testing <- All_tpm2 %>% subset(rownames(All_tpm2) %in% allGeneSetList[["MTb.TB.Phenotypes.TopGeneSets"]][["human_sputum: top 25 genes"]]) # Guess this doesn't need to be a matrix

pheatmap(testing, scale = "row")



# Grab the columns needed to give colors
Color_annotation_df <- pipeSummary_3 %>%
  filter(SampleID %in% colnames(testing)) %>%
  select(SampleID, Type2) %>%
  column_to_rownames("SampleID")

# Reorder annotation rows to match columns of tpm file
Color_annotation_df <- Color_annotation_df[colnames(testing), , drop = FALSE]

# Define the colors
my_annotation_colors <- list(
  Type2 = c("W0 sputum (cure)" = "#0072B2",
            "W0 sputum (relapse)" = "red", 
            "Caseum mimic" = "green4",
            "Marmoset" = "#6A3D9A", 
            "Rabbit" = "#E69F00", 
            "Broth" = "#999999")
)

pheatmap(testing, 
         annotation_col = annotation_df, 
         annotation_colors = my_annotation_colors,
         scale = "row")





###########################################################
######################## ALL DATA #########################

pheatmap(my_tpm , 
         # annotation_col = my_pipeSummary["Week"], 
         # annotation_row = gene_annot["Product"],
         annotation_colors = my_annotation_colors,
         scale = "row")

my_tpm_2 <- my_tpm[rowSums(my_tpm == 0) != ncol(my_tpm), ]

my_tpm_2_matrix <- my_tpm_2 %>% 
  # rename("W0_250754" = "S_250754",
  #        "W0_355466" = "S_355466",
  #        "W0_503557" = "S_503557",
  #        "W2_503937" = "S_503937",
  #        "W2_575533" = "S_575533_MtbrRNA",
  #        "W2_577208" = "S_577208") %>%
  as.matrix()

testing2 <- testing %>% 
  as.matrix()


pheatmap(testing2, 
         annotation_col = my_pipeSummary["Week"], 
         scale = "row",
         cutree_rows = 5,
         cutree_cols = 5)


###########################################################
############### ALL DATA WITH CLUSTERING ##################
# https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/

pheatmap(my_tpm_2_matrix, 
         # annotation_row = my_gene_col, 
         annotation_col = my_pipeSummary["Week"],
         scale = "row",
         cutree_rows = 5,
         cutree_cols = 5)

# Try to get good row annotations based on MTb functional group
# Start with MTb.TB.Phenotypes.AllGeneSets
# Convert to dataframe
Gene_Category <- do.call(rbind, lapply(names(Walter2015GeneSets), function(category) {
  data.frame(Gene = Walter2015GeneSets[[category]], Category = category, stringsAsFactors = FALSE)
})) %>% 
  filter(!Category %in% c("Cluster A", "Cluster B", "Cluster D", "Cluster E")) %>% 
  distinct(Gene, .keep_all = TRUE) %>% # Keep only the first gene occurance
  column_to_rownames(var = "Gene")


pheatmap(my_tpm_2_matrix, 
         # annotation_row = Gene_Category, 
         fontsize_row = 1,
         annotation_col = my_pipeSummary["Week"],
         annotation_colors = my_annotation_colors,
         scale = "row",
         cutree_rows = 8,
         cutree_cols = 5)





# Pull out what the clusters are: 
# use silent = TRUE to suppress the plot
my_heatmap <- pheatmap(my_tpm_2_matrix, 
                       annotation_row = Gene_Category, 
                       annotation_col = my_pipeSummary["Week"],
                       annotation_colors = my_annotation_colors,
                       scale = "row",
                       cutree_rows = 7,
                       cutree_cols = 5,
                       silent = TRUE)
# Extract the row clustering information
row_clusters <- cutree(my_heatmap$tree_row, k = 7)
# Convert to a data frame for easier handling
row_cluster_df <- data.frame(Gene = names(row_clusters), Cluster = row_clusters)
# View the first few rows
head(row_cluster_df)


###########################################################
############ W0 AND BROTH WITH CLUSTERING #################

my_tpm_3_matrix <- my_tpm_2 %>% select(-c("S_503937", "S_577208", "S_575533_MtbrRNA")) %>%
  as.matrix()

Gene_Category <- do.call(rbind, lapply(names(MTb.TB.Phenotypes.AllGeneSets), function(category) {
  data.frame(Gene = MTb.TB.Phenotypes.AllGeneSets[[category]], Category = category, stringsAsFactors = FALSE)
})) %>% 
  # filter(!Category %in% c("Cluster A", "Cluster B", "Cluster D", "Cluster E")) %>% 
  distinct(Gene, .keep_all = TRUE) %>% # Keep only the first gene occurance
  column_to_rownames(var = "Gene")

pheatmap(my_tpm_3_matrix, 
         annotation_row = Gene_Category, 
         fontsize_row = 1,
         annotation_col = my_pipeSummary["Week"],
         annotation_colors = my_annotation_colors,
         scale = "row",
         cutree_rows = 6,
         cutree_cols = 2)











###########################################################
################### TESTING FOR SHINY #####################


allGeneSetList[["MTb.TB.Phenotypes.TopGeneSets"]][["microaerophilic: top 25 genes"]]
my_data <- my_tpm %>% subset(rownames(my_tpm) %in% allGeneSetList[["MTb.TB.Phenotypes.TopGeneSets"]][["human_sputum: top 25 genes"]])
p <- pheatmap(my_data, 
              annotation_col = my_pipeSummary["Week"], 
              annotation_colors = my_annotation_colors,
              scale = "row")
p
heatmap(as.matrix(my_data))


selected_genes <- c("Rv0081", "Rv0494", "Rv2011c", "Rv1473A")
my_data <- my_tpm[rownames(my_tpm) %in% selected_genes, , drop = FALSE]
pheatmap(my_data, 
              annotation_col = my_pipeSummary["Week"], 
              annotation_row = gene_annot["Product"],  # Conditional annotation
              annotation_colors = my_annotation_colors,
              scale = "row", 
              fontsize = 18)



