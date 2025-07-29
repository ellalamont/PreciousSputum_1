# PCA on all the biological samples passing filter
# E. Lamont
# 7/26/25

# Look into ggbiplot for more PCA stuff??
# https://cran.r-project.org/web/packages/ggbiplot/readme/README.html

# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

# Two options in base R, prcomp() and princomp()
# prcomp() is preferred according to the website above



source("Import_data.R") # To get All_tpm


# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=14),
        # legend.title = element_text(size = 14),
        legend.title = element_blank(),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=9), 
        plot.margin = margin(10, 10, 10, 20)# ,
        # panel.background = element_rect(fill='transparent'),
        # plot.background = element_rect(fill='transparent', color=NA),
        # legend.background = element_rect(fill='transparent'),
        # legend.box.background = element_blank()
  )



###########################################################
################ ALL PASSING FILTER PCA ###################
# Passing filter is >1,000,000 genomic reads and >80% genes with at least 10 reads

# Convert gene column to rownames
All_tpm2 <- All_tpm %>% column_to_rownames(var = "X")

# Transform the data
All_tpm_t <- as.data.frame(t(All_tpm2))

# Remove columns that are all zero so the scale works for prcomp
All_tpm_t2 <- All_tpm_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(All_tpm_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 23.3% of variance
summary_PCA[2,1] # PC2 explains 8.2% of variance
summary_PCA[3,1] # PC3 explains 7.1% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, All_pipeSummary, by = "SampleID", )

PCA_tpm_1 <- my_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2, fill = Type, shape = Type)) + 
  geom_point(aes(fill = Type, shape = Type), size = 5, alpha = 0.7, stroke = 0.8) + 
  scale_fill_manual(values=c(`Week 0 sputum` = "#0072B2", `Caseum mimic` = "green4", `Broth`= "#999999", `Marmoset` = "#6A3D9A", `Rabbit` = "#E69F00")) +  
  scale_shape_manual(values=c(`Week 0 sputum` = 21, `Caseum mimic` = 22, `Broth`= 23, `Marmoset` = 24, `Rabbit` = 25)) + 
  # geom_text_repel(aes(label = Week), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "PCA: >1M reads and >80% genes with at least 10 reads",
       subtitle = "TPM",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
PCA_tpm_1
ggsave(PCA_tpm_1,
       file = paste0("TPM_AllFilterSamples_1.pdf"),
       path = "Figures/PCA",
       width = 8, height = 5, units = "in")

# Add colour for relapse
PCA_tpm_2 <- my_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2, fill = Type, shape = Type, text = SampleID)) + 
  geom_point(aes(fill = Type2, shape = Type2), size = 5, alpha = 0.7, stroke = 0.8) + 
  scale_fill_manual(values=c(`W0 sputum (cure)` = "#0072B2",  `W0 sputum (relapse)` = "red", `Caseum mimic` = "green4", `Broth`= "#999999", `Marmoset` = "#6A3D9A", `Rabbit` = "#E69F00")) +  
  scale_shape_manual(values=c(`W0 sputum (cure)` = 21, `W0 sputum (relapse)` = 21,  `Caseum mimic` = 22, `Broth`= 23, `Marmoset` = 24, `Rabbit` = 25)) + 
  # geom_text_repel(aes(label = Week), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "PCA: >1M reads and >80% genes with at least 10 reads",
       subtitle = "TPM",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
PCA_tpm_2
ggsave(PCA_tpm_2,
       file = paste0("TPM_AllFilterSamples_2.pdf"),
       path = "Figures/PCA",
       width = 8, height = 5, units = "in")

# 3D plot
# https://plotly.com/r/pca-visualization/
PCA_3D <- plot_ly(my_PCA_df, x = ~PC1, y = ~PC2, z = ~PC3,
                  type = "scatter3d", mode = "markers",
                  color = ~Type2# , 
                  # colors = c12,
                  # text = ~Replicate
)
PCA_3D



###########################################################
############## PCA JUST SPUTUM AND BROTH ##################

tpm_subset <- All_tpm %>% select("X", all_of(SputumSampleList), all_of(BrothSampleList))

# Convert gene column to rownames
my_tpm <- tpm_subset %>% column_to_rownames(var = "X")

# Transform the data
my_tpm_t <- as.data.frame(t(my_tpm))

# Remove columns that are all zero so the scale works for prcomp
my_tpm_t2 <- my_tpm_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(my_tpm_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 10.9% of variance
summary_PCA[2,1] # PC2 explains 8.4% of variance
summary_PCA[3,1] # PC3 explains 7.9% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, All_pipeSummary, by = "SampleID", )

PCA_tpm_1 <- my_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2, fill = Type, shape = Type)) + 
  geom_point(aes(fill = Type, shape = Type), size = 5, alpha = 0.7, stroke = 0.8) + 
  scale_fill_manual(values=c(`Week 0 sputum` = "#0072B2", `Caseum mimic` = "green4", `Broth`= "#999999", `Marmoset` = "#6A3D9A", `Rabbit` = "#E69F00")) +  
  scale_shape_manual(values=c(`Week 0 sputum` = 21, `Caseum mimic` = 22, `Broth`= 23, `Marmoset` = 24, `Rabbit` = 25)) + 
  # geom_text_repel(aes(label = Week), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "PCA: >1M reads and >80% genes with at least 10 reads",
       subtitle = "TPM",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
PCA_tpm_1
ggsave(PCA_tpm_1,
       file = paste0("TPM_Sputum.Broth_1.pdf"),
       path = "Figures/PCA",
       width = 8, height = 5, units = "in")

# Add colour for relapse
PCA_tpm_2 <- my_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2, fill = Type, shape = Type, text = SampleID)) + 
  geom_point(aes(fill = Type2, shape = Type2), size = 5, alpha = 0.7, stroke = 0.8) + 
  scale_fill_manual(values=c(`W0 sputum (cure)` = "#0072B2",  `W0 sputum (relapse)` = "red", `Caseum mimic` = "green4", `Broth`= "#999999", `Marmoset` = "#6A3D9A", `Rabbit` = "#E69F00")) +  
  scale_shape_manual(values=c(`W0 sputum (cure)` = 21, `W0 sputum (relapse)` = 21,  `Caseum mimic` = 22, `Broth`= 23, `Marmoset` = 24, `Rabbit` = 25)) + 
  # geom_text_repel(aes(label = Week), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "PCA: >1M reads and >80% genes with at least 10 reads",
       subtitle = "TPM",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
PCA_tpm_2
ggsave(PCA_tpm_2,
       file = paste0("TPM_Sputum.Broth_2.pdf"),
       path = "Figures/PCA",
       width = 8, height = 5, units = "in")
ggplotly(fig2_PC1vsPC2)

# 3D plot
# https://plotly.com/r/pca-visualization/
PCA_3D <- plot_ly(my_PCA_df, x = ~PC1, y = ~PC2, z = ~PC3,
                  type = "scatter3d", mode = "markers",
                  color = ~Type2# , 
                  # colors = c12,
                  # text = ~Replicate
)
PCA_3D


###########################################################
################### PCA JUST SPUTUM #######################

tpm_subset <- All_tpm %>% select("X", all_of(SputumSampleList))

# Convert gene column to rownames
my_tpm <- tpm_subset %>% column_to_rownames(var = "X")

# Transform the data
my_tpm_t <- as.data.frame(t(my_tpm))

# Remove columns that are all zero so the scale works for prcomp
my_tpm_t2 <- my_tpm_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(my_tpm_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 11.7% of variance
summary_PCA[2,1] # PC2 explains 8.9% of variance
summary_PCA[3,1] # PC3 explains 8.1% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, All_pipeSummary, by = "SampleID", )

PCA_tpm_1 <- my_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2, fill = Type, shape = Type)) + 
  geom_point(aes(fill = Type, shape = Type), size = 5, alpha = 0.7, stroke = 0.8) + 
  scale_fill_manual(values=c(`Week 0 sputum` = "#0072B2", `Caseum mimic` = "green4", `Broth`= "#999999", `Marmoset` = "#6A3D9A", `Rabbit` = "#E69F00")) +  
  scale_shape_manual(values=c(`Week 0 sputum` = 21, `Caseum mimic` = 22, `Broth`= 23, `Marmoset` = 24, `Rabbit` = 25)) + 
  # geom_text_repel(aes(label = Week), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "PCA: >1M reads and >80% genes with at least 10 reads",
       subtitle = "TPM",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
PCA_tpm_1
ggsave(PCA_tpm_1,
       file = paste0("TPM_Sputum_1.pdf"),
       path = "Figures/PCA",
       width = 8, height = 5, units = "in")

# Add colour for relapse
PCA_tpm_2 <- my_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2, fill = Type, shape = Type, text = SampleID)) + 
  geom_point(aes(fill = Type2, shape = Type2), size = 5, alpha = 0.7, stroke = 0.8) + 
  scale_fill_manual(values=c(`W0 sputum (cure)` = "#0072B2",  `W0 sputum (relapse)` = "red", `Caseum mimic` = "green4", `Broth`= "#999999", `Marmoset` = "#6A3D9A", `Rabbit` = "#E69F00")) +  
  scale_shape_manual(values=c(`W0 sputum (cure)` = 21, `W0 sputum (relapse)` = 21,  `Caseum mimic` = 22, `Broth`= 23, `Marmoset` = 24, `Rabbit` = 25)) + 
  geom_text_repel(aes(label = Patient), size= 2.5, box.padding = 0.4, max.overlaps = Inf) + 
  labs(title = "PCA: >1M reads and >80% genes with at least 10 reads",
       subtitle = "TPM",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
PCA_tpm_2
ggsave(PCA_tpm_2,
       file = paste0("TPM_Sputum_2.pdf"),
       path = "Figures/PCA",
       width = 8, height = 5, units = "in")
ggplotly(fig2_PC1vsPC2)

# 3D plot
# https://plotly.com/r/pca-visualization/
PCA_3D <- plot_ly(my_PCA_df, x = ~PC1, y = ~PC2, z = ~PC3,
                  type = "scatter3d", mode = "markers",
                  color = ~Type2# , 
                  # colors = c12,
                  # text = ~Replicate
)
PCA_3D




###########################################################
############## ALL SAMPLES ONLY Rv# GENES #################
# See how the PCA changes when I only include Rv## genes

# Remove all the non Rv genes
All_tpm_filtered <- All_tpm %>% 
  filter(grepl("^Rv[0-9]+[A-Za-z]?$", X))

# Convert gene column to rownames
All_tpm2 <- All_tpm_filtered %>% column_to_rownames(var = "X")

# Transform the data
All_tpm_t <- as.data.frame(t(All_tpm2))

# Remove columns that are all zero so the scale works for prcomp
All_tpm_t2 <- All_tpm_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(All_tpm_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 22.9% of variance
summary_PCA[2,1] # PC2 explains 8.5% of variance
summary_PCA[3,1] # PC3 explains 7.4% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, All_pipeSummary, by = "SampleID", )

PCA_tpm_1 <- my_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2, fill = Type, shape = Type)) + 
  geom_point(aes(fill = Type, shape = Type), size = 5, alpha = 0.7, stroke = 0.8) + 
  scale_fill_manual(values=c(`Week 0 sputum` = "#0072B2", `Caseum mimic` = "green4", `Broth`= "#999999", `Marmoset` = "#6A3D9A", `Rabbit` = "#E69F00")) +  
  scale_shape_manual(values=c(`Week 0 sputum` = 21, `Caseum mimic` = 22, `Broth`= 23, `Marmoset` = 24, `Rabbit` = 25)) + 
  # geom_text_repel(aes(label = Week), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "PCA: >1M reads and >80% genes with at least 10 reads",
       subtitle = "TPM, only Rv## Gene! Doesn't look different",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
PCA_tpm_1
ggsave(PCA_tpm_1,
       file = paste0("TPM_AllFilterSamples_3.pdf"),
       path = "Figures/PCA",
       width = 8, height = 5, units = "in")



