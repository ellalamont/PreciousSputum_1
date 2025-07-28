# Batch corrected PCA on all the biological samples passing filter 
# E. Lamont
# 7/26/25

# Look into ggbiplot for more PCA stuff??
# https://cran.r-project.org/web/packages/ggbiplot/readme/README.html

# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

# Two options in base R, prcomp() and princomp()
# prcomp() is preferred according to the website above

source("Import_data.R") # To get All_RawReads


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
###################### PROCESS DATA #######################

pipeSummary2 <- All_pipeSummary %>% filter(SampleID %in% c("X", GoodSampleList, "H37Ra_Broth_4_S7", "H37Ra_Broth_5_S8", "H37Ra_Broth_6_S9"))

# Add Batch designation to the metadata:
pipeSummary2 <- pipeSummary2 %>%
  mutate(Batch = case_when(
    Run == "ProbeTest3" ~ 1,
    Run == "ProbeTest5" ~ 2,
    Run == "PredictTB_Run1" ~ 3
  ))

All_RawReads_2 <- All_RawReads %>%
  column_to_rownames(var = "X") %>%
  as.matrix()

# Need to make sure the order of the SampleIDs matches
pipeSummary2 <- pipeSummary2 %>%
  arrange(match(SampleID, colnames(All_RawReads_2)))


###########################################################
#################### BATCH CORRECTION #####################
# This basically all taken from Mark and I'm not sure what's happening
# https://academic.oup.com/nargab/article/2/3/lqaa078/5909519?login=true

batch <- pipeSummary2$Batch # Check this and make sure it is in the correct order!!
counts_corrected <- ComBat_seq(All_RawReads_2, batch = batch)


###########################################################
########################### PCA ###########################
# Just do the PCA the way I have been doing it

# Convert the batch corrected counts to counts per million (cpm)
my_cpm <- cpm(counts_corrected)

# transform the data 
my_cpm_t <- as.data.frame(t(my_cpm))

# Remove columns that are all zero so the scale works for prcomp
my_cpm_t <- my_cpm_t %>% select_if(colSums(.) != 0) # Down to 4494 genes 

# Make the actual PCA
my_PCA_cpm <- prcomp(my_cpm_t, scale = TRUE)

# See the % Variance explained
summary(my_PCA_cpm)
summary_PCA_cpm <- format(round(as.data.frame(summary(my_PCA_cpm)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA_cpm[1,1] # PC1 explains 21.1% of variance
summary_PCA_cpm[2,1] # PC2 explains 8.6% of variance
summary_PCA_cpm[3,1] # PC3 explains 7.8% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_cpm_df <- as.data.frame(my_PCA_cpm$x[, 1:3]) # Extract the first 3 PCs
my_PCA_cpm_df <- data.frame(SampleID = row.names(my_PCA_cpm_df), my_PCA_cpm_df)
my_PCA_cpm_df <- merge(my_PCA_cpm_df, All_pipeSummary, by = "SampleID", )

PCA_BatchCorrected <- my_PCA_cpm_df %>% 
  ggplot(aes(x = PC1, y = PC2, fill = Type, shape = Type)) + 
  geom_point(aes(fill = Type, shape = Type), size = 5, alpha = 0.7, stroke = 0.8) + 
  scale_fill_manual(values=c(`Week 0 sputum` = "#0072B2", `Caseum mimic` = "green4", `Broth`= "#999999", `Marmoset` = "#6A3D9A", `Rabbit` = "#E69F00")) +  
  scale_shape_manual(values=c(`Week 0 sputum` = 21, `Caseum mimic` = 22, `Broth`= 23, `Marmoset` = 24, `Rabbit` = 25)) + 
  # geom_text_repel(aes(label = Week), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "PCA: >1M reads and >80% genes with at least 10 reads",
       subtitle = "Batch corrected -> CPM",
       x = paste0("PC1: ", summary_PCA_cpm[1,1], "%"),
       y = paste0("PC2: ", summary_PCA_cpm[2,1], "%")) +
  my_plot_themes
PCA_BatchCorrected
ggsave(PCA_BatchCorrected,
       file = paste0("CPM_BatchCorrected_AllFilterSamples_1.pdf"),
       path = "Figures/PCA",
       width = 8, height = 5, units = "in")

# Add colour for relapse
PCA_BatchCorrected_2 <- my_PCA_cpm_df %>% 
  ggplot(aes(x = PC1, y = PC2, fill = Type, shape = Type, text = SampleID)) + 
  geom_point(aes(fill = Type2, shape = Type2), size = 5, alpha = 0.7, stroke = 0.8) + 
  scale_fill_manual(values=c(`W0 sputum (cure)` = "#0072B2",  `W0 sputum (relapse)` = "red", `Caseum mimic` = "green4", `Broth`= "#999999", `Marmoset` = "#6A3D9A", `Rabbit` = "#E69F00")) +  
  scale_shape_manual(values=c(`W0 sputum (cure)` = 21, `W0 sputum (relapse)` = 21,  `Caseum mimic` = 22, `Broth`= 23, `Marmoset` = 24, `Rabbit` = 25)) + 
  # geom_text_repel(aes(label = Week), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "PCA: >1M reads and >80% genes with at least 10 reads",
       subtitle = "Batch corrected -> CPM",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
PCA_BatchCorrected_2
ggsave(PCA_BatchCorrected_2,
       file = paste0("CPM_BatchCorrected_AllFilterSamples_2.pdf"),
       path = "Figures/PCA",
       width = 8, height = 5, units = "in")
ggplotly(fig2_PC1vsPC2)









