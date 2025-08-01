# Compare the same sample between different runs

source("Import_data.R") # To get Run1_tpm and ProbeTest5_tpm_subset


# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        # axis.text.x = element_text(angle = 45, size=14, vjust=1, hjust=1),
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
################ COLLECT DATA OF INTEREST #################

Run1_THP1Spike <- Run1_tpm %>% select(X, THP1_1e6_1_S67)
# ProbeTest5_tpm_subset

merged_THP1Spike <- merge(Run1_THP1Spike, ProbeTest5_tpm_subset, all = T) %>%
  rename(ProbeTest5_THP1.Ra1e6 = THP1_1e6_1a_S28,
         PredictTBRun1_THP1.Ra1e6 = THP1_1e6_1_S67, 
         Gene = X)

# Log10 transform the data
merged_THP1Spike_Log10 <- merged_THP1Spike %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values

# Remove all the non Rv genes and see how much better it gets
merged_THP1Spike_Log10_filtered <- merged_THP1Spike_Log10 %>% 
  filter(grepl("^Rv[0-9]+[A-Za-z]?$", Gene))

###########################################################
############# GGCORRPLOT OF THP1 1e6 SPIKED ###############

# Using all the genes
Sample1 <- "ProbeTest5_THP1.Ra1e6" # Captured
Sample2 <- "PredictTBRun1_THP1.Ra1e6" # Not Captured
ScatterCorr <- merged_THP1Spike_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("THP1 spiked with 1e6 on two separate runs"),
       subtitle = "Pearson correlation",
       x = paste0("Log2(TPM+1) ", Sample1), y = paste0("Log2(TPM+1) ", Sample2)) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = paste0("THP1Spiked1e6_CompareRuns_AllGenes.pdf"),
       path = "Figures/THP1Spiked_RunCompare",
       width = 7, height = 5, units = "in")

# Using only the Rv genes
Sample1 <- "ProbeTest5_THP1.Ra1e6" # Captured
Sample2 <- "PredictTBRun1_THP1.Ra1e6" # Not Captured
ScatterCorr <- merged_THP1Spike_Log10_filtered %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("THP1 spiked with 1e6 on two separate runs; only Rv genes"),
       subtitle = "Pearson correlation",
       x = paste0("Log2(TPM+1) ", Sample1), y = paste0("Log2(TPM+1) ", Sample2)) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
ggsave(ScatterCorr,
       file = paste0("THP1Spiked1e6_CompareRuns_RvGenes.pdf"),
       path = "Figures/THP1Spiked_RunCompare",
       width = 7, height = 5, units = "in")
