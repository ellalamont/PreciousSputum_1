# Look at where the hits are in the no RT controls
# E. Lamont
# 7/31/25

# If hits are random across different samples, likely DNA contamination
# If hits are always at specific places, likely misalignment
# Should Also pull in data from THP1 only samples from ProbeTest 5


source("Import_data.R")

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "null",legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
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

# Get THP1 controls from ProbeTest5_tpm
THP1_tpm <- ProbeTest5_tpm %>% select("X", contains("control"))

# get noRT reads from Run1_tpm
noRT_tpm <- Run1_tpm %>% select("X", contains("noRT"))



###########################################################
################### GRAPH OF noRT SAMPLES #################

noRT_fig <- noRT_tpm %>%
  # head(20) %>% 
  pivot_longer(cols = starts_with("W0"), names_to = "Sample", values_to = "TPM") %>% 
  mutate(TPM_log10 = log10(TPM + 1)) %>% 
  ggplot(aes(x = X, y = TPM_log10, fill = Sample), text = X) + # for log10 graph
  # ggplot(aes(x = X, y = TPM, fill = Sample), text = X) + # for normal graph
  geom_bar(position="stack", stat="identity") +  # position = stack for raw reads, position = fill for propotion
  scale_fill_manual(values = cbPalette4) + 
  labs(title = "Genes hits for 4 noRT sputum samples",
       subtitle = "big gene is Rv2503c", 
       y = "log10(TPM + 1)",
       # y = "TPM", 
       x = "Gene (alphabetically ordered)") + 
  scale_y_continuous(expand = c(0,0)) + 
  my_plot_themes
# noRT_fig
# ggplotly(noRT_fig)
ggsave(noRT_fig,
       file = paste0("noRT_TPMLog10.pdf"),
       path = "Figures/Hits_in_noRT",
       width = 20, height = 5, units = "in")
# ggsave(noRT_fig,
#        file = paste0("noRT_TPM.pdf"),
#        path = "Figures/Hits_in_noRT",
#        width = 20, height = 5, units = "in")

###########################################################
################### GRAPH OF THP1 SAMPLES #################

THP1_fig <- THP1_tpm %>%
  # head(20) %>% 
  pivot_longer(cols = starts_with("THP1"), names_to = "Sample", values_to = "TPM") %>% 
  mutate(TPM_log10 = log10(TPM + 1)) %>% 
  ggplot(aes(x = X, y = TPM_log10, fill = Sample), text = X) + # for log10 graph
  # ggplot(aes(x = X, y = TPM, fill = Sample), text = X) + # for normal graph
  geom_bar(position="stack", stat="identity") +  # position = stack for raw reads, position = fill for propotion
  scale_fill_manual(values = cbPalette4) + 
  labs(title = "Genes hits for 4 THP1 controls (not spiked)",
       subtitle = "big gene is Rv1398c", 
       y = "log10(TPM + 1)",
       # y = "TPM", 
       x = "Gene (alphabetically ordered)") + 
  scale_y_continuous(expand = c(0,0)) + 
  my_plot_themes
# THP1_fig
# ggplotly(THP1_fig)
# ggsave(THP1_fig,
#        file = paste0("THP1_fig_TPMLog10.pdf"),
#        path = "Figures/Hits_in_noRT",
#        width = 20, height = 5, units = "in")
# ggsave(THP1_fig,
#        file = paste0("THP1Control_TPM.pdf"),
#        path = "Figures/Hits_in_noRT",
#        width = 20, height = 5, units = "in")






