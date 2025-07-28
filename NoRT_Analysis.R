# Sequenced some sputum sample with no RT (reverse transcriptase). Compare these to the samples sequenced normally
# E. Lamont
# 7/28/25

source("Import_data.R") # To get All_pipeSummary

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
####################### N_Genomic #########################

RTvsNot_N_Genomic <- All_pipeSummary %>%
  filter(Patient %in% c("P_12008", "P_12019", "P_12070", "P_13045")) %>% 
  ggplot(aes(x = RT, y = N_Genomic)) + 
  geom_point(size = 4, alpha = 0.8, stroke = 0.8, color = "black") + 
  geom_line(aes(group = Patient), color = "black", size = 0.5, linetype = "dashed") + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 3, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  scale_y_continuous(limits = c(0,12000000), breaks = seq(0, 12000000, 2000000)) +
  labs(title = "W0 Sputum with or without reverse transcriptase",
       subtitle = NULL, 
       x = NULL, 
       y = "# reads aligning to H37Rv") + 
  # scale_x_discrete(labels = c("Yes" = "With RT",
                              # "No" = "Without RT")) + 
  my_plot_themes
RTvsNot_N_Genomic
ggsave(RTvsNot_N_Genomic,
       file = paste0("RTvsNot_N_Genomic.pdf"),
       path = "Figures/NoRT_Analysis",
       width = 7, height = 5, units = "in")


###########################################################
####################### P_Genomic #########################

RTvsNot_P_Genomic <- All_pipeSummary %>%
  filter(Patient %in% c("P_12008", "P_12019", "P_12070", "P_13045")) %>% 
  ggplot(aes(x = RT, y = P_Genomic)) + 
  geom_point(size = 4, alpha = 0.8, stroke = 0.8, color = "black") + 
  geom_line(aes(group = Patient), color = "black", size = 0.5, linetype = "dashed") + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 3, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) + 
  labs(title = "W0 Sputum with or without reverse transcriptase",
       subtitle = NULL, 
       x = NULL, 
       y = "# reads aligning to H37Rv") + 
  # scale_x_discrete(labels = c("Yes" = "With RT",
  # "No" = "Without RT")) + 
  my_plot_themes
RTvsNot_P_Genomic
ggsave(RTvsNot_P_Genomic,
       file = paste0("RTvsNot_P_Genomic.pdf"),
       path = "Figures/NoRT_Analysis",
       width = 7, height = 5, units = "in")

###########################################################
#################### AtLeast.10.Reads #####################

RTvsNot_TenReads <- All_pipeSummary %>% 
  filter(Patient %in% c("P_12008", "P_12019", "P_12070", "P_13045")) %>% 
  ggplot(aes(x = RT, y = AtLeast.10.Reads)) + 
  # geom_point(aes(fill = Patient), shape = 21, alpha = 0.8, size = 4) + 
  geom_point(size = 4, alpha = 0.8, stroke = 0.8, color = "black") + 
  geom_line(aes(group = Patient), color = "black", size = 0.5, linetype = "dashed") + 
  geom_hline(yintercept = 4499*0.8, linetype = "dashed", alpha = 0.5) + 
  annotate("text", x = 0.9, y = 4499*0.8, label = "80%", hjust = 4, vjust = -0.5, color = "black") + 
  geom_hline(yintercept = 4499*0.5, linetype = "dashed", alpha = 0.5) + 
  annotate("text", x = 0.9, y = 4499*0.5, label = "50%", hjust = 4, vjust = -0.5, color = "black") + 
  labs(title = "W0 Sputum with or without reverse transcriptase: Genes with >= 10 reads aligning for all sample types",
       subtitle = NULL, 
       x = NULL, 
       y = "# of genes with at least 10 reads aligning") + 
  scale_y_continuous(limits = c(0,4500), breaks = seq(0, 4500, 500)) + 
  # scale_x_continuous(trans = "log10") + 
  my_plot_themes
RTvsNot_TenReads
ggsave(RTvsNot_TenReads,
       file = paste0("RTvsNot_TenReads.pdf"),
       path = "Figures/NoRT_Analysis",
       width = 7, height = 5, units = "in")



