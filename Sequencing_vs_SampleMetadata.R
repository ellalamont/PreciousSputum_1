# Compare Sputum sequencing results with the same metadata

# source("Import_data.R")
source("Import_SampleMetadata.R")

# Week 2 has been excluded from all below

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

# Stop scientific notation
# options(scipen = 999) 
options(scipen = 0) # To revert back to default

my_regression_line <- stat_poly_line(method = "lm", se = TRUE, level = 0.95, color = "black", alpha = 0.3)
my_regression_equations <- stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                                          after_stat(rr.label),
                                                          after_stat(p.value.label), 
                                                          sep = "*\", \"*")))

###########################################################
####################### MERGE METADATA ####################

# Starting with All_pipeSummary and my_metadata

# my_metadata is already a subset of the columns that I want to work on now

my_metadata_2 <- my_metadata %>% 
  filter(Patient %in% All_pipeSummary$Patient) %>%
  filter(Visit == "Day 0")

pipeSummary_2 <- All_pipeSummary %>%
  filter(Type == "Week 0 sputum") %>%
  select(SampleID, total_RNA_ng, Week, RawReads, N_RiboClear, P_RiboClear, N_Genomic, P_Genomic, N_NoHit, P_NoHit, AtLeast.10.Reads, AtLeast.100.Reads, Run, Type, Type2, SampleID2, Patient, Outcome)


merged_metadata <- merge(my_metadata_2, pipeSummary_2, by = "Patient", all = T)
  


###########################################################
#################### N+P GENOMIC TO Ct ####################

# N_GENOMIC
ctVsReads <- merged_metadata %>% 
  filter(XpertCT_wk0 != "NA") %>% 
  mutate(XpertCT_wk0 = as.numeric(XpertCT_wk0)) %>% 
  ggplot(aes(x = XpertCT_wk0, y = N_Genomic)) +
  geom_point(size = 4, alpha = 0.8, stroke = 0.8, color = "black") + 
  # geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  scale_y_continuous(limits = c(-3000000,12000000), breaks = seq(-3000000,12000000, 2000000), labels = scales::scientific ) +
  scale_x_continuous(limits = c(12,31), breaks = seq(12,31,2), expand = c(0,0)) +
  labs(title = "W0 Sputum Run 1: Ct value vs number reads aligned to Mtb",
       subtitle = NULL,
       y = "# reads aligning to Mtb genome") + 
  my_plot_themes + 
  stat_poly_line(method = "lm", se = TRUE, level = 0.95, color = "grey23", alpha = 0.25) + 
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label),
                                 after_stat(p.value.label),
                                 sep = "*\", \"*")),
               label.x = "right", label.y = "top", parse = T)
ctVsReads # + my_regression_line + my_regression_equations
ggsave(ctVsReads,
       file = paste0("ctVsReads_1.pdf"),
       path = "Figures/Sputum_Sequencing_vs_SampleMetadata",
       width = 7, height = 5, units = "in")


# P_GENOMIC
ctVsPercent <- merged_metadata %>% 
  filter(XpertCT_wk0 != "NA") %>% 
  mutate(XpertCT_wk0 = as.numeric(XpertCT_wk0)) %>% 
  ggplot(aes(x = XpertCT_wk0, y = P_Genomic)) +
  geom_point(size = 4, alpha = 0.8, stroke = 0.8, color = "black") + 
  # geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  geom_hline(yintercept = 0, linetype = "solid", alpha = 1) + 
  scale_y_continuous(limits = c(-30,100), breaks = seq(-30, 100, 10)) +
  scale_x_continuous(limits = c(12,31), breaks = seq(12,31,2), expand = c(0,0)) +
  labs(title = "W0 Sputum Run 1: Ct value vs percent reads aligned to Mtb",
       subtitle = NULL,
       y = "% reads aligning to Mtb genome") + 
  my_plot_themes + 
  stat_poly_line(method = "lm", se = TRUE, level = 0.95, color = "grey23", alpha = 0.25) + 
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label),
                                 after_stat(p.value.label),
                                 sep = "*\", \"*")),
               label.x = "right", label.y = "top", parse = T)
ctVsPercent # + my_regression_line + my_regression_equations
ggsave(ctVsPercent,
       file = paste0("ctVsPercent_1.pdf"),
       path = "Figures/Sputum_Sequencing_vs_SampleMetadata",
       width = 7, height = 5, units = "in")

###########################################################
################### N+P GENOMIC TO TTD ####################

# N_GENOMIC
ttdVsReads <- merged_metadata %>% 
  filter(TTD != "NA") %>% 
  mutate(TTD = as.numeric(TTD)) %>% 
  ggplot(aes(x = TTD, y = N_Genomic)) +
  geom_point(size = 4, alpha = 0.8, stroke = 0.8, color = "black") + 
  # geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  scale_y_continuous(limits = c(-8000000,12000000), breaks = seq(-8000000,12000000, 2000000), labels = scales::scientific ) +
  scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2), expand = c(0,0)) +
  labs(title = "W0 Sputum Run 1: TTD vs number reads aligned to Mtb",
       subtitle = NULL,
       y = "# reads aligning to Mtb genome") + 
  my_plot_themes + 
  stat_poly_line(method = "lm", se = TRUE, level = 0.95, color = "grey23", alpha = 0.25) + 
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label),
                                 after_stat(p.value.label),
                                 sep = "*\", \"*")),
               label.x = "right", label.y = "top", parse = T)
ttdVsReads # + my_regression_line + my_regression_equations
ggsave(ttdVsReads,
       file = paste0("ttdVsReads_1.pdf"),
       path = "Figures/Sputum_Sequencing_vs_SampleMetadata",
       width = 7, height = 5, units = "in")

# P_GENOMIC
ttdVsPercent <- merged_metadata %>% 
  filter(TTD != "NA") %>% 
  mutate(TTD = as.numeric(TTD)) %>% 
  ggplot(aes(x = TTD, y = P_Genomic)) +
  geom_point(size = 4, alpha = 0.8, stroke = 0.8, color = "black") + 
  # geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  geom_hline(yintercept = 0, linetype = "solid", alpha = 1) + 
  scale_y_continuous(limits = c(-60,100), breaks = seq(-60, 100, 10)) +
  scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2), expand = c(0,0)) +
  labs(title = "W0 Sputum Run 1: TTD vs percent reads aligned to Mtb",
       subtitle = NULL,
       y = "percent reads aligning to Mtb genome") + 
  my_plot_themes + 
  stat_poly_line(method = "lm", se = TRUE, level = 0.95, color = "grey23", alpha = 0.25) + 
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label),
                                 after_stat(p.value.label),
                                 sep = "*\", \"*")),
               label.x = "right", label.y = "top", parse = T)
ttdVsPercent # + my_regression_line + my_regression_equations
ggsave(ttdVsPercent,
       file = paste0("ttdVsPercent_1.pdf"),
       path = "Figures/Sputum_Sequencing_vs_SampleMetadata",
       width = 7, height = 5, units = "in")



###########################################################
############### N+P GENOMIC TO Total RNA ng ###############

# N_GENOMIC
ngVsReads <- merged_metadata %>% 
  filter(total_RNA_ng != "NA") %>% 
  mutate(total_RNA_ng = as.numeric(total_RNA_ng)) %>% 
  ggplot(aes(x = total_RNA_ng, y = N_Genomic)) +
  geom_point(size = 4, alpha = 0.8, stroke = 0.8, color = "black") + 
  # geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  scale_y_continuous(limits = c(0,12000000), breaks = seq(0,12000000, 2000000), labels = scales::scientific ) +
  # scale_x_continuous(limits = c(190,30000), breaks = c(200, 1000, 10000, 30000), trans = "log10") +
  scale_x_log10() + 
  labs(title = "W0 Sputum Run 1: total RNA ng vs number reads aligned to Mtb",
       subtitle = NULL,
       x = "Log10(total RNA ng)",
       y = "# reads aligning to Mtb genome") + 
  my_plot_themes + 
  stat_poly_line(method = "lm", se = TRUE, level = 0.95, color = "grey23", alpha = 0.25) + 
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label),
                                 after_stat(p.value.label),
                                 sep = "*\", \"*")),
               label.x = "right", label.y = "top", parse = T)
ngVsReads # + my_regression_line + my_regression_equations
# ggplotly(ngVsReads)
ggsave(ngVsReads,
       file = paste0("ngVsReads_1.pdf"),
       path = "Figures/Sputum_Sequencing_vs_SampleMetadata",
       width = 7, height = 5, units = "in")

# P_GENOMIC
ngVsPercent <- merged_metadata %>% 
  filter(total_RNA_ng != "NA") %>% 
  mutate(total_RNA_ng = as.numeric(total_RNA_ng)) %>% 
  ggplot(aes(x = total_RNA_ng, y = P_Genomic)) +
  geom_point(size = 4, alpha = 0.8, stroke = 0.8, color = "black") + 
  # geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) +
  # scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2), expand = c(0,0)) +
  scale_x_log10() + 
  labs(title = "W0 Sputum Run 1: total RNA ng vs percent reads aligned to Mtb",
       subtitle = NULL,
       x = "Log10(total RNA ng)",
       y = "percent reads aligning to Mtb genome") + 
  my_plot_themes + 
  stat_poly_line(method = "lm", se = TRUE, level = 0.95, color = "grey23", alpha = 0.25) + 
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label),
                                 after_stat(p.value.label),
                                 sep = "*\", \"*")),
               label.x = "right", label.y = "top", parse = T)
ngVsPercent # + my_regression_line + my_regression_equations
ggsave(ngVsPercent,
       file = paste0("ngVsPercent_1.pdf"),
       path = "Figures/Sputum_Sequencing_vs_SampleMetadata",
       width = 7, height = 5, units = "in")





