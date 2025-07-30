# Compare the reads to the CFUs of the rabbit, marmoset, and mimic samples
# E. Lamont
# 7/30/25


source("Import.data.R") # To get All_pipeSummary


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
###################### CFU vs GENOMIC #####################

# N_GENOMIC
CFU.Reads <- All_pipeSummary %>% 
  filter(Type %in% c("Marmoset", "Caseum mimic", "Rabbit")) %>% 
  ggplot(aes(x = CFU_per_g.or.mL, y = N_Genomic)) +
  geom_point(aes(fill = Type, shape = Type), size = 4, alpha = 0.8, stroke = 0.8, color = "black") + 
  scale_fill_manual(values=c(`Caseum mimic` = "green4", `Marmoset` = "#6A3D9A", `Rabbit` = "#E69F00")) +  
  scale_shape_manual(values=c(`Caseum mimic` = 22, `Marmoset` = 24, `Rabbit` = 25)) + 
  # geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  scale_y_continuous(limits = c(0,10000000), breaks = seq(0, 10000000, 1000000)) +
  # scale_x_continuous(limits = c(12,31), breaks = seq(12,31,2), expand = c(0,0)) +
  scale_x_continuous(trans = "log10") + 
  labs(title = "Marmoset, rabbit, caseum mimic: CFU/g (/mL for mimic) vs number reads aligned to Mtb",
       subtitle = "Marmoset data frome ProbeTest3; two rabbits are replicates from same lesion",
       x = "log10(CFU/g or /mL)", 
       y = "# reads aligning to Mtb genome") + 
  my_plot_themes
CFU.Reads # + my_regression_line + my_regression_equations
ggsave(CFU.Reads,
       file = paste0("CFU.Reads_1.pdf"),
       path = "Figures/Rabbit.Marm.Mimc",
       width = 8, height = 5, units = "in")

# P_GENOMIC
CFU.Percent <- All_pipeSummary %>% 
  filter(Type %in% c("Marmoset", "Caseum mimic", "Rabbit")) %>% 
  ggplot(aes(x = CFU_per_g.or.mL, y = P_Genomic)) +
  geom_point(aes(fill = Type, shape = Type), size = 4, alpha = 0.8, stroke = 0.8, color = "black") + 
  scale_fill_manual(values=c(`Caseum mimic` = "green4", `Marmoset` = "#6A3D9A", `Rabbit` = "#E69F00")) +  
  scale_shape_manual(values=c(`Caseum mimic` = 22, `Marmoset` = 24, `Rabbit` = 25)) + 
  # geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) +
  # scale_x_continuous(limits = c(12,31), breaks = seq(12,31,2), expand = c(0,0)) +
  scale_x_continuous(trans = "log10") + 
  labs(title = "Marmoset, rabbit, caseum mimic: CFU/g (/mL for mimic) vs percent reads aligned to Mtb",
       subtitle = "Marmoset data frome ProbeTest3; two rabbits are replicates from same lesion",
       x = "log10(CFU/g or /mL)", 
       y = "% reads aligning to Mtb genome") + 
  my_plot_themes
CFU.Percent # + my_regression_line + my_regression_equations
ggsave(CFU.Percent,
       file = paste0("CFU.Percent_1.pdf"),
       path = "Figures/Rabbit.Marm.Mimc",
       width = 8, height = 5, units = "in")


###########################################################
################## CFU vs AT LEAST 10 READS ###############

CFU.10reads <- All_pipeSummary %>% 
  filter(Type %in% c("Marmoset", "Caseum mimic", "Rabbit")) %>% 
  ggplot(aes(x = CFU_per_g.or.mL, y = AtLeast.10.Reads)) + 
  geom_point(aes(fill = Type, shape = Type), size = 4, alpha = 0.8, stroke = 0.8, color = "black") + 
  scale_fill_manual(values=c(`Caseum mimic` = "green4", `Marmoset` = "#6A3D9A", `Rabbit` = "#E69F00")) +  
  scale_shape_manual(values=c(`Caseum mimic` = 22, `Marmoset` = 24, `Rabbit` = 25)) + 
  geom_hline(yintercept = 4499*0.8, linetype = "dashed", alpha = 0.5) + 
  annotate("text", x = 50000, y = 4499*0.8, label = "80%", hjust = 0.5, vjust = -0.5, color = "black") + 
  geom_hline(yintercept = 4499*0.5, linetype = "dashed", alpha = 0.5) + 
  annotate("text", x = 50000, y = 4499*0.5, label = "50%", hjust = 0.5, vjust = -0.5, color = "black") + 
  labs(title = "Marmoset, rabbit, caseum mimic: Genes with >= 10 reads aligning",
       subtitle = "Marmoset data frome ProbeTest3; two rabbits are replicates from same lesion", 
       x = "log10(CFU/g or /mL)", 
       y = "# of genes with at least 10 reads aligning") + 
  scale_y_continuous(limits = c(0,4500), breaks = seq(0, 4500, 500)) + 
  scale_x_continuous(trans = "log10") +
  my_plot_themes
CFU.10reads
ggsave(CFU.10reads,
       file = paste0("CFU.10reads_1.pdf"),
       path = "Figures/Rabbit.Marm.Mimc",
       width = 8, height = 5, units = "in")


