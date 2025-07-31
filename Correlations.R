# Correlation plots for all
# E. Lamont 
# 7/30/25

# http://www.sthda.com/english/wiki/ggcorrplot-visualization-of-a-correlation-matrix-using-ggplot2

source("Import_data.R") # to get All_tpm, which just keeps the samples I am interested in that pass filter


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
###################### PROCESS DATA #######################

# Need to remove the gene column or it won't work
All_tpm2 <- All_tpm %>% column_to_rownames("X")

# Log10 transform
All_tpm2_Log10 <- All_tpm2 %>%
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values

# pairs(All_tpm_numeric)


###########################################################
############### RARITY ALL GRAPHS TOGETHER ################

# https://borisleroy.com/en/2013/06/09/correlation-plots-in-r/
# install.packages("Rarity")
library(Rarity)

# DO NOT RUN THIS!!! FILE IS TOO BIG!!!!
# Pearson
# Make a new folder for the figure first
# pdf("Figures/Correlations/rarity_PearsonLog10_AllSamples_v1.pdf", width = 20, height = 20)
# corPlot(All_tpm2_Log10, method = "pearson",
#         title = "Pearson Correlation Log10(TPM+1)")
# dev.off()









