# Authors: Ashley Ives 
# Purpose: Takes output of script "0a" and generates Figure S2: Histogram of relative standard deviations based on label-free quantification data. Facets data into control and cytokine treated samples.

library(dplyr)
library(tibble)
library(MSnSet.utils)
library(proteomicsCV)
library(tidyverse)

# 1. Load raw intensity data, convert to log2 scale and perform median normalization

nm = load("msnset_humanislet_int_notnormalized.RData")

exprs(m) <- log2(exprs(m))

norm_coeffs <- apply(exprs(m), 2, median, na.rm = T) #calculates the median of each column in the expression matrix while ignoring any NA values.
exprs(m) <- sweep(exprs(m), 2, norm_coeffs, FUN = "-") #applies normalization to msnset via sweep, because it's log transformed us Fun = "_"

exprs(m) <- 2^(exprs(m))

data <- as.data.frame((exprs(m)))

# Subset the columns based on their names containing "Treated" or "Untreated"
treated <- data[, grepl("Treated", colnames(data))]
untreated <- data[, grepl("Nottreated", colnames(data))]

# Plotting of RSDs, split into cytokine treated and control samples ------------

treatedplot <- treated %>%
  rowwise() %>%
  mutate(FIrsd = sd(c_across(everything()) * 100, na.rm = TRUE) / mean(c_across(everything()), na.rm = TRUE)) %>%
  ggplot()+
  aes(x = FIrsd)+
  geom_histogram(binwidth = 2, color = "black", fill = "white", position = 'dodge') +
  geom_vline(aes(xintercept=median(FIrsd, na.rm = TRUE)), color = "red", linetype = "dashed", size = 2) +
  # scale_x_continuous(limits = c(0,100), expand = c(0,1))+
  # scale_y_continuous(expand = c(0,1))+
  theme_bw(base_size = 16) +
  theme( panel.border = element_blank(),
         panel.background = element_rect(fill= 'white'),
         axis.text.y=element_text(color = 'black', size = 16),
         axis.text.x=element_text(color='black',size = 16),
         panel.grid.minor = element_blank()) +
  xlab("Relative Standard Deviation (%)")+
  ggtitle("Treated islets, N = 6")+
  xlim(0,200)
treatedplot

untreatedplot <- untreated %>%
  rowwise() %>%
  mutate(FIrsd = sd(c_across(everything()) * 100, na.rm = TRUE) / mean(c_across(everything()), na.rm = TRUE)) %>%
  ggplot()+
  aes(x = FIrsd)+
  geom_histogram(binwidth = 2, color = "black", fill = "white", position = 'dodge') +
  geom_vline(aes(xintercept=median(FIrsd, na.rm = TRUE)), color = "red", linetype = "dashed", size = 2) +
  # scale_x_continuous(limits = c(0,100), expand = c(0,1))+
  # scale_y_continuous(expand = c(0,1))+
  theme_bw(base_size = 16) +
  theme( panel.border = element_blank(),
         panel.background = element_rect(fill= 'white'),
         axis.text.y=element_text(color = 'black', size = 16),
         axis.text.x=element_text(color='black',size = 16),
         panel.grid.minor = element_blank()) +
  xlab("Relative Standard Deviation (%)")+
  ggtitle("Untreated islets, N = 6")+
  xlim(0,200)
untreatedplot

library(patchwork)

combined_plot <- (untreatedplot/treatedplot)+
  plot_annotation(tag_levels = list('A'))& 
  theme(plot.margin = unit(c(0.001, 0.001, 0.001, 0.001), "cm"))
combined_plot

ggsave(plot = combined_plot, filename= "FigureS_RSD.png", scale=1,
       width = 170,
       height = 170,
       dpi = 800,
       units = c("mm"))


  