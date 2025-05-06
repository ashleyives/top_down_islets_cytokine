
library(dplyr)
library(ggpubr)
library(ggrepel)
library(tidyverse)
library(grid) #for weird volcano plot labels 
library(gprofiler2) #gost

colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Data is made in Tyler_DEA.R 
load("res_forvolcano.RData")

res <- res %>%
  filter(notna > 4) #must be in at least 4 of 12 total samples

#what is downregulated 
down <- res %>%
  filter(P.Value < 0.05) %>%
  filter(logFC < 0) 

#what is upregulated 
up <- res %>%
  filter(P.Value < 0.05) %>%
  filter(logFC > 0)

#gost of up and down regulated terms, background should be observed genes, not entire proteome  
#takes hgnc terms, correct background is not human genome but all observed proteins in sample
gostres <- gprofiler2::gost(query = unique(up$Gene),
                            organism = "hsapiens",
                            ordered_query = F, 
                            multi_query = FALSE, significant = FALSE, exclude_iea = F, 
                            measure_underrepresentation = FALSE, evcodes = TRUE, 
                            user_threshold = 0.05, correction_method = "g_SCS", 
                            domain_scope = "annotated", 
                            custom_bg = unique(res$Gene), 
                            numeric_ns = "", sources = c("GO"), as_short_link = FALSE)

gostres$result %>% View()

gostres <- gprofiler2::gost(query = unique(down$Gene),
                            organism = "hsapiens",
                            ordered_query = F, 
                            multi_query = FALSE, significant = FALSE, exclude_iea = F, 
                            measure_underrepresentation = FALSE, evcodes = TRUE, 
                            user_threshold = 0.05, correction_method = "g_SCS", 
                            domain_scope = "annotated", 
                            custom_bg = unique(res$Gene), 
                            numeric_ns = "", sources = c("GO"), as_short_link = FALSE)

gostres$result %>% View()

minpoint <- 2
maxpoint <- 6

xpositionlab <- 2.4 

#p values for all entries 
panelp <- res %>%
  # mutate(significant = case_when(
  #   P.Value > 0.05 & (logFC < -1 | logFC > 1) ~ "log2FC",
  #   P.Value < 0.05 & logFC <= -1 ~ "Downregulated",
  #   P.Value < 0.05 & logFC >= 1 ~ "Significant after adjustment",
  #   P.Value > 0.05 & logFC >= -1 & logFC <= 1 ~ "NS",
  #   TRUE ~ "Other" # Adding a fallback category for debugging or handling unexpected cases
  # )) %>%
  mutate(significant = case_when(
    P.Value > 0.05 & (logFC < -1 | logFC > 1) ~ "black",
    P.Value < 0.05 & logFC <= -1 ~ "red3",
    P.Value < 0.05 & logFC >= 1 ~ "darkgreen",
    P.Value > 0.05 & logFC >= -1 & logFC <= 1 ~ "black",
    TRUE ~ "black" # Adding a fallback category for debugging or handling unexpected cases
  )) %>%
  mutate(MissingPercent = 100 * (1 - MissingPercent)) %>% # Reverse the scale and convert to percentage
  # ggplot(aes(x = logFC, y = -log10(adj.P.Val), color = significant, size = MissingPercent)) +
  ggplot(aes(x = logFC, y = -log10(P.Value), color = significant, size = MissingPercent)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray", size = 1) +
  geom_vline(xintercept = -1, linetype = "dashed", color = "gray", size = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray", size = 1) +
  geom_point(alpha = 0.5) +
  theme_minimal(base_size = 18) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.key.size = unit(0.5, "cm"),  # Shrinking the size of the legend keys
    legend.text = element_text(size = 10),  # Shrinking the size of the legend text
    legend.title = element_text(size = 12),   # Shrink the size of the legend title
    plot.margin = unit(c(1, 1, 3.5, 1), "lines") # Add margin space at the bottom
  ) +
  scale_size_continuous(
    range = c(maxpoint, minpoint), # Adjust the size range for points
    labels = c("0%", "20%", "40%", "60%", "80%", "100%"), # Ensure labels cover the percentage range
    breaks = c(0, 20, 40, 60, 80, 100) # Breaks corresponding to percentages
  ) +
  guides(
    size = guide_legend(title = "Missing Data (%)"),
    color = guide_legend(override.aes = list(size = 6), title = "") # Increase the size of points in color legend
  ) +
  labs(x = expression(Log[2] * " fold change"), y = expression(-Log[10] * "P")) +
  scale_color_manual(values = c("black" = "black", "red3" = "red3", "darkgreen" = "darkgreen"),
                     labels = c(
                       "black" = "Not significant",
                       "red3" = "Downregulated",
                       "darkgreen" = "Upregulated"
                     )) +
  annotate("text", x = 5, y = 0.2, label = "Up in Treatment", color = "black", size = 5, hjust = 1) +
  annotate("text", x = -5, y = 0.2, label = "Down in Treatment", color = "black", size = 5, hjust = 0)+
  xlim(-5,5)+
  geom_text_repel(aes(label = ifelse(-log10(adj.P.Val) > -log10(0.05), Gene, NA)), size = 5, nudge_x = 0.2, nudge_y = -0.05, hjust = 0, vjust = 1, show.legend  = F)+
  ggtitle("All Proteforms")
panelp

#only gcg with categories based on fragment range 
panelgcg <- res %>%
  filter(Gene == "GCG") %>%
  mutate(segment = case_when(
    between(firstAA, 18, 52) & between(lastAA, 18, 52) ~ "Glicentin-related \n polypeptide",
    between(firstAA, 53, 91) & between(lastAA, 53, 91) ~ "Oxyntomodulin",
    between(firstAA, 92, 178) & between(lastAA, 92, 178) ~ "Major proglucagon \n fragment or GLP-1",
    TRUE ~ "Other"
  )) %>%
  # filter(segment == "Other")
  mutate(MissingPercent = 100 * (1 - MissingPercent)) %>% # Reverse the scale and convert to percentage
  ggplot(aes(x = logFC, y = -log10(P.Value), color = segment, size = MissingPercent)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray", size = 1) +
  geom_vline(xintercept = -1, linetype = "dashed", color = "gray", size = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray", size = 1) +
  geom_point(alpha = 0.5) +
  theme_minimal(base_size = 18) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.key.size = unit(0.5, "cm"),  # Shrinking the size of the legend keys
    legend.text = element_text(size = 10),  # Shrinking the size of the legend text
    legend.title = element_text(size = 12),   # Shrink the size of the legend title
    plot.margin = unit(c(1, 1, 3.5, 1), "lines") # Add margin space at the bottom
  ) +
  scale_size_continuous(
    range = c(maxpoint, minpoint),  # Adjust the size range for points
    labels = c("0%", "20%", "40%", "60%", "80%", "100%"), # Ensure labels cover the percentage range
    breaks = c(0, 20, 40, 60, 80, 100) # Breaks corresponding to percentages
  ) +
  scale_color_manual(values = colorBlindGrey8) +
  guides(
    size = guide_legend(title = "Missing Data (%)"),
    color = guide_legend(override.aes = list(size = 6), title = "GCG Fragment") # Increase the size of points in color legend
  ) +
  labs(x = expression(Log[2] * " fold change"), y = expression(-Log[10] * "P")) +
  annotate("text", x = xpositionlab, y = 0.2, label = "Up in Treatment", color = "black", size = 5, hjust = 1) +
  annotate("text", x = -xpositionlab, y = 0.2, label = "Down in Treatment", color = "black", size = 5, hjust = 0)+
  xlim(-2.5,2.5)+
  ggtitle("GCG Proteforms")
panelgcg

#how many proteoforms go up or down per region 
res %>%
  filter(Gene == "GCG") %>%
  mutate(segment = case_when(
    between(firstAA, 18, 52) & between(lastAA, 18, 52) ~ "Glicentin-related polypeptide",
    between(firstAA, 53, 91) & between(lastAA, 53, 91) ~ "Oxyntomodulin",
    between(firstAA, 92, 178) & between(lastAA, 92, 178) ~ "Major proglucagon fragment/GLP",
    TRUE ~ "Other"
  )) %>%
  # filter(segment == "GLP") %>%
  # filter(segment == "Oxyntomodulin") %>%
  filter(segment == "Glicentin") %>%
  filter(logFC > 0) %>%
  nrow()
  

#only CHGA with categories based on fragment range 
panelchga <- res %>%
  filter(Gene == "CHGA") %>%
  mutate(segment = case_when(
    between(firstAA, 1, 131) & between(lastAA, 1, 131) ~ "Vasostatin1/2",
    between(firstAA, 132, 225) & between(lastAA, 132, 225) ~ "EA-92",
    # between(firstAA, 228, 275) & between(lastAA, 228, 275) ~ "228-275",
    between(firstAA, 358, 390) & between(lastAA, 358, 390) ~ "LF-19/Catestatin",
    between(firstAA, 413, 457) & between(lastAA, 413, 457) ~ "GR-44/ER-37",
    between(firstAA, 316, 339) & between(lastAA, 316, 339) ~ "SS-18",
    TRUE ~ "Other"
  )) %>%
  mutate(MissingPercent = 100 * (1 - MissingPercent)) %>% # Reverse the scale and convert to percentage
  ggplot(aes(x = logFC, y = -log10(P.Value), color = segment, size = MissingPercent)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray", size = 1) +
  geom_vline(xintercept = -1, linetype = "dashed", color = "gray", size = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray", size = 1) +
  geom_point(alpha = 0.5) +
  theme_minimal(base_size = 18) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.key.size = unit(0.5, "cm"),  # Shrinking the size of the legend keys
    legend.text = element_text(size = 10),  # Shrinking the size of the legend text
    legend.title = element_text(size = 12),   # Shrink the size of the legend title
    plot.margin = unit(c(1, 1, 3.5, 1), "lines") # Add margin space at the bottom
  ) +
  scale_size_continuous(
    range = c(maxpoint, minpoint),  # Adjust the size range for points
    labels = c("0%", "20%", "40%", "60%", "80%", "100%"), # Ensure labels cover the percentage range
    breaks = c(0, 20, 40, 60, 80, 100) # Breaks corresponding to percentages
  ) +
  scale_color_manual(values = colorBlindGrey8) +
  guides(
    size = guide_legend(title = "Missing Data (%)"),
    color = guide_legend(override.aes = list(size = 6), title = "CHGA Fragment") # Increase the size of points in color legend
  ) +
  labs(x = expression(Log[2] * " fold change"), y = expression(-Log[10] * "P")) +
  annotate("text", x = xpositionlab, y = 0.2, label = "Up in Treatment", color = "black", size = 5, hjust = 1) +
  annotate("text", x = -xpositionlab, y = 0.2, label = "Down in Treatment", color = "black", size = 5, hjust = 0)+
  xlim(-2.5,2.5)+
  ggtitle("CHGA Proteforms")
panelchga

res %>%
  filter(Gene == "CHGA") %>%
  mutate(segment = case_when(
    between(firstAA, 1, 131) & between(lastAA, 1, 131) ~ "Vasostatin1/2",
    between(firstAA, 132, 225) & between(lastAA, 132, 225) ~ "EA-92",
    between(firstAA, 228, 275) & between(lastAA, 228, 275) ~ "228-275",
    between(firstAA, 358, 390) & between(lastAA, 358, 390) ~ "LF-19/Catestatin",
    between(firstAA, 413, 457) & between(lastAA, 413, 457) ~ "GR-44/ER-37",
    between(firstAA, 316, 339) & between(lastAA, 316, 339) ~ "SS-18",
    TRUE ~ "Other"
  )) %>%
  filter(segment == "Vasostatin1/2") %>%
  # filter(segment == "SS-18") %>%
  # filter(segment == "LF-19/Catestatin") %>%
  filter(logFC > 0) %>%
  nrow()


#only ins 
panelins <- res %>%
  filter(Gene == "INS") %>%
  mutate(segment = case_when(
    between(firstAA, 25, 56) & between(lastAA, 25, 56) ~ "B Chain",
    between(firstAA, 57, 89) & between(lastAA, 57, 89) ~ "C Peptide",
    between(firstAA, 90, 110) & between(lastAA, 90, 110) ~ "A Chain",
    TRUE ~ "Other"
  )) %>%
  mutate(MissingPercent = 100 * (1 - MissingPercent)) %>% # Reverse the scale and convert to percentage
  ggplot(aes(x = logFC, y = -log10(P.Value), color = segment, size = MissingPercent)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray", size = 1) +
  geom_vline(xintercept = -1, linetype = "dashed", color = "gray", size = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray", size = 1) +
  geom_point(alpha = 0.5) +
  theme_minimal(base_size = 18) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.key.size = unit(0.5, "cm"),  # Shrinking the size of the legend keys
    legend.text = element_text(size = 10),  # Shrinking the size of the legend text
    legend.title = element_text(size = 12),   # Shrink the size of the legend title
    plot.margin = unit(c(1, 1, 3.5, 1), "lines") # Add margin space at the bottom
  ) +
  scale_size_continuous(
    range = c(maxpoint, minpoint),  # Adjust the size range for points
    labels = c("0%", "20%", "40%", "60%", "80%", "100%"), # Ensure labels cover the percentage range
    breaks = c(0, 20, 40, 60, 80, 100) # Breaks corresponding to percentages
  ) +
  scale_color_manual(values = colorBlindGrey8) +
  guides(
    color = guide_legend(order=1, override.aes = list(size = 6), title = "INS Fragment"),
    size = guide_legend(order=2,title = "Missing Data (%)")
    # Increase the size of points in color legend
  ) +
  labs(x = expression(Log[2] * " fold change"), y = expression(-Log[10] * "P")) +
  annotate("text", x = xpositionlab, y = 0.2, label = "Up in Treatment", color = "black", size = 5, hjust = 1) +
  annotate("text", x = -xpositionlab, y = 0.2, label = "Down in Treatment", color = "black", size = 5, hjust = 0)+
  xlim(-2.5,2.5)+
  ggtitle("INS Proteforms")
panelins

library(patchwork)

combined_plot <- (panelp + panelins) / (panelgcg + panelchga)+
  plot_annotation(tag_levels = list('A'))& 
  theme(plot.margin = unit(c(0.001, 0.001, 0.001, 0.001), "cm"))
combined_plot

ggsave(plot = combined_plot, filename= "Figure3_volcanoplots.png", 
       scale=2.2, 
       width = 170,
       height = 150,
       dpi = 800,
       units = c("mm"))



#also check SST, no clear pattern here
panelsst <- res %>%
  filter(Gene == "SST") %>%
  mutate(segment = case_when(
    between(firstAA, 89, 116) & between(lastAA, 89, 116) ~ "Somatostatin",
    between(firstAA, 1, 88) & between(lastAA, 1, 88) ~ "Pro",
    TRUE ~ "Other"
  )) %>%
  mutate(MissingPercent = 100 * (1 - MissingPercent)) %>% # Reverse the scale and convert to percentage
  ggplot(aes(x = logFC, y = -log10(P.Value), color = segment, size = MissingPercent)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray", size = 1) +
  geom_vline(xintercept = -1, linetype = "dashed", color = "gray", size = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray", size = 1) +
  geom_point(alpha = 0.5) +
  theme_minimal(base_size = 18) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.key.size = unit(0.5, "cm"),  # Shrinking the size of the legend keys
    legend.text = element_text(size = 10),  # Shrinking the size of the legend text
    legend.title = element_text(size = 12),   # Shrink the size of the legend title
    plot.margin = unit(c(1, 1, 3.5, 1), "lines") # Add margin space at the bottom
  ) +
  scale_size_continuous(
    range = c(maxpoint, minpoint),  # Adjust the size range for points
    labels = c("0%", "20%", "40%", "60%", "80%", "100%"), # Ensure labels cover the percentage range
    breaks = c(0, 20, 40, 60, 80, 100) # Breaks corresponding to percentages
  ) +
  scale_color_manual(values = colorBlindGrey8) +
  guides(
    color = guide_legend(order=1, override.aes = list(size = 6), title = "SST Fragment"),
    size = guide_legend(order=2,title = "Missing Data (%)")
    # Increase the size of points in color legend
  ) +
  labs(x = expression(Log[2] * " fold change"), y = expression(-Log[10] * "P")) +
  # annotate("text", x = 4.5, y = 5.5, label = "Up in Treatment", color = "black", size = 5, hjust = 1) +
  # annotate("text", x = -4.5, y = 5.5, label = "Down in Treatment", color = "black", size = 5, hjust = 0)+
  xlim(-2.5,2.5)+
  ggtitle("SST Proteforms")
panelsst

