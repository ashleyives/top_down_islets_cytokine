# Authors: Ashley Ives
# Purpose: Takes output of script "1_DEA" and generates Figure 3: Volcano plots of all and curated proteoforms.

library(dplyr)
library(ggpubr)
library(ggrepel)
library(tidyverse)
library(grid) #for weird volcano plot labels 

# 0. Set plot defaults----------------------------------------------------------

colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

minpoint <- 2
maxpoint <- 6

xpositionlab <- 2.4 

# 1. Load data generated in "1_DEA.R" and convert to usable format.------------- 

load("res_forvolcano.RData")

res <- res %>%
  filter(notna > 4) #must be in at least 4 of 12 total samples

export <- res %>%
  dplyr::select(-., -mods, -t, -B, -AveExpr, -df.total, -notna, -MissingPercent, -contrast, -pcGroup, -protLength, -PF, -proteoform_id)

write_xlsx(export, "supplementary_table_logfc.xlsx") #for supplementary 

# 2. Generate Figure 3A. Volcano plot of all quantifiable proteoforms.---------- 

panelp <- res %>%
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

# 3. Generate Figure 3C. Volcano plot of all quantifiable glucagon proteoforms.
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

# 4. Generate Figure 3D. Volcano plot of all quantifiable chromograninA proteoforms.
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

# 5. Generate Figure 3B. Volcano plot of all quantifiable insulin proteoforms.

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

# 6. Combine all plots and save-------------------------------------------------

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

