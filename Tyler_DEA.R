library(dplyr)
library(ggpubr)
library(ggrepel)
library(tidyverse)
library(grid) #for weird volcano plot labels 

# source("human_islet/limmaFit.R")
# source("human_islet/limmaDEA.R")
source("limmaFit.R")
source("limmaDEA.R")

# setwd("C:/Users/sage980/OneDrive - PNNL/JMF_ANI_Collabs/human_islet")
#switch with file specification or "switch directories" in PNNL ADDINS

# original object, called m 
# load("msnset_humanislet_log2center_august2024_forTyler.RData")

#==========================================================================
#calc LogFC 

#mlc 
load("msnset_humanislet_int_log2center_oct2024_forTyler_wmodanno.RData")

m <-  mlc 

m <- m[, with(pData(m), order(Pair, Treated))]

fData(m) <- fData(m) %>%
  dplyr::rename(n_unexpected_modifications = `#unexpected modifications`)

## Keep features that have complete data for at least two pairs ----
pair_design <- model.matrix(~ 0 + m$Pair)

mat_nonmiss <- apply(!is.na(exprs(m)), 2, as.integer)
rownames(mat_nonmiss) <- rownames(m)
n_nonmiss <- mat_nonmiss %*% pair_design

keep <- apply(n_nonmiss == 2L, 1, function(xi) sum(xi) >= 2L)
table(keep)
#numbers changed from when i oringially tried this 
# FALSE  TRUE
#   597   829

m2 <- m[keep, ]


## Normalize data ----
m2 <- MSnSet.utils::normalizeByGlob(m2, method = "once")

boxplot(exprs(m2)) 
# Different samples seem to have different variances


## MDS plot ----

# Color by pairs
idx <- match(m2$Pair, sort(unique(m2$Pair)))

col <- c("black", "red", "green", "blue", "orange", "purple")
col <- col[idx]

labels <- c(m2$Pair)

#this is screwed up now, plotting full dataset name instead of letter 
pdf(file = "humanislet_processed_MSnSet_MDS_plot_2.pdf",
    width = 6, height = 5)
plotMDS(exprs(m2), col = col, labels=labels)
dev.off()

## Differential analysis ----

# Create fitted object. The block argument makes this a paired comparison
fit <- limmaFit(object = m2,
                model.str = "~ 0 + Treated",
                contrasts = c("TreatedTreated - TreatedNottreated"),
                block = "Pair", plot = TRUE,
                trend = TRUE, robust = TRUE)

# DEA results table
res <- limmaDEA(fit = fit) %>%
  arrange(P.Value) %>%
  mutate(contrast = gsub("Treated(?! )", "", contrast, perl = TRUE))

# save(res, file="Tyler_DEA_results.RData")
# 
# write.table(res, file = "Tyler_DEA_results.txt", 
#             quote = FALSE, row.names = FALSE, sep = "\t")

library(dplyr)

# Identify list columns
list_columns <- sapply(res, is.list)

# Convert list columns to character vectors to make them CSV-friendly
res <- res %>%
  mutate(across(where(is.list), ~ sapply(., toString)))

#add data completeness
library(proDA)

x_final <- as.data.frame(m) %>%
  t()

rowwise_na_summary <- apply(x_final, 1, function(x) sum(!is.na(x))) %>%
  as.data.frame() %>%
  rownames_to_column(var = "PF")

rowwise_na_summary$notna <- rowwise_na_summary$.

rowwise_na_summary <- rowwise_na_summary %>%
  mutate(MissingPercent = (notna/ncol(x_final)))

#it's NOTna/total number, so 1.0 means PF is in all samples 
res <- res %>%
  left_join(rowwise_na_summary, by="PF")

# # Write to CSV after conversion
# write.csv(res, "humanislets_DEA_results_2.csv", row.names = FALSE)

# Check shape of p-value histogram
pdf(file = "humanislet_DEA_p-value_histogram_2.pdf",
    height = 5, width = 6)
hist(res$P.Value, breaks = seq(0, 1, 0.05),
     xlab = "P-Value",
     main = NA)
dev.off()

# Count number of significant proteoforms at 10% FDR
table(res$adj.P.Val < 0.1)
# FALSE  TRUE
#   825     4

# Due to the small number of samples, and the appearance of the MDS plot, the
# FDR threshold was relaxed from the usual 5% to 10%.

#==========================================================================
#plotting

library(dplyr)
library(ggpubr)
library(ggrepel)
library(tidyverse)
library(grid) #for weird volcano plot labels 

# save(res, file="res_forvolcano.RData")

load("res_forvolcano.RData")

write.csv(res, file= "logFC_TDislet.csv")

res %>%
  filter(P.Value < 0.05) %>%
  filter(logFC > 0) %>%
  group_by(Gene) %>%
  tally() %>%
  View()

res %>%
  filter(P.Value < 0.05) %>%
  filter(logFC < 0) %>%
  group_by(Gene) %>%
  tally() %>%
  View()

res %>%
  filter(Gene == "INS") %>%
  ggplot(aes(x = logFC, y = -log10(P.Value), size = MissingPercent)) +
  geom_point()+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray", size = 1) 

offset <- 0

res %>%
  filter(Gene == "INS") %>%
  # filter(P.Value < 0.05) %>% 
  filter(logFC > 1 | logFC < -1) %>%
  mutate(Significant = ifelse(P.Value < 0.05, "Yes", "No")) %>%
  filter(MissingPercent > 0.7) %>%
  ggplot(aes(x = firstAA, xend = lastAA, y = -log10(P.Value), yend = -log10(P.Value), color=Significant)) + 
  geom_segment(size=2) +
  # geom_text(aes(label = paste0(Proteoform, featureName), y = -log10(P.Value)+0.01), size = 3, vjust = 0) +  # Adjust y value and size as needed
  # scale_color_gradient2(low = "blue", mid = "black", high = "red", midpoint = 0) + 
  labs(x = "Amino Acid Sequence", y = "-log10(P.value)") +
  theme_minimal()+
  annotate("rect", xmin=1-offset, xmax=24-offset, ymin=-Inf, ymax=Inf, alpha=0.2, fill="gray") +
  annotate("rect", xmin=25-offset, xmax=54-offset, ymin=-Inf, ymax=Inf, alpha=0.2, fill="blue") +
  annotate("rect", xmin=57-offset, xmax=87-offset, ymin=-Inf, ymax=Inf, alpha=0.2, fill="red") +
  annotate("rect", xmin=90-offset, xmax=110-offset, ymin=-Inf, ymax=Inf, alpha=0.2, fill="green") 

res %>%
  filter(Gene == "GCG") %>%
  # filter(P.Value < 0.05) %>% 
  # filter(logFC > 1 | logFC < -1) %>%
  mutate(Significant = ifelse(P.Value < 0.05, "Yes", "No")) %>%
  filter(MissingPercent > 0.7) %>%
  ggplot(aes(x = firstAA, xend = lastAA, y = -log10(P.Value), yend = -log10(P.Value), color=Significant)) + 
  geom_segment(size=2) +
  geom_text(aes(label = paste0(Proteoform, featureName), y = -log10(P.Value)+0.01), size = 3, vjust = 0) +  # Adjust y value and size as needed
  # scale_color_gradient2(low = "blue", mid = "black", high = "red", midpoint = 0) + 
  labs(x = "Amino Acid Sequence", y = "-log10(P.value)") +
  theme_minimal()+
  annotate("rect", xmin=1-offset, xmax=20-offset, ymin=-Inf, ymax=Inf, alpha=0.2, fill="gray") + #sp
  annotate("rect", xmin=21-offset, xmax=50-offset, ymin=-Inf, ymax=Inf, alpha=0.2, fill="blue") + #glicentin-related polypeptide 
  annotate("rect", xmin=53-offset, xmax=81-offset, ymin=-Inf, ymax=Inf, alpha=0.2, fill="red") + #glucagon
  annotate("rect", xmin=53-offset, xmax=89-offset, ymin=-Inf, ymax=Inf, alpha=0.2, fill="red") + #Oxyntomodulin
  annotate("rect", xmin=92-offset, xmax=128-offset, ymin=-Inf, ymax=Inf, alpha=0.2, fill="green") + #GLP1
  annotate("rect", xmin=146-offset, xmax=178-offset, ymin=-Inf, ymax=Inf, alpha=0.2, fill="purple") #GLP2

#======================================================================
#for manuscript 
volcanoplot <- res %>%
  mutate(significant = case_when(
    adj.P.Val > 0.05 & (logFC < -1 | logFC > 1) ~ "log2FC",
    adj.P.Val < 0.05 & logFC <= -1 ~ "Downregulated",
    adj.P.Val < 0.05 & logFC >= 1 ~ "Significant after adjustment",
    adj.P.Val > 0.05 & logFC >= -1 & logFC <= 1 ~ "NS",
    TRUE ~ "Other" # Adding a fallback category for debugging or handling unexpected cases
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
    plot.margin = unit(c(1, 1, 3.5, 1), "lines") # Add margin space at the bottom
  ) +
  scale_size_continuous(
    range = c(10, 2), # Adjust the size range for points
    labels = c("0%", "20%", "40%", "60%", "80%", "100%"), # Ensure labels cover the percentage range
    breaks = c(0, 20, 40, 60, 80, 100) # Breaks corresponding to percentages
  ) +
  guides(
    size = guide_legend(title = "Missing Data (%)"),
    color = guide_legend(override.aes = list(size = 6), title = "") # Increase the size of points in color legend
  ) +
  labs(x = expression(Log[2] * " fold change"), y = expression(-Log[10] * "P")) +
  scale_color_manual(values = c("NS" = "black", "log2FC" = "darkgreen", "Significant after adjustment" = "red", "Downregulated" = "red"),
                     labels = c(
                       "NS" = "Not significant",
                       "log2FC" = expression(Log[2] * " fold change"),
                       "Upregulated" = "Upregulated",
                       "Downregulated" = "Downregulated"
                     )) +
  annotate("text", x = 4.5, y = 2.5, label = "Up in Treatment", color = "black", size = 5, hjust = 1) +
  annotate("text", x = -4.5, y = 2.5, label = "Down in Treatment", color = "black", size = 5, hjust = 0)+
  xlim(-5,5)+
  geom_text_repel(aes(label = ifelse(-log10(adj.P.Val) > -log10(0.05), Gene, NA)), size = 5, nudge_x = 0.2, nudge_y = -0.05, hjust = 0, vjust = 1, show.legend  = F)
volcanoplot
# geom_text_repel(aes(label = ifelse(-log10(adj.P.Val) > -log10(0.05), Gene, NA)), size = 5, nudge_x = -0.1, nudge_y = -0.1, hjust = 1, vjust = 1)
# geom_text_repel(aes(label = ifelse(-log10(adj.P.Val) > -log10(0.05), MissingPercent, "")), size = 5) #add repel to points if needed 
# geom_text_repel(aes(label = ifelse(-log10(adj.P.Val) > -log10(0.05), featureName, "")), size = 5) #add repel to points if needed 

# ggsave(plot = volcanoplot, filename= "paired_volcano.png", scale=2.2, width = 5,
#        height = 3.5,
#        units = c("in"))

load("x_recluster_oct2024_int.RData")

#this is either methylation 
#could also be F-> Cys with carbidomethyl (Phe->CamCys)
x_recluster %>%
  mutate(PF = paste(Gene, pcGroup, sep="_")) %>%
  filter(PF == "GCG_106") %>%
  # filter(PF == "HMGN2_8") %>%
  View()

#looking at all forms of GLP1 
x_recluster %>%
  mutate(PF = paste(Gene, pcGroup, sep="_")) %>%
  filter(Gene == "GCG") %>%
  filter(firstAA > 88 & firstAA <132) %>%
  filter(lastAA > 88 & lastAA <132) %>%
  group_by(PF) %>%
  slice_min(`E-value`) %>%
  View()

#boxplot of methylated GLP1, m2 is globalized to norm in this script, data is log2centered before hand 
exprs(m2) %>%
  as.data.frame() %>%
  rownames_to_column(var="PF") %>%
  filter(PF == "GCG_106") %>%
  # filter(PF == "HMGN2_8") %>%
  pivot_longer(-PF, names_to = "sample_name",
               values_to = "value") %>%
  separate(sample_name, into=c("TD", "islet", "Letter","Treated", "Time", "Biorep"), remove = FALSE) %>%
  mutate(Pair = case_when(Letter == "A" ~ "Pair1",
                          Letter == "H" ~ "Pair1",
                          Letter == "B" ~ "Pair2",
                          Letter == "D" ~ "Pair2",
                          Letter == "C" ~ "Pair3",
                          Letter == "N" ~ "Pair3",
                          Letter == "E" ~ "Pair4",
                          Letter == "J"~ "Pair4",
                          Letter == "F" ~ "Pair5",
                          Letter == "K" ~ "Pair5",
                          Letter == "G" ~ "Pair6",
                          Letter == "M" ~ "Pair6",
                          Letter == "I" ~ "Pair7",
                          Letter == "L" ~ "Pair7")) %>%
  filter(!is.na(value)) %>%
  # View()
  ggplot()+
  geom_boxplot(aes(x=Treated, y=value, color=Treated))+
  geom_point(aes(x=Treated, y=value,shape=Pair, size=5))+
  theme_minimal(base_size=24)+
  labs(x="Treatment", y="log2center, \nnorm_global", main="methyl-GLP1", shape = "Patient")
  
#there are intensity values for HMGN2_8 in files without ms2 scans, so match_features has been doing work here, do we need to sanity check this? 
x_recluster %>%
  mutate(PF = paste(Gene, pcGroup, sep="_")) %>%
  filter(PF == "HMGN2_8") %>%
  View()

#there 8 distinct PFs of HMNG2
x_recluster %>%
  mutate(PF = paste(Gene, pcGroup, sep="_")) %>%
  filter(Gene == "HMGN2") %>%
  group_by(PF) %>%
  slice_min(`E-value`) %>%
  View()

# HMGN only 
plot <- res %>%
  mutate(MissingPercent = 100 * (1 - MissingPercent)) %>% # Reverse the scale and convert to percentage
  filter(Gene == "HMGN2") %>%
  ggplot(aes(x=logFC, y=-log10(adj.P.Val), size = MissingPercent))+
  geom_hline(yintercept = -log10(0.05))+
  geom_point(alpha =0.5)+
  theme_minimal(base_size = 18)+
  ggtitle("HMGN2")+
  geom_text_repel(aes(label = paste(firstAA, lastAA, MissingPercent, sep = "-")), size = 5)
plot
  # geom_text_repel(aes(label = ifelse(-log10(adj.P.Val) > -log10(0.2), featureName, "")), size = 5)

# res %>%
#   filter(Gene == "IAPP") %>%
#   ggplot(aes(x=logFC, y=-log10(adj.P.Val)))+
#   geom_hline(yintercept = -log10(0.05))+
#   geom_point(color="blue", size=4, alpha =0.5)+
#   theme_minimal(base_size = 18)+
#   ggtitle("IAPP")+
#   geom_text_repel(aes(label = ifelse(-log10(adj.P.Val) > -log10(0.2), featureName, "")), size = 5) 
# 
# res %>%
#   filter(Gene == "GCG") %>%
#   ggplot(aes(x=logFC, y=-log10(adj.P.Val)))+
#   geom_hline(yintercept = -log10(0.05))+
#   geom_point(color="red", size=4, alpha =0.5)+
#   theme_minimal(base_size = 18)+
#   ggtitle("GCG")+
#   geom_text_repel(aes(label = ifelse(-log10(adj.P.Val) > -log10(0.2), featureName, "")), size = 5) 
# 
# res %>%
#   filter(Gene == "CHGA") %>%
#   ggplot(aes(x=logFC, y=-log10(adj.P.Val)))+
#   geom_hline(yintercept = -log10(0.05))+
#   geom_point(color="green", size=4, alpha =0.5)+
#   theme_minimal(base_size = 18)+
#   ggtitle("CHGA")+
#   geom_text_repel(aes(label = ifelse(-log10(adj.P.Val) > -log10(0.2), featureName, "")), size = 5) 





