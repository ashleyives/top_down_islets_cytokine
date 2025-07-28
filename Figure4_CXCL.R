# Authors: Ashley Ives
# Purpose: Takes output of script "0a" and generates Figure 4: Heatmap and table summary of observed CXCL proteoforms. 
# Performs hypergeometric probability distribution test to determine which unique observations are statistically significant.

library(tidyverse)
library(Seurat)
library(viridis)
library(gt)

# 1. Load data generated in "0a_loading_intensity_data.R" and convert to usable format.  

load("msnset_humanislet_int_notnormalized.RData")

exprs(m) <- log2(exprs(m))

norm_coeffs <- apply(exprs(m), 2, median, na.rm = T) #calculates the median of each column in the expression matrix while ignoring any NA values.
exprs(m) <- sweep(exprs(m), 2, norm_coeffs, FUN = "-") #applies normalization to msnset via sweep, because it's log transformed us Fun = "_"

exprs(m) <- 2^(exprs(m))

x <- exprs(m) %>%
  as.data.frame() %>%
  rownames_to_column(var = "PF") %>%
  pivot_longer(-PF,
               values_to = "Intensity",
               names_to = "SampleID") %>%
  mutate(Treatment = case_when(grepl("Treated", SampleID) ~ "Treated",
                               grepl("Nottreated", SampleID) ~"Untreated")) %>%
  mutate(Pair = case_when(grepl("_A_", SampleID) ~ "Pair1",
                          grepl("_H_", SampleID) ~ "Pair1",
                          grepl("_B_", SampleID) ~ "Pair2",
                          grepl("_D_", SampleID) ~ "Pair2",
                          grepl("_C_", SampleID) ~ "Pair3",
                          grepl("_N_", SampleID) ~ "Pair3",
                          grepl("_E_", SampleID) ~ "Pair4",
                          grepl("_J_", SampleID) ~ "Pair4",
                          grepl("_F_", SampleID) ~ "Pair5",
                          grepl("_K_", SampleID) ~ "Pair5",
                          grepl("_G_", SampleID) ~ "Pair6",
                          grepl("_M_", SampleID) ~ "Pair6",
                          grepl("_I_", SampleID) ~ "Pair7",
                          grepl("_L_", SampleID) ~ "Pair7"))

SampleID_meta <- x %>%
  distinct(SampleID, Treatment, Pair)

# 2. Perform hypergeometric probability distribution----------------------------

uniquePFs <- x %>%
  filter(!is.na(Intensity)) %>%
  group_by(PF) %>%
  add_count(name = "TotalN") %>%
  group_by(Treatment, PF) %>%
  add_count(name = "GroupN") %>%
  ungroup() %>%
  group_by(Treatment) %>%
  mutate(GroupSize = max(GroupN)) %>%
  filter(GroupSize == GroupN) %>%
  mutate(Check1 = TotalN - GroupSize) %>%
  ungroup() %>%
  mutate(Prob0 = (factorial(nrow(SampleID_meta)-GroupSize))/(factorial(Check1)*factorial((nrow(SampleID_meta)-GroupSize)-Check1))) %>%
  mutate(Prob1 = (factorial(nrow(SampleID_meta)))/(factorial(TotalN)*factorial(nrow(SampleID_meta)-TotalN))) %>%
  mutate(Probability = Prob0/Prob1) %>%
  group_by(Treatment, PF, Probability) %>%
  summarize(Avg = mean(Intensity))  %>%
  ungroup() %>%
  filter(Probability < 0.01) 

meta <- fData(m) %>%
  dplyr::select(firstAA, lastAA, PF, Proteoform, Gene, UniProtAcc)%>%
  mutate(isSig = ifelse(PF %in% uniquePFs$PF, TRUE, FALSE))

pair_mapping <- c('Pair1' = 'A', 'Pair2' = 'B', 'Pair3' = 'C', 'Pair5' = 'D', 'Pair6' = 'E', 'Pair7' = 'F')

# Define the specific order for PF values
specific_order <- c("CXCL11_2", "CXCL11_1", "CXCL10_8", "CXCL10_7", "CXCL10_6", 
                    "CXCL10_4", "CXCL10_3", "CXCL10_2", "CXCL9_1",  "CXCL1_1")

# 3. Generate Figure 4. Heatmap of CXCL proteoforms.---------------------------- 

plot <- x %>%
  filter(grepl("CXCL", PF)) %>%
  mutate(Patient = recode(Pair, !!!pair_mapping)) %>%
  mutate(Sample = paste(Treatment, Patient, sep  = "_")) %>%
  left_join(meta) %>%
  mutate(PF = factor(PF, levels = specific_order)) %>%  
  ggplot()+
  # aes(x = Sample, y = fct_reorder(PF, Intensity), fill = log2(Intensity))+
  aes(x = Sample, y = (PF), fill = log2(Intensity))+
  geom_tile()+
  geom_text(aes(label = ifelse(isSig & Treatment == "Treated", "*", "")), color = "white", size = 10) +
  scale_fill_viridis(na.value = "gray") + 
  # scale_fill_distiller(palette = "RdPu",
  #                      na.value = "black")+
  theme_minimal(base_size = 16)+
  theme(axis.text.y=element_text(color = 'black', size = 12),
        axis.text.x=element_text(angle = 90, vjust = 0.5, hjust= 0.5, color='black',size = 12))+
  xlab("Treatment_Patient")+
  ylab("Proteoform Identifier")
plot

table <- meta %>%
  filter(grepl("CXCL", PF)) %>%
  dplyr::select(PF, Gene, UniProtAcc, firstAA, lastAA, isSig) %>%
  mutate(Significance = ifelse(isSig, "*", "")) %>%
  dplyr::select(-isSig) %>%  # Remove the isSig column as it's not needed in the table
  mutate(PF = factor(PF, levels = specific_order)) %>%  
  arrange(desc(PF))  %>% # Arrange the data based on the ordered PF factor levels
  gt() %>%
  # tab_header(
  #   title = "Loadings and Custom Labels"
  # ) %>%
  cols_label(
    PF = md("**Identifier**"),
    Gene = md("**Gene**"),
    UniProtAcc = md("**Accession**"),
    firstAA = md("**firstAA**"),
    lastAA = md("**lastAA**"),
    Significance = md("**Significant**") # Add label for significance
  )

library(patchwork)

combplot <- (plot/table)+
  plot_layout(heights = unit(c(2, 1.5), "null")) + # Adjust the relative heights
  plot_annotation(tag_levels = 'A')  &
  theme(plot.tag = element_text(size = 24))  # Adjust the text size
combplot

ggsave(plot = combplot, filename= "FigureS_CXCL.png", scale=1.5,
       width = 170,
       height = 170,
       dpi = 800,
       units = c("mm"))
