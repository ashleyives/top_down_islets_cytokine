library(tidyverse)
library(TopPICR) 
library(MSnbase)
library(MSnSet.utils)
library(PNNL.DMS.utils)
library(ggpubr)
library(ggplot2)
library(ggnewscale)
library(scales) #colorblind palette
library(ggthemes)

textsize <- 18 
linesize <- 1
linealpha <- 1

# Assign colorblind-friendly palette to a variable
colorblind_palette <- colorblind_pal()(8)

#SST 
###########################################
###########################################
###########################################
###########################################
###########################################
###########################################


load("msnset_humanislet_sc_JMFclustering_wmodanno.RData")

x <- exprs(m) %>%
  as.data.frame()

meta <- m@featureData@data

meta_short <- meta %>%
  mutate(proteoform_id = PF) %>%
  dplyr::select(Gene, proteoform_id, UniProtAcc, mass)

x_long <- x %>%
  rownames_to_column(var = "proteoform_id") %>%
  remove_rownames() %>%
  pivot_longer(-proteoform_id, names_to = "SubjectID",
               values_to = "sc") %>%
  inner_join(.,meta_short)

plot <- fData(m) %>%
  dplyr::mutate(cleanSeq = gsub("\\[.+?\\]",
                                "", Proteoform), cleanSeq = gsub("\\(|\\)", "", cleanSeq),
                cleanSeq = sub("^[A-Z]?\\.(.*)\\.[A-Z]?", "\\1", cleanSeq)) %>%
  filter(!str_detect(Proteoform, "-57.02")) %>% #remove a missed alkylation, which is about same mass as Arg->Val
  filter(!str_detect(Proteoform, "150.0")) %>% #remove DTT adduct 
  # filter(!str_detect(mods, "NH3")) %>% 
  filter(Gene == "SST") %>%
  arrange(desc(count)) %>%  # Sorts the data in descending order based on spectral_counts
  slice_head(n = 30)   # Selects the top 20 rows

# plot %>%
#   select(-mods) %>%
#   mutate(feature_name = PF) %>%
#   select(count, firstAA, lastAA, Proteoform, feature_name, PF, Gene, pcGroup,cleanSeq) %>%
#   write.csv(file="unadjusted_forSSTplot.csv")

#find scans with lowest E value for each pfr filtered in plot above 
load("x_recluster_oct2024_sc.RData")

plot2 <- read.csv("adjusted_forSSTplot.csv")

#set offset for SP, can change back to 0 if you want to rescale plots 
# offset <- 24
offset <- 0

labels <- plot2 %>%
  distinct(feature_name ,.keep_all = TRUE) %>%
  mutate(is_modified = if_else(mod_str == "", FALSE, TRUE)) %>% 
  mutate(firstAA = firstAA-offset) %>%
  mutate(lastAA =lastAA-offset) %>%
  mutate(feature_name2 = paste0("SST","(",firstAA,"-", lastAA,")", if_else(is_modified  == TRUE,"*",""))) %>%
  mutate(feature_name3 = fct_rev(feature_name)) %>%
  select(feature_name, feature_name2) %>%
  arrange((feature_name))

#need to make a new df so you can reference it with "$", be careful plotting labels 
labelsinplot <- plot2 %>%
  arrange(desc(lastAA), desc(firstAA)) %>%
  mutate(feature_name = ordered(feature_name, levels= unique(feature_name))) %>%
  mutate(count_variable = paste("Abundance of SST")) %>%
  left_join(labels)

p0 <- plot2 %>%
  arrange(desc(lastAA), desc(firstAA)) %>%
  mutate(feature_name = ordered(feature_name, levels= unique(feature_name))) %>%
  mutate(count_variable = paste("Abundance of SST")) %>%
  ggplot() +
  aes(y=feature_name, x=count) +
  #before position dodge it was stacking counts multiple times because of different tests, so values looked way higher 
  geom_bar(stat="identity", color = NA, fill = "black", position = "dodge") +
  scale_x_reverse(labels = scales::scientific, breaks = c(0, 1e2, 1e3, 2e3)) +
  scale_fill_gradientn(#colors = MSnSet.utils::hot2.colors(50),
    colors = terrain.colors(50),
    name = expression(paste("log"[10],"(count)"))) +
  # guides(size = guide_legend(title = expression(paste("-log"[10],"(adj p-value)")))) +
  facet_grid(~ count_variable, scales = "free_x", drop = T) + #seems to just be adding header to top, also allows for alignment at end? 
  theme_bw(base_size=16) +
  theme(axis.text.y.right = element_text(hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "left",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_blank(), #removes feature labels
        axis.ticks.y = element_blank(), #removes feature labels
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        # panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background =element_rect(fill="white"),
        strip.text = element_text(face="bold"))+
  labs(y="Proteoform Identifier")+
  scale_y_discrete(name = NULL, expand = expansion(0.02),
                   labels= labelsinplot$feature_name2)
#labels= ABlabels$feature_name2)
p0

temp <-  plot2 %>%
  arrange(desc(lastAA), desc(firstAA)) %>%
  mutate(feature_name = ordered(feature_name, levels= unique(feature_name))) %>%
  mutate(is_modified = mod_str != "") %>%
  mutate(mod_str2 = case_when(mod_str == "" ~ "None", TRUE ~ mod_str)) %>%
  mutate(modification = ordered(.$mod_str2, levels = unique(.$mod_str2))) %>%
  mutate(feature_name = ordered(feature_name, levels= unique(feature_name))) %>% # ordering proteoforms according to N- and C- termini
  mutate(ymin = as.numeric(feature_name) - 0.2,
         ymax = as.numeric(feature_name) + 0.2) 

labelsinplot2 <- temp %>%
  arrange(desc(lastAA), desc(firstAA)) %>%
  mutate(feature_name = ordered(feature_name, levels= unique(feature_name))) %>%
  left_join(labels)

#for abeta   
# temp <- temp %>%
#   mutate(refAA = case_when(terminus == "firstAA" ~ 1, TRUE ~ 40))  

mod_cols <- c("grey60",'#e41a1c','#377eb8','#4daf4a','#ff7f00','#ffff33','#a65628','#f781bf', '#00ACC1', '#673AB7')

library(ggplot2)
library(gridExtra)
library(grid)

rect_data <- data.frame(
  xmin = c(1, 25, 89) - offset,
  xmax = c(24, 88, 116) - offset,
  fill_label = c("SP", "Propeptide", "Somastatin"),
  fill_color  = c("gray", "blue", "green")
)

# Ensure that fill_label in rect_data is a factor to maintain order
rect_data$fill_label <- factor(rect_data$fill_label, levels = rect_data$fill_label)

rect_data$"Fragments of SST" <- rect_data$fill_label 

#need to adjust to SP length so plot starts at -SP
p3 <- temp %>%
  mutate(facet_var = "Proteoforms of SST (P61278)") %>%
  mutate(
    previousAA = str_extract(Proteoform, "^[A-Z](?=\\.)"),  # Extract letter before the first "."
    nextAA = str_extract(Proteoform, "(?<=\\.)([A-Z])$")    # Extract letter after the last "."
  ) %>%
  mutate(firstAA = firstAA - offset) %>%
  mutate(lastAA = lastAA - offset) %>%
  ggplot() +
  # Add the first set of rectangles using rect_data
  # geom_rect(data = rect_data, 
  #           aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = `Fragments of SST`), 
  #           alpha = 0.2, inherit.aes = FALSE, show.legend = TRUE) +
  geom_vline(data = rect_data, aes(xintercept = xmin, color = `Fragments of SST`), size = linesize, alpha = linealpha, linetype = "dashed", inherit.aes = FALSE, show.legend = TRUE) +
  geom_vline(data = rect_data, aes(xintercept = xmax, color = `Fragments of SST`), size = linesize, alpha = linealpha, linetype = "dashed", inherit.aes = FALSE, show.legend = FALSE) +
  scale_colour_colorblind()+
  new_scale_fill() + # Add new scale for the next fill aesthetic
  new_scale_color() + # Add new scale for the next color aesthetic
  geom_rect(aes(ymin = ymin, ymax = ymax,
                xmin = firstAA, xmax = lastAA,
                fill = modification), 
            stat = "identity", color = "black", size = 0.1) +
  facet_grid(cols = vars(facet_var)) +
  # First fill scale for rect_data
  scale_fill_manual(
    name = "Fragments of SST",
    values = setNames(rect_data$fill_color, rect_data$`Fragments of SST`),
    guide = guide_legend(order = 1)
  ) +
  # Second fill scale for main data
  scale_fill_manual(
    name = "Modifications of SST",
    values = mod_cols,
    guide = guide_legend(order = 2)
  ) +
  scale_y_continuous(breaks = seq_along(temp$feature_name),
                     labels = temp$feature_name,
                     expand = expansion(0.015), position = "right") +
  annotation_custom(grob = textGrob(label = "Residue",
                                    gp = gpar(fontsize = textsize),
                                    y = -30, default.units = "pt")) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = textsize) +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"),
        legend.position = "right",
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  scale_y_continuous(breaks = seq_along(temp$feature_name),
                     labels = labelsinplot2$feature_name2,
                     expand = expansion(0.015), position = "right")

p3

library(patchwork)
combined_plot_sst <- p0 + p3 +
  plot_layout(widths = c(0.5, 1),
              guides = "collect")
combined_plot_sst

#CHGA 
###########################################
###########################################
###########################################
###########################################
###########################################
###########################################
###########################################


load("msnset_humanislet_sc_JMFclustering_wmodanno.RData")

x <- exprs(m) %>%
  as.data.frame()

meta <- m@featureData@data

meta_short <- meta %>%
  mutate(proteoform_id = PF) %>%
  dplyr::select(Gene, proteoform_id, UniProtAcc, mass)

x_long <- x %>%
  rownames_to_column(var = "proteoform_id") %>%
  remove_rownames() %>%
  pivot_longer(-proteoform_id, names_to = "SubjectID",
               values_to = "sc") %>%
  inner_join(.,meta_short)

#all rows of n should be 12 if 0s are filled correctly 
x_long %>%
  group_by(proteoform_id) %>%
  tally() %>%
  View()

#summarize by column sc across all 12 samples 
x_long %>%
  group_by(proteoform_id) %>%
  summarize(meansc = mean(sc),
            sdsc = sd(sc)) %>%
  View()

plot <- fData(m) %>%
  dplyr::mutate(cleanSeq = gsub("\\[.+?\\]",
                                "", Proteoform), cleanSeq = gsub("\\(|\\)", "", cleanSeq),
                cleanSeq = sub("^[A-Z]?\\.(.*)\\.[A-Z]?", "\\1", cleanSeq)) %>%
  filter(!str_detect(Proteoform, "-57.02")) %>% #remove a missed alkylation, which is about same mass as Arg->Val
  filter(!str_detect(Proteoform, "150.0")) %>% #remove DTT adduct 
  # filter(!str_detect(mods, "NH3")) %>% 
  filter(Gene == "CHGA") %>%
  arrange(desc(count)) %>%  # Sorts the data in descending order based on spectral_counts
  slice_head(n = 30)   # Selects the top 20 rows

# plot %>%
#   select(-mods) %>%
#   mutate(feature_name = PF) %>%
#   select(count, firstAA, lastAA, Proteoform, feature_name, PF, Gene, pcGroup,cleanSeq) %>%
#   write.csv(file="unadjusted_forchgaplot_2.csv")
# write.csv(file="unadjusted_forchgaplot.csv")
#save this csv and manually edit to make a plot below 

load("x_recluster_oct2024_sc.RData")

#find scans with lowest E value for each pfr filtered in plot above 
#remember that lowest e value scan might not match exact pfr pulled from fData
x_recluster %>%
  mutate(PF = paste(Gene, pcGroup, sep="_")) %>%
  filter(PF %in% plot$PF) %>%
  group_by(PF) %>%
  slice_min(`E-value`) %>%
  View()

x_recluster %>%
  mutate(PF = paste(Gene, pcGroup, sep="_")) %>%
  filter(PF %in% plot$PF) %>%
  filter(str_detect(Proteoform, "Phosph")) %>%
  # group_by(PF) %>%
  # slice_min(`E-value`) %>%
  # filter(`#unexpected modifications` == 1) %>%
  #filter(PF == "GCG_164") %>%
  View()


x_recluster %>%
  mutate(PF = paste(Gene, pcGroup, sep="_")) %>%
  filter(PF %in% plot$PF) %>%
  group_by(PF) %>%
  slice_min(`E-value`) %>%
  View()

#polished plot 

# plot2 <- read.csv("adjusted_forchgaplot.csv")
plot2 <- read.csv("adjusted_forchgaplot_2.csv")

#set offset for SP, can change back to 0 if you want to rescale plots 
# offset <- 18
offset <- 0

labels <- plot2 %>%
  distinct(feature_name ,.keep_all = TRUE) %>%
  mutate(is_modified = if_else(mod_str == "", FALSE, TRUE)) %>% 
  mutate(firstAA = firstAA-offset) %>%
  mutate(lastAA =lastAA-offset) %>%
  mutate(feature_name2 = paste0("Chga","(",firstAA,"-", lastAA,")", if_else(is_modified  == TRUE,"*",""))) %>%
  mutate(feature_name3 = fct_rev(feature_name)) %>%
  select(feature_name, feature_name2) %>%
  arrange((feature_name))

#need to make a new df so you can reference it with "$", be careful plotting labels 
labelsinplot <- plot2 %>%
  arrange(desc(lastAA), desc(firstAA)) %>%
  mutate(feature_name = ordered(feature_name, levels= unique(feature_name))) %>%
  mutate(count_variable = paste("Abundance of CHGA")) %>%
  left_join(labels)

p0 <- plot2 %>%
  arrange(desc(lastAA), desc(firstAA)) %>%
  mutate(feature_name = ordered(feature_name, levels= unique(feature_name))) %>%
  mutate(count_variable = paste("Abundance of CHGA")) %>%
  ggplot() +
  aes(y=feature_name, x=count) +
  #before position dodge it was stacking counts multiple times because of different tests, so values looked way higher 
  geom_bar(stat="identity", color = NA, fill = "black", position = "dodge") +
  scale_x_reverse(labels = scales::scientific, breaks = c(0, 1e1, 3e1, 6e1, 2e2)) +
  scale_fill_gradientn(#colors = MSnSet.utils::hot2.colors(50),
    colors = terrain.colors(50),
    name = expression(paste("log"[10],"(count)"))) +
  # guides(size = guide_legend(title = expression(paste("-log"[10],"(adj p-value)")))) +
  facet_grid(~ count_variable, scales = "free_x", drop = T) + #seems to just be adding header to top, also allows for alignment at end? 
  theme_bw(base_size=16) +
  theme(axis.text.y.right = element_text(hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "left",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_blank(), #removes feature labels
        axis.ticks.y = element_blank(), #removes feature labels
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        # panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background =element_rect(fill="white"),
        strip.text = element_text(face="bold"))+
  labs(y="Proteoform Identifier")+
  scale_y_discrete(name = NULL, expand = expansion(0.02),
                   labels= labelsinplot$feature_name2)
#labels= ABlabels$feature_name2)
p0

temp <-  plot2 %>%
  arrange(desc(lastAA), desc(firstAA)) %>%
  mutate(feature_name = ordered(feature_name, levels= unique(feature_name))) %>%
  mutate(is_modified = mod_str != "") %>%
  mutate(mod_str2 = case_when(mod_str == "" ~ "None", TRUE ~ mod_str)) %>%
  mutate(modification = ordered(.$mod_str2, levels = unique(.$mod_str2))) %>%
  mutate(feature_name = ordered(feature_name, levels= unique(feature_name))) %>% # ordering proteoforms according to N- and C- termini
  mutate(ymin = as.numeric(feature_name) - 0.2,
         ymax = as.numeric(feature_name) + 0.2) 

labelsinplot2 <- temp %>%
  arrange(desc(lastAA), desc(firstAA)) %>%
  mutate(feature_name = ordered(feature_name, levels= unique(feature_name))) %>%
  left_join(labels)

#for abeta   
# temp <- temp %>%
#   mutate(refAA = case_when(terminus == "firstAA" ~ 1, TRUE ~ 40))  

mod_cols <- c("grey60",'#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999', '#00ACC1', '#673AB7')

library(ggplot2)
library(gridExtra)
library(grid)

rect_data <- data.frame(
  xmin = c(1, 19, 134, 228, 272, 322, 358, 413) - offset,
  xmax = c(18, 94, 225, 260, 319, 339, 376, 456) - offset,
  fill_label = c("SP", "Vasostatin-1", "EA-92", "ES-43", "Pancreastatin", "SS-18", "LF-19", "GR-44"),
  fill_color  = c("gray", "blue", "red", "green", "purple", "black", "yellow", "orange")
)

# Ensure that fill_label in rect_data is a factor to maintain order
rect_data$fill_label <- factor(rect_data$fill_label, levels = rect_data$fill_label)

rect_data$"Fragments of CHGA" <- rect_data$fill_label 

#need to adjust to SP length so plot starts at -SP
p3 <- temp %>%
  mutate(facet_var = "Proteoforms of CHGA (P10645)") %>%
  mutate(
    previousAA = str_extract(Proteoform, "^[A-Z](?=\\.)"),  # Extract letter before the first "."
    nextAA = str_extract(Proteoform, "(?<=\\.)([A-Z])$")    # Extract letter after the last "."
  ) %>%
  mutate(firstAA = firstAA - offset) %>%
  mutate(lastAA = lastAA - offset) %>%
  ggplot() +
  # Add the first set of rectangles using rect_data
  # geom_rect(data = rect_data, 
  #           aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = `Fragments of CHGA`), 
  #           alpha = 0.2, inherit.aes = FALSE, show.legend = TRUE) +
  geom_vline(data = rect_data, aes(xintercept = xmin, color = `Fragments of CHGA`), size = linesize, alpha = linealpha, linetype = "dashed", inherit.aes = FALSE, show.legend = TRUE) +
  geom_vline(data = rect_data, aes(xintercept = xmax, color = `Fragments of CHGA`), size = linesize, alpha = linealpha, linetype = "dashed", inherit.aes = FALSE, show.legend = FALSE) +
  scale_colour_colorblind()+
  new_scale_fill() + # Add new scale for the next fill aesthetic
  new_scale_color() + # Add new scale for the next color aesthetic
  geom_rect(aes(ymin = ymin, ymax = ymax,
                xmin = firstAA, xmax = lastAA,
                fill = modification), 
            stat = "identity", color = "black", size = 0.1) +
  facet_grid(cols = vars(facet_var)) +
  # First fill scale for rect_data
  scale_fill_manual(
    name = "Fragments of CHGA",
    values = setNames(rect_data$fill_color, rect_data$`Fragments of CHGA`),
    guide = guide_legend(order = 1)
  ) +
  # Second fill scale for main data
  scale_fill_manual(
    name = "Modifications of CHGA",
    values = mod_cols,
    guide = guide_legend(order = 2)
  ) +
  scale_y_continuous(breaks = seq_along(temp$feature_name),
                     labels = temp$feature_name,
                     expand = expansion(0.015), position = "right") +
  annotation_custom(grob = textGrob(label = "Residue",
                                    gp = gpar(fontsize = textsize),
                                    y = -30, default.units = "pt")) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = textsize) +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"),
        legend.position = "right",
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  scale_y_continuous(breaks = seq_along(temp$feature_name),
                     labels = labelsinplot2$feature_name2,
                     expand = expansion(0.015), position = "right")

p3

library(patchwork)
combined_plot_chga <- p0 + p3 +
  plot_layout(widths = c(0.5, 1),
              guides = "collect")
combined_plot_chga


# #combine and save plots 
# combined_plot <- combined_plot_sst / combined_plot_chga+
#   plot_annotation(tag_levels = list('A'))& 
#   theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
# combined_plot
# 
# ggsave(plot = combined_plot, filename= "FigureS_sst_chga.png", 
#        scale=2.2, 
#        width = 170,
#        height = 170,
#        dpi = 800,
#        units = c("mm"))
# 


#CHGB
###########################################
###########################################
###########################################
###########################################
###########################################
###########################################

load("msnset_humanislet_sc_JMFclustering_wmodanno.RData")

x <- exprs(m) %>%
  as.data.frame()

meta <- m@featureData@data

meta_short <- meta %>%
  mutate(proteoform_id = PF) %>%
  dplyr::select(Gene, proteoform_id, UniProtAcc, mass)

x_long <- x %>%
  rownames_to_column(var = "proteoform_id") %>%
  remove_rownames() %>%
  pivot_longer(-proteoform_id, names_to = "SubjectID",
               values_to = "sc") %>%
  inner_join(.,meta_short)

#all rows of n should be 12 if 0s are filled correctly 
x_long %>%
  group_by(proteoform_id) %>%
  tally() %>%
  View()

#summarize by column sc across all 12 samples 
x_long %>%
  group_by(proteoform_id) %>%
  summarize(meansc = mean(sc),
            sdsc = sd(sc)) %>%
  View()

# plot <- fData(m) %>%
#   dplyr::mutate(cleanSeq = gsub("\\[.+?\\]",
#                                 "", Proteoform), cleanSeq = gsub("\\(|\\)", "", cleanSeq),
#                 cleanSeq = sub("^[A-Z]?\\.(.*)\\.[A-Z]?", "\\1", cleanSeq)) %>%
#   filter(!str_detect(Proteoform, "-57.02")) %>% #remove a missed alkylation, which is about same mass as Arg->Val
#   filter(!str_detect(Proteoform, "150.0")) %>% #remove DTT adduct 
#   # filter(!str_detect(mods, "NH3")) %>% 
#   filter(Gene == "CHGB") %>%
#   arrange(desc(count)) %>%  # Sorts the data in descending order based on spectral_counts
#   slice_head(n = 30)   # Selects the top 20 rows
# 
# plot %>%
#   select(-mods) %>%
#   mutate(feature_name = PF) %>%
#   select(count, firstAA, lastAA, Proteoform, feature_name, PF, Gene, pcGroup,cleanSeq) %>%
#   write.csv(file="unadjusted_forchgbplot.csv")

# save this csv and manually edit to make a plot below

load("x_recluster_oct2024_sc.RData")

#polished plot 

plot2 <- read.csv("adjusted_forchgbplot.csv")

#set offset for SP, can change back to 0 if you want to rescale plots 
# offset <- 18
offset <- 0

labels <- plot2 %>%
  distinct(feature_name ,.keep_all = TRUE) %>%
  mutate(is_modified = if_else(mod_str == "", FALSE, TRUE)) %>% 
  mutate(firstAA = firstAA-offset) %>%
  mutate(lastAA =lastAA-offset) %>%
  mutate(feature_name2 = paste0("Chgb","(",firstAA,"-", lastAA,")", if_else(is_modified  == TRUE,"*",""))) %>%
  mutate(feature_name3 = fct_rev(feature_name)) %>%
  select(feature_name, feature_name2) %>%
  arrange((feature_name))

#need to make a new df so you can reference it with "$", be careful plotting labels 
labelsinplot <- plot2 %>%
  arrange(desc(lastAA), desc(firstAA)) %>%
  mutate(feature_name = ordered(feature_name, levels= unique(feature_name))) %>%
  mutate(count_variable = paste("Abundance of CHGB")) %>%
  left_join(labels)

p0 <- plot2 %>%
  arrange(desc(lastAA), desc(firstAA)) %>%
  mutate(feature_name = ordered(feature_name, levels= unique(feature_name))) %>%
  mutate(count_variable = paste("Abundance of CHGB")) %>%
  ggplot() +
  aes(y=feature_name, x=count) +
  #before position dodge it was stacking counts multiple times because of different tests, so values looked way higher 
  geom_bar(stat="identity", color = NA, fill = "black", position = "dodge") +
  scale_x_reverse(labels = scales::scientific, breaks = c(0, 1e1, 3e1, 6e1, 2e2)) +
  scale_fill_gradientn(#colors = MSnSet.utils::hot2.colors(50),
    colors = terrain.colors(50),
    name = expression(paste("log"[10],"(count)"))) +
  # guides(size = guide_legend(title = expression(paste("-log"[10],"(adj p-value)")))) +
  facet_grid(~ count_variable, scales = "free_x", drop = T) + #seems to just be adding header to top, also allows for alignment at end? 
  theme_bw(base_size=16) +
  theme(axis.text.y.right = element_text(hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "left",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_blank(), #removes feature labels
        axis.ticks.y = element_blank(), #removes feature labels
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        # panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background =element_rect(fill="white"),
        strip.text = element_text(face="bold"))+
  labs(y="Proteoform Identifier")+
  scale_y_discrete(name = NULL, expand = expansion(0.02),
                   labels= labelsinplot$feature_name2)
#labels= ABlabels$feature_name2)
p0

temp <-  plot2 %>%
  arrange(desc(lastAA), desc(firstAA)) %>%
  mutate(feature_name = ordered(feature_name, levels= unique(feature_name))) %>%
  mutate(is_modified = mod_str != "") %>%
  mutate(mod_str2 = case_when(mod_str == "" ~ "None", TRUE ~ mod_str)) %>%
  mutate(modification = ordered(.$mod_str2, levels = unique(.$mod_str2))) %>%
  mutate(feature_name = ordered(feature_name, levels= unique(feature_name))) %>% # ordering proteoforms according to N- and C- termini
  mutate(ymin = as.numeric(feature_name) - 0.2,
         ymax = as.numeric(feature_name) + 0.2) 

labelsinplot2 <- temp %>%
  arrange(desc(lastAA), desc(firstAA)) %>%
  mutate(feature_name = ordered(feature_name, levels= unique(feature_name))) %>%
  left_join(labels)

#for abeta   
# temp <- temp %>%
#   mutate(refAA = case_when(terminus == "firstAA" ~ 1, TRUE ~ 40))  

mod_cols <- c('#e41a1c',"grey60",'#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999', '#00ACC1', '#673AB7')

library(ggplot2)
library(gridExtra)
library(grid)

rect_data <- data.frame(
  xmin = c(1, 440, 575, 617) - offset,
  xmax = c(20, 513, 585, 673) - offset,
  fill_label = c("SP", "GAWK", "PE-11", "CCB"),
  fill_color  = c("gray", "blue", "red", "green")
)

# Ensure that fill_label in rect_data is a factor to maintain order
rect_data$fill_label <- factor(rect_data$fill_label, levels = rect_data$fill_label)

rect_data$"Fragments of CHGB" <- rect_data$fill_label 

#need to adjust to SP length so plot starts at -SP
p3 <- temp %>%
  mutate(facet_var = "Proteoforms of CHGB (P05060)") %>%
  mutate(
    previousAA = str_extract(Proteoform, "^[A-Z](?=\\.)"),  # Extract letter before the first "."
    nextAA = str_extract(Proteoform, "(?<=\\.)([A-Z])$")    # Extract letter after the last "."
  ) %>%
  mutate(firstAA = firstAA - offset) %>%
  mutate(lastAA = lastAA - offset) %>%
  ggplot() +
  # Add the first set of rectangles using rect_data
  # geom_rect(data = rect_data, 
  #           aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = `Fragments of CHGB`), 
  #           alpha = 0.2, inherit.aes = FALSE, show.legend = TRUE) +
  geom_vline(data = rect_data, aes(xintercept = xmin, color = `Fragments of CHGB`), size = linesize, alpha = linealpha, linetype = "dashed", inherit.aes = FALSE, show.legend = TRUE) +
  geom_vline(data = rect_data, aes(xintercept = xmax, color = `Fragments of CHGB`), size = linesize, alpha = linealpha, linetype = "dashed", inherit.aes = FALSE, show.legend = FALSE) +
  scale_colour_colorblind()+
  new_scale_fill() + # Add new scale for the next fill aesthetic
  new_scale_color() + # Add new scale for the next color aesthetic
  geom_rect(aes(ymin = ymin, ymax = ymax,
                xmin = firstAA, xmax = lastAA,
                fill = modification), 
            stat = "identity", color = "black", size = 0.1) +
  facet_grid(cols = vars(facet_var)) +
  # First fill scale for rect_data
  # scale_fill_manual(
  #   name = "Fragments of CHGB",
  #   values = setNames(rect_data$fill_color, rect_data$`Fragments of CHGB`),
  #   guide = guide_legend(order = 1)
  # ) +
  # Second fill scale for main data
  scale_fill_manual(
    name = "Modifications of CHGB",
    values = mod_cols,
    guide = guide_legend(order = 2)
  ) +
  scale_y_continuous(breaks = seq_along(temp$feature_name),
                     labels = temp$feature_name,
                     expand = expansion(0.015), position = "right") +
  annotation_custom(grob = textGrob(label = "Residue",
                                    gp = gpar(fontsize = textsize),
                                    y = -30, default.units = "pt")) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = textsize) +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"),
        legend.position = "right",
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  scale_y_continuous(breaks = seq_along(temp$feature_name),
                     labels = labelsinplot2$feature_name2,
                     expand = expansion(0.015), position = "right")

p3

#CHGB plot 
library(patchwork)
combined_plot_chgb <- p0 + p3 +
  plot_layout(widths = c(0.5, 1),
              guides = "collect")

#combine plots and save 
###########################################

combined_plot <- combined_plot_chga / combined_plot_chgb+
  plot_annotation(tag_levels = list('A'))& 
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
combined_plot

ggsave(plot = combined_plot, filename= "FigureS_chga_chgb.png", 
       scale=2.2, 
       width = 170,
       height = 170,
       dpi = 800,
       units = c("mm"))


ss_plot <- combined_plot_sst+
  plot_annotation(tag_levels = list('A'))& 
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
ss_plot

ggsave(plot = ss_plot, filename= "FigureS_sst.png", 
       scale=2.2, 
       width = 170,
       height = 85,
       dpi = 800,
       units = c("mm"))



