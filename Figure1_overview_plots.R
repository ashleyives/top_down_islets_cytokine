
library(tidyverse)
library(TopPICR) 
library(MSnbase)
library(MSnSet.utils)
library(PNNL.DMS.utils)

library(tidyverse)
library(EnhancedVolcano)
library(EnhancedVolcano)
library(ggplot2)
library(ggrepel)

#cs_xanno
# nm = load("results_5570_unadjusted.RData")

load("x_meta.RData")
load("msnset_humanislet_log2int.RData")

textsize <- 18 

x <- exprs(m) %>%
  as.data.frame()

x_meta <- x_meta %>%
  mutate(proteoform_id = PF)

meta <- m@featureData@data

meta_short <- meta %>%
  left_join(x_meta) %>%
  dplyr::select(Gene, proteoform_id, UniProtAcc, mass)

x_long <- x %>%
  rownames_to_column(var = "proteoform_id") %>%
  remove_rownames() %>%
  pivot_longer(-proteoform_id, names_to = "SubjectID",
               values_to = "Intensity") %>%
  inner_join(.,meta_short)

#histogram of masses as SI figure 
mean_mass <- mean(x_meta$mass, na.rm = TRUE)
median_mass <- median(x_meta$mass, na.rm = TRUE)

# Creating the histogram of masses 
hist <- ggplot(x_meta, aes(x = mass)) + 
  geom_histogram(binwidth = 500, color = "black", fill = "#999999", alpha = 0.7) +
  geom_vline(aes(xintercept = mean_mass, color = "Mean"), linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = median_mass, color = "Median"), linetype = "dashed", size = 1) +
  scale_color_manual(name = "Statistics", values = c(Mean = "red", Median = "blue")) +
  labs(title = "Histogram of monoisotopic masses", x = "Monoisotopic mass (Da)", y = "Frequency") +
  theme_bw(base_size = textsize) +
  theme( panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())+
  scale_x_continuous(expand=c(0, 0))+
  scale_y_continuous(expand=c(0, 0), limits = c(0,250))

ggsave(plot = hist, filename= "FigureS_mass_hist.png", scale=1,
       width = 170,
       height = 100,
       dpi = 800,
       units = c("mm"))



#data completeness insulin only 
holes_pfr_ins <- x_long %>%
  filter(str_detect(proteoform_id, "INS")) %>%
  filter(!is.na(Intensity)) %>%
  group_by(proteoform_id) %>%
  tally() %>%
  mutate(n = n/length(unique(x_long$SubjectID))) %>%
  mutate(Group = case_when(between(n, 0, 0.25) ~ "25-0%",
                           between(n, 0.25, 0.5) ~ "50-25%",
                           between(n, 0.5, 0.75) ~ "75-50%",
                           between(n, 0.75, 1) ~ "100-75%")) %>%
  group_by(Group) %>%
  tally() 


#data completeness
holes_pfr <- x_long %>%
  filter(!is.na(Intensity)) %>%
  group_by(proteoform_id) %>%
  tally() %>%
  mutate(n = n/length(unique(x_long$SubjectID)) ) %>%
  mutate(Group = case_when(between(n, 0, 0.25) ~ "25-0%",
                           between(n, 0.25, 0.5) ~ "50-25%",
                           between(n, 0.5, 0.75) ~ "75-50%",
                           between(n, 0.75, 1) ~ "100-75%")) %>%
  group_by(Group) %>%
  tally() %>%
  ggplot()+
  aes(x= fct_relevel(Group,"100-75%",
                     "75-50%", "50-25%", "25-0%" ),
      y = n, fill = fct_relevel(Group,"100-75%",
                                "75-50%", "50-25%", "25-0%"))+
  geom_bar(stat= "identity")+
  theme_bw(base_size = textsize) +
  theme(legend.position = "none",
        panel.background = element_rect(fill= 'white'),
        axis.text.y=element_text(color = 'black', size = textsize),
        axis.text.x=element_text(angle = 45,
                                 vjust = 1,
                                 hjust= 1,
                                 color='black',
                                 size = textsize),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        text=element_text(family="Helvetica"),
        axis.line = element_line())+
  scale_fill_brewer(palette = "Blues", direction = -1)+
  xlab("Data Completeness (Proteoform)")+
  ylab("Total proteoforms(n)")+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) #force to start at 0 
holes_pfr

#data completeness with insulin marked 
# holes_pfr_withins <- x_long %>%
#   filter(!is.na(Intensity)) %>%
#   group_by(proteoform_id) %>%
#   tally() %>%
#   mutate(n = n/length(unique(x_long$SubjectID)) ) %>%
#   mutate(Group = case_when(between(n, 0, 0.25) ~ "25-0%",
#                            between(n, 0.25, 0.5) ~ "50-25%",
#                            between(n, 0.5, 0.75) ~ "75-50%",
#                            between(n, 0.75, 1) ~ "100-75%")) %>%
#   group_by(Group) %>%
#   tally() %>%
#   left_join(holes_pfr_ins, by = "Group", suffix = c("_all", "_ins")) %>%
#   ggplot() +
#   aes(x = fct_relevel(Group, "100-75%", "75-50%", "50-25%", "25-0%")) +
#   geom_bar(aes(y = n_all, fill = fct_relevel(Group, "100-75%", "75-50%", "50-25%", "25-0%")), 
#            stat = "identity") +
#   geom_segment(aes(x = as.numeric(fct_relevel(Group, "100-75%", "75-50%", "50-25%", "25-0%")) - 0.4, 
#                    xend = as.numeric(fct_relevel(Group, "100-75%", "75-50%", "50-25%", "25-0%")) + 0.4, 
#                    y = n_ins, 
#                    yend = n_ins), 
#                colour = "black", size = 1.2, linetype = "dashed") +
#   # geom_text(aes(x = Group, y = n_ins + 50, label = "Insulin"), # Add text annotation here
#   #           colour = "black", size = 6)+
#   theme_bw(base_size = textsize) +
#   theme(legend.position = "none",
#         panel.background = element_rect(fill = 'white'),
#         axis.text.y = element_text(color = 'black', size = 22),
#         axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, color = 'black', size = 22),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.border = element_blank(),
#         text = element_text(family = "Helvetica"),
#         axis.line = element_line()) +
#   scale_fill_brewer(palette = "Blues", direction = -1) +
#   xlab("Data Completeness (Proteoform)") +
#   ylab("Total proteoforms(n)") +
#   scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) 
# holes_pfr_withins


#these are averages 
pf_count <- x_long %>%
  filter(!is.na(Intensity)) %>%
  group_by(SubjectID) %>%
  tally() %>%
  summarize(mean = mean(n),
            SD = sd(n)) %>%
  mutate(Type = "Proteoforms")

#these are averages 
pf_count_ins <- x_long %>%
  filter(str_detect(proteoform_id, "INS")) %>%
  filter(!is.na(Intensity)) %>%
  group_by(SubjectID) %>%
  tally() %>%
  summarize(mean = mean(n),
            SD = sd(n)) %>%
  mutate(Type = "Proteoforms")

#these are averages
gene_count <- x_long %>%
  filter(!is.na(Intensity)) %>%
  distinct(Gene, SubjectID) %>%
  group_by(SubjectID) %>%
  tally() %>%
  summarize(mean = mean(n),
            SD = sd(n)) %>%
  mutate(Type = "Genes")

combo <- full_join(pf_count, gene_count)

proteoforms_position <- which(levels(factor(combo$Type)) == "Proteoforms")

panela <- combo %>%
  ggplot()+
  aes(x = Type, y = mean, fill = Type)+
  geom_bar(stat = "identity")+
  theme_bw(base_size = textsize) +
  theme(legend.position = "none",
        panel.background = element_rect(fill= 'white'),
        axis.text.y=element_text(color = 'black', size = textsize),
        axis.text.x=element_text(angle = 45,
                                 vjust = 1, 
                                 hjust= 1, 
                                 color='black',
                                 size = textsize),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        text=element_text(family="Helvetica"),
        axis.line = element_line())+
  xlab("")+
  scale_fill_manual(values = c("gray", "lightblue"))+
  geom_errorbar( aes(x=Type,
                     ymin=mean-SD, ymax=mean+SD),
                 width=0.4, colour="black", alpha=1, size=1.2)+
  ylab("Mean Observations (n)")+
  # geom_point(aes(x = "Proteoforms", y = pf_count_ins$mean),  shape = 16, size = 4, colour = "black") +#number of ins proteoforms
  # geom_errorbar(aes(x = "Proteoforms", ymin = pf_count_ins$mean - pf_count_ins$SD, ymax = pf_count_ins$mean + pf_count_ins$SD), 
  #               width = 0.2, colour = "black", size = 1) +
  # geom_text(aes(x = "Proteoforms", y = pf_count_ins$mean + 100, label = "Insulin"), # Add text annotation here
  #           colour = "black", size = 6)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1000), breaks = seq(0, 1000, by = 250)) #force to start at 0 
  # geom_segment(aes(x = proteoforms_position - 0.4, 
  #                  xend = proteoforms_position + 0.4, 
  #                  y = pf_count_ins$mean, yend = pf_count_ins$mean),
  #              colour = "black", size = 1.2, linetype = "dashed") 
panela
#ggsave(plot = panela, filename = "C:/Users/ives435/OneDrive - PNNL/data_sandbox/Top_Down/ROSMAP103_4kDa/figures_manuscript/barplot_genes_pfs.png", width = 5, height = 7)  

unique_genes <- x_long %>%
  filter(!is.na(Intensity)) %>%
  distinct(Gene) %>% 
  tally() 

unique_pfr <- x_long %>%
  filter(!is.na(Intensity)) %>%
  distinct(proteoform_id) %>% 
  tally() 


unique_pfr_ins <- x_long %>%
  filter(str_detect(proteoform_id, "INS")) %>%
  filter(!is.na(Intensity)) %>%
  distinct(proteoform_id) %>% 
  tally() 

unique_acc <- x_long %>%
  filter(!is.na(Intensity)) %>%
  distinct(UniProtAcc) %>% 
  tally() 

proteoforms_position <- which(levels(factor(combo$Type)) == "Proteoforms")

panelb <- combo %>%
  mutate(mean = c(as.numeric(unique_pfr),as.numeric(unique_genes))) %>%
  ggplot()+
  aes(x = Type, y = mean, fill = Type)+
  geom_bar(stat = "identity")+
  theme_bw(base_size = textsize) +
  theme(legend.position = "none",
        panel.background = element_rect(fill= 'white'),
        axis.text.y=element_text(color = 'black', size = textsize),
        axis.text.x=element_text(angle = 45,
                                 vjust = 1, 
                                 hjust= 1, 
                                 color='black',
                                 size = textsize),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        text=element_text(family="Helvetica"),
        axis.line = element_line())+
  xlab("")+
  scale_fill_manual(values = c("gray", "lightblue"))+
  ylab("Distinct Observations (n)")+
  # geom_segment(aes(x = proteoforms_position - 0.4, 
  #                  xend = proteoforms_position + 0.4, 
  #                  y = unique_pfr_ins$n, yend = unique_pfr_ins$n),
  #              colour = "black", size = 1.2, linetype = "dashed") +
  # geom_text(aes(x = "Proteoforms", y = unique_pfr_ins$n + 100, label = "Insulin"), # Add text annotation here
  #           colour = "black", size = 6)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) #force to start at 0 
panelb
#ggsave(plot = panelb, filename = "C:/Users/ives435/OneDrive - PNNL/data_sandbox/Top_Down/ROSMAP103_4kDa/figures_manuscript/unique_obs.png", width = 5, height = 7)  

library(patchwork)

# combined_plot <- (panela + panelb) / holes_pfr_withins + 
#   plot_layout(widths = c(1, 1, 2))+ 
#   plot_annotation(tag_levels = 'A')& 
#   theme(
#         plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"))
#   # theme(plot.tag.position = c(0, 1),
#   #       plot.tag = element_text(size = 8, hjust = 0, vjust = 0))

# combined_plot <- panela + panelb + holes_pfr +
#   plot_layout(widths = c(1, 1, 2))+
#   plot_annotation(tag_levels = 'A')&
#   theme(
#         plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"))
#   # theme(plot.tag.position = c(0, 1),
#   #       plot.tag = element_text(size = 8, hjust = 0, vjust = 0))

# Create a blank plot
blank_plot <- ggplot() + 
  theme_void(base_size = textsize) +  # Ensure it's completely empty
  labs(tag = "A")  # Add a tag

# Combine plots
combined_plot <- (blank_plot / (panela + panelb + holes_pfr + plot_layout(widths = c(1, 1, 2)))) +
  plot_layout(heights = c(2, 1)) +  # Specify heights to make the blank plot x as tall
  plot_annotation(tag_levels = list('A')) &  # Use custom tags so that blank plot tag remains 'A'
  theme(
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
  )

# Print the combined plot
print(combined_plot)

#manually saved 
# ggsave(plot = combined_plot, filename= "Figure1_overview.png", scale=2.2, 
#        width = 170,
#        height = 130,
#        dpi = 800,
#        units = c("mm"))
