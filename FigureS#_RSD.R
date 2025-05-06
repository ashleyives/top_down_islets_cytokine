library(dplyr)
library(tibble)
library(MSnSet.utils)
library(proteomicsCV)
library(tidyverse)

#not log transformed data 

# nm = load("msnset_humanislet_int_notnormalized.RData")
# 
# mcontrol <- m[,(str_detect(sampleNames(m), "Nottreated"))] 
#   
# mtreatment <- m[,(str_detect(sampleNames(m), "Treated"))] 
# 
# cvscontrol <- protCV(data.frame(exprs(mcontrol))) %>%
#   data.frame() %>%
#   mutate(Treatment = "Control") %>%
#   mutate(Normalization = "Unnormalized")
# 
# cvstreatment <- protCV(data.frame(exprs(mtreatment))) %>%
#   data.frame() %>%
#   mutate(Treatment = "Treated")%>%
#   mutate(Normalization = "Unnormalized")
# 
# norm_coeffs <- apply(exprs(m), 2, median, na.rm = T) #calculates the median of each column in the expression matrix while ignoring any NA values.
# exprs(m) <- sweep(exprs(m), 2, norm_coeffs, FUN = "/") #applies normalization to msnset via sweep
# 
# mcontrol <- m[,(str_detect(sampleNames(m), "Nottreated"))] 
# 
# mtreatment <- m[,(str_detect(sampleNames(m), "Treated"))] 
# 
# cvscontrol2 <- protCV(data.frame(exprs(mcontrol))) %>%
#   data.frame() %>%
#   mutate(Treatment = "Control") %>%
#   mutate(Normalization = "Median normalized")
# 
# cvstreatment2 <- protCV(data.frame(exprs(mtreatment))) %>%
#   data.frame() %>%
#   mutate(Treatment = "Treated")%>%
#   mutate(Normalization = "Median normalized")
# 
# plot_data <- rbind(cvscontrol, cvstreatment, cvscontrol2, cvstreatment2)
# 
# ggplot(plot_data, aes(x = Treatment, y = ., fill = Treatment)) +
#   geom_violin() +
#   stat_summary(fun.data = function(y) data.frame(y = median(y), label = paste0("Median: ", round(median(y), 2))),
#                geom = "point", shape = 21, size = 2, color = "black", fill = "red") +
#   stat_summary(fun.data = function(y) data.frame(y = median(y), label = round(median(y), 2)),
#                geom = "text", vjust = -0.5) +
#   labs(x = "Treatment", y = "RSD %") +
#   theme_bw(base_size = 24)+
#   facet_wrap(~Normalization)

#trying RSD via James method, distribution looks the same 
nm = load("msnset_humanislet_int_notnormalized.RData")

exprs(m) <- log2(exprs(m))

norm_coeffs <- apply(exprs(m), 2, median, na.rm = T) #calculates the median of each column in the expression matrix while ignoring any NA values.
exprs(m) <- sweep(exprs(m), 2, norm_coeffs, FUN = "-") #applies normalization to msnset via sweep, because it's log transformed us Fun = "_"

exprs(m) <- 2^(exprs(m))

data <- as.data.frame((exprs(m)))

# Subset the columns based on their names containing "Treated" or "Untreated"
treated <- data[, grepl("Treated", colnames(data))]
untreated <- data[, grepl("Nottreated", colnames(data))]

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


  