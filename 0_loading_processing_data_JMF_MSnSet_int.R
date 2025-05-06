library(tidyverse)
library(TopPICR) 
library(MSnbase)
library(MSnSet.utils)
library(PNNL.DMS.utils)

if("toppic_output_humanislet.RData" %in% list.files()){
  load("toppic_output_humanislet.RData")
}else{
  toppic_output <- read_TopPIC_DMS(5570)
  save(toppic_output, file = "toppic_output_humanislet.RData")
}


ids <- toppic_output$ms2identifications %>%
  filter(!str_detect(Dataset, "TopDown_HumanIslet_E_Nottreated_24hr_1_Aragorn_21Mar24_BMEB-21-09-03")) %>%
  filter(!str_detect(Dataset, "TopDown_HumanIslet_J_Treated_24hr_2_Aragorn_21Mar24_BMEB-21-09-03"))


feat <- toppic_output$ms1features %>%
  filter(!str_detect(Dataset, "TopDown_HumanIslet_E_Nottreated_24hr_1_Aragorn_21Mar24_BMEB-21-09-03")) %>%
  filter(!str_detect(Dataset, "TopDown_HumanIslet_J_Treated_24hr_2_Aragorn_21Mar24_BMEB-21-09-03"))

# remove NA annotations. We can't handle them at the FDR filter tuning step
# Those are non-uniprot IDs, like short ORFs and contaminants
ids <- ids %>%
  dplyr::filter(!is.na(AnnType))

ids <- ids %>%
  mutate(CV = as.character(str_sub(Dataset,-2,-1))) %>%
  mutate(Dataset = as.character(str_sub(Dataset,1,-5))) 
feat <- feat %>%
  mutate(CV = as.character(str_sub(Dataset,-2,-1))) %>%
  mutate(Dataset = as.character(str_sub(Dataset,1,-5))) 


#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
# Identified feature steps
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!

# Remove erroneous genes ---------------
# What does this do? In case a proteoform assigned to multiple genes,
# this function selects the gene with the lowest E-value.
x_err <- rm_false_gene(ids) #rename !!!


# getting protLength (parent protein length)
# No need to match. Just need to link accessions.
# 1. read fasta and 
# 2. link by `Protein accession` or UniProtAcc

fst <- Biostrings::readAAStringSet(
  file.path(r"(\\gigasax\DMS_FASTA_File_Archive\Dynamic\Forward\)",
            "ID_008274_9D2E95FB.fasta"))

add_protein_length <- function(x, fst_obj){
  fst_df <- data.frame(`Protein accession` = names(fst_obj),
                       protLength = BiocGenerics::width(fst_obj), 
                       row.names = NULL, 
                       check.names = FALSE) %>%
    mutate(`Protein accession` = word(`Protein accession`))
  x <- inner_join(x, fst_df)
  return(x)
}


x_err <- add_protein_length(x_err, fst)

compute_fdr(x_err)

# looks like we don't need this step anymore
x_aug <- x_err

# Filter by cleanSeq count ---------------

# really minimal filter
# x_fbc <- filter_by_count(x = x_aug,
#                          count_within = c("Dataset", "Scan(s)"),
#                          count = "cleanSeq",
#                          threshold = 2) # better to switch to 3 and drop the E-value filter altogether.

# FDR control ---------------

# First find the E-value threshold for each of the three annotation types.
# FDR is calculated at Gene level.
the_cutoff <- find_evalue_cutoff(x = x_aug,
                                 fdr_threshold = 0.01)

# Apply the actual FDR control/filter with the threshold values from above.
x_fdr <- apply_evalue_cutoff(x = x_aug,
                             e_vals = the_cutoff)

compute_fdr(x_fdr)


# Proteoform inference ---------------
x_ipf <- infer_prot(x_fdr)


###in progress? 
x_inferred <- set_pf_level(x_ipf)

compute_fdr(x_ipf)

# Align retention time ---------------

# Create the model that will be used to align each retention time. The model is
# created between a reference data set and all other data sets.
the_model <- form_model(
  x = x_inferred,
  ref_ds = find_ref_ds(x = x_inferred),
  control = loess.control(surface = "direct"), # Use direct to avoid NAs.
  span = 0.5,
  family = "symmetric")

# Align retention times according to the model created previously.
x_art <- align_rt(
  x = x_inferred,
  model = the_model,
  var_name = "Retention time"
)


# Recalibrate the mass ---------------

# Calculate the error between the `Precursor mass` and the `Adjusted precursor
# mass`. This acts as the model for recalibrating the mass. A reference data set
# is not used when calculating the mass error model.
x_error <- calc_error(x = x_art,
                      ref_ds = find_ref_ds(x = x_inferred))

x_error$rt_sd <- 300
# Recalibrate the mass according the the errors computed previously.
x_rcm <- recalibrate_mass(x = x_art,
                          errors = x_error,
                          var_name = "Precursor mass") 
# Cluster ---------------
# Key step. But doesn't take too long.
x_cluster <- cluster(x = x_rcm,
                     errors = x_error,
                     method = "single",
                     height = 3,
                     min_size = 2)


# Group clusters ---------------
x_recluster <- create_pcg(x = x_cluster,
                          errors = x_error,
                          n_mme_sd = 5,
                          n_rt_sd = 3,
                          #ppm_cutoff = 2.1,
                          n_Da = 3)

save(x_recluster, file="x_recluster_oct2024_int.RData")

x_meta <- create_mdata(x = x_recluster,
                       errors = x_error,
                       n_mme_sd = 5,
                       n_rt_sd = 3
) %>%
  mutate(PF = paste(Gene, pcGroup, sep = "_"))

save(x_meta, file="x_meta.RData")


# Align unidentified features ---------------

# Align the unidentified feature retention times with the model created by the
# identified feature retention times.
feat_art <- align_rt(
  x = feat,
  model = the_model,
  var_name = "Time_apex")

# Recalibrate unidentified feature mass ---------------

# Recalibrate the unidentified feature masses with the model created by the
# identified feature masses.
feat_rcm <- recalibrate_mass(
  x = feat_art,
  errors = x_error,
  var_name = "Mass")


feat_rtv <- match_features(ms2 = x_recluster,
                           ms1 = feat_rcm,
                           errors = x_error,
                           n_mme_sd = 5,
                           n_rt_sd = 3
)

x <- feat_rtv %>%
  mutate(feature_name = paste(Gene, pcGroup, CV, sep = "_"),
         sample_name = Dataset) %>%
  dplyr::select(-Dataset)

x_expr <- x %>%
  pivot_wider(id_cols = "feature_name",
              names_from = "sample_name",
              values_from = "Intensity") %>%
  as.data.frame() %>%
  {rownames(.) <- .$feature_name;.} %>%
  dplyr::select(-feature_name) %>%
  as.matrix()

# x_expr <- log2(x_expr)
#uncheck if you want log2 

# features
x_feat <- x %>%
  group_by(Gene,pcGroup, CV, feature_name) %>%
  summarize(median_intensity = median(Intensity),
            count = n(),
            .groups = "keep") %>%
  ungroup() %>%
  left_join(x_meta, by=c("Gene","pcGroup")) %>%
  mutate(proteoform_id = paste(Gene, pcGroup, sep="_")) %>%
  as.data.frame() %>%
  {rownames(.) <- .$feature_name;.}

x_pheno <- x %>%
  separate(sample_name, into=c("TD", "islet", "Letter","Treated", "Time", "Biorep"), remove = FALSE) %>%
  distinct(sample_name, Letter, Treated,Time, Biorep) %>%
  as.data.frame() %>%
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
  {rownames(.) <- .$sample_name;.} 

m <- MSnSet(x_expr, x_feat[rownames(x_expr),], x_pheno[colnames(x_expr),])

image_msnset(m)

save(x_feat, file="featuredata_beforerollup.RData")

#if log2
# m <- rrollup(m, "proteoform_id", rollFun = "-", verbose = FALSE)
# save(m, file="msnset_humanislet_log2int.RData")

#if normal scale 
m <- rrollup(m, "proteoform_id", rollFun = "/", verbose = FALSE)

#' recover proteoform info, it all get's lost in rollup step 
x_meta <- x_meta %>%
  mutate(proteoform_id = paste(Gene, pcGroup, sep="_")) %>%
  as.data.frame() %>%
  {rownames(.) <- .$proteoform_id;.}
fData(m) <- x_meta[featureNames(m),]

save(m, file="msnset_humanislet_int_notnormalized.RData")

mlc <- log2_zero_center(m)

save(mlc, file="msnset_humanislet_int_log2center_oct2024_forTyler.RData")

###annotate with mods 

x <- fData(mlc)

#these lines reformat Nterm acetyl from new version ([Acetyl]-aa to old version (aa)[Acetyl ])
x <- x %>% 
  mutate(Proteoform = case_when(str_detect(Proteoform, "\\[Acetyl]-") ~ paste0(substr(gsub(pattern = "\\[Acetyl]-", x = Proteoform, replacement = ""), 1,2), "(", substr(gsub(pattern = "\\[Acetyl]-", x = Proteoform, replacement = ""), 3,3), ")[Acetyl]", substr(gsub(pattern = "\\[Acetyl]-", x = Proteoform, replacement = ""),4,nchar(gsub(pattern = "\\[Acetyl]-", x = Proteoform, replacement = ""))))
                                , !str_detect(Proteoform, "\\[Acetyl]-") ~ Proteoform)) 

unimods <- TopPICR::create_mod_data(
  mod_file = "TopPIC_Dynamic_Mods.txt",
  mod_path = ".",
  use_unimod = TRUE
)

annotate_Nterm_acetyls <- function(x, nterm_tol = 3, acetyl_id = "Acetyl"){
  
  temp_mod_names <- "mods"
  if("mod_names" %in% names(x$mods[[1]]))
    temp_mod_names <- "mod_names"
  
  for(i in 1:nrow(x)){
    if(x[i,"firstAA"] <= nterm_tol){
      y <- x[i,"mods"][[1]]
      if(length(y$mods) > 0){
        if(y$mods[1] == acetyl_id & y$mods_left_border[1] == 1){
          y$mods[1] <- paste("N-", acetyl_id, sep="")
          y[[temp_mod_names]][1] <- y$mods[1]
          posi_str <- map2_chr(y$mods_left_border, y$mods_right_border, paste, sep="-")
          y$mods_str <- paste(map2_chr(y[[temp_mod_names]], posi_str, paste, sep="@"), collapse=", ")
          x[i,"mods"][[1]] <- list(y)
        }
      }
    }
  }
  return(x)
}


x <- x %>%
  mutate(mods = map(Proteoform, TopPICR:::extract_mods))

mass_annotation_table <- TopPICR:::get_mass_annotation_table(x, unimods, 0.3, 0.1)

x <- x %>%
  mutate(mods = map(mods, TopPICR:::annotate_masses, mass_annotation_table, matching_tol = .Machine$double.eps))

x <- annotate_Nterm_acetyls(as.data.frame(x), nterm_tol = 3, acetyl_id = "Acetyl")
rownames(x) <- x$proteoform_id

####need to switchh from m to m1 if int or sc 

fData(mlc) <- x

save(mlc, file= "msnset_humanislet_int_log2center_oct2024_forTyler_wmodanno.RData")



# load("msnset_humanislet_log2int.RData")

#old volcano plots pipeline before Tyler did paired data analysis, does work but data isn't used in paper 

# library(proDA)
# 
# x_final <- as.data.frame(m) %>%
#   tibble::rownames_to_column(var = "Dataset") %>%
#   filter(Dataset != "TopDown_HumanIslet_E_Nottreated_24hr_1_Aragorn_21Mar24_BMEB-21-09-03") %>%
#   filter(Dataset != "TopDown_HumanIslet_J_Treated_24hr_2_Aragorn_21Mar24_BMEB-21-09-03") %>%
#   column_to_rownames(var = "Dataset") %>%
#   t() %>%
#   median_normalization() 
# 
# rowwise_na_summary <- apply(x_final, 1, function(x) sum(!is.na(x))) %>%
#   as.data.frame() %>%
#   rownames_to_column(var = "PF")
# 
# rowwise_na_summary$notna <- rowwise_na_summary$.
# 
# rowwise_na_summary <- rowwise_na_summary %>%
#   mutate(MissingPercent = (notna/ncol(x_final)))
# 
# x_final[,1:12] <- 2^(x_final[,1:12])
# 
# x_final2 <- x_final %>%
#   as.data.frame() %>%
#   rownames_to_column(var = "PF") %>%
#   mutate(name = PF) %>%
#   mutate(ID = PF) %>%
#   group_by(name) %>%
#   add_count(name  = "n") %>%
#   filter(n == 1) %>%
#   dplyr::select(-n) 
# 
# LFQ_columns <- grep("_",colnames(x_final2)) # get LFQ column numbers
# 
# 
# 
# #for some reason this spot fails but if you rerun it ~3 times it works?????
# sampleIDs <- colnames(x_final2)[LFQ_columns] %>%
#   as.data.frame() %>%
#   dplyr::rename(Dataset = ".")  %>%
#   filter(grepl("_", Dataset)) %>%
#   mutate(Group = case_when(grepl("Nottreated", Dataset) ~ "Untreated",
#                            grepl("Treated", Dataset) ~ "Treated")) %>%
#   mutate(Pair = case_when(grepl("_A_", Dataset) ~ "Pair1",
#                            grepl("_H_", Dataset) ~ "Pair1",
#                            grepl("_B_", Dataset) ~ "Pair2",
#                            grepl("_D_", Dataset) ~ "Pair2",
#                            grepl("_C_", Dataset) ~ "Pair3",
#                            grepl("_N_", Dataset) ~ "Pair3",
#                            grepl("_E_", Dataset) ~ "Pair4",
#                            grepl("_J_", Dataset) ~ "Pair4",
#                            grepl("_F_", Dataset) ~ "Pair5",
#                            grepl("_K_", Dataset) ~ "Pair5",
#                            grepl("_G_", Dataset) ~ "Pair6",
#                            grepl("_M_", Dataset) ~ "Pair6",
#                            grepl("_I_", Dataset) ~ "Pair7",
#                            grepl("_L_", Dataset) ~ "Pair7"))
# 
# experimental_design <- sampleIDs %>%
#   mutate(label = Dataset) %>%
#   mutate(condition = Group) %>%
#   group_by(condition) %>%
#   # mutate(replicate = dense_rank(desc(label))) %>%
#   mutate(replicate = dense_rank((label))) %>%
#   ungroup() %>%
#   dplyr::select(label, condition, replicate, Pair) %>%
#   arrange(label)
# library(DEP)
# 
# data_se <- make_se(x_final2, LFQ_columns, experimental_design) 
# 
# data_diff <- DEP::test_diff(data_se, type = "all") 
# 
# 
# ##########!!!! Note that DEP adjusted p.values are NOT CORRECT. The below code
# ######### is meant to properly adjust with BH correction
# 
# data_diff@elementMetadata@listData <- data_diff@elementMetadata@listData %>%
#   as.data.frame() %>%
#   mutate_at(vars(matches("p.val")), ~ p.adjust(.x, method = "BH")) %>%
#   dplyr::select(-contains("p.adj"),
#                 -contains("CI.L"),
#                 -contains("CI.R")) %>%
#   rename_with(.fn = ~ str_replace(., "p.val", "p.adj"),
#               .cols = contains("p.val")) %>%
#   as.list() 
# 
# dep <- add_rejections(data_diff, alpha = 0.05, lfc = 1)
# 
# 
# data_results <- as.data.frame(dep@elementMetadata@listData)
# 
# 
# library(EnhancedVolcano)
# 
# dep <- add_rejections(data_diff, alpha = 0.05, lfc = 1)
# 
# 
# data_results <- as.data.frame(dep@elementMetadata@listData)
# 
# 
# library(EnhancedVolcano)
# 
# cs_x <- data_results %>%
#   mutate(adj.pvalue = Treated_vs_Untreated_p.adj) %>%
#   mutate(log2FC = Treated_vs_Untreated_diff) %>%
#   dplyr::select(PF, adj.pvalue, log2FC) %>%
#   filter(!is.na(adj.pvalue)) 
# 
# rownames(cs_x) <- cs_x$PF 
# 
# cs_x <- cs_x %>%
#   inner_join(x_meta) %>%
#   inner_join(rowwise_na_summary)
# 
# cs_x <- cs_x %>%
#   filter(notna > 2)
# 
# ##########filter for proteoforms that are in at least 3 out of 13 samples 
# 
# 
# ########this section annotates mods, treats cs_x like fData(m)
# 
# cs_xanno <- cs_x
# 
# #these lines reformat Nterm acetyl from new version ([Acetyl]-aa to old version (aa)[Acetyl ])
# cs_xanno <- cs_xanno %>% 
#   mutate(Proteoform = case_when(str_detect(Proteoform, "\\[Acetyl]-") ~ paste0(substr(gsub(pattern = "\\[Acetyl]-", x = Proteoform, replacement = ""), 1,2), "(", substr(gsub(pattern = "\\[Acetyl]-", x = Proteoform, replacement = ""), 3,3), ")[Acetyl]", substr(gsub(pattern = "\\[Acetyl]-", x = Proteoform, replacement = ""),4,nchar(gsub(pattern = "\\[Acetyl]-", x = Proteoform, replacement = ""))))
#                                 , !str_detect(Proteoform, "\\[Acetyl]-") ~ Proteoform)) 
# 
# unimods <- TopPICR::create_mod_data(
#   mod_file = "TopPIC_Dynamic_Mods.txt",
#   mod_path = ".",
#   use_unimod = TRUE
# )
# 
# annotate_Nterm_acetyls <- function(x, nterm_tol = 3, acetyl_id = "Acetyl"){
#   
#   temp_mod_names <- "mods"
#   if("mod_names" %in% names(x$mods[[1]]))
#     temp_mod_names <- "mod_names"
#   
#   for(i in 1:nrow(x)){
#     if(x[i,"firstAA"] <= nterm_tol){
#       y <- x[i,"mods"][[1]]
#       if(length(y$mods) > 0){
#         if(y$mods[1] == acetyl_id & y$mods_left_border[1] == 1){
#           y$mods[1] <- paste("N-", acetyl_id, sep="")
#           y[[temp_mod_names]][1] <- y$mods[1]
#           posi_str <- map2_chr(y$mods_left_border, y$mods_right_border, paste, sep="-")
#           y$mods_str <- paste(map2_chr(y[[temp_mod_names]], posi_str, paste, sep="@"), collapse=", ")
#           x[i,"mods"][[1]] <- list(y)
#         }
#       }
#     }
#   }
#   return(x)
# }
# 
# 
# cs_xanno <- cs_xanno %>%
#   mutate(mods = map(Proteoform, TopPICR:::extract_mods))
# 
# mass_annotation_table <- TopPICR:::get_mass_annotation_table(cs_xanno, unimods, 0.3, 0.1)
# 
# cs_xanno <- cs_xanno %>%
#   mutate(mods = map(mods, TopPICR:::annotate_masses, mass_annotation_table, matching_tol = .Machine$double.eps))
# 
# cs_xanno <- annotate_Nterm_acetyls(as.data.frame(cs_xanno), nterm_tol = 3, acetyl_id = "Acetyl")
# rownames(cs_xanno) <- cs_xanno$proteoform_id
# 
# save(sampleIDs, cs_xanno, x_final, file="results_5570.RData")
# 
# ############## Here are a few candidate proteoforms that are found exclusively
# ###### in one condition or the other
# x_final_long <- x_final2 %>%
#   ungroup() %>%
#   dplyr::select(-name, -ID) %>%
#   pivot_longer(-PF, names_to = "Dataset",
#                values_to = "Intensity") %>%
#   full_join(., sampleIDs) %>%
#   filter(!is.na(Intensity)) %>%
#   group_by(PF, Group) %>%
#   add_count(name = "n") %>%
#   group_by(PF) %>%
#   add_count(name = "n2") %>%
#   ungroup() %>%
#   mutate(n3 = n2-n) %>%
#   distinct(PF, Group, n,n2,n3) %>%
#   filter(n >=5 & n3 == 0) %>%
#   inner_join(., x_meta)
