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

save(x_recluster, file="x_recluster_oct2024_sc.RData")

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


x <- x_recluster %>%
  mutate(feature_name = paste(Gene, pcGroup, CV, sep = "_"),
         sample_name = Dataset) %>%
  dplyr::select(-Dataset)

#uncheck if you want log2 

#counts all PSMs across CVs, does not require rrollup 
#not derived from x, derived from x_recluster, but
x_expr <- x %>%
  mutate(feature_name = paste(Gene, pcGroup, sep="_")) %>% #added this line, not sure if this works correctly 
  select(feature_name, sample_name, `Scan(s)`) %>%
  pivot_wider(values_from = `Scan(s)`,
              names_from = sample_name,
              values_fn = length,
              values_fill = 0) %>%
  as.data.frame() %>%
  {rownames(.) <- .$feature_name;.} %>%
  select(-feature_name) %>%
  as.matrix()

# x_expr <- log2(x_expr)
#uncheck if you want log2 

# features
x_feat <- x %>%
  group_by(Gene, pcGroup) %>% #feature name contains CV as well so it messes up rownames
  dplyr::summarize(count = n(), .groups = "keep") %>%
  ungroup() %>%
  left_join(x_meta, by=c("Gene","pcGroup")) %>%
  mutate(proteoform_id = paste(Gene, pcGroup, sep="_")) %>%
  as.data.frame() %>%
  column_to_rownames(var ="proteoform_id")
  # {rownames(.) <- .$feature_name;.}

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

save(m, file="msnset_humanislet_sc_JMFclustering.RData")

###annotate with mods 

x <- fData(m)

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
# rownames(x) <- x$proteoform_id

####need to switchh from m to m1 if int or sc 

fData(m) <- x

save(m, file= "msnset_humanislet_sc_JMFclustering_wmodanno.RData")


