# Load panel sequencing data and extract the gene list
msk505 <- readLines("GENIE_synapse_16.1_all/data_gene_panel_MSK-IMPACT505.txt")
gene_list_msk505 <- msk505[3]
msk505 <- unlist(strsplit(gsub("gene_list:\\s*", "", gene_list_msk505), "\\s+")) 

# Simulate artificial WES samples
# Calculate total mutations per sample
total_mutations_sbs <- colSums(mut_matrix_sigma_filtered) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample")

# Add a Status column that will show the ERCC2 status (MUT or WT)
total_mutations_sbs <- total_mutations_sbs %>%
  mutate(Status = ifelse(Sample %in% ercc2_helicase_samples, "ERCC2-MUT", "ERCC2-WT"))
colnames(total_mutations_sbs) <- c("Sample", "TotalMutations", "Status")

# Plot the distirbution
ggplot(total_mutations_sbs, aes(x = Status, y = TotalMutations, fill = Status)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Total Mutations Across Samples (SBS) - TCGA-BLCA WES",
    x = "Sample Status",
    y = "Total Number of Mutations",
    fill = "Status"
  ) +
  scale_fill_manual(values = c("ERCC2-MUT" = "skyblue", "ERCC2-WT" = "lightpink"))

# Summary statistics for each group
summary(total_mutations_sbs %>% filter(Status == "ERCC2-MUT"))
summary(total_mutations_sbs %>% filter(Status == "ERCC2-WT"))

# Calculate total INDEL mutations per sample
total_mutations_indels <- colSums(indels_count) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample")

# Add a Status column that will show the ERCC2 status (MUT or WT)
total_mutations_indels <- total_mutations_indels %>%
  mutate(Status = ifelse(Sample %in% ercc2_helicase_samples, "ERCC2-MUT", "ERCC2-WT"))
colnames(total_mutations_indels) <- c("Sample", "TotalMutations", "Status")

# PLot the distribution 
ggplot(total_mutations_indels, aes(x = Status, y = TotalMutations, fill = Status)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Total Mutations Across Samples (INDEL) - TCGA-BLCA WES",
    x = "Sample Status",
    y = "Total Number of Mutations",
    fill = "Status"
  ) +
  scale_fill_manual(values = c("ERCC2-MUT" = "skyblue", "ERCC2-WT" = "lightpink"))

# Summary statistics 
summary(total_mutations_indels %>% filter(Status == "ERCC2-MUT"))
summary(total_mutations_indels %>% filter(Status == "ERCC2-WT"))

# Filter maf for MSK IMPACT genes
# Separate ERCC2 MUT and WT samples
ercc2_mut_samples <- maf %>% 
  filter(Tumor_Sample_Barcode %in% ercc2_helicase_samples) %>% 
  pull(Tumor_Sample_Barcode) %>% 
  unique() # 34 samples

ercc2_wt_samples <- maf %>% 
  filter(!Tumor_Sample_Barcode %in% ercc2_helicase_samples) %>% 
  pull(Tumor_Sample_Barcode) %>% 
  unique() # 380 samples

# Separate mutation types: SBS and indels
maf_sbs <- maf %>% filter(Variant_Type == "SNP")
maf_indel <- maf %>% filter(Variant_Type %in% c("INS", "DEL"))

# Mutation ranges for ERCC2 MUT and WT determiend from the distributions 
mut_ranges <- list(
  SBS = list(MUT = c(200, 900), WT = c(80, 550)),  # SBS ranges
  Indel = list(MUT = c(3, 21), WT = c(2, 10))   # Indel ranges
)

# Function to simulate WES samples
set.seed(123)

simulate_wes_samples <- function(sample_count, sbs_range, indel_range, sbs_data, indel_data, sample_prefix) {
  simulated_samples <- list()
  
  for (i in 1:sample_count) {
    # Sample the number of SBS mutations for this sample
    total_sbs <- sample(seq(sbs_range[1], sbs_range[2]), 1)
    
    # Sample the number of INDEL mutations for this sample
    total_indel <- sample(seq(indel_range[1], indel_range[2]), 1)
    
    # Sample mutations from the respective datasets
    sampled_sbs <- sbs_data %>%
      sample_n(total_sbs, replace = TRUE) %>% 
      mutate(Tumor_Sample_Barcode = paste0(sample_prefix, "_Simulated_Sample_", i))
    
    sampled_indel <- indel_data %>%
      sample_n(total_indel, replace = TRUE) %>% 
      mutate(Tumor_Sample_Barcode = paste0(sample_prefix, "_Simulated_Sample_", i))
    
    # Combine SBS and Indel mutations for this sample
    simulated_sample <- bind_rows(sampled_sbs, sampled_indel)
    simulated_samples[[i]] <- simulated_sample
  }
  
  # Combine all simulated samples
  simulated_samples <- bind_rows(simulated_samples)
  return(simulated_samples)
}

# Simulate 250 ERCC2 MUT samples with both SBS and Indels
simulated_mut <- simulate_wes_samples(
  sample_count = 250,
  sbs_range = mut_ranges$SBS$MUT,
  indel_range = mut_ranges$Indel$MUT,
  sbs_data = maf_sbs %>% filter(Tumor_Sample_Barcode %in% ercc2_mut_samples),
  indel_data = maf_indel %>% filter(Tumor_Sample_Barcode %in% ercc2_mut_samples),
  sample_prefix = "MUT"
)

# Simulate 250 ERCC2 WT samples with both SBS and Indels
simulated_wt <- simulate_wes_samples(
  sample_count = 250,
  sbs_range = mut_ranges$SBS$WT,
  indel_range = mut_ranges$Indel$WT,
  sbs_data = maf_sbs %>% filter(Tumor_Sample_Barcode %in% ercc2_wt_samples),
  indel_data = maf_indel %>% filter(Tumor_Sample_Barcode %in% ercc2_wt_samples),
  sample_prefix = "WT"
)

# Combine simulated data
simulated_wes_data <- bind_rows(simulated_mut, simulated_wt)

# Filter for only MSK-IMPACT505 genes
simulated_wes_msk505 <- simulated_wes_data %>% filter(Hugo_Symbol %in% msk505)
