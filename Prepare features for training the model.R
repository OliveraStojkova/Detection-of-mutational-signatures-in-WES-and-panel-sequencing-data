# Prepare features for MVA

# Define a function to calculate features for a single sample
calculate_features_for_sample <- function(sample_idx, 
                                          sbs5_cosine_similarity_matrix, 
                                          id8_cosine_similarity_matrix, 
                                          sbs_likelihood_matrix, 
                                          indel_likelihood_matrix,
                                          sbs_exposures, 
                                          id_exposures, 
                                          likelihood_ratio_vector_sbs,
                                          likelihood_ratio_vector_id,
                                          sbs_counts, 
                                          indel_counts) {
  
  # Check if sample exists in SBS counts, if not add NA
  if (sample_idx %in% colnames(sbs_counts)) {
    cosine_sbs5 <- sbs5_cosine_similarity_matrix[sample_idx, ]  
    likelihood_sbs5_positive <- sbs_likelihood_matrix[2, sample_idx]  
    likelihood_sbs5_negative <- sbs_likelihood_matrix[1, sample_idx] 
    sbs5_exposure <- sbs_exposures[sample_idx, "SBS5"] 
    likelihood_ratio_sbs <- likelihood_ratio_vector_sbs[sample_idx, ]  
    total_sbs_count <- sum(sbs_counts[, sample_idx]) 
  } else {
    cosine_sbs5 <- NA
    likelihood_sbs5_positive <- NA
    likelihood_sbs5_negative <- NA
    sbs5_exposure <- NA
    likelihood_ratio_sbs <- rep(NA, ncol(likelihood_ratio_vector_sbs))
    total_sbs_count <- NA
  }
  
  # Check if sample exists in Indel counts, if not add NA
  if (sample_idx %in% colnames(indel_counts)) {
    cosine_id8 <- id8_cosine_similarity_matrix[sample_idx, ]  
    likelihood_indel_c1 <- indel_likelihood_matrix[1, sample_idx]  
    likelihood_indel_c2 <- indel_likelihood_matrix[2, sample_idx]
    likelihood_indel_c3 <- indel_likelihood_matrix[3, sample_idx]
    likelihood_indel_c4 <- indel_likelihood_matrix[4, sample_idx]
    likelihood_indel_c5 <- indel_likelihood_matrix[5, sample_idx]
    
    insertion_rows <- grep("insertion", rownames(indel_counts))
    total_ins_count <- sum(indel_counts[insertion_rows, sample_idx])
    deletion_rows <- grep("deletion", rownames(indel_counts))
    total_del_count <- sum(indel_counts[deletion_rows, sample_idx])  
    
    id8_exposure <- id_exposures[sample_idx, "ID8"] 
    likelihood_ratio_id <- likelihood_ratio_vector_id[sample_idx, ]
  } else { 
    cosine_id8 <- NA
    likelihood_indel_c1 <- NA 
    likelihood_indel_c2 <- NA
    likelihood_indel_c3 <- NA
    likelihood_indel_c4 <- NA
    likelihood_indel_c5 <- NA
    total_ins_count <- NA
    total_del_count <- NA
    id8_exposure <- NA
    likelihood_ratio_id <- rep(NA, ncol(likelihood_ratio_vector_id))
  }
  
  # Combine features into a data frame for the sample (rows = number of samples)
  features <- data.frame(
    cosine_sbs5 = cosine_sbs5,  
    cosine_id8 = cosine_id8,
    likelihood_sbs5_positive = likelihood_sbs5_positive,
    likelihood_sbs5_negative = likelihood_sbs5_negative,
    likelihood_indel_c1 = likelihood_indel_c1,
    likelihood_indel_c2 = likelihood_indel_c2,
    likelihood_indel_c3 = likelihood_indel_c3,
    likelihood_indel_c4 = likelihood_indel_c4,
    likelihood_indel_c5 = likelihood_indel_c5,
    sbs5_exposure = sbs5_exposure,
    id8_exposure = id8_exposure,
    likelihood_ratio_sbs = likelihood_ratio_sbs,
    likelihood_ratio_id = likelihood_ratio_id,
    total_sbs_count = total_sbs_count,
    total_ins_count = total_ins_count,
    total_del_count = total_del_count
  )
  
  return(features)
}


# Loop over all samples in artificially simulated WES samples
all_samples <- union(colnames(simulated_wes_indel_mut_matrix), colnames(simulated_wes_sbs_mut_matrix))  # Include all unique samples - 500 smaples total
all_features <- list()  

for (sample_idx in all_samples) {
  sample_features <- calculate_features_for_sample(
    sample_idx = sample_idx, 
    sbs5_cosine_similarity_matrix = cosine_similarity_sbs,
    id8_cosine_similarity_matrix = cosine_similarity_indel,
    sbs_likelihood_matrix = likelihood_posterior_sbs$likelihoods,
    indel_likelihood_matrix = likelihood_posterior_indel$likelihoods,
    sbs_exposures = iterative_exposures_sbs$exposures,
    id_exposures = iterative_exposures_indel$exposures,
    likelihood_ratio_vector_sbs = nnls_likelihood_ratio_sbs,
    likelihood_ratio_vector_id = nnls_likelihood_ratio_indel, 
    sbs_counts = simulated_wes_sbs_mut_matrix,
    indel_counts = simulated_wes_indel_mut_matrix
  )
  all_features[[sample_idx]] <- sample_features
}

# Combine all features into one data frame
features_df <- as.data.frame(do.call(rbind, all_features))

# Add sample names
features_df$sample_names <- all_samples

# Define lbeles (ERCC2-MUT vs ERCC2-WT) for training the model
features_df$label <- ifelse(
  grepl("^MUT_", features_df$sample_names), "ercc2-MUT",
  ifelse(grepl("^WT_", features_df$sample_names), "ercc2-WT", NA)
)

# Change the label to binary 
features_df$label <- ifelse(features_df$label == "ercc2-MUT", 1, 0)

# Change any NaN values to NA
features_df[] <- lapply(features_df, function(x) ifelse(is.nan(x), NA, x)) # A few NA value in cosine similarity for id