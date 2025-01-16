### Implement the model (GBC-3) on the MSK-IMPACT505 dataset

# Get all necessary features
# Cosine similarity
cos_sim_msk_sbs5 <- cos_sim_matrix(msk505_sbs_mut_matrix, as.matrix(cosmic_signature_names[, "SBS5"]))
cos_sim_msk_id8 <- cos_sim_matrix(msk505_indel_mut_matrix, as.matrix(cosmic_signature_indel[, "ID8"]))

# Likelihood sbs and indel
msk505_sbs_likelihood <- calculate_likelihood_and_posterior(S_t_matrix = msk505_sbs_mut_matrix, D_i = D_i_sbs, P_D_i = P_D_i_sbs)
colnames(msk505_sbs_likelihood$likelihoods) <- colnames(msk505_sbs_mut_matrix)

msk505_indel_likelihood <- calculate_likelihood_and_posterior(S_t_matrix = msk505_indel_mut_matrix, D_i = D_i_indel, P_D_i = P_D_i_indel)
colnames(msk505_indel_likelihood$likelihoods) <- colnames(msk505_indel_mut_matrix)

# Iterative exposures - nnls
msk505_iterative_exposures_sbs <- iterative_nnls(cosmic_signature_names[, bladder_tissue_sbs], msk505_sbs_mut_matrix, max_signatures = 5)
msk505_iterative_exposures_indel <- iterative_nnls(cosmic_signature_indel[, bladder_tissue_indel], msk505_indel_mut_matrix, max_signatures = 5)

# Likelihood ratio + alternative solutions
alternative_solutions_sbs_msk505 <- get_alternative_solutions(msk505_sbs_mut_matrix, cosmic_signature_names[, bladder_tissue_sbs])
alternative_solutions_indel_msk505 <- get_alternative_solutions(msk505_indel_mut_matrix, cosmic_signature_indel[, bladder_tissue_indel])

nnls_likelihood_ratio_sbs_msk505 <- calculate_likelihood_nnls(cosmic_signature_names[, bladder_tissue_sbs], msk505_sbs_mut_matrix, alternative_solutions_sbs_msk505)
nnls_likelihood_ratio_indel_msk505 <- calculate_likelihood_nnls(cosmic_signature_indel[, bladder_tissue_indel], msk505_indel_mut_matrix, alternative_solutions_indel_msk505)

# Make features dataframe
all_samples_msk505 <- union(colnames(msk505_sbs_mut_matrix), colnames(msk505_indel_mut_matrix)) # 912 samples total
all_features_msk <- list()  

for (sample_idx in all_samples_msk505) {
  sample_features <- calculate_features_for_sample(
    sample_idx = sample_idx, 
    sbs5_cosine_similarity_matrix = cos_sim_msk_sbs5,
    id8_cosine_similarity_matrix = cos_sim_msk_id8,
    sbs_likelihood_matrix = msk505_sbs_likelihood$likelihoods,
    indel_likelihood_matrix = msk505_indel_likelihood$likelihoods,
    sbs_exposures = msk505_iterative_exposures_sbs$exposures,
    id_exposures = msk505_iterative_exposures_indel$exposures,
    likelihood_ratio_vector_sbs = nnls_likelihood_ratio_sbs_msk505,
    likelihood_ratio_vector_id = nnls_likelihood_ratio_indel_msk505, 
    sbs_counts = msk505_sbs_mut_matrix,
    indel_counts = msk505_indel_mut_matrix
  )
  all_features_msk[[sample_idx]] <- sample_features
}

# Combine all features into one data frame
features_df_msk505 <- as.data.frame(do.call(rbind, all_features_msk)) # 912 samples

# Predict scores and labels with gbc-3 model
msk_scores_gbc3 <- predict(
  gbm_model_3, 
  newdata = features_df_msk505, 
  n.trees = best_trees, 
  type = "response"
)

predicted_labels_features_df_msk505 <- ifelse(msk_scores_gbc3 >= optimal_cutpoint, 1, 0)
# Gives the number of samples that have a label 1 (ERCC2-MUT) and 0 (ERCC2-WT)
table(predicted_labels_features_df_msk505)

# Predict scores with gbbc-2
msk_scores_gbc2 <- predict(
  gbm_model_2, 
  newdata = features_df_msk505, 
  n.trees = best_trees, 
  type = "response"
)

# Plot GBC2 vs GBC3
scatter_data_msk505 <- data.frame(
  GBC2 = msk_scores_gbc2,
  GBC3 = msk_scores_gbc3
)

ggplot(scatter_data_msk505, aes(x = GBC2, y = GBC3)) +
  geom_point(alpha = 0.8, size = 1, color = "black") +
  geom_vline(xintercept = optimal_cutpoint, linetype = "dashed", color = "black", size =0.5) + 
  geom_hline(yintercept = optimal_cutpoint, linetype = "dashed", color = "black", size = 0.5) + 
  theme_minimal() +
  labs(
    x = "GBC-2 Scores",
    y = "GBC-3 Scores",
    title = "GBC-2 vs GBC-3 Scores for MSK-IMPACT505"
  )
