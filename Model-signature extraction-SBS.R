# Filter sbs mutational matrix for > 50 mutations
sbs_mutation_matrix_50_sigma <- sbs_mutation_matrix[, colSums(sbs_mutation_matrix) >= 50] 

# Define MSI-positive and POLE-positive samples
pole_msi_samples <- c("TCGA-DK-A6AW-01A-11D-A30E-08", "TCGA-XF-A8HG-01A-11D-A364-08", "TCGA-ZF-AA4W-01A-12D-A38G-08")
is_pole_msi <- colnames(sbs_mutation_matrix_50_sigma) %in% pole_msi_samples

# Remove the samples from the mutational count matrix
mut_matrix_sigma_filtered <- sbs_mutation_matrix_50_sigma[, !is_pole_msi] # 385 samples left

# NLLS decomposition of the mutational matrix and the cosmic catalog signatures (subsetted only for the signatures found in bladder cancer) - fitting 
nnls_fitting_1 <- fit_to_signatures_strict(mut_matrix_sigma_filtered, cosmic_signature_names[, bladder_tissue_sbs], method = "best_subset", max_delta = 0.004)

# Plot of absolute contribution absolute and relative
plot_contribution(nnls_fitting_1$fit_res$contribution[, 1:15], mode = 'absolute', signatures = cosmic_signature_names[, bladder_tissue_sbs], coord_flip = TRUE)
plot_contribution(nnls_fitting_1$fit_res$contribution[, 1:15], mode = 'relative', coord_flip = TRUE)

# Check for samples with contribution < 0.05
samples_with_low_contribution <- apply(nnls_fitting_1$fit_res$contribution, 2, function(col) any(col <= 0.05))
samples_with_low_contr_index <- which(samples_with_low_contribution) 

# Try and do the fitting to the same set of signatures except without SBS29 and SBS40, most of the samples don't have those signatures
bladder_tissue_sbs_subset <- c("SBS1", "SBS2", "SBS5", "SBS13")
nnls_fitting_2 <- fit_to_signatures_strict(mut_matrix_sigma_filtered, cosmic_signature_names[, bladder_tissue_sbs_subset], method = "best_subset", max_delta = 0.004)

plot_contribution(nnls_fitting_2$fit_res$contribution[, 1:15], mode = 'absolute', signatures = cosmic_signature_names[, bladder_tissue_sbs_subset], coord_flip = TRUE)
plot_contribution(nnls_fitting_2$fit_res$contribution[, 1:15], mode = 'relative', coord_flip = TRUE)

# Check for samples with low contribution from signatures
samples_with_low_contribution_2 <- apply(nnls_fitting_2$fit_res$contribution, 2, function(col) any(col <= 0.05))
samples_with_low_contr_index_2 <- which(samples_with_low_contribution_2) 

# Plot original vs reconstructed for the final profile between the mutational matrix 
plot_original_vs_reconstructed(mut_matrix_sigma_filtered[, 1:15], nnls_fitting_2$fit_res$reconstructed[, 1:15])