# Pseudocount, because fit_to_signatures gives error otherwise
indel <- indels_count+0.001 

# Remove the samples from the mutational count matrix
is_pole_msi_indel <- colnames(indel) %in% pole_msi_samples
indel <- indel[, !is_pole_msi_indel] # 348 samples remaining
indels_count <- indels_count[, !is_pole_msi_indel]

# NNLS fit
cosmic_catalog_subset_indel <- cosmic_signature_indel[, bladder_tissue_indel] 
nnls_fitting_indel <- fit_to_signatures_strict(indel, cosmic_catalog_subset_indel, method = "best_subset", max_delta = 0.004) 

# Plot of absolute contribution absolute and relative
plot_contribution(nnls_fitting_indel$fit_res$contribution[, 1:15], mode = 'absolute', signatures = cosmic_catalog_subset, coord_flip = TRUE)
plot_contribution(nnls_fitting_indel$fit_res$contribution[, 1:15], mode = 'relative', coord_flip = TRUE)

# Check for samples with contribution < 0.05
samples_with_low_contribution_indel <- apply(nnls_fitting_indel$fit_res$contribution, 2, function(col) any(col <= 0.05))
samples_with_low_contr_index_indel <- which(samples_with_low_contribution_indel) # all samples

# Original vs reconstructed
plot_original_vs_reconstructed(indels_count, nnls_fitting_indel$fit_res$reconstructed)

# Remove id1 and id2 and redo fitting, these 2 only had contribution in 1 sample
bladder_tissue_indel_s <- c("ID4", "ID5","ID8","ID10")
cosmic_catalog_subset_indel_s <- cosmic_signature_indel[, bladder_tissue_indel_s] 

nnls_fitting_indel_2 <- fit_to_signatures_strict(indel, cosmic_catalog_subset_indel_s, method = "best_subset", max_delta = 0.004) 

# Plot of absolute contribution absolute and relative
plot_contribution(nnls_fitting_indel_2$fit_res$contribution[, 1:15], mode = 'absolute', signatures = cosmic_catalog_subset_indel_s, coord_flip = TRUE)
plot_contribution(nnls_fitting_indel_2$fit_res$contribution[, 1:15], mode = 'relative', coord_flip = TRUE)

# Check for samples with contribution < 0.05
samples_with_low_contribution_indel_2 <- apply(nnls_fitting_indel_2$fit_res$contribution, 2, function(col) any(col <= 0.05))
samples_with_low_contr_index_indel_2 <- which(samples_with_low_contribution_indel_2) 

# Plot original vs. reconstructed
plot_original_vs_reconstructed(indels_count[, 1:15], nnls_fitting_indel_2$fit_res$reconstructed[, 1:15])
