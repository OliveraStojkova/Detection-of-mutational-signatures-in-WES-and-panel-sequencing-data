# Get known indel signatures
indel_cosmic_signatures = get_known_signatures(muttype = 'indel',
                                               source = 'COSMIC_v3.2') 

# Subset the signatures present in bladder and do the fitting
bladder_tissue_indel <- c("ID1", "ID2", "ID3", "ID4", "ID5", "ID8", "ID9", "ID10")

indel_fit_res_subset <- fit_to_signatures(indels_count, indel_cosmic_signatures[, bladder_tissue_indel])

# Visualization for the first 15 samples
plot_contribution(indel_fit_res$contribution[,1:15],
                  coord_flip = TRUE,
                  mode = "relative"
)

# Visualization for the first 15 samples - Subset
plot_contribution(indel_fit_res_subset$contribution[,1:15],
                  coord_flip = TRUE,
                  mode = "relative"
)

# Visualization ERCC2 vs other
indel_fit_ercc2 <- indel_fit_res_subset$contribution[, colnames(indel_fit_res_subset$contribution) %in% ercc2_mutated_samples, drop = FALSE]
indel_fit_no_ercc2 <- indel_fit_res_subset$contribution[, !(colnames(indel_fit_res_subset$contribution) %in% ercc2_mutated_samples), drop = FALSE]

# ERCC2
plot_contribution(indel_fit_ercc2[,1:15],
                  coord_flip = TRUE,
                  mode = "relative"
)

# NO ERCC2 mutation
plot_contribution(indel_fit_no_ercc2[,1:15],
                  coord_flip = TRUE,
                  mode = "relative"
)

# For ERCC2 - helicase vs. non-helicase
sample_ids_indel <- colnames(indel_fit_ercc2)
short_ids_indel <- substr(sample_ids_indel, 1, 12)
colnames(indel_fit_ercc2) <- short_ids_indel

is_helicase_indel <- colnames(indel_fit_ercc2) %in% helicase_domain_ercc2_samples

plot_contribution(indel_fit_ercc2[, is_helicase_indel][, 1:15],
                  coord_flip = TRUE,
                  mode = "absolute"
)

plot_contribution(indel_fit_ercc2[, !is_helicase_indel],
                  coord_flip = TRUE,
                  mode = "relative"
)

# Original vs reconstructed plot 
plot_original_vs_reconstructed(indels_count[,1:20], indel_fit_res_subset$reconstructed[, 1:20], 
                               y_intercept = 0.95)

