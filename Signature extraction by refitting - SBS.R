# Get known COSMIC signatures
cosmic_signatures = get_known_signatures(
  source = 'COSMIC_v3.2', genome = 'GRCh38') 

# Subset signatures only in bladder tissue
bladder_tissue_sbs <- c("SBS1", "SBS2", "SBS5", "SBS8", "SBS13", "SBS29", "SBS40")

# Fitting - normalized (>=50 mutations) +subset
sbs_fit_normalized <- fit_to_signatures(normalized_sbs_matrix, cosmic_signatures[, bladder_tissue_sbs])

# Visualization for the first 15 samples - Normalized matrix + subset
plot_contribution(sbs_fit_normalized$contribution[,1:15],
                  coord_flip = TRUE,
                  mode = "relative"
)

# Visualization for the first 15 samples - ERCC2 vs. other (Normalized + subset)
sbs_fit_normalized_ercc2 <- sbs_fit_normalized$contribution[, colnames(sbs_fit_normalized$contribution) %in% ercc2_mutated_samples, drop = FALSE]
sbs_fit_normalized_no_ercc2 <- sbs_fit_normalized$contribution[, !(colnames(sbs_fit_normalized$contribution) %in% ercc2_mutated_samples), drop = FALSE]

# ERCC2
plot_contribution(sbs_fit_normalized_ercc2[,1:15],
                  coord_flip = TRUE,
                  mode = "relative", # Mode can also be absolute
) 

# NO ERCC2 mutation
plot_contribution(sbs_fit_normalized_no_ercc2[,1:15],
                  coord_flip = TRUE,
                  mode = "relative"
)

# For those seen in ERCC2, whihc ones are in the helicase domain
sample_ids <- colnames(sbs_fit_normalized_ercc2)
short_ids <- substr(sample_ids, 1, 12)

is_helicase <- short_ids %in% helicase_domain_ercc2_samples 
new_colnames <- paste0(short_ids)
colnames(sbs_fit_normalized_ercc2) <- new_colnames

# Split data into helicase and non-helicase
helicase_data <- sbs_fit_normalized_ercc2[, is_helicase]
non_helicase_data <- sbs_fit_normalized_ercc2[, !is_helicase]

# Plot ERCC2 helicase samples, can be compared to the one without ERCC2 mutation at all, or just the ERCC2 mutated samples that are not in helicase domain
plot_contribution(helicase_data[, 1:15],
                  coord_flip = TRUE,
                  mode = "relative"
)

plot_contribution(non_helicase_data,
                  coord_flip = TRUE,
                  mode = "relative"
)

# Original vs reconstructed plot - normalized and subset
plot_original_vs_reconstructed(normalized_sbs_matrix[,1:20], sbs_fit_normalized$reconstructed[,1:20],
                               y_intercept = 0.95)