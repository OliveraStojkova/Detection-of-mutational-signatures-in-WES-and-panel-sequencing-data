# Retrieve INDELs
indels_grl <- get_mut_type(maf_grl, type = "indel")

# Set the genome (specify)
GenomeInfoDb::genome(indels_grl) = 'hg38'

# Indel context - excludes sampels thathave no indels = less samples
indels_seq_context <- get_indel_context(indels_grl, ref_genome = ref_genome)

# Indel count matrix
indels_count <- count_indel_contexts(indels_seq_context)  # This matrix will be used later to extract signatures

dim(indels_count) # 351 samples

# Separate based on the presence of ERCC2 mutation
indel_ercc2_matrix <- indels_count[, colnames(indels_count) %in% ercc2_mutated_samples, drop = FALSE]
indel_no_ercc2_matrix <- indels_count[, !(colnames(indels_count) %in% ercc2_mutated_samples), drop = FALSE]

# Plot indel context with counts (samples 1:3)
plot_indel_contexts(indels_count[, 1:3], condensed = TRUE)
plot_indel_contexts(indel_ercc2_matrix[, 1:3], condensed = TRUE)
plot_indel_contexts(indel_no_ercc2_matrix[, 1:3], condensed = TRUE)

# Plot ERCC2 helicase domain vs not mutated in ERCC2
indel_non_helicase <- indels_count
colnames(indel_non_helicase) <- substr(colnames(indels_count), 1, 12)
indel_non_helicase <- indel_non_helicase[, !(colnames(indel_non_helicase) %in% helicase_domain_ercc2_samples), drop = FALSE]

plot_indel_contexts(indel_ercc2_matrix[, is_helicase_indel_counts][, 5:7], condensed = TRUE)
plot_indel_contexts(indel_non_helicase[, 5:7], condensed = TRUE)

# Plot only the main contexts, without taking the number of repeat units or microhomology length into account (samples 1-4)
plot_main_indel_contexts(indel_ercc2_matrix[, is_helicase_indel_counts][, 5:8])
plot_main_indel_contexts(indel_non_helicase[, 5:8])