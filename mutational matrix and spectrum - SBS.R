# Retrieve SBSs
sbs_grl <- get_mut_type(maf_grl, type = "snv")

# Set the genome (specify)
GenomeInfoDb::genome(sbs_grl) = 'hg38'

# Get mutation matrix for SBS
sbs_mutation_matrix <- mut_matrix(sbs_grl, ref_genome = ref_genome) # 414 samples

# Separate based on the presence of ERCC2
sbs_ercc2_matrix <- sbs_mutation_matrix[, ercc2_mutated_samples, drop = FALSE]
sbs_no_ercc2_matrix <- sbs_mutation_matrix[, !(colnames(sbs_mutation_matrix) %in% ercc2_mutated_samples), drop = FALSE]

# Plot the 96 profile - for 3 samples
plot_96_profile(sbs_mutation_matrix[, 1:3]))
plot_96_profile(sbs_ercc2_matrix[, 1:3]))
plot_96_profile(sbs_no_ercc2_matrix[, 1:3])

# ERCC2 mutated in helicase, vs not mutated ERCC2
sbs_mutation_matrix1 <- sbs_mutation_matrix
colnames(sbs_mutation_matrix1) <- substr(colnames(sbs_mutation_matrix), 1, 12) 

# Contains samples not in ERCC2 helicase, so no mutated ERCC2
sbs_mutation_matrix1 <- sbs_mutation_matrix1[, !(colnames(sbs_mutation_matrix1) %in% helicase_domain_ercc2_samples), drop = FALSE]

# Plot and compare - ERCC2 helicase mutation vs. no mutation in ERCC2
plot_96_profile(sbs_ercc2_matrix[, is_helicase][, 1:2] )
plot_96_profile(sbs_mutation_matrix1[, 1:2])

