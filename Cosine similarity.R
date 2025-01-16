# Cosine similarity with SBS5 and ID8 signature

cosine_similarity_sbs <- cos_sim_matrix(simulated_wes_sbs_mut_matrix, as.matrix(cosmic_signature_names[, "SBS5"]))

cosine_similarity_indel <- cos_sim_matrix(simulated_wes_indel_mut_matrix, as.matrix(cosmic_signature_indel[, "ID8"]))