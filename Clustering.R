# Clustering - separate for SBS and INDEL
sbs_contr_mat <-nnls_fitting_2$fit_res$contribution 
indel_contr_mat <- nnls_fitting_indel_2$fit_res$contribution

# Cluster SBS
# Normalize columns so each column sums to 1
sbs_contr_normalized <- apply(sbs_contr_mat, 2, function(x) x / sum(x))

# Transpose the matrix
transposed_sbs_matrix <- t(sbs_contr_normalized)

# Compute distance matrix (euclidean)
dist_matrix_sbs <- dist(transposed_sbs_matrix)  

# Hierarchical clustering
hc_sbs <- hclust(dist_matrix_sbs, method = "average")
plot(hc_sbs, labels = FALSE, main = "Hierarchical Clustering Dendrogram")

# Determine the optimal number of clusters using the elbow method
wcss_sbs <- sapply(1:10, function(k) {
  cluster_assignments <- cutree(hc_sbs, k = k) # Cluster assignment
  sum(sapply(unique(cluster_assignments), function(cluster) {
    cluster_points <- clustering_data[cluster_assignments == cluster, , drop = FALSE]
    # Handle single-row clusters
    cluster_center <- if (nrow(cluster_points) > 1) {
      colMeans(cluster_points)
    } else {
      cluster_points
    }
    sum(rowSums((cluster_points - cluster_center)^2))
  }))
})

# Plot elbow plot
plot(1:10, wcss_sbs, type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of Clusters", ylab = "Within-Cluster Sum of Squares (WCSS)",
     main = "Elbow Method for Optimal Clusters")

# Assign clusters
k <- 2  
clusters_sbs <- cutree(hc_sbs, k)

# Add cluster information
sbs_mat_with_clusters <- as.data.frame(transposed_sbs_matrix)
sbs_mat_with_clusters$cluster <- clusters_sbs

# Make short id to match the short id in helicsase domain to see to which cluster the samples that are in the helicase domain belong
helicase_s_sbs <- sbs_mat_with_clusters[
  grepl(paste(helicase_domain_ercc2_samples, collapse = "|"), rownames(sbs_mat_with_clusters)), 
] # 34 samples, 29 are in cluster 2

# Calculate average contribution per cluster
average_sbs_mat <- aggregate(. ~ cluster, data = sbs_mat_with_clusters, FUN = mean)

# Normalize rows so each cluster's contributions sum to 1
average_sbs_normalized_mat <- t(apply(average_sbs_mat[, -1], 1, function(x) x / sum(x))) # SBS5 in cluster 1 is 0 and in cluster 2 is 0.45 

melted_data_sbs <- melt(sbs_contr_normalized)
colnames(melted_data_sbs) <- c("Signature", "Sample", "Contribution")

# Add cluster assignments to the melted data
melted_data_sbs$Cluster <- clusters_sbs[melted_data_sbs$Sample]

# Plot of relative contribution of each signature to each sample with clusters (continuous)
ggplot(melted_data_sbs, aes(x = Sample, y = Contribution, fill = Signature)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~Cluster, scales = "free_x", space = "free_x") +  # Group samples by clusters
  theme_minimal() +
  labs(x = "Samples", y = "Contribution", fill = "Signature") +
  scale_fill_brewer(palette = "Set3") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.major.x = element_blank()
  )

# Make box plot of the contribution of each signature in each cluster
ggplot(melted_data_sbs, aes(x = Signature, y = Contribution, fill = factor(Cluster))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.8)) +
  theme_minimal() +
  labs(
    x = "SBS Signature",
    y = "Contribution",
    fill = "Cluster",
    title = "SBS Signature Contributions by Cluster"
  ) +
  scale_fill_brewer(palette = "Set3") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

# Cluster INDELS
indel_contr_normalized <- apply(indel_contr_mat, 2, function(x) x / sum(x))

# Transpose the matrix
transposed_indel_matrix <- t(indel_contr_normalized)

# Compute distance matrix (Euclidean)
dist_matrix_indel <- dist(transposed_indel_matrix)  

# Hierarchical clustering
hc_indel <- hclust(dist_matrix_indel, method = "average")

plot(hc_indel, labels = FALSE, main = "Hierarchical Clustering Dendrogram")

# Determine the optimal number of clusters using the elbow method
wcss_indel <- sapply(1:10, function(k) {
  cluster_assignments <- cutree(hc_indel, k = k) # Cluster assignment
  sum(sapply(unique(cluster_assignments), function(cluster) {
    cluster_points <- clustering_data[cluster_assignments == cluster, , drop = FALSE]
    # Handle single-row clusters
    cluster_center <- if (nrow(cluster_points) > 1) {
      colMeans(cluster_points)
    } else {
      cluster_points
    }
    sum(rowSums((cluster_points - cluster_center)^2))
  }))
})

# Plot elbow plot
plot(1:10, wcss_indel, type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of Clusters", ylab = "Within-Cluster Sum of Squares (WCSS)",
     main = "Elbow Method for Optimal Clusters")

# Assign clusters
k <- 5   
clusters_indel <- cutree(hc_indel, k)

# Add cluster information
indel_mat_with_clusters <- as.data.frame(transposed_indel_matrix)
indel_mat_with_clusters$cluster <- clusters_indel

helicase_s_indel <- indel_mat_with_clusters[
  grepl(paste(helicase_domain_ercc2_samples, collapse = "|"), rownames(indel_mat_with_clusters)), 
] # 33 samples, 2 are in cluster 2, 1 in cluster 4, the other 30 in cluster 1

# Calculate average contribution per cluster
average_indel_mat <- aggregate(. ~ cluster, data = indel_mat_with_clusters, FUN = mean)

# Normalize rows so each cluster's contributions sum to 1
average_indel_normalized_mat <- t(apply(average_indel_mat[, -1], 1, function(x) x / sum(x))) # ID8 in all clusters, if i do it with all id signatures without removing, i have one cluster without id8

melted_data_indel <- melt(indel_contr_normalized)
colnames(melted_data_indel) <- c("Signature", "Sample", "Contribution")

# Add cluster assignments to the melted data
melted_data_indel$Cluster <- clusters_indel[melted_data_indel$Sample]

# Plot of relative contribution of each signature to each sample with clusters (continuous)
ggplot(melted_data_indel, aes(x = Sample, y = Contribution, fill = Signature)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~Cluster, scales = "free_x", space = "free_x") +  # Group samples by clusters
  theme_minimal() +
  labs(x = "Samples", y = "Contribution", fill = "Signature") +
  scale_fill_brewer(palette = "Set3") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.major.x = element_blank()
  )

# box plot
ggplot(melted_data_indel, aes(x = Signature, y = Contribution, fill = factor(Cluster))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.8)) +
  theme_minimal() +
  labs(
    x = "INDEL Signature",
    y = "Contribution",
    fill = "Cluster",
    title = "INDEL Signature Contributions by Cluster"
  ) +
  scale_fill_brewer(palette = "Set3") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )
