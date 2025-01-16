# Likelihood and posterior calculation: The likelihood is the probability of observing the mutational spectrum St, given a specific cluster Di. 
# It is calculated directly using the cluster distributions and is often represented as a product of probabilities (or the sum of logs in the log-likelihood).
# The posterior probability, i.e the probability of a cluster Di being the origin of St, after considering the prior probabilities and normalizing by the total 
# probability of St across all clusters. The function below does both. 
# Input data - mutational spectrum matrix (S_t_matrix); cluster distributions (D_i) and cluster proportions (number of samples in current cluster/ total number of samples)(P_D_i)

# SBS clusters and proportions (number of samples in current cluster/total number of samples) - calcuated based on the clusters determined 
D_i_sbs <- list(
  cluster1 = average_sbs_normalized_mat[1, ],
  cluster2 = average_sbs_normalized_mat[2, ]
)

c1_sbs <- melted_data_sbs%>%filter(Cluster == 1)
c1_sbs <- unique(c1_sbs$Sample)
n_c1_sbs <- length(c1_sbs) # 88

c2_sbs <- melted_data_sbs%>%filter(Cluster == 2)
c2_sbs <- unique(c2_sbs$Sample)
n_c2_sbs <- length(c2_sbs) # 297

n_sbs_total <- n_c1_sbs + n_c2_sbs # 385

P_D_1_sbs <- n_c1_sbs/n_sbs_total
P_D_2_sbs <- n_c2_sbs/n_sbs_total

P_D_i_sbs <- c(P_D_1_sbs, P_D_2_sbs)

# INDEL clusters and proportions (number of samples in current cluster/total number of samples)
D_i_indel <- list(
  cluster1 = average_indel_normalized_mat[1, ],
  cluster2 = average_indel_normalized_mat[2, ],
  cluster3 = average_indel_normalized_mat[3, ],
  cluster4 = average_indel_normalized_mat[4, ],
  cluster5 = average_indel_normalized_mat[5, ]
)

c1_indel <- melted_data_indel%>%filter(Cluster == 1)
c1_indel <- unique(c1_indel$Sample)
n_c1_indel <- length(c1_indel) # 208

c2_indel <- melted_data_indel%>%filter(Cluster == 2)
c2_indel <- unique(c2_indel$Sample)
n_c2_indel <- length(c2_indel) # 36

c3_indel <- melted_data_indel%>%filter(Cluster == 3)
c3_indel <- unique(c3_indel$Sample)
n_c3_indel <- length(c3_indel) # 33

c4_indel <- melted_data_indel%>%filter(Cluster == 4)
c4_indel <- unique(c4_indel$Sample)
n_c4_indel <- length(c4_indel) # 37

c5_indel <- melted_data_indel%>%filter(Cluster == 5)
c5_indel <- unique(c5_indel$Sample)
n_c5_indel <- length(c5_indel) # 34

n_indel_total <- n_c1_indel + n_c2_indel + n_c3_indel + n_c4_indel + n_c5_indel # 348

P_D_1_indel <- n_c1_indel/n_indel_total
P_D_2_indel <- n_c2_indel/n_indel_total
P_D_3_indel <- n_c3_indel/n_indel_total
P_D_4_indel <- n_c4_indel/n_indel_total
P_D_5_indel <- n_c5_indel/n_indel_total

P_D_i_indel <- c(P_D_1_indel, P_D_2_indel, P_D_3_indel,P_D_4_indel,P_D_5_indel)

# Function to calculate the likelihood and the posterior
calculate_likelihood_and_posterior <- function(S_t_matrix, D_i, P_D_i) {
  n_samples <- ncol(S_t_matrix)
  n_clusters <- length(D_i) 
  
  likelihoods <- matrix(0, nrow = n_clusters, ncol = n_samples)
  posteriors <- matrix(0, nrow = n_clusters, ncol = n_samples)
  
  for (sample_idx in 1:n_samples) {
    S_t <- S_t_matrix[, sample_idx] # Gets the mutational spectrum for the current sample
    
    # Calculate log-likelihoods for each cluster
    log_P_S_t_given_D <- sapply(1:n_clusters, function(cluster_idx) {
      # Use the log of the probability distribution for the current cluster
      cluster_D_i <- D_i[[cluster_idx]]
      
      # Avoid zero probabilities by replacing them with a small value
      cluster_D_i_safe <- ifelse(cluster_D_i == 0, 1e-10, cluster_D_i)
      
      # Sum the log of probabilities for observed mutations
      sum_log_probabilities <- sum(S_t * log(cluster_D_i_safe))
      
      return(sum_log_probabilities) 
    })
    
    # Store the log-likelihoods
    likelihoods[, sample_idx] <- exp(log_P_S_t_given_D) # Convert back to likelihood space
    
    # Calculate the log of the total probability of observed spectrum using log-sum-exp
    max_log_likelihood <- max(log_P_S_t_given_D) 
    log_P_S_t <- log(sum(exp(log_P_S_t_given_D - max_log_likelihood))) + max_log_likelihood
    
    # Calculate posterior probabilities P(Di|St) in log-space
    log_P_D_given_S <- log_P_S_t_given_D + log(P_D_i) - log_P_S_t
    
    # Convert log posterior probabilities back to regular space 
    posteriors[, sample_idx] <- exp(log_P_D_given_S)
    
    # Normalize to ensure the sum of posteriors equals 1 
    posteriors[, sample_idx] <- posteriors[, sample_idx] / sum(posteriors[, sample_idx])
  }
  
  return(list(
    likelihoods = likelihoods,
    posteriors = posteriors
  ))
}

# Apply the direct likelihood calculation to SBS and INDEL clusters
likelihood_posterior_sbs <- calculate_likelihood_and_posterior(S_t_matrix = simulated_wes_sbs_mut_matrix, D_i = D_i_sbs, P_D_i = P_D_i_sbs)
likelihood_posterior_indel <- calculate_likelihood_and_posterior(S_t_matrix = simulated_wes_indel_mut_matrix, D_i = D_i_indel, P_D_i = P_D_i_indel)
