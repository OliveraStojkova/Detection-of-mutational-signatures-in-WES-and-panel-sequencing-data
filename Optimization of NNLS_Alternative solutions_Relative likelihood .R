# Optimization calculation of exposures with NNLS to get the best signature composition
# Get all possible signature compositions for each sample
# Relative likelihood calculation (best/alternativ solutions), 

# Function for iterative NNLS decomposition, max_signatures = 5 for WES and panels, 6 for WGS
iterative_nnls <- function(cosmic_signatures, mutational_spectrum, max_signatures = 5) {
  
  # Iterate through all samples
  exposures <- matrix(0, ncol = ncol(cosmic_signatures), nrow = ncol(mutational_spectrum))
  colnames(exposures) <- colnames(cosmic_signatures)
  rownames(exposures) <- colnames(mutational_spectrum)
  
  errors <- numeric(ncol(mutational_spectrum))
  best_combinations <- vector("list", ncol(mutational_spectrum))
  
  # Iterate through all samples
  for (sample_idx in 1:ncol(mutational_spectrum)) {
    sample_spectrum <- mutational_spectrum[, sample_idx]
    best_exposure <- NULL
    best_error <- Inf
    best_signatures <- NULL
    
    # Increases the number of signatures with each iteration
    for (num_signatures in 2:max_signatures) {
      # Generate all combinations of signatures
      signature_combinations <- combn(colnames(cosmic_signatures), num_signatures)
      
      for (combo_idx in 1:ncol(signature_combinations)) {
        # Get the current combination of signatures
        current_signatures <- signature_combinations[, combo_idx]
        subset_signatures <- cosmic_signatures[, current_signatures, drop = FALSE]
        
        # Perform NNLS decomposition
        fit <- nnls(subset_signatures, sample_spectrum)
        reconstructed_spectrum <- subset_signatures %*% coef(fit)
        
        # Calculate reconstruction error
        error <- sum((sample_spectrum - reconstructed_spectrum)^2)
        
        # Update the best solution if the error is reduced
        if (error < best_error) {
          best_error <- error
          best_exposure <- coef(fit)
          best_signatures <- current_signatures
        }
      }
    }
    
    # Save the best exposure, error, and combination for the current sample
    exposures[sample_idx, best_signatures] <- best_exposure
    errors[sample_idx] <- best_error
    best_combinations[[sample_idx]] <- best_signatures
  }
  
  # Combine results into a list
  result <- list(
    exposures = exposures,
    errors = errors,
    best_combinations = best_combinations
  )
  
  return(result)
}

# Implement the function on SBS and INDEL count matrices
iterative_exposures_sbs <- iterative_nnls(cosmic_signature_names[, bladder_tissue_sbs], simulated_wes_sbs_mut_matrix, max_signatures = 5)
iterative_exposures_indel <- iterative_nnls(cosmic_signature_indel[, bladder_tissue_indel], simulated_wes_indel_mut_matrix, max_signatures = 5)

# Get SBS5 and ID8 exposures for common samples
sbs5_exposure <- iterative_exposures_sbs$exposures[, "SBS5"]
id8_exposure <- iterative_exposures_indel$exposures[, "ID8"]

# Function to generate all solutions
library(gtools)

get_alternative_solutions <- function(mutational_spectrum, signatures_matrix) {
  n_signatures <- ncol(signatures_matrix)
  n_samples <- ncol(mutational_spectrum)
  all_alternative_solutions <- vector("list", n_samples)
  
  # Generate combinations of 5 signatures at a time
  signature_combinations <- combinations(n_signatures, 5)
  
  # Loop through each sample
  for (sample_idx in 1:n_samples) {
    mutational_spectrum_sample <- mutational_spectrum[, sample_idx]
    sample_solutions <- list()
    
    # Loop through each combination of 5 signatures
    for (combination_idx in 1:nrow(signature_combinations)) {
      # Select the combination of signatures
      selected_signatures <- signature_combinations[combination_idx, ]
      selected_signatures_matrix <- signatures_matrix[, selected_signatures, drop = FALSE]
      
      # Recalculate exposures using NNLS
      nnls_result <- nnls(selected_signatures_matrix, mutational_spectrum_sample)
      
      # Store alternative solution
      sample_solutions[[combination_idx]] <- list(
        signatures = colnames(selected_signatures_matrix),
        exposures = nnls_result$x
      )
    }
    
    # Store all solutions for the current sample
    all_alternative_solutions[[sample_idx]] <- sample_solutions
  }
  
  return(all_alternative_solutions)
}

# Implement the function on SBS and INDEL count matrices 
alternative_solutions_sbs <- get_alternative_solutions(simulated_wes_sbs_mut_matrix, cosmic_signature_names[, bladder_tissue_sbs])
alternative_solutions_indel <- get_alternative_solutions(simulated_wes_indel_mut_matrix, cosmic_signature_indel[, bladder_tissue_indel])

# Calculate the likelihood of a given NNLS decomposition compared to other possible decompositions
calculate_likelihood_nnls <- function(cosmic_signatures, mutational_spectrum, alternative_solutions) {
  likelihoods <- numeric(ncol(mutational_spectrum))
  sample_names <- colnames(mutational_spectrum)  
  
  # Loop through each sample
  for (sample_idx in 1:ncol(mutational_spectrum)) {
    sample_spectrum <- mutational_spectrum[, sample_idx]
    
    # Get all alternative solutions for the current sample
    sample_solutions <- alternative_solutions[[sample_idx]]
    
    # Best solution is the one with the smallest residual error (NNLS fit)
    best_solution_idx <- which.min(sapply(sample_solutions, function(x) sum((sample_spectrum - (cosmic_signatures[, x$signatures] %*% x$exposures))^2)))
    
    # Extract the best solution signatures and exposures
    best_solution <- sample_solutions[[best_solution_idx]]
    best_signatures <- best_solution$signatures
    exposures <- best_solution$exposures
    
    # Reconstruct spectrum with the best solution - M
    selected_signatures <- cosmic_signatures[, best_signatures, drop = FALSE]
    M <- selected_signatures %*% exposures
    likelihood_with <- sum(dpois(sample_spectrum, M, log = TRUE))  # Log likelihood for M
    
    # Calculate likelihood for each possible alternative solution
    likelihood_withouts <- sapply(sample_solutions, function(alt_solution) {
      alternative_signatures <- alt_solution$signatures
      alt_exposures <- alt_solution$exposures
      
      alt_selected_signatures <- cosmic_signatures[, alternative_signatures, drop = FALSE]
      M_prime <- alt_selected_signatures %*% alt_exposures
      sum(dpois(sample_spectrum, M_prime, log = TRUE))  # Log likelihood for M'
    })
    
    # Relative likelihood: compare best solution with others
    likelihoods[sample_idx] <- exp(likelihood_with - max(likelihood_withouts))  # Ratio of likelihoods
  }
  
  # Combine sample names and likelihoods into a matrix
  likelihood_results <- matrix(
    likelihoods, 
    ncol = 1, 
    dimnames = list(sample_names, "likelihood")
  )
  
  return(likelihood_results)
}

# Implement the function to the SBS and INDEL count matrices
nnls_likelihood_ratio_sbs <- calculate_likelihood_nnls(cosmic_signature_names[, bladder_tissue_sbs], simulated_wes_sbs_mut_matrix, alternative_solutions_sbs)
nnls_likelihood_ratio_indel <- calculate_likelihood_nnls(cosmic_signature_indel[, bladder_tissue_indel], simulated_wes_indel_mut_matrix, alternative_solutions_indel)
