### Train the model on the artificially simulated WES samples
# 1) Train GBC-1 on the whole features_df with CV, to get the best number of trees
# 2) Make a train and test dataset from the features_df and train GBC-2 on the train set, without CV using the number of trees determiend in 1)
# 3) Train GBC-3 on the whole dataset again without CV, using th enumber of trees detemrmined in 1)

# Load libraries for gradient boosting model
library(gbm3)
library(caret)

# Train GBC-1 with cross-validation to determine optimal parameters
set.seed(123)
gbm_model_1 <- gbm(
  label ~ .,
  data = features_df %>% select(-sample_names),
  distribution = "bernoulli",
  n.trees = 10000,         # Upper limit for CV
  interaction.depth = 3,   # Depth of each tree
  shrinkage = 0.01,        # Learning rate
  bag.fraction = 0.2,      # Bagging fraction
  cv.folds = 5,            # 5-fold cross-validation
  verbose = FALSE
)

# Optimal number of trees based on CV
best_trees <- gbm.perf(gbm_model_1, method = "cv", plot.it = TRUE)

# Stratified splitting based on ERCC2-MUT/WT
set.seed(123)
train_index <- createDataPartition(features_df$label, p = 0.8, list = FALSE, times = 1)
train_data <- features_df[train_index, ] # 400 samples
test_data <- features_df[-train_index, ] # 100 samples, 50 wt - 50 mut

# Train GBC-2 on the training set
gbm_model_2 <- gbm(
  label ~ .,
  data = train_data %>% select(-sample_names),
  distribution = "bernoulli",
  n.trees = best_trees,    # Optimal number of trees from CV
  interaction.depth = 3,
  shrinkage = 0.01,
  bag.fraction = 0.2
)

# Evaluate GBC-2 on the test set
test_preds <- predict(gbm_model_2, newdata = test_data, n.trees = best_trees, type = "response")

# Train GBC-3 on the entire dataset
set.seed(123)
gbm_model_3 <- gbm(
  label ~ .,
  data = features_df %>% select(-sample_names),
  distribution = "bernoulli",
  n.trees = best_trees,    
  interaction.depth = 3,
  shrinkage = 0.01,
  bag.fraction = 0.2
)

# Predict probabilities for the entire dataset
scores_gbc3 <- predict(
  gbm_model_3, 
  newdata = features_df %>% select(-sample_names), 
  n.trees = best_trees, 
  type = "response"
)

# Find threshold for assigning labels after prediction
predicted_scores_whole_datagbc3 <- scores_gbc3
true_labels <- features_df$label

# Create a data frame for use with cutpointr
scores_labels_data <- data.frame(scores = predicted_scores_whole_datagbc3, labels = true_labels)

# Use cutpointr to calculate the optimal cutoff
result <- cutpointr(scores_labels_data, x = scores, class = labels, 
                    method = maximize_metric, metric = youden)

plot(result)

# Optimal cutpoint
optimal_cutpoint <- result$optimal_cutpoint

# Assign labels to predictions based on optimal cutpoint 
predicted_labels_features_df <- ifelse(predicted_scores_whole_datagbc3 >= optimal_cutpoint, 1, 0)

# Evaluate performance
evaluate_gbm3 <- confusionMatrix(as.factor(predicted_labels_features_df), as.factor(true_labels))

# Extract feature importance
feature_importance_gbm3 <- summary(gbm_model_3, plotit = TRUE)

### Pots and metrics for the model

# Predict scores for GBC-1 and GBC-3 on the test dataset
scores_gbc1_test <- predict(gbm_model_1, newdata = test_data, n.trees = best_trees, type = "response")
scores_gbc2_test <- test_preds 
scores_gbc3_test <- predict(gbm_model_3, newdata = test_data, n.trees = best_trees, type = "response")

# Combine predictions into a data frame
scatter_data <- data.frame(
  GBC1 = scores_gbc1_test,
  GBC2 = scores_gbc2_test,
  GBC3 = scores_gbc3_test,
  Label = test_data$label
)

# Map labels to ERCC2-WT and ERCC2-MUT
scatter_data <- scatter_data %>%
  mutate(Label_Desc = ifelse(Label == 1, "ERCC2-MUT", "ERCC2-WT"))

scatter_data <- scatter_data %>%
  mutate(
    Region = case_when(
      GBC3 >= optimal_cutpoint & GBC2 < optimal_cutpoint ~ "Region-A",
      GBC3 >= optimal_cutpoint & GBC2 >= optimal_cutpoint ~ "Region-B",
      GBC3 < optimal_cutpoint & GBC2 < optimal_cutpoint ~ "Region-C",
      GBC3 < optimal_cutpoint & GBC2 >= optimal_cutpoint ~ "Region-D"
    )
  )

# Scatterplot: GBC-1 vs GBC-2
ggplot(scatter_data, aes(x = GBC1, y = GBC2, color = Label_Desc)) +
  geom_point() +
  geom_vline(xintercept = optimal_cutpoint, linetype = "dashed", color = "black", size = 0.5) +
  geom_hline(yintercept = optimal_cutpoint, linetype = "dashed", color = "black", size = 0.5) +
  theme_minimal() +
  labs(
    x = "GBC-1 Scores",
    y = "GBC-2 Scores",
    color = "Label",
    title = "GBC-1 vs GBC-2 Scores"
  ) +
  scale_color_manual(values = c("ERCC2-MUT" = "blue", "ERCC2-WT" = "red"))

# Scatterplot: GBC-2 vs GBC-3 with regions
ggplot(scatter_data, aes(x = GBC2, y = GBC3, color = Label_Desc)) +
  geom_point() +
  geom_vline(xintercept = optimal_cutpoint, linetype = "dashed", color = "black", size = 0.5) +
  geom_hline(yintercept = optimal_cutpoint, linetype = "dashed", color = "black", size = 0.5) +
  theme_minimal() +
  labs(
    x = "GBC-2 Scores",
    y = "GBC-3 Scores",
    color = "Region",
    title = "GBC-2 vs GBC-3 Scores"
  ) +
  scale_color_manual(values = c("ERCC2-MUT" = "blue", "ERCC2-WT" = "red")) +
  annotate("text", x = 0.3, y = 0.8, label = "Region-B", color = "black", size = 4, hjust = 2) +
  annotate("text", x = 0.3, y = 0.3, label = "Region-C", color = "black", size = 4, hjust = 2) +
  annotate("text", x = 0.8, y = 0.8, label = "Region-A", color = "black", size = 4, hjust = 2) +
  annotate("text", x = 0.8, y = 0.3, label = "Region-D", color = "black", size = 4, hjust = 2)

