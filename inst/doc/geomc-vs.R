## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = FALSE,
  tidy = FALSE,
  comment = "#>",
  progress = FALSE,
  echo = TRUE,
  dev = 'pdf',
  dpi = 120,
  fig.width = 6,
  fig.height = 4,
  out.width = "100%"
)

## ----setup--------------------------------------------------------------------
library(geommc)

## ----basic-usage, eval=FALSE--------------------------------------------------
# result <- geomc.vs(
#   X = X,              # n \times p predictor matrix
#   y = y               # n-vector response
# )

## ----data-1-------------------------------------------------------------------
# Set parameters
n <- 50          # Sample size
p <- 100         # Number of predictors
nonzero <- 3     # Number of true predictors
trueidx <- 1:3   # Indices of true predictors
nonzero.value <- 4  # Effect size

# True regression coefficients
TrueBeta <- numeric(p)
TrueBeta[trueidx] <- nonzero.value

# Generate correlated predictors
rho <- 0.5
set.seed(3)
xone <- matrix(rnorm(n * p), n, p)
X <- sqrt(1 - rho) * xone + sqrt(rho) * rnorm(n)

# Generate response
intercept <- 0.5
sigma_true <- 1
y <- intercept + X %*% TrueBeta + rnorm(n, 0, sigma_true)

cat("Data dimensions:\n")
cat("  n (observations):", n, "\n")
cat("  p (predictors):", p, "\n")
cat("  Number of true non-zero predictors:", nonzero, "\n")
cat("  True coefficient indices:", paste(trueidx, collapse = ", "), "\n")
cat("  True coefficient values:", nonzero.value, "\n")

## ----run-vs-1-----------------------------------------------------------------
# Run geomc.vs
set.seed(123)
result <- geomc.vs(
  X = X,
  y = y,
  model.summary = TRUE # Compute model summaries
)

## ----results-structure--------------------------------------------------------
names(result)

## ----mip-1--------------------------------------------------------------------
# Top 10 variables by MIP
top_vars <- order(result$mip, decreasing = TRUE)[1:10]

cat("Top 10 variables by Marginal Inclusion Probability:\n\n")
mip_df <- data.frame(
  Variable = top_vars,
  MIP = round(result$mip[top_vars], 4),
  True = ifelse(top_vars %in% trueidx, "Yes", "No")
)
print(mip_df, row.names = FALSE)

## ----plot-mip-1, fig.width=6, fig.height=4------------------------------------
par(mfrow = c(1, 1))

# Plot MIPs
plot(1:p, result$mip, 
     type = "h", lwd = 2, col = "gray",
     xlab = "Variable Index", 
     ylab = "MIP",
     main = "Marginal Inclusion Probabilities",
     ylim = c(0, 1))

# Highlight true variables
points(trueidx, result$mip[trueidx], 
       col = "red", pch = 19, cex = 2)

# Add threshold line
abline(h = 0.5, col = "blue", lty = 2, lwd = 2)

legend("topright", 
       legend = c("True variables", "MIP > 0.5 threshold"),
       col = c("red", "blue"), 
       pch = c(19, NA), 
       lty = c(NA, 2),
       lwd = 2)

## ----median-model-1-----------------------------------------------------------
cat("Median Probability Model:\n")
cat("  Number of variables selected:", length(result$median.model), "\n")
cat("  Variable indices:", paste(result$median.model, collapse = ", "), "\n\n")

# Check if we recovered the true model
correctly_identified <- all(trueidx %in% result$median.model)
false_positives <- setdiff(result$median.model, trueidx)
false_negatives <- setdiff(trueidx, result$median.model)

cat("Model Recovery:\n")
cat("  All true variables identified:", correctly_identified, "\n")
cat("  Number of false positives:", length(false_positives), "\n")
cat("  Number of false negatives:", length(false_negatives), "\n")

if (length(false_positives) > 0) {
  cat("  False positive indices:", paste(false_positives, collapse = ", "), "\n")
}

## ----coef-estimates-1---------------------------------------------------------
# Compare estimated vs true coefficients for true variables
cat("Coefficient Estimates (True Variables):\n\n")
coef_comparison <- data.frame(
  Variable = trueidx,
  True_Beta = TrueBeta[trueidx],
  Posterior_Mean = round(result$beta.mean[(1+trueidx)], 3),
  Beta_WAM = round(result$beta.wam[(1+trueidx)], 3),
  MIP = round(result$mip[trueidx], 4)
)
print(coef_comparison)
intercept_comparison<- data.frame(
  True_intercept =intercept,
  Posterior_Mean = round(result$beta.mean[1], 3),
  Beta_WAM = round(result$beta.wam[1], 3)
)
print(intercept_comparison)

## ----model-space, fig.width=7, fig.height=5-----------------------------------
# Model sizes visited
model_sizes <- apply(result$samples, 1, sum)

hist(model_sizes, breaks = 20, 
     main = "Distribution of Model Sizes Visited",
     xlab = "Number of Variables in Model",
     col = "lightblue", border = "white")
abline(v = nonzero, col = "red", lwd = 3, lty = 2)
legend("topright", 
       legend = "True model size",
       col = "red", lty = 2, lwd = 2)

cat("\nModel Size Summary:\n")
cat("  Mean model size:", round(mean(model_sizes), 2), "\n")
cat("  Median model size:", median(model_sizes), "\n")
cat("  True model size:", nonzero, "\n")

## ----generate-sparse-data-----------------------------------------------------
# Problem dimensions
n <- 80        # Number of documents
p <- 200       # Vocabulary size (number of words)
sparsity <- 0.95  # 95% of entries are zero

# True important words (only 5 words matter)
true_words <- c(10, 25, 50, 100, 150)
true_effects <- c(2.5, -2.0, 3.0, -1.8, 2.2)

# Create sparse design matrix (document-term matrix)
# Most words don't appear in most documents
X_dense <- matrix(0, n, p)

# Add word counts (Poisson distributed when word appears)
set.seed(3)
for (i in 1:n) {
  # Each document has about 5% of words present
  present_words <- sample(1:p, size = round(p * (1 - sparsity)))
  X_dense[i, present_words] <- rpois(length(present_words), lambda = 2)
}

# Ensure true words appear with higher frequency
for (word in true_words) {
  appears_in <- sample(1:n, size = round(0.7 * n))
  X_dense[appears_in, word] <- rpois(length(appears_in), lambda = 3)
}

# Remove zero-variance columns
col_vars <- apply(X_dense, 2, var)
nonzero_cols <- which(col_vars > 0)

cat("Original number of columns:", p, "\n")
cat("Columns with zero variance:", sum(col_vars == 0), "\n")
cat("Columns retained:", length(nonzero_cols), "\n\n")

# Keep mapping from new to old indices
old_to_new <- rep(NA, p)
old_to_new[nonzero_cols] <- 1:length(nonzero_cols)

# Update X_dense to only include non-zero variance columns
X_dense <- X_dense[, nonzero_cols, drop = FALSE]

# Update true_words indices to reflect new column positions
true_words_old <- true_words
true_words <- match(true_words_old, nonzero_cols)
true_words <- true_words[!is.na(true_words)]  # Keep only those that survived

cat("Original true words:", paste(true_words_old, collapse = ", "), "\n")
cat("New true words indices:", paste(true_words, collapse = ", "), "\n")
cat("True words retained:", length(true_words), "out of", length(true_words_old), "\n\n")

# Update true_effects to match retained true words
true_effects <- true_effects[!is.na(match(true_words_old, nonzero_cols))]

# Update p to reflect actual number of columns
p <- ncol(X_dense)
cat("Updated p:", p, "\n\n")

library(Matrix)
# Convert to sparse matrix (dgCMatrix)
X_sparse <- Matrix(X_dense, sparse = TRUE)

# Check sparsity
actual_sparsity <- 1 - (nnzero(X_sparse) / (n * p))
cat("Actual sparsity:", round(actual_sparsity * 100, 1), "%\n")
cat("Non-zero elements:", nnzero(X_sparse), "\n")
cat("Matrix class:", class(X_sparse), "\n")

## ----visualize-sparse, fig.width=6, fig.height=4------------------------------
# Visualize the sparse matrix pattern
image(
    x = 1:min(100, p),
    y = 1:min(50, n),
    z = as.matrix(t(X_sparse[1:min(50, n), 1:min(100, p)])),
    main = "Sparse Matrix Pattern (50 docs × 100 words)",
    ylab = "Documents",
    xlab = "Words",
    col = colorRampPalette(c("white", "blue", "darkblue"))(100)
)

## ----memory-comparison--------------------------------------------------------
# Compare memory usage
size_sparse <- object.size(X_sparse)
size_dense <- object.size(X_dense)

cat("Memory usage:\n")
cat("  Sparse matrix:", format(size_sparse, units = "Kb"), "\n")
cat("  Dense matrix:", format(size_dense, units = "Kb"), "\n")
cat("  Reduction:", 
    round(100 * (1 - as.numeric(size_sparse) / as.numeric(size_dense)), 1), 
    "%\n")

## ----generate-response--------------------------------------------------------
# Create coefficient vector (now with updated p)
beta_true <- numeric(p)
beta_true[true_words] <- true_effects

# Generate response (document sentiment/category)
set.seed(3)
linear_predictor <- X_sparse %*% beta_true
y <- as.numeric(linear_predictor) + rnorm(n, sd = 1.5)

cat("Response variable summary:\n")
print(summary(y))

## ----run-geomc-vs-------------------------------------------------------------
# Run variable selection
set.seed(123)
result <- geomc.vs(
  X = X_sparse,         # Sparse dgCMatrix
  y = y,
  model.summary = TRUE # Compute model summaries (default: FALSE)
)

## ----samples-explained--------------------------------------------------------
# Dimensions
cat("Sample matrix dimensions:", dim(result$samples), "\n")
cat("(rows = iterations, columns = variables)\n\n")

# First few iterations
cat("First 5 iterations (first 10 variables):\n")
print(result$samples[1:5, 1:20])

# Each row is a model
cat("\nIteration 1 includes variables:", which(result$samples[1, ] == 1), "\n")
cat("Iteration 2 includes variables:", which(result$samples[2, ] == 1), "\n")

# The proportion of accepted proposals
cat("\nAcceptance rate:", round(result$acceptance.rate, 3), "\n")

## ----logpost, fig.width=6, fig.height=4---------------------------------------
# Trace plot of log posterior
plot(result$log.p, type = "l",
     xlab = "Iteration", ylab = "Log Posterior",
     main = "Trace of Log Posterior Probability")

cat("\nLog posterior summary:\n")
cat("  Min:", round(min(result$log.p), 2), "\n")
cat("  Max:", round(max(result$log.p), 2), "\n")
cat("  Mean:", round(mean(result$log.p), 2), "\n")

## ----mip-analysis-------------------------------------------------------------
# Identify top words by MIP
top_n <- 20
top_indices <- order(result$mip, decreasing = TRUE)[1:top_n]

cat("Top", top_n, "words by Marginal Inclusion Probability:\n\n")
cat("(Note: Word indices are in the reduced matrix space)\n\n")

mip_table <- data.frame(
  Rank = 1:top_n,
  Word_Index = top_indices,
  MIP = round(result$mip[top_indices], 4),
  True_Word = ifelse(top_indices %in% true_words, "Yes", "No"),
  True_Effect = ifelse(
    top_indices %in% true_words,
    sprintf("%.2f", beta_true[top_indices]),
    ""
  )
)
print(mip_table, row.names = FALSE)

# Check recovery
cat("\nRecovery of true words:\n")
for (i in seq_along(true_words)) {
  word_idx <- true_words[i]
  cat(sprintf("  Word %3d (original %3d): MIP = %.4f, Rank = %2d\n",
              word_idx,
              true_words_old[i],
              result$mip[word_idx],
              which(order(result$mip, decreasing = TRUE) == word_idx)))
}

## ----plot-mips, fig.width=7, fig.height=6-------------------------------------
par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))

# Plot 1: All MIPs
plot(
  1:p, result$mip,
  type = "h", col = "gray70", lwd = 1,
  xlab = "Word Index (in reduced matrix)",
  ylab = "MIP",
  main = "Marginal Inclusion Probabilities (All Words)",
  ylim = c(0, 1)
)
points(true_words, result$mip[true_words], 
       col = "red", pch = 19, cex = 2)
abline(h = 0.5, col = "blue", lty = 2, lwd = 2)
legend("topright",
       legend = c("True word", "MIP > 0.5"),
       col = c("red", "blue"),
       pch = c(19, NA),
       lty = c(NA, 2),
       lwd = 2)


## ----plot-wmips, fig.width=7, fig.height=6------------------------------------
par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))

# Plot 2: All WMIPs
plot(
  1:p, result$wmip,
  type = "h", col = "gray70", lwd = 1,
  xlab = "Word Index (in reduced matrix)",
  ylab = "WMIP",
  main = "Weighted Marginal Inclusion Probabilities (All Words)",
  ylim = c(0, 1)
)
points(true_words, result$wmip[true_words], 
       col = "red", pch = 19, cex = 2)
abline(h = 0.5, col = "blue", lty = 2, lwd = 2)
legend("topright",
       legend = c("True word", "WMIP > 0.5"),
       col = c("red", "blue"),
       pch = c(19, NA),
       lty = c(NA, 2),
       lwd = 2)


## ----model-performance--------------------------------------------------------
# Median probability model
selected_words <- result$median.model

cat("Median Probability Model:\n")
cat("  Words selected:", length(selected_words), "\n")
cat("  Word indices (in reduced matrix):", paste(selected_words, collapse = ", "), "\n\n")

# Confusion matrix
true_positive <- sum(true_words %in% selected_words)
false_positive <- length(selected_words) - true_positive
false_negative <- length(true_words) - true_positive
true_negative <- p - length(true_words) - false_positive

cat("Classification Performance:\n")
cat("  True Positives:", true_positive, "/", length(true_words), "\n")
cat("  False Positives:", false_positive, "\n")
cat("  False Negatives:", false_negative, "\n")
cat("  Sensitivity:", round(true_positive / length(true_words), 3), "\n")
cat("  Precision:", 
    round(true_positive / max(1, length(selected_words)), 3), "\n")

## ----wam-model-performance----------------------------------------------------
# Weighted average model
selected_words <- result$wam

cat("Weighted Average Model:\n")
cat("  Words selected:", length(selected_words), "\n")
cat("  Word indices (in reduced matrix):", paste(selected_words, collapse = ", "), "\n\n")

# Confusion matrix
true_positive <- sum(true_words %in% selected_words)
false_positive <- length(selected_words) - true_positive
false_negative <- length(true_words) - true_positive
true_negative <- p - length(true_words) - false_positive

cat("Classification Performance:\n")
cat("  True Positives:", true_positive, "/", length(true_words), "\n")
cat("  False Positives:", false_positive, "\n")
cat("  False Negatives:", false_negative, "\n")
cat("  Sensitivity:", round(true_positive / length(true_words), 3), "\n")
cat("  Precision:", 
    round(true_positive / max(1, length(selected_words)), 3), "\n")



## ----coefficient-estimates----------------------------------------------------
cat("Coefficient Estimates for True Words:\n\n")

coef_table <- data.frame(
  Word_New = true_words,
  Word_Original = true_words_old[1:length(true_words)],
  True_Effect = true_effects,
  Post_Mean = round(result$beta.mean[(1+true_words)], 3),
  Beta_WAM = round(result$beta.wam[(1+true_words)], 3),
  MIP = round(result$mip[true_words], 4),
  WMIP = round(result$wmip[true_words], 4),
  Selected = ifelse(true_words %in% selected_words, "Yes", "")
)
print(coef_table, row.names = FALSE)

## ----sparse-model-space, fig.width=7, fig.height=6----------------------------
# Sizes of visited models
model_sizes <- apply(result$samples, 1, sum)

hist(
  model_sizes,
  breaks = 30,
  main = "Distribution of Model Sizes Visited",
  xlab = "Number of Words in Model",
  ylab = "Frequency",
  col = "lightblue",
  border = "white"
)
abline(v = length(true_words), col = "red", lwd = 3, lty = 2)
abline(v = median(model_sizes), col = "blue", lwd = 2, lty = 2)
legend("topleft",
       legend = c(
         paste("True size =", length(true_words)),
         paste("Median =", median(model_sizes))
       ),
       col = c("red", "blue"),
       lty = 2,
       lwd = 2)

cat("\nModel size summary:\n")
cat("  Mean:", round(mean(model_sizes), 2), "\n")
cat("  Median:", median(model_sizes), "\n")
cat("  Range:", range(model_sizes), "\n")

## ----run-geomc-vs-asym--------------------------------------------------------
# Run variable selection
set.seed(123)
result_asym <- geomc.vs(
  X = X_sparse,         # Sparse dgCMatrix
  y = y,
  symm = FALSE,
  move.prob = c(0.4,0.4,0.2),
  model.summary = TRUE # Compute model summaries (default: FALSE)
)

## ----asym-analysis------------------------------------------------------------
cat("\nAcceptance rate:", round(result_asym$acceptance.rate, 3), "\n")
# Identify top words by MIP
top_n <- 20
top_indices <- order(result_asym$mip, decreasing = TRUE)[1:top_n]

cat("Top", top_n, "words by Marginal Inclusion Probability:\n\n")
cat("(Note: Word indices are in the reduced matrix space)\n\n")

mip_table <- data.frame(
  Rank = 1:top_n,
  Word_Index = top_indices,
  MIP = round(result_asym$mip[top_indices], 4),
  True_Word = ifelse(top_indices %in% true_words, "Yes", "No"),
  True_Effect = ifelse(
    top_indices %in% true_words,
    sprintf("%.2f", beta_true[top_indices]),
    ""
  )
)
print(mip_table, row.names = FALSE)

# Check recovery
cat("\nRecovery of true words:\n")
for (i in seq_along(true_words)) {
  word_idx <- true_words[i]
  cat(sprintf("  Word %3d (original %3d): MIP = %.4f, Rank = %2d\n",
              word_idx,
              true_words_old[i],
              result_asym$mip[word_idx],
              which(order(result_asym$mip, decreasing = TRUE) == word_idx)))
}

## ----asym-plot-mips, fig.width=7, fig.height=6--------------------------------
par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))

# Plot 1: All MIPs
plot(
  1:p, result_asym$mip,
  type = "h", col = "gray70", lwd = 1,
  xlab = "Word Index (in reduced matrix)",
  ylab = "MIP",
  main = "Marginal Inclusion Probabilities (All Words)",
  ylim = c(0, 1)
)
points(true_words, result_asym$mip[true_words], 
       col = "red", pch = 19, cex = 2)
abline(h = 0.5, col = "blue", lty = 2, lwd = 2)
legend("topright",
       legend = c("True words", "MIP > 0.5"),
       col = c("red", "blue"),
       pch = c(19, NA),
       lty = c(NA, 2),
       lwd = 2)

## ----asym-model-performance---------------------------------------------------
# Median probability model
selected_words <- result_asym$median.model

cat("Median Probability Model:\n")
cat("  Words selected:", length(selected_words), "\n")
cat("  Word indices (in reduced matrix):", paste(selected_words, collapse = ", "), "\n\n")

# Confusion matrix
true_positive <- sum(true_words %in% selected_words)
false_positive <- length(selected_words) - true_positive
false_negative <- length(true_words) - true_positive
true_negative <- p - length(true_words) - false_positive

cat("Classification Performance:\n")
cat("  True Positives:", true_positive, "/", length(true_words), "\n")
cat("  False Positives:", false_positive, "\n")
cat("  False Negatives:", false_negative, "\n")
cat("  Sensitivity:", round(true_positive / length(true_words), 3), "\n")
cat("  Precision:", 
    round(true_positive / max(1, length(selected_words)), 3), "\n")

## ----asym-coefficient-estimates-----------------------------------------------
cat("Coefficient Estimates for True Words:\n\n")

coef_table <- data.frame(
  Word_New = true_words,
  Word_Original = true_words_old[1:length(true_words)],
  True_Effect = true_effects,
  Post_Mean = round(result_asym$beta.mean[(1+true_words)], 3),
  Beta_WAM = round(result_asym$beta.wam[(1+true_words)], 3),
  MIP = round(result_asym$mip[true_words], 4),
  WMIP = round(result_asym$wmip[true_words], 4),
  Selected = ifelse(true_words %in% selected_words, "Yes", "")
)
print(coef_table, row.names = FALSE)

## ----asym-model-space, fig.width=7, fig.height=6------------------------------
# Sizes of visited models
model_sizes <- apply(result_asym$samples, 1, sum)

hist(
  model_sizes,
  breaks = 30,
  main = "Distribution of Model Sizes Visited",
  xlab = "Number of Words in Model",
  ylab = "Frequency",
  col = "lightblue",
  border = "white"
)
abline(v = length(true_words), col = "red", lwd = 3, lty = 2)
abline(v = median(model_sizes), col = "blue", lwd = 2, lty = 2)
legend("topleft",
       legend = c(
         paste("True size =", length(true_words)),
         paste("Median =", median(model_sizes))
       ),
       col = c("red", "blue"),
       lty = 2,
       lwd = 2)

cat("\nModel size summary:\n")
cat("  Mean:", round(mean(model_sizes), 2), "\n")
cat("  Median:", median(model_sizes), "\n")
cat("  Range:", range(model_sizes), "\n")

## ----tuning, eval=FALSE-------------------------------------------------------
# result <- geomc.vs(
#   X = X,
#   y = y,
#   initial = c(1,2),         # Initial model (the set of active variables)
#   n.iter = 100,        # Number of MCMC iterations
#   burnin = 2,           # Burn-in period (default: 1)
#   eps = 0.5,            # Perturbation parameter (default: 0.5)
#   symm = FALSE,          # Use symmetric proposals (default: TRUE)
#   move.prob = c(.3, .3, .4),  # Probabilities for add/delete/swap moves
#   lam0 = 0,             # Prior parameter (default: 0)
#   a0 = 0,               # Prior shape parameter (default: 0)
#   b0 = 0,               # Prior scale parameter (default: 0)
#   lam = n / p^2,        # Prior parameter (default: n/p^2)
#   w = sqrt(n) / p,      # Prior inclusion probability (default: sqrt(n)/p)
#   model.summary = TRUE, # Compute model summaries (default: FALSE)
#   model.threshold = 0.5, # Threshold for median model (default: 0.5)
#   show.progress = TRUE # Show progress during sampling
# )

