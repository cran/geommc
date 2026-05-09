## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = FALSE,
  tidy = FALSE,
  comment = "#>",
  progress = FALSE
)

## ----setup--------------------------------------------------------------------
# Install using install.packages('geommc')
library(geommc)

## ----mvnorm-example-----------------------------------------------------------
# Define the log-target density
log_target_mvnorm <- function(x, target.mean, target.Sigma) {
  d  <- length(x)
  xc <- x - target.mean
  Q  <- solve(target.Sigma)
  -0.5 * drop(t(xc) %*% Q %*% xc)
}

# Target distribution parameters
target_mean <- c(1, -2)
target_Sigma <- matrix(c(1.5, 0.7, 0.7, 2.0), 2, 2)

# Run geomc with default settings
set.seed(3)
result <- geomc(
  logp = log_target_mvnorm,
  initial = c(0, 0),
  n.iter = 1000,
  target.mean = target_mean,
  target.Sigma = target_Sigma
)

## ----examine-results----------------------------------------------------------
names(result)

## ----results-summary----------------------------------------------------------
# Sample means (should be close to target_mean)
cat("Sample means:\n")
print(colMeans(result$samples))
cat("\nTarget means:\n")
print(target_mean)

# Sample covariance (should be close to target_Sigma)
cat("\nSample covariance:\n")
print(cov(result$samples))
cat("\nTarget covariance:\n")
print(target_Sigma)

## ----plots-mvnorm, fig.width=7, fig.height=6----------------------------------
par(mfrow = c(3, 2))

# Trace plots
plot(result$samples[, 1], type = "l", 
     main = "Trace plot",
     ylab = expression(y[1]), xlab = "Iteration")
abline(h = target_mean[1], col = "red", lty = 2, lwd = 2)

plot(result$samples[, 2], type = "l", 
     main = "Trace plot",
     ylab = expression(y[2]), xlab = "Iteration")
abline(h = target_mean[2], col = "red", lty = 2, lwd = 2)

# Density plots
plot(density(result$samples[, 1]), 
     main = "Density",
     xlab = expression(y[1]))
abline(v = target_mean[1], col = "red", lty = 2, lwd = 2)

plot(density(result$samples[, 2]), 
     main = "Density",
     xlab = expression(y[2]))
abline(v = target_mean[2], col = "red", lty = 2, lwd = 2)

## ----acf-plot, fig.width=6, fig.height=4--------------------------------------
# Autocorrelation plot of samples
par(mfrow = c(1, 2))
acf(result$samples[, 1], main = "Autocorrelation plot", 
    xlab =    expression(y[1]), ylab = "Autocorrelation")
acf(result$samples[, 2], main = "Autocorrelation plot", 
    xlab = expression(y[2]), ylab = "Autocorrelation")

## ----scatter-plot, fig.width=6, fig.height=6----------------------------------
# Scatter plot of samples
plot(result$samples[, 1], result$samples[, 2],
     pch = 16, col = rgb(0, 0, 0, 0.3),
     xlab = expression(y[1]), ylab = expression(y[2]),
     main = "MCMC Samples")
points(target_mean[1], target_mean[2], col = "red", pch = 4, cex = 2, lwd = 3)
legend("topright", legend = "True mean", col = "red", pch = 4, lwd = 2)

## ----normal-inference---------------------------------------------------------
# Generate data
set.seed(42)
true_mu <- 5
true_sigma <- 2
n <- 100
w <- rnorm(n, mean = true_mu, sd = true_sigma)

cat("True parameters:\n")
cat("  mu =", true_mu, "\n")
cat("  sigma =", true_sigma, "\n")
cat("  sigma^2 =", true_sigma^2, "\n\n")

cat("Sample statistics:\n")
cat("  mean =", round(mean(w), 3), "\n")
cat("  sd =", round(sd(w), 3), "\n")

# Define log-posterior
log_target <- function(par, x, mu0, tau0, alpha0, beta0) {
  mu <- par[1]
  sigma2 <- par[2]
  
  # Check constraint
  if (sigma2 <= 0) return(-Inf)
  
  n <- length(x)
  SSE <- sum((x - mu)^2)
  val <- -(n/2) * log(sigma2) - SSE / (2 * sigma2)
  val <- val - (mu - mu0)^2 / (2 * tau0^2)
  val <- val - (alpha0 + 1) * log(sigma2) - beta0 / sigma2
  
  return(val)
}

# Hyperparameters (weakly informative priors)
mu0 <- 0      # prior mean for mu
tau0 <- 10    # prior sd for mu
alpha0 <- 2.01  # shape parameter for inverse-gamma
beta0 <- 1.01   # scale parameter for inverse-gamma

# Run geomc
set.seed(3)
result <- geomc(
  logp = log_target,
  initial = c(0, 1),  # Starting at mu=0, sigma^2=1
  n.iter = 1000,
  x = w,
  mu0 = mu0,
  tau0 = tau0,
  alpha0 = alpha0,
  beta0 = beta0
)
burnin<-50

## ----trace-plots-normal, fig.width=7, fig.height=6----------------------------
par(mfrow = c(1, 2))

# Trace plots
plot(result$samples[, 1], type = "l", 
     main = expression(paste("Trace plot: ", mu)),
     ylab = expression(mu), xlab = "Iteration")
abline(v = burnin, col = "gray", lty = 2)
abline(h = true_mu, col = "red", lty = 2, lwd = 2)

plot(result$samples[, 2], type = "l", 
     main = expression(paste("Trace plot: ", sigma^2)),
     ylab = expression(sigma^2), xlab = "Iteration")
abline(v = burnin, col = "gray", lty = 2)
abline(h = true_sigma^2, col = "red", lty = 2, lwd = 2)

## ----posterior-summaries------------------------------------------------------
burnin <- 50
# Remove burn-in (first 500 samples)
samples_post_burnin <- result$samples[-(1:burnin), ]

# Posterior mean estimates
mu_post_mean <- mean(samples_post_burnin[, 1])
sigma2_post_mean <- mean(samples_post_burnin[, 2])
sigma_post_mean <- mean(sqrt(samples_post_burnin[, 2]))

cat("Posterior mean estimates:\n")
cat("  mu:", round(mu_post_mean, 3), "\n")
cat("  sigma^2:", round(sigma2_post_mean, 3), "\n")
cat("  sigma:", round(sigma_post_mean, 3), "\n\n")

# 95% credible intervals
mu_ci <- quantile(samples_post_burnin[, 1], c(0.025, 0.975))
sigma2_ci <- quantile(samples_post_burnin[, 2], c(0.025, 0.975))

cat("95% Credible Intervals:\n")
cat("  mu: [", round(mu_ci[1], 3), ",", round(mu_ci[2], 3), "]\n")
cat("  sigma^2: [", round(sigma2_ci[1], 3), ",", round(sigma2_ci[2], 3), "]\n")

## ----posterior-plots, fig.width=7, fig.height=6-------------------------------
layout(matrix(c(1,2,3), nrow = 1), widths = c(4, 4, 2))
par(mar = c(4, 4, 2, 1))

plot(density(samples_post_burnin[,1]), main=expression(mu), xlab=expression(mu))
abline(v = true_mu, col="red", lty=2, lwd=2)
abline(v = mu_post_mean, col="blue", lty=2, lwd=2)

plot(density(samples_post_burnin[,2]), main=expression(sigma^2), xlab=expression(sigma^2))
abline(v = true_sigma^2, col="red", lty=2, lwd=2)
abline(v = sigma2_post_mean, col="blue", lty=2, lwd=2)

par(mar = c(0,0,0,0))
plot.new()
legend("center", legend=c("True value","Posterior mean est."),
       col=c("red","blue"), lty=2, lwd=2, bty="n")

## ----output-structure---------------------------------------------------------
#- **samples**: Matrix of MCMC samples (rows = iterations, columns = parameters)
head(result$samples)
#- **log.p**: Log-density values at each sample
# Trace plot of log-density values
plot(result$log.p[-(1:burnin)], type = "l", 
     main = "Trace plot of log-density values",
     ylab = "log.p", xlab = "Iteration")
#- **acceptance.rate**: Proportion of accepted proposals
cat("\nAcceptance rate:", round(result$acceptance.rate, 3), "\n")
#- **var.base**: Variance of the random walk base density used
cat("Variance of the random walk base=", result$var.base, "\n")
#- **mean.ap.tar**: Mean of the approximate target density
cat("Mean of the approximate target:\n")
print(result$mean.ap.tar)
#- **var.ap.tar**: Variance of the approximate target density
cat("\nVariance of the approximate target:\n")
print(result$var.ap.tar)
#- **model.case**: Which model configuration was used
cat("model configuration:\n")
print(result$model.case)
#- **gaus**: Whether Gaussian densities for both base and ap.tar were assumed
cat("Whether Gaussian densities for both base and ap.tar were assumed=", result$gaus, "\n")
#- **ind**: Whether independence was assumed for both base and ap.tar
cat("Whether independence was assumed for both base and ap.tar=", result$ind, "\n")

## ----mixture-default----------------------------------------------------------
# Define the log-target
log_target_mixture <- function(y) {
  log(0.5 * dnorm(y) + 0.5 * dnorm(y, mean = 10, sd = 1.4))
}

# Run geomc with default settings
set.seed(3)
result <- geomc(
  logp = list(log.target = log_target_mixture),
  initial = 0,
  n.iter = 1000
)

## ----mixture-default-plot, fig.width=7, fig.height=8--------------------------
par(mfrow = c(2, 1))

# Trace plot
plot(result$samples, type = "l", 
     main = "Trace plot: Mixture of Normals",
     ylab = "y", xlab = "Iteration")
abline(h = c(0, 10), col = c("blue", "red"), lty = 2)

# Histogram with true density overlay
hist(result$samples, breaks = 50, freq = FALSE,
     main = "Posterior samples vs True density",
     xlab = "y", col = "lightblue", border = "white")

# Add true density
y_grid <- seq(-5, 15, length.out = 1000)
true_density <- 0.5 * dnorm(y_grid) + 0.5 * dnorm(y_grid, mean = 10, sd = 1.4)
lines(y_grid, true_density, col = "red", lwd = 2)
legend("topright", legend = "True density", col = "red", lwd = 2)

## ----mixture-custom-base------------------------------------------------------
set.seed(123)
result_custom <- geomc(
  logp = list(
    log.target = log_target_mixture,
    mean.base = function(x) x,           # Centered at current state
    var.base = function(x) 4            # Variance = 4
  ),
  initial = 0,
  n.iter = 1000
)
#Note the variance of the random walk base density of the default setting
cat("Variance of the default geomc =", result$var.base, "\n")
cat("The acceptance rate of the default geomc:\n")
result$acceptance.rate
cat("The acceptance rate of geomc with the custom base density:\n")
result_custom$acceptance.rate

## ----mixture-custom-plot, fig.width=7, fig.height=4---------------------------
par(mfrow = c(1, 2))

hist(result$samples, breaks = 50, freq = FALSE,
     main = "Default settings",
     xlab = "y", col = "lightblue", border = "white", ylim = c(0, 0.25))
lines(y_grid, true_density, col = "red", lwd = 2)

hist(result_custom$samples, breaks = 50, freq = FALSE,
     main = "Custom random walk base",
     xlab = "y", col = "lightgreen", border = "white", ylim = c(0, 0.25))
lines(y_grid, true_density, col = "red", lwd = 2)

## ----mixture-informed---------------------------------------------------------
set.seed(123)
result_informed <- geomc(
  logp = list(
    log.target = log_target_mixture,
    # Specify the two densities g_1 and g_2
    mean.ap.tar = function(x) c(0, 10),  # Means of the two densities
    var.ap.tar = function(x) c(1, 1.4^2)  # Variances of the two densities
  ),
  initial = 0,
  n.iter = 1000
)

## ----mixture-informed-plot, fig.width=7, fig.height=8-------------------------
par(mfrow = c(3, 1))

plot(result$samples, type = "l", main = "Default settings", ylab = "y", 
     xlab = "Iteration")
abline(h = c(0, 10), col = c("blue", "red"), lty = 2)

plot(result_custom$samples, type = "l", main = "Custom random walk", 
     ylab = "y", xlab = "Iteration")
abline(h = c(0, 10), col = c("blue", "red"), lty = 2)

plot(result_informed$samples, type = "l", main = "Informed approximate targets", 
     ylab = "y", xlab = "Iteration")
abline(h = c(0, 10), col = c("blue", "red"), lty = 2)

## ----mixture-informed-apbase--------------------------------------------------
set.seed(123)
result_informed_apbase <- geomc(
  logp = list(
    log.target = log_target_mixture,
    mean.base = function(x) x, # specify the base mean
    var.base = function(x) 4, # specify the base variance
    # Specify the two densities g_1 and g_2
    mean.ap.tar = function(x) c(0, 10),  # Means of the two densities
    var.ap.tar = function(x) c(1, 1.4^2)  # Variances of the two densities
      ),
  initial = 0,
  n.iter = 1000
)
# Compare the values of model.case
cat("The model configuration under default setting=\n")
print(result$model.case)
cat("The model configuration under result_custom\n")
print(result_custom$model.case)
cat("The model configuration under result_informed\n")
print(result_informed$model.case)
cat("The model configuration under result_informed_apbase\n")
print(result_informed_apbase$model.case)

## ----mixture-informed-apbase-plot, fig.width=7, fig.height=8------------------
par(mfrow = c(3, 1))

plot(result$samples, type = "l", main = "Default settings", ylab = "y", 
     xlab = "Iteration")
abline(h = c(0, 10), col = c("blue", "red"), lty = 2)

plot(result_informed$samples, type = "l", main = "Informed approximate targets", 
     ylab = "y", xlab = "Iteration")
abline(h = c(0, 10), col = c("blue", "red"), lty = 2)
plot(result_informed_apbase$samples, type = "l", 
     main = "Informed approximate targets and custom base", 
     ylab = "y", xlab = "Iteration")
abline(h = c(0, 10), col = c("blue", "red"), lty = 2)


## ----mixture-eps--------------------------------------------------------------
set.seed(123)
result_informed_apbase_eps_large <- geomc(
  logp = list(
    log.target = log_target_mixture,
    mean.base = function(x) x,
    var.base = function(x) 4,
    # Specify the two densities g_1 and g_2
    mean.ap.tar = function(x) c(0, 10),  # Means of the two densities
    var.ap.tar = function(x) c(1, 1.4^2)  # Variances of the two densities
      ),
  initial = 0,
  n.iter = 1000,
  eps = 0.9
)

set.seed(123)
result_informed_apbase_eps_small <- geomc(
  logp = list(
    log.target = log_target_mixture,
    mean.base = function(x) x,
    var.base = function(x) 4,
    # Specify the two densities g_1 and g_2
    mean.ap.tar = function(x) c(0, 10),  # Means of the two densities
    var.ap.tar = function(x) c(1, 1.4^2)  # Variances of the two densities
      ),
  initial = 0,
  n.iter = 1000,
  eps = 0.01
)

set.seed(123)
result_rwm_mix <- geommc:::rwm(
  log_target_mixture,
  initial = 0,
  n_iter = 1000,
  sig = 4,
  return_sample = TRUE
)
cat("Geometric MCMC:\n")
cat("With default eps:\n")
cat("  Acceptance rate:", round(result_informed_apbase$acceptance.rate, 3), "\n")
cat("With eps =0.9:\n")
cat("  Acceptance rate:", round(result_informed_apbase_eps_large$acceptance.rate, 3), "\n")
cat("With eps=0.01:\n")
cat("  Acceptance rate:", round(result_informed_apbase_eps_small$acceptance.rate, 3), "\n")
cat("Random Walk Metropolis:\n")
cat("  Acceptance rate:", round(result_rwm_mix$acceptance_rate, 3), "\n")

## ----mixture-informed-apbase-eps-plot, fig.width=7, fig.height=11-------------
par(mfrow = c(2, 1))
plot(result_informed_apbase$samples, type = "l", 
     main = "Informed approximate targets and custom base (Default eps)", 
     ylab = "y", xlab = "Iteration")
abline(h = c(0, 10), col = c("blue", "red"), lty = 2)
plot(result_informed_apbase_eps_large$samples, type = "l", 
     main = "Informed approximate targets and custom base (eps=0.9)", 
     ylab = "y", xlab = "Iteration")
abline(h = c(0, 10), col = c("blue", "red"), lty = 2)
plot(result_informed_apbase_eps_small$samples, type = "l", 
     main = "Informed approximate targets and custom base (eps=0.01)", 
     ylab = "y", xlab = "Iteration")
abline(h = c(0, 10), col = c("blue", "red"), lty = 2)
plot(result_rwm_mix$samples, type = "l", 
     main = "Random walk Matropolis", 
     ylab = "y", xlab = "Iteration")
abline(h = c(0, 10), col = c("blue", "red"), lty = 2)


## ----mixture-informed_another-------------------------------------------------
# Sampler that draws from the mixture density g
samp.ap.tar_mixture <- function(x, kk = 1) {
  component <- sample.int(2, 1, prob = c(0.5, 0.5))
  if (component == 1) {
    return(rnorm(1))
  } else {
    return(10 + 1.4 * rnorm(1))
  }
}
set.seed(123)
result_informed_another <- geomc(
  logp = list(
    log.target = log_target_mixture,
    dens.ap.tar=function(y,x) 0.5*dnorm(y)+0.5*dnorm(y,mean=10,sd=1.4),
    samp.ap.tar=samp.ap.tar_mixture
  ),
  initial = 0,
  n.iter = 1000,
  gaus = FALSE # Not using Gaussian assumption
)

## ----mixture-informed-another-plot, fig.width=7, fig.height=8-----------------
par(mfrow = c(3, 1))

plot(result$samples, type = "l", main = "Default settings", 
     ylab = "y", xlab = "Iteration")
abline(h = c(0, 10), col = c("blue", "red"), lty = 2)

plot(result_informed$samples, type = "l", 
     main = "Informed normal approximate targets", 
     ylab = "y", xlab = "Iteration")
abline(h = c(0, 10), col = c("blue", "red"), lty = 2)
plot(result_informed_another$samples, type = "l", 
     main = "Informed mixture approximate target", 
     ylab = "y", xlab = "Iteration")
abline(h = c(0, 10), col = c("blue", "red"), lty = 2)

## ----mixture-importance-------------------------------------------------------
set.seed(123)
result_imp <- geomc(
  logp = list(
    log.target = log_target_mixture,
    dens.base = function(y, x) dnorm(y, mean = x, sd = 2),
    samp.base = function(x) x + 2 * rnorm(1),
    dens.ap.tar = function(y, x) 0.5 * dnorm(y) + 0.5 * dnorm(y, mean = 10, sd = 1.4),
    samp.ap.tar = samp.ap.tar_mixture
  ),
  initial = 0,
  n.iter = 1000,
  gaus = FALSE,  # Not using Gaussian assumption
  imp = list(
    enabled = TRUE,
    n.samp = 50,
    samp.base = FALSE  # Draw samples from approximate target density
  ),
  show.progress = FALSE
)

## ----importance-plot, fig.width=7, fig.height=8-------------------------------
par(mfrow = c(3, 1))

plot(result_informed$samples, type = "l", 
     main = "Informed normal approximate target", 
     ylab = "y", xlab = "Iteration")
abline(h = c(0, 10), col = c("blue", "red"), lty = 2)
plot(result_informed_another$samples, type = "l", 
     main = "Informed mixture approximate target (numerical integration)",
     ylab = "y", xlab = "Iteration")
abline(h = c(0, 10), col = c("blue", "red"), lty = 2)
plot(result_imp$samples, type = "l", 
     main = "Informed mixture approximate target (importance sampling)", 
     ylab = "y", xlab = "Iteration")
abline(h = c(0, 10), col = c("blue", "red"), lty = 2)

## ----bivariate-mixture--------------------------------------------------------
# Define log-target for bivariate mixture
log_target_mvnorm_mix <- function(x, mean1, Sigma1, mean2, Sigma2) {
  ldens1 <- geommc:::ldens_mvnorm(x, mean1, Sigma1)
  ldens2 <- geommc:::ldens_mvnorm(x, mean2, Sigma2)
  return(log(0.5 * exp(ldens1) + 0.5 * exp(ldens2)))
}

# Parameters for the two components
mean1 <- c(0, 0)
Sigma1 <- diag(2)
mean2 <- c(10, 10)
Sigma2 <- 2 * diag(2)

# Run geomc with informed approximate targets
set.seed(123)
result_biv <- geomc(
  logp = list(
    log.target = log_target_mvnorm_mix,
    mean.base = function(x) x,
    var.base = function(x) 2 * diag(2),
    mean.ap.tar = list(c(0, 0), c(10, 10)),
    var.ap.tar = list(diag(2), 2 * diag(2))
  ),
  initial = c(5, 5),  # Start between the modes
  n.iter = 1000,
  show.progress = FALSE,
  mean1 = mean1,
  Sigma1 = Sigma1,
  mean2 = mean2,
  Sigma2 = Sigma2
)

cat("Acceptance rate:", round(result_biv$acceptance.rate, 3), "\n")
cat("\nSample means:\n")
print(colMeans(result_biv$samples))
cat("\nExpected mean: (5, 5)\n")

## ----bivariate-plot, fig.width=8, fig.height=8--------------------------------
par(mfrow = c(2, 2))

# Trace plots
plot(result_biv$samples[, 1], type = "l", 
     main = "Trace plot: y1",
     ylab = expression(y[1]), xlab = "Iteration")
abline(h = 5, col = "red", lty = 2)

plot(result_biv$samples[, 2], type = "l", 
     main = "Trace plot: y2",
     ylab = expression(y[2]), xlab = "Iteration")
abline(h = 5, col = "red", lty = 2)

# Marginal densities
plot(density(result_biv$samples[, 1]), 
     main = "Marginal density: y1",
     xlab = expression(y[1]))

plot(density(result_biv$samples[, 2]), 
     main = "Marginal density: y2",
     xlab = expression(y[2]))

## ----bivariate-scatter, fig.width=7, fig.height=7-----------------------------
# Scatter plot
plot(result_biv$samples[, 1], result_biv$samples[, 2],
     pch = 16, col = rgb(0, 0, 0, 0.3),
     xlab = expression(y[1]), ylab = expression(y[2]),
     main = "Bivariate Mixture Samples")
points(mean1[1], mean1[2], col = "blue", pch = 4, cex = 3, lwd = 3)
points(mean2[1], mean2[2], col = "red", pch = 4, cex = 3, lwd = 3)
legend("topright", 
       legend = c("Mode 1", "Mode 2"),
       col = c("blue", "red"), pch = 4, lwd = 2)

## ----bivariate-mixture-diffuse------------------------------------------------
# Run geomc with diffuse g
set.seed(123)
result_biv_diffuse <- geomc(
  logp = list(
    log.target = log_target_mvnorm_mix,
    mean.ap.tar = list(c(0, 0)),
    var.ap.tar = list(400*diag(2))
  ),
  initial = c(5, 5),  # Start between the modes
  n.iter = 1000,
  mean1 = mean1,
  Sigma1 = Sigma1,
  mean2 = mean2,
  Sigma2 = Sigma2
)

cat("Acceptance rate:", round(result_biv_diffuse$acceptance.rate, 3), "\n")
cat("\nSample means:\n")
print(colMeans(result_biv_diffuse$samples))
cat("\nExpected mean: (5, 5)\n")

## ----bivariate-plot-diffuse, fig.width=8, fig.height=8------------------------
par(mfrow = c(2, 2))

# Trace plots
plot(result_biv_diffuse$samples[, 1], type = "l", 
     main = "Trace plot: y1",
     ylab = expression(y[1]), xlab = "Iteration")
abline(h = 5, col = "red", lty = 2)

plot(result_biv_diffuse$samples[, 2], type = "l", 
     main = "Trace plot: y2",
     ylab = expression(y[2]), xlab = "Iteration")
abline(h = 5, col = "red", lty = 2)

# Marginal densities
plot(density(result_biv_diffuse$samples[, 1]), 
     main = "Marginal density: y1",
     xlab = expression(y[1]))

plot(density(result_biv_diffuse$samples[, 2]), 
     main = "Marginal density: y2",
     xlab = expression(y[2]))

## ----rwm-comparison-----------------------------------------------------------
# Run random walk Metropolis
set.seed(123)
result_rwm <- geommc:::rwm(
  log_target = function(x) {
    log_target_mvnorm_mix(x, mean1, Sigma1, mean2, Sigma2)
  },
  initial = c(5, 5),
  n_iter = 1000,
  sig = 2,
  return_sample = TRUE
)

cat("Random Walk Metropolis:\n")
cat("  Acceptance rate:", round(result_rwm$acceptance_rate, 3), "\n")
cat("  Sample means:", round(colMeans(result_rwm$samples), 3), "\n\n")

cat("Geometric MCMC:\n")
cat("  Acceptance rate:", round(result_biv$acceptance.rate, 3), "\n")
cat("  Sample means:", round(colMeans(result_biv$samples), 3), "\n")

## ----rwm-plot, fig.width=7, fig.height=7--------------------------------------
par(mfrow = c(1, 2))

# Random Walk Metropolis
plot(result_rwm$samples[, 1], result_rwm$samples[, 2],
     pch = 16, col = rgb(0, 0, 0, 0.3),
     xlab = expression(y[1]), ylab = expression(y[2]),
     main = "Random Walk Metropolis",
     xlim = c(-5, 15), ylim = c(-5, 15))
points(mean1[1], mean1[2], col = "blue", pch = 4, cex = 3, lwd = 3)
points(mean2[1], mean2[2], col = "red", pch = 4, cex = 3, lwd = 3)

# Geometric MCMC
plot(result_biv$samples[, 1], result_biv$samples[, 2],
     pch = 16, col = rgb(0, 0, 0, 0.3),
     xlab = expression(y[1]), ylab = expression(y[2]),
     main = "Geometric MCMC",
     xlim = c(-5, 15), ylim = c(-5, 15))
points(mean1[1], mean1[2], col = "blue", pch = 4, cex = 3, lwd = 3)
points(mean2[1], mean2[2], col = "red", pch = 4, cex = 3, lwd = 3)

## ----normal-inference-MALA----------------------------------------------------
# Generate data
set.seed(42)
true_mu <- 5
true_sigma <- 2
n <- 100
w <- rnorm(n, mean = true_mu, sd = true_sigma)

# Define log-posterior
# Log posterior for theta = (mu, eta = log(sigma2))
log_target_normal <- function(theta, x, mu0, tau0, alpha0, beta0) {
  mu <- theta[1]
  eta <- theta[2]
  sigma2 <- exp(eta)
  
  n <- length(x)
  SSE <- sum((x - mu)^2)
  
  val <- -(n/2)*eta - SSE/(2*sigma2)
  val <- val - (mu - mu0)^2/(2*tau0^2)
  val <- val - (alpha0 + 1)*eta - beta0/sigma2
  
  return(val)
}

# Gradient of log posterior wrt (mu, eta)
grad_log_target_normal <- function(theta, x, mu0, tau0, alpha0, beta0) {
  mu <- theta[1]
  eta <- theta[2]
  sigma2 <- exp(eta)
  
  n <- length(x)
  SSE <- sum((x - mu)^2)
  
  # d/dmu
  grad_mu <- (sum(x) - n*mu)/sigma2 - (mu - mu0)/tau0^2
  
  # d/deta
  grad_eta <- -(n/2) - (alpha0 + 1)              
  grad_eta <- grad_eta + SSE/(2*sigma2) + beta0/sigma2    
  return(c(grad_mu, grad_eta))
}

# Hyperparameters set as before (weakly informative priors)
mu0 <- 0      # prior mean for mu
tau0 <- 10    # prior sd for mu
alpha0 <- 2.01  # shape parameter for inverse-gamma
beta0 <- 1.01   # scale parameter for inverse-gamma

#step-size for MALA proposal
step <- 0.001
mala_mean<- function(theta) 
  theta+(step/2)*grad_log_target_normal(theta, w, mu0, tau0, alpha0, beta0)
# Run geomc
set.seed(123)
result <- geomc(
  logp = list(log.target=log_target_normal,
              mean.base= function(theta)
                mala_mean(theta),
              var.base=function(theta) step*diag(2)),
  initial = c(0, log(1)),  # Starting at mu=0, eta= log(1)
  n.iter = 1000,
  x = w,
  mu0 = mu0,
  tau0 = tau0,
  alpha0 = alpha0,
  beta0 = beta0
)

## ----trace-plots-normal-mala, fig.width=7, fig.height=6-----------------------
par(mfrow = c(1, 2))

# Trace plots
plot(result$samples[, 1], type = "l", 
     main = expression(paste("Trace plot: ", mu)),
     ylab = expression(mu), xlab = "Iteration")
abline(v = burnin, col = "gray", lty = 2)
abline(h = true_mu, col = "red", lty = 2, lwd = 2)

plot(exp(result$samples[, 2]), type = "l", 
     main = expression(paste("Trace plot: ", sigma^2)),
     ylab = expression(sigma^2), xlab = "Iteration")
abline(v = burnin, col = "gray", lty = 2)
abline(h = true_sigma^2, col = "red", lty = 2, lwd = 2)

## ----binomial-rw--------------------------------------------------------------
# Binomial parameters
size <- 5 #m=size
true_prob <- 0.3 #p=0.3

# Define discrete random walk base density
dens.base <- function(y, x){
  if (x > 0 && x < size) {
    if (y == x - 1 || y == x + 1) {
      return(0.5)
    } else {
      return(0)
    }
  }
  if (x == 0) {
    if (y == 0 || y == 1) {
      return(0.5)
    } else {
      return(0)
    }
  }
  if (x == size) {
    if (y == size - 1 || y == size) {
      return(0.5)
    } else {
      return(0)
    }
  }
  return(0)
}

samp.base<- function(x) {
  if (x > 0 && x < size) {
    return(x + sample(c(-1, 1), 1))
  }
  if (x == 0) {
    return(sample(c(0, 1), 1))
  }
  if (x == size) {
    return(sample(c(size - 1, size), 1))
  }
} 

# Define approximate target as the discrete uniform distribution 
dens.ap.tar <- function(y, x) 1 / (size + 1)
samp.ap.tar <- function(x, kk = 1) sample(0:size, 1)

# Define the Bhattacharyya coefficients \langle \sqrt{f(\cdot|x)}, #\sqrt{g(\cdot|x)} \rangle
bhat.coef=function(x) sqrt(2 / (size + 1))

set.seed(123)
result_binom_rw <- geomc(
  logp = list(
    log.target = function(y) dbinom(y, size, true_prob, log = TRUE),
    dens.base = dens.base,
    samp.base = samp.base,
    dens.ap.tar = dens.ap.tar,
    samp.ap.tar = samp.ap.tar,
    bhat.coef=bhat.coef
  ),
  initial = 2,
  n.iter = 1000,
  gaus = FALSE,  # Not Gaussian
)

cat("Acceptance rate:", round(result_binom_rw$acceptance.rate, 3), "\n")

## ----binomial-plot-rw, fig.width=7, fig.height=5------------------------------
# Observed frequencies
obs_freq_rw <- table(factor(result_binom_rw$samples, levels = 0:size))
obs_prop_rw <- obs_freq_rw / sum(obs_freq_rw)

# True probabilities
true_prob_vals <- dbinom(0:size, size, true_prob)

# Plot comparison
barplot(rbind(obs_prop_rw, true_prob_vals),
        beside = TRUE,
        names.arg = 0:size,
        col = c("lightblue", "red"),
        main = "Binomial Distribution: Samples vs Truth",
        xlab = "Value",
        ylab = "Probability",
        legend.text = c("Sampled", "True"),
        args.legend = list(x = "topright"))

## ----binomial-table-rw--------------------------------------------------------
# Detailed comparison
comparison <- data.frame(
  Value = 0:size,
  Sampled = as.numeric(obs_prop_rw),
  True = true_prob_vals,
  Difference = as.numeric(obs_prop_rw) - true_prob_vals
)
print(round(comparison, 4))

## ----binomial-rrw-------------------------------------------------------------
# Binomial parameters
size <- 5
true_prob <- 0.3

# Define discrete reflecting random walk base density
dens.base <- function(y, x){
  if (x > 0 && x < size) {
    if (y == x - 1 || y == x + 1) {
      return(0.5)
    } else {
      return(0)
    }
  }
  if (x == 0) {
    if (y == 1) {
      return(1)
    } else {
      return(0)
    }
  }
  if (x == size) {
    if (y == size - 1) {
      return(1)
    } else {
      return(0)
    }
  }
  return(0)
}

samp.base<- function(x) {
  if (x == 0) {
    return(1)
  }
  if (x == size) {
    return(size - 1)
  }
  return(x + sample(c(-1, 1), 1))
}


# Define approximate target as the discrete uniform distribution 
dens.ap.tar <- function(y, x) 1 / (size + 1)
samp.ap.tar <- function(x, kk = 1) sample(0:size, 1)

# Define the Bhattacharyya coefficients \langle \sqrt{f(\cdot|x)}, #\sqrt{g(\cdot|x)} \rangle
bhat.coef<- function(x) {
    ifelse(x == 0 || x == size, 1, sqrt(2)) / sqrt(size + 1)
}

set.seed(123)
result_binom_rrw <- geomc(
  logp = list(
    log.target = function(y) dbinom(y, size, true_prob, log = TRUE),
    dens.base = dens.base,
    samp.base = samp.base,
    dens.ap.tar = dens.ap.tar,
    samp.ap.tar = samp.ap.tar,
    bhat.coef = bhat.coef
  ),
  initial = 2,
  n.iter = 1000,
  gaus = FALSE,  # Not Gaussian
)

cat("Acceptance rate:", round(result_binom_rrw$acceptance.rate, 3), "\n")

## ----binomial-plot-rrw, fig.width=7, fig.height=5-----------------------------
# Observed frequencies
obs_freq_rrw <- table(factor(result_binom_rrw$samples, levels = 0:size))
obs_prop_rrw <- obs_freq_rrw / sum(obs_freq_rrw)

# True probabilities
true_prob_vals <- dbinom(0:size, size, true_prob)

# Plot comparison
barplot(rbind(obs_prop_rrw, true_prob_vals),
        beside = TRUE,
        names.arg = 0:size,
        col = c("lightblue", "red"),
        main = "Binomial Distribution: Samples vs Truth",
        xlab = "Value",
        ylab = "Probability",
        legend.text = c("Sampled", "True"),
        args.legend = list(x = "topright"))

## ----binomial-table-rrw-------------------------------------------------------
# Detailed comparison
comparison <- data.frame(
  Value = 0:size,
  Sampled = as.numeric(obs_prop_rrw),
  True = true_prob_vals,
  Difference = as.numeric(obs_prop_rrw) - true_prob_vals
)
print(round(comparison, 4))

## ----binomial-----------------------------------------------------------------
# Binomial parameters
size <- 5
true_prob <- 0.3

# Define discrete uniform base density
dens.base <- function(y, x) 1 / (size + 1)
samp.base <- function(x) sample(0:size, 1)

# Define approximate target (another binomial with different probability)
dens.ap.tar <- function(y, x) dbinom(y, size, 0.7) #p_1=0.7
samp.ap.tar <- function(x, kk = 1) rbinom(1, size, 0.7)
set.seed(123)
result_binom <- geomc(
  logp = list(
    log.target = function(y) dbinom(y, size, true_prob, log = TRUE),
    dens.base = dens.base,
    samp.base = samp.base,
    dens.ap.tar = dens.ap.tar,
    samp.ap.tar = samp.ap.tar
  ),
  initial = 2,
  n.iter = 1000,
  ind = TRUE,  # Base and approximate target don't depend on current state
  gaus = FALSE,  # Not Gaussian
  imp = list(
    enabled = TRUE,
    n.samp = 1000,
    samp.base = TRUE #samples from the base density is used in the  importance sampling 
  )
)

cat("Acceptance rate:", round(result_binom$acceptance.rate, 3), "\n")

## ----binomial-plot, fig.width=7, fig.height=5---------------------------------
# Observed frequencies
obs_freq<- table(factor(result_binom$samples, levels = 0:size))
obs_prop <- obs_freq / sum(obs_freq)

# True probabilities
true_prob_vals <- dbinom(0:size, size, true_prob)

# Plot comparison
barplot(rbind(obs_prop, true_prob_vals),
        beside = TRUE,
        names.arg = 0:size,
        col = c("lightblue", "red"),
        main = "Binomial Distribution: Samples vs Truth",
        xlab = "Value",
        ylab = "Probability",
        legend.text = c("Sampled", "True"),
        args.legend = list(x = "topright"))

## ----binomial-table-----------------------------------------------------------
# Detailed comparison
comparison <- data.frame(
  Value = 0:size,
  Sampled = as.numeric(obs_prop),
  True = true_prob_vals,
  Difference = as.numeric(obs_prop) - true_prob_vals
)
print(round(comparison, 4))

