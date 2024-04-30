# load data
setwd("C:/Users/Mayin/Documents/1GRADUATE/1. Study/2. 24Spring/5224 Bayesian Statistics/5224_Project/Bayesian-Statistics-Project")
source("slice_pred_plot.r")
data_ori <- read.csv("data.csv", header = TRUE)
n <- nrow(data)

# Define a function for min-max normalization
min_max_normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

# Apply the function to each numeric column in the data frame
data<- data.frame(lapply(data_ori, function(x) {
  if (is.numeric(x)) min_max_normalize(x) else x
}))
# cor of NGP and MHG
cor(data$NGP, data$MHG)
########################
# OLS
########################
# Prepare an empty list for storing models by group and a data frame for estimates
group_list <- list()
ols_estimate <- data.frame()
res_var <- list()

# Loop over unique cultivars
for (cultivar in unique(data$Cultivar)) {
  group <- data[data$Cultivar == cultivar,]
  group_list[[cultivar]] <- group
  
  # Fit a linear model using all predictors
  fit <- lm(GY ~ NGP + MHG + NGP * MHG + I(MHG^2), data = group)
  res_var[[cultivar]] <- summary(fit)$sigma^2

  # Extract coefficients and bind them to the ols_estimate data frame
  ols_estimate <- rbind(ols_estimate, data.frame(
    cultivar = cultivar,
    beta_0 = coef(fit)[[1]],  # Intercept
    beta_NGP = coef(fit)[["NGP"]],
    beta_MHG = coef(fit)[["MHG"]],
    beta_NGP_MHG = coef(fit)[["NGP:MHG"]],
    beta_MHG2 = coef(fit)[["I(MHG^2)"]]
  ))
}

# Display the estimates
ols_estimate
theta  <- colMeans(ols_estimate[ , 2:6])
Sigma <- cov(ols_estimate[ , 2:6])
sigma2 <- mean(unlist(res_var))

theta
Sigma
sigma2

# plot OLS 
slice_pred_plot(data, ols_estimate, slice = "NGP", method = "OLS")
# ggsave("my_plot.png", slice_pred_plot(data, ols_estimate, slice = "NGP", method = "OLS"))
slice_pred_plot(data, ols_estimate, slice = "MHG", method = "OLS")
ggsave("2.3OLS_M.png", slice_pred_plot(data, ols_estimate, slice = "MHG", method = "OLS"))

########################
# Bayesian
########################
library(MASS)
library(MCMCpack)

#============================
# weak infromative prior
#============================

theta  <- colMeans(ols_estimate[ , 2:6])
Sigma <- cov(ols_estimate[ , 2:6])
sigma2 <- mean(unlist(res_var))

# number of groups
m <- 40
nj <- 8
p <- 5

## prior parameters
mu0_w <- c(rep(0, p))
Gamma0_w <- Sigma * 10
S0_w <- Sigma / 100
eta0_w <- p
nu0_w <- 10
sigma02_w <- sigma2

## starting values
mu_w <- c(rep(0, p))
S_w <- diag(p)
beta1_w <- ols_estimate[ , "beta_0"]
beta2_w <- ols_estimate[ , "beta_NGP"]
beta3_w <- ols_estimate[ , "beta_MHG"]
beta4_w <- ols_estimate[ , "beta_NGP_MHG"]
beta5_w <- ols_estimate[ , "beta_MHG2"]

## Gibbs sampler: store
NN <- 5000
burnin <- 500

MU_w <- matrix(0, p, NN)
SIGMA_w <- matrix(0, p^2, NN)
S2_w <- NULL
BETA1_w <- matrix(0, m, NN)
BETA2_w <- matrix(0, m, NN)
BETA3_w <- matrix(0, m, NN)
BETA4_w <- matrix(0, m, NN)
BETA5_w <- matrix(0, m, NN)

set.seed(5224)

for(s in 1:(NN+burnin)) {
  ##calculate the SSR
  ssr_w <- numeric(length(m))
  for (i in 1:m) {
    group <- group_list[[i]]
    predicted <- beta1_w[i] + beta2_w[i] * group$NGP + beta3_w[i] * group$MHG + beta4_w[i] * group$NGP * group$MHG + beta5_w[i] * group$MHG^2
    residuals <- group$GY - predicted
    ssr_w[i] <- sum(residuals^2)
  }

  ##update s2
  s2_w <- 1 / rgamma(1, (nu0_w + m * nj )/2,
               (nu0_w *sigma02_w +sum(ssr_w)) / 2)

  ##update mu (theta)
  col_mean_beta <- colMeans(cbind(beta1_w, beta2_w, beta3_w, beta4_w, beta5_w))
  var.mu<-  solve(solve(Gamma0_w) + m * solve(S_w))
  mean.mu<- var.mu %*% (solve(Gamma0_w) %*% mu0_w + m * solve(S_w) %*% col_mean_beta)
  mu_w <- mvrnorm(1, mean.mu, var.mu)
  
  ##update Sigma (S)
  S_theta <- matrix(0, p, p)
  for (j in 1:m){
    beta_j <- c(beta1_w[j], beta2_w[j], beta3_w[j], beta4_w[j], beta5_w[j])
    S_theta <- S_theta + (beta_j - mu_w) %*% t(beta_j - mu_w)
  }
  var.S <- solve(S0_w + S_theta)
  S_w <- solve(rwish(eta0_w + m, var.S))

  ## update beta
  lambda <- 0
  for (j in 1:m){
    group <- group_list[[j]]
    X <- cbind(1, group$NGP, group$MHG, group$NGP * group$MHG, group$MHG^2)
    var.beta<-  solve(solve(S_w + lambda * diag(nrow(S_w))) + t(X) %*% X / s2_w + lambda * diag(nrow(S_w)))
    mean.beta<- var.beta %*% (solve(S_w + lambda * diag(nrow(S_w))) %*% mu_w + t(X) %*% group$GY / s2_w)
    beta <- mvrnorm(1, mean.beta, var.beta)
    if(s > burnin){
      BETA1_w[j,s-burnin] <- beta[1]
      BETA2_w[j,s-burnin] <- beta[2]
      BETA3_w[j,s-burnin] <- beta[3]
      BETA4_w[j,s-burnin] <- beta[4]
      BETA5_w[j,s-burnin] <- beta[5]
    }
  }

  # save the parameter values
  if(s > burnin){
    MU_w[,s-burnin] <- mu_w
    SIGMA_w[,s-burnin] <- as.vector(S_w)
    S2_w<-c(S2_w,s2_w) 
  }
}       

bayes_estimate_w <- data.frame(
  cultivar = unique(data$Cultivar)[1:40],
  beta_0 = apply(BETA1_w, 1, mean),
  beta_NGP = apply(BETA2_w, 1, mean),
  beta_MHG = apply(BETA3_w, 1, mean),
  beta_NGP_MHG = apply(BETA4_w, 1, mean),
  beta_MHG2 = apply(BETA5_w, 1, mean)
)

# display bayes_estimate_w
bayes_estimate_w

# plot for bayes_estimate_w
slice_pred_plot(data, bayes_estimate_w, slice = "NGP", method = "Bayes(Wwak Informative Prior)")
slice_pred_plot(data, bayes_estimate_w, slice = "MHG", method = "Bayes(Wwak Informative Prior)")

#-- Diagnose
## trace plot & ESS
# theta
ess_mu_w <- numeric(p)
par(mfrow=c(2,3))
for (i in 1:5){
  plot(MU_w[i,], type="l", ylab=paste("mu_",i), xlab="iteration")
  ess_mu_w[i] <- effectiveSize(MU_w[i,])
}
ess_mu_w
# sigma
par(mfrow=c(1,1))
plot(S2_w, type="l", ylab="sigma2", xlab="iteration")
effectiveSize(S2_w)
# Sigma
ess_Sigma_w <- numeric(p^2)
par(mfrow=c(5,5))
for (i in 1:25){
  plot(SIGMA_w[i,], type="l", ylab=paste("Sigma_",i), xlab="iteration")
  ess_Sigma_w[i] <- effectiveSize(SIGMA_w[i,])
}
ess_Sigma_w
# beta (shoule have 40 plots, just show the first 3)
par(mfrow=c(5,3))
for(i in 1:3){
  plot(BETA1_w[i,], type="l", ylab=paste("beta1_",i), xlab="iteration")
  plot(BETA2_w[i,], type="l", ylab=paste("beta2_",i), xlab="iteration")
  plot(BETA3_w[i,], type="l", ylab=paste("beta3_",i), xlab="iteration")
  plot(BETA4_w[i,], type="l", ylab=paste("beta4_",i), xlab="iteration")
  plot(BETA5_w[i,], type="l", ylab=paste("beta5_",i), xlab="iteration")
}
ess_beta1_w <- numeric(m)
ess_beta2_w <- numeric(m)
ess_beta3_w <- numeric(m)
ess_beta4_w <- numeric(m)
ess_beta5_w <- numeric(m)
for (i in 1:m){
  ess_beta1_w[i] <- effectiveSize(BETA1_w[i,])
  ess_beta2_w[i] <- effectiveSize(BETA2_w[i,])
  ess_beta3_w[i] <- effectiveSize(BETA3_w[i,])
  ess_beta4_w[i] <- effectiveSize(BETA4_w[i,])
  ess_beta5_w[i] <- effectiveSize(BETA5_w[i,])
}
ess_beta_w_df <- data.frame(
  ess_beta1_w = ess_beta1_w,
  ess_beta2_w = ess_beta2_w,
  ess_beta3_w = ess_beta3_w,
  ess_beta4_w = ess_beta4_w,
  ess_beta5_w = ess_beta5_w
)
ess_beta_w_df


#============================
# infromative prior
#============================
## prior parameters
mu0_i <- theta
Gamma0_i <- S0_i <- Sigma
eta0_i <- 2*p
nu0_i <- 2
sigma02_i <- sigma2

# initial points
mu_i <- theta
S_i <- Sigma
beta1_i <- ols_estimate[ , "beta_0"]
beta2_i <- ols_estimate[ , "beta_NGP"]
beta3_i <- ols_estimate[ , "beta_MHG"]
beta4_i <- ols_estimate[ , "beta_NGP_MHG"]
beta5_i <- ols_estimate[ , "beta_MHG2"]

## Gibbs sampler: store
NN <- 5000
burnin <- 500

MU_i <- matrix(0, p, NN)
SIGMA_i <- matrix(0, p^2, NN)
S2_i <- NULL
BETA1_i <- matrix(0, m, NN)
BETA2_i <- matrix(0, m, NN)
BETA3_i <- matrix(0, m, NN)
BETA4_i <- matrix(0, m, NN)
BETA5_i <- matrix(0, m, NN)

set.seed(5224)

for(s in 1:(NN+burnin)) {
  ##calculate the SSR
  ssr_i <- numeric(length(m))
  for (i in 1:m) {
    group <- group_list[[i]]
    predicted <- beta1_i[i] + beta2_i[i] * group$NGP + beta3_i[i] * group$MHG + beta4_i[i] * group$NGP * group$MHG + beta5_i[i] * group$MHG^2
    residuals <- group$GY - predicted
    ssr_i[i] <- sum(residuals^2)
  }

  ##update s2
  s2_i <- 1 / rgamma(1, (nu0_i + m * nj )/2,
               (nu0_i *sigma02_i +sum(ssr_i)) / 2)

  ##update mu (theta)
  col_mean_beta <- colMeans(cbind(beta1_i, beta2_i, beta3_i, beta4_i, beta5_i))
  var.mu<-  solve(solve(Gamma0_i) + m * solve(S_i))
  mean.mu<- var.mu %*% (solve(Gamma0_i) %*% mu0_i + m * solve(S_i) %*% col_mean_beta)
  mu_i <- mvrnorm(1, mean.mu, var.mu)
  
  ##update Sigma (S)
  S_theta <- matrix(0, p, p)
  for (j in 1:m){
    beta_j <- c(beta1_i[j], beta2_i[j], beta3_i[j], beta4_i[j], beta5_i[j])
    S_theta <- S_theta + (beta_j - mu_i) %*% t(beta_j - mu_i)
  }
  var.S <- solve(S0_i + S_theta)
  S_i <- solve(rwish(eta0_i + m, var.S))

  ## update beta
  lambda <- 1e-6
  for (j in 1:m){
    group <- group_list[[j]]
    X <- cbind(1, group$NGP, group$MHG, group$NGP * group$MHG, group$MHG^2)
    var.beta<-  solve(solve(S_i + lambda * diag(nrow(S_i))) + t(X) %*% X / s2_i + lambda * diag(nrow(S_i)))
    mean.beta<- var.beta %*% (solve(S_i + lambda * diag(nrow(S_i))) %*% mu_i + t(X) %*% group$GY / s2_i)
    beta <- mvrnorm(1, mean.beta, var.beta)
    if(s > burnin){
      BETA1_i[j,s-burnin] <- beta[1]
      BETA2_i[j,s-burnin] <- beta[2]
      BETA3_i[j,s-burnin] <- beta[3]
      BETA4_i[j,s-burnin] <- beta[4]
      BETA5_i[j,s-burnin] <- beta[5]
    }
  }

  # save the parameter values
  if(s > burnin){
    MU_i[,s-burnin] <- mu_i
    SIGMA_i[,s-burnin] <- as.vector(S_i)
    S2_i<-c(S2_i,s2_i) 
  }
}       

bayes_estimate_i <- data.frame(
  cultivar = unique(data$Cultivar)[1:40],
  beta_0 = apply(BETA1_i, 1, mean),
  beta_NGP = apply(BETA2_i, 1, mean),
  beta_MHG = apply(BETA3_i, 1, mean),
  beta_NGP_MHG = apply(BETA4_i, 1, mean),
  beta_MHG2 = apply(BETA5_i, 1, mean)
)

# display bayes_estimate_i
bayes_estimate_i

# plot for bayes_estimate_i
slice_pred_plot(data, bayes_estimate_i, slice = "NGP", method = "Bayes(Wwak Informative Prior)")
slice_pred_plot(data, bayes_estimate_i, slice = "MHG", method = "Bayes(Wwak Informative Prior)")

#-- Diagnose
## trace plot & ESS
# theta
ess_mu_i <- numeric(p)
par(mfrow=c(2,3))
for (i in 1:5){
  plot(MU_i[i,], type="l", ylab=paste("mu_",i), xlab="iteration")
  ess_mu_i[i] <- effectiveSize(MU_i[i,])
}
ess_mu_i
# sigma
par(mfrow=c(1,1))
plot(S2_i, type="l", ylab="sigma2", xlab="iteration")
effectiveSize(S2_i)
# Sigma
ess_Sigma_i <- numeric(p^2)
par(mfrow=c(5,5))
for (i in 1:25){
  plot(SIGMA_i[i,], type="l", ylab=paste("Sigma_",i), xlab="iteration")
  ess_Sigma_i[i] <- effectiveSize(SIGMA_i[i,])
}
ess_Sigma_i
# beta (shoule have 40 plots, just show the first 3)
par(mfrow=c(5,3))
for(i in 1:3){
  plot(BETA1_i[i,], type="l", ylab=paste("beta1_",i), xlab="iteration")
  plot(BETA2_i[i,], type="l", ylab=paste("beta2_",i), xlab="iteration")
  plot(BETA3_i[i,], type="l", ylab=paste("beta3_",i), xlab="iteration")
  plot(BETA4_i[i,], type="l", ylab=paste("beta4_",i), xlab="iteration")
  plot(BETA5_i[i,], type="l", ylab=paste("beta5_",i), xlab="iteration")
}
ess_beta1_i <- numeric(m)
ess_beta2_i <- numeric(m)
ess_beta3_i <- numeric(m)
ess_beta4_i <- numeric(m)
ess_beta5_i <- numeric(m)
for (i in 1:m){
  ess_beta1_i[i] <- effectiveSize(BETA1_i[i,])
  ess_beta2_i[i] <- effectiveSize(BETA2_i[i,])
  ess_beta3_i[i] <- effectiveSize(BETA3_i[i,])
  ess_beta4_i[i] <- effectiveSize(BETA4_i[i,])
  ess_beta5_i[i] <- effectiveSize(BETA5_i[i,])
}
ess_beta_i_df <- data.frame(
  ess_beta1_i = ess_beta1_i,
  ess_beta2_i = ess_beta2_i,
  ess_beta3_i = ess_beta3_i,
  ess_beta4_i = ess_beta4_i,
  ess_beta5_i = ess_beta5_i
)
ess_beta_i_df


#####################################
# Analysis: Result from different Prior
####################################
library(ggplot2)
library(gridExtra)
# Comparison of prior and posterior: theta
## weak informative prior
theta_ci_w <- data.frame() # 95 CI
theta_prior_w <- mvrnorm(NN, mu0_w, Gamma0_w)
theta_w_plots <- list()
for (i in 1:p) {
  theta_prior_w_i <- theta_prior_w[, i]
  theta_post_w_i <- MU_w[i, ]
  pp <- ggplot() +
    geom_density(aes(x = theta_prior_w_i), fill = "blue", alpha = 0.5) +
    geom_density(aes(x = theta_post_w_i), fill = "red", alpha = 0.5) +
    labs(x = paste0("theta_", i), y = "Density", title = paste0(i, "th Element of Theta"))
  theta_w_plots[[i]] <- pp

  ci <- quantile(theta_post_w_i, probs = c(0.025, 0.975))
  theta_ci_w <- rbind(theta_ci_w, data.frame(
    row = i,
    lb = ci[1],
    ub = ci[2]
  ))
}
theta_ci_w
grid.arrange(grobs = theta_w_plots, nrow=5)
# ggsave("6.1theta_w.png", grid.arrange(grobs = theta_w_plots, nrow=5))
## informative prior
theta_ci_i <- data.frame() # 95 CI
theta_prior_i <- mvrnorm(NN, mu0_i, Gamma0_i)
theta_i_plots <- list()
for (i in 1:5) {
  theta_prior_i_i <- theta_prior_i[, i]
  theta_post_i_i <- MU_i[i, ]
  pp <- ggplot() +
    geom_density(aes(x = theta_prior_i_i), fill = "blue", alpha = 0.5) +
    geom_density(aes(x = theta_post_i_i), fill = "red", alpha = 0.5) +
    labs(x = paste0("theta_", i), y = "Density", title = paste0(i, "th Element of Theta"))
  theta_i_plots[[i]] <- pp

  ci <- quantile(theta_post_i_i, probs = c(0.025, 0.975))
  theta_ci_i <- rbind(theta_ci_i, data.frame(
    row = i,
    lb = ci[1],
    ub = ci[2]
  ))
}
theta_ci_i
grid.arrange(grobs = theta_i_plots, nrow=5)
ggsave("6.2theta_i.png", grid.arrange(grobs = theta_i_plots, nrow=5))

# Comparison of prior and posterior: sigma2
sigma2_prior_w <- 1/ rgamma(NN, nu0_w/2, nu0_w * sigma02_w / 2)
sigma2_prior_i <- 1/ rgamma(NN, nu0_i/2, nu0_i * sigma02_i / 2)
p1 <- ggplot() +
  geom_density(aes(x = sigma2_prior_w), fill = "blue", alpha = 0.5) +
  geom_density(aes(x = S2_w), fill = "red", alpha = 0.5) +
  xlim(min(S2_w), max(S2_w)+0.001) + 
  labs(x = "sigma2", y = "Density", title = "weak informative prior: Sigma2")
p2 <- ggplot() +
  geom_density(aes(x = sigma2_prior_i), fill = "blue", alpha = 0.5) +
  geom_density(aes(x = S2_i), fill = "red", alpha = 0.5) +
  xlim(min(S2_i), max(S2_i)+0.001) +
  labs(x = "sigma2", y = "Density", title = "informative prior: Sigma2")
grid.arrange(p1, p2, nrow=2)

# Comparison of prior and posterior: Beta
############# CI overlap
beta_ci_plot <- function(beta_w, beta_i, title){
  ci_beta_w <- data.frame()
  ci_beta_i <- data.frame()
  for (i in 1:m){
    ci_w <- c(quantile(beta_w[i, ], 0.025), quantile(beta_w[i, ], 0.975))
    ci_i <- c(quantile(beta_i[i, ], 0.025), quantile(beta_i[i, ], 0.975))
    ci_beta_w <- rbind(ci_beta_w, data.frame(
      group = i,
      lb = ci_w[1],
      ub = ci_w[2]
    ))
    ci_beta_i <- rbind(ci_beta_i, data.frame(
      group = i,
      lb = ci_i[1],
      ub = ci_i[2]
    ))}
  ggplot() +
    geom_rect(data = ci_beta_w, aes(xmin = group+0.05 , xmax = group + 0.95, ymin = lb, ymax = ub), fill = "blue", alpha = 0.7) +
    geom_rect(data = ci_beta_i, aes(xmin = group+0.05, xmax = group + 0.95, ymin = lb, ymax = ub), fill = "yellow", alpha = 0.6) +
    geom_hline(yintercept = mean(beta_w), linetype = "dashed", color = "blue", size = 1) +
    geom_hline(yintercept = mean(beta_i), linetype = "dashed", color = "yellow", size = 1) +
    labs(x = "Group", y = "Interval") +
    ggtitle(title) +
    theme_minimal()
}

beta_ci_plot(BETA1_w, BETA1_i, "Comparison 95% Posterior Interval of Beta_1")
beta_ci_plot(BETA2_w, BETA2_i, "Comparison 95% Posterior Interval of Beta_2")
beta_ci_plot(BETA3_w, BETA3_i, "Comparison 95% Posterior Interval of Beta_3")
beta_ci_plot(BETA4_w, BETA4_i, "Comparison 95% Posterior Interval of Beta_4")
beta_ci_plot(BETA5_w, BETA5_i, "Comparison 95% Posterior Interval of Beta_5")
## 结论：就算prior不一样，结果类似： vary

# Comparison of prior and posterior: SIGMA
## weak informative prior
Sigma_prior <- NULL
SSS <- solve(var(ols_estimate[,2:6]))
for (i in 1:NN) {
  sample <- solve(rwish(eta0_w, solve(S0_i)))
  Sigma_prior <- rbind(Sigma_prior, c(sample))
}
## posterior
Sigma_post <- SIGMA_w

## plot + CI
SIGMA_w_ci <- data.frame()
plots <- list()
p <- 5
for (i in 1:p){
  for (j in 1:p){
    index <- p*(i-1)+j
    sigma_prior <- Sigma_prior[, index]
    sigma_post <- Sigma_post[index, ]
    
    df <- data.frame(value = c(sigma_prior, sigma_post),
                     group = rep(c("Prior", "Post"), each = length(sigma_prior)))
    
    pp <- ggplot(df, aes(x = value, fill = group)) +
      geom_density(alpha = 0.4) +
      labs(x = paste("sigma_",i,j), y = "Density", 
           title = paste(i, j, "Element of Sigma")) +
      xlim(min(sigma_post), max(sigma_post)) +
      theme_minimal() +
      theme(legend.position = "none")
    
    plots[[index]] <- pp

    ci_w <- c(quantile(sigma_post, 0.025), quantile(sigma_post, 0.975))
    SIGMA_w_ci <- rbind(SIGMA_w_ci, data.frame(
      row = index,
      lb = ci_w[1],
      ub = ci_w[2]
    ))
  }
}
grid.arrange(grobs = plots, nrow = 5)
SIGMA_w_ci
## informative
Sigma_prior_i <- NULL
for (i in 1:NN) {
  sample <- solve(rwish(eta0_i, solve(S0_i)))
  Sigma_prior_i <- rbind(Sigma_prior_i, c(sample))
}
## posterior
Sigma_post_i <- SIGMA_i
## plot
plots <- list()
p <- 5
for (i in 1:p){
  for (j in 1:p){
    index <- p*(i-1)+j
    sigma_prior_i <- Sigma_prior_i[, index]
    sigma_post_i <- Sigma_post_i[index, ]
    
    df <- data.frame(value = c(sigma_prior_i, sigma_post_i),
                     group = rep(c("Prior", "Post"), each = length(sigma_prior_i)))
    
    pp <- ggplot(df, aes(x = value, fill = group)) +
      geom_density(alpha = 0.4) +
      labs(x = paste("sigma_",i,j), y = "Density", 
           title = paste(i, j, "Element of Sigma")) +
      xlim(min(sigma_post_i), max(sigma_post_i)) +
      theme_minimal() +
      theme(legend.position = "none")
    plots[[index]] <- pp
  }
}
grid.arrange(grobs = plots, nrow = 5)


#####################################
# Relation between NGP and MHG (over randam cultivar)
####################################
## grid
NGP_ori_range <- range(data_ori$NGP)
MHG_ori_range <- range(data_ori$MHG)
GY_ori_range <- range(data_ori$GY)

NGP_grid <- seq(0, 1, 0.01)
MHG_grid <- seq(0, 1, 0.01)

# weak informative prior
pred_w <- data.frame()
for (j in 1:100){
  print(j)
  for (k in 1:100){
    for (i in 1:m){
      beta_w <- bayes_estimate_w[i, ]
      pred_each <- NULL
      pred <- beta_w$beta_0 + beta_w$beta_NGP * NGP_grid[j] + beta_w$beta_MHG * MHG_grid[k] + beta_w$beta_NGP_MHG * NGP_grid[j] * MHG_grid[k] + beta_w$beta_MHG2 * MHG_grid[k]^2
      pred_each <- c(pred_each, pred)
    }
    pred_w <- rbind(pred_w, data.frame(NGP = NGP_grid[j], MHG = MHG_grid[k], predict_mean = mean(pred_each)))
  }
}
## de-normalize
pred_w$NGP_ori <- pred_w$NGP * (NGP_ori_range[2] - NGP_ori_range[1]) + NGP_ori_range[1]
pred_w$MHG_ori <- pred_w$MHG * (MHG_ori_range[2] - MHG_ori_range[1]) + MHG_ori_range[1]
pred_w$predict_mean_ori <- pred_w$predict_mean * (GY_ori_range[2] - GY_ori_range[1]) + GY_ori_range[1]
pred_w
## 3D plot
library(rgl)
par(mfrow = c(1, 1))
blue_palette <- colorRampPalette(c("lightblue", "darkblue")) # color palette
with(pred_w, plot3d(NGP, MHG, predict_mean, xlab = "NGP", ylab = "MHG", zlab = "Predict Mean", col = blue_palette(length(predict_mean)), alpha = 0.5))
## plot: contour + heatmap
ggplot(pred_w, aes(x = NGP_ori, y = MHG_ori, z = predict_mean_ori, fill = predict_mean_ori)) +
  geom_tile() +
  geom_contour(aes(colour = after_stat(level)), bins = 10) +
  scale_fill_gradientn(colors = blue_palette(100)) +
  labs(x = "NGP", y = "MHG", fill = "Predict Mean") +
  ggtitle("Predict Mean over NGP and MHG (weak informative prior)") +
  theme_minimal()
## Find the row with the maximum prediction
pred_w[which.max(pred_w$predict_mean),]
## 95 CI
max_pred_w_ci <- data.frame()
# NGP_max <- pred_w[which.max(pred_w$predict_mean), "NGP"]
NGP_max <- 1
MHG_max <- 0
for (i in 1:m){
  beta_w <- bayes_estimate_w[i, ]
  max_pred <- beta_w$beta_0 + beta_w$beta_NGP * NGP_max + beta_w$beta_MHG * MHG_max + beta_w$beta_NGP_MHG * NGP_max * MHG_max + beta_w$beta_MHG2 * MHG_max^2
  max_pred_ori <- max_pred * (GY_ori_range[2] - GY_ori_range[1]) + GY_ori_range[1]
  max_pred_w_ci <- rbind(max_pred_w_ci, data.frame(
    group = i,
    predict = max_pred_ori
  ))
}
quantile(max_pred_w_ci$predict, probs = c(0.4, 0.6))
mean(max_pred_w_ci$predict)
ggplot()+
  geom_density(data = max_pred_w_ci, aes(x = predict), fill = "blue", alpha = 0.5) +
  labs(x = "Predict", y = "Density", title = "Predict Density at Maximum Point (weak informative prior)") +
  theme_minimal()
# informative prior
pred_i <- data.frame()
for (j in 1:100){
  print(j)
  for (k in 1:100){
    for (i in 1:m){
      beta_i <- bayes_estimate_i[i, ]
      pred_each <- NULL
      pred <- beta_i$beta_0 + beta_i$beta_NGP * NGP_grid[j] + beta_i$beta_MHG * MHG_grid[k] + beta_i$beta_NGP_MHG * NGP_grid[j] * MHG_grid[k] + beta_i$beta_MHG2 * MHG_grid[k]^2
      pred_each <- c(pred_each, pred)
    }
    pred_i <- rbind(pred_i, data.frame(NGP = NGP_grid[j], MHG = MHG_grid[k], predict_mean = mean(pred_each)))
  }
}
## de-normalize
pred_i$NGP_ori <- pred_i$NGP * (NGP_ori_range[2] - NGP_ori_range[1]) + NGP_ori_range[1]
pred_i$MHG_ori <- pred_i$MHG * (MHG_ori_range[2] - MHG_ori_range[1]) + MHG_ori_range[1]
pred_i$predict_mean_ori <- pred_i$predict_mean * (GY_ori_range[2] - GY_ori_range[1]) + GY_ori_range[1]
pred_i
## 3D plot
with(pred_i, plot3d(NGP, MHG, predict_mean, xlab = "NGP", ylab = "MHG", zlab = "Predict Mean", col = blue_palette(length(predict_mean)), alpha = 0.5))
## plot: contour + heatmap
ggplot(pred_i, aes(x = NGP_ori, y = MHG_ori, z = predict_mean_ori, fill = predict_mean_ori)) +
  geom_tile() +
  geom_contour(aes(colour = after_stat(level)), bins = 10) +
  scale_fill_gradientn(colors = blue_palette(100)) +
  labs(x = "NGP", y = "MHG", fill = "Predict Mean") +
  ggtitle("Predict Mean over NGP and MHG (informative prior)") +
  theme_minimal()
## Find the row with the maximum prediction
pred_i[which.max(pred_i$predict_mean),]
