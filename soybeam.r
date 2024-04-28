# load data
setwd("C:/Users/Mayin/Documents/1GRADUATE/1. Study/2. 24Spring/5224 Bayesian Statistics/5224_Project/Bayesian-Statistics-Project")
data <- read.csv("data.csv", header = TRUE)
colnames(dat)

qqnorm(dat$GY)
qqline(dat$GY)

setwd("/Users/yuhanliu/Desktop")
# Load necessary libraries
library(dplyr)

fit <- lm(GY ~  MHG, data = group)

## OLS

# Prepare an empty list for storing models by group and a data frame for estimates
group_list <- list()
ols_estimate <- data.frame()
res_var <- list()

# Loop over unique cultivars
for (cultivar in unique(data$Cultivar)) {
  group <- data[data$Cultivar == cultivar,]
  group_list[[cultivar]] <- group
  
  # Fit a linear model using all predictors
  fit <- lm(GY ~ Season + PH + IFP + NLP + NGP + NGL + NS + MHG, data = group)
  res_var[[cultivar]] <- summary(fit)$sigma^2

  # Extract coefficients and bind them to the ols_estimate data frame
  ols_estimate <- rbind(ols_estimate, data.frame(
    cultivar = cultivar,
    beta_0 = coef(fit)[[1]],  # Intercept
    beta_Season = coef(fit)[["Season"]],
    beta_PH = coef(fit)[["PH"]],
    beta_IFP = coef(fit)[["IFP"]],
    beta_NLP = coef(fit)[["NLP"]],
    beta_NGP = coef(fit)[["NGP"]],
    beta_NGL = coef(fit)[["NGL"]],
    beta_NS = coef(fit)[["NS"]],
    beta_MHG = coef(fit)[["MHG"]]
  ))
}

# Display the estimates
ols_estimate
