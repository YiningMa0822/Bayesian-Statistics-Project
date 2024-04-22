#--------------------------------------------------------Yining's Laptop
setwd("C:/Users/Mayin/Documents/1GRADUATE/1. Study/2. 24Spring/5224 Bayesian Statistics/5224_Project/Bayesian-Statistics-Project/")


# ------------------------------------------------------ Data wrangling
data <- read.csv('student-por.csv')
data['sum'] <- rowSums(data[, c('G1', 'G2', 'G3')])
colnames(data)
school_binary <- ifelse(data$school == "GP", 1, 0)
sex_binary <- ifelse(data$sex == "F", 1, 0)
address_binary <- ifelse(data$address == "U", 1, 0)
famsize_binary <- ifelse(data$famsize == "GT3", 1, 0)
Pstatus_binary <- ifelse(data$Pstatus == "T", 1, 0)
Mjob_onehot <- model.matrix(~ Mjob - 1, data)
Fjob_onehot <- model.matrix(~ Fjob - 1, data)
reason_onehot <- model.matrix(~ reason - 1, data)
guardian_onehot <- model.matrix(~ guardian - 1, data)
schoolsup_binary <- ifelse(data$schoolsup == "yes", 1, 0)
famsup_binary <- ifelse(data$famsup == "yes", 1, 0)
paid_binary <- ifelse(data$paid == "yes", 1, 0)
activities_binary <- ifelse(data$activities == "yes", 1, 0)
nursery_binary <- ifelse(data$nursery == "yes", 1, 0)
higher_binary <- ifelse(data$higher == "yes", 1, 0)
internet_binary <- ifelse(data$internet == "yes", 1, 0)
romantic_binary <- ifelse(data$romantic == "yes", 1, 0)

encode_data <- cbind(school_binary, sex_binary, data$age, address_binary, 
famsize_binary, Pstatus_binary, data$Medu, data$Fedu, Mjob_onehot, Fjob_onehot, 
reason_onehot, guardian_onehot, data$traveltime, data$studytime, data$failures, 
schoolsup_binary, famsize_binary, paid_binary, activities_binary, nursery_binary, 
higher_binary, internet_binary, romantic_binary, data$famrel, data$freetime, 
data$goout, data$Dalc, data$Walc, data$health, data$absences)

encode_data <- apply(encode_data, 2, function(col) as.numeric(as.character(col)))
dim(encode_data)
# ------------------------------------------------------ Bayesian Factor Analysis
## Determine the number of factors
cov_encode_data <- cov(encode_data)

# Eigenvalues
eigen(cov_encode_data)$values
# Load the necessary package

# Scree plot to select number of factors:
plot(eigen(cov_encode_data)$values, xlab = 'Eigenvalue Number', ylab = 'Eigenvalue Size', main = 'Scree Graph', type = 'b', xaxt = 'n')
axis(1, at = seq(1, 10, by = 1))
abline(h = 1)
# We decide to use 3 factors: elbow at 3

# Meaningful names for the output
rownames(cov_encode_data) = colnames(cov_encode_data)

library(psych)
library(caret)

# Now you can use the findCorrelation function
correlation_matrix <- cor(cov_encode_data)
highly_correlated <- findCorrelation(correlation_matrix, cutoff = 0.90)
print(highly_correlated)

# Delete highly correlated variables
cov_data <- cov_encode_data[-highly_correlated, -highly_correlated]

# Convert cov_encode_data to a matrix
cov_data <- as.matrix(cov_data)
dim(cov_data)
# # factanal function uses the MLE method:
# (fit = factanal(factors = 2, covmat = cov_data))

# Try a different factor analysis method
fit <- psych::fa(r = cov_data, nfactors = 3)

print(fit, digits = 2, cutoff=.3, sort = TRUE)
# plot factor 1 by factor 2 
load = fit$loadings[,1:2] 
load
plot(load,type="n") # set up plot
text(load, cex=.7) # add variable names
# Plot the first two factors

