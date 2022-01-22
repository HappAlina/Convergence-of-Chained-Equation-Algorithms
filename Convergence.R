# Convergence of chained equation algorithms

library(mi)
library(mice)
library(VIM)
library(dplyr)
library(miceadds)


# Data generation

set.seed(1)
n <- 1000
number_variables <- 3
mu <- runif(number_variables, -5, 10) # warum hier nicht rnorm, beta ist ja normalverteilt

v <- rWishart(3, 3, matrix(0.5, number_variables, number_variables)) # not sure what to do here
# warum nicht inverse gamma, da SIgma ja glaub ich invers gamma verteilt ist

y <- rnorm(n, mu[1], v[1])
X1 <- 18 - 0.7 * y + rnorm(n, mu[2], v[2])
X2 <- 5 + 0.5 * y + 0.7 * X1 + rnorm(n, mu[3], v[3])

data_orig <- data.frame(y, X1, X2)

## original values
(means_orig <- apply(data_orig, 2, FUN = mean))
(vars_orig <- apply(data_orig, 2, FUN = var))
(cov_orig <- cov(data_orig))
(cor_orig <- cor(data_orig))

# create missing values

## mcar
percent_miss <- 0.3
data_mcar <- data_orig

mcar_values1 <- sample(nrow(data_orig), round(percent_miss * n))
mcar_values2 <- sample(nrow(data_orig), round(percent_miss * n))
mcar_values3 <- sample(nrow(data_orig), round(percent_miss * n))

is.na(data_mcar$y[mcar_values1]) <- T
is.na(data_mcar$X1[mcar_values2]) <- T
is.na(data_mcar$X2[mcar_values3]) <- T

##ALTERNATIVE MIT WAHRSCHEINLICHKEIT
# mis_x6 <- sample(
#   x = c(TRUE, FALSE),
#   size = n,
#   replace = TRUE,
#   prob = c(0.2, 0.8)
# ) & mis_x5

## graphics of missing pattern

aggr(data_mcar, numbers = T, prop = F)
md.pattern(data_mcar, plot = T)

matrixplot(data_mcar, sortby = "y")
matrixplot(data_mcar, sortby = "X1")
matrixplot(data_mcar, sortby = "X2")

marginplot(data_mcar[, c("y", "X1")])
marginplot(data_mcar[, c("y", "X2")])
marginplot(data_mcar[, c("X1", "X2")])

# imputation

## mice
mice_mcar <- mice(data_mcar)
summary(mice_mcar)

means_mcar <- mice_mcar %>% complete() %>% colMeans()
vars_mcar <- apply(mice_mcar %>% complete(), 2, FUN = var)
cov_mcar <- mice_mcar %>% complete() %>% cov()
cor_mcar <- mice_mcar %>% complete() %>% cor(., use = "pairwise.complete.obs")

# This returns the difference of these values. As there are multiple numbers,
# I only return the highest, cause it is the most problematic (or is it? 
# They might be dependent on the original variance of the variables)
max(means_mcar - means_orig)
max(vars_mcar - vars_orig)
max(cov_mcar - cov_orig)
max(abs(cor_mcar) - abs(cor_orig))

# Rhat
Rhat.mice(mice_mcar)







# MAR SACHEN VON PAUL (Session 1)
mar_values1 <- order(dat$tvtot + rnorm(n, 0, 0.5 * sd(dat$tvtot)),
                     decreasing = TRUE)[1:(round(0.3 * n))]


ind_mar  <- 1 - rbinom(n, 1, logistic(y[, 1]))   # depends on y1
prop.table(table(ind_mar))


dat_mar_p <- dat

# model missingness via linear regression model
# depending on x_1, x_2 and random error
z_miss_mar_p <- 0.5 + 2 * x_1 - 0.7 * x_2 + rnorm(n1, 0, 3)

plot(density(z_miss_mar_p))
abline(v = quantile(z_miss_mar_p, p_mis))

# MNAR SACHEN VON PAUL (Session 1)
ind_mnar <- 1 - rbinom(n, 1, logistic(y[, 2]))   # depends on y2
prop.table(table(ind_mnar))

# Siehe auch library(missMethods)

dat_mnar_p <- dat

# model missingness by logistic regression (probit):
z_miss_mnar_p <-
  0.5 + 1 * x_1 - 0.7 * x_2 - 5 * x_3 + rnorm(n1, 0, 3)  # the missing of a value now also depends on x_3 itself
plot(density(z_miss_mnar_p))
abline(v = quantile(z_miss_mnar_p, p_mis))







# SWISS CHEESE PATTERN VON PAUL:
x_df <- x_df_org

mis_x4 <- rbinom(n, 1, 0.2) %>% as.logical()
mis_x5 <- rbinom(n, 1, 0.4) %>% as.logical()
mis_x6 <- rbinom(n, 1, 0.6) %>% as.logical()

x_df[cbind(matrix(rep(FALSE, n * 3), n, 3), mis_x4, mis_x5, mis_x6)] <-
  NA

aggr(x_df)
