
# Convergence of chained equation algorithms

rm(list=ls())


library(mi)
library(mice)
library(VIM)
library(dplyr)
library(miceadds)


set.seed(1)
R_mean_10 = list(len=100)
R_mean_30 = list(len=100)
R_mean_50 = list(len=100)
R_mean_80 = list(len=100)


R_var_10 = list(len=100)
R_var_30 = list(len=100)
R_var_50 = list(len=100)
R_var_80 = list(len=100)
# Data generation

for (i in 1:20){
  
  n <- 1000
  number_variables <- 3
  mu <- runif(number_variables, -5, 10)
  v <- rWishart(3, 3, matrix(0.5, number_variables, number_variables)) # not sure what to do here
  
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
  for (p_miss in c(0.1, 0.3, 0.5, 0.8)){
    percent_miss <- p_miss
    data_mcar <- data_orig
    
    mcar_values1 <- sample(nrow(data_orig), round(percent_miss * n))
    mcar_values2 <- sample(nrow(data_orig), round(percent_miss * n))
    mcar_values3 <- sample(nrow(data_orig), round(percent_miss * n))
    
    is.na(data_mcar$y[mcar_values1]) <- T
    is.na(data_mcar$X1[mcar_values2]) <- T
    is.na(data_mcar$X2[mcar_values3]) <- T
    
    ## graphics of missing pattern
    
    #aggr(data_mcar, numbers = T, prop = F)
    #md.pattern(data_mcar, plot = T)
    
    #matrixplot(data_mcar, sortby = "y")
    #matrixplot(data_mcar, sortby = "X1")
    #matrixplot(data_mcar, sortby = "X2")
    
    #marginplot(data_mcar[, c("y", "X1")])
    #marginplot(data_mcar[, c("y", "X2")])
    #marginplot(data_mcar[, c("X1", "X2")])
    
    # imputation
    
    ## mice
    mice_mcar <- mice(data_mcar, maxit=5)
    summary(mice_mcar)
    
    means_mcar <- mice_mcar %>% complete() %>% colMeans()
    vars_mcar <- apply(mice_mcar %>% complete(), 2, FUN = var)
    cov_mcar <- mice_mcar %>% complete() %>% cov()
    cor_mcar <- mice_mcar %>% complete() %>% cor(., use = "pairwise.complete.obs")
    
    # This returns the difference of these values. As there are multiple numbers,
    # I only return the highest, cause it is the most problematic (or is it? 
    # They might be dependent on the original variance of the variables)
    #max(means_mcar - means_orig)
    #max(vars_mcar - vars_orig)
    #max(cov_mcar - cov_orig)
    #max(abs(cor_mcar) - abs(cor_orig))
    
    # Rhat
    
    if (p_miss == 0.1){
      R_mean_10[i] = Rhat.mice(mice_mcar)[3]
      R_var_10[i] = Rhat.mice(mice_mcar)[4]
    }
    if (p_miss == 0.3){
      R_mean_30[i] = Rhat.mice(mice_mcar)[3]
      R_var_30[i] = Rhat.mice(mice_mcar)[4]
    }
    if (p_miss == 0.5){
      R_mean_50[i] = Rhat.mice(mice_mcar)[3]
      R_var_50[i] = Rhat.mice(mice_mcar)[4]
    }
    if (p_miss == 0.8){
      R_mean_80[i] = Rhat.mice(mice_mcar)[3]
      R_var_80[i] = Rhat.mice(mice_mcar)[4]
    }
    
    
    
    # R_mean[i] = Rhat.mice(mice_mcar)[3]
    #  R_var[i] = Rhat.mice(mice_mcar)[4]
  }
}

# how many means have converged (3 variables x 100 iterations) 
# i have changed the number of iterations for computation time so it isnt at 100 yet
lapply(R_mean_10,function(x) x[which(x>1.1)]) %>%  unlist() %>% length()
lapply(R_mean_30,function(x) x[which(x>1.1)]) %>%  unlist() %>% length()
lapply(R_mean_50,function(x) x[which(x>1.1)]) %>%  unlist() %>% length()
lapply(R_mean_80,function(x) x[which(x>1.1)]) %>%  unlist() %>% length()

# how many variances have converges (3 variables x 100 iterations)
lapply(R_var_10,function(x) x[which(x>1.1)]) %>%  unlist() %>% length()
lapply(R_var_30,function(x) x[which(x>1.1)]) %>%  unlist() %>% length()
lapply(R_var_50,function(x) x[which(x>1.1)]) %>%  unlist() %>% length()
lapply(R_var_80,function(x) x[which(x>1.1)]) %>%  unlist() %>% length()

#show the values that are above 1.1
print(as.data.frame(t(lapply(R_mean_50,function(x) x[which(x>1.1)]))))



##############################################################################
################# MAR ########################################################
##############################################################################





## mar
percent_miss <- 0.1
data_mar <- data_orig

#remove the lowest quantile 
percent_miss <- 0.1
data_mar <- data_orig
z_miss_mar_p <- 0.5 + 2 * X1 - 0.7 * X2 + rnorm(n, 0, 3)
mis_mar_p <- z_miss_mar_p < quantile(z_miss_mar_p, 0.1)
data_mar$y[mis_mar_p] <- NA

summary(data_mar$y)
summary(data_orig$y)

# question: should we set more than one variable to mnar and the other ones?







# MNAR 
data_mnar <- data_orig
z_miss_mnar_p <-
  0.5 + 1 * X1 - 0.7 * X2 - 5 * y + rnorm(n, 0, 3)
mis_mnar_p <- z_miss_mnar_p < quantile(z_miss_mnar_p, 0.1)
data_mnar$y[mis_mnar_p] <- NA
summary(data_mnar$y)
summary(data_orig$y)


#a





# SWISS CHEESE PATTERN VON PAUL:
x_df <- x_df_org

mis_x4 <- rbinom(n, 1, 0.2) %>% as.logical()
mis_x5 <- rbinom(n, 1, 0.4) %>% as.logical()
mis_x6 <- rbinom(n, 1, 0.6) %>% as.logical()

x_df[cbind(matrix(rep(FALSE, n * 3), n, 3), mis_x4, mis_x5, mis_x6)] <-
  NA

aggr(x_df)