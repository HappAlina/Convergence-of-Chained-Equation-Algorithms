### Convergence of chained equation algorithms

rm(list=ls())

library(ggplot2)
library(mice)
#library(VIM)
library(dplyr)
#library(mi)
library(miceadds)


#Lists for Rhat-Estimator, that will be calculated in the for loop
R_mean_30 = list(len=100) #mcar
R_mean_70 = list(len=100)
R_mean_mar_30 = list(len=100)
R_mean_mar_70 = list(len=100)
R_mean_flex = list(len=100)
R_var_30 = list(len=100)
R_var_70 = list(len=100)
R_var_mar_30 = list(len=100)
R_var_mar_70 = list(len=100)
R_var_flex = list(len=100)
R_Rsq_30 = list(len=100)
R_Rsq_70 = list(len=100)
R_Rsq_mar_30 = list(len=100)
R_Rsq_mar_70 = list(len=100)
R_Rsq_flex = list(len=100)

# The following function belongs to the source code of the Rhat.mice function of the 
# miceadds package and will be used later on in our code. See:
# https://rdrr.io/cran/miceadds/src/R/Rhat1.R
# https://rdrr.io/cran/miceadds/src/R/Rhat.mice.R
# We will later implement a Rhat for the R-squared to check if the variables'
# interdependency converged. Therefore we cannot use the Rhat.mice function,
# which can only calculate Rhats for mean and variance
# The function takes a matrix with the iterations (n) as rows and the multiple 
# imputations (m) as columns 

Rhat1 <- function(mat)
{
  m <- ncol(mat)
  n <- nrow(mat)
  b <- apply(mat,2,mean)
  B <- sum((b-mean(mat))^2)*n/(m-1)
  w <- apply(mat,2, stats::var)
  W <- mean(w)
  s2hat <- (n-1)/n*W + B/n
  Vhat <- s2hat + B/m/n
  covWB <- n /m * (stats::cov(w,b^2)-2*mean(b)*stats::cov(w,b))
  varV <- (n-1)^2 / n^2 * stats::var(w)/m +
    (m+1)^2 / m^2 / n^2 * 2*B^2/(m-1) +
    2 * (m-1)*(n-1)/m/n^2 * covWB
  df <- 2 * Vhat^2 / varV
  R <- sqrt((df+3) * Vhat / (df+1) / W)
  return(R)
}

# The matrices that will later be put into the function Rhat1 
maty <- matrix(nrow=5, ncol=5)
matX1 <- matrix(nrow=5, ncol=5)
matX2 <- matrix(nrow=5, ncol=5)
matX3 <- matrix(nrow=5, ncol=5)
matX4 <- matrix(nrow=5, ncol=5)
rhat_sq <- vector(length=5) # temporary saved value of rhat for R²


#Reproducibility
set.seed(2)

##### Data generation #####

for (i in 1:100){
  
  n <- 1000
  number_variables <- 6
  mu <- runif(number_variables, -5, 10)
  v <- runif(number_variables, 0.1, 15)
  
  y <- rnorm(n, mu[1], v[1])
  X1 <- 18 - 0.7 * y + rnorm(n, mu[2], v[2])
  X2 <- 5 + 0.5 * y + 0.7 * X1 + rnorm(n, mu[3], v[3])
  X3 <- 1 - 0.3 * y - 1 * X1 + 0.2 * X2 + rnorm(n, mu[4], v[4])
  X4 <-13 + 0.9 * y + 0.4 * X1 - 0.6 * X2 - 0.5 * X3 + rnorm(n, mu[5], v[5])
  X5 <- 6 + 0.1 * y - 1.3 * X1 + 0.75 * X2 - 0.2 * X3 + 0.8 * X4 + rnorm(n, mu[6], v[6])
  data_orig <- data.frame(y, X1, X2, X3, X4, X5)
  
  
  ##### mcar #####
  ##create missing values
  #we test two missing proportions, 30% and 70%
  for (percent_miss in c(0.3, 0.7)){
    data_mcar <- data_orig
    
    mcar_values1 <- sample(nrow(data_orig), round(percent_miss * n))
    mcar_values2 <- sample(nrow(data_orig), round(percent_miss * n))
    mcar_values3 <- sample(nrow(data_orig), round(percent_miss * n))
    mcar_values4 <- sample(nrow(data_orig), round(percent_miss * n))
    mcar_values5 <- sample(nrow(data_orig), round(percent_miss * n))
    
    is.na(data_mcar$y[mcar_values1]) <- T
    is.na(data_mcar$X1[mcar_values2]) <- T
    is.na(data_mcar$X2[mcar_values3]) <- T
    is.na(data_mcar$X3[mcar_values4]) <- T
    is.na(data_mcar$X4[mcar_values5]) <- T
    #x5 stays complete

    ##### mar #####
    ##create missing values
    data_mar <- data_orig
    #numbers for dependencies were chosen randomly
    z_miss_mar_py <- 0.5 + 2 * X1 - 0.7 * X2 + 0.4 * X3 - 0.8 * X4 + 0.1 * X5 + rnorm(n, 0, 3)
    mis_mar_py <- z_miss_mar_py < quantile(z_miss_mar_py, percent_miss)
    
    #x5 stays complete
    z_miss_mar_pX1 <- 0.3 + 2.5 * y - 0.5 * X2 + 0.6 * X3 - 0.9 * X4 + 0.2 * X5 + rnorm(n, 0, 3)
    mis_mar_pX1 <- z_miss_mar_pX1 < quantile(z_miss_mar_pX1, percent_miss)
    
    z_miss_mar_pX2 <- 1.5 + 1.2 * y - 0.5 * X1 + 0.7 * X3 - 0.2 * X4 + 0.7 * X5 + rnorm(n, 0, 3)
    mis_mar_pX2 <- z_miss_mar_pX2 < quantile(z_miss_mar_pX2, percent_miss)
    
    z_miss_mar_pX3 <- 0.7 + 0.2 * y - 0.8 * X1 + 0.3 * X2 + 0.1 * X4 + 0.5 * X5 + rnorm(n, 0, 3)
    mis_mar_pX3 <- z_miss_mar_pX3 < quantile(z_miss_mar_pX3, percent_miss)
    
    z_miss_mar_pX4 <- 0.4 + 0.9 * y + 0.6 * X1 + 0.7 * X2 - 0.9 * X3 + 0.3 * X5 + rnorm(n, 0, 3)
    mis_mar_pX4 <- z_miss_mar_pX4 < quantile(z_miss_mar_pX4, percent_miss)
    
    #setting the lowest quantile NA
    data_mar$y[mis_mar_py] <- NA
    data_mar$X1[mis_mar_pX1] <- NA
    data_mar$X2[mis_mar_pX2] <- NA
    data_mar$X3[mis_mar_pX3] <- NA
    data_mar$X4[mis_mar_pX4] <- NA

    
    ##### Imputation #####
    
    #we use the mice package
    #Imputation for MCAR
    mice_mcar <- mice(data_mcar, maxit=5, printFlag = F, seed = 2)
   
    #Imputation for MAR
    mice_mar <- mice(data_mar, maxit=5, printFlag = F, seed = 2)
    
    
    ##### Rhat for R² #####
    #For each iteration and model (m) we need the imputed dataset. 
    #Therefore maxit is first set to 1 (then 2, ... until 5). 
    #For each iteration the default amount of models (m=5) will be created. 
    #For the 5x5=25 datasets a linear model is calculated and the 
    #corresponding R² is extracted in the sub-loop. 
    #The result is a 5x5 matrix of R²s (maty, matX1...). 
    #Now the function presented in the beginning of the script Rhat1 is used 
    #to get the Rhat-value of this matrix, which is then saved as rhat_sq.
    #We each build a linear model for y, X1, X2, X3 and 
    #X4 as the dependent variable.

    for (impute_data in c(data_mcar, data_mar)){
      
      for (iter in 1:5){
        mice_iter <- mice(data_mcar, maxit=iter, printFlag = F, seed = 2) # simulate the 5 iterations
        modely <- lapply(1:mice_iter$m, function(m){ # do for all 5 m
          lm(y ~ .,data = mice::complete(mice_iter, action = m))})
        modelX1 <- lapply(1:mice_iter$m, function(m){
          lm(X1 ~ .,data = mice::complete(mice_iter, action = m))})
        modelX2 <- lapply(1:mice_iter$m, function(m){
          lm(X2 ~ .,data = mice::complete(mice_iter, action = m))})
        modelX3 <- lapply(1:mice_iter$m, function(m){
          lm(X3 ~ .,data = mice::complete(mice_iter, action = m))})
        modelX4 <- lapply(1:mice_iter$m, function(m){
          lm(X4 ~ .,data = mice::complete(mice_iter, action = m))})
        for (m in 1:5){
          maty[iter, m] <- summary(modely[[m]])$r.squared
          matX1[iter, m] <- summary(modelX1[[m]])$r.squared
          matX2[iter, m] <- summary(modelX2[[m]])$r.squared
          matX3[iter, m] <- summary(modelX3[[m]])$r.squared
          matX4[iter, m] <- summary(modelX4[[m]])$r.squared
        }
      }
      rhat_sq[1] <- Rhat1(maty)
      rhat_sq[2] <- Rhat1(matX1)
      rhat_sq[3] <- Rhat1(matX2)
      rhat_sq[4] <- Rhat1(matX3)
      rhat_sq[5] <- Rhat1(matX4)
    }
      
      ##### Rhat for mean and variance from package miceadds #####
      if (percent_miss == 0.3){
        R_mean_30[i] = Rhat.mice(mice_mcar)[3] #col 3 is the Rhat of the mean
        R_var_30[i] = Rhat.mice(mice_mcar)[4] #col 4 is the Rhat of the variance
        R_Rsq_30[[i]] = rhat_sq #saving the calculated Rhats for R²
      }
      if (percent_miss == 0.7){
        R_mean_70[i] = Rhat.mice(mice_mcar)[3]
        R_var_70[i] = Rhat.mice(mice_mcar)[4]
        R_Rsq_70[[i]] = rhat_sq
      }
      if (percent_miss == 0.3){
        R_mean_mar_30[i] = Rhat.mice(mice_mar)[3]
        R_var_mar_30[i] = Rhat.mice(mice_mar)[4]
        R_Rsq_mar_30[[i]] = rhat_sq
      }
      if (percent_miss == 0.7){
        R_mean_mar_70[i] = Rhat.mice(mice_mar)[3]
        R_var_mar_70[i] = Rhat.mice(mice_mar)[4]
        R_Rsq_mar_70[[i]] = rhat_sq
      }
  }
  
  
  ##### Flex Pattern #####
  
  #When different questionnaires are given to two groups of people
  #with only a few overlapping questions (e.g. demographics),
  #there will be a missing pattern by design. 
  #This pattern will be simulated here.
  
  data_flex <- data_orig
  
  #Since the order of the observations is random, no sampling is necessary.
  #We build two groups: Group 1 is observed for y, X1 and X2 
  #and group 2 is observed for X3 and X4.
  #X5 is observed in both groups.
  #Only 10 % of the data are complete cases. 
  
  is.na(data_flex$y[1:540]) <- T
  is.na(data_flex$X1[1:540]) <- T
  is.na(data_flex$X2[1:540]) <- T
  is.na(data_flex$X3[541:900]) <- T
  is.na(data_flex$X4[541:900]) <- T
  
  #Imputation for Flex Pattern
  mice_flex <- mice(data_flex, maxit=5, printFlag = F, seed = 2)
  
  R_mean_flex[i] <-  Rhat.mice(mice_flex)[3]
  R_var_flex[i] <-  Rhat.mice(mice_flex)[4]
  
  # Rhat for R²
  #The same code as for MCAR and MAR is used but we have to repeat it here
  #because this pattern does not differentiate between 30 % and 70 %
  
  for (iter in 1:5){
    mice_iter <- mice(data_mcar, maxit=iter, printFlag = F, seed = 2) # simulate the 5 iterations
    modely <- lapply(1:mice_iter$m, function(m){ # do for all 5 m
      lm(y ~ .,data = mice::complete(mice_iter, action = m))})
    modelX1 <- lapply(1:mice_iter$m, function(m){
      lm(X1 ~ .,data = mice::complete(mice_iter, action = m))})
    modelX2 <- lapply(1:mice_iter$m, function(m){
      lm(X2 ~ .,data = mice::complete(mice_iter, action = m))})
    modelX3 <- lapply(1:mice_iter$m, function(m){
      lm(X3 ~ .,data = mice::complete(mice_iter, action = m))})
    modelX4 <- lapply(1:mice_iter$m, function(m){
      lm(X4 ~ .,data = mice::complete(mice_iter, action = m))})
    for (m in 1:5){
      maty[iter, m] <- summary(modely[[m]])$r.squared
      matX1[iter, m] <- summary(modelX1[[m]])$r.squared
      matX2[iter, m] <- summary(modelX2[[m]])$r.squared
      matX3[iter, m] <- summary(modelX3[[m]])$r.squared
      matX4[iter, m] <- summary(modelX4[[m]])$r.squared
    }
  }
  rhat_sq[1] <- Rhat1(maty)
  rhat_sq[2] <- Rhat1(matX1)
  rhat_sq[3] <- Rhat1(matX2)
  rhat_sq[4] <- Rhat1(matX3)
  rhat_sq[5] <- Rhat1(matX4)
  
  R_Rsq_flex[[i]] <-  rhat_sq
  
  print(i) # print the iteration to be aware how long it will still take
  print("iteration")
  
}

##### Calculate the share of converged estimators #####

# This returns the share of the variables that have converged (have an Rhat<1.1)
# It is calculated by taking the absolute number of converged variables and 
# dividing it by the total number of tested variables i*5
# i = the amount of loops/ the amount of generated datasets that are tested
# 5 = the 5 variables that have missing values 

#MCAR
con_mean_mcar_30 <- (lapply(R_mean_30,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)
con_mean_mcar_70 <- (lapply(R_mean_70,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)
con_var_mcar_30 <- (lapply(R_var_30,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)
con_var_mcar_70 <- (lapply(R_var_70,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)
con_Rsq_mcar_30 <- (lapply(R_Rsq_30,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)
con_Rsq_mcar_70 <- (lapply(R_Rsq_70,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)

#MAR
con_mean_mar_30 <- (lapply(R_mean_mar_30,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)
con_mean_mar_70 <-(lapply(R_mean_mar_70,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)
con_var_mar_30 <-(lapply(R_var_mar_30,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)
con_var_mar_70 <-(lapply(R_var_mar_70,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)
con_Rsq_mar_30 <-(lapply(R_Rsq_mar_30,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)
con_Rsq_mar_70 <-(lapply(R_Rsq_mar_70,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)

#Flex Pattern
con_mean_flex <- (lapply(R_mean_flex,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)
con_var_flex <- (lapply(R_var_flex,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)
con_Rsq_flex <- (lapply(R_Rsq_flex,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)

##### Graphical representation #####
#Collect data for graph
#Make a data frame from the share of the converged estimators calculated before
missmech <- c(rep("MCAR30", 3), rep("MCAR70", 3), rep("MAR30", 3), rep("MAR70", 3), rep("Flex-Pattern", 3) )
est <- rep(c("Mean" , "Variance" , "R²"), 5)
prop <- c(con_mean_mcar_30, con_var_mcar_30, con_Rsq_mcar_30, con_mean_mcar_70, con_var_mcar_70, con_Rsq_mcar_70, con_mean_mar_30, con_var_mar_30, con_Rsq_mar_30, con_mean_mar_70, con_var_mar_70, con_Rsq_mar_70, con_mean_flex, con_var_flex, con_Rsq_flex)
plotdata <- data.frame(missmech, est, prop)

# Plot
ggplot(plotdata, aes(fill=est, y=prop, x=missmech)) + 
  geom_bar(position="dodge", stat="identity") +
  xlab("Missing Mechanisms & Pattern") +
  ylab("Share of Converged Estimators") +
  scale_fill_manual(name="Estimators", values = c("Mean" = "lightblue", "R²" = "lightcoral", "Variance" = "grey")) +
  theme_bw()
