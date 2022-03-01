## Convergence of chained equation algorithms

rm(list=ls())

library(ggplot2)
library(mice)
library(VIM)
library(dplyr)
library(mi)
library(miceadds)


#Listen für Rhat-Schätzer
R_mean_30 = list(len=100) #mcar
R_mean_70 = list(len=100)
R_mean_mar_30 = list(len=100)
R_mean_mar_70 = list(len=100)
R_mean_my = list(len=100)
R_var_30 = list(len=100)
R_var_70 = list(len=100)
R_var_mar_30 = list(len=100)
R_var_mar_70 = list(len=100)
R_var_my = list(len=100)
R_Rsq_30 = list(len=100)
R_Rsq_70 = list(len=100)
R_Rsq_mar_30 = list(len=100)
R_Rsq_mar_70 = list(len=100)
R_Rsq_my = list(len=100)

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

# The matrix that will later be put into the function Rhat1 will be called mat1
mat1 <- matrix(nrow=5, ncol=5)
rhat_sq <- vector(length=5) # temporary saved value of rhat for R²


#Reproduzierbarkeit
set.seed(1)

# Data generation

for (i in 1:2){
  
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
  
  ## original values
  #(means_orig <- apply(data_orig, 2, FUN = mean))
  #(vars_orig <- apply(data_orig, 2, FUN = var))
  #(cov_orig <- cov(data_orig))
  #(cor_orig <- cor(data_orig))
  
  # create missing values
  
  ## mcar
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
    #x5 bleibt vollständig
    ## graphics of missing pattern
    
    #aggr(data_mcar, numbers = T, prop = F)
    #md.pattern(data_mcar, plot = T)
    
    #matrixplot(data_mcar, sortby = "x")
    #matrixplot(data_mcar, sortby = "X1")
    #matrixplot(data_mcar, sortby = "X2")
    
    #marginplot(data_mcar[, c("y", "X1")])
    #marginplot(data_mcar[, c("y", "X2")])
    #marginplot(data_mcar[, c("X1", "X2")])
    
    ## mar
    
    #remove the lowest quantile
    data_mar <- data_orig
    z_miss_mar_py <- 0.5 + 2 * X1 - 0.7 * X2 + 0.4 * X3 - 0.8 * X4 + 0.1 * X5 + rnorm(n, 0, 3)
    mis_mar_py <- z_miss_mar_py < quantile(z_miss_mar_py, percent_miss)
    
    #für alle Variablen außer x5 durchführen
    z_miss_mar_pX1 <- 0.3 + 2.5 * y - 0.5 * X2 + 0.6 * X3 - 0.9 * X4 + 0.2 * X5 + rnorm(n, 0, 3)
    mis_mar_pX1 <- z_miss_mar_pX1 < quantile(z_miss_mar_pX1, percent_miss)
    
    z_miss_mar_pX2 <- 1.5 + 1.2 * y - 0.5 * X1 + 0.7 * X3 - 0.2 * X4 + 0.7 * X5 + rnorm(n, 0, 3)
    mis_mar_pX2 <- z_miss_mar_pX2 < quantile(z_miss_mar_pX2, percent_miss)
    
    z_miss_mar_pX3 <- 0.7 + 0.2 * y - 0.8 * X1 + 0.3 * X2 + 0.1 * X4 + 0.5 * X5 + rnorm(n, 0, 3)
    mis_mar_pX3 <- z_miss_mar_pX3 < quantile(z_miss_mar_pX3, percent_miss)
    
    z_miss_mar_pX4 <- 0.4 + 0.9 * y + 0.6 * X1 + 0.7 * X2 - 0.9 * X3 + 0.3 * X5 + rnorm(n, 0, 3)
    mis_mar_pX4 <- z_miss_mar_pX4 < quantile(z_miss_mar_pX4, percent_miss)
    
    #NA setzen der Beobachtungen
    data_mar$y[mis_mar_py] <- NA
    data_mar$X1[mis_mar_pX1] <- NA
    data_mar$X2[mis_mar_pX2] <- NA
    data_mar$X3[mis_mar_pX3] <- NA
    data_mar$X4[mis_mar_pX4] <- NA
    #Muster mal ansehen
    #Anteile überprüfen, ob 0.3 hinkommt.
    #summary(data_mar$y)
    #summary(data_orig$y)
    
    
    # Imputation
    
    ## mice
    #mice mcar
    mice_mcar <- mice(data_mcar, maxit=5)
    summary(mice_mcar)
    
    #Funktioniert nicht
    #means_mcar <- mice_mcar %>% complete() %>% colMeans()
    #vars_mcar <- apply(mice_mcar %>% complete(), 2, FUN = var)
    #cov_mcar <- mice_mcar %>% complete() %>% cov()
    #cor_mcar <- mice_mcar %>% complete() %>% cor(., use = "pairwise.complete.obs")
    
    # This returns the difference of these values. As there are multiple numbers,
    # I only return the highest, cause it is the most problematic (or is it? 
    # They might be dependent on the original variance of the variables)
    #max(means_mcar - means_orig)
    #max(vars_mcar - vars_orig)
    #max(cov_mcar - cov_orig)
    #max(abs(cor_mcar) - abs(cor_orig))
    
    #mice mar
    mice_mar <- mice(data_mar, maxit=5)
    summary(mice_mar)
    
    #means_mar <- mice_mar %>% complete() %>% colMeans()
    #vars_mar <- apply(mice_mar %>% complete(), 2, FUN = var)
    #cov_mar <- mice_mar %>% complete() %>% cov()
    #cor_mar <- mice_mar %>% complete() %>% cor(., use = "pairwise.complete.obs")
    
    
    ##########################################################################
    ####################### Rhat for R² ######################################
    for (impute_data in c(data_mcar, data_mar)){
      
      
      # for the y variable as depvar
      for (iter in 1:5){
        mice_iter <- mice(data_mcar, maxit=iter) # simulate the first iteration
        models <- lapply(1:mice_iter$m, function(m){
          lm(y ~ .,
             data = mice::complete(mice_iter, action = m))}) # action iterates between the m
        
        for (m in 1:5){
          mat1[iter, m] <- summary(models[[m]])$r.squared
        }
      }
      rhat_sq[1] <- Rhat1(mat1)
      
      # for variable x1 as depvar
      for (iter in 1:5){
        mice_iter <- mice(data_mcar, maxit=iter) # simulate the first iteration
        models <- lapply(1:mice_iter$m, function(m){
          lm(X1 ~ .,
             data = mice::complete(mice_iter, action = m))}) # action iterates between the m
        
        for (m in 1:5){
          mat1[iter, m] <- summary(models[[m]])$r.squared
        }
      }
      rhat_sq[2] <- Rhat1(mat1)
      
      # for variable x2 as depvar
      for (iter in 1:5){
        mice_iter <- mice(data_mcar, maxit=iter) # simulate the first iteration
        models <- lapply(1:mice_iter$m, function(m){
          lm(X2 ~ .,
             data = mice::complete(mice_iter, action = m))}) # action iterates between the m
        
        for (m in 1:5){
          mat1[iter, m] <- summary(models[[m]])$r.squared
        }
      }
      rhat_sq[3] <- Rhat1(mat1)
      
      # for variable x3 as depvar
      for (iter in 1:5){
        mice_iter <- mice(data_mcar, maxit=iter) # simulate the first iteration
        models <- lapply(1:mice_iter$m, function(m){
          lm(X3 ~ .,
             data = mice::complete(mice_iter, action = m))}) # action iterates between the m
        
        for (m in 1:5){
          mat1[iter, m] <- summary(models[[m]])$r.squared
        }
      }
      rhat_sq[4] <- Rhat1(mat1)  
      
      # for variable x4 as depvar
      for (iter in 1:5){
        mice_iter <- mice(data_mcar, maxit=iter) # simulate the first iteration
        models <- lapply(1:mice_iter$m, function(m){
          lm(X4 ~ .,
             data = mice::complete(mice_iter, action = m))}) # action iterates between the m
        
        for (m in 1:5){
          mat1[iter, m] <- summary(models[[m]])$r.squared
        }
      }
      rhat_sq[5] <- Rhat1(mat1)
      
      # Rhat
      if (percent_miss == 0.3){
        R_mean_30[i] = Rhat.mice(mice_mcar)[3]
        R_var_30[i] = Rhat.mice(mice_mcar)[4]
        R_Rsq_30[[i]] = rhat_sq
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
  }
  #Hr. Meinfelders Muster
  data_my <- data_orig
  #Da die Daten auf Zufallszügen basieren, ist hier kein richtiges Sampling notwendig
  #10% der Daten als CC
  #Gleiche große Anzahl an Missings in den Variablen, nicht realisisch aber ist es problematisch?
  #Funktioniert das so? Ansonsten wie in Zeilen 102-106 anpassen
  is.na(data_my$y[1:540]) <- T
  is.na(data_my$X1[1:540]) <- T
  is.na(data_my$X2[1:540]) <- T
  is.na(data_my$X3[541:900]) <- T
  is.na(data_my$X4[541:900]) <- T
  
  #mice "my"
  mice_my <- mice(data_my, maxit=5)
  summary(mice_my)
  
  R_mean_my[i] <-  Rhat.mice(mice_my)[3]
  R_var_my[i] <-  Rhat.mice(mice_my)[4]
  
  
  
  #############################################################################
  ############################# Rhat of R² ####################################
  # for the y variable as depvar
  for (iter in 1:5){
    mice_iter <- mice(data_my, maxit=iter) # simulate the first iteration
    models <- lapply(1:mice_iter$m, function(m){
      lm(y ~ .,
         data = mice::complete(mice_iter, action = m))}) # action iterates between the m
    
    for (m in 1:5){
      mat1[iter, m] <- summary(models[[m]])$r.squared
    }
  }
  rhat_sq[1] <- Rhat1(mat1)
  
  # for variable x1 as depvar
  for (iter in 1:5){
    mice_iter <- mice(data_my, maxit=iter) # simulate the first iteration
    models <- lapply(1:mice_iter$m, function(m){
      lm(X1 ~ .,
         data = mice::complete(mice_iter, action = m))}) # action iterates between the m
    
    for (m in 1:5){
      mat1[iter, m] <- summary(models[[m]])$r.squared
    }
  }
  rhat_sq[2] <- Rhat1(mat1)
  
  # for variable x2 as depvar
  for (iter in 1:5){
    mice_iter <- mice(data_my, maxit=iter) # simulate the first iteration
    models <- lapply(1:mice_iter$m, function(m){
      lm(X2 ~ .,
         data = mice::complete(mice_iter, action = m))}) # action iterates between the m
    
    for (m in 1:5){
      mat1[iter, m] <- summary(models[[m]])$r.squared
    }
  }
  rhat_sq[3] <- Rhat1(mat1)
  
  # for variable x3 as depvar
  for (iter in 1:5){
    mice_iter <- mice(data_my, maxit=iter) # simulate the first iteration
    models <- lapply(1:mice_iter$m, function(m){
      lm(X3 ~ .,
         data = mice::complete(mice_iter, action = m))}) # action iterates between the m
    
    for (m in 1:5){
      mat1[iter, m] <- summary(models[[m]])$r.squared
    }
  }
  rhat_sq[4] <- Rhat1(mat1)
  
  # for variable x4 as depvar
  for (iter in 1:5){
    mice_iter <- mice(data_my, maxit=iter) # simulate the first iteration
    models <- lapply(1:mice_iter$m, function(m){
      lm(X4 ~ .,
         data = mice::complete(mice_iter, action = m))}) # action iterates between the m
    
    for (m in 1:5){
      mat1[iter, m] <- summary(models[[m]])$r.squared
    }
  }
  rhat_sq[5] <- Rhat1(mat1)
  R_Rsq_my[[i]] <-  rhat_sq
  
  
}

# This returns the share of the variables that have converged (have an Rhat<1.1)
# It is calculated by taking the absolute number of converged variables and 
# dividing it by the total number of tested variables i*5
# i = the amount of loops/ the amount of generated datasets that are tested
# 5 = the 5 variables that have missing values 

#Anteil an konvergierten Werten
con_mean_mcar_30 <- (lapply(R_mean_30,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)
con_mean_mcar_70 <- (lapply(R_mean_70,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)
con_var_mcar_30 <- (lapply(R_var_30,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)
con_var_mcar_70 <- (lapply(R_var_70,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)
con_Rsq_mcar_30 <- (lapply(R_Rsq_30,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)
con_Rsq_mcar_70 <- (lapply(R_Rsq_70,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)

# how many variances have converges (6 variables x 100 iterations)


#show the values that are under 1.1
print(as.data.frame(t(lapply(R_mean_30,function(x) x[which(x<1.1)]))))

###MAR konvergierte Anteile
con_mean_mar_30 <- (lapply(R_mean_mar_30,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)
con_mean_mar_70 <-(lapply(R_mean_mar_70,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)
con_var_mar_30 <-(lapply(R_var_mar_30,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)
con_var_mar_70 <-(lapply(R_var_mar_70,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)
con_Rsq_mar_30 <-(lapply(R_Rsq_mar_30,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)
con_Rsq_mar_70 <-(lapply(R_Rsq_mar_70,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)

###Hr. Meinfelders Muster konvergierte Anteile
con_mean_my <- (lapply(R_mean_my,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)
con_var_my <- (lapply(R_var_my,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)
con_Rsq_my <- (lapply(R_Rsq_my,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)

###Collect data for sample
missmech <- c(rep("MCAR30", 3), rep("MCAR70", 3), rep("MAR30", 3), rep("MAR70", 3), rep("MY_PATTERN", 3) )
est <- rep(c("Mean" , "Variance" , "R^2"), 5)
prop <- c(con_mean_mcar_30, con_var_mcar_30, con_Rsq_mcar_30, con_mean_mcar_70, con_var_mcar_70, con_Rsq_mcar_70, con_mean_mar_30, con_var_mar_30, con_Rsq_mar_30, con_mean_mar_70, con_var_mar_70, con_Rsq_mar_70, con_mean_my, con_var_my, con_Rsq_my)
plotdata <- data.frame(missmech, est, prop)

# Plot
ggplot(plotdata, aes(fill=est, y=prop, x=missmech)) + 
  geom_bar(position="dodge", stat="identity") +
  xlab("Missing Mechanisms") +
  ylab("Proportion of Converged Estimators") +
  scale_fill_discrete(name= "Estimators") +
  ggtitle("Proportions of Converged Estimators by Pattern") +
  theme_bw()
