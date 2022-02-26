
## Convergence of chained equation algorithms

rm(list=ls())


library(mice)
library(VIM)
library(dplyr)
library(mi)
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
  number_variables <- 6
  mu <- runif(number_variables, -5, 10) 
  v <- runif(number_variables, 0, 15)
  
  y <- rnorm(n, mu[1], v[1])
  x1 <- 18 - 0.7 * y + rnorm(n, mu[2], v[2])
  x2 <- 5 + 0.5 * y + 0.7 * x1 + rnorm(n, mu[3], v[3])
  x3 <- 1 - 0.3 * y - 1 * x1 + 0.2 * x2 + rnorm(n, mu[4], v[4])
  x4 <-13 + 0.9 * y + 0.4 * x1 - 0.6 * x2 - 0.5 * x3 + rnorm(n, mu[5], v[5])
  x5 <- 6 + 0.1 * y - 1.3 * x1 + 0.75 * x2 - 0.2 * x3 + 0.8 * x4 + rnorm(n, mu[6], v[6])
  data_orig <- data.frame(y, x1, x2, x3, x4, x5)
  
  ## original values
  (means_orig <- apply(data_orig, 2, FUN = mean))
  (vars_orig <- apply(data_orig, 2, FUN = var))
  (cov_orig <- cov(data_orig))
  (cor_orig <- cor(data_orig))
  
  # create missing values
  
  ## mcar
  for (p_miss in c(0.3, 0.7)){
    percent_miss <- p_miss
    data_mcar <- data_orig
    
    mcar_values1 <- sample(nrow(data_orig), round(percent_miss * n))
    mcar_values2 <- sample(nrow(data_orig), round(percent_miss * n))
    mcar_values3 <- sample(nrow(data_orig), round(percent_miss * n))
    mcar_values4 <- sample(nrow(data_orig), round(percent_miss * n))
    mcar_values5 <- sample(nrow(data_orig), round(percent_miss * n))
    
    is.na(data_mcar$y[mcar_values1]) <- T
    is.na(data_mcar$x1[mcar_values2]) <- T
    is.na(data_mcar$x2[mcar_values3]) <- T
    is.na(data_mcar$x3[mcar_values4]) <- T
    is.na(data_mcar$x4[mcar_values5]) <- T
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
    data_mar <- data_orig
    
    #remove the lowest quantile
    data_mar <- data_orig
    z_miss_mar_px <- 0.5 + 2 * x1 - 0.7 * x2 + 0.4 * x3 - 0.8 * x4 + 0.1 * x5 + rnorm(n, 0, 3)
    mis_mar_px <- z_miss_mar_px < quantile(z_miss_mar_px, p_miss)
    
    #für alle Variablen außer x5 durchführen
    z_miss_mar_px1 <- 0.3 + 2.5 * y - 0.5 * x2 + 0.6 * x3 - 0.9 * x4 + 0.2 * x5 + rnorm(n, 0, 3)
    mis_mar_px1 <- z_miss_mar_px1 < quantile(z_miss_mar_px1, p_miss)

    z_miss_mar_px2 <- 1.5 + 1.2 * y - 0.5 * x1 + 0.7 * x3 - 0.2 * x4 + 0.7 * x5 + rnorm(n, 0, 3)
    mis_mar_px2 <- z_miss_mar_px2 < quantile(z_miss_mar_px2, p_miss)
    
    z_miss_mar_px3 <- 0.7 + 0.2 * y - 0.8 * x1 + 0.3 * x2 + 0.1 * x4 + 0.5 * x5 + rnorm(n, 0, 3)
    mis_mar_px3 <- z_miss_mar_px3 < quantile(z_miss_mar_px3, p_miss)
    
    z_miss_mar_px4 <- 0.4 + 0.9 * y + 0.6 * x1 + 0.7 * x2 - 0.9 * x3 + 0.3 * x5 + rnorm(n, 0, 3)
    mis_mar_px4 <- z_miss_mar_px4 < quantile(z_miss_mar_px4, p_miss)
    
    #NA setzen der Beobachtungen
    data_mar$y[mis_mar_px] <- NA
    data_mar$x1[mis_mar_px1] <- NA
    data_mar$x2[mis_mar_px2] <- NA
    data_mar$x3[mis_mar_px3] <- NA
    data_mar$x4[mis_mar_px4] <- NA
    #Muster mal ansehen
    #Anteile überprüfen, ob 0.3 hinkommt.
    summary(data_mar$y)
    summary(data_orig$y)
    
    #Hr. Meinfelders Muster
    data_my <- data_orig
    #Da die Daten auf Zufallszügen basieren, ist hier kein richtiges Sampling notwendig
    #10% der Daten als CC
    #Gleiche große Anzahl an Missings in den Variablen, nicht realisisch aber ist es problematisch?
    #Funktioniert das so?! Ansonsten wie in Zeilen 102-106 anpassen
    is.na(data_my$y[1:180]) <- T
    is.na(data_my$x1[181:360]) <- T
    is.na(data_my$x2[361:540]) <- T
    is.na(data_my$x3[541:720]) <- T
    is.na(data_my$x4[721:900]) <- T
    
    # Imputation
    
    ## mice
    #mice mcar
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
    
    #mice mar
    mice_mar <- mice(data_mar, maxit=5)
    summary(mice_mar)
    
    means_mar <- mice_mar %>% complete() %>% colMeans()
    vars_mar <- apply(mice_mar %>% complete(), 2, FUN = var)
    cov_mar <- mice_mar %>% complete() %>% cov()
    cor_mar <- mice_mar %>% complete() %>% cor(., use = "pairwise.complete.obs")
    
    #mice "my"
    mice_my <- mice(data_my, maxit=5)
    summary(mice_my)
    
    means_my <- mice_my %>% complete() %>% colMeans()
    vars_my <- apply(mice_my %>% complete(), 2, FUN = var)
    cov_my <- mice_my %>% complete() %>% cov()
    cor_my <- mice_my %>% complete() %>% cor(., use = "pairwise.complete.obs")
    
    # Rhat
    
    if (p_miss == 0.3){
      R_mean_30[i] = Rhat.mice(mice_mcar)[3]
      R_var_30[i] = Rhat.mice(mice_mcar)[4]
    }
    if (p_miss == 0.7){
      R_mean_70[i] = Rhat.mice(mice_mcar)[3]
      R_var_70[i] = Rhat.mice(mice_mcar)[4]
    }
    #mit mice_my und mice_mar hinzufügen
    
    
    # R_mean[i] = Rhat.mice(mice_mcar)[3]
    #  R_var[i] = Rhat.mice(mice_mcar)[4]
  }
}

# how many means have converged (3 variables x 100 iterations) 
# i have changed the number of iterations for computation time so it isnt at 100 yet
lapply(R_mean_30,function(x) x[which(x>1.1)]) %>%  unlist() %>% length()
lapply(R_mean_70,function(x) x[which(x>1.1)]) %>%  unlist() %>% length()

# how many variances have converges (3 variables x 100 iterations)
lapply(R_var_30,function(x) x[which(x>1.1)]) %>%  unlist() %>% length()
lapply(R_var_70,function(x) x[which(x>1.1)]) %>%  unlist() %>% length()

#show the values that are above 1.1
print(as.data.frame(t(lapply(R_mean_50,function(x) x[which(x>1.1)]))))


###MNAR ist nicht relevant, oder? Kann somit entfernt werden?! Sowie das Swiss Cheese von Paul?
# MNAR 
data_mnar <- data_orig
z_miss_mnar_p <-
  0.5 + 1 * x1 - 0.7 * x2 - 5 * y + rnorm(n, 0, 3)
mis_mnar_p <- z_miss_mnar_p < quantile(z_miss_mnar_p, 0.1)
data_mnar$y[mis_mnar_p] <- NA
summary(data_mnar$y)
summary(data_orig$y)






# SWISS CHEESE PATTERN VON PAUL:
x_df <- x_df_org

mis_x4 <- rbinom(n, 1, 0.2) %>% as.logical()
mis_x5 <- rbinom(n, 1, 0.4) %>% as.logical()
mis_x6 <- rbinom(n, 1, 0.6) %>% as.logical()

x_df[cbind(matrix(rep(FALSE, n * 3), n, 3), mis_x4, mis_x5, mis_x6)] <-
  NA

aggr(x_df)