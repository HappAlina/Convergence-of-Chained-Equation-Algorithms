## Convergence of chained equation algorithms

rm(list=ls())


library(mice)
library(VIM)
library(dplyr)
library(mi)
library(miceadds)

#Reproduzierbarkeit
set.seed(1)
#Listen für Schätzer
R_mean_30 = list(len=100)
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

# Data generation

for (i in 1:10){
  
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
    summary(data_mar$y)
    summary(data_orig$y)
    
    
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

    # Rhat
    #aktuell werden alle Variablen übertragen, soll das so oder macht das nicht gerade die Probleme?
    if (percent_miss == 0.3){
      R_mean_30[i] = Rhat.mice(mice_mcar)[3]
      R_var_30[i] = Rhat.mice(mice_mcar)[4]
    }
    if (percent_miss == 0.7){
      R_mean_70[i] = Rhat.mice(mice_mcar)[3]
      R_var_70[i] = Rhat.mice(mice_mcar)[4]
    }
    if (percent_miss == 0.3){
      R_mean_mar_30[i] = Rhat.mice(mice_mar)[3]
      R_var_mar_30[i] = Rhat.mice(mice_mar)[4]
    }
    if (percent_miss == 0.7){
      R_mean_mar_70[i] = Rhat.mice(mice_mar)[3]
      R_var_mar_70[i] = Rhat.mice(mice_mar)[4]
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
  
  R_mean_my[i] = Rhat.mice(mice_my)[3]
  R_var_my[i] = Rhat.mice(mice_my)[4]
}

# how many means have converged (6 variables x 100 iterations) 
# i have changed the number of iterations for computation time so it isnt at 100 yet
(lapply(R_mean_30,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)
#Anteil an konvergierten Werten
(lapply(R_mean_70,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)

# how many variances have converges (6 variables x 100 iterations)
(lapply(R_var_30,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)
(lapply(R_var_70,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)

#show the values that are under 1.1
print(as.data.frame(t(lapply(R_mean_30,function(x) x[which(x<1.1)]))))

###MAR konvergierte Anteile
(lapply(R_mean_mar_30,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)
(lapply(R_mean_mar_70,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)
(lapply(R_var_mar_30,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)
(lapply(R_var_mar_70,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)

###Hr. Meinfelders Muster konvergierte Anteile
(lapply(R_mean_my,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)
(lapply(R_var_my,function(x) x[which(x<1.1)]) %>%  unlist() %>% length())/(i*5)

