library(foreach)
library(doParallel)
library(doRNG)
library(Rfast)

simulation <- function(N=200, dist="bern", par1 = c(0,0),  par2 = NULL, burnin=3, ar=NULL, nsim= 10^4, method="ER", measure="sd", signlevel=0.05){
  
  #Generate Result Data Frame
  colnames <-  c()
  for(i in 1:ncol(par1)){
    colnames <- c(colnames, paste("p1", i, sep = "", collapse = NULL))
  }
  if(!is.null(par2)){
    for(i in 1:ncol(par2)){
      colnames <- c(colnames, paste("p2", i, sep = "", collapse = NULL))
    }
  }
  for(j in 1:2){
    colnames <- c(colnames, paste( "n", j-1, sep = "", collapse = NULL))
    colnames <- c(colnames, paste("Var(n", j-1, ")",sep = "", collapse = NULL))
  }
  colnames <- c(colnames, "EMR", "Var_EMR",  "Z", "Z_A", "Z_F","BM", "Per_Superior", "Var_Superior")
  if(!is.null(ar)){
    colnames <- c(colnames, ar)
    colnames <- c(colnames, paste("BIAS(", ar, ")",sep = "", collapse = NULL))
  }
  result <- data.frame(matrix(ncol = length(colnames), nrow = max(nrow(par1),nrow(par2))))
  colnames(result) <- colnames
  
  cl <- parallel::makeCluster(4)  #increase according to the limitations of your cluster if you run it on the cluster
  doParallel::registerDoParallel(cl)
  
  ########## Theoretical Variance Calculation ############################
  if(dist=="bern"){
      mean <- par1
      fun <- function(x){x*(1-x)}
      variance <-  data.frame(lapply(par1,fun))
      theta <- par1[,2]*(1-par1[,1]) + 0.5*(par1[,1]*par1[,2]+(1-par1[,1])*(1-par1[,2]))
      psi <- 0.5*theta +0.25
  }
  if(dist=="norm"){
      mean <- par1
      fun <- function(x){x^2}
      variance <-  data.frame(lapply(par2,fun))
      xn <- (mean[,1] - mean[,2])/sqrt(variance[,1]+variance[,2])
      theta <-  1 - dnorm(xn)
      psi <- 0.5*theta +0.25
  }
  if(dist=="expon"){
      fun <- function(x){1/x}
      mean <- data.frame(lapply(par1,fun))
      fun <- function(x){1/(x^2)}
      variance <-  data.frame(lapply(par1,fun))
      theta <-  par1[,1]/(par1[,1]+par1[,2])
      psi <- 0.5*theta +0.25
  }
  if(dist=="beta"){
      fun <- function(x){5/(x+5)}
      mean <- data.frame(lapply(par1,fun))
      fun <- function(x){5*x/((5+x)^2+(5+x+1))}
      variance <-  data.frame(lapply(par1,fun))
      theta <- c(0.5,0.5631895, 0.6317465, 0.704089, 0.7777655, 0.813794,0.8486965,
                 0.881401, 0.91132, 0.937633, 0.959517, 0.976518, 0.988456, 0.991853,
                  0.994524, 0.995621, 0.997357, 0.998557, 0.998991, 0.999326, 0.999584)
      psi <- 0.5*theta +0.25
      psi <- c(rep(0.5, length(theta)),theta[-1])
  }
  if(dist=="Lickert"){
      ## Variance needed to be estimated through simulations in extra Code
      mean <-  read.csv("Lickert_Mean.csv", header = TRUE)
      mean <- mean[,2:3]
      variance <-  read.csv("Lickert_Variance.csv", header = TRUE)
      variance <- variance[,2:3]
      theta <- c(0.5, 0.544, 0.593, 0.6446757, 0.696642, 0.7469708, 0.796875, 0.818026,
                  0.840608, 0.852698, 0.8653635, 0.878488, 0.892124, 0.906136, 0.92029,
                  0.93435, 0.9479362, 0.9606685, 0.9720625, 0.9817345, 0.9893505,
                  0.9947455, 0.998026, 0.999626)
      psi <- 0.5*theta +0.25
      psi <- c(rep(0.5, length(psi)),psi[-1])
  }
  ############## ER ######################################################
  if(method == "ER" ){
    result_part_emp <- data.frame(matrix(ncol = length(colnames), nrow = 0))
    result.ER <- foreach(i = 1:length(par1[,1]), .combine='rbind', .packages=c("BSDA", "rankFD")) %dorng% { #%dorng% 
  
      result_part <-result_part_emp
      source('ER.R')
      p <- as.numeric(par1[i,])
      start <- length(p)
      result_part[1,1:start] <- p
      if(is.null(par2)){
          sim.ER <- sim_ER(N=N, dist=dist, par1 = p,  burnin= burnin, nsim=nsim, measure=measure)
      } else{
        sig <- as.numeric(par2[i,])
        start <- start+ length(sig)
        result_part[1,1:start] <- c(p, sig)
        sim.ER <- sim_ER(N=N, dist=dist, par1 = p, burnin= burnin, par2=sig, nsim=nsim)
      }
      
      result_part[1,(start+1)] <- mean(sim.ER[,7]); result_part[1,(start+2)] <- var(sim.ER[,8]) 
      result_part[1,(start+3)] <- mean(sim.ER[,8]); result_part[1,(start+4)] <- var(sim.ER[,8]) 
          
      result_part[1,(start+5)] <- mean(sim.ER[,1])
      result_part[1,(start+6)] <- var(sim.ER[,1])
          
      reject_Z <- sum(sim.ER[,2] < signlevel)
      reject_Z_A <- sum(sim.ER[,3] < signlevel)
      reject_Z_F <- sum(sim.ER[,4] < signlevel)
      reject_BM <- sum(sim.ER[,5] < signlevel)
          
      result_part[1,(start+7)] <- reject_Z/nsim
      result_part[1,(start+8)] <- reject_Z_A/nsim
      result_part[1,(start+9)] <- reject_Z_F/nsim
      result_part[1,(start+10)] <- reject_BM/nsim
          
      result_part[1,(start+11)] <- mean(sim.ER[,6])
      result_part[1,(start+12)] <- var(sim.ER[,6])
          
      result_part
        
    } #%do for no parallel
    result[,1:(length(colnames))] <- result.ER
  }
  if(method == "FR" ){
    result_part_emp <- data.frame(matrix(ncol = length(colnames), nrow = 0))
    result.FR <- foreach(i = 1:length(par1[,1]), .combine='rbind', .packages=c("BSDA", "rankFD")) %dorng% { #%dorng% 
      
      result_part <-result_part_emp
      source('FR.R')
      p <- as.numeric(par1[i,])
      start <- length(p)
      result_part[1,1:start] <- p
      if(is.null(par2)){
        sim.FR <- sim_FR(N=N, dist=dist, par1 = p,  FR=ar, burnin= burnin, nsim=nsim, measure=measure)
      } else{
        sig <- as.numeric(par2[i,])
        start <- start+ length(sig)
        result_part[1,1:start] <- c(p, sig)
        sim.FR <- sim_FR(N=N, dist=dist, par1 = p, burnin= burnin,FR=ar, par2=sig, nsim=nsim)
      }
      
      result_part[1,(start+1)] <- mean(sim.FR[,7]); result_part[1,(start+2)] <- var(sim.FR[,8]) 
      result_part[1,(start+3)] <- mean(sim.FR[,8]); result_part[1,(start+4)] <- var(sim.FR[,8]) 
      
      result_part[1,(start+5)] <- mean(sim.FR[,1])
      result_part[1,(start+6)] <- var(sim.FR[,1])
      
      reject_Z <- sum(sim.FR[,2] < signlevel)
      reject_Z_A <- sum(sim.FR[,3] < signlevel)
      reject_Z_F <- sum(sim.FR[,4] < signlevel)
      reject_BM <- sum(sim.FR[,5] < signlevel)
      
      result_part[1,(start+7)] <- reject_Z/nsim
      result_part[1,(start+8)] <- reject_Z_A/nsim
      result_part[1,(start+9)] <- reject_Z_F/nsim
      result_part[1,(start+10)] <- reject_BM/nsim
      
      result_part[1,(start+11)] <- mean(sim.FR[,6])
      result_part[1,(start+12)] <- var(sim.FR[,6])
      
      result_part
      
    } #%do for no parallel
    result[,1:(length(colnames))] <- result.FR
  }
  if( method == "ERADE"){
    result_part_emp <- data.frame(matrix(ncol = length(colnames), nrow = 0))
    result.ERADE <- foreach(i = 1:length(par1[,1]), .combine='rbind', .packages=c("BSDA", "rankFD")) %dorng% { #%do for no parallel
      result_part <-result_part_emp
      source('ERADE.R')
      p <- as.numeric(par1[i,])
      mu <- as.numeric(mean[i,])
      var <- as.numeric(variance[i,])
      if(dist=="norm") {
        std <- as.numeric(par2[i,])
        start <- length(c(p,std))
        sim.ERADE <- sim_ERADE(N=N, dist=dist, par1 = p, par2 = std,  nsim=nsim, ar=ar, burnin=burnin)
        result_part[1,1:start] <- c(p,std)
      }else{
        start <- length(p)
        sim.ERADE <- sim_ERADE(N=N, dist=dist, par1 = p,   measure= measure, nsim=nsim, ar=ar, burnin=burnin)
        result_part[1,1:start] <- p
      }
          
      result_part[1,(start+1)] <- mean(sim.ERADE[,7]); result_part[1,(start+2)] <- var(sim.ERADE[,8]) 
      result_part[1,(start+3)] <- mean(sim.ERADE[,8]); result_part[1,(start+4)] <- var(sim.ERADE[,8]) 
          
      result_part[1,(start+5)] <- mean(sim.ERADE[,1])
      result_part[1,(start+6)] <- var(sim.ERADE[,1])
          
      reject_Z <- sum(sim.ERADE[,2] < signlevel)
      reject_Z_A <- sum(sim.ERADE[,3] < signlevel)
      reject_Z_F <- sum(sim.ERADE[,4] < signlevel)
      reject_BM <- sum(sim.ERADE[,5] < signlevel)
          
      result_part[1,(start+7)] <- reject_Z/nsim
      result_part[1,(start+8)] <- reject_Z_A/nsim
      result_part[1,(start+9)] <- reject_Z_F/nsim
      result_part[1,(start+10)] <- reject_BM/nsim
          
      result_part[1,(start+11)] <- mean(sim.ERADE[,6])
      result_part[1,(start+12)] <- var(sim.ERADE[,6])
      if(ar=="RSHIR_Z1"){
        result_part[1,(start+13)] <- 1- sqrt(p[1])/(sqrt(p[1]) + sqrt(p[2]))
      } 
      if(ar=="Neyman_Z1"){
        result_part[1,(start+13)] <- 1 - var[1]/(var[1]+var[2])
      }
      if(ar=="Neyman_Z0"){
        result_part[1,(start+13)] <- var[1]/(var[1]+var[2])
      }
      result_part[1,(start+14)] <- mean(sim.ERADE[,7]) - result_part[1,(start+13)]
      
      result_part
    }    
    result[,1:(length(colnames))] <- result.ERADE 
  }
    
  
  #Print and save Results
  if(method=="ERADE"){
    filename <-  paste(method, ar, N, burnin,dist,"All_Table_t1.csv",sep="_")
  } else {
    filename <-  paste(method, N, burnin, dist,"All_Table_t1.csv",sep="_")
  }
  
  write.csv(round(result,4), filename, row.names=TRUE)
  print(round(result,4))
  

  # Print Results
  rounded_data <- round(result[,c(1,2,9, 10,11,3,4,5,6)], 4)
  # Round the 7th row of 'result' and multiply by N
  rounded_N_value <- round(N * result[7], 1)
  # Combine the two results into one table (data frame)
  final_table <- cbind(rounded_data, rounded_N_value)
  final_table[,3] <- round(final_table[,3]*100,1)
  final_table[,4] <- round(final_table[,4]*100,1)
  final_table[,5] <- round(final_table[,5]*100,1)
  # Rename the last column as "ENS"
  colnames(final_table)[ncol(final_table)] <- "ENS"
  # Print the combined table with renamed column
  print(final_table)
  
  stopCluster(cl)
}

# Generate Parameter Data Frame
# Bernoulli
p0 <- c(0.05, 0.05)
p1 <- c(0.05, 0.3)

par1_Ber <- data.frame(p0,p1)

#Number of Simulations 
nsim = 10^4
n=60

# Simulations for Table 1
simulation(N=n, dist="bern", par1 = par1_Ber, nsim=nsim,  burnin=round(n/2), method = "ER", measure = "sd", signlevel = 0.05) 
simulation(N=n, dist="bern", par1 = par1_Ber, nsim=nsim,  burnin=6, method = "FR", ar=0.6667, measure = "sd", signlevel = 0.05) 
simulation(N=n, dist="bern", par1 = par1_Ber, nsim=nsim, burnin=6, method = "ERADE", ar="Neyman_Z1", measure = "sd", signlevel = 0.05) #10^4

