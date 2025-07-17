library(BSDA)
library(rankFD)

source('Data_Generator.R')

sim_ER <-  function(N=200, dist="bern",  par1 = c(0,0), par2 = NULL, measure="sd", burnin=3, nsim=10^4, alpha=0.05, one.sided = FALSE){
  
  if(dist=="bern"){
    out = matrix(nrow = nsim, ncol = 8)
  }
  else {
    out = matrix(nrow = nsim, ncol = 6)
  }
    for(i in 1:nsim){
      out[i,] = two_arm_ER(N=N, dist=dist,  par1 =par1, par2 = par2, measure=measure,burnin = burnin, Z_corrected=FALSE)
    }
  return(out)
}

two_arm_ER <- function(N=200, dist="bern",  par1 = c(0,0), par2 = NULL, measure="sd", burnin=3, Z_corrected=FALSE){
  # Generate Samples 
  n1 = min(max(rbinom(1, N, 0.5),burnin), N-burnin)
  n0 = N - n1
  x0 <- generate_data(N=n0, dist=dist, parV=c(par1[1],par2[1])) 
  x1 <- generate_data(N=n1, dist=dist, parV=c(par1[2],par2[2]))
  
  # Response
  s0 <- sum(x0)
  s1 <- sum(x1)
  response <- (s0+s1)/N
  
  #Imbalance Measure
  sup <- superior(dist=dist,par1 = par1, par2 = par2)
  if(sup==2){ #2 means that no arm is theoretically supirior 
    per.sup <- NA
  } else {
    if(sup==0){ 
      per.sup <- n0/N
    }
    if(sup==1){ 
      per.sup <- n1/N
    }
  }
  
  ##### Z - Test & BM-Test
  source('WaldTest.R')
  ##### Z - Test 
  if(dist=="bern"){
    Z_P <-wald.test.binary(c(x0),c(x1), measure = measure)
  } else {
    Z_P <- wald.test(c(x0),c(x1))
  }
    # One-Point-Sample-Correction
    sigx <- sd(x0)
    sigy <- sd(x1)
    if(sigx==0 && sigy==0){
      BM_P <- Z_P[1]
    } else {
      #BM <- lawstat::brunner.munzel.test(c(x0,0,1),c(x1,0,1))
      BM <- lawstat::brunner.munzel.test(c(x0),c(x1))
      if(BM$statistic == Inf){
        BM_P <- 0
      } else if(BM$statistic == -Inf) {
        BM_P <-  1
      } else {
        BM_P <- BM$p.value
      }
    }
    n <- c(n0,n1)/N
  return(c(response, Z_P, BM_P, per.sup, n))
}

# Scenario 1: Type-I Error (No true difference between arms in means)
# dist="norm", par1 = c(mean_arm0, mean_arm1), par2=c(sd_arm0, sd_arm1)
# Here, par1 = c(0,0) means both arms have a mean of 0.
# par2=c(0.5,2) means arm0 has SD 0.5, arm1 has SD 2.
cat("--- Scenario 1: Type-I Error Assessment (Normal Distribution, No Mean Difference) ---\n")
result.ER_scenario1 <- sim_ER(N=350, dist="norm", par1 = c(0,0), par2=c(0.5,2), nsim=10^4, measure="sd")

cat("1. Expected Mean Response (EMR):\n")
cat(sprintf("  Overall Mean Response: %.4f\n", mean(result.ER_scenario1[,1]))) # Column 1 is 'response'

cat("\n2. Type-I Error Rate (alpha = 0.05) for Z-Test:\n")
cat(sprintf("  Z-Test Type-I Error Rate: %.4f\n", sum(result.ER_scenario1[,2] < 0.05) / 10^4)) # Column 2 is Z_P_value

cat("\n3. Allocation towards Experimental Arm:\n")
cat(sprintf("  Mean Proportion in Experimental Arm: %.4f\n", mean(result.ER_scenario1[,6]))) # Column 6 is n1/N (proportion in experimental arm)

# Scenario 2: Power Analysis (True difference between arms in means)
# Here, par1 = c(0,0.5) means arm0 has mean 0, arm1 has mean 0.5 (a true difference).
# par2=c(0.5,2) means arm0 has SD 0.5, arm1 has SD 2.
cat("\n--- Scenario 2: Power Assessment (Normal Distribution, Mean Difference Exists) ---\n")
result.ER_scenario2 <- sim_ER(N=350, dist="norm", par1 = c(0,0.5), par2=c(0.5,2), nsim=10^4, measure="sd")

cat("1. Expected Mean Response (EMR):\n")
cat(sprintf("  Overall Mean Response: %.4f\n", mean(result.ER_scenario2[,1])))

cat("\n2. Power (alpha = 0.05) for Z-Test:\n")
cat(sprintf("  Z-Test Power: %.4f\n", sum(result.ER_scenario2[,2] < 0.05) / 10^4))

cat("\n3. Allocation towards Experimental Arm:\n")
cat(sprintf("  Mean Proportion in Experimental Arm: %.4f\n", mean(result.ER_scenario2[,6])))