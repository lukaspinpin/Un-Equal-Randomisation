library(BSDA)
library(rankFD)

source('Data_Generator.R')

sim_FR <-  function(N=200, dist="bern",  par1 = c(0,0), par2 = NULL, measure="sd", FR=0.5, burnin=3, nsim=10^4, alpha=0.05, one.sided = FALSE){
  
    if(dist=="bern"){
      out = matrix(nrow = nsim, ncol = 8)
    }
    else {
      out = matrix(nrow = nsim, ncol = 6)
    }
    for(i in 1:nsim){
      out[i,] = two_arm_FR(N=N, dist=dist,  par1 =par1, par2 = par2, measure=measure, FR=FR, burnin = burnin, Z_corrected=FALSE)
    }
  return(out)
}

two_arm_FR <- function(N=200, dist="bern",  par1 = c(0,0), par2 = NULL, measure="sd", FR=0.5, burnin=3, Z_corrected=FALSE){
  # Generate Samples 
  n1 = round(N*FR)
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
    Z_P <- wald.test(x0,x1)
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

# --- Simulation Scenarios and Output ---

# Scenario Group 1: Type-I Error (No true difference in means, par1 = c(0,0))
# Arm0 (control): Mean 0, SD 0.5
# Arm1 (experimental): Mean 0, SD 2

cat("--- Fixed Randomization (FR) Simulations: Type-I Error Assessments ---\n")
cat(" (Common Mean = 0 for both arms, SDs: Control=0.5, Experimental=2)\n\n")

# Scenario 1.1: FR = 0.66
FR_val <- 0.66
cat(sprintf("--- Scenario: FR = %.2f (Type-I Error) ---\n", FR_val))
result.FR <- sim_FR(N=350, dist="norm", par1 = c(0,0), par2=c(0.5,2), FR=FR_val, nsim=10^4, measure="sd")
cat("1. Expected Mean Response (EMR):\n")
cat(sprintf("  Overall Mean Response: %.4f\n", mean(result.FR[,1])))
cat("2. Type-I Error Rate (alpha = 0.05) for Z-Test:\n")
cat(sprintf("  Z-Test Type-I Error Rate: %.4f\n", sum(result.FR[,2] < 0.05)/10^4))
cat("3. Allocation towards Experimental Arm:\n")
cat(sprintf("  Mean Proportion in Experimental Arm: %.4f\n\n", mean(result.FR[,6])))

# Scenario 1.2: FR = 0.75
FR_val <- 0.75
cat(sprintf("--- Scenario: FR = %.2f (Type-I Error) ---\n", FR_val))
result.FR <- sim_FR(N=350, dist="norm", par1 = c(0,0), par2=c(0.5,2), FR=FR_val, nsim=10^4, measure="sd")
cat("1. Expected Mean Response (EMR):\n")
cat(sprintf("  Overall Mean Response: %.4f\n", mean(result.FR[,1])))
cat("2. Type-I Error Rate (alpha = 0.05) for Z-Test:\n")
cat(sprintf("  Z-Test Type-I Error Rate: %.4f\n", sum(result.FR[,2] < 0.05)/10^4))
cat("3. Allocation towards Experimental Arm:\n")
cat(sprintf("  Mean Proportion in Experimental Arm: %.4f\n\n", mean(result.FR[,6])))

# Scenario 1.3: FR = 0.80
FR_val <- 0.80
cat(sprintf("--- Scenario: FR = %.2f (Type-I Error) ---\n", FR_val))
result.FR <- sim_FR(N=350, dist="norm", par1 = c(0,0), par2=c(0.5,2), FR=FR_val, nsim=10^4, measure="sd")
cat("1. Expected Mean Response (EMR):\n")
cat(sprintf("  Overall Mean Response: %.4f\n", mean(result.FR[,1])))
cat("2. Type-I Error Rate (alpha = 0.05) for Z-Test:\n")
cat(sprintf("  Z-Test Type-I Error Rate: %.4f\n", sum(result.FR[,2] < 0.05)/10^4))
cat("3. Allocation towards Experimental Arm:\n")
cat(sprintf("  Mean Proportion in Experimental Arm: %.4f\n\n", mean(result.FR[,6])))

# Scenario Group 2: Power (True difference in means, par1 = c(0,0.5))
# Arm0 (control): Mean 0, SD 0.5
# Arm1 (experimental): Mean 0.5, SD 2

cat("\n--- Fixed Randomization (FR) Simulations: Power Assessments ---\n")
cat(" (Control Mean=0, Experimental Mean=0.5; SDs: Control=0.5, Experimental=2)\n\n")

# Scenario 2.1: FR = 0.66
FR_val <- 0.66
cat(sprintf("--- Scenario: FR = %.2f (Power) ---\n", FR_val))
result.FR <- sim_FR(N=350, dist="norm", par1 = c(0,0.5), par2=c(0.5,2), FR=FR_val, nsim=10^4, measure="sd")
cat("1. Expected Mean Response (EMR):\n")
cat(sprintf("  Overall Mean Response: %.4f\n", mean(result.FR[,1])))
cat("2. Power (alpha = 0.05) for Z-Test:\n")
cat(sprintf("  Z-Test Power: %.4f\n", sum(result.FR[,2] < 0.05)/10^4))
cat("3. Allocation towards Experimental Arm:\n")
cat(sprintf("  Mean Proportion in Experimental Arm: %.4f\n\n", mean(result.FR[,6])))

# Scenario 2.2: FR = 0.75
FR_val <- 0.75
cat(sprintf("--- Scenario: FR = %.2f (Power) ---\n", FR_val))
result.FR <- sim_FR(N=350, dist="norm", par1 = c(0,0.5), par2=c(0.5,2), FR=FR_val, nsim=10^4, measure="sd")
cat("1. Expected Mean Response (EMR):\n")
cat(sprintf("  Overall Mean Response: %.4f\n", mean(result.FR[,1])))
cat("2. Power (alpha = 0.05) for Z-Test:\n")
cat(sprintf("  Z-Test Power: %.4f\n", sum(result.FR[,2] < 0.05)/10^4))
cat("3. Allocation towards Experimental Arm:\n")
cat(sprintf("  Mean Proportion in Experimental Arm: %.4f\n\n", mean(result.FR[,6])))

# Scenario 2.3: FR = 0.80
FR_val <- 0.80
cat(sprintf("--- Scenario: FR = %.2f (Power) ---\n", FR_val))
result.FR <- sim_FR(N=350, dist="norm", par1 = c(0,0.5), par2=c(0.5,2), FR=FR_val, nsim=10^4, measure="sd")
cat("1. Expected Mean Response (EMR):\n")
cat(sprintf("  Overall Mean Response: %.4f\n", mean(result.FR[,1])))
cat("2. Power (alpha = 0.05) for Z-Test:\n")
cat(sprintf("  Z-Test Power: %.4f\n", sum(result.FR[,2] < 0.05)/10^4))
cat("3. Allocation towards Experimental Arm:\n")
cat(sprintf("  Mean Proportion in Experimental Arm: %.4f\n\n", mean(result.FR[,6])))