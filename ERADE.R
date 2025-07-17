###############################################################################
library(BSDA)

source('Data_Generator.R')

sim_ERADE <-  function(N=200, dist="bern", par1 = c(0,0), par2 = NULL, measure="sd", burnin=3, ar="RSHIR_Z0",nsim=10^4, alpha=0.05, one.sided = FALSE){
  
    if(dist=="bern"){
      out = matrix(nrow = nsim, ncol = 8)
    }
    else {
      out = matrix(nrow = nsim, ncol = 6)
    }
    for(i in 1:nsim){
      if(i==1){
        out[i,] = two_arm_ERADE(N=N, dist=dist,  par1 =par1, par2 = par2, measure=measure, burnin = burnin,ar=ar, first=FALSE)
      } else {
        out[i,] = two_arm_ERADE(N=N, dist=dist,  par1 =par1, par2 = par2, measure=measure, burnin = burnin,ar=ar)
      }
    }
  return(out)
}

###############################################################################

two_arm_ERADE = function(N=200, dist="bern",  par1 = c(0,0), par2 = NULL, measure="sd", burnin=3, ar="RSHIR_Z0", alphaE=0.5, Z_corrected=FALSE, first=FALSE){
  
  A = rep(NA, N)   # Vector of allocations
  X = rep(NA, N)   # Vector of responses
  est <- rep(NA, N) 
  
  #Burnin
  if(burnin>0){
    A[1:burnin] = 1
    X[1:burnin] = generate_data(N=burnin, dist=dist, parV=c(par1[1],par2[1]))
    A[(burnin+1):(2*burnin)] = 2
    X[(burnin+1):(2*burnin)] = generate_data(N=burnin, dist=dist, parV=c(par1[2],par2[2]))
    est[1:(2*burnin)] <- 0.5
  }
  
  for (i in (2*burnin+1):N){
    n1 = sum(A==1, na.rm = TRUE)
    n2 = sum(A==2, na.rm = TRUE)
    
    if(ar=="RSHIR_Z1"){
      if(dist=="bern") {
        s1 = sum(X[A==1], na.rm = TRUE)
        p1.hat = (s1)/(n1)
        s2 = sum(X[A==2], na.rm = TRUE)
        p2.hat = (s2)/(n2)
        if(measure=="sd"){
          #if(p1.hat*(1-p1.hat)== 0 || p2.hat*(1-p2.hat) == 0){ # Uncomment and comment our line below to sample with equal variacne while one of the variance estimators is equals to zero
          if((sqrt(p1.hat) + sqrt(p2.hat)) == 0){
            rho1.hat = 0.5
          } else {
            rho1.hat = sqrt(p1.hat)/(sqrt(p1.hat) + sqrt(p2.hat))
          }
        } else if(measure=="lrr"){
          if((sqrt(p1.hat)*(1-p2.hat) + sqrt(p2.hat)*(1-p1.hat)) == 0){
            rho1.hat = 0.5
          } else {
            rho1.hat = 1- sqrt(p2.hat)*(1-p1.hat) / (sqrt(p1.hat)*(1-p2.hat) + sqrt(p2.hat)*(1-p1.hat))
          }
        }
        p = 1-rho1.hat
      }
      if(dist=="norm"){
        sig1 <- sqrt(var(X[A==1], na.rm = TRUE))
        sig2 <- sqrt(var(X[A==2], na.rm = TRUE))
        if(sig1+sig2==0){
          rho1.hat = 0.5
        } else {
          p = sig2/(sig1+sig2)
        }
      }
    }
    if(ar=="RSHIR_Z0"){
      if(dist=="bern") {
        s1 = sum(X[A==1], na.rm = TRUE)
        p1.hat = (s1)/(n1)
        s2 = sum(X[A==2], na.rm = TRUE)
        p2.hat = (s2)/(n2)
        if(measure=="sd"){
          if(p1.hat*(1-p1.hat)== 0 || p2.hat*(1-p2.hat) == 0){
            p = 0.5
          } else {
            p = find_root(p1.hat, p2.hat)
          }
        } 
      }
    }
    if(ar=="Neyman_Z0"){
      if(measure=="sd"){
        sig1_N <- sd(X[A==1], na.rm = TRUE)
        sig2_N <- sd(X[A==2], na.rm = TRUE)
        #if(sig1_N == 0 || sig2_N==0){  #### burn-in
          if(sig1_N+sig2_N==0){
          p <- 0.5
        } else {
          p <- sig1_N/(sig1_N+sig2_N)
        }
      }
    }
    if(ar=="Neyman_Z1"){
      if(measure=="sd"){
        sig1_N <- sd(X[A==1], na.rm = TRUE)
        sig2_N <- sd(X[A==2], na.rm = TRUE)
        #if(sig1_N == 0 || sig2_N==0){  #### burn-in
        if(sig1_N+sig2_N==0){
          p <- 0.5
        } else {
          p <- sig2_N/(sig1_N+sig2_N)
        }
      } else if(measure=="lrr"){
        s1 = sum(X[A==1], na.rm = TRUE)
        p1.hat = (s1)/(n1)
        s2 = sum(X[A==2], na.rm = TRUE)
        p2.hat = (s2)/(n2)
        if(sqrt((1-p1.hat)*p2.hat) +  sqrt((1-p2.hat)*p1.hat)==0){
          p <- 0.5
        } else {
          p <- sqrt((1-p1.hat)*p2.hat) / (sqrt((1-p1.hat)*p2.hat) + sqrt((1-p2.hat)*p1.hat))
          
        }
      } else if(measure=="lor"){
        s1 = sum(X[A==1], na.rm = TRUE)
        p1.hat = (s1)/(n1)
        s2 = sum(X[A==2], na.rm = TRUE)
        p2.hat = (s2)/(n2)
        if(sqrt((1-p1.hat)*p1.hat) +  sqrt((1-p2.hat)*p2.hat)==0){
          p <- 0.5
        } else {
          p <- sqrt((1-p1.hat)*p1.hat) / (sqrt((1-p1.hat)*p1.hat) + sqrt((1-p2.hat)*p2.hat))
        }
      }
      est[i] <- p
    }
    
    if(p==0){p <- 1/N}
    if(p==1){p <- 1-1/N}
    if(is.na(p)){
      print(c(p1.hat,p2.hat))
    }
    est[i] <- p
    alloc.prop = n1/i
    
    phi1 = ERADE(1-p, alloc.prop, alphaE) #USE AUC here instead 
    A[i] = rbinom(1, 1, 1-phi1) + 1
    X[i]  <- generate_data(N=1, dist=dist, parV=c(par1[A[i]],par2[A[i]])) 
    
  }
  
  if(first){
    filename <-  paste(ar, N, dist, par1[1], par1[2], "Trial1.pdf",sep="_")
    pdf(filename, height = 4, width = 8)
    # Set plot parameters
    plot_colors <- "#0072B2"  # blue color for line
    title <- "Temporal Behavior of Estimator"  # plot title
    x_label <- "Patient"  # x-axis label
    y_label <- "Estimator"  # y-axis label
    
    # Create plot with custom formatting
    plot(1:N, est, type="l", col=plot_colors, lwd=2, ylim=c(0,1), xlab=x_label, ylab=y_label, main=title)
    grid()  # add gridlines to plot
    box()  # add border to plot
    dev.off() 
  }
  
  # Response
  x0 <- X[A==1]
  x1 <- X[A==2]

  n0 <- length(x0)
  n1 <- length(x1)
  s0 <- sum(x0)
  s1 <- sum(x1)
  response <- (s0+s1)/N
  
  #Imbalance Measure
  sup <- superior(dist=dist,par1 = par1, par2 = par2)
  if(sup==2){ #2 means that no arm is theoretically supirior 
    #   imb <- NA
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
  sigx <- sd(x0); sigy <- sd(x1)
  if(sigx==0 && sigy==0){
    BM_P <- Z_P[1]
  } else {
    BM <- lawstat::brunner.munzel.test(c(x0),c(x1))
    if(BM$statistic == Inf){
      BM_P <- 0
    } else if(BM$statistic == -Inf) {
      BM_P <-  1
    } else {
      BM_P <- BM$p.value
    }
  }
  #estimator <- p
  n <- c(n0,n1)/N
  return(c(response, Z_P, BM_P, per.sup, n))
}


###############################################################################
# Help Functions 

#Function to Calculate ERADE probability based on allcoation and optimal proportion 
ERADE = function(rho.hat, alloc.prop, alpha=0.5){
  if(alloc.prop > rho.hat){
    p = alpha*rho.hat
  } else if(alloc.prop < rho.hat){
    p = 1 - alpha*(1 - rho.hat)
  } else{
    p = rho.hat
  }
  
  return(p)
}


# Function to Calculate RSHIR_Z0
find_root <- function(p0, p1) {
  # Check that inputs are within (0,1)
  if (p0 <= 0 || p0 >= 1 || p1 <= 0 || p1 >= 1) {
    stop("p0 and p1 must be in the range (0,1)")
  }
  
  # Define the equation to solve
  equation <- function(rho) {
    term1 <- (p0 - p1) * ((p0 * (1 - p0 + rho * p0)) / rho + (p1 - rho * p1^2) / (1 - rho) - 2 * p0 * p1)
    term2 <- (1 - p0 + rho * p0 - rho * p1) * ((p1 * (1 - p1)) / (1 - rho)^2 - (p0 * (1 - p0)) / rho^2)
    term1 + term2
  }
  
  # Solve for rho numerically using uniroot
  tryCatch({
    solution <- uniroot(equation, lower = 1e-6, upper = 1 - 1e-6)$root
    return(solution)
  }, error = function(e) {
    message("No solution found within the interval (0,1)")
    return(NA)
  })
}



##### Table 2 Simulations (ERADE) ################

# Scenario 1: Type-I Error (No true difference in means, par1 = c(0,0))
# Arm0 (control): Mean 0, SD 0.5
# Arm1 (experimental): Mean 0, SD 2
cat("--- Adaptive Randomization (ERADE - RSHIR_Z1): Type-I Error Assessment ---\n")
cat(" (Common Mean = 0 for both arms, SDs: Control=0.5, Experimental=2)\n\n")

# Run ERADE for Type-I Error (par1 = c(0,0))
result.ERADE_TypeI <- sim_ERADE(N=350, dist="norm", par1 = c(0,0), par2=c(0.5,2), burnin=35, ar="RSHIR_Z1", nsim=10^4, measure="sd")

cat("1. Expected Mean Response (EMR):\n")
cat(sprintf("  Overall Mean Response: %.4f\n", mean(result.ERADE_TypeI[,1])))

cat("\n2. Type-I Error Rate (alpha = 0.05) for Z-Test:\n")
cat(sprintf("  Z-Test Type-I Error Rate: %.4f\n", sum(result.ERADE_TypeI[,2] < 0.05)/10^4))

cat("\n3. Allocation towards Experimental Arm:\n")
cat(sprintf("  Mean Proportion in Experimental Arm: %.4f\n\n", mean(result.ERADE_TypeI[,6])))

# Scenario 2: Power (True difference in means, par1 = c(0,0.5))
# Arm0 (control): Mean 0, SD 0.5
# Arm1 (experimental): Mean 0.5, SD 2
cat("\n--- Adaptive Randomization (ERADE - RSHIR_Z1): Power Assessment ---\n")
cat(" (Control Mean=0, Experimental Mean=0.5; SDs: Control=0.5, Experimental=2)\n\n")

# Run ERADE for Power (par1 = c(0,0.5))
result.ERADE_Power <- sim_ERADE(N=350, dist="norm", par1 = c(0,0.5), par2=c(0.5,2), burnin=35, ar="RSHIR_Z1", nsim=10^4, measure="sd")

cat("1. Expected Mean Response (EMR):\n")
cat(sprintf("  Overall Mean Response: %.4f\n", mean(result.ERADE_Power[,1])))

cat("\n2. Power (alpha = 0.05) for Z-Test:\n")
cat(sprintf("  Z-Test Power: %.4f\n", sum(result.ERADE_Power[,2] < 0.05)/10^4))

cat("\n3. Allocation towards Experimental Arm:\n")
cat(sprintf("  Mean Proportion in Experimental Arm: %.4f\n\n", mean(result.ERADE_Power[,6])))