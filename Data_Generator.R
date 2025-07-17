generate_data <- function(N=200, dist="bern", parV = 0.5){
  # Function to generate data based on specified distribution
  # Args:
  #   N: Number of samples to generate (default = 200)
  #   dist: Distribution type ("bern", "norm", "Lickert", "Expon", "Beta")
  #   parV: Parameters for the chosen distribution
  
  if(dist=="bern"){
    # Generate Bernoulli distributed data with probability parV
    return(rbinom(N,1, parV))
  }
  if(dist=="norm"){
    # Generate Normal distributed data with given mean and standard deviation
    return(rnorm(n=N, mean=parV[1], sd=parV[2]))
  }
  if(dist=="Lickert"){
    # Placeholder for Likert-scale data generation (not implemented)
    return(0)
  }
  if(dist=="Expon"){
    # Placeholder for Exponential distribution data generation (not implemented)
    return(0)
  }
  if(dist=="Beta"){
    # Placeholder for Beta distribution data generation (not implemented)
    return(0)
  }
}

superior <- function(dist="bern", par1=c(0,0), par2=FALSE){
  # Function to identify the superior arm given a distribution and parameter sets
  # Args:
  #   dist: Distribution type ("bern", "norm", "Lickert", "Beta", "Expon")
  #   par1: First set of parameters (vector)
  #   par2: Second set of parameters (vector or FALSE if not applicable)
  # Returns:
  #   1 if experimental arm is superior, 0 if control is superior, 2 if equal
  
  if(dist=="bern"){
    if(par1[1] < par1[2]){
      return(1)
    }
    if(par1[1] > par1[2]){
      return(0)
    }
    return(2)
  }
  if(dist=="norm"){
    if(par1[1] < par1[2]){
      return(1)
    }
    if(par1[1] > par1[2]){
      return(0)
    }
    return(2)
  }
  if(dist=="Lickert" || dist=="Beta"){
    # Compute expected values for comparison
    exp1 <- (par1[1] + par2[1]) / par1[1]
    exp2 <- (par1[2] + par2[2]) / par1[2]
    
    if(exp1[1] > exp1[2]){
      return(1)
    }
    if(exp1[1] < exp1[2]){
      return(0)
    }
    return(2)
  }
  if(dist=="Expon"){
    if(par1[1] > par1[2]){
      return(1)
    }
    if(par1[1] < par1[2]){
      return(0)
    }
    return(2)
  }
}
