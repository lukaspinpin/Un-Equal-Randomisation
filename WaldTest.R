wald.test.binary <- function(x0,x1, measure="sd", signlevel=0.05){
  # Function for Binary Wald test
  
  # Estimators
  # MLEs
  p0 <- mean(x0); p1 <- mean(x1)
  q0 <- 1-p0 ; q1 <- 1-p1
  sdx0 <- p0*q0; sdx1 <- p1*(1-p1)
  n0 <- length(x0); n1 <- length(x1); n <- n0+n1
  
  # For pooled Variance
  pbar <- (p0*n0+p1*n1)/n; qbar <- 1-pbar
  
  # Z1 is Wald test, Z2 is score test and Z3 is Wald with Agresti & Caffo Correction
  if(measure=="sd"){
    if(pbar*qbar==0){
      if(p0==p1){Z2 = 0}
      if(p0<p1){Z2 = -Inf}
      if(p0>p1){Z2 = Inf}
    }
    else{
      Z2 <- (p0-p1)/sqrt(n*pbar*qbar/(n0*n1)) 
    }
    if((sdx0==0 && sdx1==0)){
      if(p0==p1){Z = 0}
      if(p0<p1){Z = -Inf}
      if(p0>p1){Z = Inf}
    } else {
      Z <- (p0-p1)/sqrt(sdx0/n0+sdx1/n1)
    }
    # Agresti and Caffo adjusted estimators
    p0 <- mean(c(x0,0,1)); p1 <- mean(c(x1,0,1))
    q0 <- 1-p0 ; q1 <- 1-p1
    sdx0 <- p0*q0; sdx1 <- p1*(1-p1)
    n0 <- length(c(x0,0,1)); n1 <- length(c(x1,0,1)); n <- n0+n1
    if((sdx0==0 && sdx1==0)){
      if(p0==p1){Z3 = 0}
      if(p0<p1){Z3 = -Inf}
      if(p0>p1){Z3 = Inf}
    } else {
      Z3 <- (p0-p1)/sqrt(sdx0/n0+sdx1/n1)
    }
  } else if(measure=="rr") {
    if(p0==1 || (p0*q1==0 && p1*q1==0) ){ #added (p0==1 || p1==1) 
      if(p0==p1){Z = 0}
      if(p0<p1){Z = -Inf}
      if(p0>p1){Z = Inf}
    } else {
      Z <- (q1/q0-1)/sqrt((p0*q1**2)/(n0*q0**3)+p1*q1/(n1*q0**2))
    }
  } else if(measure=="or"){
    if( p0==1 || p0==0 || p1==1 || p1==0){   #odds ratio changed (p0==1 || p1==1) instead of (p0==1 && p1==1)
      if(p0==p1){Z = 0}
      if(p0<p1){Z = -Inf}
      if(p0>p1){Z = Inf}
    }
    else{
      Z <- (p0*q1/(q0*p1)-1)/sqrt((p0*q1**2)/(n0*p1**2*q0**3)+(p0**2*q1)/(n1*p1**3*q0**2))
    }
  } else if(measure=="lrr"){
    if(p0==1 || p1==1 || (p0==0 && p1==0) ){
      if(p0==p1){Z = 0}
      if(p0<p1){Z = -Inf}
      if(p0>p1){Z = Inf}
    } else {
      Z <- log(q1/q0)/ sqrt(p0/(n0*q0)+p1/(n1*q1))
    }
  } else if(measure=="lor"){
    if(p0==0 || q0==0 || p1==0 || q1==0){
      if(p0==p1){Z = 0}
      if(p0<p1){Z = -Inf}
      if(p0>p1){Z = Inf}
    } else {
      Z <- log(p0*q1/(q0*p1))/ sqrt(1/(n0*p0*q0)+1/(n1*p1*q1))
    }
  }
  return(c(2*pnorm(-1*abs(Z)),2*pnorm(-1*abs(Z3)),2*pnorm(-1*abs(Z2))))
}

# Non-binary Wald test
wald.test <- function(x0,x1, measure="sd"){
  p0 <- mean(x0); p1 <- mean(x1)
  varx0 <- var(x0); varx1 <- var(x1)
  n0 <- length(x0); n1 <- length(x1)
  
  
  if((varx0==0 && varx1==0)){
    if(p0==p1){Z = 0}
    if(p0<p1){Z = -Inf}
    if(p0>p1){Z = Inf}
  } else {
    if(measure=="sd"){
      Z <- (p0-p1)/sqrt(varx0/n0+varx1/n1)
    }
  }
  return(2*pnorm(-1*abs(Z)))
}