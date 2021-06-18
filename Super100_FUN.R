CI.GLM <- function(event, Age, expo){
  
  fit <- glm(event ~ -1 + as.factor(Age), offset = log(expo),
             family=poisson)
  
  V <- vcov(fit)
  gr <- diag(exp(fit$coef))
  tgVg <- t(gr) %*% V %*% gr
  se.h <- sqrt(diag(tgVg))
  h.hat <- exp(fit$coeff)
  h.up <- h.hat + 1.96*se.h
  h.low <- h.hat - 1.96*se.h
  
  GLM<- data.frame(h.hat, h.up, h.low, Age)
  
  return(GLM)
  
}


# Calculate hazard non- parametrically
mx.calc <- function(entry, exit, delta, interval, s.age, nx){

dati <- pyears(Surv(entry, exit, delta) ~ interval,scale=1)
event <- dati$event
expo  <- dati$pyears
mx    <- event / expo 
Age <- seq(s.age, (length(mx)*nx-nx)+s.age, by = nx)
out <- data.frame(Age, event, expo, mx)

return(out)}

# Calculate confidence intervals



CImx <- function(Age, mx, nx, n=20000, nsim=10000){
  
  # Transform rates into probabilties
  # Assume constant rate within each age interval
  qx <- 1-exp(-mx * nx)# Annual
  
  # Output table with observed probabilities
  qxMat <- matrix(0, nrow = nsim, ncol = length(qx))
  
  # Run simulation for annual probabilities
  for (i in 1:nsim) {
    # Initial population
    n1 <- n
    # Simulate number of deaths in each age interval
    for (k in 1:length(qx)) {
      # If there are still individuals alive
      if (n1 > 0) {
        # Assign a random uniform number between 0 and 1 to each individual
        pop <- runif(n1)
        # Number of deaths
        die <- length(which(pop < qx[k]))
        # Store observed probability
        qxMat[i, k] <- die/n1
        # Survivors to the next age
        n1 <- n1 - die
      }
      
    }
  }
  
  # Back transform into quarterly rates
  mxMat <- -log(1-qxMat) / nx 
  
  # 95% Confidence intervals: The values may differ slightly at each simulation
  M025 <- apply(mxMat, 2, quantile, probs = 0.025, na.rm = T)
  M975 <- apply(mxMat, 2, quantile, probs = 0.975, na.rm = T)
  
  # CI for annual probabilities
  qx2 <- 1- exp(-mx)
  qxMat2 <- 1- exp(-mxMat)
  
  Q025 <-  apply(qxMat2, 2, quantile, probs = 0.025, na.rm = T)
  Q975 <-  apply(qxMat2, 2, quantile, probs = 0.975, na.rm = T)
  
  # some other measures
  
  Hx <- cumsum(mx)
  lx <- exp(-Hx)
  
  out <- data.frame(Age, mx, mx.low = M025,
                    mx.up = M975, qx = qx2, qx.low=Q025,
                    qx.up=Q975, Hx, lx)
  
  return(out)}




CIqx <- function(Age, qx, nx=0.25, n=20000, nsim=10000){
  
  # Transform rates into probabilties
  # Assume constant rate within each age interval
 
  # Output table with observed probabilities
  qxMat <- matrix(0, nrow = nsim, ncol = length(qx))
  
  # Run simulation for annual probabilities
  for (i in 1:nsim) {
    # Initial population
    n1 <- n
    # Simulate number of deaths in each age interval
    for (k in 1:length(qx)) {
      # If there are still individuals alive
      if (n1 > 0) {
        # Assign a random uniform number between 0 and 1 to each individual
        pop <- runif(n1)
        # Number of deaths
        die <- length(which(pop < qx[k]))
        # Store observed probability
        qxMat[i, k] <- die/n1
        # Survivors to the next age
        n1 <- n1 - die
      }
      
    }
  }
  
  # Back transform into quarterly rates
  mx <- -log(1-qx)
  mxMat <- -log(1-qxMat)
  # 95% Confidence intervals: The values may differ slightly at each simulation
  M025 <- apply(mxMat, 2, quantile, probs = 0.025, na.rm = T)
  M975 <- apply(mxMat, 2, quantile, probs = 0.975, na.rm = T)
  
  # CI for annual probabilities
  qx2 <- 1- exp(-mx)
  qxMat2 <- 1- exp(-mxMat)
  
  Q025 <-  apply(qxMat2, 2, quantile, probs = 0.025, na.rm = T)
  Q975 <-  apply(qxMat2, 2, quantile, probs = 0.975, na.rm = T)
  
  # some other measures
  
  Hx <- cumsum(mx)
  lx <- exp(-Hx)
  
  out <- data.frame(Age, mx, mx.low = M025,
                    mx.up = M975, qx = qx2, qx.low=Q025,
                    qx.up=Q975, Hx, lx)
  
  return(out)}



# Break point test
breaks <-function(mx, Age, method = "PELT", penalty = "AIC"){
points<- cpt.mean(mx, method = method, penalty = penalty)@cpts
breaks <- Age[points]
out <- data.frame(breaks)
return(out)
}

# Rate of Ageing
bx <- function(Rx, Rx.low, Rx.up, Age, nrol, fn = mean, skip, nx){
  rol     <-rollapply(Rx, nrol, fn, by = skip)
  rol.low <-rollapply(Rx.low, nrol, fn, by = skip)
  rol.up  <-rollapply(Rx.up,  nrol, fn, by = skip)
  bx <-diff(rol) / nx
  bx.low <-diff(rol.low) / nx
  bx.up <-diff(rol.up) / nx
  Age <- Age[1:length(bx)]
  out<- data.frame(Age, bx, bx.low, bx.up)
  return(out)}


# SIMULATE CI for USA
simUSA <- function(mx, n, nx, nage, nsimCI = 100) {
  
  # mx      Vector of death rates
  # n       Initial sample size
  # nx      year intervals
  # nsimCI  Number of simulations to calculate Confidence intervals

  # Transfrom RATES into PROBABILITIES
  qx <- 1-exp(-mx * nx)
  
  # Matrix to store observed probabilities
  matQx <- c()
  
  # Estimate confidence intervals
  for (j in 1:nsimCI) {
    
    # Vector to store observed probabilities
    obsQx <- c()
    
    # Initial sample size
    n0 <- n
    
    # Run simulation over age groups
    for (k in 1:length(nage)) {
      # If there are still individuals alive
      if (n0 > 0) {
        # Assign a random uniform number between 0 and 1 to each individual
        pop <- runif(n0)
        # Number of deaths
        die <- length(which(pop < qx[k]))
        # Store observed probability
        obsQx <- c(obsQx, die/n0)
        # Survivors to the next age
        n0 <- n0 - die
      } else obsQx <- c(obsQx, 0)
    }
    
    # Fill matrix with results
    matQx <- rbind(matQx, obsQx)
  }
  
  # Return matrix with results of the simulations
  rownames(matQx) <- NULL
  return(matQx)
}

# Confidence intervals fro USA
ciUSA <- function(qxMat, Age, nx = 0.25) {
  
  # Back transform PROBABILITIES into quarterly RATES
  mxMat <- -log(1-qxMat) / nx 
  
  # 95% Confidence intervals: The values may differ slightly at each simulation
  M025 <- apply(mxMat, 2, quantile, probs = 0.025, na.rm = T)
  M50 <- apply(mxMat, 2, median, na.rm = T)
  Mmean <- apply(mxMat, 2, 
                 function(x) {
                   id <- which(x > quantile(x, .005, na.rm = T) &
                                 x < quantile(x, .995, na.rm = T))
                   x <- mean(x[id], na.rm = T)
                   if (is.na(x)) x <- 0
                   return(x)
                 })
  M975 <- apply(mxMat, 2, quantile, probs = 0.975, na.rm = T)
  
  # Some other measures
  Hx <- cumsum(M50)
  Sx <- exp(-Hx)
  
  # Output
  out <- data.frame(Age, mxMean = Mmean, mxMedian = M50, mx.low = M025,
                    mx.up = M975, Hx = Hx, Sx = Sx)
  return(out)
  
}

