simData <- function (n.station, noise.sd, stn.sd, z.phi, x.phi) {
  
  
  mu <- 50
  sd <- 40
  mult <- 1e3
  z <- seq(0, 250, 5) #explanatory variable (depth)
  rho <- mult*dnorm(z, mu, sd)/(pnorm(max(z), mu, sd) - pnorm(min(z), mu, sd))
  stn.re <- rnorm(n.station, mean = 0, sd = stn.sd) #station specific random effect
  
  #regular grid for 100 stations
  x <- sort(rep(sort(rep(seq(1:sqrt(n.station)), length(rho))), sqrt(n.station)))
  y <- rep(sort(rep(seq(1:sqrt(n.station)), length(rho))), sqrt(n.station))
  
  #ar process down z 
  l.obs.z <- NULL
  for (i in c(1:n.station)) {
    l.temp <- log(rho) + arima.sim(n = length(rho), model = list(ar = z.phi), sd = noise.sd) 
    l.obs.z <- c(l.obs.z,l.temp)
  }
  
  
  # isotropic ar1 process across x (for each y and z combination)
  # additive on log scale
  
  x.cor <- NULL
  
  for (j in 1:max(y)){
    x.temp <- arima.sim(n = max(x), model = list(ar = x.phi), sd = stn.sd)   
    x.cor <- c(x.cor, x.temp)
  }
  
  x.cor.v <- rep(x.cor, each=length(z))
  x.cor.x <- rep(1:max(x), each=length(z), times=max(y))
  
  x.cor.y <- rep(1:max(y), 1, each = length(z)*max(x))
  x.cor <- x.cor.v[order(x.cor.x, x.cor.y)] #sort to be the same as the rest of the data set
  
  #total observations
  l.obs <- l.obs.z + x.cor
  obs <- exp(l.obs)
  
  #data frame
  glm.spl <- data.frame(obs, l.obs, rep(z, n.station), rep(c(1:n.station), 1, each = length(z)), x, y)
  names(glm.spl) <- c("obs", "l.obs", "z", "stn", "x", "y")
  glm.spl$stn <- as.factor(glm.spl$stn)
  glm.spl$z.fact <- as.factor(as.integer(glm.spl$z))
  glm.spl$x.fact <- as.factor(as.integer(glm.spl$x))
  glm.spl$y.fact <- as.factor(as.integer(glm.spl$y))
  
  glm.spl <- glm.spl[order(glm.spl$y,glm.spl$x,glm.spl$z), ]
  
  return(glm.spl)
  
}