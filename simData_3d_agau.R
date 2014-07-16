simData <- function (n.station, noise.sd, stn.sd, z.phi, x.phi, y.phi) {
  
  
  mu <- 50
  sd <- 40
  mult <- 1e3
  z <- seq(0, 250, 5) #explanatory variable (depth)
  z.int <- rep(c(1:length(z)), n.station) #explanatory variable (depth)
  stn <- rep(c(1:n.station), 1, each = length(z))
  rho <- mult*dnorm(z, mu, sd)/(pnorm(max(z), mu, sd) - pnorm(min(z), mu, sd))
  stn.re <- rnorm(n.station, mean = 0, sd = stn.sd) #station specific random effect
  
  #regular grid for 100 stations
  x <- sort(rep(sort(rep(seq(1:sqrt(n.station)), length(rho))), sqrt(n.station)))
  y <- rep(sort(rep(seq(1:sqrt(n.station)), length(rho))), sqrt(n.station))
  
  #random noise matrix
  r.noise <- rnorm(length(x), 0, noise.sd)
  
  #correlation matrix for x correlation
  cor.ij <- matrix(0, ncol = max(x), nrow = max(x))
  for (i in 1:max(x)) {
    for (j in 1:max(x)) {
      cor.ij[i, j] <- x.phi^((i - j)^2)
    }
  }
  
  #correlation matrix for y correlation
  cor.ij.y <- matrix(0, ncol = max(y), nrow = max(y))
  for (i in 1:max(y)) {
    for (j in 1:max(y)) {
      cor.ij.y[i, j] <- y.phi^((i - j)^2)
    }
  }
  
  #combined ar component
  t.cor <- rep(0, length(r.noise))
  for (k in 2:length(z)) {  
    for (i in 1:max(y)) {
      for (j in 1:max(x)) {
        w <- which(y == i & z.int == k & x == j)
        x.cor <- matrix(rep(cor.ij[, j], max(y)), ncol = max(x), byrow = T) #correlation matrix 
        y.cor <- matrix(rep(cor.ij.y[, i], max(x)), ncol = max(x))
        z.data <- matrix(r.noise[z.int == k], ncol = max(x))
        t.cor[w] <- t.cor[y == i & z.int == (k - 1) & x == j]*z.phi + sum((x.cor*y.cor)*z.data) + rnorm(1, 0, noise.sd/2)
      }
    }
  }
  
  
  #total observations
  l.obs <- rep(log(rho), n.station) + t.cor + rep(stn.re, 1, each = length(rho))
  obs <- exp(l.obs)
  
  #data frame
  glm.spl <- data.frame(obs, l.obs, rep(z, n.station), as.factor(rep(c(1:n.station), 1, each = length(z))), x, y)
  names(glm.spl) <- c("obs", "l.obs", "z", "stn", "x", "y")
  glm.spl$z.fact <- as.factor(as.integer(glm.spl$z))
  glm.spl$x.fact <- as.factor(as.integer(glm.spl$x))
  glm.spl$y.fact <- as.factor(as.integer(glm.spl$y))
  glm.spl <- glm.spl[order(glm.spl$z, glm.spl$x), ] #sort by order of rcov structure
  
  return(glm.spl)
  
}