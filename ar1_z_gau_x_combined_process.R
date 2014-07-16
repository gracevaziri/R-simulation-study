#simulates a regular grid with anisotropic correlation in the x and y directions
library(asreml)
library(mgcv)

n.station <- 100 #specify number of stations to run
mu <- 50
sd <- 40
noise.sd <- 0.2 #noise sd for ar process
stn.sd   <- 0.65 #sd for random station effect
phi.true <- 0.55 #autocorrelation (ar1) down depths
x.phi <- 0.65 #autocorrelation across x

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

#combined ar component
t.cor <- rep(0, length(r.noise))
for (i in 1:max(y)) {
  for (k in 2:length(z)) {  
    for (j in 2:(max(x) - 1)) {
      w <- which(y == i & z.int == k & x == j)
      t.cor[w] <- t.cor[y == i & z.int == (k - 1) & x == j]*phi.true + sum(r.noise[y == i & z.int == k]*cor.ij[c(1:max(x)), j])  
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
glm.spl <- glm.spl[order(glm.spl$z, glm.spl$y), ] #sort by order of rcov structure


#---------------------------------- asreml -------------------------------------#

asreml.fit <- asreml(fixed = l.obs ~ z, random =~ spl(z) + stn, data = glm.spl, 
                     splinepoints = list(z = seq(0, 250, 25)), rcov =~ ar1(z.fact):y:gau(x.fact))
summary(asreml.fit)


vals <- matrix(c(stn.sd, noise.sd, phi.true, x.phi, round(summary(asreml.fit)$varcomp[2,2]^0.5, 2), round(summary(asreml.fit)$varcomp[3,2]^0.5, 2), round(summary(asreml.fit)$varcomp[4,2], 2), round(summary(asreml.fit)$varcomp[5,2], 2)), ncol = 2)
colnames(vals) <- c("true", "fitted")
rownames(vals) <- c("stn", "noise", "z ar1", "x gau")
vals










