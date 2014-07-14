#simulates a regular grid with anisotropic correlation in the x and y directions
library(asreml)
library(mgcv)

n.station <- 100 #specify number of stations to run
mu <- 50
sd <- 40
noise.sd <- 0.2 #noise sd for ar process
stn.sd   <- 0.65 #sd for random station effect
phi.true <- 0.6 #autocorrelation (ar1) down depths
x.phi <- 0.6 #autocorrelation across x
y.phi <- 0.6 #autocorrelation across y

mult <- 1e3
z <- seq(0, 250, 5) #explanatory variable (depth)
z.int <- rep(c(1:length(z)), n.station) #explanatory variable (depth)
stn <- rep(c(1:n.station), 1, each = length(z))
rho <- mult*dnorm(z, mu, sd)/(pnorm(max(z), mu, sd) - pnorm(min(z), mu, sd))
stn.re <- rnorm(n.station, mean = 0, sd = stn.sd) #station specific random effect

#regular grid for 100 stations
x <- sort(rep(sort(rep(seq(1:sqrt(n.station)), length(rho))), sqrt(n.station)))
y <- rep(sort(rep(seq(1:sqrt(n.station)), length(rho))), sqrt(n.station))

#white noise
r.noise <- rnorm(length(rho)*n.station, 0, noise.sd)


#ar process down z 
l.obs.z <- NULL
for (i in 1:n.station) {
  for (j in 2:length(z)){
    l.temp[j] <- r.noise[stn == i & z.int == (j - 1)]*phi.true
  }
  l.obs.z <- c(l.obs.z,l.temp)
}


#ar process across x 
x.cor <- NULL
for (i in 1:max(y)) {
  for (k in 1:length(z)) {  
    x.temp <- rep(0, max(x))
    for (j in 2:max(x)) {
      x.temp[j] <- r.noise[y == i & z.int == k & x == (j - 1)]*phi.true 
    }
    x.cor <- c(x.cor,x.temp)
  }
}

#ar process across y 
y.cor <- NULL
for (i in 1:max(x)) {
  for (k in 1:length(z)) {  
    y.temp <- rep(0, max(y))
    for (j in 2:max(y)) {
      y.temp[j] <- r.noise[x == i & z.int == k & y == (j - 1)]*phi.true 
    }
    y.cor <- c(y.cor, y.temp)
  }
}

#total observations
x.cor.x <- rep(1:max(x), length(z)*max(y))
x.cor.y <- rep(1:max(y), 1, each = length(z)*max(x))
x.cor <- x.cor[order(x.cor.x, x.cor.y)]
y.cor.x <- rep(1:max(x), 1, each = length(z)*max(y))
y.cor.y <- rep(1:max(y), length(z)*max(x))
y.cor <- y.cor[order(y.cor.x, y.cor.y)]
l.obs <- rep(log(rho), n.station) + l.obs.z + x.cor + y.cor + r.noise+ rep(stn.re, 1, each = length(rho))
obs <- exp(l.obs)

#data frame
glm.spl <- data.frame(obs, l.obs, rep(z, n.station), rep(c(1:n.station), 1, each = length(z)), x, y)
names(glm.spl) <- c("obs", "l.obs", "z", "stn", "x", "y")
glm.spl$stn <- as.factor(glm.spl$stn)
glm.spl$z.fact <- as.factor(as.integer(glm.spl$z))
glm.spl$x.fact <- as.factor(as.integer(glm.spl$x))
glm.spl$y.fact <- as.factor(as.integer(glm.spl$y))

glm.spl <- glm.spl[order(glm.spl$z, glm.spl$y), ]


#---------------------------------- asreml -------------------------------------#

asreml.fit <- asreml(fixed = l.obs ~ z, random =~ spl(z) + stn, data = glm.spl, 
                     splinepoints = list(z = seq(0, 250, 25)), rcov =~ ar1(z.fact):ar1(y.fact):ar1(x.fact))
summary(asreml.fit)


vals <- matrix(c(stn.sd, noise.sd, phi.true, x.phi, y.phi, round(summary(asreml.fit)$varcomp[2,2]^0.5, 2), round(summary(asreml.fit)$varcomp[3,2]^0.5, 2), round(summary(asreml.fit)$varcomp[4,2], 2), round(summary(asreml.fit)$varcomp[5,2], 2), round(summary(asreml.fit)$varcomp[6,2], 2)), ncol = 2)
colnames(vals) <- c("true", "fitted")
rownames(vals) <- c("stn", "noise", "z ar1", "y ar1", "x ar1")
vals













