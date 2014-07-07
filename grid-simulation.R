#gridded simulation with evenly spaced stations
library(asreml)
library(mgcv)
library(lmeSplines)

n.station <- 25 #specify number of stations to run
mu <- 50
sd <- 40
noise.sd <- 0.2 #noise sd for ar process
stn.sd   <- 0.65 #sd for random station effect

mult <- 1e3
z <- seq(0, 250, 5) #explanatory variable (depth)
rho <- mult*dnorm(z, mu, sd)/(pnorm(max(z), mu, sd) - pnorm(min(z), mu, sd))
#z.ar <- arima.sim(n = length(rho), model = list(ar = 0.8), sd = noise.sd) #ar process down z (same for each station)
stn.re <- rnorm(n.station, mean = 0, sd = stn.sd) #station specific random effect

#calculate obs, the observed data, and log(obs)
#l.obs <- rep(log(rho), n.station) + rep(z.ar, n.station) + rep(stn.re, 1, each = length(rho))
l.obs <- NULL

#ar process down z (same for each station)

for (i in c(1:n.station)) {
  
  l.temp <- log(rho) + arima.sim(n = length(rho), model = list(ar = 0.8), sd = noise.sd) +
    rep(stn.re[i],length(rho))
  l.obs <- c(l.obs,l.temp)
}

obs <- exp(l.obs)

#data frame
glm.spl <- data.frame(obs, l.obs, rep(z, n.station), rep(c(1:n.station), 1, each = length(z)), as.factor(z))
names(glm.spl) <- c("obs", "l.obs", "z", "stn", "z.fact")
glm.spl$stn <- as.factor(glm.spl$stn)

#generate locations for grid using 25 stations
x <- sort(rep(sort(rep(seq(1:sqrt(n.station)), length(rho))), sqrt(n.station)))
y <- rep(sort(rep(seq(1:sqrt(n.station)), length(rho))), sqrt(n.station))
grid.dat <- cbind(glm.spl, x, y)

#asreml with ar1 down depths
asreml.fit <- asreml(fixed = l.obs ~ z, random =~ spl(z) + stn + x + y, data = grid.dat, 
                     splinepoints = list(z = seq(0, 250, 25)), rcov =~ stn:ar1(z.fact))
summary(asreml.fit)

summary(asreml.fit)$varcomp[2,2]^0.5
summary(asreml.fit)$varcomp[5,2]^0.5
summary(asreml.fit)$varcomp[6,2]

















