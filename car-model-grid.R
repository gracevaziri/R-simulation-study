#simulates a regular grid with CAR model along x direction
library(asreml)
library(mgcv)

n.station <- 100 #specify number of stations to run
mu <- 50
sd <- 40
noise.sd <- 0.2 #noise sd for ar process
stn.sd   <- 0.65 #sd for random station effect
phi.true <- 0.6 #autocorrelation (ar1) down depths
x.phi <- 0.5 #autocorrelation (car1) across x

mult <- 1e3
z <- seq(0, 250, 5) #explanatory variable (depth)
rho <- mult*dnorm(z, mu, sd)/(pnorm(max(z), mu, sd) - pnorm(min(z), mu, sd))
stn.re <- rnorm(n.station, mean = 0, sd = stn.sd) #station specific random effect

#regular grid for 100 stations
x <- sort(rep(sort(rep(seq(1:sqrt(n.station)), length(rho))), sqrt(n.station)))
y <- rep(sort(rep(seq(1:sqrt(n.station)), length(rho))), sqrt(n.station))

#ar process down z (same for each station)
l.obs.z <- NULL
for (i in c(1:n.station)) {
  l.temp <- log(rho) + arima.sim(n = length(rho), model = list(ar = phi.true), sd = noise.sd/2, 
                                 innov=rnorm(length(rho), mean=0, sd=noise.sd/2)) +  rep(stn.re[i],length(rho)) 
  l.obs.z <- c(l.obs.z,l.temp)
}

#car process across x (same for each depth*y combination)
start.vals <- rnorm(max(x), 0, noise.sd/2) 
x.cor <- rep(0, max(x)) 

for (j in c(2:(max(x) - 1))){
  x.cor[j] <- start.vals[j] + start.vals[j - 1]*(x.phi/2) +  start.vals[j + 1]*(x.phi/2)
}
x.cor[1] <- start.vals[1] + start.vals[2]*(x.phi)
x.cor[max(x)] <- start.vals[max(x)] + start.vals[max(x) - 1]*(x.phi)

#total observations
l.obs <- l.obs.z + rep(x.cor, max(y)*length(z))

obs <- exp(l.obs)

#data frame
glm.spl <- data.frame(obs, l.obs, rep(z, n.station), rep(c(1:n.station), 1, each = length(z)), x, y)
names(glm.spl) <- c("obs", "l.obs", "z", "stn", "x", "y")
glm.spl$stn <- as.factor(glm.spl$stn)
glm.spl$z.fact <- as.factor(as.integer(glm.spl$z))
glm.spl$x.fact <- as.factor(as.integer(glm.spl$x))
glm.spl$y.fact <- as.factor(as.integer(glm.spl$y))

#---------------------------------- asreml -------------------------------------#

asreml.fit <- asreml(fixed = l.obs ~ z, random =~ spl(z) + stn, data = glm.spl, 
                     splinepoints = list(z = seq(0, 250, 25)), rcov =~ stn:ar1(z.fact) + sar(x.fact))
asreml.fit <- update(asreml.fit)
summary(asreml.fit)


vals <- matrix(c(stn.sd, noise.sd, phi.true, x.phi, summary(asreml.fit)$varcomp[2,2]^0.5, summary(asreml.fit)$varcomp[5,2]^0.5, summary(asreml.fit)$varcomp[6,2], summary(asreml.fit)$varcomp[4,2]), ncol = 2)
colnames(vals) <- c("true", "fitted")
rownames(vals) <- c("stn", "noise", "z ar1", "x car1")
vals

#variogram
v1 <- Variogram(gam.fit$lme, form =~ z.int|stn, breaks=seq(0,15,1),  robust = TRUE)
plot(v1$dist, v1$variog)
lines(y = (1 - phi^(seq(0, 15,1))), x = seq(0, 15,1), lwd = 2)






