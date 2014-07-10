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


#isotropic gaussic process across x and y
x.cor <- NULL
for (j in 1:(max(y)*length(z))){
  x.temp <- arima.sim(n = max(x), model = list(ar = x.phi), sd = noise.sd/2, innov=rnorm(max(x), mean=0, sd=noise.sd/2))   
  x.cor <- c(x.cor, x.temp)
}

#total observations
x.cor.x <- rep(1:max(x), length(z)*max(y))
x.cor.y <- rep(1:max(y), 1, each = length(z)*max(x))
x.cor <- x.cor[order(x.cor.x, x.cor.y)]
l.obs <- l.obs.z + x.cor
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
                     splinepoints = list(z = seq(0, 250, 25)), rcov =~ ar1(z.fact):y:ar1(x.fact))
summary(asreml.fit)


vals <- matrix(c(stn.sd, noise.sd, phi.true, x.phi, round(summary(asreml.fit)$varcomp[2,2]^0.5, 2), round(summary(asreml.fit)$varcomp[3,2]^0.5, 2), round(summary(asreml.fit)$varcomp[4,2], 2), round(summary(asreml.fit)$varcomp[5,2], 2)), ncol = 2)
colnames(vals) <- c("true", "fitted")
rownames(vals) <- c("stn", "noise", "z ar1", "x ar1")
vals



#variogram
v1 <- Variogram(gam.fit$lme, form =~ z.int|stn, breaks=seq(0,15,1),  robust = TRUE)
plot(v1$dist, v1$variog)
lines(y = (1 - phi^(seq(0, 15,1))), x = seq(0, 15,1), lwd = 2)











