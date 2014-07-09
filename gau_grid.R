#simulates a regular grid with anisotropic correlation in the x and y directions
library(asreml)
library(mgcv)

n.station <- 100 #specify number of stations to run
mu <- 50
sd <- 40
noise.sd <- 0.2 #noise sd for ar process
stn.sd   <- 0.65 #sd for random station effect
phi.true <- 0.6 #autocorrelation (ar1) down depths
x.phi <- 0.5 #autocorrelation across x
y.phi <- 0.5 #autocorrelation across y

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


#correlation matrix for x correlation
cor.ij <- matrix(0, ncol = max(x), nrow = max(x))
for (i in 1:max(x)) {
  for (j in 1:max(x)) {
    cor.ij[i, j] <- x.phi^((i - j)^2)
  }
}

#correlation matrix for y correlation
cory.ij <- matrix(0, ncol = max(x), nrow = max(x))
for (i in 1:max(x)) {
  for (j in 1:max(x)) {
    cory.ij[i, j] <- y.phi^((i - j)^2)
  }
}

#anisotropic gaussic process across x and y
start.vals <- matrix(rnorm(max(x)*max(y), 0, noise.sd/2), ncol = max(x))
x.cor <- rep(0, max(x)) 

for (j in c(1:max(x))){
  for (i in 1:max(y)) {
    for (k in  c(1:max(x)){
      for (l in  c(1:max(x)){
        x.cor[j] <- start.vals[j] + start.vals[i]*cor.ij[i, j] * start.vals[i]*cory.ij[i, j]       
      }
    }  
  }
}



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

glm.spl <- glm.spl[order(glm.spl$z, glm.spl$y), ]


#---------------------------------- asreml -------------------------------------#

asreml.fit <- asreml(fixed = l.obs ~ z, random =~ spl(z) + stn, data = glm.spl, 
                     splinepoints = list(z = seq(0, 250, 25)), rcov =~ ar1(z.fact):y:ar1(x.fact))
summary(asreml.fit)


vals <- matrix(c(stn.sd, noise.sd, phi.true, x.phi, round(summary(asreml.fit)$varcomp[2,2]^0.5, 2), round(summary(asreml.fit)$varcomp[3,2]^0.5, 2), round(summary(asreml.fit)$varcomp[4,2], 2), round(summary(asreml.fit)$varcomp[5,2], 2)), ncol = 2)
colnames(vals) <- c("true", "fitted")
rownames(vals) <- c("stn", "noise", "z ar1", "x car1")
vals

#variogram
v1 <- Variogram(gam.fit$lme, form =~ z.int|stn, breaks=seq(0,15,1),  robust = TRUE)
plot(v1$dist, v1$variog)
lines(y = (1 - phi^(seq(0, 15,1))), x = seq(0, 15,1), lwd = 2)











