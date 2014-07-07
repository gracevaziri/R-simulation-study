library(mgcv)


n.station <- 20 #specify number of stations to run
mu <- 50
sd <- 40
noise.sd <- 0.2 #sd of white noise
stn.sd <- rnorm(n.station, 0.5, 0.2) #sd for random station effect


mult <- 1e3
z <- seq(0, 250, 5) #explanatory variable (depth)
rho <- mult*dnorm(z, mu, sd)/(pnorm(max(z), mu, sd) - pnorm(min(z), mu, sd))
noise <- rnorm(n = length(rho), mean = 0, sd = noise.sd) #white noise
stn.re <- rnorm(length(stn.sd), mean = 0, sd = stn.sd) #station specific random effect

#calculate obs, the observed data, and log(obs)
obs <- rep(rho, n.station) + rep(noise, n.station) + rep(stn.re, 1, each = length(rho))
l.obs <- log(obs + 2)

#fit gamm with station random effect
glm.spl <- data.frame(obs, l.obs, rep(z, n.station), rep(c(1:n.station), 1, each = length(z)))
names(glm.spl) <- c("obs", "l.obs", "z", "stn")
glm.spl$stn <- as.factor(glm.spl$stn)
lme.fit <- gamm(l.obs ~ s(z),  random=list(stn =~1), data = glm.spl)
summary(lme.fit$lme)

#plot fitted values and observed data
plot(log(rho + 2), z, ylim = c(max(z), min(z)), type = 'l', xlim = c(min(na.omit(l.obs)), max(na.omit(l.obs))), xlab = "l.fluoro")
points(l.obs, rep(z, n.station), col = glm.spl$stn)
title("gamm with station random effect")
lines(fitted(lme.fit$gam)[glm.spl$stn == 1], col = "black", z, lwd = 2)

#variogram
gamm.var <- Variogram(lme.fit$lme)
plot(gamm.var)
