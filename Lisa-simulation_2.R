library(mgcv)

n.station <- 25 #specify number of stations to run
mu <- 50
sd <- 40
noise.sd <- 0.2 #sd of white noise

#  random effects variance can be set to be constant (SGC)
#  too difficult an estimation problem if each station re comes from a different distribution
#  unless you can standardise by a transformation such as log for lognormal re distribution
#  do this explicitly using HGLMs or implicitly (for just lognormal) using GLMM or LMM for log fluoro

#stn.sd <- rnorm(n.station, 0.5, 0.2) #sd for random station effect
stn.sd <- 0.2 #sd for random station effect

summary(stn.sd)

mult <- 1e3
z <- seq(0, 250, 5) #explanatory variable (depth)
rho <- mult*dnorm(z, mu, sd)/(pnorm(max(z), mu, sd) - pnorm(min(z), mu, sd))
noise <- rnorm(n = length(rho), mean = 0, sd = noise.sd) #white noise
stn.re <- rnorm(n.station, mean = 0, sd = stn.sd) #station specific random effect

#calculate obs, the observed data, and log(obs)
l.obs <- rep(log(rho), n.station) + rep(noise, n.station) + rep(stn.re, 1, each = length(rho))
obs <- exp(l.obs)

#fit gamm with station random effect
glm.spl <- data.frame(obs, l.obs, rep(z, n.station), rep(c(1:n.station), 1, each = length(z)))
names(glm.spl) <- c("obs", "l.obs", "z", "stn")
glm.spl$stn <- as.factor(glm.spl$stn)
gam.fit <- gamm(l.obs ~ s(z),  random=list(stn =~1), data = glm.spl)

summary(gam.fit$gam)
summary(gam.fit$lme)

lme.fit <- gamm(l.obs ~ s(z), random=list(stn =~1), data = glm.spl)


lme.fit.CAR <- gamm(l.obs ~ s(z, bs = "cr"),  random=list(stn =~1), correlation=corCAR1(form= ~ z | stn), data = glm.spl)

summary(lme.fit.CAR$gam)
summary(lme.fit.CAR$lme)


#plot fitted values and observed data
plot(log(rho), z, ylim = c(max(z), min(z)), type = 'l', xlim = c(min(na.omit(l.obs)), max(na.omit(l.obs))), xlab = "l.fluoro")
points(l.obs, rep(z, n.station), col = glm.spl$stn)
title("gamm with station random effect")
lines(fitted(lme.fit$gam)[glm.spl$stn == 1], col = "black", z, lwd = 2)

#variogram
gamm.var <- Variogram(lme.fit.CAR$lme, form = ~ z | stn, robust = TRUE)
plot(gamm.var)


#asreml with station random effect

z_ind <- outer(glm.spl$z, seq(0, 250, 25), FUN = "==") %*% matrix(data = rep(1, 11), ncol = 1, nrow = 11)

glm.sub <- glm.spl[z_ind == 1, ]
asreml.fit <- asreml(fixed = l.obs ~ z, random =~ spl(z) + stn, data = glm.spl, splinepoints = list(z = seq(0, 250, 25)))
summary(asreml.fit)


#lme splines
glm.spl <- data.frame(obs, l.obs, rep(z, n.station), rep(c(1:n.station), 1, each = length(z)))
names(glm.spl) <- c("obs", "l.obs", "z", "stn")
glm.spl$stn <- as.factor(glm.spl$stn)
glm.spl$all <- rep(1, nrow(glm.spl))
glm.spl$Zt  <- smspline(~ z, data = glm.spl)
lme.fit <- lme(l.obs ~ z,  random = list(all=pdIdent(~Zt - 1), stn =~1), data = glm.spl)
summary(lme.fit)

plot(log(rho),z,ylim=c(max(z),min(z)),type='l', xlim = c(min(l.obs), max(l.obs)))
points(l.obs, glm.spl$z)
lines(fitted(lme.fit)[glm.spl$stn == 1],z,col='red',lwd=2)

#----------------------------- GRIDDED DATA SET -----------------------------------#

#generate locations for grid using 25 stations
x <- sort(rep(sort(rep(seq(1:sqrt(n.station)), length(rho))), sqrt(n.station)))
y <- rep(sort(rep(seq(1:sqrt(n.station)), length(rho))), sqrt(n.station))
grid.dat <- cbind(glm.spl, x, y)

#asreml without correlation by distance
asreml.fit <- asreml(fixed = l.obs ~ z, random =~ spl(z) + stn + x + y, data = grid.dat, splinepoints = list(z = seq(0, 250, 25)))
summary(asreml.fit)

summary(asreml.fit)$varcomp[2,2]^0.5
summary(asreml.fit)$varcomp[5,2]^0.5

#--------------------------- AR PROCESS DOWN DEPTHS -------------------------------#

n.station <- 20 #specify number of stations to run
mu <- 50
sd <- 40
noise.sd <- 0.2 #noise sd for ar process
stn.sd   <- 0.5 #sd for random station effect

mult <- 1e3
z <- seq(0, 250, 5) #explanatory variable (depth)
rho <- mult*dnorm(z, mu, sd)/(pnorm(max(z), mu, sd) - pnorm(min(z), mu, noise.sd))
z.ar <- arima.sim(n = length(rho), model = list(ar = 0.8), sd = noise.sd)
stn.re <- rnorm(n.station, mean = 0, sd = stn.sd) #station specific random effect

#calculate obs, the observed data, and log(obs)
l.obs <- rep(log(rho), n.station) + rep(z.ar, n.station) + rep(stn.re, 1, each = length(rho))
obs <- exp(l.obs)

#build data frame
glm.spl <- data.frame(obs, l.obs, rep(z, n.station), rep(c(1:n.station), 1, each = length(z)), as.factor(z))
names(glm.spl) <- c("obs", "l.obs", "z", "stn", "z.fact")
glm.spl$stn <- as.factor(glm.spl$stn)

#asreml
asreml.fit <- asreml(fixed = l.obs ~ z, random =~ spl(z) + stn, data = glm.spl, 
                     splinepoints = list(z = seq(0, 250, 25)), rcov =~ stn:ar1(z.fact))
summary(asreml.fit)
summary(asreml.fit)$varcomp[2,2]^0.5
summary(asreml.fit)$varcomp[3,2]^0.5
summary(asreml.fit)$varcomp[4,2]^0.5







