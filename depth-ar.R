#--------------------------- AR PROCESS DOWN DEPTHS -------------------------------#
library(asreml)
library(mgcv)
library(lmeSplines)

n.station <- 20 #specify number of stations to run
mu <- 50
sd <- 40
noise.sd <- 0.2 #noise sd for ar process
stn.sd   <- 0.5 #sd for random station effect
phi <- 0.8

mult <- 1e3
z <- seq(0, 250, 5) #explanatory variable (depth)
rho <- mult*dnorm(z, mu, sd)/(pnorm(max(z), mu, sd) - pnorm(min(z), mu, noise.sd))
z.ar <- arima.sim(n = length(rho), model = list(ar = phi), sd = noise.sd) #ar process down z (same for each station)
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

v1 <- variogram.asreml(asreml.fit)
plot(v1$gamma[v1$stn == 1])



#gamm 
gam.fit <- gamm(l.obs ~ s(z),  random=list(stn =~1), data = glm.spl, correlation = corAR1(phi, form =~ z|stn))
summary(gam.fit$lme)
summary(gam.fit$gam)


#lme splines
glm.spl$all <- rep(1, nrow(glm.spl))
glm.spl$Zt  <- smspline(~ z, data = glm.spl)
lme.fit <- lme(l.obs ~ z,  random = list(all=pdIdent(~Zt - 1), stn =~1), data = glm.spl, correlation = corAR1(phi, form =~ z|stn))
summary(lme.fit)








