# simulates data on a regular 3D grid with autocorrelated errors
library(mgcv)

mu=50;sd=40
noise.sd=0.2

mult=1e3
a <- 2 #shape
s <- 30 #scale
z <- seq(0, 250, by = 5)
rho = mult*(1/(s^a*gamma(a))*z^(a-1)*exp(-(x/s)))
obs = rho + rnorm(n = length(rho), mean=0, sd = noise.sd)
plot(z, obs)






mu=50;sd=40
noise.sd=0.2

mult=1e3
z=seq(0,250,5)
rho=mult*dnorm(z,mu,sd)/(pnorm(max(z),mu,sd)-pnorm(min(z),mu,sd))
obs=rho+rnorm(n=length(rho),mean=0,sd=noise.sd)

lmObj=lm(obs~poly(z,3))
preds=predict(lmObj,se.fit=TRUE)
par(mfrow=c(1,2))
plot(rho,z,ylim=c(max(z),min(z)),type='l')
title('Polynomial')
points(obs,z)
lines(preds$fit,z,col='grey',lwd=2)
lines(preds$fit-2*preds$se.fit*2,z,col='grey',lty=2,lwd=2)
lines(preds$fit+2*preds$se.fit*2,z,col='grey',lty=2,lwd=2)


#gam
gam2 <- gam(log(obs + 2) ~ s(z), family=Gamma()) 
plot(log(rho + 1), z,ylim=c(max(z),min(z)),type='l')
title("gam - without noise")
points(log(obs + 1), z)
lines(gam2$fit, z, col = 'red', lwd = 2)

#asreml spline mixed model
asreml.01 <- asreml(fixed = obs ~ z, random =~ spl(z))
summary(asreml.01)
wald(asreml.01)

plot(rho, z, ylim = c(max(z), min(z)), type = 'l')
title("asreml")
points(obs, z)
lines(fitted(asreml.01), z, col = 'red', lwd = 2)

#lme with splines (same as asreml)
glm.spl <- data.frame(obs, z)
glm.spl$all <- rep(1, nrow(glm.spl))
glm.spl$Zt  <- smspline(~ z, data = glm.spl)
lme.fit <- lme(obs ~ z,  random=list(all=pdIdent(~Zt - 1)), data = glm.spl)
summary(lme.fit)

plot(rho,z,ylim=c(max(z),min(z)),type='l')
title("lme with splines")
points(obs,z)
lines(fitted(lme.fit),z,col='red',lwd=2)

#lme using three stations without noise
mu=50;sd=40
noise.sd=0

mult=1e3
z=seq(0,250,5)
rho=mult*dnorm(z,mu,sd)/(pnorm(max(z),mu,sd)-pnorm(min(z),mu,sd))
obs=rho+rnorm(n=length(rho),mean=0,sd=noise.sd)
obs2 <- rho+rnorm(n=length(rho),mean=0,sd=noise.sd) - 1
obs3 <- rho+rnorm(n=length(rho),mean=0,sd=noise.sd) + 1

glm.spl <- data.frame(c(obs, obs2, obs3), rep(z, 3), sort(rep(c(1:3), length(z))))
names(glm.spl) <- c("obs", "z", "stn")
glm.spl$stn <- as.factor(glm.spl$stn)
glm.spl$all <- rep(1, nrow(glm.spl))
glm.spl$Zt  <- smspline(~ z, data = glm.spl)
lme.fit <- lme(obs ~ z + stn,  random=list(all=pdIdent(~Zt - 1)), data = glm.spl)
summary(lme.fit)

plot(rho,z,ylim=c(max(z),min(z)),type='l', xlim = c(min(obs2), max(obs3)))
points(obs, z)
points(obs2, z, col = "red")
points(obs3, z, col = "blue")
title("spline lme - 3 stations without noise")
lines(fitted(lme.fit)[glm.spl$stn == 1],z,col='black',lwd=2)
lines(fitted(lme.fit)[glm.spl$stn == 2],z,col='red',lwd=2)
lines(fitted(lme.fit)[glm.spl$stn == 3],z,col='blue',lwd=2)




#lme using three stations with noise
mu=50;sd=40
noise.sd=0.2
noise.sd2=0.2
noise.sd3=0.2

mult=1e3
z=seq(0,250,5)
rho=mult*dnorm(z,mu,sd)/(pnorm(max(z),mu,sd)-pnorm(min(z),mu,sd))
obs=rho+rnorm(n=length(rho),mean=0,sd=noise.sd)
obs2 <- rho+rnorm(n=length(rho),mean=0,sd=noise.sd2) - 1
obs3 <- rho+rnorm(n=length(rho),mean=0,sd=noise.sd3) + 1

glm.spl <- data.frame(c(obs, obs2, obs3), rep(z, 3), sort(rep(c(1:3), length(z))))
names(glm.spl) <- c("obs", "z", "stn")
glm.spl$stn <- as.factor(glm.spl$stn)
glm.spl$all <- rep(1, nrow(glm.spl))
glm.spl$Zt  <- smspline(~ z, data = glm.spl)
lme.fit <- lme(obs ~ z + stn,  random=list(all=pdIdent(~Zt - 1)), data = glm.spl)
summary(lme.fit)

plot(rho,z,ylim=c(max(z),min(z)),type='l', xlim = c(min(obs2), max(obs3)))
points(obs, z)
points(obs2, z, col = "red")
points(obs3, z, col = "blue")
title("spline lme - 3 stations with noise")
lines(fitted(lme.fit)[glm.spl$stn == 1],z,col='black',lwd=2)
lines(fitted(lme.fit)[glm.spl$stn == 2],z,col='red',lwd=2)
lines(fitted(lme.fit)[glm.spl$stn == 3],z,col='blue',lwd=2)



#gam using three stations without noise
mu=50;sd=40
noise.sd=0

mult=1e3
z=seq(0,250,5)
rho=mult*dnorm(z,mu,sd)/(pnorm(max(z),mu,sd)-pnorm(min(z),mu,sd))
obs=rho+rnorm(n=length(rho),mean=0,sd=noise.sd)
obs2 <- rho+rnorm(n=length(rho),mean=0,sd=noise.sd) - 1
obs3 <- rho+rnorm(n=length(rho),mean=0,sd=noise.sd) + 1

glm.spl <- data.frame(c(obs, obs2, obs3), rep(z, 3), sort(rep(c(1:3), length(z))))
names(glm.spl) <- c("obs", "z", "stn")
glm.spl$stn <- as.factor(glm.spl$stn)
gam <- gam(obs + 1 ~ s(z) + stn, data = glm.spl)


plot(rho,z,ylim=c(max(z),min(z)),type='l', xlim = c(min(obs2), max(obs3)))
points(obs, z)
points(obs2, z, col = "red")
points(obs3, z, col = "blue")
title("gam - 3 stations without noise")
lines(gam$fit[glm.spl$stn == 1] - 1,z,col='black',lwd=2)
lines(gam$fit[glm.spl$stn == 2] - 1,z,col='red',lwd=2)
lines(gam$fit[glm.spl$stn == 3] - 1,z,col='blue',lwd=2)


#gam using three stations with noise
mu=50;sd=40
noise.sd=0.2

mult=1e3
z=seq(0,250,5)
rho=mult*dnorm(z,mu,sd)/(pnorm(max(z),mu,sd)-pnorm(min(z),mu,sd))
obs=rho+rnorm(n=length(rho),mean=0,sd=noise.sd)
obs2 <- rho+rnorm(n=length(rho),mean=0,sd=noise.sd) - 1
obs3 <- rho+rnorm(n=length(rho),mean=0,sd=noise.sd) + 1

glm.spl <- data.frame(c(obs, obs2, obs3), rep(z, 3), sort(rep(c(1:3), length(z))))
names(glm.spl) <- c("obs", "z", "stn")
glm.spl$stn <- as.factor(glm.spl$stn)
gam <- gam(obs + 1 ~ s(z) + stn, data = glm.spl)


plot(rho,z,ylim=c(max(z),min(z)),type='l', xlim = c(min(obs2), max(obs3)))
points(obs, z)
points(obs2, z, col = "red")
points(obs3, z, col = "blue")
title("gam - 3 stations with noise")
lines(gam$fit[glm.spl$stn == 1] - 1,z,col='black',lwd=2)
lines(gam$fit[glm.spl$stn == 2] - 1,z,col='red',lwd=2)
lines(gam$fit[glm.spl$stn == 3] - 1,z,col='blue',lwd=2)


# single station with correlated residuals
mu=50;sd=40
noise.sd=0.2

mult <- 1e3
z <- seq(0,250,5)
rho <- mult*dnorm(z,mu,sd)/(pnorm(max(z),mu,sd)-pnorm(min(z),mu,sd))
err <- arima.sim(length(rho), model = list(ar = 0.8), mean = 0, sd = noise.sd)
obs <- rho + err


glm.spl <- data.frame(obs, z)
glm.spl$all <- rep(1, nrow(glm.spl))
glm.spl$Zt  <- smspline(~ z, data = glm.spl)
lme.fit <- lme(obs ~ z,  random=list(all=pdIdent(~Zt - 1)), data = glm.spl, correlation = corAR1(0.8))
summary(lme.fit)

plot(rho,z,ylim=c(max(z),min(z)),type='l')
title("lme with splines")
points(obs,z)
lines(fitted(lme.fit),z,col='red',lwd=2)

# add station random effect
n.station <- 20
mu=50;sd=40
noise.sd=0.2
stn.sd <- rnorm(n.station, 0.5, 0.2)


mult <- 1e3
z <- seq(0,250,5)
rho <- mult*dnorm(z,mu,sd)/(pnorm(max(z),mu,sd)-pnorm(min(z),mu,sd))
noise <- rnorm(n = length(rho), mean = 0, sd = noise.sd)
stn.re <- rnorm(length(stn.sd), mean = 0, sd = stn.sd)

obs <- rep(rho, n.station) + rep(noise, n.station) + rep(stn.re, 1, each = length(rho))
l.obs <- log(obs + 2)

glm.spl <- data.frame(obs, l.obs, rep(z, n.station), sort(rep(c(1:n.station), length(z))))
names(glm.spl) <- c("obs", "l.obs", "z", "stn")
glm.spl$stn <- as.factor(glm.spl$stn)
lme.fit <- gamm(l.obs ~ s(z),  random=list(stn =~1), data = glm.spl)

summary(lme.fit$lme)


plot(log(rho + 2),z,ylim=c(max(z),min(z)),type='l', xlim = c(min(na.omit(l.obs)), max(na.omit(l.obs))), xlab = "l.fluoro")
points(l.obs, rep(z, n.station), col = glm.spl$stn)
title("gamm with station random effect")
lines(fitted(lme.fit$gam)[glm.spl$stn == 1], col = "black", z,lwd=2)

#asreml with station random effect
asreml.fit <- asreml(fixed = l.obs ~ z, random =~ spl(z) + stn, data = glm.spl)
summary(asreml.fit)



#generate locations for grid using 100 stations
x <- sort(rep(sort(rep(seq(1:10), length(rho))), 10))
y <- rep(sort(rep(seq(1:10), length(rho))), 10)
grid.dat <- cbind(glm.spl, x, y)

#model with no autocorrelation by distance
gamm.fit <- gamm(l.obs ~ s(z) + s(x) + s(y),  random=list(stn =~1), data = grid.dat)
summary(gamm.fit$lme)

#error by distance
err <- arima.sim(length(rho), model = list(ar = 0.8), mean = 0, sd = 0.2)


plot(Variogram(lme.fit$lme))




