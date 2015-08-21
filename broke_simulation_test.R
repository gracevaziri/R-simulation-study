#fits asreml model developed in simulation to BROKE-West data
#uses the distance between stations rather than latitude and longitude
#used in HDR conference talk in December 2014
#author: Lisa-Marie Harrison
#date: 18/09/2014

setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/Data")
dat.cut <- read.csv(file = "procCTD.csv", header= T)
library(asreml)
library(nlme)
library(lattice)
library(mgcv)

names(dat.cut) <- c("survey", "stn", "lat", "long", "start.time", "end.time", "depth", "transmittance", "cond", "temp", "sal", "par", "oxygen", "fluoro", "x2", "ice", "wm")

#remove null values
dat.cut$sal[dat.cut$sal == -9] <- NA
dat.cut$temp[dat.cut$temp == -9] <- NA
dat.cut$par[dat.cut$par == -9] <- NA
dat.cut$fluoro[dat.cut$fluoro == -9] <- NA

#compute log transformed fluoro values
dat.cut$l.fluoro <- log(dat.cut$fluoro)
dat.cut$l.fluoro[is.nan(dat.cut$l.fluoro)] <- NA

#get latitude and longitude for each station
n.station <- length(unique(dat.cut$stn))
lat  <- dat.cut$lat[duplicated(dat.cut$stn) == FALSE]
long <- dat.cut$long[duplicated(dat.cut$stn) == FALSE]

#plot location of BROKE-West station with station number overlayed
plot(long, lat, col = "white", xlab = "longitude", ylab = "latitude")
text(long, lat, c(2:118))
title("location of BROKE-West CTD stations")


#find maximum fluorescence depth
max.depth <- 0
for (i in 1:length(unique(dat.cut$stn))) {
  max.depth[i] <- which.max(dat.cut$l.fluoro[dat.cut$stn == unique(dat.cut$stn)[i]])
}
dat.cut$max.depth <- rep(unique(dat.cut$profile.depth)[max.depth], each = 125)


#function to convert degrees to radians
deg2rad <- function(deg) {
  return(deg*pi/180)
}

#Calculates the distance between two points with radian latitude/longitude using Haversine formula (hf)
gcd.hf <- function(lat1, long1, lat2, long2) {
  R <- 6371 # Earth mean radius [km]
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  d = R * c
  return(d) # Distance in km
}

#distance of each station from station 1 in x and y directions
x <- 0
y <- 0
rad_long <- deg2rad(long)
rad_lat  <- deg2rad(lat)
top_lat <- deg2rad(max(lat))
top_long <- deg2rad(max(long))
for (i in 1:n.station) {
  x[i] <- gcd.hf(rad_lat[i], top_long, top_lat, top_long)/100    
  y[i] <- gcd.hf(top_lat, rad_long[i], top_lat, top_long)/100
}


#data frame
glm.spl <- data.frame(dat.cut$l.fluoro, dat.cut$depth, as.factor(dat.cut$stn), rep(x, 1, each = length(unique(dat.cut$depth))), rep(y, 1, each = length(unique(dat.cut$depth))), dat.cut$temp, dat.cut$par, dat.cut$sal, dat.cut$oxygen, dat.cut$ice, as.factor(dat.cut$wm))
names(glm.spl) <- c("l.obs", "z", "stn", "x", "y", "temp", "par", "sal", "oxy", "ice", "wm")
glm.spl$z.fact <- as.factor(as.integer(glm.spl$z))
glm.spl$x.fact <- as.factor(glm.spl$x)
glm.spl$y.fact <- as.factor(glm.spl$y)
glm.spl <- glm.spl[order(glm.spl$z, glm.spl$x, glm.spl$y), ] #sort by order of rcov structure
glm.spl$l.obs[glm.spl$l.obs == -Inf] <- NA

#centre and scale covariates to mean = 0 and sd = 1
#this is required if using na.method = "include" since this sets the missing values to 0
glm.spl$temp <- scale(glm.spl$temp)
glm.spl$par  <- scale(glm.spl$par)
glm.spl$sal  <- scale(glm.spl$sal)
glm.spl$oxy  <- scale(glm.spl$oxy)
glm.spl$ice  <- scale(glm.spl$ice)
glm.spl$oxy  <- scale(glm.spl$oxy)

#------------------------------- FIT ASREML MODELS -----------------------------------#

#fit asreml model
asreml.fit <- asreml(fixed = l.obs ~ z + par + temp:wm + oxy + sal, random =~ spl(z, 10) + spl(par, 10) + 
                        spl(temp, 10):wm + spl(oxy, 10) + spl(sal, 10) + stn, 
                     data = glm.spl, rcov=~ ar1(z.fact):agau(x.fact, y.fact),
                     na.method.X = "include", workspace = 50000000)
asreml.fit <- update(asreml.fit)
summary(asreml.fit)


#plot fitted against observed for all stations
lat.plot <- xyplot(glm.spl$l.obs + fitted(asreml.fit) ~ glm.spl$z | glm.spl$stn, 
                   outer = FALSE, type = "l", xlab = "depth (m)", ylab = "l.fluoro")
update(lat.plot, par.settings = simpleTheme(lwd = c(2, 1), col = c("dodgerblue", "red")))


#plot fitted against observed for 3 stations
s <- sample(c(2:118), 3, replace = F)
lat.plot <- xyplot(glm.spl$l.obs[glm.spl$stn %in% s] + fitted(asreml.fit)[glm.spl$stn %in% s] ~ glm.spl$z[glm.spl$stn %in% s] | glm.spl$stn[glm.spl$stn %in% s], 
                   outer = FALSE, type = "l", xlab = list("depth (m)", cex = 2), ylab = list("l.fluoro", cex = 2), scales = list(cex = 2), cex.axis = 2)
update(lat.plot, par.settings = simpleTheme(lwd = c(2, 1), col = c("dodgerblue", "red")), lwd = 2,  cex.axis = 2, cex.lab = 2)

#convert back to fluoro
lat.plot <- xyplot(exp(glm.spl$l.obs) + exp(fitted(asreml.fit)) ~ glm.spl$z | glm.spl$stn, 
                   outer = FALSE, type = "l", xlab = list("depth (m)", cex = 2), ylab = list("l.fluoro", cex = 2), scales = list(cex = 2), cex.axis = 2)
update(lat.plot, par.settings = simpleTheme(lwd = c(2, 1), col = c("dodgerblue", "red")), lwd = 3,  cex.axis = 2, cex.lab = 2)

#--------------------------- CHECK AUTOCORRELATION ----------------------------#

#variogram of residuals for model with and without correlation strucutre
d <- 62
gamma <- asreml.variogram(glm.spl$z[glm.spl$stn == d], z = resid(asreml.fit)[glm.spl$stn == d])$gamma
dist  <- asreml.variogram(glm.spl$z[glm.spl$stn == d], z = resid(asreml.fit)[glm.spl$stn == d])$x
plot(dist, gamma)

gamma <- asreml.variogram(glm.spl$z[glm.spl$stn == d], z = resid(asreml.null)[glm.spl$stn == d])$gamma
dist  <- asreml.variogram(glm.spl$z[glm.spl$stn == d], z = resid(asreml.null)[glm.spl$stn == d])$x
plot(dist, gamma)

#likelihood ratio test to check whether adding correlation structure works
1 - pchisq(2 * (asreml.fit$loglik - asreml.null$loglik), 1) 



#---------------------------------- FIT GAMM -----------------------------------------#

#fit gamm to compare
gamm.fit <- gamm(l.obs ~ s(z) + s(temp, by = wm) + s(par) + s(ice), random = list(stn =~ 1, x =~1, y =~1), 
                 data = glm.spl, correlation = corAR1(0.9, 1 ~ z | x | y))
summary(gamm.fit$gam)


#plot fitted against observed by station
lat.plot <- xyplot(glm.spl$l.obs + fitted(gamm.fit$gam) ~ glm.spl$z | glm.spl$stn, outer = FALSE, type = "l")
update(lat.plot, par.settings = simpleTheme(lwd = c(2, 1), col = c("dodgerblue", "red")))



#bubble plot of mean residuals by station
res <- aggregate(residuals(asreml.fit), by = list(glm.spl$stn), FUN = mean, na.rm = TRUE)$x
bubble_dat <- as.data.frame(cbind(long, lat, res))
colnames(bubble_dat) <- c("long", "lat", "res")

p1 <- ggplot(bubble_dat[bubble_dat$res < 0, ], guide = FALSE) + 
  geom_point(aes(x=long, y=lat, size=abs(res)), colour="red", fill = "red", shape = 21)+ scale_size_area(max_size = 15) +
    scale_x_continuous(name="Longitude") +
  scale_y_continuous(name="Latitude") +
  theme_bw() + 
  theme(legend.title=element_blank(), text = element_text(size=20)) 
p1 + geom_point(data = bubble_dat[bubble_dat$res >= 0, ], aes(x=long, y=lat, size=abs(res)), colour="blue", fill = "blue", shape = 21, add = TRUE)



#bubble plot of mean fluorescence by station
res <- aggregate(exp(glm.spl$l.obs), by = list(glm.spl$stn), FUN = mean, na.rm = TRUE)$x
bubble_dat <- as.data.frame(cbind(long, lat, res))
colnames(bubble_dat) <- c("long", "lat", "res")

ggplot(bubble_dat, guide = FALSE) + 
  geom_point(aes(x=long, y=lat, size=res), shape = 21)+ scale_size_area(max_size = 15) +
  scale_x_continuous(name="Longitude") +
  scale_y_continuous(name="Latitude") +
  theme_bw() + 
  theme(legend.title=element_blank(), text = element_text(size=20)) 

#-------------------------- AVERAGE PREDICTIONS -------------------------------#

par(mfrow = c(1, 2))

#temperature
pred <- predict(asreml.fit, classify = "temp")
pval <- pred$predictions$pvals["predicted.value"]$predicted.value
temp <- pred$predictions$pvals["temp"]$temp
se <- pred$predictions$pvals["standard.error"]$standard.error

logci <- pval + se%*%t(qnorm(c(0.025,0.5,0.975)))
ci <- exp(logci)
dimnames(ci)[[2]]<-c("lower95", "est", "upper95")

plot(temp, ci[, 2], xlab = "temperature (degrees celcius)", ylab = "", 
     type = "l", ylim = c(min(ci[, 1]), max(ci[, 3])), cex.lab = 2, cex.axis = 2)
points(temp, ci[, 1], type = "l", lty = 2)
points(temp, ci[, 3], type = "l", lty = 2)


#par
pred <- predict(asreml.fit, classify = "par")
pval <- pred$predictions$pvals["predicted.value"]$predicted.value
par <- pred$predictions$pvals["par"]$par
se <- pred$predictions$pvals["standard.error"]$standard.error

logci <- pval + se%*%t(qnorm(c(0.025,0.5,0.975)))
ci <- exp(logci)
dimnames(ci)[[2]]<-c("lower95", "est", "upper95")

plot(par, ci[, 2], xlab = expression("par" ~ (??E ~ m^{???2} ~ s^{???1})), ylab = "", 
     type = "l", ylim = c(min(ci[, 1]), max(ci[, 3])), cex.lab = 2, cex.axis = 2)
points(par, ci[, 1], type = "l", lty = 2)
points(par, ci[, 3], type = "l", lty = 2)


#oxygen
pred <- predict(asreml.fit, classify = "oxy")
pval <- pred$predictions$pvals["predicted.value"]$predicted.value
oxy <- pred$predictions$pvals["oxy"]$oxy
se <- pred$predictions$pvals["standard.error"]$standard.error

logci <- pval + se%*%t(qnorm(c(0.025,0.5,0.975)))
ci <- exp(logci)
dimnames(ci)[[2]]<-c("lower95", "est", "upper95")

plot(oxy, ci[, 2], xlab = expression("dissolved oxygen" ~ (??mol ~ L^{???1})), ylab = "", 
     type = "l", ylim = c(min(ci[, 1]), max(ci[, 3])), cex.lab = 2, cex.axis = 2)
points(oxy, ci[, 1], type = "l", lty = 2)
points(oxy, ci[, 3], type = "l", lty = 2)



#depth
pred <- predict(asreml.fit, classify = "z")
pval <- pred$predictions$pvals["predicted.value"]$predicted.value
z <- pred$predictions$pvals["z"]$z
se <- pred$predictions$pvals["standard.error"]$standard.error

logci <- pval + se%*%t(qnorm(c(0.025,0.5,0.975)))
ci <- exp(logci)
dimnames(ci)[[2]]<-c("lower95", "est", "upper95")

plot(z, ci[, 2], xlab = "depth (m)", ylab = "", 
     type = "l", ylim = c(min(ci[, 1]), max(ci[, 3])), cex.lab = 2, cex.axis = 2)
points(z, ci[, 1], type = "l", lty = 2)
points(z, ci[, 3], type = "l", lty = 2)


#salinity
pred <- predict(asreml.fit, classify = "sal")
pval <- pred$predictions$pvals["predicted.value"]$predicted.value
z <- pred$predictions$pvals["sal"]$sal
se <- pred$predictions$pvals["standard.error"]$standard.error

logci <- pval + se%*%t(qnorm(c(0.025,0.5,0.975)))
ci <- exp(logci)
dimnames(ci)[[2]]<-c("lower95", "est", "upper95")

plot(z, ci[, 2], xlab = "salinity (psu)", ylab = "", 
     type = "l", ylim = c(min(ci[, 1]), max(ci[, 3])), cex.lab = 2, cex.axis = 2)
points(z, ci[, 1], type = "l", lty = 2)
points(z, ci[, 3], type = "l", lty = 2)

#-------------------------- CTD VERTICAL PROFILE ------------------------------#

par(mfrow = c(1, 3))

plot(dat.cut$temp[dat.cut$stn == 2], -dat.cut$depth[dat.cut$stn == 2], type = "l",
     xlab = "temperature", ylab = "depth (m)", lwd = 2, cex.lab = 1.5)
title("Temperature")
plot(dat.cut$oxygen[dat.cut$stn == 2], -dat.cut$depth[dat.cut$stn == 2], type = "l",
     xlab = "dissolved oxygen", ylab = "", lwd = 2, cex.lab = 1.5)
title("Dissolved Oxygen")


#------------------------ CLIMATE CHANGE PREDICTIONS --------------------------#


#normal prediction for station 2
ci <- matrix(NA, nrow = 125, ncol = 3)
for (i in 1:length(glm.spl$z[glm.spl$stn == 2])){
  pred <- predict(asreml.fit, classify = "temp:z:stn", levels = list("stn" = 1, "temp" = glm.spl$temp[glm.spl$stn == 2][i], "z" = i))
  pval <- pred$predictions$pvals["predicted.value"]$predicted.value
  se <- pred$predictions$pvals["standard.error"]$standard.error
  
  logci <- pval + se%*%t(qnorm(c(0.025,0.5,0.975)))
  ci[i, ] <- exp(logci)
}
dimnames(ci)[[2]]<-c("lower95", "est", "upper95")


#temperature increase of 2 degrees at station 2
ci2 <- matrix(NA, nrow = 125, ncol = 3)
for (i in 1:length(glm.spl$z[glm.spl$stn == 2])){
  pred <- predict(asreml.fit, classify = "temp:z:stn", levels = list("stn" = 1, "temp" = glm.spl$temp[glm.spl$stn == 2][i] + 2, "z" = i))
  pval <- pred$predictions$pvals["predicted.value"]$predicted.value
  se <- pred$predictions$pvals["standard.error"]$standard.error
  
  logci <- pval + se%*%t(qnorm(c(0.025,0.5,0.975)))
  ci2[i, ] <- exp(logci)
}
dimnames(ci2)[[2]]<-c("lower95", "est", "upper95")

#plot before and after
plot(glm.spl$z[glm.spl$stn == 2], ci[, 2], xlab = "depth (m)", ylab = "predicted fluoro", 
     type = "l", ylim = c(min(na.omit(ci[, 1])), 1), lwd = 2)
points(glm.spl$z[glm.spl$stn == 2], ci2[, 2], type = "l", col = "red", lwd = 2)
legend(150, 1, c("Current station 2 prediction", "2 degrees celcius temperature increase"), col = c("black", "red"), lwd = 2, bty = "n")



#---------------- FIT ASREML MODEL WITHOUT CORRELATION STRUCTURE --------------#

#fit asreml model
asreml.full <- asreml(fixed = l.obs ~ z + par + temp:wm + oxy + sal, random =~ spl(z, 10) + spl(par, 10) + 
                       spl(temp, 10):wm + spl(oxy, 10) + spl(sal, 10) + stn, 
                     data = glm.spl, rcov=~ ar1(z.fact):agau(x.fact, y.fact),
                     na.method.X = "include", workspace = 50000000, aom = T)
asreml.full <- update(asreml.full)
summary(asreml.full)


#fit null asreml model
asreml.null <- asreml(fixed = l.obs ~ z  + temp:wm + oxy + sal, random =~ spl(z, 10) + 
                       spl(temp, 10):wm + spl(oxy, 10) + spl(sal, 10) + stn,
                     data = glm.spl, na.method.X = "include", workspace = 50000000, aom = T)
summary(asreml.null)

#intercept model
asreml.null <- asreml(fixed = l.obs ~ z,
                      data = glm.spl, na.method.X = "include", workspace = 50000000, aom = T)
summary(asreml.null)

#plot fitted against observed for all stations
lat.plot <- xyplot(glm.spl$l.obs + fitted(asreml.full) ~ glm.spl$z | glm.spl$stn, 
                   outer = FALSE, type = "l", xlab = "depth (m)", ylab = "l.fluoro")
update(lat.plot, par.settings = simpleTheme(lwd = c(2, 1), col = c("dodgerblue", "red")))


lat.plot <- xyplot(glm.spl$l.obs + fitted(asreml.null) ~ glm.spl$z | glm.spl$stn, 
                   outer = FALSE, type = "l", xlab = "depth (m)", ylab = "l.fluoro")
update(lat.plot, par.settings = simpleTheme(lwd = c(2, 1), col = c("dodgerblue", "red")))



plot(asreml.null$resid, asreml.full$resid, xlim = c(-4, 4), ylim = c(-4, 4))
plot(asreml.null$fitted, asreml.full$fitted)


hist(glm.spl$l.obs - residuals(asreml.full))
hist(glm.spl$l.obs - residuals(asreml.null))


#calculate R squared

ss_tot <- sum(na.omit(glm.spl$l.obs - rep(aggregate(glm.spl$l.obs, list(glm.spl$stn), FUN = mean, na.rm = T)$x, 125))^2)


ss_res_full <- sum(na.omit(asreml.full$resid)^2)
r_squared_full <- 1 - (ss_res_full/ss_tot)

ss_res_null <- sum(na.omit(asreml.null$resid)^2)
r_squared_null <- 1 - (ss_res_null/ss_tot)

dat <- cbind(glm.spl$l.obs, fitted(asreml.full))
dat <- na.omit(dat)
KL.divergence(dat[, 1], dat[, 2], k = 1)


