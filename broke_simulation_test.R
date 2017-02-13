#fits asreml model developed in simulation to BROKE-West data
#uses the distance between stations rather than latitude and longitude
#used in HDR conference talk in December 2014
#author: Lisa-Marie Harrison
#date: 18/09/2014

if (Sys.info()[4] == "SCI-6246") {
  setwd(dir = "C:/Users/43439535/Documents/Lisa/phd/Mixed models")
} else {
  setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models")
}

dat.cut <- read.csv(file = "Data/procCTD.csv", header= T)
library(asreml)
library(nlme)
library(lattice)
library(mgcv)
library(ggplot2)

names(dat.cut) <- c("survey", "stn", "lat", "long", "start.time", "end.time", "depth", "transmittance", "cond", "temp", "sal", "par", "oxygen", "fluoro", "x2", "ice", "wm")

#source functions
function_list <- c("gcdHF.R", 
                   "deg2rad.R", 
                   "depthFluoroMax.R", 
                   "distFromStn1.R", 
                   "calc_asreml_conditional_marginal_Rsquared.R",
                   "asremlAIC.R")

for (f in function_list) {
  
    source(paste("R code/R-functions-southern-ocean/", f, sep = ""))
  
}



#remove null values
dat.cut$sal[dat.cut$sal == -9] <- NA
dat.cut$temp[dat.cut$temp == -9] <- NA
dat.cut$par[dat.cut$par == -9] <- NA
dat.cut$fluoro[dat.cut$fluoro == -9] <- NA

#compute log transformed fluoro values
dat.cut$l.fluoro <- log(dat.cut$fluoro)
dat.cut$l.fluoro[is.nan(dat.cut$l.fluoro)] <- NA

#get latitude and longitude for each station
lat  <- dat.cut$lat[duplicated(dat.cut$stn) == FALSE]
long <- dat.cut$long[duplicated(dat.cut$stn) == FALSE]

#plot location of BROKE-West station with station number overlayed
plot(long, lat, col = "white", xlab = "longitude", ylab = "latitude", cex.lab = 2)
text(long, lat, sort(unique(dat.cut$stn)))
title("location of BROKE-West CTD stations", cex = 2)

#find depth of fluoro maximum
max.depth <- depthFluoroMax(dat.cut)

#distance of each station from station 1 in x and y directions
dist_stn_1 <- distFromStn1(lat, long)
x <- dist_stn_1$x
y <- dist_stn_1$y

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
glm.spl <- cbind(glm.spl[, c(1:5, 11:14)], apply(glm.spl[, c(6:10)], 2, scale))


#------------------------------- FIT ASREML MODELS -----------------------------------#

#fit asreml model
asreml.fit <- asreml(fixed = l.obs ~ z + temp:wm + oxy + sal, random =~ spl(z, 10) + 
                        spl(temp, 10):wm + spl(oxy, 10) + spl(sal, 10) + stn, 
                     data = glm.spl, rcov=~ ar1(z.fact):agau(x.fact, y.fact),
                     na.method.X = "include", workspace = 50000000)
asreml.fit <- update(asreml.fit)
summary(asreml.fit)

#assess collinearity using Cholesky decomposition of correlation matrix
#Values near 0 in the diagonal indicates variables are collinear enough to cause problems
chol(with(na.omit(glm.spl), cor(data.frame(z, temp, sal, oxy, ice))))


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

#---------------------------------- FIT GAMM -----------------------------------------#

#fit gamm to compare
gamm.fit <- gamm(l.obs ~ s(z) + s(temp, by = wm) + s(sal) + s(oxy), random = list(stn =~ 1, x =~1, y =~1), 
                 data = glm.spl, correlation = corAR1(0.9, 1 ~ z | x | y))
summary(gamm.fit$gam)

#----------------------------- PLOTS FOR PAPER -------------------------------#

#plot fitted against observed by station
lat.plot <- xyplot(glm.spl$l.obs + fitted(gamm.fit$gam) ~ glm.spl$z | glm.spl$stn, outer = FALSE, type = "l")
update(lat.plot, par.settings = simpleTheme(lwd = c(2, 1), col = c("dodgerblue", "red")))

#polygon of coastline
library(maptools)
library(mapdata)
library(rgdal)
library(raster)
library(rgeos)

shape <- readOGR("C:/Users/Lisa/Documents/phd/southern ocean/BROKE-West/polygon", "moa_2004_coastline_v1.1")

shape_ll <- spTransform(shape, CRS("+proj=longlat +datum=WGS84"))

shape_crop <- as(extent(rbind(range(bubble_dat$long), range(bubble_dat$lat))), "SpatialPolygons")

proj4string(shape_crop) <- CRS(proj4string(shape_ll))

out <- gIntersection(shape_ll, shape_crop, byid=TRUE)


#bubble plot of mean residuals by station
pdf("C:/Users/Lisa/Dropbox/uni/MEE submitted/images/Figure_2a.pdf")
res <- aggregate(residuals(asreml.fit), by = list(glm.spl$stn), FUN = mean, na.rm = TRUE)$x
bubble_dat <- as.data.frame(cbind(long, lat, res))
colnames(bubble_dat) <- c("long", "lat", "res")

p1 <- ggplot(bubble_dat[bubble_dat$res < 0, ], guide = FALSE) + 
  geom_point(aes(x=long, y=lat, size=abs(res)), colour="red", fill = "red", shape = 21)+ scale_size_area(max_size = 8) +
  geom_polygon(data=fortify(shape_ll[1, 1]), aes(x=long, y=lat), color = "black", fill = NA) +
  coord_map(xlim = range(bubble_dat$long) + c(-10, 10), ylim = range(bubble_dat$lat) + c(-2, 1)) +
  scale_x_continuous(name="Longitude") +
  scale_y_continuous(name="Latitude") +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  theme(legend.title=element_blank(), text = element_text(size=20), 
        legend.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text = element_text(color = "black"),
        panel.background = element_blank())
p1 + geom_point(data = bubble_dat[bubble_dat$res >= 0, ], aes(x=long, y=lat, size=abs(res)), colour="black", fill = "black", shape = 21)

dev.off()


#bubble plot of mean fluorescence by station
pdf("C:/Users/43439535/Dropbox/uni/MEE submitted/images/Figure_2b.pdf")
res <- aggregate(exp(glm.spl$l.obs), by = list(glm.spl$stn), FUN = mean, na.rm = TRUE)$x
bubble_dat <- as.data.frame(cbind(long, lat, res))
colnames(bubble_dat) <- c("long", "lat", "res")

ggplot(bubble_dat, guide = FALSE) + 
  geom_point(aes(x=long, y=lat, size=res), shape = 21, fill = "black")+ scale_size_area(max_size = 8) +
  geom_polygon(data=fortify(shape_ll[1, 1]), aes(x=long, y=lat), color = "black", fill = NA) +
  coord_map(xlim = range(bubble_dat$long) + c(-10, 10), ylim = range(bubble_dat$lat) + c(-2, 1)) +
  scale_x_continuous(name="Longitude") +
  scale_y_continuous(name="Latitude") +
  theme_bw() + 
  theme(legend.title=element_blank(), text = element_text(size=20), axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        panel.background = element_blank()) 
dev.off()

#-------------------------- AVERAGE PREDICTIONS -------------------------------#

pdf("C:/Users/Lisa/Dropbox/uni/MEE submitted/images/Figure_3.pdf")

par(mar=c(4.1,4.1,3.1,2.1),mfrow=c(2, 2))
par(oma = c(3, 6, 0, 0))

#temperature
pred <- predict(asreml.fit, classify = "temp", ignore = "z")
pval <- pred$predictions$pvals["predicted.value"]$predicted.value
temp <- pred$predictions$pvals["temp"]$temp
se <- pred$predictions$pvals["standard.error"]$standard.error

logci <- pval + se%*%t(qnorm(c(0.025,0.5,0.975)))
ci <- exp(logci)
dimnames(ci)[[2]]<-c("lower95", "est", "upper95")

plot(temp, ci[, 2], xlab = expression(temperature~(~degree~C)), ylab = "", bty = "n",
     type = "l", ylim = c(min(ci[, 1]), max(ci[, 3])), cex.lab = 2, cex.axis = 2)
polygon(c(temp, rev(temp)), c(ci[, 1], rev(ci[, 3])), col = rgb(0.84, 0.84, 0.84, 1), border = NA)
points(temp, ci[, 2], lwd = 2, type = "l")
rug(unique(na.omit(glm.spl$temp)), ticksize = 0.03, side = 1, lwd = 0.5)

#oxygen
pred <- predict(asreml.fit, classify = "oxy", ignore = "z")
pval <- pred$predictions$pvals["predicted.value"]$predicted.value
oxy <- pred$predictions$pvals["oxy"]$oxy
se <- pred$predictions$pvals["standard.error"]$standard.error

logci <- pval + se%*%t(qnorm(c(0.025,0.5,0.975)))
ci <- exp(logci)
dimnames(ci)[[2]]<-c("lower95", "est", "upper95")

plot(oxy, ci[, 2], xlab = expression("dissolved oxygen" ~ (mu~mol ~ L^{-1})), ylab = "", 
     type = "l", ylim = c(min(ci[, 1]), max(ci[, 3])), cex.lab = 2, cex.axis = 2, bty = "n")
polygon(c(oxy, rev(oxy)), c(ci[, 1], rev(ci[, 3])), col = rgb(0.84, 0.84, 0.84, 1), border = NA)
points(oxy, ci[, 2], lwd = 2, type = "l")
rug(unique(na.omit(glm.spl$oxy)), ticksize = 0.03, side = 1, lwd = 0.5)


#depth
pred <- predict(asreml.fit, classify = "z")
pval <- pred$predictions$pvals["predicted.value"]$predicted.value
z <- pred$predictions$pvals["z"]$z
se <- pred$predictions$pvals["standard.error"]$standard.error

logci <- pval + se%*%t(qnorm(c(0.025,0.5,0.975)))
ci <- exp(logci)
dimnames(ci)[[2]]<-c("lower95", "est", "upper95")

plot(z, ci[, 2], xlab = "depth (m)", ylab = "", bty = "n",
     type = "l", ylim = c(min(ci[, 1]), max(ci[, 3])), cex.lab = 2, cex.axis = 2)
polygon(c(z, rev(z)), c(ci[, 1], rev(ci[, 3])), col = rgb(0.84, 0.84, 0.84, 1), border = NA)
points(z, ci[, 2], lwd = 2, type = "l")
rug(unique(na.omit(glm.spl$z)), ticksize = 0.03, side = 1, lwd = 0.5)

#salinity
pred <- predict(asreml.fit, classify = "sal", ignore = "z")
pval <- pred$predictions$pvals["predicted.value"]$predicted.value
z <- pred$predictions$pvals["sal"]$sal
se <- pred$predictions$pvals["standard.error"]$standard.error

logci <- pval + se%*%t(qnorm(c(0.025,0.5,0.975)))
ci <- exp(logci)
dimnames(ci)[[2]]<-c("lower95", "est", "upper95")

plot(z, ci[, 2], xlab = "salinity (psu)", ylab = "", lwd = 2, bty = "n",
     type = "l", ylim = c(min(ci[, 1]), max(ci[, 3])), cex.lab = 2, cex.axis = 2)
polygon(c(z, rev(z)), c(ci[, 1], rev(ci[, 3])), col = rgb(0.84, 0.84, 0.84, 1), border = NA)
points(z, ci[, 2], lwd = 2, type = "l")
rug(unique(na.omit(glm.spl$sal)), ticksize = 0.03, side = 1, lwd = 0.5)

mtext(expression(hat(Fl) ~ (mu~g ~ L^{-1})), side = 2, outer = TRUE, line = 2, cex = 2)
dev.off()

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
asreml.full <- asreml(fixed = l.obs ~ z + temp:wm + oxy + sal, random =~ spl(z, 10) +  
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
asreml.int <- asreml(fixed = l.obs ~ 1,
                      data = glm.spl, na.method.X = "include", workspace = 50000000, aom = T)
summary(asreml.int)

#plot fitted against observed for all stations
lat.plot <- xyplot(glm.spl$l.obs + fitted(asreml.full) ~ glm.spl$z | glm.spl$stn, 
                   outer = FALSE, type = "l", xlab = "depth (m)", ylab = "l.fluoro")
update(lat.plot, par.settings = simpleTheme(lwd = c(2, 1), col = c("dodgerblue", "red")))


lat.plot <- xyplot(glm.spl$l.obs + fitted(asreml.null) ~ glm.spl$z | glm.spl$stn, 
                   outer = FALSE, type = "l", xlab = "depth (m)", ylab = "l.fluoro")
update(lat.plot, par.settings = simpleTheme(lwd = c(2, 1), col = c("dodgerblue", "red")))

#AIC values
source("R code/R-functions-southern-ocean/asremlAIC.R")
asremlAIC(asreml.full)
asremlAIC(asreml.null)
asremlAIC(asreml.int)

#likelihood ratio tests
1 - pchisq(2 * (asreml.null$loglik - asreml.int$loglik), 1) 
1 - pchisq(2 * (asreml.full$loglik - asreml.null$loglik), 1) 

#calculate R squared
calcRsquared(asreml.int, varStruct = F)
calcRsquared(asreml.null, rand = "stn", varStruct = F)
calcRsquared(asreml.full, rand = "stn", varStruct = T)

#RMSE at each station to show range in goodness of fit
rmse_stn <- 0
for (s in levels(glm.spl$stn)) {
  rmse_stn[as.numeric(s)] <- sqrt(mean(na.omit(glm.spl$l.obs[glm.spl$stn == s] - asreml.fit$fitted.values[glm.spl$stn == s])^2))
}
rmse_stn <- rmse_stn[-1]

par(mar = c(4, 5, 0, 0))
hist(rmse_stn, main = "", xlab = "RMSE", cex.axis = 2, cex.lab = 2, col = "grey", lwd = 2)

#variogram of residuals for model with and without correlation structure
a <- matrix(NA, nrow = 125, ncol = 118)
for (i in unique(glm.spl$stn)) {
  a[, as.numeric(i)] <- acf(residuals(asreml.null)[glm.spl$stn == i], plot=FALSE, lag.max = 125, na.action = na.pass)$acf
}
a <- a[, -1]

b <- matrix(NA, nrow = 125, ncol = 118)
for (i in unique(glm.spl$stn)) {
  b[, as.numeric(i)] <- acf(residuals(asreml.full)[glm.spl$stn == i], plot=FALSE, lag.max = 125, na.action = na.pass)$acf
}
b <- b[, -1]

pdf("C:/Users/Lisa/Dropbox/uni/MEE submitted/images/Figure_4.pdf")
par(mar = c(5, 5, 2, 2), mfrow = c(1, 1))
plot(seq(2, 250, by = 2), rowMeans(a, na.rm = T), type = "l", ylim = c(-0.5, 1), cex.axis = 2, cex.lab = 2,
     xlab = "depth lag (m)", ylab = "mean ACF", bty = "n", lwd = 2)
points(seq(2, 250, by = 2), rowMeans(b, na.rm = T), type = "l", lty = 2, lwd = 2)
legend(150, 1.1, c("null model", "full model"), lwd = 2, lty = c(1, 2), bty = "n", cex = 1.5, 
       pt.cex = 1, y.intersp = 0.8, seg.len=1.5, x.intersp=0.3)
dev.off()

#------------------------- PLOT FULL AND NULL MODEL ---------------------------#

#compares full and null model for salinity and temperature to show overfit splines
#graphs need to be rearanged in GIMP to fit journal size requirements

pdf("C:/Users/Lisa/Dropbox/uni/MEE submitted/images/Figure_5.pdf")

par(mfrow = c(2, 1))
par(oma = c(0, 5, 0, 0)) #make room for y-axis label

#salinity full
pred <- predict(asreml.full, classify = "sal", ignore = "z")
pval <- pred$predictions$pvals["predicted.value"]$predicted.value
z <- pred$predictions$pvals["sal"]$sal
se <- pred$predictions$pvals["standard.error"]$standard.error

logci <- pval + se%*%t(qnorm(c(0.025,0.5,0.975)))
ci <- exp(logci)
dimnames(ci)[[2]]<-c("lower95", "est", "upper95")

plot(z, ci[, 2], xlab = "salinity (psu)", ylab = "", lwd = 2, bty = "n",
     type = "l", ylim = c(min(ci[, 1]), max(ci[, 3]) + 0.5), cex.lab = 2, cex.axis = 2)
polygon(c(z, rev(z)), c(ci[, 1], rev(ci[, 3])), col = rgb(0.84, 0.84, 0.84, 1), border = NA)
points(z, ci[, 2], lwd = 2, type = "l")

#salinity null
pred <- predict(asreml.null, classify = "sal", ignore = "z")
pval <- pred$predictions$pvals["predicted.value"]$predicted.value
z <- pred$predictions$pvals["sal"]$sal
se <- pred$predictions$pvals["standard.error"]$standard.error

logci <- pval + se%*%t(qnorm(c(0.025,0.5,0.975)))
ci <- exp(logci)
dimnames(ci)[[2]]<-c("lower95", "est", "upper95")

points(z, ci[, 2], lwd = 2, col = "darkgrey", type = "l")
polygon(c(z, rev(z)), c(ci[, 1], rev(ci[, 3])), col = rgb(1, 0, 0, 0.3), border = NA)
points(z, ci[, 2], lwd = 2, type = "l", col = "red")

#temperature - full
pred <- predict(asreml.full, classify = "temp", ignore = "z")
pval <- pred$predictions$pvals["predicted.value"]$predicted.value
temp <- pred$predictions$pvals["temp"]$temp
se <- pred$predictions$pvals["standard.error"]$standard.error

logci <- pval + se%*%t(qnorm(c(0.025,0.5,0.975)))
ci <- exp(logci)
dimnames(ci)[[2]]<-c("lower95", "est", "upper95")

plot(temp, ci[, 2], xlab = expression(temperature~(~degree~C)), ylab = "", lwd = 2, bty = "n",
     type = "l", ylim = c(min(ci[, 1])-0.1, max(ci[, 3])), cex.lab = 2, cex.axis = 2)
polygon(c(temp, rev(temp)), c(ci[, 1], rev(ci[, 3])), col = rgb(0.84, 0.84, 0.84, 1), border = NA)
points(temp, ci[, 2], lwd = 2, type = "l")

#temperature - null
pred <- predict(asreml.null, classify = "temp", ignore = "z")
pval <- pred$predictions$pvals["predicted.value"]$predicted.value
temp <- pred$predictions$pvals["temp"]$temp
se <- pred$predictions$pvals["standard.error"]$standard.error

logci <- pval + se%*%t(qnorm(c(0.025,0.5,0.975)))
ci <- exp(logci)
dimnames(ci)[[2]]<-c("lower95", "est", "upper95")

points(temp, ci[, 2], lwd = 2, type = "l", col = "red")
polygon(c(temp, rev(temp)), c(ci[, 1], rev(ci[, 3])), col = rgb(1, 0, 0, 0.3), border = NA)
points(temp, ci[, 2], lwd = 2, type = "l", col = "red")

mtext(expression(hat(Fl) ~ (mu~g ~ L^{-1})), side = 2, outer = TRUE, line = 2, cex = 2)
dev.off()

#---------------- SUMMARY STATISTICS OF EXPLANATORY VARIABLES -----------------#

exp_var <- glm.spl[, c(2, 6:10)]

stats_tab <- matrix(round(apply(exp_var, 2, mean, na.rm = T), 2), ncol = 1)
stats_tab <- cbind(stats_tab, round(apply(exp_var, 2, sd, na.rm = T), 3))
stats_tab <- cbind(stats_tab, round(apply(exp_var, 2, min, na.rm = T), 2))
stats_tab <- cbind(stats_tab, round(apply(exp_var, 2, max, na.rm = T), 2))
colnames(stats_tab) <- c("mean", "sd", "min", "max")






