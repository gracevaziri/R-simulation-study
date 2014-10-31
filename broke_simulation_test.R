#fits asreml model developed in simulation to BROKE-West data
#uses the distance between stations rather than latitude and longitude
#author: Lisa-Marie Harrison
#date: 18/09/2014

setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/Data")
dat.cut <- read.csv(file = "rstnCTD.csv", header= T)
library(asreml)
library(nlme)
library(lattice)
library(mgcv)

#remove null values
dat.cut$sal[dat.cut$sal == -9] <- NA
dat.cut$temp[dat.cut$temp == -9] <- NA
dat.cut$par[dat.cut$par == -9] <- NA
dat.cut$fluoro[dat.cut$fluoro == -9] <- NA

#get latitude and longitude for each station
n.station <- length(unique(dat.cut$stn))
lat  <- dat.cut$lat[duplicated(dat.cut$stn) == FALSE]
long <- dat.cut$long[duplicated(dat.cut$stn) == FALSE]

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

dist_x <- matrix(0, ncol = n.station, nrow = n.station)
for (i in 1:n.station) {
  for (k in 1:n.station) {
    dist_x[i, k] <- gcd.hf(deg2rad(lat[i]), deg2rad(long[i]), deg2rad(lat[k]), deg2rad(long[i]))/100
  }
}

dist_y <- matrix(0, ncol = n.station, nrow = n.station)
for (i in 1:n.station) {
  for (k in 1:n.station) {
    dist_y[i, k] <- gcd.hf(deg2rad(lat[i]), deg2rad(long[i]), deg2rad(lat[i]), deg2rad(long[k]))/100
  }
}

#get distance of each station from station 1 in x and y directions
x <- dist_x[1, ]
y <- dist_y[1, ]


#data frame
glm.spl <- data.frame(log(dat.cut$fluoro - min(na.omit(dat.cut$fluoro))), dat.cut$profile.depth, as.factor(dat.cut$stn), rep(x, 1, each = length(unique(dat.cut$profile.depth))), rep(y, 1, each = length(unique(dat.cut$profile.depth))), dat.cut$temp, dat.cut$par, dat.cut$sal, dat.cut$oxy, dat.cut$ice, as.factor(dat.cut$water.mass))
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

#------------------------------- FIT ASREML MODELS -----------------------------------#

#fit asreml model
asreml.fit <- asreml(fixed = l.obs ~ z + par + temp:wm + ice , random =~ spl(z, 10) + spl(par, 10) + 
                        spl(temp, 10):wm +  spl(ice, 10) + stn, 
                     data = glm.spl, rcov=~ ar1(z.fact):agau(x.fact, y.fact),
                     na.method.X = "include", workspace = 50000000)
asreml.fit <- update(asreml.fit)
summary(asreml.fit)



#plot fitted against observed for all stations
lat.plot <- xyplot(glm.spl$l.obs + fitted(asreml.fit) ~ glm.spl$z | glm.spl$stn, 
                   outer = FALSE, type = "l", xlab = "depth (m)", ylab = "l.fluoro")
update(lat.plot, par.settings = simpleTheme(lwd = c(2, 1), col = c("dodgerblue", "red")))

#fit the same asreml model but without the correlation structure
fit <- asreml(fixed = l.obs ~ z + par + temp:wm + ice , random =~ spl(z, 10) + spl(par, 10) + 
                       spl(temp, 10):wm +  spl(ice, 10) + stn, 
                     data = glm.spl, na.method.X = "include", workspace = 50000000)

b.int <- fit$coefficients$fixed["(Intercept)"]
temp <- fit$coefficients$fixed[c(1:4)]
ice <- fit$coefficients$fixed["ice"]
par <- fit$coefficients$fixed["par"]
z <- fit$coefficients$fixed["z"]
wm <- as.numeric(glm.spl$wm) #water mass indicator variable

n <- fit$coefficients$random
z_spl <- n[which(substr(names(n), start = 1, stop = 5) %in% "spl(z")]
par_spl <- n[which(substr(names(n), start = 1, stop = 5) %in% "spl(p")]
ice_spl <- n[which(substr(names(n), start = 1, stop = 5) %in% "spl(i")]
temp_spl <- n[which(substr(names(n), start = 1, stop = 5) %in% "spl(t")]

pred.fixed <- b.int + z*depth + ice*i + par*light + temp[wm]*temp
pred.random <- depth*1

#----------------------- PREDICT VALUES AT OTHER STATIONS ---------------------#

#predict at each point using parameter values for other stations to test extrapolation
#using depth, par, temperature, ice and watermass, predict l.fluoro using the model
extra.dat <- read.csv(file = "ctd_data.csv", header = T)
extra.dat <- extra.dat[extra.dat$stn %in% c(3, 13, 23, 33, 43, 53, 63, 73, 83, 93), ] #use 10 stations that weren't used to fit the model

pval <- 0
se <- 0
for(i in 1:nrow(extra.dat)) {

  capture.output({
    pred <- predict(asreml.fit, classify = "wm:ice:par:temp:z", levels = list("z" = extra.dat$profile.depth[i], "temp" = extra.dat$temp[i], "ice" = extra.dat$ice[i], "par" = extra.dat$par[i], "wm" = extra.dat$water_mass[i]))
  }, file = tempfile())
  pval[i] <- pred$predictions$pvals["predicted.value"]$predicted.value
  se[i] <- pred$predictions$pvals["standard.error"]$standard.error
    
  if(i %% 100 == 0) print(i)
  
}  
#write.csv(cbind(pval, se), "pred_with_correlation.csv", row.names = F)

#compare predicted values to actual values
plot(extra.dat$l.fluoro, pval, xlab = "observed", ylab = "predicted")
title("Observed vs predicted for extra stations")

lat.plot <- xyplot(extra.dat$l.fluoro + pval ~ extra.dat$profile.depth | extra.dat$stn, 
                   outer = FALSE, type = "l", xlab = "depth (m)", ylab = "l.fluoro")
update(lat.plot, par.settings = simpleTheme(lwd = c(2, 1), col = c("dodgerblue", "red")))
title("Observed (blue) vs predicted (red) by station")


#--------------------------- CHECK AUTOCORRELATION ----------------------------#

#variogram of residuals for model with and without correlation strucutre
d <- 62
gamma <- asreml.variogram(glm.spl$z[glm.spl$stn == d], z = resid(asreml.fit)[glm.spl$stn == d])$gamma
dist  <- asreml.variogram(glm.spl$z[glm.spl$stn == d], z = resid(asreml.fit)[glm.spl$stn == d])$x
plot(dist, gamma)

gamma <- asreml.variogram(glm.spl$z[glm.spl$stn == d], z = resid(fit)[glm.spl$stn == d])$gamma
dist  <- asreml.variogram(glm.spl$z[glm.spl$stn == d], z = resid(asreml.fit)[glm.spl$stn == d])$x
plot(dist, gamma)

#likelihood ratio test to check whether adding correlation structure works
1 - pchisq(2 * (asreml.fit$loglik - fit$loglik), 1) 



#---------------------------------- FIT GAMM -----------------------------------------#

#fit gamm to compare
gamm.fit <- gamm(l.obs ~ s(z) + s(temp) + s(sal) + s(par) + s(oxy) + s(ice), random = list(stn =~ 1, x =~1, y =~1), 
                 data = glm.spl, correlation = corAR1(0.9, z ~ 1))
summary(gamm.fit$gam)


#plot fitted against observed by station
lat.plot <- xyplot(glm.spl$l.obs + fitted(gamm.fit$gam) ~ glm.spl$z | glm.spl$stn, outer = FALSE, type = "l")
update(lat.plot, par.settings = simpleTheme(lwd = c(2, 1), col = c("dodgerblue", "red")))


