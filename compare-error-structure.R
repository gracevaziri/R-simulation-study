#runs cross validation using BROKE-West data
#fits asreml model developed in simulation to BROKE-West data
#compares full model to null model without error structure
#author: Lisa-Marie Harrison
#date: 24/02/2015

setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/Data")
dat <- read.csv(file = "procCTD.csv", header= T)
library(asreml)
library(lattice)

#simplify names
names(dat) <- c("survey", "stn", "lat", "long", "start.time", "end.time", "depth", "transmittance", "cond", "temp", "sal", "par", "oxygen", "fluoro", "x2", "ice", "wm")

#remove null values
dat$sal[dat$sal == -9] <- NA
dat$temp[dat$temp == -9] <- NA
dat$par[dat$par == -9] <- NA
dat$fluoro[dat$fluoro == -9] <- NA

#compute log transformed fluoro values
dat$l.fluoro <- log(dat$fluoro)
dat$l.fluoro[is.nan(dat$l.fluoro)] <- NA

#get latitude and longitude for each station
n.station <- length(unique(dat$stn))
lat  <- dat$lat[duplicated(dat$stn) == FALSE]
long <- dat$long[duplicated(dat$stn) == FALSE]

#plot location of BROKE-West station with station number overlayed
plot(long, lat, col = "white", xlab = "longitude", ylab = "latitude")
text(long, lat, c(2:118))
title("location of BROKE-West CTD stations")

deg2rad <- function(deg) {
  #converts degrees to radians
  #input: degree coordinate
  #returns: radian coordinate 
  
  return(deg*pi/180)
}

gcd.hf <- function(lat1, long1, lat2, long2) {
  #calculates distance between two coordinates using the Haversine formula (hf)
  #input: radian latitude and longitude coordinates
  #returns: distance between coordinates in km
  
  R <- 6371 # Earth mean radius [km]
  delta.long <- (long2 - long1)
  delta.lat  <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1, sqrt(a)))
  d = R * c
  return(d) 
  
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
glm.spl <- data.frame(dat$l.fluoro, dat$depth, as.factor(dat$stn), rep(x, 1, each = length(unique(dat$depth))), rep(y, 1, each = length(unique(dat$depth))), dat$temp, dat$par, dat$sal, dat$oxygen, dat$ice, as.factor(dat$wm))
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


#fit full asreml model
asreml.full <- asreml(fixed = l.obs ~ z + par + temp:wm + oxy, random =~ spl(z, 10) + spl(par, 10) + 
                       spl(temp, 10):wm + spl(oxy, 10) + stn,
                     data = glm.spl, rcov=~ ar1(z.fact):agau(x.fact, y.fact),
                     na.method.X = "include", workspace = 50000000)
asreml.full <- update(asreml.full)
summary(asreml.full)$varcomp


#fit null asreml model
asreml.null <- asreml(fixed = l.obs ~ z + par + temp:wm + oxy, random =~ spl(z, 10) + spl(par, 10) + 
                        spl(temp, 10):wm + spl(oxy, 10) + stn, data = glm.spl,
                      na.method.X = "include", workspace = 50000000)
summary(asreml.null)$varcomp


#plot fitted values
plot(asreml.full$fitted, asreml.null$fitted)

par(mfrow = c(1, 2))
hist(asreml.full$residuals, xlim = c(-4, 4))
hist(asreml.null$residuals)

#plot standardised conditional residuals
full_cond_res <- asreml.full$aom$R[, 2]
null_cond_res <- asreml.null$aom$R[, 2]
plot(full_cond_res, null_cond_res)

plot(asreml.full)
plot(asreml.null)


#likelihood ratio test to check whether adding correlation structure works
1 - pchisq(2 * (asreml.full$loglik - asreml.null$loglik), 1) 




