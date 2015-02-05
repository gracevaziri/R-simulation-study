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


#find depth of fluorescence maximum at each station
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
asreml.fit <- asreml(fixed = l.obs ~ z + par + temp:wm + oxy, random =~ spl(z, 10) + spl(par, 10) + 
                       spl(temp, 10):wm + spl(oxy, 10) + stn, 
                     data = glm.spl, rcov=~ ar1(z.fact):agau(x.fact, y.fact),
                     na.method.X = "include", workspace = 50000000)
asreml.fit <- update(asreml.fit)
summary(asreml.fit)

