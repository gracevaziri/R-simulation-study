#simulates an irregular grid with anisotropic gaussian correlation in the x and y directions, 
#ar1 correlation in the z direction
#uses BROKE-West ctd station coordinates
library(asreml)
library(mgcv)
library(fields)
dat <- read.csv("C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/Data/stn_coordinates.csv", header = T)

#get latitude and longitude for each station
lat  <- dat$latitude
long <- dat$longitude

n.station <- length(lat)
mu <- 50
sd <- 40
noise.sd <- 0.2 #noise sd
stn.sd   <- 0.05 #sd for random station effect
z.phi  <- 0.4 #ar1 autocorrelation down z
x.phi <- 0.5 #gaussian autocorrelation across x
y.phi <- 0.4 #gaussian autocorrelation across y

mult <- 1e3
z <- seq(0, 250, 5) #explanatory variable (depth)
z.int <- rep(c(1:length(z)), n.station) #explanatory variable (depth)
stn <- rep(c(1:n.station), 1, each = length(z))
rho <- mult*dnorm(z, mu, sd)/(pnorm(max(z), mu, sd) - pnorm(min(z), mu, sd))
stn.re <- rnorm(n.station, mean = 0, sd = stn.sd) #station specific random effect

#random noise matrix
r.noise <- rnorm(length(lat)*length(z), 0, noise.sd)

#function to convert degrees to radians
deg2rad <- function(deg) {
  return(deg*pi/180)
}

#function to calculate the distance between two points with radian lat/long 
#Uses the Haversine formula
gcd.hf <- function(lat1, long1, lat2, long2) {
  
  R <- 6371 # Earth mean radius (km)
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  d = R * c
  return(d) # Returns distance (km)
  
}

#calculate distance of each station from each other station in x-direction
dist_x <- matrix(0, ncol = n.station, nrow = n.station)
for (i in 1:n.station) {
  for (k in 1:n.station) {
    dist_x[i, k] <- gcd.hf(deg2rad(lat[i]), deg2rad(long[i]), deg2rad(lat[k]), deg2rad(long[i]))/100   
  }
}

#calculate distance of each station from each other station in y-direction
dist_y <- matrix(0, ncol = n.station, nrow = n.station)
for (i in 1:n.station) {
  for (k in 1:n.station) {
    dist_y[i, k] <- gcd.hf(deg2rad(lat[i]), deg2rad(long[i]), deg2rad(lat[i]), deg2rad(long[k]))/100
  }
}

#get distance of each station from station 1 in x and y directions
x <- dist_x[1, ]
y <- dist_y[1, ]

#create distance matrix for x and y directions using distance from station 1
adist_x <- matrix(0, ncol = n.station, nrow = n.station)
for (i in 1:n.station) {
  for (j in 1:n.station) {    
    adist_x[i, j] <- (x[i] - x[j])^2    
  }
}

adist_y <- matrix(0, ncol = n.station, nrow = n.station)
for (i in 1:n.station) {
  for (j in 1:n.station) {    
    adist_y[i, j] <- (y[i] - y[j])^2     
  }
}

#change distances into a distance object
adist_x <- as.matrix(as.dist(adist_x, diag = FALSE, upper = FALSE))
adist_y <- as.matrix(as.dist(adist_y, diag = FALSE, upper = FALSE))

#create the correlation structure
omega1 <- (x.phi^adist_x ) * (y.phi^adist_y)

#calculate correlation weights, and invert weights matrix
weights <- chol(solve(omega1))
weights_inv <- solve(weights)

#form the combined correlated error component
t.cor <- rep(0, length(r.noise))
for (k in 2:length(z)) {
  z.data <- matrix(r.noise[z.int == k], ncol = 1) #random noise for all stations at depth k
  xy_error <- weights_inv %*% z.data
  for (i in 1:n.station) {
    w <- which(stn == i & z.int == k)
    t.cor[w] <- t.cor[stn == i & z.int == (k - 1)]*z.phi + xy_error[i] #first component is ar1(z), second is agau(x, y)
  }
}

#calculate the total observations
l.obs <- rep(log(rho), n.station) + t.cor + rep(stn.re, 1, each = length(rho))
obs <- exp(l.obs)


#data frame
glm.spl <- data.frame(obs, l.obs, rep(z, n.station), as.factor(rep(c(1:n.station), 1, each = length(z))), rep(x, 1, each = 51), rep(y, 1, each = 51))
names(glm.spl) <- c("obs", "l.obs", "z", "stn", "x", "y")
glm.spl$z.fact <- as.factor(as.integer(glm.spl$z))
glm.spl$x.fact <- as.factor(glm.spl$x)
glm.spl$y.fact <- as.factor(glm.spl$y)
glm.spl <- glm.spl[order(glm.spl$z, glm.spl$x, glm.spl$y), ] #sort by order of rcov structure


#---------------------------------- asreml -------------------------------------#

asreml.fit <- asreml(fixed = l.obs ~ z, random =~ spl(z) + stn, data = glm.spl, 
                     splinepoints = list(z = seq(0, 250, 25)), rcov=~ ar1(z.fact):agau(x, y))
summary(asreml.fit)


vals <- matrix(c(stn.sd, noise.sd, z.phi, x.phi, y.phi, round(summary(asreml.fit)$varcomp[2,2]^0.5, 2), round(summary(asreml.fit)$varcomp[3,2]^0.5, 2), round(summary(asreml.fit)$varcomp[4,2], 2), round(summary(asreml.fit)$varcomp[5,2], 2), round(summary(asreml.fit)$varcomp[6,2], 2)), ncol = 2)
colnames(vals) <- c("true", "fitted")
rownames(vals) <- c("stn", "noise", "z ar1", "x agau", "y agau")
vals


#plot fitted against observed
plot(l.obs[1:51])
points(fitted(asreml.fit)[glm.spl$stn == 1], col = "red")



