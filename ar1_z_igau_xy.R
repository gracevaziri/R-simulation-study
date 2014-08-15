#simulates a regular grid with anisotropic exponential correlation in the x and y directions, 
#ar1 correlation in the z direction
#author: Lisa-Marie Harrison
#date: 14/08/2014
library(asreml)
library(mgcv)
library(fields)

n.station <- 100 #specify number of stations to run
mu <- 50
sd <- 40
noise.sd <- 0.2 #noise sd for ar process
stn.sd   <- 0.05 #sd for random station effect
z.phi <- 0.4 #autocorrelation (ar1) down depths
xy.phi <- 0.7 #exponential autocorrelation across xy

mult <- 1e3
z <- seq(0, 250, 5) #explanatory variable (depth)
z.int <- rep(c(1:length(z)), n.station) #explanatory variable (depth)

stn <- rep(c(1:n.station), 1, each = length(z))
rho <- mult*dnorm(z, mu, sd)/(pnorm(max(z), mu, sd) - pnorm(min(z), mu, sd))
stn.re <- rnorm(n.station, mean = 0, sd = stn.sd) #station specific random effect

#regular grid for 100 stations
x <- sort(rep(sort(rep(seq(1:sqrt(n.station)), length(rho))), sqrt(n.station)))
y <- rep(sort(rep(seq(1:sqrt(n.station)), length(rho))), sqrt(n.station))

#random noise matrix
r.noise <- rnorm(length(x), 0, noise.sd)

cols <- 10
rows <- 10
#create dataframe (xyz format)
X <- rep(1:10, each = rows)
Y <- rep(1:10, cols)

# create a spatial autocorrelation signature
# coordinate list
coords <- data.frame(X, Y)
# distance matrix

loc <- matrix(c(1:100), ncol = 10, byrow = T)
dist <- matrix(0, ncol = 100, nrow = 100)
for (i in 1:100) {
  for (k in 1:100) {
    w_current <- which(loc == i, arr.ind = TRUE)      
    w_new <- which(loc == k, arr.ind = TRUE)
    dist[i, k] <- (w_current[1] - w_new[1])^2 + (w_current[2] - w_new[2])^2     
  }
}

# create a correlation structure (exponential)
omega1 <- xy.phi^dist
# calculate correlation weights, and invert weights matrix
weights <- chol(solve(omega1))
weights_inv <- solve(weights)

#combined correlated error component
t.cor <- rep(0, length(r.noise))
for (k in 2:length(z)) {
  z.data <- matrix(r.noise[z.int == k], ncol = max(x))
  xy_error <- matrix(weights_inv %*% matrix(z.data, ncol = 1, byrow = T), ncol = 10, byrow = T)
  for (i in 1:max(y)) {
    for (j in 1:max(x)) {
      w <- which(y == i & z.int == k & x == j)
      t.cor[w] <- t.cor[y == i & z.int == (k - 1) & x == j]*z.phi + xy_error[i, j]
    }
  }
}


#total observations
l.obs <- rep(log(rho), n.station) + t.cor + rep(stn.re, 1, each = length(rho))
obs <- exp(l.obs)

#data frame
glm.spl <- data.frame(obs, l.obs, rep(z, n.station), as.factor(rep(c(1:n.station), 1, each = length(z))), x, y)
names(glm.spl) <- c("obs", "l.obs", "z", "stn", "x", "y")
glm.spl$z.fact <- as.factor(as.integer(glm.spl$z))
glm.spl$x.fact <- as.factor(as.integer(glm.spl$x))
glm.spl$y.fact <- as.factor(as.integer(glm.spl$y))
glm.spl <- glm.spl[order(glm.spl$z, glm.spl$x), ] #sort by order of rcov structure


#---------------------------------- asreml -------------------------------------#

asreml.fit <- asreml(fixed = l.obs ~ z, random =~ spl(z) + stn, data = glm.spl, 
                     splinepoints = list(z = seq(0, 250, 25)), rcov=~ ar1(z.fact):igau(x.fact, y.fact), maxIter = 50)
summary(asreml.fit)


vals <- matrix(c(stn.sd, noise.sd, z.phi, xy.phi, round(summary(asreml.fit)$varcomp[2,2]^0.5, 2), round(summary(asreml.fit)$varcomp[3,2]^0.5, 2), round(summary(asreml.fit)$varcomp[4,2], 2), round(summary(asreml.fit)$varcomp[5,2], 2)), ncol = 2)
colnames(vals) <- c("true", "fitted")
rownames(vals) <- c("stn", "noise", "z ar1", "xy iexp")
vals


#plot fitted against observed
plot(l.obs[1:51])
points(fitted(asreml.fit)[glm.spl$stn == 1], col = "red")








