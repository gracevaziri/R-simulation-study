#runs a simulation multiple times using the same input values to check the average
#uses ar1 process
library(asreml)
source(file = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/R code/R-simulation-study/simData_2d_ar1.R")

#set all values to zero initially
stn.sd   <- 0
noise.sd <- 0
z.phi    <- 0
x.phi    <- 0

for (i in 1:200) {
  
  glm.spl <- simData(n.station = 100, noise.sd = 0.2, stn.sd = 0.1, z.phi = 0.45,
                     x.phi = 0.8)
  
  
  asreml.fit <- asreml(fixed = l.obs ~ z, data = glm.spl, 
                       splinepoints = list(z = seq(0, 250, 25)),
                       random= ~ spl(z) + idv(y.fact):ar1(x.fact), rcov =~ y.fact:x.fact:ar1(z.fact))
  
  #extract variance components
  stn.sd[i]   <- summary(asreml.fit)$varcomp[2,2]^0.5
  noise.sd[i] <- summary(asreml.fit)$varcomp[4,2]^0.5
  z.phi[i]    <- summary(asreml.fit)$varcomp[3,2]
  x.phi[i]    <- summary(asreml.fit)$varcomp[5,2]
  
}

#plot histogram of every variance component with fitted average and true values overlayed
par(mfrow = c(2, 3))
hist(stn.sd, main = "stn.sd")
abline(v = 0.1, col = "red", lwd = 2)
abline(v = mean(stn.sd), col = "blue", lwd = 2)
hist(noise.sd, main = "noise.sd", xlim = c(0.19, 0.24))
abline(v = 0.2, col = "red", lwd = 2)
abline(v = mean(noise.sd), col = "blue", lwd = 2)
hist(z.phi, main = "z.phi")
abline(v = 0.45, col = "red", lwd = 2)
abline(v = mean(z.phi), col = "blue", lwd = 2)
hist(x.phi, main = "x.phi")
abline(v = 0.8, col = "red", lwd = 2)
abline(v = mean(x.phi), col = "blue", lwd = 2)


