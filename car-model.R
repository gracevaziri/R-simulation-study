#car process down z
library(asreml)
library(mgcv)

n.station <- 2500 #specify number of stations to run
mu <- 50
sd <- 40
noise.sd <- 0.2 #noise sd for ar process
stn.sd   <- 0.7 #sd for random station effect
phi.true <- 0.6 #CAR autocorrelation down depths

mult <- 1e3
z <- seq(0, 250, 5) #explanatory variable (depth)
rho <- mult*dnorm(z, mu, sd)/(pnorm(max(z), mu, sd) - pnorm(min(z), mu, sd))
stn.re <- rnorm(n.station, mean = 0, sd = stn.sd) #station specific random effect

#CAR1 autoregressive process down depths
l.obs <- NULL
for(i in c(1:n.station)) {
  start.vals <- rnorm(length(rho), 0, 0.2)
  l.temp <- rep(0, 51)
  for (j in c(2:(length(rho) - 1))){
    l.temp[j] <- start.vals[j] + start.vals[j - 1]*(phi.true/2) +  start.vals[j + 1]*(phi.true/2) +  stn.re[i]
  }
  l.temp[1] <- start.vals[1] + start.vals[2]*(phi.true) +  stn.re[i]
  l.temp[51] <- start.vals[51] + start.vals[50]*(phi.true) +  stn.re[i]
  l.temp <- log(rho) + l.temp
  l.obs <- c(l.obs,l.temp)
}
obs <- exp(l.obs)

#data frame
glm.spl <- data.frame(obs, l.obs, rep(z, n.station), rep(c(1:n.station), 1, each = length(z)), as.factor(z))
names(glm.spl) <- c("obs", "l.obs", "z", "stn", "z.fact")
glm.spl$stn <- as.factor(glm.spl$stn)
glm.spl$z.fact <- as.factor(as.integer(glm.spl$z.fact))


#---------------------------------- asreml ---------------------------------------#

asreml.fit <- asreml(fixed = l.obs ~ z, random =~ spl(z) + stn, data = glm.spl, 
                     splinepoints = list(z = seq(0, 250, 25)), rcov =~ stn:sar(z.fact))

summary(asreml.fit)
summary(asreml.fit)$varcomp[2,2]^0.5
summary(asreml.fit)$varcomp[3,2]^0.5
summary(asreml.fit)$varcomp[4,2]

#----------------------------------- gamm ---------------------------------------#

#make variable for time series distance rather than using depths
glm.spl$z.int <- as.integer(glm.spl$z.fact)

#create lower triangular correlation matrix
mat <- matrix(0, nrow = 51, ncol = 51)
for (i in 2:nrow(mat)) {
  j <- i - 1
  mat[i, j] <- phi.true/2
}
val <- mat[lower.tri(mat, diag = FALSE)]
cs1Symm <- corSymm(value = val, form = ~ z.int | stn)
cs1Symm <- Initialize(cs1Symm, data = glm.spl)

#fit gamm
gam.fit <- gamm(l.obs ~ s(z),  random=list(stn =~1), data = glm.spl, 
                correlation = cs1Symm)

summary(gam.fit$lme)
summary(gam.fit$gam)
intervals <- intervals(gam.fit$lme)
phi <- intervals$corStruct[1,2]























