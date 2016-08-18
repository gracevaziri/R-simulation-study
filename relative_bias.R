dat <- read.csv("C:/Users/43439535/Documents/Lisa/phd/Mixed models/Data/3d_ar1_agau_irregular_grid_estimates.csv", header = T)


#true input values
par_vals <- c("stn.sd" = 0.22, "noise.sd" = 0.45, "z.phi" = 0.35,
              "x.phi" = 0.5, "y.phi" = 0.4)

#estimator means
est_mean <- apply(dat, 2, mean)

#relative bias
est_mean - par_vals

#confidence intervals


#coverage probabilities