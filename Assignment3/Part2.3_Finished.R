library(ctsmTMB)

dat <- read.csv("ex3_largecase.csv", stringsAsFactors = FALSE)

dat_ev <- subset(dat, Event_ID == 1)
t0 <- as.POSIXct(dat_ev$Timestamp[1])
dat_ev$t <- as.numeric(difftime(as.POSIXct(dat_ev$Timestamp),
                                t0,
                                units = "hours"))

scale_factor <- 1/1000  
dat_ctsm <- data.frame(
  t = dat_ev$t,
  R = dat_ev$Rainfall,       
  P = dat_ev$Pumpflow * scale_factor * 60,  
  S = dat_ev$Volume * scale_factor,
  Y = dat_ev$Volume * scale_factor
)
dat_ctsm$R <- dat_ctsm$R / max(dat_ctsm$R)

###################
##### Model #######
###################
model <- ctsmTMB$new()
model$addSystem(
  dG ~ A * R * dt - 3/K1 * G * dt  + sigma * dw1,
  dC ~ 3/K1 * G * dt - 3/K2 * C * dt
        - (1/(1+exp(-alpha*(C - beta)))) * C * dt + sigma * dw2,
  dT ~ (1/(1+exp(-alpha*(C - beta)))) * C * dt - 3/K3 * T * dt 
        + sigma * dw3,
  dS ~ (3/K3) * T * dt - P * dt + sigma * dw4
  )

model$addObs(Y ~ S)
model$setVariance(Y ~ sigma_e^2)

model$addInput(R)
model$addInput(P)


# Parameters
model$setParameter(
  A       = c(initial = 5,    lower = 0,    upper = 50),  
  K1      = c(initial = 1, lower = 1e-3, upper = 500),
  K2      = c(initial = 1, lower = 1e-3, upper = 500),
  K3      = c(initial = 1, lower = 1e-3, upper = 500),
  #K4      = c(initial = 1, lower = 1e-3, upper = 500),
  alpha = c(initial = 1, lower = 0.01, upper = 50),
  beta  = c(initial = 0.5, lower = 0.01, upper = 50),
  sigma  = c(initial = 0.1, lower = 1e-4, upper = 30),
  #sigma2  = c(initial = 0.1, lower = 1e-4, upper = 30),
  #sigma3  = c(initial = 0.1, lower = 1e-4, upper = 30),
  #sigma4  = c(initial = 0.1, lower = 1e-4, upper = 30),
  sigma_e = c(initial = 0.1, lower = 1e-4, upper = 30)
)

# Initial States
model$setInitialState(
  list(c(dat_ctsm$R[1], dat_ctsm$P[1], dat_ctsm$P[1], dat_ctsm$Y[1]), 1e-1 * diag(4))
)

# Estimate
fit <- model$estimate(dat_ctsm)
print(fit)
cat("logLik:", -fit$nll, "\n")
cat("A:", fit$par.fixed["A"], "\n")


# Print Results
print(fit)
cat("logLik:", -fit$nll, "\n")
cat("Estimated Parameters:\n")
print(fit$par.fixed)


# Predict
pred_1k <-model$predict(dat_ctsm,k.ahead=3)

# Set up a 2-row layout
par(mfrow = c(3, 1))
y_range <- range(c(dat_ctsm$y, pred_1k$states$S), na.rm = TRUE)
y_range <- c(0, 30)

# Plot 1: Rainwater over time (black line)
plot(dat_ctsm$t, dat_ctsm$R, type = "l", col = "blue",  # Use common y-axis range
     ylab = "Stormwater", xlab = "Time (hours)", 
     main = "Observed Rainwater over time (Event 1)", xaxt = "n")
axis(1, at = seq(0, max(dat_ctsm$t), by = 2), labels = seq(0, max(dat_ctsm$t), by = 2))

# Plot 2: Stormwater over time (black line)
plot(dat_ctsm$t, dat_ctsm$Y, type = "l", col = "black", 
     ylim = y_range,  # Use common y-axis range
     ylab = "Stormwater", xlab = "Time (hours)", 
     main = "Observed Stormwater over time (Event 1)", xaxt = "n")
axis(1, at = seq(0, max(dat_ctsm$t), by = 2), labels = seq(0, max(dat_ctsm$t), by = 2))

# Plot 3: Predicted Stormwater (red line)
plot(pred_1k$states$t.j, pred_1k[["observations"]][["Y"]], type = "l", col = "red", lty = 2, 
     ylim = y_range,  # Use the same y-axis range
     ylab = "Stormwater", xlab = "Time (hours)", 
     main = "Predicted Stormwater over time (Event 1)", xaxt = "n")
axis(1, at = seq(0, max(dat_ctsm$t), by = 2), labels = seq(0, max(dat_ctsm$t), by = 2))


res <- fit$residuals
plot(fit)


