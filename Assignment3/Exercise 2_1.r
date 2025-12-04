# 1. Create a user library path in your home folder
user_lib <- file.path(
  Sys.getenv("USERPROFILE"),  # e.g. "C:/Users/yourname"
  "R",
  "win-library",
  paste(R.version$major, R.version$minor, sep = ".")
)

dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)

# 2. Make this the first library in .libPaths()
.libPaths(c(user_lib, .libPaths()))

# 3. Check that the first path is NOT under Program Files anymore
.libPaths()

install.packages("ctsmTMB")

#-------------

library(ctsmTMB)

dat <- read.csv("C:/Users/snehi/Documents/Advanced TSA/Assignment 3/ex1_rainfallrunoff.csv",
                stringsAsFactors = FALSE)

#convert the timestamp to hours
dat$t <- as.numeric(difftime(as.POSIXct(dat$timestamp),
                             as.POSIXct(dat$timestamp[1]),
                             units = "hours"))

data_ctsm <- data.frame(t = dat$t,
                        U = dat$rainfall,
                        stormwater = dat$stormwater)
#Question 2.1.1
#initialise the model
model <- ctsmTMB$new()

model$addSystem(dX1 ~ (A * U - (1 / K) * X1) * dt + sigma * dw1)
model$addSystem(dX2 ~ ((1 / K) * X1) * dt + sigma * dw2)

model$addObs(stormwater ~ X2)
model$setVariance(stormwater ~ sigma_e^2)

model$addInput(U)

model$setParameter(
  A       = c(initial = 10,  lower = 1e-4, upper = 1e5),
  K       = c(initial = 2,   lower = 1e-3, upper = 1e3),
  sigma   = c(initial = 0.1, lower = 1e-6, upper = 10),
  sigma_e = c(initial = 0.1, lower = 1e-6, upper = 10)
)

init <- list(mean = c(X1 = 0, X2 = 0), cov = 1e-2 * diag(2))
model$setInitialState(initial.state = init)

fit <- model$estimate(data_ctsm)

print(fit)
cat("logLik:", -fit$nll, "\n")
cat("A:", fit$par.fixed["A"], "\n")
cat("K:", fit$par.fixed["K"], "\n")

#k-step predictions 
pred <- model$predict(data_ctsm)   
plot(pred, type = "observations", against = "stormwater")

#Question 2.1.2
#helper, build and fit model with m states
fit_linear_reservoir <- function(m, data) {
  stopifnot(m >= 1)
  model <- ctsmTMB$new()
  
  #state equations
  for (i in 1:m) {
    if (i == 1) {
      #first reservoir
      eq <- dX1 ~ (A * U - (1 / K) * X1) * dt + sigma * dw1
    } else if (i < m) {
      #intermediate reservoirs
      eq <- as.formula(
        sprintf("dX%d ~ ((1 / K) * X%d - (1 / K) * X%d) * dt + sigma * dw%d",
                i, i - 1, i, i)
      )
    } else {
      #last reservoir– no outflow term
      eq <- as.formula(
        sprintf("dX%d ~ ((1 / K) * X%d) * dt + sigma * dw%d",
                i, i - 1, i)
      )
    }
    model$addSystem(eq)
  }
  
  #observation:last state
  obs_eq <- as.formula(sprintf("stormwater ~ X%d", m))
  model$addObs(obs_eq)
  model$setVariance(stormwater ~ sigma_e^2)
  
  # input
  model$addInput(U)
  
  #parameters-same across models
  model$setParameter(
    A       = c(initial = 10,  lower = 1e-4, upper = 1e5),
    K       = c(initial = 2,   lower = 1e-3, upper = 1e3),
    sigma   = c(initial = 0.1, lower = 1e-6, upper = 10),
    sigma_e = c(initial = 0.1, lower = 1e-6, upper = 10)
  )
  
  #initial state
  state_names <- paste0("X", 1:m)
  init <- list(
    mean = setNames(rep(0, m), state_names),
    cov  = 1e-2 * diag(m)
  )
  model$setInitialState(initial.state = init)
  
  fit <- model$estimate(data)
  
  list(model = model, fit = fit)
}

#going up until n=9
n_states_vec <- 2:9

fits <- lapply(n_states_vec, function(m) fit_linear_reservoir(m, data_ctsm))

#model fit stats
fit_stats <- lapply(seq_along(fits), function(i) {
  m   <- n_states_vec[i]
  fit <- fits[[i]]$fit
  
  logLik <- -fit$nll
  k      <- length(fit$par.fixed)    # number of estimated parameters
  n      <- nrow(data_ctsm)          # sample size
  
  AIC <- -2 * logLik + 2 * k
  BIC <- -2 * logLik + log(n) * k
  
  data.frame(
    n_states = m,
    logLik   = logLik,
    k        = k,
    AIC      = AIC,
    BIC      = BIC
  )
})

model_sel_tab <- do.call(rbind, fit_stats)
print(model_sel_tab)

#now printing various params
param_tables <- lapply(seq_along(fits), function(i) {
  m   <- n_states_vec[i]
  fit <- fits[[i]]$fit
  
  est <- fit$par.fixed
  se  <- fit$sd.fixed
  
  # keep only the four main parameters (just in case par.fixed has more)
  keep <- c("A", "K", "sigma", "sigma_e")
  keep <- keep[keep %in% names(est)]
  
  est <- est[keep]
  se  <- se[keep]
  
  tval <- est / se
  pval <- 2 * (1 - pnorm(abs(tval)))
  
  data.frame(
    n_states = m,
    param    = keep,
    estimate = as.numeric(est),
    se       = as.numeric(se),
    t        = as.numeric(tval),
    p        = as.numeric(pval)
  )
})

param_tab_long <- do.call(rbind, param_tables)
print(param_tab_long)


#Question 2.1.3
# Best model is the 4-state model (m = 4)
fit4 <- fit_linear_reservoir(4, data_ctsm)$fit
fit2 <- fit_linear_reservoir(2, data_ctsm)$fit 

#Extract A and K for both models
pars_2  <- fit2$par.fixed[c("A", "K")]
pars_4  <- fit4$par.fixed[c("A", "K")]

pars_2
pars_4

#Question 2.1.4
#fit the best (4-state) model
fit4_obj <- fit_linear_reservoir(4, data_ctsm)
fit4 <- fit4_obj$fit

#correlation matrix of parameters
summary(fit4, correlation = TRUE)

#extract correlation matrix explicitly
cor_mat <- fit4$cov.fixed
sd_vec  <- fit4$sd.fixed
cor_matrix <- cor_mat / (sd_vec %o% sd_vec)
print(round(cor_matrix, 3))

residuals4 <- fit4$residuals$normalized$stormwater   # <– this is the right slot
residuals4 <- as.numeric(na.omit(residuals4))
length(residuals4)

#Q–Q plot
qqnorm(residuals4, main = "Q-Q plot of residuals (4-state model)")
qqline(residuals4, col = "red", lwd = 2)

#acf and pacf
par(mfrow = c(1, 2))
acf(residuals4, main = "ACF of residuals (4-state model)")
pacf(residuals4, main = "PACF of residuals (4-state model)")
par(mfrow = c(1, 1))

#Question 2.1.5
#fitted 4-state parameters
fit4_obj <- fit_linear_reservoir(4, data_ctsm)
fit4     <- fit4_obj$fit

A_hat      <- fit4$par.fixed["A"]
K_hat      <- fit4$par.fixed["K"]
sigma_hat  <- fit4$par.fixed["sigma"]
sigmae_hat <- fit4$par.fixed["sigma_e"]

A_hat; K_hat; sigma_hat; sigmae_hat

t  <- data_ctsm$t
U  <- data_ctsm$U            #rainfall input
n  <- length(t)
m  <- 4                      #number of states = 4 (best model)

dt <- c(0, diff(t))          #time step in hours

#simulate 4-state linear reservoir)
set.seed(123)

X     <- matrix(0, nrow = n, ncol = m)   # states X1,...,X4
Y_sim <- numeric(n)                      # simulated stormwater

for (i in 2:n) {
  dW <- rnorm(m, mean = 0, sd = sqrt(dt[i]))  # Wiener increments
  
  #first reservoir
  X[i, 1] <- X[i-1, 1] +
    (A_hat * U[i-1] - (1 / K_hat) * X[i-1, 1]) * dt[i] +
    sigma_hat * dW[1]
  
  #intermediate reservoirs
  if (m > 2) {
    for (j in 2:(m - 1)) {
      X[i, j] <- X[i-1, j] +
        ((1 / K_hat) * X[i-1, j-1] - (1 / K_hat) * X[i-1, j]) * dt[i] +
        sigma_hat * dW[j]
    }
  }
  
  #last reservoir 
  X[i, m] <- X[i-1, m] +
    ((1 / K_hat) * X[i-1, m-1]) * dt[i] +
    sigma_hat * dW[m]
  
  #observation (stormwater) from last state + obs noise
  Y_sim[i] <- X[i, m] + sigmae_hat * rnorm(1, 0, 1)
}

Y_sim[1] <- 0

#observed vs simulated stormwater
plot(t, data_ctsm$stormwater, type = "l",
     xlab = "Time [h]", ylab = "Stormwater",
     main = "Observed vs simulated stormwater (4-state SDE)")
lines(t, Y_sim, col = "blue", lwd = 2)
legend("topleft",
       legend = c("Observed", "Simulated"),
       col    = c("black", "blue"),
       lty    = 1, bty = "n")

#Residuals (obs - sim) over time
res_ts <- data_ctsm$stormwater - Y_sim

plot(t, res_ts, type = "l", col = "gray",
     xlab = "Time [h]", ylab = "Residual",
     main = "Residuals over time (observed - simulated)")
abline(h = 0, col = "red", lty = 2)

