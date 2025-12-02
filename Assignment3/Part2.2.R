library("ctsmTMB")
library("lubridate")
library("tidyverse")
library("dplyr")
library("ggplot2")
library("patchwork")

data <- read.csv('/Users/sunjiawen/Desktop/Advanced ts/02427_CE03_Rainfall_runoff_exercise_data/ex2_overflow.csv')

head(data)

# Then we need to convert data format as the question demanding.
df <- data.frame(
  t = as.numeric(difftime(data$timestamp, min(data$timestamp), units = "hours")),
  U = as.numeric(data$rainfall),
  y = as.numeric(data$stormwater)
)

# Try to plot two graphs as Part2 shows
par(mfrow = c(2, 1), mar=c(4,4,2,1))
# Runoff Plot
plot(df$t, df$y, type="l", lwd=2,
     xlab="", ylab="Runoff", main="")

# Rainfall Plot
plot(df$t, df$U, type="l", lwd=1,
     xlab="Time", ylab="Rainfall", main="")

# 2.2.1 Linear Reservoir Model (n = 2)
# Linear reservoir model with 3 states (n = 2)
m <- ctsmTMB$new()


# STATE 1
m$addSystem(
  dX1 ~ A * U * dt - (2/K) * X1 * dt + sigma * dw
)

# STATE 2
m$addSystem(
  dX2 ~ (2/K) * X1 * dt - (2/K) * X2 * dt + sigma * dw
)

# STATE 3
m$addSystem(
  dX3 ~ (2/K) * X2 * dt + sigma * dw
)

# y observes X3
m$addObs(
  y ~ X3
)

# Add measurement noise variance (sig_e)
m$setVariance(
  y ~ sigma_y^2
)

# Add Input(Rainfall)
m$addInput(U)

# Set Initial States
m$setInitialState(list(c(df$U[1], df$U[1], df$y[1]), 1e-1 * diag(3)))


# Parameter initialization
m$setParameter(
  A = c(initial = 1, lower = 0, upper = 500),
  K = c(initial = 20, lower = 1e-3, upper = 500),
  sigma = c(initial = 0.1, lower = 1e-10, upper = 30),
  sigma_y = c(initial = 0.1, lower = 1e-10, upper = 30)
)

# Predict the model
fit_lin <- m$estimate(df, method = "ekf", compile = TRUE)
pred_lin <- m$predict(df, k.ahead = 10)

summary(fit_lin)


## Visualization
s10 <- pred_lin$states %>% filter(k.ahead == 10)

ggplot() +
  geom_line(aes(x = s10$t.j, y = s10$X3, color = "10-step Prediction"),
            size = 1.1) +
  
  geom_line(aes(x = df$t, y = df$y, color = "Observed"),
            size = 0.9) +
  
  scale_color_manual(values = c(
    "Observed" = "black",
    "10-step Prediction" = "red"
  )) +
  
  labs(title = "10-step Ahead Prediction vs Observed",
       x = "Time", y = "Stormwater",
       color = "") +
  
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")


# Residuals extraction
residuals <- fit_lin$residuals$residuals$y

ggplot(data.frame(t = df$t, residuals = residuals), aes(x = t, y = residuals)) +
  geom_line() +
  labs(title = "Residuals over Time", x = "Time", y = "Residuals", color = "blue") +
  theme_minimal()

# QQ plot
qqnorm(residuals)
qqline(residuals, col = "red")

# ACF (LACF)
acf(residuals)

# PACF (PLACF)
pacf(residuals)


# 2.2.2 Introduce overflow into the system
m2 <- ctsmTMB$new()

#----------------------------------
# STATE EQUATIONS
#----------------------------------

# STATE 1
m2$addSystem(
  dX1 ~ A * U * dt 
  - (2/K) * X1 * dt 
  + sigma * dw1
)

# STATE 2
m2$addSystem(
  dX2 ~ (2/K) * X1 * dt 
  - (2/K) * X2 * dt 
  + sigma * dw2
)

# STATE 3 with overflow
m2$addSystem(
  dX3 ~ (2/K) * (1 / (1 + exp(-alpha * (X2 - beta)))) * X2 * dt
  + sigma * dw3
)

#----------------------------------
# OBSERVATION EQUATIONS
#----------------------------------

m2$addObs(
  y ~ X3
)

# Measurement noise
m2$setVariance(
  y  ~ sigma_storm^2
)

#----------------------------------
# INPUT
#----------------------------------
m2$addInput(U)

#----------------------------------
# INITIAL STATE
#----------------------------------
m2$setInitialState(list(c(df$U[1], df$U[1], df$y[1]), 1e-1 * diag(3)))

#----------------------------------
# PARAMETERS
#----------------------------------
m2$setParameter(
  A = c(initial = 1, lower = 0, upper = 500),
  K = c(initial = 1, lower = 1e-3, upper = 500),
  
  alpha = c(initial = 1, lower = 0.01, upper = 1000),
  beta  = c(initial = 0.5, lower = 0.01, upper = 1000),
  
  sigma = c(initial = 0.1, lower = 1e-10, upper = 30),
  sigma_storm  = c(initial = 0.1, lower = 1e-10, upper = 30)
)

#----------------------------------
# FIT THE MODEL
#----------------------------------
fit <- m2$estimate(df, method = "ekf", compile=T)

# Predictions
pred <- m2$predict(df, k.ahead = 10)

summary(fit)

# Extract the 1-step ahead (prior) estimates\
s <- pred$states %>% filter(k.ahead == 10)

ggplot() +
  geom_ribbon(aes(x = s$t.j,
                  ymin = s$X3 - 2*sqrt(s$var.X3),
                  ymax = s$X3 + 2*sqrt(s$var.X3)),
              fill = "grey80", alpha = 0.6) +
  
  geom_line(aes(x = s$t.j, y = s$X3, color = "10-step Prediction"),
            size = 1.1) +
  
  geom_line(aes(x = df$t, y = df$y, color = "Observed"),
            size = 0.9) +
  
  scale_color_manual(values = c(
    "Observed" = "black",
    "10-step Prediction" = "red"
  )) +
  
  labs(title = "10-step Ahead Prediction vs Observed",
       x = "Time",
       y = "Stormwater",
       color = "") +
  
  theme_minimal() +
  theme(legend.position = "top")



residuals_2 <- fit$residuals$residuals$y

ggplot(data.frame(t = df$t, residuals = residuals_2), 
       aes(x = t, y = residuals_2)) +
  geom_line() +
  labs(title = "Residuals",
       x = "Time", y = "Residual") +
  theme_minimal()

# QQ plot
qqnorm(residuals_2)
qqline(residuals_2, col = "red")

# ACF (LACF)
acf(residuals_2)

# PACF (PLACF)
pacf(residuals_2)


# report the log-likelihood and all parameters


logLik_lin <- -fit_lin$nll
logLik_lin

logLik_sig <- -fit$nll
logLik_sig



## 2.2.3
# Report the correlation matrix for the parameters.

cor_matlin <- cov2cor(fit_lin$cov.fixed)
cor_matlin
cor_matlin2 <- cov2cor(fit$cov.fixed)
cor_matlin2

fit_lin$par.fixed
fit$par.fixed



