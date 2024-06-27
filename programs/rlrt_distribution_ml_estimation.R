# density function for 0.5 * chi^2(1) + 0.5 * chi^2(2)
drlrt <- function(x) {
  0.5 * dchisq(x, df = 1) + 0.5 * dchisq(x, df = 2)
}

# likelihood function for sum of chisq distribution
log_likelihood <- function(params, data) {
  df1 <- params[1]
  df2 <- params[2]
  -sum(log(0.5 * dchisq(data, df = df1) + 0.5 * dchisq(data, df = df2)))
}

# ml estimation for sum of chisq distribution
estimate_parameters <- function(data) {
  initial_params <- c(1, 2) # initial guesses
  result <- optim(
    par = initial_params,
    fn = log_likelihood,
    data = data,
    method = "L-BFGS-B",
    lower = c(0.1, 0.1), # inf limits for degrees of freedom
    upper = c(Inf, Inf)  # sup limits for degrees of freedom
  )
  return(result$par)
}


# simulate data
set.seed(123)
data <- 0.5 * rchisq(100, df = 1) + 0.5 * rchisq(100, df = 2)

# estimate parameters
estimated_params <- estimate_parameters(data)
print(estimated_params)

# compare densities
dev.new()
hist(data, probability = TRUE, breaks = 30, col = "grey", main = "Compared distributions")
curve(drlrt(x), add = TRUE, col = "blue", lwd = 2)
curve(
  0.5 * dchisq(x,
    df = estimated_params[1]
  ) + 0.5 * dchisq(x, df = estimated_params[2]),
  add = TRUE, col = "red", lwd = 2
)
legend("topright",
  legend = c("Initial Density", "Estimated Density"),
  col = c("blue", "red"), lwd = 2
)
