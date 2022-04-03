# Metropolis algorithm
metropolis <- function(chain, theta_cur, sigma = 0.2, dfunc, space_min = -Inf, space_max = Inf, burnin = 0) {
  if (!is.function(dfunc)) stop("argument dfunc is not a function")
  theta_freq <- theta_cur # theta_freq <- NULL

  while (length(theta_freq) < chain) {
    delta_theta <- rnorm(1, mean = 0, sd = sigma)
    theta_pro <- theta_cur + delta_theta
    # theta_pro <- rnorm(1, mean = theta_cur, sd = sigma)

    if (theta_pro < space_min) {
      pmoving <- 0
      theta_pro <- space_min
    } else if (theta_pro > space_max) {
      pmoving <- 0
      theta_pro <- space_max
    } else if (dfunc(theta_cur) == 0) {
      pmoving <- 1
    } else {
      pmoving <- min(dfunc(theta_pro) / dfunc(theta_cur), 1)
    }

    if (runif(1) <= pmoving) {
      theta_freq <- append(theta_freq, theta_pro)
      theta_cur <- theta_pro
    } else {
      theta_freq <- append(theta_freq, theta_cur)
    }

  }

  return(tail(theta_freq, chain - burnin))
}

# Example
# posterior function 
# prior Beta(a = 1, b = 1)
# likelihood Bernoulli: N = 20, z = 14
posterior <- function(theta, a, b, N, z) {
  dbeta(theta, z + a, N - z + b)
}
a = 1
b = 1
N = 20 
z = 14
posterior_func <- function(theta) posterior(theta, a, b, N, z)

dposterior1 <- metropolis(chain = 5000, theta_cur = 0.1, dfunc = posterior_func, space_min = 0, space_max = 1, burnin = 2500)
dposterior2 <- metropolis(chain = 5000, theta_cur = 0.9, dfunc = posterior_func, space_min = 0, space_max = 1, burnin = 2500)
dposterior3 <- metropolis(chain = 5000, theta_cur = 0.5, dfunc = posterior_func, space_min = 0, space_max = 1, burnin = 2500)