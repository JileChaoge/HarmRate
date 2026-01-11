#' @title Estimating the Marginal Harm Rate
#'
#' @description Assume that \eqn{Y^1} and \eqn{Y^0} are potential binary
#'     outcomes for each individual. We can use \code{HarmRateRho} to estimate
#'     the harm rate, defined as \deqn{\text{HR} = \mathbb{P}(Y^1 = 0, Y^0 = 1),}
#'     which represents the proportion of individuals who would experience a worse
#'     outcome if they receive the treatment rather than the control. We estimate
#'     it given a marginal correlation coefficient \eqn{\rho} between \eqn{Y^1} and \eqn{Y^0}.
#'
#' @param A A vector of treatment.
#' @param Y A vector of outcome of interest.
#' @param mu0 The vector of outcome regression functions for potential outcome \eqn{Y^0}.
#' @param mu1 The vector of outcome regression functions for potential outcome \eqn{Y^1}.
#' @param ps The vector of propensity scores.
#' @param rho A scalar, denoting the correlation coefficient between \eqn{Y^0}
#'      and \eqn{Y^1}.
#'
#' @return A list which includes:
#' \itemize{
#'    \item \code{harm.rate.est}: the estimate of \eqn{\text{HR}}.
#'    \item \code{harm.rate.std}: the corresponding standard deviation.}
#' @export
#'
#'
#' @examples
#'
#' # generate simulated data
#' n <- 1000
#' X1 <- rnorm(n)
#' X2 <- rnorm(n)
#' A <-  rbinom(n, size = 1, prob = plogis(0.5*X1 - 0.5*X2))
#' U <-  rnorm(n)  # runif(n, -0.5, 0.5)
#' Y0 <- rbinom(n, size = 1, prob = plogis(0.5*X1 + 0.5*X2 + 3*U))
#' Y1 <- rbinom(n, size = 1, prob = plogis(0.5*X1 + 0.5*X2 + 1 + 1.5*U))
#' Y <- A*Y1 + (1-A)*Y0
#'
#' true.harm.rate <- sum((Y0 == 1) & (Y1 == 0)) / n
#'
#' # estimate the harm rate
#' X.design <- data.frame(cbind(1, X1,X2))
#' # propensity score
#' mod.ps <- glm(A ~ X1 + X2, family = binomial(link = 'logit'))
#' ps <- predict.glm(mod.ps, newdata = X.design, type = 'response')
#' # outcome regression functions
#' mod.mu1 <- glm(Y ~ X1 + X2, family = binomial(link = 'logit'),
#'                subset = A == 1)
#' mod.mu0 <- glm(Y ~ X1 + X2, family = binomial(link = 'logit'),
#'                subset = A == 0)
#' mu1 <- predict.glm(mod.mu1, newdata = X.design, type = 'response')
#' mu0 <- predict.glm(mod.mu0, newdata = X.design, type = 'response')
#'
#' rho.lst <- seq(0, 0.3, 0.05)
#' est.lst <- rep(NA, length(rho.lst))
#' std.lst <- rep(NA, length(rho.lst))
#'
#' for(i in 1:length(rho.lst)){
#'   rho <- rho.lst[i]
#'   mod <- HarmRateRho(A, Y, mu0, mu1, ps, rho = rho)
#'   est.lst[i] <-mod$harm.rate.est
#'   std.lst[i] <- mod$harm.rate.std
#' }
#'
#' est.lst
#' std.lst
#'
#'
HarmRateRho <- function(A, Y, mu0, mu1, ps, rho){

  n <- length(A)

  # step 1. construct the influence score
  phi.0 <- mu0*(1-mu1) + (1-A)*(Y-mu0)*(1-mu1)/(1-ps) - A*(Y-mu1)*mu0/ps
  phi.r <- (1-2*mu1)*sqrt(mu0*(1-mu0)/(mu1*(1-mu1)))*A*(Y-mu1)/(2*ps) +
    (1-2*mu0)*sqrt(mu1*(1-mu1)/(mu0*(1-mu0)))*(1-A)*(Y-mu0)/(2*(1-ps)) +
    sqrt(mu0*(1-mu0)*mu1*(1-mu1))

  g.eta <- mu0*(1-mu1) - rho * sqrt(mu0*(1-mu0)*mu1*(1-mu1))
  psi <- 1*(g.eta>=0) * (phi.0 - rho*phi.r)

  # step 2.point estimate and standard deviation of harm rate
  beta.hat <- mean(psi)
  std <- stats::sd(psi)/sqrt(n)

  return(list(harm.rate.est = beta.hat, harm.rate.std = std))
}

#' @references
#' Peng W., Peng D. Zhi G. and Yue L. (2025).
#' Quantifying Individual Risk for Binary Outcome,
#' \emph{arXiv:2402.10537}
