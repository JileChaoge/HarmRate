#' @title Calculating Fréchet-Hoeffding Bounds of Harm Rate Conditional on Covariates
#'
#' @description \code{HarmRateFHBound} calculates the sharp Fréchet-Hoeffding
#'     bounds of the harm rate conditioanl on covariates, defined as \deqn{\text{HR(x)} = \mathbb{P}(Y^1 = 0, Y^0 = 1\mid X=x).} These bounds are based only on the unconfoundedness assumption (\eqn{Y^1 \perp \!\!\! \perp Y^0 \mid X}).
#'
#' @param mu0 The vector of outcome regression functions for potential outcome \eqn{Y^0}.
#' @param mu1 The vector of outcome regression functions for potential outcome \eqn{Y^1}.
#'
#' @return A list which includes:
#' \itemize{
#'    \item \code{harm.rate.FH.low}: the estimated Fréchet-Hoeffding lower bound of \eqn{\text{HR}(x)}.
#'    \item \code{harm.rate.FH.up}: the estimated Fréchet-Hoeffding upper bound of \eqn{\text{HR}(x)}.}
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
#' true.marginal.harm.rate <- sum((Y0 == 1) & (Y1 == 0)) / n
#'
#' # estimate the harm rate
#' X.design <- data.frame(cbind(1, X1,X2))
#' # outcome regression functions
#' mod.mu1 <- glm(Y ~ X1 + X2, family = binomial(link = 'logit'),
#'                subset = A == 1)
#' mod.mu0 <- glm(Y ~ X1 + X2, family = binomial(link = 'logit'),
#'                subset = A == 0)
#' mu1 <- predict.glm(mod.mu1, newdata = X.design, type = 'response')
#' mu0 <- predict.glm(mod.mu0, newdata = X.design, type = 'response')
#'
#' mod <- HarmRateFHBound(mu0, mu1)
#' harm.rate.FH.low <- mod$harm.rate.FH.low
#' harm.rate.FH.up <- mod$harm.rate.FH.up
#'
#' harm.rate.FH.low
#' harm.rate.FH.up
#' mean(harm.rate.FH.low)
#' mean(harm.rate.FH.up)
#'
#'
HarmRateFHBound <- function(mu0, mu1){

  FH_low <- pmax(mu0-mu1, 0)
  FH_up <- pmin(mu0, 1-mu1)

  return(list(harm.rate.FH.low = FH_low, harm.rate.FH.up = FH_up))
}

#' @references
#' Peng W., Peng D. Zhi G. and Yue L. (2025).
#' Quantifying Individual Risk for Binary Outcome,
#' \emph{arXiv:2402.10537}
