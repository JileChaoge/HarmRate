#' @title Calculating Fréchet-Hoeffding Bounds of Rho Conditional on Covariates
#'
#' @description \code{RhoFHBound} calculates the sharp Fréchet-Hoeffding bounds
#'     of the Pearson correlation coefficient for binary outcomes, defined as \deqn{\rho(x)=\text{Corr}(Y^0,Y^1\mid X=x).}
#'     These bounds are based only on unconfoundedness assumption (\eqn{Y^1 \perp \!\!\! \perp Y^0 \mid X}).
#'
#' @param mu0 The vector of outcome regression functions for potential outcome \eqn{Y^0}.
#' @param mu1 The vector of outcome regression functions for potential outcome \eqn{Y^1}.
#'
#' @return A list which includes:
#' \itemize{
#'    \item \code{rho.FH.low}: the estimated Fréchet-Hoeffding lower bound of \eqn{\rho(x)}.
#'    \item \code{rho.FH.up}: the estimated Fréchet-Hoeffding upper bound of \eqn{\rho(x)}.}
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
#' X.design <- data.frame(cbind(1, X1,X2))
#'
#' # outcome regression functions
#' mod.mu1 <- glm(Y ~ X1 + X2, family = binomial(link = 'logit'),
#'                subset = A == 1)
#' mod.mu0 <- glm(Y ~ X1 + X2, family = binomial(link = 'logit'),
#'                subset = A == 0)
#' mu1 <- predict.glm(mod.mu1, newdata = X.design, type = 'response')
#' mu0 <- predict.glm(mod.mu0, newdata = X.design, type = 'response')
#'
#' mod <- RhoFHBound(mu0, mu1)
#' rho.FH.low <- mod$rho.FH.low
#' rho.FH.up <- mod$rho.FH.up
#'
#' rho.FH.low
#' rho.FH.up
#'
#'
RhoFHBound <- function(mu0, mu1){

  L_rho <- max(-(1-mu0)*(1-mu1),  -mu0*mu1)/sqrt(mu0*(1-mu0)*mu1*(1-mu1))
  U_rho <- min( mu0*(1-mu1),  (1-mu0)*mu1)/sqrt(mu0*(1-mu0)*mu1*(1-mu1))

  return(list(rho.FH.low = L_rho, rho.FH.up = U_rho))
}

#' @references
#' Peng W., Peng D. Zhi G. and Yue L. (2025).
#' Quantifying Individual Risk for Binary Outcome,
#' \emph{arXiv:2402.10537}
