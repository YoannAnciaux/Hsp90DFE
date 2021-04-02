#### FIT: DFE FGM ####
#' DFE of FGM from Martin & Lenormand (2015).
#'
#' Density of the distribution of fitness effects of random mutations from \href{https://doi.org/10.1111/evo.12671}{Martin
#' & Lenormand (2015)}. The density corresponds to equation (2) in the article:
#' \deqn{f_s(s, n, \lambda, s_o) = \frac{2}{\lambda}%
#' f_{\chi_{n}^2} (\frac{2 (s_o-s)}{\lambda}, \frac{2 s_o}{\lambda})}{%
#' fs(s, n, \lambda, so) = 2/\lambda d\chi[2,n](2 (so - s) / \lambda, 2 so / \lambda)}
#' The density is encoded using the functions from the package \code{\link[distr]{distr}}
#'
#' @param s Vector of real numbers. Quantiles of the density. The density is 0
#' if \eqn{s \ge so}.
#' @param n Natural number. The dimensionality, i.e. the number of phenotypic
#' dimensions (traits) under selection.
#' @param lambda Positive real number. Mutational variance per trait among random
#' mutations scaled by the strength of selection. It scales the mean fitness effect
#' of mutations.
#' @param so Positive real number. Fitness distance between the wild-type (phenotype
#' without mutations) and the optimum.
#'
#' @return The density for each quantiles (selection coefficients) in a vector
#' format. The length of the vector is equal to the length of \code{s}.
#'
#' @examples
#' s <- seq(-2, 1, 0.1)
#' fs <- dfgm_martin(s, 2, 0.1, 1)
#'
#' @source Martin, G., & Lenormand, T. (2015). The fitness effect of mutations
#' across environments: Fisher's geometrical model with multiple optima.
#' Evolution, 69(6), 1433-1447.
#'
#'@export
dfgm_martin <- function(s, n, lambda, so) {
  C <- distr::Chisq(df = n, ncp = 2 * so / lambda)
  S <- so - lambda / 2 * C
  return(distr::d(S)(s))
}

#' DFE of FGM from Tenaillon (2014).
#'
#' Density of the distribution of fitness effects of random mutations from
#' \href{https://doi.org/10.1146/annurev-ecolsys-120213-091846}{Tenaillon (2014)}
#' adapted to show the same parameters as Martin & Lenormand (2015) plus the
#' parameters \code{Q} and \code{alpha}.
#' The density corresponds to the first equation in the paragraph *"Distribution of
#' mutation effects"* in the article.
#' The density is encoded using the functions from the package \code{\link[gsl]{gsl}}
#'
#' @inheritParams dfgm_martin
#' @param alpha Positive real number. Scaling parameter in the formulation of
#' Tenaillon et al. (2007). (Identifiability problem with \code{lambda})
#' @param Q Positive real number. Parameter that modifies the concavity of the
#' fitness decline from the optimum.
#'
#' @return The density for each quantiles (selection coefficients) in a vector
#' format. The length of the vector is equal to the length of \code{s}.
#'
#' @examples
#' s <- sample(seq(-2, 1.5, 0.1))
#' fs <- dfgm_tenaillon(s, 2, 0.1, 1, 1/2, 2)
#'
#' @source
#' \itemize{
#'   \item Tenaillon, O. (2014). The utility of Fisher's geometric model in
#'   evolutionary genetics. Annual review of ecology, evolution, and systematics,
#'   45, 179-201.
#'   \item Martin, G., & Lenormand, T. (2015). The fitness effect of mutations
#'   across environments: Fisher's geometrical model with multiple optima.
#'   Evolution, 69(6), 1433-1447.
#'   \item Tenaillon, O., Silander, O. K., Uzan, J. P., & Chao, L. (2007).
#'   Quantifying organismal complexity using a population genetic approach.
#'   PloS one, 2(2), e217.
#' }
#'
#'@export
dfgm_tenaillon <- function(s, n, lambda, so, alpha, Q) {
  # The density is 0 is s > so, however the function gsl::bessel_In produces
  # NA/Nan/Inf when used with s > so. To overcome this numerical issue we compute
  # the density only for the values of s <= so and set the density to 0 for the
  # other.S
  d <- numeric(length(s))
  idx_s_infeqso <- which(s <= so)
  s <- s[idx_s_infeqso]

  # The formula from Tenaillon (2014) is separated in the for multiplicative terms
  # A B C D for readibility.
  A <- exp( -(((-s + so) / alpha)^(2/Q) + (so / alpha)^(2/Q)) / (2 * lambda) ) / (alpha * Q * lambda)
  B <- ((-s + so) / alpha)^((n/2 + 1 - Q)/Q)
  C <- (so / alpha)^((1 - n/2) / Q)
  D <- gsl::bessel_In(n = n/2 - 1, x = 1/lambda * ((-s + so) / alpha)^(1/Q) * (so / alpha)^(1/Q))
  d_s_infeqso <- A*B*C*D
  d[idx_s_infeqso] <- d_s_infeqso
  return(d)
}

#### FIT: MLE ####
#' MLE on the DFE for FGM
#'
#' Maximum likelihood estimation (mle) for the parameters of the DFE density for
#' different form of the FGM.
#'
#' @inheritParams dfgm_martin
#' @param sd_error Vector of real numbers. Standard deviation for the normal
#' distribution of the measurement error for each value of \code{s}. It can be
#' either NULL (default) in which case the loglikelihood does not account for
#' measurement error. Or it can be a vector of real number the same size as \code{s}.
#' @param model String. "Martin" performs the mle on the density of the DFE from
#' Martin & Lenormand (2015) (see \code{\link{dfgm_martin}}). "Tenaillon" performs
#' the mle on the density of the DFE from Tenaillon (2014) (see \code{\link{dfgm_tenaillon}}).
#' Default "Martin".
#' @param start Numeric vector (named). Used as starting values in \code{\link{maxLik}{maxLik}}.
#' It must contains the parameters for the chosen \code{model}:
#' \itemize{
#'   \item{"Martin"}{ c(n = #, lambda = #, so = #)}
#'   \item{"Tenaillon"}{ c(n = #, lambda = #, so = #, alpha = #, Q = #)}
#' }
#' @param method String. Maximisation method. See \code{\link{maxLik}{maxLik}}.
#' Default "NM".
#' @param ... See \code{\link{maxLik}{maxLik}} all the possible parameters
#'
#' @return Returns the output of \code{\link[maxLik]{maxLik}} for the chosen model.
#'
#' @seealso \code{\link{dfgm_martin}}, \code{\link{dfgm_tenaillon}},
#' \code{\link{maxLik}{maxLik}}.
#'
#' @examples
#' #### model : "Martin" ####
#' # Parameters
#' nsample <- 1000; n <- 4; lambda <- 0.5; so <- 1; alpha <- 1/2; Q <- 2;
#' # Simulated data
#' s <- rfgm(nsample, n, lambda, so, alpha, Q)
#' # Noise on initial parameters
#' initial_parameter <- c(n = ceiling(abs(rnorm(1, n))),
#'                        lambda = abs(rnorm(1, lambda, sd = 0.1)),
#'                        so = abs(rnorm(1, n)))
#' # Constraints on parameters (not applioed for all methods see maxLik in maxLik package)
#' consA <- rbind(c(1, 0, 0),
#'                c(0, 1, 0),
#'                c(0, 0, 1))
#' consB <- c(0, 0, -max(s))
#' # MLE
#' res <- mle_dfgm(s, model = "Martin", start = initial_parameter,
#'                 method = "NM", constraints = list(ineqA = consA, ineqB = consB))
#' # MLE with measurement error
#' s_with_error <- s + rnorm(length(s), mean = 0, sd = 0.1)
#' res <- mle_dfgm(s_with_error, sd_error = rep(0.1, length(s)), model = "Martin", start = initial_parameter,
#'                 method = "NM", constraints = list(ineqA = consA, ineqB = consB))
#'
#' #### model : "Tenaillon" ####
#' # Parameters
#' nsample <- 1000; n <- 4; lambda <- 0.5; so <- 1; alpha <- 1/2; Q <- 2;
#' # Simulated data
#' s <- rfgm(nsample, n, lambda, so, alpha, Q)
#' # Noise on initial parameters
#' initial_parameter <- c(n = ceiling(abs(rnorm(1, n))),
#'                        lambda = abs(rnorm(1, lambda, sd = 0.1)),
#'                        so = abs(rnorm(1, n)),
#'                        alpha = alpha,
#'                        Q = abs(rnorm(1, Q)))
#' # Constraints on parameters (not applioed for all methods see maxLik in maxLik package)
#' consA <- rbind(c(1, 0, 0, 0, 0),
#'                c(0, 1, 0, 0, 0),
#'                c(0, 0, 1, 0, 0),
#'                c(0, 0, 0, 1, 0),
#'                c(0, 0, 0, 0, 1))
#' consB <- c(0, 0, -max(s), 0, 0)
#' # MLE
#' res <- mle_dfgm(s, model = "Tenaillon", start = initial_parameter,
#'                 method = "NM", constraints = list(ineqA = consA, ineqB = consB))
#'
#' @source
#' \itemize{
#'   \item Tenaillon, O. (2014). The utility of Fisher's geometric model in
#'   evolutionary genetics. Annual review of ecology, evolution, and systematics,
#'   45, 179-201.
#'   \item Martin, G., & Lenormand, T. (2015). The fitness effect of mutations
#'   across environments: Fisher's geometrical model with multiple optima.
#'   Evolution, 69(6), 1433-1447.
#'  }
#'
#' @export
mle_dfgm <- function(s, sd_error = NULL, model = "Martin", start, method = "NM", ...){

  coll <- checkmate::makeAssertCollection()
  checkmate::assert_numeric(s, finite = T, any.missing = F, null.ok = F, add = coll)
  checkmate::assert_numeric(sd_error, finite = T, any.missing = F, null.ok = T, add = coll)
  checkmate::assert_subset(model, choices = c("Martin", "Tenaillon"), empty.ok = F,
                           add = coll)
  checkmate::assert_vector(start, strict = T, any.missing = F, min.len = 3, max.len = 5,
                           unique = F, names = "named", null.ok = F, add = coll)
  checkmate::reportAssertions(coll)

  # probabilities for the dfe with error
  p1 <- pnorm(1/2)-pnorm(-1/2)
  p2 <- pnorm(3/2)-pnorm(1/2)
  p3 <- pnorm(5/2)-pnorm(3/2)

  # Log-Likelihood for the model "Tenaillon"
  if (model == "Martin") {
    checkmate::assert_names(names(start), permutation.of = c("n", "lambda", "so"))
    start <- start[c("n", "lambda", "so")]

    if (is.null(sd_error)) { #loglikelihood without normal measurement error
      loglik <- function(param) {
        n <- param[1]
        lambda <- param[2]
        so <- param[3]
        if(n > 0 & lambda > 0 & so > 0) {
          ll <- log(dfgm_martin(s, n, lambda, so))
        } else {
          ll <- rep(NA, length(s))
        }
        return(ll)
      }
    } else { #loglikelihood with normal measurement error
      loglik <- function(param) {
        n <- param[1]
        lambda <- param[2]
        so <- param[3]
        if(n > 0 & lambda > 0 & so > 0) {
          ll <- log(p1 / (p1+2*p2+2*p3) * dfgm_martin(s, n, lambda, so) +
                      p2 / (p1+2*p2+2*p3) * (dfgm_martin(s - sd_error, n, lambda, so) + dfgm_martin(s + sd_error, n, lambda, so)) +
                      p3 / (p1+2*p2+2*p3) * (dfgm_martin(s - 2 * sd_error, n, lambda, so) + dfgm_martin(s + 2 * sd_error, n, lambda, so)))
        } else {
          ll <- rep(NA, length(s))
        }
        return(ll)
      }
    }
  }
  # Log-Likelihood for the model "Tenaillon"
  if (model == "Tenaillon") {
    checkmate::assert_names(names(start), permutation.of = c("n", "lambda", "so", "alpha", "Q"))
    start <- start[c("n", "lambda", "so", "alpha", "Q")]

    if (is.null(sd_error)) { #loglikelihood without normal measurement error
      loglik <- function(param) {
        n <- param[1]
        lambda <- param[2]
        so <- param[3]
        alpha <- param[4]
        Q <- param[5]
        if(n > 0 & lambda > 0 & so > 0 & alpha > 0 & Q > 0) {
          ll <- log(dfgm_tenaillon(s, n, lambda, so, alpha, Q))
        } else {
          ll <- rep(NA, length(s))
        }
        return(ll)
      }
    } else { #loglikelihood with normal measurement error
      loglik <- function(param) {
        n <- param[1]
        lambda <- param[2]
        so <- param[3]
        alpha <- param[4]
        Q <- param[5]
        if(n > 0 & lambda > 0 & so > 0 & alpha > 0 & Q > 0) {
          ll <- log(p1 * dfgm_martin(s, n, lambda, so) +
                      p2 * (dfgm_martin(s - sd_error, n, lambda, so) + dfgm_martin(s + sd_error, n, lambda, so)) +
                      p3 * (dfgm_martin(s - 2 * sd_error, n, lambda, so) + dfgm_martin(s + 2 * sd_error, n, lambda, so)))
        } else {
          ll <- rep(NA, length(s))
        }
        return(ll)
      }
    }
  }

  # MLE
  tryCatch(maxLik::maxLik(loglik, start = start, method = method, ...),
           error = function(e) {
             est_tmp <- rep(NA, length(start))
             names(est_tmp) <- names(start)
             list(estimate = est_tmp,
                  maximum = NA,
                  hessian = matrix(NA,
                                   ncol = length(start),
                                   nrow = length(start),
                                   dimnames = list(names(start),
                                                   names(start))),
                  code = -1,
                  message = e)
           })
}

#### FIT: Wrapper fit dfe ####

fit_DFE <- function(s, sd = NULL, model = "Martin", method = "NM", round_n = T,
                    bounds = NULL, iterlim = 500, nmax = ncol(s)-1) {

  coll <- checkmate::makeAssertCollection()
  checkmate::assert_vector(s, any.missing = T, all.missing = F, add = coll)
  checkmate::assert_vector(sd, any.missing = T, all.missing = F, len = length(s),
                           null.ok = T, add = coll)


  # Remove NA in s does not test for NA in sd
  idx_noNA <- which(!is.na(s))
  if(length(idx_noNA) != length(s)) warning("NAs have been removed")
  s <- s[idx_noNA]
  if(!is.null(sd)) {
    sd <- sd[idx_noNA]
  }

  # Constraints for the MLE
  bounds_default <- list(n      = c(0.9, Inf),
                         so     = c(max(sqrt(.Machine$double.eps), max(s)), Inf),
                         lambda = c(sqrt(.Machine$double.eps), Inf),
                         alpha  = c(sqrt(.Machine$double.eps), Inf),
                         Q      = c(sqrt(.Machine$double.eps), Inf))
  if(!is.null(bounds)) {
    bounds_default[names(bounds)] <- bounds
  }
  bounds <- bounds_default

  consA = c()
  consB = c()
  if(bounds$n[1] != -Inf) {
    consA = rbind(consA, c(1, 0, 0))
    consB = c(consB, -bounds$n[1])
  }
  if(bounds$n[2] != Inf) {
    consA = rbind(consA, c(-1, 0, 0))
    consB = c(consB, bounds$n[2])
  }
  if(bounds$lambda[1] != -Inf) {
    consA = rbind(consA, c(0, 1, 0))
    consB = c(consB, -bounds$lambda[1])
  }
  if(bounds$lambda[2] != Inf) {
    consA = rbind(consA, c(0, -1, 0))
    consB = c(consB, bounds$lambda[2])
  }
  if(bounds$so[1] != -Inf) {
    consA = rbind(consA, c(0, 0, 1))
    consB = c(consB, -bounds$so[1])
  }
  if(bounds$so[2] != Inf) {
    consA = rbind(consA, c(0, 0, -1))
    consB = c(consB, bounds$so[2])
  }
  if(model == "Tenaillon") {
    consA = cbind(consA, numeric(nrow(consA)), numeric(nrow(consA)))
    if(bounds$alpha[1] != -Inf) {
      consA = rbind(consA, c(0, 0, 0, 1, 0))
      consB = c(consB, -bounds$alpha[1])
    }
    if(bounds$alpha[2] != Inf) {
      consA = rbind(consA, c(0, 0, 0, -1, 0))
      consB = c(consB, bounds$alpha[2])
    }
    if(bounds$Q[1] != -Inf) {
      consA = rbind(consA, c(0, 0, 0, 0, 1))
      consB = c(consB, -bounds$Q[1])
    }
    if(bounds$Q[2] != Inf) {
      consA = rbind(consA, c(0, 0, 0, 0, -1))
      consB = c(consB, bounds$Q[2])
    }
  }


  ml <- mle_dfgm(s, sd_error = sd, model = model,
                 start = c(n = nmax, lambda = abs(mean(s)), so = max(2 * c(.Machine$double.eps, max(s)))),
                 method = method, constraints = list(ineqA = consA, ineqB = consB),
                 control = list(iterlim = iterlim))
  par_AIC <- AIC(ml)
  if(round_n) {# round n and do another MLE with n fixed.
    if(!anyNA(ml$estimate)) {
      ml$estimate["n"] <- round(ml$estimate["n"])
      ml <- mle_dfgm(s = s, sd_error = sd, model = model, start = ml$estimate,
                     method = "NM", constraints = list(ineqA = consA, ineqB = consB),
                     fixed = c("n"), control = list(iterlim = iterlim))
      # the AIC computed here is not exactly the correct AIC as the likelihood is
      # the one with n fixed but we do not consider n fixed in the number of parameters
      # because we already did an optimization to find n
      par_AIC <- 2 * (length(ml$estimate) - ml$maximum)
      attr(par_AIC, "df") <- length(ml$estimate)
    }
  }

  list(estimate = ml$estimate,
       AIC = par_AIC,
       ml = ml)
}

# nb_mut = 1e3
# n = c(A = 2, B = 2, C = 1, D = 2)
# so = c(A = 1, B = 1.1, C = 0.1, D = 5)
# lambda = c(A = 0.1, B = 0.01, C = 0.05, D = 0.2)
# alpha = c(A = 1/2, B = 1/2, C = 1/2, D = 1/2)
# Q = c(A = 2, B = 2, C = 2, D = 2)
# angular_coord = "RAND"
#
# sim = simulate_mutant(nb_mut = nb_mut, n = n, so = so, lambda = lambda, alpha = alpha, Q = Q,
#                       nmax = length(n)-2)
#
# sd = NULL
# model = "Martin"
# method = "NM"
# round_n = T
# bounds = list(n = c(0.9, ncol(sim$s) - 2 + 0.1))
# iterlim = 5000
#
# fit = sapply(X = split(sim$s, col(sim$s, as.factor = T)),
#              simplify = F,
#              USE.NAMES = T,
#              FUN = fit_DFE,
#              sd = sd, model = model, method = method, round_n = round_n, bounds = bounds, iterlim = iterlim)
