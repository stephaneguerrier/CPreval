#' @title Simulate R and R0
#' @description Simulation function for: R the number of people, out of n, that have been declared positive
#' in the random sample, and among these ones, the ones that have been declared as infected before, i.e. R0.
#' @param p         A \code{numeric} that provides the true poportion of proportion of people in this population who are positive.
#' @param pi0       A \code{numeric} that provides the proportion of people (in the whole population) who are positive, as measured through a non-random,
#' but systematic sampling (e.g. based on medical selection).
#' @param n         A \code{numeric} that provides the sample size.
#' @param alpha0    A \code{numeric} that provides the False Negative (FN) rate for the sample R0. Default value is \code{0}.
#' @param alpha0    A \code{numeric} that provides the False Negative (FN) rate for the sample R0. Default value is \code{0}.
#' @param beta      A \code{numeric} that provides the False Positive (FP) rate for the sample R. Default value is \code{0}.
#' @param beta      A \code{numeric} that provides the False Positive (FP) rate for the sample R. Default value is \code{0}.
#' @param seed      A \code{numeric} that provides the simulation seed. Default value is \code{NULL}.
#' @param ...       Additional arguments.
#' @return A \code{list} containing R, R0, pi0, n, alpha, alpha0, beta and beta0.
#' @export
#' @author Stephane Guerrier
#' @examples
#' # Samples without measurement error
#' sim_Rs(p = 3/100, pi0 = 1/100, n = 1500, seed = 18)
#'
#' # With measurement error
#' sim_Rs(p = 3/100, pi0 = 1/100, n = 1500, alpha0 = 0.01,
#' alpha = 0.01, beta0 = 0.05, beta = 0.05, seed = 18)
sim_Rs = function(p, pi0, n, alpha0 = 0, alpha = 0, beta0 = 0, beta = 0, seed = NULL,...){
  # Simulation seed
  if (!is.null(seed)){
    set.seed(seed)
  }

  # Check inputs
  if (n %% 1 != 0){
    stop("The input n must be an integer.")
  }

  if (min(p, pi0) <= 0 || max(p, pi0) >= 1 || pi0 >= p){
    stop("The inputs pi0 and p must be such that 0 < pi0 < p < 1")
  }

  if (alpha0 < 0 || alpha0 > 0.5){
    stop("The input alpha0 must be 0 <= alpha0 < 0.5.")
  }

  if (alpha < 0 || alpha > 0.5){
    stop("The input alpha must be 0 <= alpha < 0.5.")
  }

  if (beta0 < 0 || beta0 > 0.5){
    stop("The input beta0 must be 0 <= beta0 < 0.5.")
  }

  if (beta < 0 || beta > 0.5){
    stop("The input beta must be 0 <= beta < 0.5.")
  }

  # Define lambda
  lambda = p/pi0

  # Draw R2
  R2 = rbinom(n = 1, size = n, prob = p)

  # Add potential false pos.
  if (alpha > 0){
    FP2 =  rbinom(n = 1, size = n - R2, prob = alpha)
  }else{
    FP2 = 0
  }

  # Add potential false neg.
  if (beta > 0){
    FN2 = rbinom(n = 1, size = R2, prob = beta)
  }else{
    FN2 = 0
  }

  # Adjust R2
  R2 = R2 + FP2 - FN2

  # Draw R1
  R1 = rbinom(n = 1, size = R2, prob = 1/lambda)

  # Add potential false pos.
  if (alpha0 > 0){
    FP1 = rbinom(n = 1, size = R2 - R1, prob = alpha0)
  }else{
    FP1 = 0
  }

  # Add potential false neg.
  if (beta0 > 0){
    FN1 = rbinom(n = 1, size = R1, prob = beta0)
  }else{
    FN1 = 0
  }

  # Adjust R1
  R1 = R1 + FP1 - FN1

  # Output
  list(R0 = R1, R = R2, pi0 = pi0, n = n, alpha0 = alpha0, alpha = alpha,
             beta0 = beta0, beta = beta, ...)
}

