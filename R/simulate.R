#' @title Simulate R1 and R2
#' @description Simulation function for: R2 the number of people, out of n, that have been declared positive
#' in the random sample, and among these ones, the ones that have been declared as infected before, i.e. R1.
#' @param pi2       A \code{numeric} that provides the true poportion of proportion of people in this population who are positive.
#' @param pi1       A \code{numeric} that provides the proportion of people (in the whole population) who are positive, as measured through a non-random,
#' but systematic sampling (e.g. based on medical selection).
#' @param n         A \code{numeric} that provides the sample size.
#' @param alpha1    A \code{numeric} that provides the False Negative (FN) rate for the sample R1. Default value is \code{0}.
#' @param alpha2    A \code{numeric} that provides the False Negative (FN) rate for the sample R2. Default value is \code{0}.
#' @param beta1     A \code{numeric} that provides the False Positive (FP) rate for the sample R1. Default value is \code{0}.
#' @param beta2     A \code{numeric} that provides the False Positive (FP) rate for the sample R2. Default value is \code{0}.
#' @param seed      A \code{numeric} that provides the simulation seed. Default value is \code{NULL}.
#' @param ...       Additional arguments.
#' @return A \code{list} containing R1, R2, pi1, n, alpha1, alpha2, beta1 and beta2.
#' @export
#' @author Stephane Guerrier
#' @examples
#' # Samples without measurement error
#' sim_Rs(pi2 = 3/100, pi1 = 1/100, n = 1500, seed = 18)
#'
#' # With measurement error
#' sim_Rs(pi2 = 3/100, pi1 = 1/100, n = 1500, alpha1 = 0.01,
#' alpha2 = 0.01, beta1 = 0.05, beta2 = 0.05, seed = 18)
sim_Rs = function(pi2, pi1, n, alpha1 = 0, alpha2 = 0, beta1 = 0, beta2 = 0, seed = NULL,...){
  # Simulation seed
  if (!is.null(seed)){
    set.seed(seed)
  }

  # Check inputs
  if (n %% 1 != 0){
    stop("The input n must be an integer.")
  }

  if (min(pi2, pi1) <= 0 || max(pi2, pi1) >= 1 || pi1 >= pi2){
    stop("The inputs pi1 and pi2 must be such that 0 < pi1 < pi2 < 1")
  }

  if (alpha1 < 0 || alpha1 > 0.5){
    stop("The input alpha1 must be 0 <= alpha1 < 0.5.")
  }

  if (alpha2 < 0 || alpha2 > 0.5){
    stop("The input alpha2 must be 0 <= alpha2 < 0.5.")
  }

  if (beta1 < 0 || beta1 > 0.5){
    stop("The input beta1 must be 0 <= beta1 < 0.5.")
  }

  if (beta2 < 0 || beta2 > 0.5){
    stop("The input beta2 must be 0 <= beta2 < 0.5.")
  }

  # Define lambda
  lambda = pi2/pi1

  # Draw R2
  R2 = rbinom(n = 1, size = n, prob = pi2)

  # Add potential false pos.
  if (alpha2 > 0){
    FP2 =  rbinom(n = 1, size = n - R2, prob = alpha2)
  }else{
    FP2 = 0
  }

  # Add potential false neg.
  if (beta2 > 0){
    FN2 = rbinom(n = 1, size = R2, prob = beta2)
  }else{
    FN2 = 0
  }

  # Adjust R2
  R2 = R2 + FP2 - FN2

  # Draw R1
  R1 = rbinom(n = 1, size = R2, prob = 1/lambda)

  # Add potential false pos.
  if (alpha1 > 0){
    FP1 = rbinom(n = 1, size = R2 - R1, prob = alpha1)
  }else{
    FP1 = 0
  }

  # Add potential false neg.
  if (beta1 > 0){
    FN1 = rbinom(n = 1, size = R1, prob = beta1)
  }else{
    FN1 = 0
  }

  # Adjust R1
  R1 = R1 + FP1 - FN1

  # Output
  list(R1 = R1, R2 = R2, pi1 = pi1, n = n, alpha1 = alpha1, alpha2 = alpha2,
       beta1 = beta1, beta2 = beta2, ...)
}



#' @title Simulate R1 and R2 (fast version)
#' @description Faster simulation function for: R2 the number of people, out of n, that have been declared positive
#' in the random sample, and among these ones, the ones that have been declared as infected before, i.e. R1.
#' @param pi2       A \code{numeric} that provides the true poportion of proportion of people in this population who are positive.
#' @param pi1       A \code{numeric} that provides the proportion of people (in the whole population) who are positive, as measured through a non-random,
#' but systematic sampling (e.g. based on medical selection).
#' @param n         A \code{numeric} that provides the sample size.
#' @param alpha1    A \code{numeric} that provides the False Negative (FN) rate for the sample R1. Default value is \code{0}.
#' @param alpha2    A \code{numeric} that provides the False Negative (FN) rate for the sample R2. Default value is \code{0}.
#' @param beta1     A \code{numeric} that provides the False Positive (FP) rate for the sample R1. Default value is \code{0}.
#' @param beta2     A \code{numeric} that provides the False Positive (FP) rate for the sample R2. Default value is \code{0}.
#' @param seed      A \code{numeric} that provides the simulation seed. Default value is \code{NULL}.
#' @param ...       Additional arguments.
#' @return A \code{list} containing R1 and R2
#' @author Stephane Guerrier
#' @export
#' @examples
#' # Samples without measurement error
#' sim_Rs_fast(pi2 = 3/100, pi1 = 1/100, n = 1500, seed = 18)
#'
#' # With measurement error
#' sim_Rs_fast(pi2 = 3/100, pi1 = 1/100, n = 1500, alpha1 = 0.01,
#' alpha2 = 0.01, beta1 = 0.05, beta2 = 0.05, seed = 18)
sim_Rs_fast = function(pi2, pi1, n, alpha1, alpha2, beta1, beta2, seed){
  # Set seed
  set.seed(seed)

  # Define lambda
  lambda = pi2/pi1

  # Draw R2
  R2 = rbinom(n = 1, size = n, prob = pi2)

  # Add potential false pos.
  if (alpha2 > 0){
    FP2 =  rbinom(n = 1, size = n - R2, prob = alpha2)
  }else{
    FP2 = 0
  }

  # Add potential false neg.
  if (beta2 > 0){
    FN2 = rbinom(n = 1, size = R2, prob = beta2)
  }else{
    FN2 = 0
  }

  # Adjust R2
  R2 = R2 + FP2 - FN2

  # Draw R1
  R1 = rbinom(n = 1, size = R2, prob = 1/lambda)

  # Add potential false pos.
  if (alpha1 > 0){
    FP1 = rbinom(n = 1, size = R2 - R1, prob = alpha1)
  }else{
    FP1 = 0
  }

  # Add potential false neg.
  if (beta1 > 0){
    FN1 = rbinom(n = 1, size = R1, prob = beta1)
  }else{
    FN1 = 0
  }

  # Adjust R1
  R1 = R1 + FP1 - FN1

  # Output
  list(R1 = R1, R2 = R2)
}
