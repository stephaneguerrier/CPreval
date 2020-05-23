#' @title Compute proportion in the survey sample
#' @description Proportion estimated using the survey sample and confidence intervals based on the Clopper–Pearson approach.
#' @param R2        A \code{numeric} that provides the people of positive people in the sample.
#' @param n         A \code{numeric} that provides the sample size.
#' @param alpha2    A \code{numeric} that provides the False Negative (FN) rate for the sample R2. Default value is \code{0}.
#' @param beta2     A \code{numeric} that provides the False Positive (FP) rate for the sample R2. Default value is \code{0}.
#' @param gamma     A \code{numeric} that used to compute a (1 - gamma) confidence region for the proportion. Default value is \code{0.05}.
#' @param ...       Additional arguments.
#' @return A \code{list} with the structure:
#' \itemize{
#'  \item estimate: Estimated proportion
#'  \item sd: Estimated standard error of the estimator
#'  \item ci: Confidence interval
#'  \item gamma: Confidence level (i.e. 1 - gamma) for confidence interval
#' }
#' @export
#' @author Stephane Guerrier
#' @examples
#' # Samples without measurement error
#' X = sim_Rs(pi2 = 30/1000, pi1 = 10/1000, n = 1500, seed = 18)
#' survey_sample(X$R2, X$n)
#'
#' # With measurement error
#' X = sim_Rs(pi2 = 30/1000, pi1 = 10/1000, n = 1500, alpha1 = 0.01,
#' alpha2 = 0.01, beta1 = 0.05, beta2 = 0.05, seed = 18)
#' survey_sample(X$R2, X$n)
#' survey_sample(X$R2, X$n, alpha = 0.01, beta2 = 0.05)
survey_sample = function(R2, n, alpha2 = 0, beta2 = 0, gamma = 0.05, ...){
  # Compute survey proportion
  pi_bar = R2/n

  # Adjust for false positive/negative
  if (max(alpha2, beta2) > 0){
    pi_bar = pi_bar/(1 - alpha2 - beta2) - alpha2
  }

  # Estimated standard error
  sd = sqrt(R2/n*(1 - R2/n)/n)

  # Compute 1 - gamma confidence interval
  upper = (qbeta(p = 1 - gamma/2, R2 + 1, n - R2) - alpha2)/(1 - alpha2 - beta2)
  lower = (qbeta(p = gamma/2, R2, n - R2 + 1) - alpha2)/(1 - alpha2 - beta2)
  list(estimate = pi_bar, sd = sd, ci = c(lower, upper), gamma = gamma, ...)
}


mle_auxiliary = function(R1, R2, n, pi1){
  pi1*(n - R2)/(n - R1) + (R2 - R1)/(n - R1)
}

compute_bias = function(pi2, pi1, n, seed, B, alpha1, alpha2, beta1, beta2){
  inter = 0
  for (i in 1:B){
    X = sim_Rs_fast(pi2 = pi2, pi1 = pi1, n = n, alpha1 = alpha1, alpha2 = alpha2,
               beta1 = beta1, beta2 = beta2, seed = i + seed)
    inter = inter + mle_auxiliary(R1 = X$R1, R2 = X$R2, n = n, pi1 = pi1)
  }
  inter/B - pi2
}

#' @title Compute MLE based on R1 and R2
#' @description Proportion estimated using the MLE and confidence intervals based the asymptotic distribution of the estimator.
#' When measurement errors (false positive/negative) are considered, the estimator is corrected by simulations using the iterative
#' bootstrap.
#' @param R1        A \code{numeric} that provides the people of positive people in the sample that were known to be positive.
#' @param R2        A \code{numeric} that provides the people of positive people in the sample.
#' @param n         A \code{numeric} that provides the sample size.
#' @param alpha1    A \code{numeric} that provides the False Negative (FN) rate for the sample R1. Default value is \code{0}.
#' @param alpha2    A \code{numeric} that provides the False Negative (FN) rate for the sample R2. Default value is \code{0}.
#' @param beta1     A \code{numeric} that provides the False Positive (FP) rate for the sample R1. Default value is \code{0}.
#' @param beta2     A \code{numeric} that provides the False Positive (FP) rate for the sample R2. Default value is \code{0}.
#' @param gamma     A \code{numeric} that used to compute a (1 - gamma) confidence region for the proportion. Default value is \code{0.05}.
#' @param B         A \code{numeric} that corresponds to the number of simulations used in correction. This is only used in the case of
#' measurement error (i.e. \code{max(alpha1, alpha2, beta1, beta2) > 0}) in the simulation-based correction method.
#' Default value is \code{10^3}.
#' @param seed      A \code{numeric} that provides the simulation seed. This is only used in the case of measurement error (i.e.
#' \code{max(alpha1, alpha2, beta1, beta2) > 0}) in the simulation-based correction method. Default value is \code{18}.
#' @param eps       A \code{numeric} that corresponds to tolerance used to decide when to stop the simulation algorithm (iterative bootstrap).
#' This is only used in the case of measurement error (i.e. \code{max(alpha1, alpha2, beta1, beta2) > 0}) in the simulation-based correction method.
#' Default value is \code{10^(-6)}.
#' @param iter_max  A \code{numeric} that corresponds to the maximum of iteration of the simulation algorithm (iterative bootstrap).
#' This is only used in the case of measurement error (i.e. \code{max(alpha1, alpha2, beta1, beta2) > 0}) in the simulation-based correction method.
#' Default value is \code{100}.
#' @param ...       Additional arguments.
#' @return A \code{list} with the structure:
#' \itemize{
#'  \item estimate: Estimated proportion
#' }
#' @export
#' @author Stephane Guerrier
#' @examples
#' # Samples without measurement error
#' X = sim_Rs(pi2 = 3/100, pi1 = 1/100, n = 1500, seed = 18)
#' mle(R1 = X$R1, R2 = X$R2, pi1 = X$pi1, n = X$n)
#'
#' # With measurement error
#' X = sim_Rs(pi2 = 30/1000, pi1 = 10/1000, n = 1500, alpha1 = 0.01,
#' alpha2 = 0.01, beta1 = 0.05, beta2 = 0.05, seed = 18)
#' mle(R1 = X$R1, R2 = X$R2, pi1 = X$pi1, n = X$n)
#' mle(R1 = X$R1, R2 = X$R2, pi1 = X$pi1, n = X$n)
mle = function(R1, R2, pi1, n, alpha1 = 0, alpha2 = 0,
               beta1 = 0, beta2 = 0, B = 10^3, seed = 18,
               eps = 10^(-6), iter_max = 100, gamma = 0.05, ...){


  if (max(alpha1, alpha2, beta1, beta2) > 0){
    res = rep(NA, (iter_max + 1))
    res[1] = inter = mle_auxiliary(R1, R2, n, pi1)
    estimate = NULL
    for (i in 1:iter_max){
      if (is.na(res[i]) || (i > 1 && (res[i] < 0 || res[i] > 1))){
        estimate = res[i-1]
        break
      }
      res[i+1] = res[1] - compute_bias(pi2 = res[i], pi1 = pi1, n = n, seed = seed,
                                       B = B, alpha1 = alpha1, alpha2 = alpha2,
                                       beta1 = beta1, beta2 = beta2)

      if (i > 5){
        inter_old = inter
        inter = mean(res[(i-5):i])

        if (abs(inter_old - inter) < eps){
          estimate = res[i]
          break
        }
      }
    }

    if (is.null(estimate)){
      estimate = res[iter_max + 1]
    }

    # Compute asymptotic covariance matrix

    sd = NULL
    ci = NULL
  }else{
    estimate = mle_auxiliary(R1, R2, n, pi1)
    sd = sqrt(((estimate - pi1)*(1 - estimate))/(n*(1 - pi1)))
    ci = estimate + c(-1, 1)*qnorm(1 - gamma/2)*sd
  }

  list(estimate = estimate, sd = sd, ci = ci, gamma = gamma, ...)
}



#' @title Compute approximated MLE based on R1 and R2
#' @description Proportion estimated using the approximated MLE and confidence intervals based the asymptotic distribution of the estimator or
#' the Clopper–Pearson approach when alpha2 = 0.
#' @param R1        A \code{numeric} that provides the people of positive people in the sample that were known to be positive.
#' @param R2        A \code{numeric} that provides the people of positive people in the sample.
#' @param n         A \code{numeric} that provides the sample size.
#' @param alpha1    A \code{numeric} that provides the False Negative (FN) rate for the sample R1. Default value is \code{0}.
#' @param alpha2    A \code{numeric} that provides the False Negative (FN) rate for the sample R2. Default value is \code{0}.
#' @param beta1     A \code{numeric} that provides the False Positive (FP) rate for the sample R1. Default value is \code{0}.
#' @param beta2     A \code{numeric} that provides the False Positive (FP) rate for the sample R2. Default value is \code{0}.
#' @param gamma     A \code{numeric} that used to compute a (1 - gamma) confidence region for the proportion. Default value is \code{0.05}.
#' @param ...       Additional arguments.
#' @return A \code{list} with the structure:
#' \itemize{
#'  \item estimate: Estimated proportion
#' }
#' @export
#' @author Stephane Guerrier
#' @examples
#' # Samples without measurement error
#' X = sim_Rs(pi2 = 3/100, pi1 = 1/100, n = 1500, seed = 18)
#' modified_mle(R1 = X$R1, R2 = X$R2, pi1 = X$pi1, n = X$n)
#'
#' # With measurement error
#' X = sim_Rs(pi2 = 30/1000, pi1 = 10/1000, n = 1500, alpha1 = 0.01,
#' alpha2 = 0.01, beta1 = 0.05, beta2 = 0.05, seed = 18)
#' modified_mle(R1 = X$R1, R2 = X$R2, pi1 = X$pi1, n = X$n)
#' modified_mle(R1 = X$R1, R2 = X$R2, pi1 = X$pi1, n = X$n)
modified_mle = function(R1, R2, pi1, n, alpha1 = 0, alpha2 = 0, beta1 = 0, beta2 = 0, gamma = 0.05, ...){
  # Define Deltas
  Delta1 = 1 - alpha1 - beta1
  Delta2 = 1 - alpha2 - beta2

  # Define as
  a1 = 2*alpha2*Delta1*pi1
  a2 = (R2 - R1)/n + pi1*Delta1*Delta2 - alpha2*(1 - alpha2)
  a3 = 2*Delta2*(1 - alpha1)

  # Estimator
  estimate = (a2 + sqrt(a2^2 + a1*a3))/a3
  list(estimate = estimate)
}








