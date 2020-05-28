

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

pi_of_theta = function(pi2, pi1, n, seed, B, alpha1, alpha2, beta1, beta2){
  inter = 0
  for (i in 1:B){
    X = sim_Rs_fast(pi2 = pi2, pi1 = pi1, n = n, alpha1 = alpha1, alpha2 = alpha2,
                    beta1 = beta1, beta2 = beta2, seed = i + seed)
    inter = inter + mle_auxiliary(R1 = X$R1, R2 = X$R2, n = n, pi1 = pi1)
  }
  inter/B
}

neg_log_like = function(pi2, R2, R1, pi1, n, alpha1, alpha2, beta1, beta2){
  # Define Deltas
  Delta1 = 1 - alpha1 - beta1
  Delta2 = 1 - alpha2 - beta2

  # Define modified probs
  pi2_star = pi2*Delta2 + alpha2
  pi_star = pi1/pi2*Delta1 + alpha1

  # Log likelihood times (-1)
  -(R2*log(pi2_star) + (n - R2)*log(1 - pi2_star) + R1*log(pi_star) + (R2 - R1)*log(1 - pi_star))
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
#' mle(R1 = X$R1, R2 = X$R2, pi1 = X$pi1, n = X$n, alpha1 = 0.01,
#' alpha2 = 0.01, beta1 = 0.05, beta2 = 0.05)
mle = function(R1, R2, pi1, n, alpha1 = 0, alpha2 = 0,
               beta1 = 0, beta2 = 0, B = 100, seed = 18,
               eps = 10^(-6), iter_max = 15, gamma = 0.05, ...){


  if (max(alpha1, alpha2, beta1, beta2) > 0){
    # res = rep(NA, (iter_max + 1))
    # res[1] = mle_auxiliary(R1, R2, n, pi1)
    # estimate = NULL
    # for (i in 1:iter_max){
    #   if (is.na(res[i]) || (i > 1 && (res[i] < 0 || res[i] > 1))){
    #     estimate = res[i-1]
    #     break
    #   }
    #   res[i+1] = res[1] - compute_bias(pi2 = res[i], pi1 = pi1, n = n, seed = seed,
    #                                    B = B, alpha1 = alpha1, alpha2 = alpha2,
    #                                    beta1 = beta1, beta2 = beta2)
    #
    #   if (i > 4){
    #     if (min(abs(res[i] - res[1:(i-1)])) < eps){
    #       estimate = res[i]
    #       break
    #     }
    #   }
    # }
    #
    # if (is.null(estimate)){
    #   estimate = res[iter_max + 1]
    # }
    #
    # # Compute asymptotic covariance matrix
    # inter = rep(NA, 10*B)
    # for (i in 1:(10*B)){
    #   X = sim_Rs_fast(pi2 = estimate, pi1 = pi1, n = n, alpha1 = alpha1, alpha2 = alpha2,
    #                   beta1 = beta1, beta2 = beta2, seed = i + seed)
    #   inter[i] = mle_auxiliary(R1 = X$R1, R2 = X$R2, n = n, pi1 = pi1)
    # }
    # V = var(inter)
    # d1 = pi_of_theta(pi2 = estimate, pi1 = pi1, n = n, alpha1 = alpha1, alpha2 = alpha2,
    #                  beta1 = beta1, beta2 = beta2, seed = i + seed, B = 10*B)
    # d2 = pi_of_theta(pi2 = estimate + 10^(-4), pi1 = pi1, n = n, alpha1 = alpha1, alpha2 = alpha2,
    #                  beta1 = beta1, beta2 = beta2, seed = i + seed, B = 10*B)
    # B = (d2 - d1)/10^(-4)
    #
    # sd = as.numeric(sqrt(V)/B)
    # lower = estimate - qnorm(1 - gamma/2)*sd
    # upper = estimate + qnorm(1 - gamma/2)*sd
    # ci = c(lower, upper)
    estimate = optimize(neg_log_like, c(pi1, 0.9999), tol = 0.000001, R2 = R2, R1 = R1,
                        pi1 = pi1, n = n, alpha1 = alpha1, alpha2 = alpha2,
                        beta1 = beta1, beta2 = beta2)$minimum
    I_fisher = (n*(alpha2 + beta2 - 1)^2)/(alpha2*estimate - estimate - alpha2 + beta2*estimate + 1) + (n*(alpha2 + beta2 - 1)^2)/(alpha2 + estimate - alpha2*estimate - beta2*estimate) + (n*pi1^2*(alpha1 + beta1 - 1)^2*(alpha2 + estimate - alpha2*estimate - beta2*estimate))/(estimate^3*(estimate - pi1 + alpha1*pi1 - alpha1*estimate + beta1*pi1)) + (n*pi1^2*(alpha1 + beta1 - 1)^2*(alpha2 + estimate - alpha2*estimate - beta2*estimate))/(estimate^3*(pi1 - alpha1*pi1 + alpha1*estimate - beta1*pi1))
    sd = 1/sqrt(I_fisher)
    ci = estimate + c(-1, 1)*qnorm(1 - gamma/2)*sd
  }else{
    estimate = mle_auxiliary(R1, R2, n, pi1)
    sd = sqrt(((estimate - pi1)*(1 - estimate))/(n*(1 - pi1)))
    ci = estimate + c(-1, 1)*qnorm(1 - gamma/2)*sd
  }

  I2 = qbeta(p = 1 - gamma/2, R2 - R1 + 1, n - R2 + R1)
  I1 = qbeta(p = gamma/2, R2 - R1, n - R2 + R1 + 1)
  Delta1 = 1 - alpha1 - beta1
  Delta2 = 1 - alpha2 - beta2
  dlt = pi1*Delta1*(Delta2 + alpha2/estimate)
  lower = ((I1 + dlt)/(1 - alpha1) - alpha2)/Delta2
  upper = ((I2 + dlt)/(1 - alpha1) - alpha2)/Delta2
  ci2 = c(lower, upper)

  list(estimate = estimate, sd = sd, ci = ci, ci2 = ci2, gamma = gamma, ...)
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
#' modified_mle(R1 = X$R1, R2 = X$R2, pi1 = X$pi1, n = X$n, alpha1 = 0.01,
#' alpha2 = 0.01, beta1 = 0.05, beta2 = 0.05)
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

  if (alpha2 == 0){
    sd = sqrt(((estimate*(1 - alpha1) - pi1*Delta1)*(1 + pi1*(1 - beta2)*Delta1 - estimate*(1 - alpha1)*(1 - beta2)))/(n*(1 - alpha1)^2*(1 - beta2)))
    I2 = qbeta(p = 1 - gamma/2, R2 - R1 + 1, n - R2 + R1)
    I1 = qbeta(p = gamma/2, R2 - R1, n - R2 + R1 + 1)
    upper = (I2/(1 - beta2) + pi1*Delta1)/(1 - alpha1)
    lower = (I1/(1 - beta2) + pi1*Delta1)/(1 - alpha1)
  }else{
    pi2_star = estimate*Delta1 + alpha2
    pi2_corrected = pi2_star*(1 - (pi1/estimate*(1 - beta1) - (1 - pi1/estimate)*alpha1))
    var_R2_R1 = n*pi2_corrected*(1 - pi2_corrected)
    c_coef = pi1*Delta1*Delta2 - alpha2*(1 - alpha2)
    coef = 1/(n*a3)*(1 + (pi2_corrected + c_coef)/sqrt(a1*a3 + (c_coef + pi2_corrected)^2))
    sd = sqrt(coef^2*var_R2_R1)
    lower = estimate - qnorm(1 - gamma/2)*sd
    upper = estimate + qnorm(1 - gamma/2)*sd
  }
  list(estimate = estimate, sd = sd, ci = c(lower, upper), gamma = gamma, ...)
}









