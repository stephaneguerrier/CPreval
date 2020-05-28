#' @title Compute proportion in the survey sample
#' @description Proportion estimated using the survey sample and confidence intervals based on the Clopper–Pearson approach.
#' @param R        A \code{numeric} that provides the people of positive people in the sample.
#' @param n        A \code{numeric} that provides the sample size.
#' @param alpha    A \code{numeric} that provides the False Negative (FN) rate for the sample R. Default value is \code{0}.
#' @param beta0     A \code{numeric} that provides the False Positive (FP) rate for the sample R. Default value is \code{0}.
#' @param gamma    A \code{numeric} that used to compute a (1 - gamma) confidence region for the proportion. Default value is \code{0.05}.
#' @param ...      Additional arguments.
#' @return A \code{CPreval} object with the structure:
#' \itemize{
#'  \item estimate: Estimated proportion
#'  \item sd: Estimated standard error of the estimator
#'  \item ci_asy: Asymptotic confidence interval
#'  \item ci_cp: Confidence interval based on the Clopper–Pearson approach
#'  \item gamma: Confidence level (i.e. 1 - gamma) for confidence interval
#'  \item method: Estimation method (in this case sample survey)
#'  \item measurement: A vector with (alpha0, alpha, beta0, beta)
#'  \item ...: Additional parameters
#' }
#' @export
#' @author Stephane Guerrier
#' @examples
#' # Samples without measurement error
#' X = sim_Rs(p = 30/1000, pi0 = 10/1000, n = 1500, seed = 18)
#' survey_sample(X$R, X$n)
#'
#' # With measurement error
#' X = sim_Rs(p = 30/1000, pi0 = 10/1000, n = 1500, alpha0 = 0.01,
#' alpha = 0.01, beta0 = 0.05, beta = 0.05, seed = 18)
#' survey_sample(X$R, X$n)
#' survey_sample(X$R, X$n, alpha = 0.01, beta = 0.05)
survey_sample = function(R, n, alpha = 0, beta = 0, gamma = 0.05, ...){
  # Compute survey proportion
  pi_bar = R/n

  # Adjust for false positive/negative
  if (max(alpha, beta) > 0){
    pi_bar = pi_bar/(1 - alpha - beta) - alpha
  }

  # Estimated standard error
  sd = sqrt(R/n*(1 - R/n)/n)

  # Compute 1 - gamma asymptotic interval
  ci_asym = pi_bar + qnorm(1 - gamma/2)*c(-1,1)*sd

  # Compute 1 - gamma confidence interval - Clopper–Pearson approach
  upper = (qbeta(p = 1 - gamma/2, R + 1, n - R) - alpha)/(1 - alpha - beta)
  lower = (qbeta(p = gamma/2, R, n - R + 1) - alpha)/(1 - alpha - beta)
  ci_cp = c(lower, upper)

  out = list(estimate = pi_bar, sd = sd, ci_asym = ci_asym, ci_cp = ci_cp, gamma = gamma,
             method = "Survey sample", measurement = c(NA, alpha, NA, beta), ...)
  class(out) = "CPreval"
  out
}

#' @title Print sample survey
#' @description Simple print function for the sample survey estimator
#' @method print CPreval
#' @export
#' @keywords internal
#' @param x    A \code{CPreval} object
#' @param ...  Further arguments passed to or from other methods
#' @return Prints the sample survey
#' @author Stephane Guerrier
#' @examples
#' X = sim_Rs(p = 3/100, pi0 = 1/100, n = 1500, seed = 18)
#' survey_sample(X$R, X$n)
print.CPreval = function(x, ...){
  cat("Method: ")
  cat(x$method)
  cat("\n\n")
  cat("Estimated proportion: ")
  cat(sprintf("%.4f", 100*x$estimate))
  cat("%\n")
  cat("Standard error      : ")
  cat(sprintf("%.4f", 100*x$sd))
  cat("%\n\n")
  cat("Confidence intervals with gamma = ")
  cat(x$gamma)
  cat(":\n")
  cat("Asymptotic Approach: ")
  cat(sprintf("%.4f", 100*x$ci_asym[1]))
  cat("% - ")
  cat(sprintf("%.4f", 100*x$ci_asym[2]))
  cat("%\n")
  cat("Clopper–Pearson    : ")
  cat(sprintf("%.4f", 100*x$ci_cp[1]))
  cat("% - ")
  cat(sprintf("%.4f", 100*x$ci_cp[2]))
  cat("%\n\n")

  if (x$method == "Survey sample"){
    cat("Assumed measurement error: alpha = ")
    cat(100*x$measurement[2])
    cat("%, beta = ")
    cat(100*x$measurement[4])
    cat("% \n\n")
  }else{
    cat("Assumed measurement error: alpha0 = ")
    cat(100*x$measurement[1])
    cat("%, alpha = ")
    cat(100*x$measurement[2])
    cat("%, \n")
    cat("                          beta0  = ")
    cat(100*x$measurement[3])
    cat("%, beta  = ")
    cat(100*x$measurement[4])
    cat("% \n")
  }
}


#' @title Compute moment-based estimator based on R0 and R
#' @description Proportion estimated using the moment-based estimator and confidence intervals based the asymptotic distribution of the estimator as well as
#' the Clopper–Pearson approach.
#' @param R0        A \code{numeric} that provides the people of positive people in the sample that were known to be positive.
#' @param R         A \code{numeric} that provides the people of positive people in the sample.
#' @param n         A \code{numeric} that provides the sample size.
#' @param alpha0    A \code{numeric} that provides the False Negative (FN) rate for the sample R0. Default value is \code{0}.
#' @param alpha     A \code{numeric} that provides the False Negative (FN) rate for the sample R. Default value is \code{0}.
#' @param beta0     A \code{numeric} that provides the False Positive (FP) rate for the sample R0. Default value is \code{0}.
#' @param beta      A \code{numeric} that provides the False Positive (FP) rate for the sample R. Default value is \code{0}.
#' @param gamma     A \code{numeric} that used to compute a (1 - gamma) confidence region for the proportion. Default value is \code{0.05}.
#' @param ...       Additional arguments.
#' @return A \code{CPreval} object with the structure:
#' \itemize{
#'  \item estimate: Estimated proportion
#'  \item sd: Estimated standard error of the estimator
#'  \item ci_asy: Asymptotic confidence interval
#'  \item ci_cp: Confidence interval based on the Clopper–Pearson approach
#'  \item gamma: Confidence level (i.e. 1 - gamma) for confidence interval
#'  \item method: Estimation method (in this case sample survey)
#'  \item measurement: A vector with (alpha0, alpha, beta0, beta)
#'  \item ...: Additional parameters
#' }
#' @export
#' @author Stephane Guerrier
#' @examples
#' # Samples without measurement error
#' X = sim_Rs(p = 3/100, pi0 = 1/100, n = 1500, seed = 18)
#' moment_estimator(R0 = X$R0, R = X$R, pi0 = X$pi0, n = X$n)
#'
#' # With measurement error
#' X = sim_Rs(p = 30/1000, pi0 = 10/1000, n = 1500, alpha0 = 0.01,
#' alpha = 0.01, beta0 = 0.05, beta = 0.05, seed = 18)
#' moment_estimator(R0 = X$R0, R = X$R, pi0 = X$pi0, n = X$n)
#' moment_estimator(R0 = X$R0, R = X$R, pi0 = X$pi0, n = X$n, alpha0 = 0.01,
#' alpha = 0.01, beta0 = 0.05, beta = 0.05)
moment_estimator = function(R0, R, pi0, n, alpha0 = 0, alpha = 0, beta0 = 0, beta = 0, gamma = 0.05, ...){
  # Compute point estimate
  # Define Deltas
  Delta0 = 1 - alpha0 - beta0
  Delta = 1 - alpha - beta

  # Define as
  a1 = 2*alpha*Delta0*pi0
  a2 = (R - R0)/n + pi0*Delta0*Delta - alpha*(1 - alpha)
  a3 = 2*Delta*(1 - alpha0)

  # Estimator
  estimate = (a2 + sqrt(a2^2 + a1*a3))/a3

  if (alpha == 0){
    # Compute sd
    sd = sqrt(((estimate*(1 - alpha0) - pi0*Delta0)*(1 + pi0*(1 - beta0)*Delta0 - estimate*(1 - alpha0)*(1 - beta)))/(n*(1 - alpha0)^2*(1 - beta)))

    # Compute 1 - gamma asymptotic interval
    ci_asym = estimate + qnorm(1 - gamma/2)*c(-1,1)*sd

    # Compute 1 - gamma confidence interval - Clopper–Pearson approach
    I2 = qbeta(p = 1 - gamma/2, R - R0 + 1, n - R + R0)
    I1 = qbeta(p = gamma/2, R - R0, n - R + R0 + 1)
    upper = (I2/(1 - beta) + pi0*Delta0)/(1 - alpha0)
    lower = (I1/(1 - beta) + pi0*Delta0)/(1 - alpha0)
    ci_cp = c(lower, upper)
  }else{
    # Compute sd
    pi2_star = estimate*Delta0 + alpha
    pi2_corrected = pi2_star*(1 - (pi0/estimate*(1 - beta0) - (1 - pi0/estimate)*alpha0))
    var_R_R0 = n*pi2_corrected*(1 - pi2_corrected)
    c_coef = pi0*Delta0*Delta - alpha*(1 - alpha)
    coef = 1/(n*a3)*(1 + (pi2_corrected + c_coef)/sqrt(a1*a3 + (c_coef + pi2_corrected)^2))
    sd = sqrt(coef^2*var_R_R0)

    # Compute 1 - gamma asymptotic interval
    ci_asym = estimate + qnorm(1 - gamma/2)*c(-1,1)*sd

    # Compute 1 - gamma confidence interval - Clopper–Pearson approach
    I2 = qbeta(p = 1 - gamma/2, R - R0 + 1, n - R + R0)
    I1 = qbeta(p = gamma/2, R - R0, n - R + R0 + 1)
    Delta0 = 1 - alpha0 - beta0
    Delta = 1 - alpha - beta
    dlt = pi0*Delta0*(Delta + alpha/estimate)
    lower = ((I1 + dlt)/(1 - alpha0) - alpha)/Delta
    upper = ((I2 + dlt)/(1 - alpha0) - alpha)/Delta
    ci_cp = c(lower, upper)
  }

  out = list(estimate = estimate, sd = sd, ci_asym = ci_asym, ci_cp = ci_cp, gamma = gamma,
             method = "Moment estimator", measurement = c(alpha0, alpha, beta0, beta), ...)
  class(out) = "CPreval"
  out
}
