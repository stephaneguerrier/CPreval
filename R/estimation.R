#' @title Compute proportion in the survey sample
#' @description Proportion estimated using the survey sample and confidence intervals based on the Clopper–Pearson approach.
#' @param R        A \code{numeric} that provides the people of positive people in the sample.
#' @param n        A \code{numeric} that provides the sample size.
#' @param alpha    A \code{numeric} that provides the False Negative (FN) rate for the sample R. Default value is \code{0}.
#' @param beta     A \code{numeric} that provides the False Positive (FP) rate for the sample R. Default value is \code{0}.
#' @param gamma    A \code{numeric} that used to compute a (1 - gamma) confidence region for the proportion. Default value is \code{0.05}.
#' @param ...      Additional arguments.
#' @return A \code{list} with the structure:
#' \itemize{
#'  \item estimate: Estimated proportion
#'  \item sd: Estimated standard error of the estimator
#'  \item ci_asy: Asymptotic confidence interval
#'  \item ci_cp: Confidence interval based on the Clopper–Pearson approach
#'  \item gamma: Confidence level (i.e. 1 - gamma) for confidence interval
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

  out = list(estimate = pi_bar, sd = sd, ci_asym = ci_asym, ci_cp = ci_cp, gamma = gamma, ...)
  class(out) = "sample_survey"
  out
}

#' @title Print sample survey
#' @description Simple print function for the sample survey estimator
#' @method print sample_survey
#' @export
#' @keywords internal
#' @param x    A \code{sample_survey} object
#' @param ...  Further arguments passed to or from other methods
#' @return Prints the sample survey
#' @author Stephane Guerrier
#' @examples
#' X = sim_Rs(p = 30/1000, pi0 = 10/1000, n = 1500, seed = 18)
#' survey_sample(X$R, X$n)
print.sample_survey = function(x, ...){
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
  cat("%\n")
}


