#' @title Compute proportion in the survey sample
#' @description Proportion estimated using the survey sample and confidence intervals based on the Clopper–Pearson approach.
#' @param R        A \code{numeric} that provides the people of positive people in the sample.
#' @param n        A \code{numeric} that provides the sample size.
#' @param alpha    A \code{numeric} that provides the False Negative (FN) rate for the sample R. Default value is \code{0}.
#' @param beta0    A \code{numeric} that provides the False Positive (FP) rate for the sample R. Default value is \code{0}.
#' @param gamma    A \code{numeric} that used to compute a (1 - gamma) confidence region for the proportion. Default value is \code{0.05}.
#' @param ...      Additional arguments.
#' @return A \code{CPreval} object with the structure:
#' \itemize{
#'  \item estimate: Estimated proportion
#'  \item sd: Estimated standard error of the estimator
#'  \item ci_asym: Asymptotic confidence interval
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
survey_sample = function(R, n, alpha = 0, beta = 0, gamma = 0.05, pi0 = 0, simulation = FALSE, ...){
  # Compute survey proportion
  pi_bar = R/n

  # Adjust for false positive/negative
  if (max(alpha, beta) > 0){
    pi_bar = (pi_bar - alpha)/(1 - alpha - beta)
  }

  if (pi_bar < pi0 || pi_bar > 1){
    if (simulation == FALSE){
      warning("The estimator has no (reliable) solution.")
    }
    pi_bar = c(pi0, 1)[which.min(abs(pi_bar - c(pi0, 1)))]
    sd = NA
    ci_asym = NA

    if (pi_bar == pi0){
      upper = (qbeta(p = 1 - gamma/2, R + 1, n - R) - alpha)/(1 - alpha - beta)
      lower = pi0
    }else{
      upper = 1
      lower = (qbeta(p = gamma/2, R, n - R + 1) - alpha)/(1 - alpha - beta)
    }

    ci_cp = c(lower, upper)

    # Construct output
    out = list(estimate = pi_bar, sd = sd, ci_asym = ci_asym, ci_cp = ci_cp, gamma = gamma,
               method = "Survey sample", measurement = c(NA, alpha, NA, beta),
               boundary = TRUE,...)
    class(out) = "CPreval"
    return(out)
  }

  # Estimated standard error
  sd = sqrt(pi_bar*(1 - pi_bar)/n)

  # Compute 1 - gamma asymptotic interval
  ci_asym = pi_bar + qnorm(1 - gamma/2)*c(-1,1)*sd

  # Compute 1 - gamma confidence interval - Clopper–Pearson approach
  upper = (qbeta(p = 1 - gamma/2, R + 1, n - R) - alpha)/(1 - alpha - beta)
  lower = (qbeta(p = gamma/2, R, n - R + 1) - alpha)/(1 - alpha - beta)

  # Adjust CI - CP
  if (lower < pi0){
    lower = pi0
  }

  if (upper > 1){
    upper = 1
  }

  ci_cp = c(lower, upper)

  # Construct output
  out = list(estimate = pi_bar, sd = sd, ci_asym = ci_asym, ci_cp = ci_cp, gamma = gamma,
             method = "Survey sample", measurement = c(NA, alpha, NA, beta),
             boundary = FALSE,...)
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
  cat("Confidence intervals at the ")
  cat(100*(1 - x$gamma))
  cat("% level:\n")
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
    cat("                           beta0  = ")
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
#'  \item ci_asym: Asymptotic confidence interval
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
moment_estimator = function(R0, R, pi0, n, alpha0 = 0, alpha = 0, beta0 = 0, beta = 0, gamma = 0.05, simulation = FALSE, ...){
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

  # Check estimator
  if (estimate < pi0 || estimate > 1){
    if (simulation == FALSE){
      warning("The estimator has no (reliable) solution.")
    }
    estimate = c(pi0, 1)[which.min(abs(estimate - c(pi0, 1)))]
    sd = NA
    ci_asym = NA
    ci_cp = NA
    out = list(estimate = estimate, sd = sd, ci_asym = ci_asym, ci_cp = ci_cp, gamma = gamma,
               method = "Moment estimator", measurement = c(alpha0, alpha, beta0, beta),
               boundary = TRUE, ...)
    class(out) = "CPreval"
    return(out)
  }

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
             method = "Moment estimator", measurement = c(alpha0, alpha, beta0, beta),
             boundary = FALSE, ...)
  class(out) = "CPreval"
  out
}

#' @title Negative Log-Likelihood function based on R0 and R
#' @description Log-Likelihood function based on R0 and R multiplied by -1.
#' @param pi2       A \code{numeric} value for p0.
#' @param R0        A \code{numeric} that provides the people of positive people in the sample that were known to be positive.
#' @param R         A \code{numeric} that provides the people of positive people in the sample.
#' @param n         A \code{numeric} that provides the sample size.
#' @param alpha0    A \code{numeric} that provides the False Negative (FN) rate for the sample R0. Default value is \code{0}.
#' @param alpha     A \code{numeric} that provides the False Negative (FN) rate for the sample R. Default value is \code{0}.
#' @param beta0     A \code{numeric} that provides the False Positive (FP) rate for the sample R0. Default value is \code{0}.
#' @param beta      A \code{numeric} that provides the False Positive (FP) rate for the sample R. Default value is \code{0}.
#' @param ...       Additional arguments.
#' @return Negative Log-Likelihood
#' @author Stephane Guerrier
neg_log_like = function(pi2, R, R0, pi0, n, alpha0, alpha, beta0, beta, ...){
  # Define Deltas
  Delta1 = 1 - alpha0 - beta0
  Delta2 = 1 - alpha - beta

  # Define modified probs
  pi2_star = pi2*Delta2 + alpha
  pi_star = pi0/pi2*Delta1 + alpha0

  # Log likelihood times (-1)
  -(R*log(pi2_star) + (n - R)*log(1 - pi2_star) + R0*log(pi_star) + (R - R0)*log(1 - pi_star))
}

#' @title Compute MLE based on R0 and R
#' @description Proportion estimated using the MLE and confidence intervals based the asymptotic distribution of the estimator
#' as well as the Clopper–Pearson approach.
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
#'  \item ci_asym: Asymptotic confidence interval
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
#' mle(R0 = X$R0, R = X$R, pi0 = X$pi0, n = X$n)
#'
#' # With measurement error
#' X = sim_Rs(p = 30/1000, pi0 = 10/1000, n = 1500, alpha0 = 0.01,
#' alpha = 0.01, beta0 = 0.05, beta = 0.05, seed = 18)
#' mle(R0 = X$R0, R = X$R, pi0 = X$pi0, n = X$n)
#' mle(R0 = X$R0, R = X$R, pi0 = X$pi0, n = X$n, alpha0 = 0.01,
#' alpha = 0.01, beta0 = 0.05, beta = 0.05)
mle = function(R0, R, pi0, n, alpha0 = 0, alpha = 0,
               beta0 = 0, beta = 0, gamma = 0.05, simulation = FALSE, ...){


  if (max(alpha0, alpha, beta0, beta) > 0){
    eps = 10^(-5)
    estimate = optimize(neg_log_like, c(pi0 - eps, 1 - eps), tol = 0.000001, R = R, R0 = R0,
                        pi0 = pi0, n = n, alpha0 = alpha0, alpha = alpha,
                        beta0 = beta0, beta = beta)$minimum

    # Check boundary
    boundary_lower = neg_log_like(pi2 = pi0, R = R, R0 = R0, pi0 = pi0, n = n,
                 alpha0 = 0.01, alpha = 0.01, beta0 = 0.05, beta = 0.05)
    boundary_upper = neg_log_like(pi2 = 1, R = R, R0 = R0, pi0 = pi0, n = n,
                 alpha0 = 0.01, alpha = 0.01, beta0 = 0.05, beta = 0.05)
    min_eval = neg_log_like(pi2 = estimate, R = R, R0 = R0, pi0 = pi0, n = n,
                 alpha0 = 0.01, alpha = 0.01, beta0 = 0.05, beta = 0.05)

    if (min(boundary_lower, boundary_upper) <= min_eval){
      if (simulation == FALSE){
        warning("The estimator has no (reliable) solution.")
      }

      estimate = c(pi0, 1)[which.min(c(boundary_lower, boundary_upper))]

      sd = NA
      ci_asym = NA

      if (estimate == pi0){
        upper = (qbeta(p = 1 - gamma/2, R + 1, n - R) - alpha)/(1 - alpha - beta)
        lower = pi0
      }else{
        upper = 1
        lower = (qbeta(p = gamma/2, R, n - R + 1) - alpha)/(1 - alpha - beta)
      }
      ci_cp = c(lower, upper)

      out = list(estimate = estimate, sd = sd, ci_asym = ci_asym, ci_cp = ci_cp, gamma = gamma,
                 method = "MLE", measurement = c(alpha0, alpha, beta0, beta),
                 boundary = TRUE,...)
      class(out) = "CPreval"
      return(out)
    }

    I_fisher = (n*(alpha + beta - 1)^2)/(alpha*estimate - estimate - alpha + beta*estimate + 1) + (n*(alpha + beta - 1)^2)/(alpha + estimate - alpha*estimate - beta*estimate) + (n*pi0^2*(alpha0 + beta0 - 1)^2*(alpha + estimate - alpha*estimate - beta*estimate))/(estimate^3*(estimate - pi0 + alpha0*pi0 - alpha0*estimate + beta0*pi0)) + (n*pi0^2*(alpha0 + beta0 - 1)^2*(alpha + estimate - alpha*estimate - beta*estimate))/(estimate^3*(pi0 - alpha0*pi0 + alpha0*estimate - beta0*pi0))
    sd = 1/sqrt(I_fisher)
    ci_asym = estimate + c(-1, 1)*qnorm(1 - gamma/2)*sd
  }else{
    estimate = pi0*(n - R)/(n - R0) + (R - R0)/(n - R0)

    if (estimate < pi0 || estimate > 1){
      if (simulation == FALSE){
        warning("The estimator has no (reliable) solution.")
      }

      estimate = c(pi0, 1)[which.min(abs(estimate - c(pi0, 1)))]
      sd = NA
      ci_asym = NA

      if (estimate == pi0){
        upper = (qbeta(p = 1 - gamma, R + 1, n - R) - alpha)/(1 - alpha - beta)
        lower = pi0
      }else{
        upper = 1
        lower = (qbeta(p = gamma, R, n - R + 1) - alpha)/(1 - alpha - beta)
      }
      ci_cp = c(lower, upper)

      out = list(estimate = estimate, sd = sd, ci_asym = ci_asym, ci_cp = ci_cp, gamma = gamma,
                 method = "MLE", measurement = c(alpha0, alpha, beta0, beta),
                 boundary = TRUE,...)
      class(out) = "CPreval"
      return(out)
    }

    sd = sqrt(((estimate - pi0)*(1 - estimate))/(n*(1 - pi0)))
    ci_asym = estimate + c(-1, 1)*qnorm(1 - gamma/2)*sd
}

  I2 = qbeta(p = 1 - gamma/2, R - R0 + 1, n - R + R0)
  I1 = qbeta(p = gamma/2, R - R0, n - R + R0 + 1)
  Delta1 = 1 - alpha0 - beta0
  Delta2 = 1 - alpha - beta
  dlt = pi0*Delta1*(Delta2 + alpha/estimate)
  lower = ((I1 + dlt)/(1 - alpha0) - alpha)/Delta2
  upper = ((I2 + dlt)/(1 - alpha0) - alpha)/Delta2

  if (lower < pi0){
    lower = pi0
  }

  if (upper > 1){
    upper = 1
  }

  ci_cp = c(lower, upper)

  out = list(estimate = estimate, sd = sd, ci_asym = ci_asym, ci_cp = ci_cp, gamma = gamma,
             method = "MLE", measurement = c(alpha0, alpha, beta0, beta),
             boundary = FALSE, ...)
  class(out) = "CPreval"
  out
}

