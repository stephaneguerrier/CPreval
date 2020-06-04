#' @title Plot function a
#' @description Plot function a to assess the identifiability of the MLE
#' @param pi0       A \code{numeric} that provides the proportion of people (in the whole population) who are positive, as measured through a non-random,
#' but systematic sampling (e.g. based on medical selection).
#' @param alpha0    A \code{numeric} that provides the False Negative (FN) rate for the sample R0. Default value is \code{0}.
#' @param alpha0    A \code{numeric} that provides the False Negative (FN) rate for the sample R0. Default value is \code{0}.
#' @param beta      A \code{numeric} that provides the False Positive (FP) rate for the sample R. Default value is \code{0}.
#' @param beta      A \code{numeric} that provides the False Positive (FP) rate for the sample R. Default value is \code{0}.
#' @param ...       Additional arguments.
#' @return A plot of the function a
#' @export
#' @author Stephane Guerrier
#' @examples
#' plot_a(pi0 = 1/100, alpha0 = 0.005, alpha = 0.005, beta0 = 0.05, beta = 0.05)
plot_a = function(pi0, alpha, alpha0, beta, beta0, ...){
  Delta = 1 - alpha - beta
  Delta0 = 1 - alpha0 - beta0
  delta = 0.0001
  m = 100
  my_pi = p0 = seq(from = pi0 + delta, to = 1 - delta, length.out = m)
  A = matrix(NA, m, m)
  for (i in 1:m){
    for (j in 1:m){
      a1 = Delta^2*Delta0*pi0*p0[j]*my_pi[i]^2*(1 + Delta*Delta0*pi0*my_pi[i])
      a2 = Delta*Delta0*alpha^2*pi0^2*(my_pi[i]*Delta*(my_pi[i] + 2*p0[j]) - (my_pi[i] + p0[j]))
      a3 = Delta0*alpha^2*pi0^2*(Delta*(2*my_pi[i] + p0[j]) - 1 + alpha)
      a4 = Delta^2*alpha0*my_pi[i]^2*p0[j]*(2*Delta0*pi0 - my_pi[i]*(1 - alpha0))
      A[i,j] = a1 + a2 + a3 + a4
    }
  }

  plot(NA, xlim = c(pi0 + delta, 1 - delta), ylim = c(pi0 + delta, 1 - delta),
       xlab = expression(pi), ylab = expression(p[0]))
  contour(my_pi, p0, A, add = TRUE)
}

#' @export
gui = function(){
  appDir = system.file("gui", package = "CPreval")
  shiny::runApp(appDir, display.mode = "normal")
}
