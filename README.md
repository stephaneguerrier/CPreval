
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Travis-CI Build
Status](https://travis-ci.org/SMAC-Group/simts.svg?branch=master)](https://travis-ci.org/SMAC-Group/simts)
[![Licence](https://img.shields.io/badge/licence-AGPL--3.0-blue.svg)](https://opensource.org/licenses/AGPL-3.0)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-3.4.0-6666ff.svg)](https://cran.r-project.org/)
[![Last-changedate](https://img.shields.io/badge/last%20change-2020--06--01-green.svg)](https://github.com/SMAC-Group/simts)

# `CPrval` Overview

This R package provides an implementation of the estimators discussed in
*Accurate Prevalence Estimation for Emerging or Rare Infectious
Diseases* by Stéphane Guerrier, Christoph Kuzmics and Maria-Pia
Victoria-Feser (submitted manuscript available upon request). The
pacakge can be installed from GitHub as follows:

``` r
# Install devtools
install.packages("devtools")

# Install the package from GitHub
devtools::install_github("stephaneguerrier/CPreval")
```

## Example

Consider the following data:

``` r
seed = 18
n    = 1500
p0   = 3/100
pi0  = 1/100
(X = sim_Rs(p = p0, pi0 = pi0, n = n, seed = seed))
#> $R0
#> [1] 12
#> 
#> $R
#> [1] 40
#> 
#> $pi0
#> [1] 0.01
#> 
#> $n
#> [1] 1500
#> 
#> $alpha0
#> [1] 0
#> 
#> $alpha
#> [1] 0
#> 
#> $beta0
#> [1] 0
#> 
#> $beta
#> [1] 0
```

Estimators:

``` r
survey_sample(X$R, X$n)
#> Method: Survey sample
#> 
#> Estimated proportion: 2.6667%
#> Standard error      : 0.4160%
#> 
#> Confidence intervals with gamma = 0.05:
#> Asymptotic Approach: 1.8514% - 3.4820%
#> Clopper–Pearson    : 1.9118% - 3.6137%
#> 
#> Assumed measurement error: alpha = 0%, beta = 0%
moment_estimator(R0 = X$R0, R = X$R, pi0 = X$pi0, n = X$n)
#> Method: Moment estimator
#> 
#> Estimated proportion: 2.8667%
#> Standard error      : 0.3495%
#> 
#> Confidence intervals with gamma = 0.05:
#> Asymptotic Approach: 2.1817% - 3.5516%
#> Clopper–Pearson    : 2.2439% - 3.6866%
#> 
#> Assumed measurement error: alpha0 = 0%, alpha = 0%, 
#>                            beta0  = 0%, beta  = 0%
mle(R0 = X$R0, R = X$R, pi0 = X$pi0, n = X$n)
#> Method: MLE
#> 
#> Estimated proportion: 2.8629%
#> Standard error      : 0.3491%
#> 
#> Confidence intervals with gamma = 0.05:
#> Asymptotic Approach: 2.1787% - 3.5471%
#> Clopper–Pearson    : 2.2439% - 3.6866%
#> 
#> Assumed measurement error: alpha0 = 0%, alpha = 0%, 
#>                            beta0  = 0%, beta  = 0%
```

With measurement error:

``` r
alpha0 = 1/100
alpha  = 1/100
beta0  = 6/100
beta   = 6/100
(X = sim_Rs(p = p0, pi0 = pi0, n = n, seed = seed,
           alpha = alpha, alpha0 = alpha0, 
           beta = beta, beta0 = beta0))
#> $R0
#> [1] 15
#> 
#> $R
#> [1] 51
#> 
#> $pi0
#> [1] 0.01
#> 
#> $n
#> [1] 1500
#> 
#> $alpha0
#> [1] 0.01
#> 
#> $alpha
#> [1] 0.01
#> 
#> $beta0
#> [1] 0.06
#> 
#> $beta
#> [1] 0.06
```

Without taking into account measurement error:

``` r
survey_sample(X$R, X$n)
#> Method: Survey sample
#> 
#> Estimated proportion: 3.4000%
#> Standard error      : 0.4679%
#> 
#> Confidence intervals with gamma = 0.05:
#> Asymptotic Approach: 2.4829% - 4.3171%
#> Clopper–Pearson    : 2.5418% - 4.4463%
#> 
#> Assumed measurement error: alpha = 0%, beta = 0%
moment_estimator(R0 = X$R0, R = X$R, pi0 = X$pi0, n = X$n)
#> Method: Moment estimator
#> 
#> Estimated proportion: 3.4000%
#> Standard error      : 0.3952%
#> 
#> Confidence intervals with gamma = 0.05:
#> Asymptotic Approach: 2.6255% - 4.1745%
#> Clopper–Pearson    : 2.6865% - 4.3072%
#> 
#> Assumed measurement error: alpha0 = 0%, alpha = 0%, 
#>                            beta0  = 0%, beta  = 0%
mle(R0 = X$R0, R = X$R, pi0 = X$pi0, n = X$n)
#> Method: MLE
#> 
#> Estimated proportion: 3.4000%
#> Standard error      : 0.3951%
#> 
#> Confidence intervals with gamma = 0.05:
#> Asymptotic Approach: 2.6256% - 4.1744%
#> Clopper–Pearson    : 2.6865% - 4.3072%
#> 
#> Assumed measurement error: alpha0 = 0%, alpha = 0%, 
#>                            beta0  = 0%, beta  = 0%
```

Taking into account measurement error:

``` r
survey_sample(X$R, X$n, alpha = alpha, beta = beta)
#> Method: Survey sample
#> 
#> Estimated proportion: 2.6559%
#> Standard error      : 0.4679%
#> 
#> Confidence intervals with gamma = 0.05:
#> Asymptotic Approach: 1.7388% - 3.5730%
#> Clopper–Pearson    : 1.6578% - 3.7057%
#> 
#> Assumed measurement error: alpha = 1%, beta = 6%
moment_estimator(R0 = X$R0, R = X$R, pi0 = X$pi0, n = X$n, 
                 alpha = alpha, alpha0 = alpha0, 
                 beta = beta, beta0 = beta0)
#> Method: Moment estimator
#> 
#> Estimated proportion: 2.8280%
#> Standard error      : 0.3860%
#> 
#> Confidence intervals with gamma = 0.05:
#> Asymptotic Approach: 2.0714% - 3.5847%
#> Clopper–Pearson    : 2.0530% - 3.8134%
#> 
#> Assumed measurement error: alpha0 = 1%, alpha = 1%, 
#>                            beta0  = 6%, beta  = 6%
mle(R0 = X$R0, R = X$R, pi0 = X$pi0, n = X$n, 
                 alpha = alpha, alpha0 = alpha0, 
                 beta = beta, beta0 = beta0)
#> Method: MLE
#> 
#> Estimated proportion: 2.8658%
#> Standard error      : 0.3823%
#> 
#> Confidence intervals with gamma = 0.05:
#> Asymptotic Approach: 2.1165% - 3.6152%
#> Clopper–Pearson    : 2.0483% - 3.8087%
#> 
#> Assumed measurement error: alpha0 = 1%, alpha = 1%, 
#>                            beta0  = 6%, beta  = 6%
```

## License

The license this source code is released under is the GNU AFFERO GENERAL
PUBLIC LICENSE (AGPL) v3.0. Please see the LICENSE file for full text.
Otherwise, please consult [TLDR
Legal](https://tldrlegal.com/license/gnu-affero-general-public-license-v3-\(agpl-3.0\))
or [GNU](https://www.gnu.org/licenses/agpl-3.0.en.html) which will
provide a synopsis of the restrictions placed upon the code.
