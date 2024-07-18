library('dplyr')
library('tidyr')
library('ggplot2')
library('numDeriv')
library('reshape2')

## TODO
## 1. Use code from ACTIV 4 HT to simulate power of trial with no interim analyses
## 2. Calculate power of same trial using MLE method

## two-sided power against null of zero
power_mle <- function(n, par, var, n_sde, alp=0.05) {
  z <- qnorm(1-alp/2)
  s <- sqrt(var * n_sde/n)
  1-pnorm(z, par/s)+pnorm(-z, par/s)
}

fit <- glm(Survived ~ Age + Sex + Class, family=binomial(), data=Titanic_df)
summary(fit)

# log-likelihood function
logLikFun <- function(params, model) {
  # design matrix
  mmx <- model.matrix(model)
  # linear predictor and GLM parameter
  eta <- mmx %*% params
  mu  <- model$family$linkinv(eta)
  # log-likelihood
  -0.5*fit$family$aic(y=model$y, mu=mu, wt=1)
}

## verify log-likelihood calculation
logLik(fit)
logLikFun(coef(fit), fit)

## vcov() is based on negative inverse Hessian
hess <- numDeriv::hessian(logLikFun, coef(fit), model=fit)
sqrt(diag(-solve(hess)))
sqrt(diag(vcov(fit)))

## vcov() is used to calculate p-values, i.e., is assumed to be
## the variance-covariance under the null hypothesis
coef(fit) / sqrt(diag(vcov(fit)))  ## z-value
2*(1-pnorm(abs(coef(fit) / sqrt(diag(vcov(fit)))))) ## p-value
summary(fit)
