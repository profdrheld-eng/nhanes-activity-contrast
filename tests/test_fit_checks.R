#!/usr/bin/env Rscript
# Synthetic failure detection, including finite coefficients from a separated Cox fit.
suppressPackageStartupMessages({library(survey);library(survival)})
arg <- sub('^--file=','',commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))])
root <- dirname(dirname(normalizePath(arg)))
source(file.path(root,'R','fit_checks.R'))
expect_error <- function(expr,pattern) {
  err <- tryCatch({force(expr);NULL},error=identity)
  stopifnot(inherits(err,'error'),grepl(pattern,conditionMessage(err),ignore.case=TRUE))
}
coef.synthetic_check <- function(object,...) object$b
vcov.synthetic_check <- function(object,...) object$v
valid <- structure(list(b=c(x=.2),v=matrix(.03,1,1)),class='synthetic_check')
stopifnot(identical(checked_model(valid),valid))
expect_error(checked_model({warning('Ran out of iterations and did not converge');valid}),'did not converge')
expect_error(checked_model({warning('Loglik converged before variable 1; coefficient may be infinite');valid}),'infinite')
expect_error(checked_model({warning('Iteration limit exceeded');valid}),'iteration limit')
bad <- valid;bad$b[1] <- Inf
expect_error(checked_model(bad),'non-finite')
bad <- valid;bad$v[1,1] <- NaN
expect_error(checked_model(bad),'non-finite')
bad <- valid;bad$converged <- FALSE
expect_error(checked_model(bad),'converge')
bad <- valid;bad$fail <- 'numerical failure'
expect_error(checked_model(bad),'failure')
# Unrelated warnings must remain observable to an outer caller.
seen <- character()
kept <- withCallingHandlers(checked_model({warning('Benign diagnostic notice');valid}),
  warning=function(w){seen <<- c(seen,conditionMessage(w));invokeRestart('muffleWarning')})
stopifnot(identical(kept,valid),identical(seen,'Benign diagnostic notice'))
# Monotone likelihood produces a finite numerical coefficient and covariance,
# but no finite maximum-likelihood coefficient. The warning must reject the fit.
d <- data.frame(time=1:40,event=rep(c(1,0),each=20),z=rep(c(1,0),each=20),
  weight=1,psu=rep(1:20,each=2),stratum=rep(1:10,each=4))
des <- svydesign(ids=~psu,strata=~stratum,weights=~weight,data=d)
expect_error(checked_model(svycoxph(Surv(time,event)~z,design=des,method='efron')),'coefficient may be infinite')
# A routine estimable fit must retain its exact coefficients and covariance.
d$event <- rep(c(1,0,1,1),10);d$z <- rep(c(0,1),20)
des <- svydesign(ids=~psu,strata=~stratum,weights=~weight,data=d)
f <- svycoxph(Surv(time,event)~z,design=des,method='efron')
g <- checked_model(svycoxph(Surv(time,event)~z,design=des,method='efron'))
stopifnot(identical(coef(f),coef(g)),identical(vcov(f),vcov(g)))
cat('PASS: checked fitting rejects non-estimability and numerical failure, preserves valid fits and unrelated warnings\n')
