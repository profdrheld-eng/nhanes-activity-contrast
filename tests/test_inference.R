#!/usr/bin/env Rscript
# Independent synthetic contrast checks; no study participants or model fits.
suppressPackageStartupMessages(library(survey))
arg <- sub('^--file=','',commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))])
root <- dirname(dirname(normalizePath(arg)))
source(file.path(root,'R','inference.R'))
coef.synthetic_fit <- function(object,...) object$b
vcov.synthetic_fit <- function(object,...) object$v
f <- structure(list(b=c(hip=-0.3,wrist=0.1),v=matrix(c(.04,.005,.005,.02),2,
  dimnames=list(c('hip','wrist'),c('hip','wrist'))),df.residual=2),class='synthetic_fit')
x <- data.frame(psu=1:8,stratum=rep(1:4,each=2),weight=1,death=c(1,0,0,0,0,0,0,0))
des <- svydesign(ids=~psu,strata=~stratum,weights=~weight,data=x)
r <- extract(f,des,c(hip=-1,wrist=1),'test','Wrist/Hip slope ratio')
main <- r[r$inference=='full_design_t',]
# The contrast is 0.4; its variance is .04+.02-2*.005=.05.
stopifnot(abs(main$coefficient-.4)<1e-14,abs(main$se-sqrt(.05))<1e-14)
stopifnot(abs(main$estimate-1.49182469764127)<1e-13,main$design_df==4,main$n==8,main$deaths==1)
# t(.975,4)=2.7764451051977987, retained independently as a tabulated constant.
stopifnot(abs(main$ci_low-exp(.4-2.7764451051977987*sqrt(.05)))<1e-12)
stopifnot(abs(main$ci_high-exp(.4+2.7764451051977987*sqrt(.05)))<1e-12)
stopifnot(identical(r$denominator_df,c(4,2,Inf)))
a <- extract(f,des,c(hip=-1,wrist=1),'test','Linear difference',FALSE)
stopifnot(all(abs(a$estimate-.4)<1e-14))
j <- joint(f,des,c('hip','wrist'),'test','two coefficients')
# Manual 2x2 inverse: b' V^-1 b = (.02*.09 + 2*.005*.03 + .04*.01)/(.04*.02-.005^2).
expected <- (.02*.09 + 2*.005*.03 + .04*.01)/(.04*.02-.005^2)
stopifnot(abs(j$F-expected/2)<1e-13,j$numerator_df==2,j$design_df==4)
cat('PASS: synthetic coefficient contrasts, exponentiation, covariance, df, intervals and joint Wald test\n')
# Ambiguous/malformed contrast requests must fail, never silently overwrite a term.
reject <- function(expr) stopifnot(inherits(tryCatch(force(expr),error=identity),'error'))
reject(extract(f,des,c(hip=1,hip=-1),'test','duplicate'))
reject(extract(f,des,c(1,-1),'test','unnamed'))
reject(extract(f,des,numeric(0),'test','empty'))
reject(extract(f,des,c(missing=1),'test','unknown'))
reject(extract(f,des,c(hip=NA_real_),'test','nonfinite'))
reject(joint(f,des,c('hip','hip'),'test','duplicate'))
# Named covariance axes may be permuted; inference must follow coefficient names.
permuted <- f;permuted$v <- f$v[c('wrist','hip'),c('wrist','hip')]
stopifnot(identical(extract(permuted,des,c(hip=1),'test','single'),
                    extract(f,des,c(hip=1),'test','single')))
cat('PASS: malformed contrasts rejected and covariance aligned by name\n')
