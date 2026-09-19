#!/usr/bin/env Rscript
arg <- sub('^--file=','',commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))])
root <- dirname(dirname(normalizePath(arg)))
source(file.path(root,'R/response_balance.R'))
r <- balance_row(c(0,2),c(1,3),c(3,1),c(TRUE,TRUE),'example')
stopifnot(abs(r$target_value-1.5)<1e-14,abs(r$respondent_after-.5)<1e-14,
          abs(r$smd_before)<1e-14,abs(r$smd_after+1/sqrt(.75))<1e-14)
b <- balance_row(c(0,1),c(1,3),c(3,1),c(TRUE,TRUE),'example','yes')
stopifnot(abs(b$smd_after+.5/sqrt(.1875))<1e-14)
constant <- balance_row(c(1,1),c(1,3),c(3,1),c(TRUE,TRUE),'constant')
stopifnot(is.na(constant$smd_before),is.na(constant$smd_after))
cat('PASS: response balance uses target population SD and target Bernoulli SD\n')
