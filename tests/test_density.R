#!/usr/bin/env Rscript
# Exercise only the binning function, without loading study data or plotting.
arg <- sub('^--file=','',commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))])
root <- dirname(dirname(normalizePath(arg)))
expressions <- parse(file.path(root,'R','plot_density.R'))
selected <- Filter(function(e)is.call(e) && identical(e[[1]],as.name('<-')) &&
  identical(e[[2]],as.name('weighted_bins')),as.list(expressions))
stopifnot(length(selected)==1)
eval(selected[[1]])
# Gaps on both axes must retain their true coordinates, including empty bins.
edges <- 0:5
x <- c(1.5,3.5,1.5);y <- c(.5,2.5,2.5);w <- c(2,7,3)
expected <- matrix(0,5,5);expected[2,1] <- 2;expected[4,3] <- 7;expected[2,3] <- 3
actual <- weighted_bins(x,y,w,edges)
stopifnot(identical(actual,expected),sum(actual)==sum(w))
# Intervals are left-closed; outer endpoints and off-scale values enter edge bins.
x <- c(-1,0,1,5,6);y <- c(6,5,1,0,-1);w <- c(1,2,4,8,16)
expected <- matrix(0,5,5);expected[1,5] <- 3;expected[2,2] <- 4;expected[5,1] <- 24
actual <- weighted_bins(x,y,w,edges)
stopifnot(identical(actual,expected),sum(actual)==sum(w))
# All-in-one-cell and zero-weight observations must not compress either axis.
actual <- weighted_bins(c(3.1,3.2),c(2.1,2.2),c(0,9),edges)
expected <- matrix(0,5,5);expected[4,3] <- 9
stopifnot(identical(actual,expected))
cat('PASS: full-grid density placement, gaps, unequal weights, boundaries and mass conservation\n')
