arg <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))])
source(file.path(dirname(normalizePath(arg)),'../R/reporting_time_split.R'))
x <- data.frame(time=c(2,5,7,8,NA),death=c(1,1,1,0,NA),
  primary_complete_case=c(TRUE,TRUE,TRUE,TRUE,FALSE),Delta=c(1,2,3,4,NA))
y <- reporting_time_split(x)
z <- y[y$primary_complete_case,]
stopifnot(nrow(z)==6,sum(z$death)==3,
  sum(z$time-z$tstart)==22,
  all(z$late_delta[z$tstart==0]==0),
  identical(z$record_id[z$tstart==5],c(3L,4L)),
  sum(z$death[z$time<=5])==2,
  sum(z$death[z$tstart>=5])==1,
  nrow(y[!y$primary_complete_case,])==1)
bad <- x;bad$time[1] <- 0
stopifnot(inherits(try(reporting_time_split(bad),silent=TRUE),'try-error'))
cat('Time-split boundary, person-time, events and excluded-record checks passed\n')
