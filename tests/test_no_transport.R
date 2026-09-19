arg <- sub('^--file=','',commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))])
root <- dirname(dirname(normalizePath(arg)))
wm <- function(x,w)sum(x*w)/sum(w)
source(file.path(root,'R','score_helpers.R'))
source(file.path(root,'R','no_transport_scores.R'))
x <- data.frame(SEQN=1:5,reference_domain=c(TRUE,TRUE,TRUE,TRUE,FALSE),WTMEC2YR=c(1,2,3,4,1000),
 SR_home_30d=c(30,60,30,90,99999),SR_moderate_30d=0,SR_vigorous_30d=c(15,0,15,0,99999),
 SR_transport_30d=c(900,0,0,60,NA),D=c(-1,0,1,2,NA),primary_complete_case=c(TRUE,TRUE,FALSE,TRUE,FALSE))
a <- no_transport_scores(x)
stopifnot(identical(a$SEQN,x$SEQN),identical(a$primary_complete_case,x$primary_complete_case),
 identical(a$reference_domain,x$reference_domain),identical(a$D,x$D))
stopifnot(identical(a$SR_no_transport[1:4],c(2,2,2,3)),is.na(a$S_no_transport[5]))
# First three tied observations carry 6/10 of reference weight: midpoint .3;
# last observation carries 4/10: midpoint .8. Excluded extreme row has no effect.
stopifnot(max(abs(a$S_no_transport[1:4]-qnorm(c(.3,.3,.3,.8))))<1e-14)
for(col in c('Delta_no_transport','Level_no_transport')) {
 z<-a[[col]][1:4];w<-x$WTMEC2YR[1:4]
 stopifnot(abs(wm(z,w))<1e-14,abs(wm(z^2,w)-1)<1e-14)
}
y<-x;y$SR_transport_30d<-c(0,9000,1000,NA,0)
cols<-c("SR_no_transport","S_no_transport","Delta_no_transport","Level_no_transport")
stopifnot(identical(no_transport_scores(y)[cols],a[cols])) # transport data cannot enter arithmetic
reject<-function(expr)stopifnot(inherits(tryCatch(force(expr),error=identity),'error'))
y<-x;y$SR_home_30d[2]<-NA;reject(no_transport_scores(y))
y<-x;y$WTMEC2YR[2]<-0;reject(no_transport_scores(y))
cat('PASS: no-transport arithmetic, weighted ties, fixed domain, sample preservation and invalid inputs\n')
