# Retain the full survey base, then subset the principal domain after splitting.
# Artificial times for excluded records only preserve survey-design information;
# these records never enter any risk set or model fit.
reporting_time_split <- function(x, cut=5) {
  stopifnot(length(cut)==1L,is.finite(cut),cut>0,
    !anyNA(x$primary_complete_case))
  keep <- x$primary_complete_case
  stopifnot(all(is.finite(x$time[keep])),all(x$time[keep]>0),
    all(x$death[keep] %in% c(0,1)))
  x$record_id <- seq_len(nrow(x))
  x$time[!keep] <- 1
  x$death[!keep] <- 0
  Surv <- survival::Surv
  y <- survival::survSplit(Surv(time,death)~.,data=x,
    cut=cut,episode='time_period')
  y$late_delta <- y$Delta*as.integer(y$tstart>=cut)
  y
}
