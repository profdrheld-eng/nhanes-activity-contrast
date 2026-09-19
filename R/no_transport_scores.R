# One hip cycle at a time; retain the original reference and analysis domains.
# Sources score_helpers.R; requires wm(). Never infer a new eligibility mask.
no_transport_scores <- function(x) {
  ref <- x$reference_domain
  stopifnot(is.logical(ref),!anyNA(ref),sum(ref)>1)
  components <- c('SR_home_30d','SR_moderate_30d','SR_vigorous_30d')
  raw <- as.matrix(x[ref,components])
  w <- x$WTMEC2YR[ref]
  stopifnot(all(is.finite(raw)),all(raw>=0),all(is.finite(w)),all(w>0),
            all(is.finite(x$D[ref])))
  x$SR_no_transport <- (x$SR_home_30d+x$SR_moderate_30d+2*x$SR_vigorous_30d)/30
  x$S_no_transport <- x$Delta_no_transport <- x$Level_no_transport <- NA_real_
  s <- midrank(log1p(x$SR_no_transport[ref]),w)
  delta <- stand(s-x$D[ref],w);level <- stand((s+x$D[ref])/2,w)
  stopifnot(all(is.finite(delta)),all(is.finite(level)))
  x$S_no_transport[ref] <- s
  x$Delta_no_transport[ref] <- delta
  x$Level_no_transport[ref] <- level
  x
}
