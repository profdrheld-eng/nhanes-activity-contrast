# Evaluate a model fit and reject numerical failures before exporting inference.
# Cox fits can return finite coefficients and covariance despite monotone
# likelihood, so checking only finite output and fit$fail is insufficient.
checked_model <- function(expr) {
  fit <- withCallingHandlers(force(expr),warning=function(w) {
    message <- conditionMessage(w)
    fatal <- grepl(paste(c('did not converge','failed to converge',
      'convergence fail','coefficient.*infinite','infinite.*coefficient',
      'iteration limit','ran out of iterations'),collapse='|'),message,
      ignore.case=TRUE)
    if(fatal)stop(paste('Model fit rejected:',message),call.=FALSE)
    # Do not muffle other warnings: retain the model's ordinary warning behavior.
  })
  if(!is.null(fit$fail))stop(paste('Model fit failure:',fit$fail),call.=FALSE)
  if(!is.null(fit$converged) && !isTRUE(fit$converged))
    stop('Model did not converge',call.=FALSE)
  b <- coef(fit);v <- vcov(fit)
  if(!length(b) || any(!is.finite(b)) || any(!is.finite(v)))
    stop('Model has missing or non-finite coefficients/covariance',call.=FALSE)
  fit
}
