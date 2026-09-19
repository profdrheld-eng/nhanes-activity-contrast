#!/usr/bin/env Rscript
# Additional sex-interaction disease exclusions retained in Tables S15 and S29.
# Other retained sensitivities reuse the identically specified fresh main fits.
arg <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))])
HERE <- dirname(normalizePath(arg))
source(file.path(HERE,'model_common.R'))
targets <- file.path(OUT,c('retained_sensitivity_results.csv','retained_model_support.csv'))
if(any(file.exists(targets))) stop('Refusing to overwrite completed or partial stage outputs')
for(name in c('2005-2006','wrist')) {
  x <- groups[[name]]
  fit_model(x,paste(name,'interaction','exclude_CVD_cancer',sep='__'),
            interaction=TRUE,keep=x$cvd_cancer=='no')
}
# Table S15 retains an earlier analysis that reranked a 12-significant-digit
# cohort export. Reproduce that documented serialization here from fresh data.
# Current Table S12 uses full-precision inputs and remains a separate analysis.
source(file.path(HERE,'score_helpers.R'))
x0 <- groups[['2005-2006']]
for(v in names(x0)[vapply(x0,is.numeric,logical(1))]) {
  finite <- is.finite(x0[[v]])
  x0[[v]][finite] <- as.numeric(sprintf('%.12g',x0[[v]][finite]))
}
x0$weight <- x0$WTMEC2YR
for(variant in c('valid_days_3','valid_days_5','calibration_1_or_2','nonwear_90')) {
  x <- x0
  metric <- x[[paste0('PAX_log1p_counts_per_wear_minute__',variant)]]
  eligible <- flag(x[[paste0('PAX_eligible_participant__',variant)]])
  ref <- x$adult_domain & flag(x$mortality_followup_valid) & flag(x$SR_complete) &
    eligible & is.finite(metric) & is.finite(x$SR_day)
  w <- x$weight[ref];s <- midrank(log1p(x$SR_day[ref]),w);dv <- midrank(metric[ref],w)
  x$Delta <- NA_real_;x$Level <- NA_real_
  x$Delta[ref] <- stand(s-dv,w);x$Level[ref] <- stand((s+dv)/2,w)
  x$primary_complete_case <- ref & complete.cases(x[,c('age','sex','race4','education3','PIR','smoking3')])
  fit_model(x,paste('2005-2006','legacy_interaction',variant,sep='__'),interaction=TRUE)
}
write.csv(do.call(rbind,results),file.path(OUT,'retained_sensitivity_results.csv'),row.names=FALSE)
write.csv(do.call(rbind,support),file.path(OUT,'retained_model_support.csv'),row.names=FALSE)
