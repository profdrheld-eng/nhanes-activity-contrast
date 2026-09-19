#!/usr/bin/env Rscript
# Wrist response-weighting balance relative to the complete-predictor target
# population. Standardizers use target SDs, not pooled two-group SDs.
balance_row <- function(x,base,adjusted,respondent,variable,level='continuous') {
  stopifnot(length(x)==length(base),length(x)==length(adjusted),
            length(x)==length(respondent),is.logical(respondent),!anyNA(respondent),
            all(is.finite(x)),all(is.finite(base)&base>0),any(respondent),
            all(is.finite(adjusted[respondent])&adjusted[respondent]>0))
  mean_target <- sum(x*base)/sum(base)
  standardizer <- if(level=='continuous') sqrt(sum(base*(x-mean_target)^2)/sum(base)) else
    sqrt(mean_target*(1-mean_target))
  before <- sum(x[respondent]*base[respondent])/sum(base[respondent])
  after <- sum(x[respondent]*adjusted[respondent])/sum(adjusted[respondent])
  data.frame(variable=variable,level=level,target_value=mean_target,
             respondent_before=before,respondent_after=after,
             smd_before=if(standardizer>0)(before-mean_target)/standardizer else NA_real_,
             smd_after=if(standardizer>0)(after-mean_target)/standardizer else NA_real_)
}

if(sys.nframe()==0) {
  args <- commandArgs(trailingOnly=TRUE)
  if(length(args)!=1)stop('Supply exactly one local work directory')
  out <- file.path(normalizePath(args[1],mustWork=TRUE),'response/wrist')
  if(any(file.exists(file.path(out,c('balance.csv','balance_summary.csv')))))stop('Refusing to overwrite balance outputs')
  # The freshly fitted survey model retains its exact response-model domain.
  fit <- readRDS(file.path(out,'response_model.rds'))
  d <- fit$survey.design$variables
  w <- read.csv(gzfile(file.path(out,'pax_response_weights.csv.gz')))
  stopifnot(!anyNA(d$SEQN),!anyDuplicated(d$SEQN),!anyNA(w$SEQN),!anyDuplicated(w$SEQN),
            all(d$response_model_domain),all(w$SEQN %in% d$SEQN))
  respondent <- d$PAX_response==1
  adjusted <- w$combined_response_weight[match(d$SEQN,w$SEQN)]
  stopifnot(sum(respondent)==nrow(w),all(is.finite(adjusted)==respondent))
  rows <- lapply(c('age','PIR','BMI'),function(v)balance_row(d[[v]],d$base_weight,adjusted,respondent,v))
  for(v in c('sex','race4','education3','smoking3','cycle')) {
    for(level in levels(d[[v]]))rows[[length(rows)+1]] <- balance_row(
      as.numeric(d[[v]]==level),d$base_weight,adjusted,respondent,v,level)
  }
  balance <- do.call(rbind,rows)
  balance$absolute_smd_reduction <- abs(balance$smd_before)-abs(balance$smd_after)
  write.csv(balance,file.path(out,'balance.csv'),row.names=FALSE)
  write.csv(data.frame(target_n=nrow(d),respondent_n=sum(respondent),
    maximum_absolute_smd_before=max(abs(balance$smd_before),na.rm=TRUE),
    maximum_absolute_smd_after=max(abs(balance$smd_after),na.rm=TRUE)),
    file.path(out,'balance_summary.csv'),row.names=FALSE)
}
