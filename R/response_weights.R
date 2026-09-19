#!/usr/bin/env Rscript
fit_script <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))])
source(file.path(dirname(normalizePath(fit_script)),'fit_checks.R'))
# Refit stabilized response weights from the newly constructed cohorts.
# No weight trimming. Mortality outcome and follow-up duration are not predictors;
# availability of valid linkage/follow-up remains part of the eligibility gate.
# Subsequent outcome models condition on these fitted weights (no full-pipeline SE).
options(stringsAsFactors=FALSE, survey.lonely.psu='adjust')
suppressPackageStartupMessages(library(survey))
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<1 || length(args)>2) stop('Usage: response_weights.R WORK_DIR [hip|wrist]')
WORK <- normalizePath(args[1],mustWork=TRUE)
eras <- if(length(args)==2) args[2] else c('hip','wrist')
if(!all(eras %in% c('hip','wrist'))) stop('Unknown response-weighting era')
flag <- function(x) !is.na(x) & tolower(as.character(x)) %in% c('true','1','yes')
for(era in eras){
  cycles <- if(era=='hip') '2005-2006' else c('2011-2012','2013-2014')
  d <- do.call(rbind,lapply(cycles,function(cy) read.csv(gzfile(file.path(WORK,cy,'cohort','analytic_cohort.csv.gz')))))
  if(anyNA(d$SEQN) || anyDuplicated(d$SEQN)) stop('Participant keys must be unique')
  for(v in c('adult_domain','mortality_followup_valid','SR_complete','PAX_primary_eligible')) d[[v]] <- flag(d[[v]])
  d$sex <- factor(d$sex,levels=c('male','female'))
  d$race4 <- factor(d$race_ethnicity4,levels=c('Non-Hispanic White','Hispanic','Non-Hispanic Black','Other'))
  d$education3 <- factor(d$education3,levels=c('more than high school','high school/GED','under high school'))
  d$smoking3 <- factor(d$smoking3,levels=c('never','former','current'))
  predictors <- c('age','sex','race4','education3','PIR','BMI','smoking3')
  d$cycle <- factor(d$cycle,levels=cycles)
  if(era=='wrist') predictors <- c(predictors,'cycle')
  d$response_path <- d$adult_domain & d$mortality_followup_valid & d$SR_complete
  d$response_model_domain <- d$response_path & complete.cases(d[,predictors])
  d$PAX_response <- as.integer(d$PAX_primary_eligible)
  d$base_weight <- d$WTMEC2YR/length(cycles)
  if(era=='hip'){
    des <- svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~base_weight,nest=TRUE,data=d)
  }else{
    d$design_stratum <- interaction(d$cycle,d$SDMVSTRA,drop=TRUE)
    d$design_psu <- interaction(d$cycle,d$SDMVSTRA,d$SDMVPSU,drop=TRUE)
    des <- svydesign(ids=~design_psu,strata=~design_stratum,weights=~base_weight,nest=TRUE,data=d)
  }
  model <- checked_model(svyglm(reformulate(predictors,response='PAX_response'),
    design=subset(des,response_model_domain),family=quasibinomial(link='logit')))
  if(!isTRUE(model$converged) || any(!is.finite(coef(model)))) stop('Response model failed to converge')
  d$response_probability <- NA_real_
  d$response_probability[d$response_model_domain] <- as.numeric(predict(model,type='response'))
  p <- d$response_probability[d$response_model_domain]
  if(any(!is.finite(p) | p<=0 | p>=1)) stop('Response probability outside (0,1)')
  domain <- d$response_model_domain
  rate <- sum(d$base_weight[domain]*d$PAX_response[domain])/sum(d$base_weight[domain])
  respondent <- domain & d$PAX_response==1
  d$response_multiplier <- NA_real_
  d$response_multiplier[respondent] <- rate/d$response_probability[respondent]
  d$combined_response_weight <- d$base_weight*d$response_multiplier
  out <- file.path(WORK,'response',era)
  if(dir.exists(out)) stop('Refusing to overwrite response-weight outputs')
  dir.create(out,recursive=TRUE)
  # Compatibility names describe the same combined examination-response weight.
  if(era=='hip'){
    x <- d[domain,c('SEQN','PAX_response','response_probability','response_multiplier','base_weight','combined_response_weight','SDMVSTRA','SDMVPSU')]
    names(x)[names(x)=='base_weight'] <- 'mec_weight'
    names(x)[names(x)=='combined_response_weight'] <- 'combined_weight'
  }else{
    x <- d[respondent,c('SEQN','cycle','response_probability','response_multiplier','base_weight','combined_response_weight')]
    names(x)[names(x)=='base_weight'] <- 'pooled_weight_era'
  }
  con <- gzfile(file.path(out,'pax_response_weights.csv.gz'),'wb')
  write.csv(x,con,row.names=FALSE);close(con)
  write.csv(data.frame(term=names(coef(model)),coefficient=unname(coef(model)),
    design_based_se=sqrt(diag(vcov(model)))),file.path(out,'model_coefficients.csv'),row.names=FALSE)
  write.csv(data.frame(era=era,response_model_n=sum(domain),respondents=sum(respondent),
    design_df=degf(subset(des,response_model_domain)),response_rate=rate,
    min_probability=min(p),max_probability=max(p),min_multiplier=min(d$response_multiplier[respondent]),
    max_multiplier=max(d$response_multiplier[respondent]),weights_trimmed=FALSE),
    file.path(out,'diagnostics.csv'),row.names=FALSE)
  saveRDS(model,file.path(out,'response_model.rds'))
  writeLines(capture.output(sessionInfo()),file.path(out,'session_info.txt'))
  message('Response weights reconstructed: ',era)
}
