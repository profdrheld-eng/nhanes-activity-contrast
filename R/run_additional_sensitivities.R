#!/usr/bin/env Rscript
arg <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))])
HERE <- dirname(normalizePath(arg))
source(file.path(HERE,'model_common.R'))
targets <- file.path(OUT,c('additional_sensitivity_results.csv','additional_model_support.csv'))
if(any(file.exists(targets))) stop('Refusing to overwrite completed or partial stage outputs')
source(file.path(HERE,'score_helpers.R'))
for(name in c('2003-2004','2005-2006','hip','2011-2012','2013-2014','wrist')){
  original<-groups[[name]];iship<-all(original$era=='hip')
  variants<-if(iship)c('valid_days_3','valid_days_5','nonwear_90','calibration_1_or_2') else c('intensity','one_day')
  for(variant in variants){
    x<-original
    if(iship){
      for(cy in unique(x$cycle)){
        z<-x$cycle==cy
        eligible<-flag(x[[paste0('PAX_eligible_participant__',variant)]])
        metric<-x[[paste0('PAX_log1p_counts_per_wear_minute__',variant)]]
        ref<-z & flag(x$adult_domain)&flag(x$mortality_followup_valid)&flag(x$SR_complete)&eligible&is.finite(metric)
        w<-x$weight[ref];s<-midrank(log1p(x$SR_day[ref]),w);dv<-midrank(metric[ref],w)
        x$Delta[z]<-NA_real_;x$Level[z]<-NA_real_
        x$Delta[ref]<-stand(s-dv,w);x$Level[ref]<-stand((s+dv)/2,w)
      }
      x$primary_complete_case<-complete.cases(x[,c('Delta','Level','age','sex','race4','education3','PIR','smoking3')])
    } else {
      x$Delta<-x[[paste0('Delta_',variant)]];x$Level<-x[[paste0('Level_',variant)]]
      x$primary_complete_case<-flag(x[[paste0(variant,'_complete_case')]])
    }
    for(inter in c(FALSE,TRUE))fit_model(x,paste(name,if(inter)'interaction' else 'common',variant,sep='__'),inter)
  }
}
# Use stabilized response weights refitted in the explicit preceding stage.
# Outcome-model variances condition on these fitted weights, as in the paper.
for(era in c('hip','wrist')){
  if(era=='hip'){
    x<-groups[['2005-2006']]
    wf<-file.path(WORK,'response/hip/pax_response_weights.csv.gz')
    w<-read.csv(gzfile(wf))
  }else{
    x<-groups[['wrist']];wf<-file.path(WORK,'response/wrist/pax_response_weights.csv.gz')
    w<-read.csv(gzfile(wf))
  }
  # Both original files contain the combined stabilized examination weight.
  weight_col<-if(era=='hip')'combined_weight' else 'combined_response_weight'
  stopifnot(weight_col %in% names(w))
  if(length(weight_col)!=1)stop(paste('Identify exact response weight column in',wf))
  neww<-w[[weight_col]][match(x$SEQN,w$SEQN)]
  x$primary_complete_case<-x$primary_complete_case & is.finite(neww)&neww>0
  # Nonrespondents remain outside the analysis domain; retain valid base weights there.
  x$weight[x$primary_complete_case]<-neww[x$primary_complete_case]
  for(inter in c(FALSE,TRUE))fit_model(x,paste(if(era=='hip')'2005-2006' else 'wrist',if(inter)'interaction' else 'common','response_weighted',sep='__'),inter)
}
for(era in c('hip','wrist')){
  x<-groups[[era]]
  fit_model(x,paste(era,'age_quadratic',sep='__'),extra='+ I(age^2)')
  des<-subset(design(x),primary_complete_case)
  f<-checked_model(svycoxph(Surv(time,death)~Delta+Level+age+race4+education3+PIR+smoking3+strata(cycle,sex),design=des,method='efron'))
  id<-paste(era,'sex_stratified',sep='__')
  results[[id]]<-extract(f,des,c(Delta=1),id,'Common contrast')
  support[[id]]<-data.frame(model=id,n=nrow(des$variables),deaths=sum(des$variables$death),strata=length(unique(des$variables$design_stratum)),
    PSUs=length(unique(des$variables$design_psu)),design_df=degf(des),parameters=length(coef(f)),residual_df=resdf(f),formula=paste(deparse(formula(f)),collapse=' '))
}
write.csv(do.call(rbind,results),file.path(OUT,'additional_sensitivity_results.csv'),row.names=FALSE)
write.csv(do.call(rbind,support),file.path(OUT,'additional_model_support.csv'),row.names=FALSE)
message('Additional sensitivities complete')
