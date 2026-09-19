#!/usr/bin/env Rscript
fit_script <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))])
source(file.path(dirname(normalizePath(fit_script)),'fit_checks.R'))
# Diagnostics deliberately distinguish survey step-time checks from the
# conventional weighted Schoenfeld and cluster-robust exact log-time checks.
options(stringsAsFactors=FALSE,survey.lonely.psu='adjust')
suppressPackageStartupMessages({library(survey);library(survival)})
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1) stop('Supply exactly one local work directory')
WORK <- normalizePath(args[1],mustWork=TRUE)
OUT <- file.path(WORK,'models')
targets <- file.path(OUT,c('ph_diagnostics.csv','influence.csv','diagnostic_support.csv','time_varying.csv'))
if(any(file.exists(targets))) stop('Refusing to overwrite completed or partial stage outputs')
fits <- readRDS(file.path(OUT,'primary_fits.rds'))
ph<-list();influence<-list();checks<-list();time<-list()
for(era in c('hip','wrist')) for(type in c('common','interaction')) {
  id<-paste(era,type,'primary',sep='__');obj<-fits[[id]];d<-obj$domain$variables
  d$analysis_weight<-d$weight/mean(d$weight);ff<-formula(obj$fit)
  cf<-checked_model(coxph(ff,data=d,weights=analysis_weight,ties='efron',robust=TRUE,cluster=design_psu,x=TRUE,y=TRUE,model=TRUE))
  stopifnot(max(abs(coef(cf)-coef(obj$fit)))<1e-8)
  z<-cox.zph(cf,transform='km')$table
  ph[[id]]<-data.frame(model=id,term=rownames(z),z,row.names=NULL)
  db<-residuals(cf,type='dfbeta',weighted=TRUE)
  influence[[id]]<-data.frame(model=id,term=names(coef(cf)),max_abs_weighted_dfbeta=apply(abs(db),2,max))
  X<-cf$x
  legacy<-scale(X)*sqrt(d$analysis_weight)
  mu<-colSums(X*d$analysis_weight)/sum(d$analysis_weight)
  centered<-sweep(X,2,mu,'-');s<-sqrt(colSums(centered^2*d$analysis_weight)/sum(d$analysis_weight))
  WX<-sweep(centered,2,s,'/')*sqrt(d$analysis_weight)
  singular<-svd(WX,nu=0,nv=0)$d
  checks[[id]]<-data.frame(model=id,n=nrow(d),deaths=sum(d$death),design_df=degf(obj$domain),
    coefficients_match_svycoxph=max(abs(coef(cf)-coef(obj$fit))),
    legacy_gram_condition_number=kappa(crossprod(legacy),exact=TRUE),
    weighted_matrix_condition_number=max(singular)/min(singular),
    max_abs_weighted_dfbeta=max(abs(db)))
  # Exact step change after five years with survey covariance, not a time transform approximation.
  split<-survSplit(Surv(time,death)~.,data=d,cut=5,episode='period')
  split$late_delta<-split$Delta*as.integer(split$tstart>=5)
  split$late_sex<-as.integer(split$sex=='female')*as.integer(split$tstart>=5)
  split$late_interaction<-split$late_delta*as.integer(split$sex=='female')
  des<-svydesign(ids=~design_psu,strata=~design_stratum,weights=~weight,nest=TRUE,data=split)
  tf<-update(ff,Surv(tstart,time,death)~. + late_delta + late_sex)
  if(type=='interaction')tf<-update(tf,.~.+late_interaction)
  f<-checked_model(svycoxph(tf,design=des,method='efron'))
  for(term in c('late_delta','late_sex',if(type=='interaction')'late_interaction')){
    b<-coef(f)[term];se<-sqrt(vcov(f)[term,term]);dd<-degf(des)
    time[[paste(id,term)]]<-data.frame(model=id,method='survey step change after 5 years',term=term,coefficient=b,se=se,
      ratio=exp(b),ci_low=exp(b-qt(.975,dd)*se),ci_high=exp(b+qt(.975,dd)*se),p=2*pt(-abs(b/se),dd),denominator_df=dd,
      n=nrow(d),deaths=sum(d$death),n_after5=sum(d$time>5),deaths_after5=sum(d$death[d$time>5]))
  }
  # Exact log-time diagnostic uses cluster-robust, not full survey, inference.
  # This reproduces the historical diagnostic approach and is labelled accordingly.
  d$sex_female<-as.integer(d$sex=='female')
  tf<-update(ff,.~.+tt(Delta)+tt(sex_female))
  f<-checked_model(coxph(tf,data=d,weights=analysis_weight,ties='efron',robust=TRUE,cluster=design_psu,
           tt=function(x,t,...)x*log(pmax(t,1/365.25))))
  for(term in c('tt(Delta)','tt(sex_female)')){
    b<-coef(f)[term];se<-sqrt(vcov(f)[term,term])
    time[[paste(id,term)]]<-data.frame(model=id,method='weighted cluster-robust exact log time (not full survey)',term=term,coefficient=b,se=se,
      ratio=exp(b),ci_low=exp(b-qnorm(.975)*se),ci_high=exp(b+qnorm(.975)*se),p=2*pnorm(-abs(b/se)),denominator_df=Inf,
      n=nrow(d),deaths=sum(d$death),n_after5=sum(d$time>5),deaths_after5=sum(d$death[d$time>5]))
  }
  rm(f);gc();message('Diagnostics: ',id)
  write.csv(do.call(rbind,ph),file.path(OUT,'ph_diagnostics.csv'),row.names=FALSE)
  write.csv(do.call(rbind,influence),file.path(OUT,'influence.csv'),row.names=FALSE)
  write.csv(do.call(rbind,checks),file.path(OUT,'diagnostic_support.csv'),row.names=FALSE)
  write.csv(do.call(rbind,time),file.path(OUT,'time_varying.csv'),row.names=FALSE)
}
