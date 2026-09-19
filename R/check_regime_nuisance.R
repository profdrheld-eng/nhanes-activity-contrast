#!/usr/bin/env Rscript
arg <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))])
HERE <- dirname(normalizePath(arg))
source(file.path(HERE,'model_common.R'))
targets <- file.path(OUT,c('regime_nuisance_results.csv','regime_nuisance_tests.csv'))
if(any(file.exists(targets))) stop('Refusing to overwrite completed or partial stage outputs')
rr<-list();jj<-list()
for(h in c(Inf,5,8))for(spec in c('sex_era','all_nuisance_era')){
  x<-d;x$time<-pmin(x$follow_up_years,h);x$death<-as.integer(x$MORTSTAT==1 & x$follow_up_years<=h)
  des<-subset(design(x),primary_complete_case)
  rhs<-if(spec=='sex_era')'Delta+sex+Delta:era+sex:era+Level+age+race4+education3+PIR+smoking3+strata(cycle)' else
    '(Delta+sex+Level+age+race4+education3+PIR+smoking3)*era-era+strata(cycle)'
  f<-checked_model(svycoxph(as.formula(paste('Surv(time,death)~',rhs)),design=des,method='efron'))
  term<-names(coef(f))[grepl('Delta:erawrist|erawrist:Delta',names(coef(f)))];stopifnot(length(term)==1)
  id<-paste(spec,h,sep='__');rr[[id]]<-extract(f,des,setNames(1,term),id,'Wrist/Hip slope ratio')
  jj[[id]]<-joint(f,des,term,id,'Regime heterogeneity')
}
write.csv(do.call(rbind,rr),file.path(OUT,'regime_nuisance_results.csv'),row.names=FALSE)
write.csv(do.call(rbind,jj),file.path(OUT,'regime_nuisance_tests.csv'),row.names=FALSE)
