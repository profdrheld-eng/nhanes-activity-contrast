#!/usr/bin/env Rscript
arg <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))])
HERE <- dirname(normalizePath(arg))
source(file.path(HERE,'model_common.R'))
targets <- file.path(OUT,c('age_sex_overlay.csv','age_sex_overlay_models.rds'))
if(any(file.exists(targets))) stop('Refusing to overwrite completed or partial stage outputs')
curves <- list(); models <- list()
original <- read.csv(file.path(OUT,'age_curves.csv'))
for (era_name in c('hip','wrist')) {
 des <- subset(design(d[d$era==era_name,]),age_domain)
 common <- checked_model(svyglm(Delta ~ Level + age10 + I(age10^2) + sex + race4 + education3 + PIR + cycle,design=des))
 fit <- checked_model(svyglm(Delta ~ Level + (age10 + I(age10^2))*sex + race4 + education3 + PIR + cycle,design=des))
 stopifnot(nobs(common)==nobs(fit),all(is.finite(coef(fit))))
 baseline <- original[original$model==paste0(era_name,'__age__primary'),]
 for (age in 20:80) {
  z<-(age-50)/10
  expected<-coef(common)['age10']*z+coef(common)['I(age10^2)']*z^2
  stopifnot(abs(expected-baseline$estimate[baseline$age==age])<1e-9)
  for(s in c('Male','Female')) {
   L<-setNames(rep(0,length(coef(fit))),names(coef(fit)))
   L['age10']<-z;L['I(age10^2)']<-z^2
   if(s=='Female'){L['age10:sexfemale']<-z;L['I(age10^2):sexfemale']<-z^2}
   stopifnot(length(L)==length(coef(fit)))
   est<-sum(L*coef(fit));se<-sqrt(max(0,drop(t(L)%*%vcov(fit)%*%L)))
   crit<-qt(.975,degf(des))
   curves[[length(curves)+1]]<-data.frame(era=era_name,sex=s,age=age,estimate=est,ci_low=est-crit*se,ci_high=est+crit*se,n=nobs(fit),design_df=degf(des))
  }
 }
 models[[era_name]]<-fit
}
write.csv(do.call(rbind,curves),file.path(OUT,'age_sex_overlay.csv'),row.names=FALSE)
saveRDS(models,file.path(OUT,'age_sex_overlay_models.rds'))
message('PASS: common curves reproduced at all 122 ages; matched samples; 244 sex-specific predictions exported')
