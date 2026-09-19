#!/usr/bin/env Rscript
arg <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))])
HERE <- dirname(normalizePath(arg))
source(file.path(HERE,'model_common.R'))
if(file.exists(file.path(OUT,'age_support_sensitivity.csv'))) stop('Refusing to overwrite stage outputs')
rows<-list()
for(name in names(groups)) for(rule in c('historical','cap80','under80','harmonized_disease_domain')){
  x<-groups[[name]]
  if(rule=='cap80')x$age10<-(pmin(x$age,80)-50)/10
  if(rule=='under80')x$age_domain<-x$age_domain & x$age<80
  if(rule=='harmonized_disease_domain')x$age_domain<-x$primary_complete_case & is.finite(x$BMI)&!is.na(x$disease)
  des<-subset(design(x),age_domain)
  f<-checked_model(svyglm(as.formula(paste('Delta~Level+age10+I(age10^2)+sex+race4+education3+PIR',if(name%in%c('hip','wrist'))'+cycle' else '')),design=des))
  for(a in if(rule=='under80')75 else c(75,80)){
    z<-(a-50)/10;id<-paste(name,rule,a,sep='__')
    rows[[id]]<-extract(f,des,c(age10=z,'I(age10^2)'=z*z),id,paste('Age',a,'versus 50'),FALSE)
  }
}
write.csv(do.call(rbind,rows),file.path(OUT,'age_support_sensitivity.csv'),row.names=FALSE)
