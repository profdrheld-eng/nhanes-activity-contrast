#!/usr/bin/env Rscript
# Same-case unadjusted, Level-only and principal-adjusted mortality reporting.
arg <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))])
HERE <- dirname(normalizePath(arg))
source(file.path(HERE,'model_common.R'))
report_dir <- file.path(WORK,'reporting')
if(!file.exists(file.path(report_dir,'descriptive_manifest.json')))stop('Run prepare_reporting.py first')
targets <- file.path(report_dir,c('mortality_adjustment.csv','mortality_support.csv','adjustment_fits.rds','session_info.txt'))
if(any(file.exists(targets)))stop('Refusing to overwrite reporting model outputs')
original <- read.csv(file.path(OUT,'mortality_results.csv'))
rows <- list(); info <- list(); models <- list()
for(name in names(groups)) {
  x <- groups[[name]]
  # Keep the full examined survey base before domain subsetting, as in principal models.
  des <- subset(design(x),primary_complete_case)
  expected_n <- sum(x$primary_complete_case)
  for(adjustment in c('unadjusted','level_only','principal')) {
    rhs <- switch(adjustment,unadjusted='Delta',level_only='Delta + Level',
      principal='Delta + sex + Level + age + race4 + education3 + PIR + smoking3')
    if(length(unique(x$cycle))>1)rhs<-paste(rhs,'+ strata(cycle)')
    formula <- as.formula(paste('Surv(time,death) ~',rhs))
    fit <- checked_model(svycoxph(formula,design=des,method='efron',x=TRUE))
    stopifnot(is.null(fit$na.action),fit$n==expected_n,fit$nevent==sum(des$variables$death))
    id<-paste(name,adjustment,sep='__')
    row<-extract(fit,des,c(Delta=1),id,'Signed contrast')
    row$group<-name;row$adjustment<-adjustment;rows[[id]]<-row
    info[[id]]<-data.frame(model=id,n=expected_n,deaths=sum(des$variables$death),
      design_df=degf(des),formula=paste(deparse(formula),collapse=' '))
    models[[id]]<-list(fit=fit,domain=des)
    if(adjustment=='principal') {
      old<-original[original$model==paste(name,'common','primary',sep='__'),]
      stopifnot(nrow(old)==3)
      for(conv in row$inference) {
        a<-row[row$inference==conv,];b<-old[old$inference==conv,]
        for(v in c('coefficient','se','estimate','ci_low','ci_high','p'))stopifnot(abs(a[[v]]-b[[v]])<1e-8)
        stopifnot(a$n==b$n,a$deaths==b$deaths,a$design_df==b$design_df)
      }
    }
  }
  message('Same-case reporting and principal-model parity passed: ',name)
}
write.csv(do.call(rbind,rows),targets[1],row.names=FALSE)
write.csv(do.call(rbind,info),targets[2],row.names=FALSE)
saveRDS(models,targets[3])
writeLines(capture.output(sessionInfo()),targets[4])
