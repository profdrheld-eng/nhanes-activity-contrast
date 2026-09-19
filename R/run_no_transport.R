#!/usr/bin/env Rscript
# Post-hoc hip-score sensitivity, preserving original reference and outcome samples.
arg <- sub('^--file=','',commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))])
HERE <- dirname(normalizePath(arg))
source(file.path(HERE,'model_common.R'))
names_out <- c('no_transport_mortality.csv','no_transport_age.csv','no_transport_support.csv',
 'no_transport_joint_tests.csv','no_transport_score_summary.csv','no_transport_scores.csv.gz')
if(any(file.exists(file.path(OUT,names_out))))stop('Refusing to overwrite completed or partial stage outputs')
source(file.path(HERE,'score_helpers.R'))
source(file.path(HERE,'no_transport_scores.R'))
frames <- list(); summaries <- list(); age_support <- list()
wcorr <- function(x,y,w) {
  x<-x-wm(x,w);y<-y-wm(y,w)
  wm(x*y,w)/sqrt(wm(x*x,w)*wm(y*y,w))
}
for(cy in c('2003-2004','2005-2006')) {
  original <- groups[[cy]];x <- no_transport_scores(original);ref <- x$reference_domain
  stopifnot(identical(x[names(original)],original),
    all(is.finite(x$Delta_no_transport[x$primary_complete_case])),
    all(is.finite(x$Level_no_transport[x$age_domain])))
  w<-x$WTMEC2YR[ref]
  # Reconstruct the original device ranks to ensure the retained D is appropriate.
  stopifnot(max(abs(midrank(x$device_log[ref],w)-x$D[ref]))<1e-9)
  summaries[[cy]]<-data.frame(cycle=cy,reference_n=sum(ref),mortality_n=sum(x$primary_complete_case),
    deaths=sum(x$death[x$primary_complete_case]),age_n=sum(x$age_domain),
    original_index_mean=wm(x$SR_day[ref],w),no_transport_index_mean=wm(x$SR_no_transport[ref],w),
    self_report_rank_correlation=wcorr(x$S[ref],x$S_no_transport[ref],w),
    contrast_correlation=wcorr(x$Delta[ref],x$Delta_no_transport[ref],w),
    common_component_correlation=wcorr(x$Level[ref],x$Level_no_transport[ref],w))
  frames[[cy]]<-x
}
frames[['hip']]<-rbind(frames[['2003-2004']],frames[['2005-2006']])
for(name in names(frames)) {
  x<-frames[[name]]
  for(variant in c('original','no_transport')) {
    delta<-if(variant=='original')'Delta' else 'Delta_no_transport'
    level<-if(variant=='original')'Level' else 'Level_no_transport'
    for(inter in c(FALSE,TRUE)) {
      id<-paste(name,if(inter)'interaction' else 'common',variant,sep='__')
      fitted<-fit_model(x,id,inter,delta=delta,level=level)
      stopifnot(identical(as.numeric(fitted$domain$variables$SEQN),as.numeric(x$SEQN[x$primary_complete_case])))
    }
    des<-subset(design(x),age_domain)
    ff<-as.formula(paste(delta,'~',level,'+ age10 + I(age10^2) + sex + race4 + education3 + PIR',if(name=='hip')'+ cycle' else ''))
    f<-checked_model(svyglm(ff,design=des));id<-paste(name,'age',variant,sep='__')
    stopifnot(identical(as.numeric(des$variables$SEQN),as.numeric(x$SEQN[x$age_domain])))
    age_results[[id]]<-extract(f,des,c(age10=3,'I(age10^2)'=9),id,'Age 80 versus 50',FALSE)
    tests[[paste(id,'global')]]<-joint(f,des,c('age10','I(age10^2)'),id,'Global age')
    tests[[paste(id,'quadratic')]]<-joint(f,des,'I(age10^2)',id,'Quadratic age')
    age_support[[id]]<-data.frame(model=id,n=nrow(des$variables),deaths=sum(des$variables$death),
      strata=length(unique(des$variables$design_stratum)),PSUs=length(unique(des$variables$design_psu)),
      design_df=degf(des),parameters=length(coef(f)),residual_df=resdf(f),formula=paste(deparse(ff),collapse=' '))
  }
}
write.csv(do.call(rbind,results),file.path(OUT,names_out[1]),row.names=FALSE)
write.csv(do.call(rbind,age_results),file.path(OUT,names_out[2]),row.names=FALSE)
write.csv(do.call(rbind,c(support,age_support)),file.path(OUT,names_out[3]),row.names=FALSE)
write.csv(do.call(rbind,tests),file.path(OUT,names_out[4]),row.names=FALSE)
write.csv(do.call(rbind,summaries),file.path(OUT,names_out[5]),row.names=FALSE)
# This local participant-level audit export must never enter the public archive.
cols<-c('cycle','SEQN','reference_domain','primary_complete_case','age_domain','WTMEC2YR',
 'SR_day','SR_no_transport','S','D','Delta','Level','S_no_transport','Delta_no_transport','Level_no_transport')
con<-gzfile(file.path(OUT,names_out[6]),'wt');write.csv(frames[['hip']][cols],con,row.names=FALSE);close(con)
message('PASS: no-transport sensitivity completed on unchanged original reference and outcome domains')
