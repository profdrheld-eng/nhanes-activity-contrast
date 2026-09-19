#!/usr/bin/env Rscript
arg <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))])
HERE <- dirname(normalizePath(arg))
source(file.path(HERE,'model_common.R'))
targets <- file.path(OUT,c('mortality_results.csv','model_support.csv','joint_tests.csv','spline_curves.csv','spline_support.csv','age_results.csv','age_curves.csv','primary_fits.rds','session_info.txt'))
if(any(file.exists(targets))) stop('Refusing to overwrite completed or partial stage outputs')
for(name in names(groups)) {
  x<-groups[[name]]
  for(inter in c(FALSE,TRUE)) {
    type<-if(inter)'interaction' else 'common'
    id<-paste(name,type,'primary',sep='__');fits[[id]]<-fit_model(x,id,inter)
    for(lag in c(2,3,5)) fit_model(x,paste(name,type,paste0('exclude_deaths_',lag,'y'),sep='__'),inter,lag=lag)
    for(h in c(5,8)) fit_model(x,paste(name,type,paste0('censor_',h,'y'),sep='__'),inter,horizon=h)
    for(alt in c('equal','zlog','zraw')) fit_model(x,paste(name,type,alt,sep='__'),inter,delta=paste0('Delta_',alt),level=paste0('Level_',alt))
  }
  message('Core and sensitivity models finished: ',name)
  if(name %in% c('hip','wrist')) {
    for(lag in c(3,5)) fit_model(x,paste(name,'common',paste0('landmark_',lag,'y'),sep='__'),lag=lag,landmark=TRUE)
    fit_model(x,paste(name,'quintiles',sep='__'),delta='quintile')
    matched<-is.finite(x$BMI)&!is.na(x$disease)
    fit_model(x,paste(name,'matched_core',sep='__'),keep=matched)
    fit_model(x,paste(name,'matched_BMI',sep='__'),keep=matched,extra='+ BMI')
    fit_model(x,paste(name,'matched_BMI_disease',sep='__'),keep=matched,extra='+ BMI + disease')
    fit_model(x,paste(name,'exclude_CVD_cancer',sep='__'),keep=x$cvd_cancer=='no')
    # Fixed three-knot restricted cubic spline: linear column plus one nonlinear column.
    primary<-x[x$primary_complete_case,];kn<-wq(primary$Delta,primary$weight,c(.1,.5,.9))
    basis<-function(z) ((pmax(z-kn[1],0)^3)-pmax(z-kn[2],0)^3*(kn[3]-kn[1])/(kn[3]-kn[2])+pmax(z-kn[3],0)^3*(kn[2]-kn[1])/(kn[3]-kn[2]))/(kn[3]-kn[1])^2
    x$nonlinear<-basis(x$Delta)
    for(inter in c(FALSE,TRUE)){
      id<-paste(name,if(inter)'interaction' else 'common','spline',sep='__')
      out<-fit_model(x,id,interaction=inter,extra=if(inter)'+ nonlinear * sex' else '+ nonlinear')
      f<-out$fit;des<-out$domain
      terms<-names(coef(f))[grepl('Delta|nonlinear',names(coef(f)))]
      nt<-terms[grepl('nonlinear',terms)]
      tests[[paste(id,'global')]]<-joint(f,des,terms,id,'Global contrast')
      tests[[paste(id,'nonlinear')]]<-joint(f,des,nt,id,'Deviation from linearity')
      grid<-seq(wq(primary$Delta,primary$weight,.01),wq(primary$Delta,primary$weight,.99),length.out=151)
      for(sex in if(inter)c('Male','Female') else 'Common') {
        rows<-lapply(grid,function(g){
          L<-c(Delta=g,nonlinear=basis(g)-basis(0))
          if(sex=='Female'){
            i1<-names(coef(f))[grepl('Delta:sexfemale|sexfemale:Delta',names(coef(f)))]
            i2<-names(coef(f))[grepl('nonlinear:sexfemale|sexfemale:nonlinear',names(coef(f)))]
            L<-c(L,setNames(g,i1),setNames(basis(g)-basis(0),i2))
          }
          rr<-extract(f,des,L,id,sex);rr<-rr[rr$inference=='full_design_t',];rr$Delta<-g;rr
        });curves[[paste(id,sex)]]<-do.call(rbind,rows)
      }
    }
    spline_support[[name]]<-data.frame(era=name,p01=wq(primary$Delta,primary$weight,.01),p05=wq(primary$Delta,primary$weight,.05),
      k10=kn[1],k50=kn[2],k90=kn[3],p95=wq(primary$Delta,primary$weight,.95),p99=wq(primary$Delta,primary$weight,.99),
      low_tail_n=sum(primary$Delta<kn[1]),high_tail_n=sum(primary$Delta>kn[3]),
      low_tail_deaths=sum(primary$death[primary$Delta<kn[1]]),high_tail_deaths=sum(primary$death[primary$Delta>kn[3]]))
  }
  for(alt in c('primary','equal','zlog','zraw')) {
    delta<-if(alt=='primary')'Delta' else paste0('Delta_',alt)
    level<-if(alt=='primary')'Level' else paste0('Level_',alt)
    des<-subset(design(x),age_domain)
    ff<-as.formula(paste(delta,'~',level,'+ age10 + I(age10^2) + sex + race4 + education3 + PIR',if(name%in%c('hip','wrist'))'+ cycle' else ''))
    f<-checked_model(svyglm(ff,design=des));id<-paste(name,'age',alt,sep='__')
    age_results[[id]]<-extract(f,des,c(age10=3,'I(age10^2)'=9),id,'Age 80 versus 50',transform=FALSE)
    tests[[paste(id,'age')]]<-joint(f,des,c('age10','I(age10^2)'),id,'Global age')
    tests[[paste(id,'age2')]]<-joint(f,des,'I(age10^2)',id,'Quadratic age')
    support[[id]]<-data.frame(model=id,n=nrow(des$variables),deaths=sum(des$variables$death),
      strata=length(unique(des$variables$design_stratum)),PSUs=length(unique(des$variables$design_psu)),design_df=degf(des),
      parameters=length(coef(f)),residual_df=resdf(f),formula=paste(deparse(ff),collapse=' '))
    age_curves[[id]]<-do.call(rbind,lapply(seq(20,80,1),function(age){a<-(age-50)/10;rr<-extract(f,des,c(age10=a,'I(age10^2)'=a*a),id,'Age versus 50',FALSE);rr<-rr[rr$inference=='full_design_t',];rr$age<-age;rr}))
  }
}

for(h in c(Inf,5,8)) {
  message('Regime heterogeneity horizon: ',h)
  x<-d;x$time<-pmin(x$follow_up_years,h);x$death<-as.integer(x$MORTSTAT==1&x$follow_up_years<=h)
  for(inter in c(FALSE,TRUE)) {
    des<-subset(design(x),primary_complete_case)
    rhs<-if(inter)'Delta * sex + Delta:era + sex:era + Delta:sex:era' else 'Delta + sex + Delta:era'
    ff<-as.formula(paste('Surv(time,death)~',rhs,'+ Level + age + race4 + education3 + PIR + smoking3 + strata(cycle)'))
    f<-checked_model(svycoxph(ff,design=des,method='efron'));id<-paste('era_heterogeneity',if(inter)'interaction' else 'common',h,sep='__')
    terms<-names(coef(f))[grepl('Delta.*era|era.*Delta',names(coef(f)))]
    tests[[id]]<-joint(f,des,terms,id,'Regime heterogeneity')
    if(!inter)results[[id]]<-extract(f,des,setNames(1,terms),id,'Wrist/Hip slope ratio')
  }
}
for(name in c('hip','wrist')){
  message('Cycle heterogeneity: ',name)
  x<-groups[[name]];des<-subset(design(x),primary_complete_case)
  f<-checked_model(svycoxph(Surv(time,death)~Delta + sex + Delta:cycle + Level + age + race4 + education3 + PIR + smoking3 + strata(cycle),design=des,method='efron'))
  terms<-names(coef(f))[grepl('Delta.*cycle|cycle.*Delta',names(coef(f)))]
  tests[[paste(name,'cycle_heterogeneity')]]<-joint(f,des,terms,paste(name,'cycle_heterogeneity'),'Cycle heterogeneity')
}
write.csv(do.call(rbind,results),file.path(OUT,'mortality_results.csv'),row.names=FALSE)
write.csv(do.call(rbind,support),file.path(OUT,'model_support.csv'),row.names=FALSE)
write.csv(do.call(rbind,tests),file.path(OUT,'joint_tests.csv'),row.names=FALSE)
write.csv(do.call(rbind,curves),file.path(OUT,'spline_curves.csv'),row.names=FALSE)
write.csv(do.call(rbind,spline_support),file.path(OUT,'spline_support.csv'),row.names=FALSE)
write.csv(do.call(rbind,age_results),file.path(OUT,'age_results.csv'),row.names=FALSE)
write.csv(do.call(rbind,age_curves),file.path(OUT,'age_curves.csv'),row.names=FALSE)
saveRDS(fits,file.path(OUT,'primary_fits.rds'))
writeLines(capture.output(sessionInfo()),file.path(OUT,'session_info.txt'))
message('Revision models finished')
