#!/usr/bin/env Rscript
# Additional post-hoc comparison of age relationships, not an isolated device effect.
# Original cycle-standardized scores are retained in every restricted domain.
args <- commandArgs(TRUE)
if(length(args)!=1)stop('Usage: Rscript R/run_age_comparison.R /path/to/work')
WORK <- normalizePath(args[1],mustWork=TRUE)
script <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))])
ROOT <- dirname(dirname(normalizePath(script)))
if(WORK==ROOT || startsWith(WORK,paste0(ROOT,.Platform$file.sep)))
  stop('Work directory must be outside the code repository')
OUT <- file.path(WORK,'age_comparison')
if(file.exists(OUT))stop('Refusing to overwrite age-comparison outputs')
stopifnot(file.exists(file.path(WORK,'models/revision_cohort.csv.gz')),
          file.exists(file.path(WORK,'models/age_results.csv')))
if(!dir.create(OUT))stop('Cannot create age-comparison output directory')
options(stringsAsFactors=FALSE,survey.lonely.psu='adjust')
suppressPackageStartupMessages(library(survey))
d<-read.csv(gzfile(file.path(WORK,'models/revision_cohort.csv.gz')))
flag<-function(x)!is.na(x)&tolower(as.character(x))%in%c('true','1')
d$primary_complete_case<-flag(d$primary_complete_case)
d$cycle<-factor(d$cycle,levels=c('2003-2004','2005-2006','2011-2012','2013-2014'))
d$era<-factor(d$era,levels=c('hip','wrist'));d$sex<-factor(d$sex,levels=c('male','female'))
d$race4<-factor(d$race_ethnicity4,levels=c('Non-Hispanic White','Hispanic','Non-Hispanic Black','Other'))
d$education3<-factor(d$education3,levels=c('more than high school','high school/GED','under high school'))
d$PIR<-as.numeric(d$PIR)
for(v in c('diabetes_history','CVD_history','cancer_history'))d[[v]]<-factor(d[[v]],levels=c('no','yes'))
d$disease<-with(d,ifelse(diabetes_history=='yes'|CVD_history=='yes'|cancer_history=='yes','yes',ifelse(diabetes_history=='no'&CVD_history=='no'&cancer_history=='no','no',NA)))
d$cvd_cancer<-with(d,ifelse(CVD_history=='yes'|cancer_history=='yes','yes',ifelse(CVD_history=='no'&cancer_history=='no','no',NA)))
d$historical<-with(d,primary_complete_case & is.finite(BMI)&ifelse(era=='hip',!is.na(disease),!is.na(cvd_cancer)))
d$harmonized<-with(d,primary_complete_case & is.finite(BMI)&!is.na(disease))
# Construct the survey design before subsetting; preserve all PSUs/strata.
d$stratum<-interaction(d$cycle,d$SDMVSTRA,drop=TRUE);d$psu<-interaction(d$cycle,d$SDMVSTRA,d$SDMVPSU,drop=TRUE);d$weight<-d$WTMEC2YR/2
stopifnot(!anyNA(d$weight),all(d$weight>0),!anyDuplicated(paste(d$cycle,d$SEQN)))
make_design<-function(x)svydesign(ids=~psu,strata=~stratum,weights=~weight,nest=TRUE,data=x)
checkfit<-function(f){stopifnot(isTRUE(f$converged),all(is.finite(coef(f))),all(is.finite(vcov(f))),f$rank==length(coef(f)));f}
# Reproduce the previously reported age analyses before calculating new tests.
baseline<-read.csv(file.path(WORK,'models/age_results.csv'));reproduction<-list()
d$age10<-(d$age-50)/10;d$age2<-d$age10^2
for(era in levels(d$era)){
 x<-droplevels(d[d$era==era,]);des<-subset(make_design(x),historical)
 f<-checkfit(svyglm(Delta~age10+age2+Level+sex+race4+education3+PIR+cycle,design=des))
 estimate<-sum(coef(f)[c('age10','age2')]*c(3,9));old<-baseline[baseline$model==paste0(era,'__age__primary')&baseline$inference=='full_design_t',]
 stopifnot(nrow(old)==1,abs(estimate-old$estimate)<1e-10,nrow(des$variables)==old$n)
 reproduction[[era]]<-data.frame(era=era,n=nrow(des$variables),estimate=estimate,previous=old$estimate)
}
write.csv(do.call(rbind,reproduction),file.path(OUT,'baseline_reproduction.csv'),row.names=FALSE)
joints<-list();support<-list();contrasts<-list();curves<-list();coefficients<-list();fits<-list();verification<-list()
basef<-~age10+age2+Level+sex+race4+education3+PIR
# Fully interact nuisance terms, with cycle intercepts absorbing the era intercept.
for(id in c('historical_under80','harmonized_under80','harmonized_cap80')){
 x<-d;if(id=='harmonized_cap80')x$age10<-(pmin(x$age,80)-50)/10
 x$age2<-x$age10^2;x$keep<-if(id=='historical_under80')x$historical else x$harmonized
 x$keep<-x$keep & x$age>=20 & if(id=='harmonized_cap80')TRUE else x$age<80
 x$keep[is.na(x$keep)]<-FALSE
 mm<-model.matrix(basef,model.frame(basef,x,na.action=na.pass))[, -1,drop=FALSE];stopifnot(nrow(mm)==nrow(x))
 colnames(mm)<-make.names(colnames(mm));gn<-paste0('g_',colnames(mm));wn<-paste0('w_',colnames(mm))
 for(j in seq_len(ncol(mm))){x[[gn[j]]]<-mm[,j];x[[wn[j]]]<-mm[,j]*as.integer(x$era=='wrist')}
 full<-make_design(x);des<-subset(full,keep);ff<-reformulate(c('cycle',gn,wn),response='Delta')
 fit<-checkfit(svyglm(ff,design=des));b<-coef(fit);v<-vcov(fit);terms<-c('w_age10','w_age2');q<-2
 beta<-b[terms];V<-v[terms,terms];stopifnot(qr(V)$rank==2)
 F<-drop(t(beta)%*%solve(V,beta))/q;ddf<-degf(des);rdf<-fit$df.residual
 for(conv in c('full_design','residual','asymptotic')){
  df<-switch(conv,full_design=ddf,residual=rdf,asymptotic=Inf)
  test<-regTermTest(fit,~w_age10+w_age2,df=df)
  pv<-if(is.infinite(df))pchisq(F*q,q,lower.tail=FALSE) else pf(F,q,df,lower.tail=FALSE)
  stopifnot(abs(pv-as.numeric(test$p))<1e-12)
  joints[[length(joints)+1]]<-data.frame(model=id,inference=conv,F=F,numerator_df=q,denominator_df=df,p=pv)
 }
 # Independent period fits must agree with the full-interaction parameterization.
 separate<-list()
 for(era in levels(x$era)){
  xd<-droplevels(x[x$era==era,]);sd<-subset(make_design(xd),keep)
  separate[[era]]<-checkfit(svyglm(Delta~age10+age2+Level+sex+race4+education3+PIR+cycle,design=sd))
  support[[length(support)+1]]<-data.frame(model=id,era=era,n=nrow(sd$variables),min_age=min(sd$variables$age),max_recorded_age=max(sd$variables$age),strata=length(unique(sd$variables$stratum)),psu=length(unique(sd$variables$psu)),design_df=degf(sd))
 }
 diffb<-coef(separate$wrist)[c('age10','age2')]-coef(separate$hip)[c('age10','age2')]
 diffv<-vcov(separate$wrist)[c('age10','age2'),c('age10','age2')]+vcov(separate$hip)[c('age10','age2'),c('age10','age2')]
 stopifnot(max(abs(diffb-beta))<1e-10,max(abs(diffv-V))<1e-10)
 verification[[id]]<-data.frame(model=id,coefficient_error=max(abs(diffb-beta)),covariance_error=max(abs(diffv-V)),design_df=ddf,residual_df=rdf,parameters=length(b),weighted_design_condition=kappa(model.matrix(fit)*sqrt(weights(des))))
 scalar<-function(L,df){est<-sum(L*b);se<-sqrt(max(0,drop(t(L)%*%v%*%L)));crit<-if(is.infinite(df))qnorm(.975) else qt(.975,df);data.frame(estimate=est,se=se,ci_low=est-crit*se,ci_high=est+crit*se,p=if(se==0)NA else if(is.infinite(df))2*pnorm(-abs(est/se)) else 2*pt(-abs(est/se),df))}
 for(age in 20:79){
  a<-(age-50)/10
  for(type in c('hip','wrist','wrist_minus_hip')){
   L<-setNames(rep(0,length(b)),names(b));if(type!='wrist_minus_hip')L[c('g_age10','g_age2')]<-c(a,a*a)
   if(type!='hip')L[c('w_age10','w_age2')]<-c(a,a*a)
   curves[[length(curves)+1]]<-cbind(data.frame(model=id,age=age,type=type),scalar(L,ddf))
   if(age==75)for(conv in c('full_design','residual','asymptotic')){
    df<-switch(conv,full_design=ddf,residual=rdf,asymptotic=Inf)
    contrasts[[length(contrasts)+1]]<-cbind(data.frame(model=id,age=age,reference=50,type=type,inference=conv,denominator_df=df),scalar(L,df))
   }
  }
 }
 coefficients[[id]]<-data.frame(model=id,term=names(b),coefficient=unname(b),se=sqrt(diag(v)))
 fits[[id]]<-fit
}
for(name in c('support','contrasts','curves','coefficients','verification'))write.csv(do.call(rbind,get(name)),file.path(OUT,paste0(name,'.csv')),row.names=FALSE)
write.csv(do.call(rbind,joints),file.path(OUT,'joint_tests.csv'),row.names=FALSE)
saveRDS(fits,file.path(OUT,'models.rds'));writeLines(capture.output(sessionInfo()),file.path(OUT,'session_info.txt'))
print(do.call(rbind,joints));print(do.call(rbind,support));print(subset(do.call(rbind,contrasts),inference=='full_design'))

# A partial or failed fit must never be accepted by the table exporter.
writeLines('3 age-comparison models complete',file.path(OUT,'COMPLETE.txt'))
