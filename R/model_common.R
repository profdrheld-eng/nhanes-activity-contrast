source(file.path(HERE,'fit_checks.R'))
# Shared current model specification, data domains, and contrast helpers.
# Sourced explicitly by each analysis stage; never parse/evaluate a script prefix.
options(stringsAsFactors=FALSE, survey.lonely.psu='adjust')
suppressPackageStartupMessages({library(survey);library(survival)})
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1) stop('Supply exactly one local work directory')
WORK <- normalizePath(args[1],mustWork=TRUE)
OUT <- file.path(WORK,'models')
d <- read.csv(gzfile(file.path(OUT,'revision_cohort.csv.gz')))
flag <- function(x) !is.na(x)&tolower(as.character(x)) %in% c('true','1')
for (v in c('adult_domain','reference_domain','primary_complete_case')) d[[v]] <- flag(d[[v]])
d$cycle <- factor(d$cycle,levels=c('2003-2004','2005-2006','2011-2012','2013-2014'))
d$era <- factor(d$era,levels=c('hip','wrist'))
d$sex <- factor(d$sex,levels=c('male','female'))
d$race4 <- factor(d$race_ethnicity4,levels=c('Non-Hispanic White','Hispanic','Non-Hispanic Black','Other'))
d$education3 <- factor(d$education3,levels=c('more than high school','high school/GED','under high school'))
d$smoking3 <- factor(d$smoking3,levels=c('never','former','current'))
for (v in c('CVD_history','cancer_history','diabetes_history')) d[[v]] <- factor(d[[v]],levels=c('no','yes'))
d$quintile <- relevel(factor(d$quintile,levels=1:5),ref='3')
d$death <- d$MORTSTAT;d$time <- d$follow_up_years
d$design_stratum <- interaction(d$cycle,d$SDMVSTRA,drop=TRUE)
d$design_psu <- interaction(d$cycle,d$SDMVSTRA,d$SDMVPSU,drop=TRUE)
d$weight <- d$WTMEC2YR/2
d$age10 <- (d$age-50)/10
d$disease <- ifelse(d$CVD_history=='yes'|d$cancer_history=='yes'|d$diabetes_history=='yes','yes',
                    ifelse(d$CVD_history=='no'&d$cancer_history=='no'&d$diabetes_history=='no','no',NA))
d$disease <- factor(d$disease,levels=c('no','yes'))
d$cvd_cancer <- ifelse(d$CVD_history=='yes'|d$cancer_history=='yes','yes',
                       ifelse(d$CVD_history=='no'&d$cancer_history=='no','no',NA))
# Reproduce historical age-domain differences explicitly; also test a harmonized domain.
d$age_domain <- d$primary_complete_case & is.finite(d$BMI) &
  ifelse(d$era=='hip',!is.na(d$disease),!is.na(d$cvd_cancer))
design <- function(x) svydesign(ids=~design_psu,strata=~design_stratum,weights=~weight,nest=TRUE,data=droplevels(x))
wm <- function(x,w) sum(x*w)/sum(w)
wq <- function(x,w,p) {o<-order(x);x<-x[o];w<-w[o];sapply(p,function(pp)x[which(cumsum(w)/sum(w)>=pp)[1]])}
results <- list(); support <- list(); tests <- list(); curves <- list(); age_results <- list(); age_curves <- list(); spline_support <- list(); fits <- list()
source(file.path(HERE,'inference.R'))
fit_model <- function(x,id,interaction=FALSE,delta='Delta',level='Level',extra='',keep=NULL,horizon=Inf,lag=0,landmark=FALSE){
  x$time<-pmin(x$follow_up_years,horizon);x$death<-as.integer(x$MORTSTAT==1 & x$follow_up_years<=horizon)
  x$keep<-x$primary_complete_case
  if(lag>0) {
    if(landmark) {x$keep<-x$keep & x$follow_up_years>lag;x$time<-x$time-lag}
    else x$keep<-x$keep & !(x$MORTSTAT==1 & x$follow_up_years<=lag)
  }
  if(!is.null(keep)) x$keep<-x$keep & keep
  x$keep[is.na(x$keep)]<-FALSE
  des<-subset(design(x),keep)
  rhs<-paste(delta,if(interaction)'* sex' else '+ sex','+',level,'+ age + race4 + education3 + PIR + smoking3',extra)
  if(length(unique(x$cycle))>1) rhs<-paste(rhs,'+ strata(cycle)')
  f<-checked_model(svycoxph(as.formula(paste('Surv(time,death)~',rhs)),design=des,method='efron',x=TRUE))
  stopifnot(all(is.finite(coef(f))),all(is.finite(vcov(f))),is.null(f$fail))
  support[[id]] <<- data.frame(model=id,n=nrow(des$variables),deaths=sum(des$variables$death),
    strata=length(unique(des$variables$design_stratum)),PSUs=length(unique(des$variables$design_psu)),
    design_df=degf(des),parameters=length(coef(f)),residual_df=resdf(f),formula=paste(deparse(formula(f)),collapse=' '))
  if(delta=='quintile') {
    for(term in names(coef(f))[grepl('^quintile',names(coef(f)))]) results[[paste(id,term)]]<<-extract(f,des,setNames(1,term),id,term)
    tests[[paste(id,'quintiles')]]<<-joint(f,des,names(coef(f))[grepl('^quintile',names(coef(f)))],id,'quintile joint')
  } else if(interaction){
    int<-names(coef(f))[grepl(paste0('^',delta,':sexfemale$|^sexfemale:',delta,'$'),names(coef(f)))]
    stopifnot(length(int)==1)
    results[[paste(id,'male')]]<<-extract(f,des,setNames(1,delta),id,'Male')
    results[[paste(id,'female')]]<<-extract(f,des,setNames(c(1,1),c(delta,int)),id,'Female')
    results[[paste(id,'interaction')]]<<-extract(f,des,setNames(1,int),id,'Female/Male slope ratio')
  } else results[[id]]<<-extract(f,des,setNames(1,delta),id,'Common contrast')
  list(fit=f,domain=des)
}

groups<-c(setNames(lapply(levels(d$cycle),function(cy)d[d$cycle==cy,]),levels(d$cycle)),
          list(hip=d[d$era=='hip',],wrist=d[d$era=='wrist',]))
