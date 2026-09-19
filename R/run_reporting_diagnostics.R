#!/usr/bin/env Rscript
# Post-hoc diagnostics for all 18 same-case S36 models. No model selection.
arg <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))])
HERE <- dirname(normalizePath(arg))
source(file.path(HERE,'model_common.R'))
source(file.path(HERE,'reporting_time_split.R'))
destination <- file.path(WORK,'reporting_diagnostics')
if(dir.exists(destination))stop('Refusing to overwrite diagnostic outputs')
dir.create(destination)
saved <- readRDS(file.path(WORK,'reporting','adjustment_fits.rds'))
stopifnot(length(saved)==18)
ph <- list();steps <- list();checks <- list();plots <- list()
for(id in names(saved)) {
  obj <- saved[[id]];parts <- strsplit(id,'__',fixed=TRUE)[[1]]
  name <- parts[1];adjustment <- parts[2];x <- groups[[name]]
  des <- subset(design(x),primary_complete_case)
  d <- des$variables;ff <- formula(obj$fit)
  stopifnot(identical(as.numeric(d$SEQN),as.numeric(obj$domain$variables$SEQN)))
  # Same normalized weights as svycoxph; zph remains a conventional diagnostic,
  # not a design-based test, even though the fitted model has robust covariance.
  d$analysis_weight <- d$weight/mean(d$weight)
  cf <- checked_model(coxph(ff,data=d,weights=analysis_weight,
    ties='efron',robust=TRUE,cluster=design_psu,x=TRUE,y=TRUE,model=TRUE))
  difference <- max(abs(coef(cf)-coef(obj$fit)))
  stopifnot(difference<1e-8)
  z <- cox.zph(cf,transform='km')
  ph[[id]] <- data.frame(model=id,group=name,adjustment=adjustment,
    term=rownames(z$table),z$table,row.names=NULL)
  png(file.path(destination,paste0(id,'_schoenfeld.png')),width=1800,height=1200,res=150)
  plot(z,var='Delta',main=paste(name,adjustment,'contrast diagnostic'))
  abline(h=coef(cf)['Delta'],lty=2,col='grey40')
  dev.off()
  split <- reporting_time_split(x)
  sdes <- subset(design(split),primary_complete_case)
  ns <- nrow(sdes$variables)
  # Splitting without time interactions must reproduce coefficients AND survey
  # covariance, proving that weights, PSUs, strata and domains were preserved.
  sf <- update(ff,Surv(tstart,time,death)~.)
  base <- checked_model(svycoxph(sf,design=sdes,method='efron'))
  db <- max(abs(coef(base)-coef(obj$fit)))
  dv <- max(abs(vcov(base)-vcov(obj$fit)))
  stopifnot(db<1e-8,dv<1e-8,degf(sdes)==degf(des),
    sum(sdes$variables$death)==sum(d$death),
    abs(sum(sdes$variables$time-sdes$variables$tstart)-sum(d$time))<1e-8)
  # Focused check: only Delta gets a step interaction, preserving each S36
  # adjustment set. Covariate time variation is screened separately above.
  tf <- update(sf,.~.+late_delta)
  fit <- checked_model(svycoxph(tf,design=sdes,method='efron'))
  contrasts <- list(early=c(Delta=1),late=c(Delta=1,late_delta=1),
    late_to_early=c(late_delta=1))
  for(label in names(contrasts)) {
    row <- extract(fit,sdes,contrasts[[label]],id,label)
    row <- row[row$inference=='full_design_t',]
    row$n <- nrow(d);row$group <- name;row$adjustment <- adjustment
    row$interval_rows <- ns;row$cut_years <- 5
    steps[[paste(id,label)]] <- row
  }
  checks[[id]] <- data.frame(model=id,n=nrow(d),deaths=sum(d$death),
    n_after5=sum(d$time>5),deaths_early=sum(d$death[d$time<=5]),
    deaths_after5=sum(d$death[d$time>5]),design_df=degf(des),
    interval_rows=ns,person_years=sum(d$time),
    coxph_coefficient_difference=difference,
    split_coefficient_difference=db,split_covariance_difference=dv)
  message('Reporting diagnostics passed: ',id)
}
write.csv(do.call(rbind,ph),file.path(destination,'schoenfeld.csv'),row.names=FALSE)
write.csv(do.call(rbind,steps),file.path(destination,'step_time.csv'),row.names=FALSE)
write.csv(do.call(rbind,checks),file.path(destination,'support.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(destination,'session_info.txt'))
writeLines('18 model diagnostics complete',file.path(destination,'COMPLETE.txt'))
