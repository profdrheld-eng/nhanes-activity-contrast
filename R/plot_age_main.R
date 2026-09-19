#!/usr/bin/env Rscript
# Main Figure 2 from freshly fitted common and sex-specific age curves.
# Base-R fonts avoid dependence on a proprietary Word font installation.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1) stop('Supply exactly one local work directory')
WORK <- normalizePath(args[1],mustWork=TRUE)
OUT <- file.path(WORK,'displays');dir.create(OUT,showWarnings=FALSE)
target <- file.path(OUT,'Figure_2_age_overlay.png')
if(file.exists(target)) stop('Refusing to overwrite Figure 2')
common <- read.csv(file.path(WORK,'models/age_curves.csv'))
overlay <- read.csv(file.path(WORK,'models/age_sex_overlay.csv'))
png(target,width=3400,height=1700,res=300,type='cairo',bg='white')
par(mfrow=c(1,2),mar=c(5.4,5.8,4.5,1.5),font=2,font.axis=2,font.lab=2,font.main=2,las=1)
for(era in c('hip','wrist')) {
  x <- common[common$model==paste0(era,'__age__primary'),]
  x <- x[order(x$age),]
  stopifnot(nrow(x)==61,all(x$inference=='full_design_t'),all(x$age==20:80))
  stopifnot(abs(x$estimate[x$age==50])<1e-10)
  title <- if(era=='hip') 'A  2003-2006 (hip)' else 'B  2011-2014 (wrist)'
  plot(NA,xlim=c(20,80),ylim=c(-.5,2),xaxs='i',yaxs='i',bty='l',
       xlab='Recorded age (years)',ylab='',main=paste0(title,'\nAdjusted age relationship | n = ',format(x$n[1],big.mark=',')))
  mtext('Contrast difference vs age 50 (SD)',side=2,line=4.2,las=0,font=2,cex=.95)
  abline(h=0,lty=2,col='grey60');abline(v=50,lty=3,col='grey75')
  colors <- c(Common='#333333',Male='#5AA9E6',Female='#E888A5')
  for(group in names(colors)) {
    y <- if(group=='Common') x else overlay[overlay$era==era & overlay$sex==group,]
    y <- y[order(y$age),]
    stopifnot(nrow(y)==61,all(y$age==20:80),all(is.finite(as.matrix(y[c('estimate','ci_low','ci_high')]))),
              min(y$ci_low)>=-.5,max(y$ci_high)<=2,abs(y$estimate[y$age==50])<1e-10)
    polygon(c(y$age,rev(y$age)),c(y$ci_low,rev(y$ci_high)),col=adjustcolor(colors[group],.14),border=NA)
    lines(y$age,y$estimate,col=colors[group],lwd=3)
  }
  legend('topleft',c('Common, sex-adjusted','Male','Female'),col=colors,lwd=3,bty='n',cex=.8)
}
dev.off()
