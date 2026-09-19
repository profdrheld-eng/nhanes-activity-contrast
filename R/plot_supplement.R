#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1) stop('Supply exactly one local work directory')
WORK <- normalizePath(args[1],mustWork=TRUE)
OUT <- file.path(WORK,'models');FIG <- file.path(WORK,'displays')
dir.create(FIG,showWarnings=FALSE)
targets <- c('mortality_spline_common.png','mortality_spline_interaction.png','age_curves_alternative_definitions.png')
if(any(file.exists(file.path(FIG,targets)))) stop('Refusing to overwrite supplementary figures')
cur<-read.csv(file.path(OUT,'spline_curves.csv'));sp<-read.csv(file.path(OUT,'spline_support.csv'));tests<-read.csv(file.path(OUT,'joint_tests.csv'))
for(type in c('common','interaction')){
  png(file.path(FIG,paste0('mortality_spline_',type,'.png')),width=3600,height=1550,res=300,type='cairo',bg='white')
  par(mfrow=c(1,2),mar=c(5.5,5.5,4.5,1.5),font=2,font.axis=2,font.lab=2,font.main=2,las=1)
  for(era in c('hip','wrist')){
    id<-paste(era,type,'spline',sep='__');x<-cur[cur$model==id,];s<-sp[sp$era==era,]
    ylim<-range(c(x$ci_low,x$ci_high));p<-tests$p_full[tests$model==id&tests$test=='Deviation from linearity']
    plot(NA,xlim=c(s$p01,s$p99),ylim=ylim,log='y',xaxs='i',xlab='',ylab='',bty='l',
         main=paste0(if(era=='hip')'A  2003–2006 (hip)' else 'B  2011–2014 (wrist)','\nNonlinearity p = ',formatC(p,digits=3,format='f')))
    mtext('Relative self-report–device contrast (SD)',side=1,line=3.4,font=2,cex=.92)
    mtext('Mortality hazard ratio (reference = 0)',side=2,line=4,font=2,las=0,cex=.92)
    rect(s$p01,ylim[1],s$p05,ylim[2],col=adjustcolor('grey60',.15),border=NA)
    rect(s$p95,ylim[1],s$p99,ylim[2],col=adjustcolor('grey60',.15),border=NA)
    abline(h=1,lty=2,col='grey35');abline(v=0,lty=3,col='grey65')
    groups<-if(type=='common')'Common' else c('Male','Female');cols<-if(type=='common')'#427A9B' else c('#86B6DE','#E899AF')
    for(i in seq_along(groups)){
      y<-x[x$contrast==groups[i],]
      polygon(c(y$Delta,rev(y$Delta)),c(y$ci_low,rev(y$ci_high)),col=adjustcolor(cols[i],.23),border=NA)
      lines(y$Delta,y$estimate,lwd=3,col=cols[i])
    }
    legend('topleft',legend=groups,col=cols,lwd=3,bty='n',cex=.9)
  }
  dev.off()
}
a<-read.csv(file.path(OUT,'age_curves.csv'))
png(file.path(FIG,'age_curves_alternative_definitions.png'),width=3400,height=2600,res=300,type='cairo',bg='white')
par(mfrow=c(2,2),mar=c(5,5.5,2.5,1),font=2,font.axis=2,font.lab=2,las=1)
for(cy in c('2003-2004','2005-2006','2011-2012','2013-2014')){
  x<-a[startsWith(a$model,paste0(cy,'__age__')),];yp<-x[x$model==paste0(cy,'__age__primary'),]
  plot(NA,xlim=c(20,80),ylim=range(c(x$estimate,yp$ci_low,yp$ci_high)),xaxs='i',xlab='Recorded age (years)',ylab='',main=cy,bty='l')
  mtext('Contrast difference vs age 50 (SD)',side=2,line=4,font=2,las=0,cex=.9)
  abline(h=0,lty=2,col='grey60')
  cols<-c('#333333','#548FB2','#C27B91','#67856C');alts<-c('primary','equal','zlog','zraw')
  for(i in seq_along(alts)){
    y<-x[x$model==paste0(cy,'__age__',alts[i]),]
    if(i==1)polygon(c(y$age,rev(y$age)),c(y$ci_low,rev(y$ci_high)),col=adjustcolor(cols[i],.12),border=NA)
    lines(y$age,y$estimate,col=cols[i],lty=i,lwd=2.5)
  }
  if(cy=='2003-2004')legend('topleft',c('Rank, vigorous ×2','Rank, vigorous ×1','Log-score z contrast','Raw-score z contrast'),col=cols,lty=1:4,lwd=2.5,bty='n',cex=.72)
}
dev.off()
