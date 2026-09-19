#!/usr/bin/env Rscript
# Figure S1: mutually exclusive terminal dispositions, shared count-width scale.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1) stop('Supply exactly one local work directory')
WORK <- normalizePath(args[1],mustWork=TRUE)
OUT <- file.path(WORK,'displays');dir.create(OUT,showWarnings=FALSE)
target <- file.path(OUT,'participant_flow_sankey.png')
targets <- file.path(OUT,c('participant_flow_sankey.png','participant_flow_checks.csv'))
if(any(file.exists(targets))) stop('Refusing to overwrite completed or partial stage outputs')
rows <- read.csv(file.path(WORK,'models/participant_flow.csv'))
cycles <- c('2003-2004','2005-2006','2011-2012','2013-2014')
value <- function(cycle,stage,field='remaining') {
  x <- rows[rows$cycle==cycle & rows$stage==stage,field]
  stopifnot(length(x)==1,is.finite(x),x>=0,x==as.integer(x));x
}
checks <- list()
check <- function(label,ok) {
  if(!isTRUE(ok)) stop(label)
  checks[[length(checks)+1]] <<- data.frame(check=label,passed=TRUE)
}
png(target,width=3600,height=2980,res=300,type='cairo',bg='white')
par(mar=rep(0,4),xaxs='i',yaxs='i')
plot.new();plot.window(xlim=c(0,3600),ylim=c(2980,0))
label <- function(x,y,s,size=34,bold=FALSE) {
  cex <- size/50;font <- if(bold)2 else 1
  check(paste('Text bounds',s),x>=0 && y>=0 && x+max(strwidth(s,cex=cex,font=font))<3600 &&
        y+abs(strheight(s,cex=cex,font=font))<2980)
  text(x,y,s,adj=c(0,1),cex=cex,font=font)
}
ribbon <- function(x0,y0,x1,y1,height,color) {
  t <- seq(0,1,length.out=161)
  x <- x0+(x1-x0)*t;y <- y0+(y1-y0)*(3*t*t-2*t*t*t)
  polygon(c(x,rev(x)),c(y,rev(y+height)),col=color,border=NA)
}
label(90,35,'Participant selection across four NHANES cycles',64,TRUE)
label(90,118,'Examined adults aged 20+ years | Unweighted counts | Flow widths proportional to participant numbers',34)
scale <- .13
for(i in seq_along(cycles)) {
  cy <- cycles[i];ox <- 90+((i-1)%%2)*1780;oy <- 220+((i-1)%/%2)*1190
  colors <- if(i<=2)c('#719CE8','#A8C2F1','#E2EBFA') else c('#EE9487','#F4B8AD','#FBE4DE')
  n <- value(cy,0);final <- value(cy,9)
  losses <- c(value(cy,3,'excluded')+value(cy,5,'excluded'),value(cy,6,'excluded'),value(cy,7,'excluded'),value(cy,9,'excluded'))
  check(paste(cy,'partition'),final+sum(losses)==n)
  for(s in 1:9) check(paste(cy,'sequential stage',s),value(cy,s-1)-value(cy,s)==value(cy,s,'excluded'))
  check(paste(cy,'zero-loss stages'),all(vapply(c(1,2,4,8),function(s)value(cy,s,'excluded')==0,logical(1))))
  label(ox,oy,paste(LETTERS[i],cy),49,TRUE)
  label(ox,oy+63,if(i<=2)'Hip-worn monitor' else 'Wrist-worn monitor',33)
  counts <- c(final,losses)
  names <- c('Complete case','Linkage / follow-up','Incomplete self-report','<4 valid device days*','Missing covariates')
  dest <- c(oy+155,oy+155+cumsum(counts[1:4]*scale+85))
  x0 <- ox+290;x1 <- ox+950;start <- oy+235;cursor <- start
  for(j in 1:5) {
    height <- counts[j]*scale;dy <- dest[j]
    ribbon(x0+20,cursor,x1,dy,height,if(j==1)colors[2] else colors[3])
    rect(x1,dy,x1+20,dy+height,col=if(j==1)colors[1] else colors[2],border=NA)
    label(x1+44,dy+max(0,(height-78)/2),paste0(names[j],'\nn = ',format(counts[j],big.mark=',')),35,j==1)
    cursor <- cursor+height
  }
  rect(x0,start,x0+20,start+n*scale,col='#25282C',border=NA)
  label(ox+5,start+n*scale/2-55,paste0('Examined\nn = ',format(n,big.mark=',')),38,TRUE)
  label(ox+335,oy+1060,paste0('Deaths in complete cases: ',format(value(cy,10),big.mark=',')),31)
}
label(90,2800,'* Includes participants without usable monitor recordings. Exclusions are sequential and mutually exclusive.',30)
details <- vapply(cycles,function(cy)paste(value(cy,3,'excluded'),value(cy,5,'excluded'),sep=' / '),'')
label(90,2850,paste0('Linkage / follow-up exclusions by cycle: ',paste(details,collapse=', '),'. Other eligibility checks excluded no additional participants.'),28)
label(90,2900,'NHANES: National Health and Nutrition Examination Survey. Equal flow-width scale across panels; terminal outcomes, not follow-up trajectories.',28)
dev.off()
write.csv(do.call(rbind,checks),file.path(OUT,'participant_flow_checks.csv'),row.names=FALSE)
