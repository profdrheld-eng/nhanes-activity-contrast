# Weighted tied midranks and population standardization for processing variants.
midrank<-function(x,w){
  o<-order(x,method='radix');xs<-x[o];ws<-w[o]
  starts<-c(1,which(diff(xs)!=0)+1);ends<-c(starts[-1]-1,length(x))
  p<-numeric(length(x));total<-sum(w);cum<-0
  for(i in seq_along(starts)){ii<-starts[i]:ends[i];tw<-sum(ws[ii]);p[ii]<-(cum+tw/2)/total;cum<-cum+tw}
  z<-numeric(length(x));z[o]<-qnorm(pmin(.9999,pmax(.0001,p)));z
}
stand<-function(x,w){m<-wm(x,w);(x-m)/sqrt(wm((x-m)^2,w))}
