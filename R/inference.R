# Named linear contrasts and survey design-based reference distributions.
# Inputs: fitted coefficients/covariance and the fitted survey domain.
# Returns all three documented reference conventions without selection.
resdf <- function(f) if(!is.null(f$degf.resid)) f$degf.resid else f$df.residual
extract <- function(f,des,L,id,label,transform=TRUE) {
  b<-coef(f);v<-vcov(f)
  if(!length(L) || is.null(names(L)) || anyNA(names(L)) ||
     any(!nzchar(names(L))) || anyDuplicated(names(L)) ||
     any(!names(L) %in% names(b)) || any(!is.finite(L)))
    stop('Contrast requires unique known coefficient names and finite weights')
  v<-v[names(b),names(b),drop=FALSE]
  a<-setNames(rep(0,length(b)),names(b));a[names(L)]<-L
  est<-sum(a*b);se<-sqrt(drop(t(a)%*%v%*%a));df<-degf(des)
  stopifnot(is.finite(est),is.finite(se),se>=0)
  do.call(rbind,lapply(c('full_design_t','residual_t','normal'),function(conv){
    dd<-switch(conv,full_design_t=df,residual_t=resdf(f),normal=Inf)
    crit<-if(is.infinite(dd))qnorm(.975) else if(dd>0)qt(.975,dd) else NA_real_
    p<-if(se==0)NA_real_ else if(is.infinite(dd))2*pnorm(-abs(est/se)) else if(dd>0)2*pt(-abs(est/se),dd) else NA_real_
    trans<-if(transform)exp else identity
    data.frame(model=id,contrast=label,inference=conv,n=nrow(des$variables),deaths=sum(des$variables$death),
               coefficient=est,se=se,estimate=trans(est),ci_low=trans(est-crit*se),ci_high=trans(est+crit*se),p=p,
               design_df=df,residual_df=resdf(f),denominator_df=dd,parameters=length(b))
  }))
}
joint <- function(f,des,terms,id,label){
  if(!length(terms) || anyNA(terms) || anyDuplicated(terms) ||
     any(!terms %in% names(coef(f))))
    stop('Joint test requires unique known coefficient names')
  b<-coef(f)[terms];v<-vcov(f)[terms,terms,drop=FALSE];q<-length(b)
  stopifnot(q>0,qr(v)$rank==q)
  W<-drop(t(b)%*%solve(v,b));F<-W/q
  data.frame(model=id,test=label,numerator_df=q,F=F,design_df=degf(des),residual_df=resdf(f),
             p_full=pf(F,q,degf(des),lower.tail=FALSE),
             p_residual=if(resdf(f)>0)pf(F,q,resdf(f),lower.tail=FALSE) else NA_real_,p_normal=pchisq(W,q,lower.tail=FALSE))
}
