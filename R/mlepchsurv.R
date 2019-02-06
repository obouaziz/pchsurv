#' Maximum likelihood estimation from the piecewise constant hazard model
#'
#' Given time to event data and a set of cuts return the maximum likelihood estimator of the log-hazard from the piecewise constant hazard model.
#' @param time observed time data.
#' @param status status of the time data. TRUE for true time and FALSE for censored time.
#' @param cuts the sequence of cuts. Default to NULL which corresponds to the exponential model.
#' @param weights an optional weight sequence. Default to NULL.
#' @param CI should the confidence intervals be computed? Default to TRUE.
#' @param alphaCI the value of alpha for 1-alpha confidence intervals. Default to 0.05.
#' @param logtransf should the confidence intervals be computed using the log-transformtation? Default to TRUE.
#' @details The maximum likelihood estimator is computed from the two exhaustive statistics O_k=sum_i Delta_i 1(c_{(k-1)}<T_i<=c_k) and R_k=sum_i (min(c_k,T_i)-c_{(k-1)}) 1(T_i>=c_{(k-1)}).
#' It is equal to O_k/R_k on each cut (c_{(k-1)},c_{k}]. They are called A and B in what follows.
#' @return
#' \tabular{lll}{
#' \code{a} \tab  \code{ } \tab the estimated log-hazard\cr
#' \tab \tab \cr
#' \code{hazard} \tab \code{ }  \tab the estimated hazard\cr
#' \tab \tab \cr
#' \code{A} \tab  \code{ } \tab the number of observed events between each cut\cr
#' \tab \tab \cr
#' \code{B} \tab \code{ }  \tab the total time at risk between each cut\cr
#' \tab \tab \cr
#' \code{CIleft} \tab  \code{ } \tab the left confidence intervals\cr
#' \tab \tab \cr
#' \code{CIright} \tab \code{ }  \tab the right confidence intervals\cr
#' }
#'
#' @keywords pchsurv
#' @family pchsurv functions
#' @export
#' @examples
#' n=400
#' cuts=c(20,40,50,70)
#' alpha=c(0,0.05,0.1,0.2,0.4)/10
#' time=rsurv(n,cuts,alpha) #generate true data from the pch model
#' censoring=runif(n,min=70,max=90)
#' time=pmin(time,censoring) #observed times
#' delta=time<censoring #gives 62% of observed data on average
#' mlepchsurv(time,delta,cuts)

mlepchsurv <- function(time,status,cuts=NULL,weights=NULL,CI=TRUE,alphaCI=0.05,logtransf=TRUE) {UseMethod("mlepchsurv")}
#' @export
mlepchsurv.default=function(time,status,cuts=NULL,weights=NULL,CI=TRUE,alphaCI=0.05,logtransf=TRUE) {
  if(length(time)!=length(status)) {stop("time and status should have same length!")}
  if (all(diff(cuts)>0)==FALSE |  all(cuts==cummax(cuts))==FALSE) {
    stop("cuts must be an increasing sequence with no ex-aequo!")}
  if (is.null(weights)) {weights=rep(1,length(time))}
  k=length(cuts)+1
  J=cut(time,breaks=c(0,cuts,Inf),labels=1:k,right=FALSE)
  cuts0=c(0,cuts)
  A=sapply(split(weights*status,J),sum)
  aux=sapply(split(weights,J),sum)[-1]
  B=sapply(split(weights*(time-cuts0[J]),J),sum)
  B[-k]=B[-k]+rev(cumsum(rev(aux)))*diff(cuts0)
  a=log(A/B)
  hazard=NULL
  if (length(cuts)>0) hazard=stepfun(cuts,exp(a))
  CIleft=CIright=rep(NA,k)
  if (CI==TRUE){
  if (sum(a==-Inf)!=0){CIleft[a==-Inf]<-0;CIright[a==-Inf]<-0}
  if (logtransf==TRUE){
    CIleft[a!=-Inf]=exp(a[a!=-Inf]-qnorm(1-alphaCI/2)/sqrt(A[a!=-Inf]))#log-transform CI
    CIright[a!=-Inf]=exp(a[a!=-Inf]+qnorm(1-alphaCI/2)/sqrt(A[a!=-Inf]))#log-transform CI
  } else {
    CIleft[a!=-Inf]=exp(a[a!=-Inf])-qnorm(1-alphaCI/2)*sqrt(A[a!=-Inf])/B
    CIright[a!=-Inf]=exp(a[a!=-Inf])+qnorm(1-alphaCI/2)*sqrt(A[a!=-Inf])/B
  }
  if (sum(CIright==Inf)!=0) {warning("some of the CIs could not be computed!")}
  }
  out<-list(a=a,hazard=hazard,A=A,B=B,CIleft=CIleft,CIright=CIright,Range=range(time),cuts=cuts,call=match.call())
  class(out)<-"pch"
  out
}
#' @export
print.pch <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nMaximum Likelihood Estimate:\n\n")
  if (sum(is.na(x$CIleft))==length(x$a)){
    matres=data.frame(convertPCH(x$cuts),cbind(round(x$a,4),round(exp(x$a),4)))
    colnames(matres)<-c("cuts","log-haz","haz")
  } else {
    matres=data.frame(convertPCH(x$cuts),cbind(round(x$a,4),round(exp(x$a),4),round(x$CIleft,4),round(x$CIright,4)))
    colnames(matres)<-c("cuts","log-haz","haz","low CI", "up CI")
  }
  print(matres, row.names = FALSE)
}
#' @export
plot.pch=function(fit,CI=FALSE,xlab="Time",ylab="Hazard",main="Hazard function",lwd=1,...)
{
  Range=c(0,fit$Range[2])
  # else {
  #   RangeL=Range[1];RangeR=Range[2]
  # }
  # seqtime=c(RangeL,fit$cuts,RangeR)
  #seq(RangeL,RangeR,by=0.5)
  if (CI==TRUE)
  {
    hazardL=stepfun(fit$cuts,fit$CIleft)
    hazardR=stepfun(fit$cuts,fit$CIright)
    plot(Range,range(fit$hazard(c(0,max(fit$cuts,Range[2]))),hazardL(c(0,fit$cuts,max(fit$cuts,Range[2]))),hazardR(c(0,fit$cuts,max(fit$cuts,Range[2])))),t='n',xlab=xlab,ylab=ylab,main=main,...)
    lines(c(0,fit$cuts,max(fit$cuts,Range[2])),fit$hazard(c(0,fit$cuts,max(fit$cuts,Range[2]))),type="s",lwd=lwd,...)
    lines(c(0,fit$cuts,max(fit$cuts,Range[2])),hazardL(c(0,fit$cuts,max(fit$cuts,Range[2]))),type="s",lwd=lwd,lty=2,...)
    lines(c(0,fit$cuts,max(fit$cuts,Range[2])),hazardR(c(0,fit$cuts,max(fit$cuts,Range[2]))),type="s",lwd=lwd,lty=2,...)
  } else {
    #plot(c(0,fit$cuts,max(fit$cuts,fit$time)),hazard(c(0,fit$cuts,max(fit$cuts,fit$time))),type="s",lwd=lwd,xlab=xlab,ylab=ylab,main=main,...)
    plot(Range,range(fit$hazard(c(0,fit$cuts,max(fit$cuts,Range[2])))),t='n',xlab=xlab,ylab=ylab,main=main,...)
    lines(c(0,fit$cuts,max(fit$cuts,Range[2])),fit$hazard(c(0,fit$cuts,max(fit$cuts,Range[2]))),type="s",lwd=lwd,...)
  }
}
#' @export
lines.pch=function(fit,Range=NA,CI=FALSE,lwd=1,...)
{
  if (is.na(Range)==TRUE){
    Range=c(0,fit$Range[2])
  }# else {
  #   RangeL=Range[1];RangeR=Range[2]
  # }
  # seqtime=c(RangeL,fit$cuts,RangeR)
  #seq(RangeL,RangeR,by=0.5)
  if (CI==TRUE)
  {
    hazardL=stepfun(fit$cuts,fit$CIleft)
    hazardR=stepfun(fit$cuts,fit$CIright)
    lines(c(0,fit$cuts,Range[2]),fit$hazard(c(0,fit$cuts,Range[2])),type="s",lwd=lwd,...)
    lines(c(0,fit$cuts,Range[2]),hazardL(c(0,fit$cuts,Range[2])),type="s",lwd=lwd,lty=2,...)
    lines(c(0,fit$cuts,Range[2]),hazardR(c(0,fit$cuts,Range[2])),type="s",lwd=lwd,lty=2,...)
  } else {
    lines(c(0,fit$cuts,Range[2]),fit$hazard(c(0,fit$cuts,Range[2])),type="s",lwd=lwd,...)
  }
}
#' @export
convertPCH<-function(cuts){
  #cuts0=c(0,cuts)
  k=length(cuts)
  newcuts=0
  vect=rep(NA,k)
  for (i in 1:k)
  {
    vect[i]=paste("[",newcuts,",",cuts[i],")",sep="")
    newcuts=cuts[i]
  }
  vect=c(vect,paste("[",newcuts,",","Inf",")",sep=""))
  vect
}

