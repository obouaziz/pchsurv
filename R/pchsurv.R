#' Estimated survival function from the piecewise constant hazard model
#'
#' Given time to event data and a set of cuts return the survival estimator from maximum likelihood estimation in the piecewise constant hazard model
#' @param time the observed time data.
#' @param status status of the time data. TRUE for true time and FALSE for censored time.
#' @param cuts a sequence of cuts. Default to NULL which corresponds to the exponential model.
#' @param seqtime a time sequence.
#' @param CI should the confidence intervals be computed? Default to TRUE.
#' @param alphaCI the value of alpha for 1-alpha confidence intervals. Default to 0.05.
#' @details The survival estimator is computed from the maximum likelihood estimator in the piecewise constant hazard model. It is computed for each value of \code{seqtime}.
#' @return
#' \tabular{lll}{
#' \code{surv} \tab  \code{ } \tab the estimated survival function\cr
#' \tab \tab \cr
#' \code{CIleft} \tab  \code{ } \tab the left confidence intervals of the survival function\cr
#' \tab \tab \cr
#' \code{CIright} \tab \code{ }  \tab the right confidence intervals of the survival function\cr
#' }
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
#' seqtime=seq(0,100,by=10)
#' result=pchsurv(time,delta,cuts,seqtime)
#' result
#' plot(result,CI=TRUE)

pchsurv <- function(time,status,cuts,seqtime=NULL,CI=TRUE,alphaCI=0.05) {UseMethod("pchsurv")}
#' @export
pchsurv.default=function(time,status,cuts,seqtime=NULL,CI=TRUE,alphaCI=0.05) {
  if(length(time)!=length(status)) {stop("time and status should have same length!")}
  if (all(diff(cuts)>0)==FALSE |  all(cuts==cummax(cuts))==FALSE) {
    stop("cuts must be an increasing sequence with no ex-aequo!")}
  if (is.null(seqtime)==TRUE) {seqtime=seq(0,max(time),length.out=1000)}
  estime=mlepchsurv(time,status,cuts)
  alpha=exp(estime$a)
  A=estime$A
  if (length(cuts)==0) {
    surv=exp(-seqtime*exp(estime$a))
    CIleft=(surv)^(exp(qnorm(1-alphaCI/2)/sqrt(A)))
    CIright=(surv)^(exp(-qnorm(1-alphaCI/2)/sqrt(A)))
  } else {
    k=length(cuts)+1
    I=k-apply(matrix(rep(seqtime,k-1),byrow=TRUE,nrow=k-1)<cuts,2,sum)
    cuts0=c(0,cuts)
    cumhaz=alpha[I]*(seqtime-cuts0[I])+c(0,cumsum(alpha*c(diff(cuts0),0))[-k])[I]
    Abis=A;Abis[alpha==0]=1#to avoid dividing by 0
    cumhaz2=(alpha[I]*(seqtime-cuts0[I]))^2/Abis[I]+c(0,cumsum((alpha*c(diff(cuts0),0))^2/Abis)[-k])[I]
    CIterm=exp(log(cumhaz2)-2*log(cumhaz))
    surv=exp(-cumhaz)
    CIleft=(surv)^(exp(qnorm(1-alphaCI/2)*sqrt(CIterm)))
    CIright=(surv)^(exp(-qnorm(1-alphaCI/2)*sqrt(CIterm)))
  }
  out<-list(surv=surv,seqtime=seqtime,CIleft=CIleft,CIright=CIright,call=match.call())
  class(out)<-"surv"
  out
}
#' @export
print.surv <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
    cat("\nSurvival Estimate:\n\n")
    if (is.null(x$CIleft)==TRUE & is.null(x$CIright)==TRUE){
      matres=data.frame(round(x$seqtime,2),round(x$surv,4))
      colnames(matres)<-c("time","survival")
    } else {
      matres=data.frame(round(x$seqtime,2),round(x$surv,4),round(x$CIleft,4),round(x$CIright,4))
      colnames(matres)<-c("time","survival","low CI", "up CI")
    }
    print(matres, row.names = FALSE)
}
#' @export
plot.surv=function(fit,CI=FALSE,xlab="Time",ylab="Survival",main="Survival function",lwd=1,...)
{
  if (CI==TRUE)
  {
    plot(range(fit$seqtime),range(c(fit$CIleft[is.finite(fit$CIleft)],fit$CIright[is.finite(fit$CIright)])),t='n',xlab=xlab,ylab=ylab,main=main,...)
    lines(fit$seqtime,fit$CIleft,type="l",lty=2,lwd=lwd,...)
    lines(fit$seqtime,fit$CIright,type="l",lty=2,lwd=lwd,...)
    lines(fit$seqtime,fit$surv,type="l",lwd=lwd,...)
  } else {
    plot(fit$seqtime,fit$surv,type="l",xlab=xlab,ylab=ylab,main=main,lwd=lwd,...)
  }
}
#' @export
lines.surv=function(fit,CI=FALSE,lwd=1,...)
{
  lines(fit$seqtime,fit$surv,type="l",lwd=lwd,...)
  if (CI==TRUE)
  {
    lines(fit$seqtime,fit$CIleft,type="l",lwd=lwd,lty=2,...)
    lines(fit$seqtime,fit$CIright,type="l",lwd=lwd,lty=2,...)
  }
}

