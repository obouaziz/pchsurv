#' The adaptive ridge and ridge estimates from the piecewise constant hazard model
#'
#' Given time to event data, a set of cuts and a sequence of penalty values, returns the adaptive ridge or the ridge estimator of the log-hazard.
#' @param time the observed time data.
#' @param status status of the time data. TRUE for true time and FALSE for censored time.
#' @param cuts a sequence of cuts. Default to NULL which corresponds to the exponential model.
#' @param w the w sequence values for the adaptive ridge algorithm at the initialization step. Default to 1. Should be of the size of the cuts.
#' @param pen a sequence of penalty values.
#' @param a the log-hazard value at the initialization step. Default to the unpenalized log-hazard estimator.
#' @param tol the tolerance parameter for convergence of the algorithm. Default to 1e-7.
#' @param itermax the maximum number of iterations. If the algorithm has not converged before itermax iterations then the algorithm exits the program. Default to 1e+5.
#' @param weights an optional weight sequence. Default to 1. Should be of the size of time and status.
#' @param AR do you want to compute the adaptive-ridge? If FALSE compute the ridge estimator. Default to TRUE.
#' @param CI should the confidence intervals be computed? Default to TRUE.
#' @param alphaCI the value of alpha for 1-alpha confidence intervals. Default to 0.05.
#' @param logtransf should the confidence intervals be computed using the log-transformtation? Default to TRUE.
#' @details The adaptive-ridge algorithm is run for the sequence of penalty values \code{pen}. The final adaptive-ridge estimator \code{final.a} is found
#' by minimising a Bayesian Information Criteria (BIC). The \code{res} list returned from the function contains all the
#' estimators (adaptive-ridge or ridge) for each penalty value. Write \code{res[[pos]]} to get the estimator for the
#' penalty corresponding to \code{pen[pos]}. Use the option AR to choose between the ridge or adaptive-ridge estimator.
#' For the ridge estimator only res is not null.Confidence intervals are computed for the final adaptive-ridge hazard estimator.
#' @return
#' \tabular{lll}{
#' \code{res} \tab  \code{ } \tab the list of adaptive-ridge estimators for each penalty value\cr
#' \tab \tab \cr
#' \code{final.a} \tab \code{ }  \tab the final log-hazard that minimizes the BIC\cr
#' \tab \tab \cr
#' \code{final.cuts} \tab  \code{ } \tab the final cut vector found from the BIC\cr
#' \tab \tab \cr
#' \code{bic} \tab \code{ }  \tab the minimal value of the BIC\cr
#' \tab \tab \cr
#' \code{pos} \tab  \code{ } \tab the position of the penalty that minimizes the BIC\cr
#' \tab \tab \cr
#' \code{CIleft} \tab  \code{ } \tab the left confidence intervals for the final hazard\cr
#' \tab \tab \cr
#' \code{CIright} \tab \code{ }  \tab the right confidence intervals for the final hazard\cr
#' }
#' @keywords pchsurv
#' @family pchsurv functions
#' @export
#' @examples
#' set.seed(45)
#' n=400
#' cuts=c(20,40,50,70)
#' alpha=c(0,0.05,0.1,0.2,0.4)/10
#' time=rsurv(n,cuts,alpha) #generate true data from the pch model
#' censoring=runif(n,min=70,max=90)
#' time=pmin(time,censoring) #observed times
#' delta=time<censoring #gives 62% of observed data on average
#'
#' #The adaptive ridge estimator
#' fit=arpchsurv(time,delta,verbose=TRUE,cuts=1:100,CI=TRUE)
#'
#' fit$final.cuts
#' exp(fit$final.a)
#'
#' plot(fit)
#' plot(fit,pos=fit$pos,CI=TRUE)
#' points(c(0,c(20,40,50,70),max(time)),c(0,alpha),type="S",lwd=2,xlim=c(10,80),lty=2,col="red") #the true hazard function
#'
#' fitsurv=pchsurv(time,delta,cuts=fit$final.cuts)
#' plot(fitsurv,CI=TRUE)
#' seqtime=seq(0,100,by=0.1)
#' lines(seqtime,exp(-pchcumhaz(seqtime,cuts,alpha)),type="l",col="red") #the true survival function
#'
#' #The ridge estimator
#'fitR=arpchsurv(time,delta,verbose=TRUE,cuts=1:100,AR=FALSE)
#'plot(fitR)
#'plot(fitR,pos=50,main="")
#'lines(fitR,pos=80,col="blue")
#'points(c(0,c(20,40,50,70),max(time)),c(0,alpha),type="S",lwd=2,xlim=c(10,80),lty=2,col="red") #the true hazard function

arpchsurv <- function(time,status,cuts=NULL,weights=rep(1.0,length(time)),pen=exp(seq(log(0.1),log(1000),length=100)),a=rep(0,length(cuts)+1),AR=TRUE,...) {UseMethod("arpchsurv")}
#' @export
arpchsurv.default=function(time,status,cuts=NULL,weights=rep(1.0,length(time)),pen=exp(seq(log(0.1),log(1000),length=100)),
                   verbose=TRUE,beta=1e-5,a=rep(0,length(cuts)+1),tol=1e-7,itermax=100000,
                    w=rep(1.0,length(cuts)),
                    AR=TRUE,CI=TRUE,alphaCI=0.05,logtransf=TRUE) {
  if(length(time)!=length(status)) {stop("time and status should have same length!")}
  # initialization
  if (is.null(cuts)) cuts=seq(0,max(time[status==1],length.out=200))#unique(sort(time[status]))
  if (all(diff(cuts)>0)==FALSE |  all(cuts==cummax(cuts))==FALSE) {
    stop("cuts must be an increasing sequence with no ex-aequo!")}
  k=length(cuts)+1
  J=cut(time,breaks=c(0,cuts,Inf),labels=1:k,right=FALSE)
  pos=1
  res=vector(length(pen),mode="list")
  CIleft=CIright=NULL
  # main loop
  if (AR==TRUE){
    for (iter in 1:itermax) {
      old.a=a;
      a=solvear.pchsurv(time,status,cuts,w,pen[pos],a,J,weights=weights)
      # update weights
      w=1/(diff(a)^2+beta^2)
      # test for convergence
      if ((max(abs(a-old.a)/abs(a))<tol)|(max(abs(a))<tol)) {
        # save model
        sel=w*diff(a)^2
        fit=mlepchsurv(time,status,cuts[sel>0.99],weights)
        if (is.null(fit$hazard)) {
          exact.a=rep(fit$a,length(a))
        } else {
          exact.a=log(fit$hazard(c(0,cuts)))
        }
        bic=-2*loglik.pchsurv(time,status,cuts,exact.a)+(sum(sel>0.99)+1)*log(sum(weights))
        res[[pos]]=list(sel=sel,w=w,a=a,fit=fit,bic=bic,exact.a=exact.a)
        # verbose
        if (verbose) cat("iter=",iter," log(pen)=",signif(log(pen[pos]),3)," dim=",sum(sel>0.99),"\n",sep="")
        # next penalty
        pos=pos+1;
      }
      if (pos>length(pen)) break
    }
    opt.bic=sapply(res,function(z)z$bic)#the bic value for each penalty
    final.pos=which.min(opt.bic)
    final.a=res[[final.pos]]$fit$a
    #final.a=res[[final.pos]]$exact.a[res[[final.pos]]$sel>0.99]
    final.cuts=res[[final.pos]]$fit$cuts#cuts[#sel>0.99]
    if (CI==TRUE)
    {
      CIleft=CIright=rep(NA,length(final.a))
      if (sum(final.a==-Inf)!=0){CIleft[final.a==-Inf]<-0;CIright[final.a==-Inf]<-0}
      A=res[[final.pos]]$fit$A
      B=res[[final.pos]]$fit$B
      if (logtransf==TRUE){
        CIleft[final.a!=-Inf]=exp(final.a[final.a!=-Inf]-qnorm(1-alphaCI/2)/sqrt(A[final.a!=-Inf]))#log-transform CI
        CIright[final.a!=-Inf]=exp(final.a[final.a!=-Inf]+qnorm(1-alphaCI/2)/sqrt(A[final.a!=-Inf]))#log-transform CI
      } else {
        CIleft[final.a!=-Inf]=exp(final.a[final.a!=-Inf])-qnorm(1-alphaCI/2)*sqrt(A[final.a!=-Inf])/B[final.a!=-Inf]
        CIright[final.a!=-Inf]=exp(final.a[final.a!=-Inf])+qnorm(1-alphaCI/2)*sqrt(A[final.a!=-Inf])/B[final.a!=-Inf]
      }
    }
  } else {
    for (iter in 1:itermax) {
      old.a=a;
      a=solvear.pchsurv(time,status,cuts,w,pen[pos],a,J,weights=weights)
      res[[pos]]=a
      # verbose
      if (verbose) cat("iter=",iter," log(pen)=",signif(log(pen[pos]),3),"\n",sep="")
      # next penalty
      pos=pos+1;
      if (pos>length(pen)) break
    }
    opt.bic=NULL
    final.pos=NULL
    final.a=NULL#this function does not select a final estimator for the ridge! A CV-criterion needs to be implemented next.
    final.cuts=NULL
  }
  # break if all penalties have been processed

  out<-list(time=time,status=status,pen=pen,res=res,iter=iter,cuts=cuts,AR=AR,final.a=final.a,final.cuts=final.cuts,bic=opt.bic,pos=final.pos,CIleft=CIleft,CIright=CIright,call=match.call())
  class(out)<-"arpchsurv"
  out
}
#' @export
print.arpchsurv <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  if (x$AR==TRUE){
  cat("\nCuts:\n")
  print(x$final.cuts)
  cat("\nAdaptive Ridge Estimate:\n\n")
  if (is.null(x$CIleft)==TRUE & is.null(x$CIright)==TRUE){
    matres=data.frame(convertPCH(x$final.cuts),cbind(round(x$final.a,4),round(exp(x$final.a),4)))
    colnames(matres)<-c("cuts","log-haz","haz")
  } else {
    matres=data.frame(convertPCH(x$final.cuts),cbind(round(x$final.a,4),round(exp(x$final.a),4),round(x$CIleft,4),round(x$CIright,4)))
    colnames(matres)<-c("cuts","log-haz","haz","low CI", "up CI")
  }
  print(matres, row.names = FALSE)} else {
    cat("There is no print for the ridge estimator, use $res[[pos]] to display the estimator for a given penalty")
  }
}
#' @export
plot.arpchsurv=function(fit,CI=FALSE,pos=NULL,xlab="Time",ylab="Hazard",main="pen",lwd=1,COL=FALSE,seqtime=NULL,...) {
  if (is.null(seqtime)){seqtime=seq(0,max(fit$time),length.out=1000)}
  if (is.null(pos)) {
    if(fit$AR==TRUE){
      #for the adaptive ridge
      a=sapply(fit$res,function(z)z$exact.a)
      a=matrix(a,ncol=length(fit$pen))
      plot(range(fit$pen),range(a[is.finite(a)]),log="x",t='n',xlab="Penalty",ylab="Parameter value (on the log-scale)",...)
      a[a==-Inf]<-min(a[is.finite(a)])-20
      for (j in 1:nrow(a)){
        if (COL==TRUE){
          lines(fit$pen,a[j,],col=j,lwd=2)#,lty=j
        } else {lines(fit$pen,a[j,],lty=j,lwd=2)}
      }
      bic=sapply(fit$res,function(z)z$bic);
      abline(v=fit$pen[which.min(bic)],lty=2)
    } else {
      #for the ridge
      a=matrix(unlist(fit$res),ncol=length(fit$pen))
      plot(range(fit$pen),range(a[is.finite(a)]),log="x",t='n',xlab="Penalty",ylab="Parameter value (on the log-scale)")
      for (j in 1:nrow(a))
        if (COL==TRUE){
          lines(fit$pen,a[j,],col=j,lwd=2)#,col=j
        } else {lines(fit$pen,a[j,],lty=j,lwd=2)}
    }
  } else {
    if (pos>length(fit$pen) | pos<=0) {stop("wrong value for pos!")}
    if (fit$AR==TRUE){
      hazard=stepfun(fit$cuts,exp(fit$res[[pos]]$exact.a))#for the adaptive ridge#stepfun(fit$cuts,exp(fit$res[[pos]]$exact.a))#for the adaptive ridge
      #time=seq(0,max(fit$time),length=100)
      if (main=="pen"){main=paste0("log(pen)=",signif(log(fit$pen[pos]),3))}
      if (CI==TRUE & pos==fit$pos){
        hazardL=stepfun(fit$final.cuts,(fit$CIleft))
        hazardR=stepfun(fit$final.cuts,(fit$CIright))
        plot(c(0,range(fit$time)[2]),range(c(hazardL(fit$time)[is.finite(hazardL(fit$time))],hazardR(fit$time)[is.finite(hazardR(fit$time))])),t='n',xlab=xlab,ylab=ylab,main=main,...)
        lines(c(0,fit$final.cuts,max(fit$final.cuts,fit$time)),hazardL(c(0,fit$final.cuts,max(fit$final.cuts,fit$time))),type="s",lty=2,...)
        lines(c(0,fit$final.cuts,max(fit$final.cuts,fit$time)),hazardR(c(0,fit$final.cuts,max(fit$final.cuts,fit$time))),type="s",lty=2,...)
        lines(c(0,fit$final.cuts,max(fit$final.cuts,fit$time)),hazard(c(0,fit$final.cuts,max(fit$final.cuts,fit$time))),type="s",...)
        #lines(time,hazardL(time),type="S",lty=2,...)
        #lines(time,hazardR(time),type="S",lty=2,...)
        #lines(time,hazard(time),type="S",lwd=lwd,...)
      } else {
        plot(c(0,fit$cuts,max(fit$cuts,fit$time)),hazard(c(0,fit$cuts,max(fit$cuts,fit$time))),type="s",lwd=lwd,xlab=xlab,ylab=ylab,main=main,...)
      }
    }
    else {
      if (main=="pen"){main=paste0("log(pen)=",signif(log(fit$pen[pos]),3))}
      hazard=stepfun(fit$cuts,exp(fit$res[[pos]]))#for the ridge
      plot(seqtime,hazard(seqtime),type="s",lwd=lwd,xlab=xlab,ylab=ylab,main=main,...)
    }
    #if (CI==TRUE & fit$AR==TRUE & pos==fit$pos){
    #  hazardL=stepfun(fit$final.cuts,(fit$CIleft))
    #  hazardR=stepfun(fit$final.cuts,(fit$CIright))
    #  lines(time,hazardL(time),type="S",lty=2,...)
    #  lines(time,hazardR(time),type="S",lty=2,...)
    #}
  }
}
#' @export
lines.arpchsurv=function(fit,pos=NULL,lwd=1,CI=FALSE,seqtime=NULL,...) {
  if (is.null(seqtime)){seqtime=seq(0,max(fit$time),length.out=1000)}
  if (is.null(pos)) {
    stop("position for the penalty must be specified!")
  } else {
    if (pos>length(fit$pen) | pos<=0) {stop("wrong value for pos!")}
    if (fit$AR==TRUE){
      hazard=stepfun(fit$cuts,exp(fit$res[[pos]]$exact.a))#for the adaptive ridge#stepfun(fit$cuts,exp(fit$res[[pos]]$exact.a))#for the adaptive ridge
      #time=seq(0,max(fit$time),length=100)
      lines(c(0,fit$cuts,max(fit$cuts,fit$time)),hazard(c(0,fit$cuts,max(fit$cuts,fit$time))),type="s",lwd=lwd,...)
      #lines(time,hazard(time),type="S",lwd=lwd,...)
      if (CI==TRUE & pos==fit$pos){
        hazardL=stepfun(fit$final.cuts,(fit$CIleft))
        hazardR=stepfun(fit$final.cuts,(fit$CIright))
        lines(c(0,fit$final.cuts,max(fit$time)),hazardL(c(0,fit$final.cuts,max(fit$time))),type="s",lty=2,lwd=lwd,...)
        lines(c(0,fit$final.cuts,max(fit$time)),hazardR(c(0,fit$final.cuts,max(fit$time))),type="s",lty=2,lwd=lwd,...)
        #lines(time,hazardL(time),type="s",lty=2,...)
        #lines(time,hazardR(time),type="s",lty=2,...)
      }
    } else {
      hazard=stepfun(fit$cuts,exp(fit$res[[pos]]))#for the ridge
      lines(seqtime,hazard(seqtime),type="s",lwd=lwd,...)
    }
  }
}

