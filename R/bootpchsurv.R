#' Bootstrap estimator of the adaptive piecewise constant hazard model
#'
#' Given time to event data and a set of cuts return the survival estimator from maximum likelihood estimation in the piecewise constant hazard model
#' @param time the observed time data.
#' @param status status of the time data. TRUE for true time and FALSE for censored time.
#' @param cuts a sequence of cuts. Default to NULL which corresponds to the exponential model.
#' @param seqtime a time sequence where the survival function is estimated.
#' @param M the number of bootrapped samples.
#' @param w the w sequence values for the adaptive ridge algorithm at the initialization step. Default to 1. Should be of the size of the cuts.
#' @param pen a sequence of penalty values.
#' @param a the log-hazard value at the initialization step. Default to the unpenalized log-hazard estimator.
#' @param tol the tolerance parameter for convergence of the algorithm. Default to 1e-7.
#' @param itermax the maximum number of iterations. If the algorithm has not converged before itermax iterations then the algorithm exits the program. Default to 1e+5.
#' @details The bootstrap procedure is performed by sampling on the initial sample. A total of \code{M} bootstrap samples are generated.
#' A new estimator is constructed by taking the median of all survival estimate for all bootrapp samples. Confidence intervals are computed
#' using the bootstrap samples.
#' @return
#' \tabular{lll}{
#' \code{Result} \tab  \code{ } \tab a matrix containing the cumulative hazard estimates for each bootstrap sample\cr
#' \tab \tab \cr
#' \code{medSurv} \tab \code{ }  \tab the median estimator obtained from all bootstrap samples\cr
#' \tab \tab \cr
#' \code{CIleft} \tab  \code{ } \tab the left confidence intervals obtained from the bootstrap samples\cr
#' \tab \tab \cr
#' \code{CIright} \tab \code{ }  \tab the right confidence intervals obtained from the bootstrap samples\cr
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
#' Result=bootpchsurv(time,delta,cuts=1:100,M=20)
#' plot(Result)
#'
#' #The adaptive ridge estimator
#' fit=arpchsurv(time,delta,verbose=TRUE,cuts=1:100,CI=TRUE)
#' fitsurv=pchsurv(time,delta,cuts=fit$final.cuts)
#' lines(fitsurv,CI=TRUE,col="blue") #the adaptive-ridge survival estimator with direct confidence intervals
#' seqtime=seq(0,100,by=0.1)
#' lines(seqtime,exp(-pchcumhaz(seqtime,cuts,alpha)),type="l",col="red") #the true survival function
#'
#'
bootpchsurv <- function(time,status,cuts,seqtime=NULL,M=100,verbose=TRUE,tol=1e-7,itermax=100000,
                        pen=exp(seq(log(0.1),log(1000),length=100)),w=rep(1.0,length(cuts)),
                        a=rep(0,length(cuts)+1)) {UseMethod("bootpchsurv")}
#' @export
bootpchsurv.default=function(time,status,cuts,seqtime=NULL,M=100,verbose=TRUE,tol=1e-7,itermax=100000,
                     pen=exp(seq(log(0.1),log(1000),length=100)),w=rep(1.0,length(cuts)),
                     a=rep(0,length(cuts)+1)){
  n=length(time)
  if (is.null(seqtime)) {seqtime=seq(0,max(time),length.out=100)}
  Result<-matrix(rep(0,length(seqtime)*M),nrow=length(seqtime))
  for (i in 1:M)
  {
    newid<-sample(1:n,n, replace=TRUE)
    newtime<-time[newid];newstatus<-status[newid]
    newfit=arpchsurv(newtime,newstatus,cuts=cuts,pen=pen,verbose=FALSE,CI=FALSE)
    Result[,i]<-pchcumhaz(seqtime,cuts,exp(newfit$res[[newfit$pos]]$exact.a))
    if (verbose==TRUE){
      cat("iter=",i,sep="")
    }
  }
  leftCI<-apply(Result,1,function(z)quantile(z,0.025))
  rightCI<-apply(Result,1,function(z)quantile(z,0.975))
  medSurv<-apply(Result,1,function(z)quantile(z,0.5))
  out<-list(time=time,status=status,cuts=cuts,Result=Result,seqtime=seqtime,leftCI=leftCI,rightCI=rightCI,medSurv=medSurv,call=match.call())
  class(out)<-"bootpch"
  out
}
#' @export
print.bootpch <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nSurvival Estimate:\n\n")
  matres=data.frame(round(x$seqtime,2),round(x$medSurv,4),round(x$leftCI,4),round(x$rightCI,4))
  colnames(matres)<-c("time","survival","low CI", "up CI")
  print(matres, row.names = FALSE)
}
#' @export
plot.bootpch<-function(bootfit,...)
{
  plot(c(0,max(bootfit$time)),c(0,1),t='n',xlab="Time",ylab="Survival estimate",...)
  lines(bootfit$seqtime,exp(-bootfit$medSurv),type="l",...)
  lines(bootfit$seqtime,exp(-bootfit$leftCI),type="l",lty=4,...)
  lines(bootfit$seqtime,exp(-bootfit$rightCI),type="l",lty=4,...)
}
