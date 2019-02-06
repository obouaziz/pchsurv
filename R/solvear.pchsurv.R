#' Penalized maximum likelihood estimation from the piecewise constant hazard model
#'
#' Given time to event data, a set of cuts and a penalty value, returns the penalized maximum likelihood estimator of the log-hazard using the adaptive-ridge penalty.
#' The penalty is iteratively updated until convergence of the algorithm.
#' @param time the observed time data.
#' @param status status of the time data. TRUE for true time and FALSE for censored time.
#' @param cuts a sequence of cuts. Default to NULL which corresponds to the exponential model.
#' @param weights the weight values at the initialization step. Default to 1.
#' @param pen the penalty value
#' @param a the log-hazard value at the initialization step. Default to the unpenalized log-hazard estimator.
#' @param tol the tolerance parameter for convergence of the algorithm. Default to 1e-7.
#' @param itermax the maximum number of iterations. If the algorithm has not converged before itermax iterations then the algorithm exits the program.
#' @param weights an optional weight sequence. Default to 1.
#' @keywords pchsurv
#' @export
#' @examples
#' n=400
#' cuts=c(20,40,50,70)
#' alpha=c(0,0.05,0.1,0.2,0.4)/10
#' time=rsurv(n,cuts,alpha) #generate true data from the pch model
#' censoring=runif(n,min=70,max=90)
#' time=pmin(time,censoring) #observed times
#' delta=time<censoring #gives 62% of observed data on average
#' result1=solvear.pchsurv(time,delta,cuts,pen=2)
#' result2=solvear.pchsurv(time,delta,cuts,pen=0)
#' result3=mlepchsurv(time,delta,cuts) #compare exp(result2) with exp(result3$a)!

solvear.pchsurv=function(time,status,cuts,w=rep(1.0,length(cuts)),pen=2.0,a=NULL,J=NULL,
                          tol=1e-7,itermax=50,weights=rep(1.0,length(time))) {
  k=length(cuts)+1
  if (is.null(J)) J=cut(time,breaks=c(0,cuts,Inf),labels=1:k,right=FALSE)
  # exhaustive statistics
  cuts0=c(0,cuts)
  A=sapply(split(weights*status,J),sum)
  aux=sapply(split(weights,J),sum)[-1]
  B=sapply(split(weights*(time-cuts0[J]),J),sum)
  B[-k]=B[-k]+rev(cumsum(rev(aux)))*diff(cuts0)
  # initialize a (if necessary)
  if (is.null(a)) a=log(A/B)
  a[a==-Inf]=-500
  D=matrix(0,k,2)
  for (iter in 1:itermax) {
    old.a=a
    da=diff(a)
    b=A-exp(a)*B+pen*c(w*da,0)-pen*c(0,w*da)
    #From me#if (sum(is.na(b))==length(b)) { warning("convergence problem for b!");a=NaN;break; }
    #From me#else {
    # D0
    D[,1]=-exp(a)*B-pen*c(w,0)-pen*c(0,w)
    # D1
    D[-k,1+1]=pen*w
    # NR update with bandsolve
    a=a-bandsolve::bandsolve(D,b)
    #if (is.na(max(abs(a-old.a)/abs(a)))) break else {
    if (max(abs(a-old.a)/abs(a))<tol) break;
  }
  return(a)
}
