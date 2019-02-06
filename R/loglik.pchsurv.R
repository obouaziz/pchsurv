#' Likelihood function of a piecewise constant hazard model
#'
#' Given time to event data, a set of cuts and a log-hazard value, compute the maximum likelihood function
#' @param time the observed time data.
#' @param status status of the time data. TRUE for true time and FALSE for censored time.
#' @param cuts a sequence of cuts. Default to NULL which corresponds to the exponential model.
#' @param a sequence of log-hazard values. Should be of length equal to length(cuts)+1.
#' @param weights an optional weight sequence. Default to NULL.
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
#' result=mlepchsurv(time,delta,cuts)
#' loglik.pchsurv(time,delta,cuts,result$a)


loglik.pchsurv=function(time,status,cuts,a,weights=NULL) {
  if (all(diff(cuts)>0)==FALSE |  all(cuts==cummax(cuts))==FALSE) {
    stop("cuts must be an increasing sequence with no ex-aequo!")}
  if (is.null(weights)) weights=rep(1,length(time))
  k=length(a)
  if (length(cuts)!=(k-1)) stop("error: length(cuts) must be equal to length(alpha)-1 !")
  J=cut(time,breaks=c(0,cuts,Inf),labels=1:k,right=FALSE)
  cuts0=c(0,cuts)
  cumhaz=exp(a)[J]*(time-cuts0[J])+c(0,cumsum(exp(a)*c(diff(cuts0),0))[-k])[J]
  return(sum(weights*((status==1)*a[J]-cumhaz)))
}
