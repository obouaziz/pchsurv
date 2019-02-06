#' Simulation of survival data
#'
#' Given a set of cuts and hazard values, generates a survival data sample from a piecewise-constant hazard model
#' @param n the sample size.
#' @param cuts the sequence of cuts. Default to NULL which corresponds to the exponential model.
#' @param alpha the value of the hazard function between each cut. Should be of length equal to length(cuts)+1.
#' @keywords pchsurv
#' @family pchsurv functions
#' @return A time vector of size \code{n} drawn from a piecewise constant hazard model with hazard equal to \code{alpha} between \code{cuts}.
#' @export
#' @examples
#' n=400
#' cuts=c(20,40,50,70)
#' alpha=c(0,0.05,0.1,0.2,0.4)/10
#' time=rsurv(n,cuts,alpha) #generate true data from the pch model
#' censoring=runif(n,min=70,max=90)
#' time=pmin(time,censoring) #observed times
#' delta=time<censoring #gives 62% of observed data on average

rsurv=function(n,cuts=NULL,alpha) {
  if (all(diff(cuts)>0)==FALSE |  all(cuts==cummax(cuts))==FALSE) {
    stop("cuts must be an increasing sequence with no ex-aequo!")}
  u=runif(n);
  k=length(alpha);
  if (length(cuts)!=(k-1)) stop("error: length(cuts) must be equal to length(alpha)-1 !")
  cuts0=c(0,cuts)
  if (k>1) {
    thresh=exp(-cumsum(alpha[-k]*diff(cuts0)))
    if (n<=200) {
      seg=apply(matrix(rep(thresh,n),byrow=T,nrow=n)>u,1,sum)+1
    } else {
      seg=rep(NA,n)
      for (i in 1:n)
        seg[i]=sum(thresh>u[i])+1
    }
  } else {
    seg=rep(1,n)
  }
  res=cuts0[seg]-(log(u)+cumsum(c(0,alpha[-k]*diff(cuts0)))[seg])/alpha[seg]
  return(res)
}
