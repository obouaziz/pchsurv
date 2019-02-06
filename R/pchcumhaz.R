#' Computation of the cumulative hazard function
#'
#' Given a sequence of time, a set of cuts and hazard values, compute the cumulative hazard function from a piecewise-constant hazard model
#' @param seqtime the time sequence.
#' @param cuts the sequence of cuts.
#' @param alpha the value of the hazard function between each cut. Should be of length equal to length(cuts)+1.
#' @keywords pchsurv
#' @export
#' @examples
#' cuts=c(20,40,50,70)
#' alpha=c(0,0.05,0.1,0.2,0.4)/10
#' seqtime=seq(0,100,by=10)
#' pchcumhaz(seqtime,cuts,alpha)
pchcumhaz=function(seqtime,cuts,alpha) {
  if (length(cuts)==0) return(time*alpha[1])
  k=length(alpha)
  I=k-apply(matrix(rep(seqtime,k-1),byrow=TRUE,nrow=k-1)<cuts,2,sum)
  cuts0=c(0,cuts)
  return(alpha[I]*(seqtime-cuts0[I])+c(0,cumsum(alpha*c(diff(cuts0),0))[-k])[I])
}
