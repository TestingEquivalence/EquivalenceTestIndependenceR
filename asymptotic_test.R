source("distance.R")

asympt_stdev<-function(p,derivative){
  vec = derivative
  vnsq_1  = sum(p*vec*vec)
  
  k=length(p)
  vnsq_2=0
  for (j1 in 1:k)
    for (j2 in 1:k)
      vnsq_2 = vnsq_2 + vec[j1] * vec[j2] * p[j1] * p[j2]
  
  
  vnsq  = vnsq_1 - vnsq_2
  return (sqrt(vnsq))
}


#' The asymptotic test is based on the asymptotic distribution of the test statistic. 
#' The test statistic is scaled Euclidian distance between the counting frequencies
#' and the product measure of the marginal distributions.
#' Therefore the asymptotic test need some sufficiently large number of the observations
#' in any cell of the contingency table.
#' It should be used carefully because the test is approximate 
#' and may be anti conservative at some points. 
#' In order to obtain a conservative test reducing of alpha  (usually halving) or
#' slight shrinkage of the tolerance parameter epsilon may be appropriate. 
#' \code{asymptotic_test_absolute} asymptotic test for approximate row column independence
#' in two way contingency tables. 
#' The test statistic is scaled Euclidian distance between the counting frequencies
#' and the product measure of the marginal distributions.
#' @param tab contingency table containing the counts of events
#' @param alpha significance level
#' @return tests returns the minimum tolerance parameter epsilon,
#' for which the approximate independence can be shown

asymptotic_test_absolute<-function(tab, alpha){
  n=sum(tab)
  tab=tab/n
  vtab=as.vector(t(tab))
  der=derivative_distance_abs(tab)
  vder=as.vector(t(der))
  
  vol = asympt_stdev(vtab,vder) / sqrt(n)
  qt=qnorm(1-alpha,0,1)
  t= distance_abs(tab)
  eps = t + qt*vol
  eps=sqrt(eps)
  return(eps)
}

#' The asymptotic test is based on the asymptotic distribution of the test statistic. 
#' The test statistic is scaled Euclidian norm of the relative deviations 
#' between the counting frequencies and the product measure of the marginal distributions.
#' Therefore the asymptotic test need some sufficiently large number of the observations
#' in any cell of the contingency table.
#' It should be used carefully because the test is approximate 
#' and may be anti conservative at some points. 
#' In order to obtain a conservative test reducing of alpha  (usually halving) or
#' slight shrinkage of the tolerance parameter epsilon may be appropriate. 
#' \code{asymptotic_test_relative} asymptotic test for approximate row column independence
#' in two way contingency tables. 
#' The test statistic is scaled Euclidian norm of the relative deviations 
#' between the counting frequencies and the product measure of the marginal distributions.
#' @param tab contingency table containing the counts of events
#' @param alpha significance level
#' @return tests returns the minimum tolerance parameter epsilon,
#' for which the approximate independence can be shown

asymptotic_test_relative<-function(tab, alpha){
  n=sum(tab)
  tab=tab/n
  vtab=as.vector(t(tab))
  der=derivative_distance_rel(tab)
  vder=as.vector(t(der))
  
  vol = asympt_stdev(vtab,vder) / sqrt(n)
  qt=qnorm(1-alpha,0,1)
  t= distance_rel(tab)
  eps = t + qt*vol
  eps=sqrt(eps/(nrow(tab)*ncol(tab)))
  return(eps)
}