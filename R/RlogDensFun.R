RlogDensFun <- function(u, theta, delta)
{
  ## The density function
  TC1 <- (1-(1-u)^theta) # # FIXME: Numerically instable if theta -> Inf,  then TC1-> 1
  TC2.log <- (-1+theta)*log(1-u) # TC2 = (1-u)^(-1+theta)

  L5 <- (rowSums(TC1^(-delta)) - 1)
  L6 <- 1-L5^(-1/delta) # FIXME: log(L6)->Inf when u->1,  v->1

  logCplDensObs <- ((-1-delta)*rowSums(log(TC1))+
                      rowSums(TC2.log) -
                      2*(1+delta)/delta*log(L5)+
                       (-2+1/theta)*log(L6)+
                      log(-1+theta+L5^(1/delta)*L6*(1+delta)*theta))

  out.log <- matrix(logCplDensObs)
  return (logCplDensObs)
}
