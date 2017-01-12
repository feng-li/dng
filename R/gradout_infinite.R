gradFun4delta_infinite <- function(u, delta, theta, redo.idx)
{
  require("Rmpfr", quietly = TRUE)
  precBits <- 1024
  ## MPFR class used for u, theta,  delta
  gradout.redoMPFR <- gradFun4delta(u = mpfr(u[redo.idx, , drop = FALSE], precBits = precBits),
                                    theta = mpfr(theta[redo.idx], precBits = precBits),
                                    delta = mpfr(delta[redo.idx], precBits = precBits))

  if(!is.finite(gradout.redoMPFR))
    warning("MPFR used with insufficient ", precBits,
            " precBits in BB7 gradient for ", parCaller)
  return (gradout.redoMPFR);

}

#############################################################

gradFun4theta_infinite <- function(u, theta, delta,redo.idx)
{
  require("Rmpfr", quietly = TRUE)
  precBits <- 1024
  ## MPFR class used for u, theta,  delta
  gradout.redoMPFR <- gradFun4theta(u = mpfr(u[redo.idx, , drop = FALSE], precBits = precBits),
                                    theta = mpfr(theta[redo.idx], precBits = precBits),
                                    delta = mpfr(delta[redo.idx], precBits = precBits))
  if(!is.finite(gradout.redoMPFR))
    warning("MPFR used with insufficient ", precBits,
            " precBits in BB7 gradient for ", parCaller)
  return (gradout.redoMPFR);
}
