logDensFun_infinite <- function(u, theta, delta, redo.idx)
{
  require("Rmpfr")
  precBits <- 1024
  ## MPFR class used for u, theta,  delta
  out.logredoMPFR <- logDensFun(u = mpfr(u[redo.idx, , drop = FALSE], precBits = precBits),
                                theta = mpfr(theta[redo.idx], precBits = precBits),
                                delta = mpfr(delta[redo.idx], precBits = precBits))
  out.logredo <- as.numeric(out.logredoMPFR)

  if(!is.finite(out.logredo))
    warning("MPFR used with insufficient ", precBits, " precBits in BB7 density.")

  return(out.logredo)
}
