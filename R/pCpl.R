pCpl <- function(u, CplNM, parCpl, log = FALSE)
{
  if(tolower(CplNM) == "bb7")
  {
    p <- 2
    theta <- parCpl[["theta"]]
    delta <- parCpl["delta"]
    u1 <- u[, 1]
    u2 <- u[, 2]
    percentile <- (1 - (1 - ((1 - (1 - u1)^theta)^(-delta) +
                               (1 - (1 -  u2)^theta)^(-delta) - 1)^(-1/delta))^(1/theta))
    out <- matrix(percentile)
  }

  else if(tolower(CplNM) == "gumbel")
{# Joe 1997. p.142 Family B6
  delta <- as.vector(parCpl[["delta"]])
  u.tilde <- -log(u)
  percentile.log <- -rowSums(u.tilde^delta)^(1/delta)
  out.log <- matrix(percentile.log)

  if(log)
  {
    out <- out.log
  }
  else
  {
    out <- exp(out.log)
  }
  }

  else
  {
    stop("Given copula is not implemented.")
  }

  return(out)
}
