\name{psplitn}
\alias{psplitn}
\title{
Distribution function of Split-normal distribution
}
\description{
Density, distribution function, quantile function and random generation for the split normal distribution.
}
\usage{
psplitn(q, mu, sigma, lmd)
}

\arguments{
  \item{q}{
vector of quantiles.
}
  \item{mu}{
vector of location parameter. (The mode of the density)
}
  \item{sigma}{
vector of standard deviations.
}
  \item{lmd}{
vector of skewness parameters (>0). If is 1, reduced to symmetric student t distribution.
}
}
\details{
Computing the CDF of the Asymmetric Normal distribution. The random variable y follows a split-normal distribution, y~N(\eqn{\mu}, \eqn{\sigma}, \eqn{\lmd}).
}
\value{
Computing the CDF of the Asymmetric Normal distribution. The random variable y follows a split-normal distribution, y~N(\eqn{\mu}, \eqn{\sigma}, \eqn{\lmd}).
}
\references{
Li, F., Villani, M., & Kohn, R. (2009). Flexible modeling of conditional distributions using smooth mixtures of asymmetric student t densities. Journal of Statistical Planning & Inference, 140(12), 3638-3654.
}
\author{
Feng Li, Jiayue Zeng
}

\seealso{
\code{\link{dsplitn}()} for the split-normal distribution.
}

\examples{
q <- c(0.25,0.5,0.75)
mu <- c(0,1,2)
sigma <- c(0.5,1,2)
lmd <- c(1,2,3)
dsplitn0 = psplitn(q,mu,sigma,lmd)

}
\keyword{distribution}
\keyword{asymmetric normal}
