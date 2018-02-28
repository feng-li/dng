\name{dsplitn}
\alias{dsplitn}

\title{
split-normal distribution
}
\description{
Density function for the split normal distribution.
}
\usage{
dsplitn(x, mu, phi, sigma, logarithm)
}

\arguments{
  \item{x}{
vector of quantiles.
}
  \item{mu}{
vector of location parameter. (The mode of the density)
}
  \item{phi}{
vector of scale parameters (>0).
}
  \item{sigma}{
vector of standard deviations.
}
  \item{logarithm}{
logical; if TRUE, probabilities p are given as log(p).
}
}
\details{
Computing the PDF of the Asymmetric Normal distribution. The random variable y follows a split-normal distribution, y~N(\eqn{\mu}, \eqn{\sigma}, \eqn{\lmd}).
}
\value{
Computing the PDF of the Asymmetric Normal distribution. The random variable y follows a split-normal distribution, y~N(\eqn{\mu}, \eqn{\sigma}, \eqn{\lmd}).
}
\references{
Li, F., Villani, M., & Kohn, R. (2009). Flexible modeling of conditional distributions using smooth mixtures of asymmetric student t densities. Journal of Statistical Planning & Inference, 140(12), 3638-3654.
}
\author{
Feng Li, Jiayue Zeng
}

\seealso{
\code{\link{psplitn}()} for the split-normal distribution.
}

\examples{
x <- c(0.25,0.5,0.75)
mu <- c(0,1,2)
phi <- c(0.5,1,2)
sigma <- c(0.5,1,2)
dsplitn0 = dsplitn(x, mu, phi, sigma, TRUE)
}

\keyword{distribution}
\keyword{asymmetric normal}
