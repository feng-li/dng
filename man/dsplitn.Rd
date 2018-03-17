\name{dsplitn}

\alias{dsplitn}

\title{
  split-normal distribution
}
\description{
  Density function for the split normal distribution.
}
\usage{
dsplitn(x, mu, sigma, lmd,logarithm)
}

\arguments{
  \item{x}{
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
  \item{logarithm}{
    logical; if TRUE, probabilities p are given as log(p).
  }
}
\details{
  Computing the PDF of the Asymmetric Normal distribution. The random variable y follows a split-normal distribution, y~N(\eqn{\mu}, \eqn{\sigma}, \eqn{\lmd}), which has density:
  \deqn{1/(1+\lambda)\sigma \sqrt(2/\pi) exp{-(y-\mu)*2/2\sigma^2}, if y<=\mu}
  \deqn{1/(1+\lambda)\sigma \sqrt(2/\pi) exp{-(y-\mu)*2/2\sigma^2 \lambda^2}, if y>\mu}
  where \eqn{\sigma>0} and \eqn{\lambda>0}. The Split-normal distribution reduce to normal distribution when \eqn{\lmd=1}.
}
\value{
  Computing the PDF of the Asymmetric Normal distribution. The random variable y follows a split-normal distribution, y~N(\eqn{\mu}, \eqn{\sigma}, \eqn{\lmd}).
}
\references{
  Li, F., Villani, M., & Kohn, R. (2010). Flexible modeling of conditional distributions using smooth mixtures of asymmetric student t densities. Journal of Statistical Planning & Inference, 140(12), 3638-3654.
  Villani, M., & Larsson, R. (2006) The Multivariate Split Normal Distribution and Asymmetric Principal Components Analysis. Sveriges Riksbank Working Paper Series, No. 175.
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
lmd <- c(0.5,1,2)
sigma <- c(0.5,1,2)
dsplitn0 = dsplitn(x, mu, sigma,lmd, TRUE)
}

\keyword{distribution}
\keyword{asymmetric normal}
