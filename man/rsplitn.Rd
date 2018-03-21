\name{rsplitn}
\alias{rsplitn}

\title{
random generation for Split-Normal distribution
}
\description{
Density, distribution function, quantile function and random generation for the split-normal distribution.
}
\usage{
rsplitn(n, mu, sigma, lmd)
}

\arguments{
  \item{n}{
number of observations. If length(n) > 1, the length is taken to be the number required.
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
The random variable \eqn{y} follows a split-normal distribution, which has density
\deqn{F(y) = 1/(1+\lambda)\sigma \sqrt(2/\pi) exp(-(y-\mu)^2)/2\sigma^2   ,if  y<=\mu} ,
\deqn{F(y) = 1/(1+\lambda)\sigma \sqrt(2/\pi) exp(-(y-\mu)^2)/2\sigma ^2 \lambda ^2  ,if  y>\mu}

where  \eqn{\sigma >0} and \eqn{\lambda >0}. The Split-normal distribution reduces to normal distribution when \eqn{\lambda = 1}
}
\value{
\code{rsplitt} generates random deviates.
Invalid arguments will result in return value NaN, with a warning.

The length of the result is determined by \code{n} for \code{rsplitn}.

The numerical arguments other than n are recycled to the length of the result. Only the first elements of the logical arguments are used.
}
\references{
Li, F., Villani, M., & Kohn, R. (2009). Flexible modeling of conditional distributions using smooth mixtures of asymmetric student t densities. Journal of Statistical Planning & Inference, 140(12), 3638-3654.
}
\author{
Feng Li, Jiayue Zeng
}


\seealso{
For split-normal distribution, also see:
\code{\link{dsplitn}()}
\code{\link{psplitn}()}
\code{\link{qsplitn}()}
\code{\link{split-n}()}
}


\examples{

n = 10
mu <- c(0,1,2)
sigma <- c(0.5,1,2)
lmd <- c(1,2,3)
rsplitn0 = rsplitn(n,mu,sigma,lmd)


}


\keyword{distribution}
\keyword{asymmetric normal}


