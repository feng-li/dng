\name{qsplitn}
\alias{qsplitn}

\title{
quantile function of split-normal distribution
}
\description{
Density, distribution function, quantile function and random generation for the split-normal distribution.
}
\usage{
qsplitn(p,mu,sigma,lmd)
}

\arguments{
  \item{p}{
vector of probability.
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
Gives the quantile functio Asymmetric Normal distribution. The random variable y follows a split-normal distribution, y~N(\eqn{\mu}, \eqn{\sigma}, \eqn{\lmd}).
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
\code{\link{rsplitn}()}
\code{\link{split-n}()}
}

\examples{

p <- c(0.25,0.5,0.75)
mu <- c(0,1,2)
sigma <- c(0.5,1,2)
lmd <- c(1,2,3)
qsplitn0 = qsplitn(p,mu,sigma,lmd)


}


\keyword{distribution}
\keyword{asymmetric normal}
