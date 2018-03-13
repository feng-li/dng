\name{split-normal}
\alias{split_mean}
\alias{splitn_var}
\alias{splitn_skewness}
\alias{splitn_kurtosis}
\alias{dng_splitn_mean}
\alias{dng_splitn_var}
\alias{dng_splitn_skewness}
\alias{dng_splitn_kurtosis}
\title{Moments of the Split Normal distribution}

\description{
Compute the mean, variance, skewness and kurtosis for the split normal distribution.
}
\usage{
splitn_mean(mu, sigma, lmd)
splitn_var(sigma, lmd)
splitn_skewness(sigma, lmd)
splitn_kurtosis(lmd)
}

\arguments{
  \item{mu}{
vector of location parameter. (The mode of the density)
}
  \item{sigma}{
vector of standard deviations.
}
  \item{lmd}{
vector of skewness parameters (>0). If is 1, reduce to normal distribution.
}
}

\details{
The random variable y follows a split-normal distribution, which has density:
\deqn{1/(1+\lambda)\sigma \sqrt(2/\pi) exp{-(y-\mu)*2/2\sigma^2}, if y<=\mu}
\deqn{1/(1+\lambda)\sigma \sqrt(2/\pi) exp{-(y-\mu)*2/2\sigma^2 \lambda^2}, if y>\mu}
where \eqn{\sigma>0} and \eqn{\lambda>0}. The Split-normal distribution reduce to normal distribution when \eqn{\lmd=1}.

}

\value{
\code{splitn_mean} gives the mean.
\code{splitn_var} gives the variance.
\code{splitn_skewness} gives the skewness.
\code{splitn_kurtosis} gives the kurtosis.
(\code{splitn_mean}, \code{splitn_var},\code{splitn_skeness} and \code{splitn_kurtosis} are all vectors.
}
\references{
Li, F., Villani, M., & Kohn, R. (2009). Flexible modeling of conditional distributions using smooth mixtures of asymmetric student t densities. Journal of Statistical Planning & Inference, 140(12), 3638-3654.
and
Working Paper The Multivariate Split Normal Distribution and Asymmetric Principal Components Analysis Villani, Mattias; Larsson, Rolf Sveriges Riksbank Working Paper Series, No. 175 Central Bank of Sweden, Stockholm.
}
\author{
Feng Li, Jiayue Zeng
}


\seealso{
\code{\link{psplitn}()} \code{\link{dsplitn}()} and for the split-normal distribution.
}
\examples{
mu <- c(0,1,2)
sigma <- c(0.5,1,2)
lmd <- c(1,2,3)

mean0 <- splitn_mean(mu, sigma, lmd)
var0 <- splitn_var(sigma, lmd)
skewness0 <- splitn_skewness(sigma, lmd)
kurtosis0 <- splitn_kurtosis(lmd)
}

\keyword{distribution}
\keyword{asymmetric normal}
