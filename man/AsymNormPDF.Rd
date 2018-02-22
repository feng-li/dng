\name{AsymNormPDF}
\alias{AsymNormPDF}

\title{
Split-Normal Distribution
}
\description{
Computing the PDF or CDF (cum=1) of the Asymmetric Normal distribution. Mean and standard deviation of the distribution is also computed.
}
\usage{
AsymNormPDF(y, mu, sigma, lmd, cum, sumStat)
}

\arguments{
  \item{y}{
n-by-1 Point where the density is evaluated.
}
  \item{mu}{
vector of location parameter. (The mode of the density)
}
  \item{sigma}{
vector of standard deviations.
}
  \item{lmd}{
vector of skewness parameters (> 0).
}
  \item{cum}{
bool If cum=TRUE, then CDF is returned, otherwise the PDF is returned.
}
  \item{sumStat}{
bool If sumStat=TRUE, then mean and standard deviation is computed, otherwise won't compute it.
}
}
\details{
Computing the PDF or CDF (cum=1) of the Asymmetric Normal distribution. Mean and standard deviation of the distribution is also computed.
The random variable y follows a split-normal distribution, y~N(\eqn{\mu}, \eqn{\sigma}, \eqn{\lmd}).
}
\value{
\code{AsymNormPDF} gives the density: vector with PDF/CDF evaluations,
densMean: vector with means,
and densStd: vector with standard deviations.
Invalid arguments will result in return value NaN, with a warning.
}
\references{
Li, F., Villani, M., & Kohn, R. (2009). Flexible modeling of conditional distributions using smooth mixtures of asymmetric student t densities. Journal of Statistical Planning & Inference, 140(12), 3638-3654.
}
\author{
Feng Li, Jiayue Zeng
}

\seealso{
\code{\link{dsplitt}()}, \code{\link{psplitt}()}, \code{\link{qsplitt}()} and \code{\link{rsplitt}()} for the split-t distribution.
}
\examples{
y <- c(0.25,0.5,0.75)
mu <- c(0,1,2)
sigma <- rep(10,3)
lmd <- c(1,2,3)

asy <- AsymNormPDF(y, mu, sigma, lmd, cum = TRUE, sumStat= TRUE)

}

\keyword{distribution}
\keyword{asymmetric normal}
