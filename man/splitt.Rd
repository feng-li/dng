\name{splitt}
\alias{dsplitt}
\alias{psplitt}
\alias{qsplitt}
\alias{rsplitt}
\title{
Split student-t distribution
}

\description{
Density, distribution function, quantile function and random generation for the split student-t distribution.
}

\usage{
dsplitt(x, mu, df, phi, lmd, log0) %??log0 or log
psplitt(q, mu, df, phi, lmd)
qsplitt(p, mu, df, phi, lmd)
rsplitt(n, mu, df, phi, lmd)
}

\arguments{
  \item{x, q}{
  vector of quantiles.
  }
  \item{p}{
  vector of probabilities.
}
  \item{n}{
  number of observations. If length(n) > 1, the length is taken to be the number required.
  }
  \item{mu}{
  vector of modes.
}
  \item{df}{
  degrees of freedom (> 0, maybe non-integer). df = Inf is allowed.
}
  \item{phi}{
  vector of scale parameters.
}
  \item{lmd}{
  vector of skewness parameters. If is 1, reduced to symmetric student t distribution.
}
  \item{log}{
  logical; if TRUE, probabilities p are given as log(p).
  }
}

\details{
The random variable y follows a split-t distribution with  \eqn{\nu}>0 degrees of freedom, y~t(\eqn{\mu}, \eqn{\phi}, \eqn{\lambda}, \eqn{\nu}), if its density function is of the form
%% FIXME add non-LaTeX alternative \deqn{LaTeX}{non-LaTex}
\deqn{C K(\mu, \phi, \nu,)I(y\leq\mu) + C K(\mu, \lambda\phi, \nu)I(y>\mu),   }
where,
\deqn{K(\mu, \phi, \nu,) =[\nu/(\nu+(y-\mu)^(2) /\phi^{2})]^((\nu+1)/2)   }
is the kernel of a student \eqn{t} density with variance \eqn{\phi^{2}\nu/(\nu-2)}
and
\deqn{c = 2[(1+\lambda)\phi (\sqrt\nu) Beta(\nu/2,1/2)]^(-1)  }is the normalization constant.

}

\value{
\code{dsplitt} {gives the density.}
\code{psplitt} {gives the distribution function.}
\code{qsplitt} {gives the quantile function.}
\code{rsplitt} {generates random deviates.}

Invalid arguments will result in return value NaN, with a warning.

The length of the result is determined by n for rt, and is the maximum of the lengths of the numerical arguments for the other functions.

The numerical arguments other than n are recycled to the length of the result. Only the first elements of the logical arguments are used.

}

\references{
Li, F., Villani, M., & Kohn, R. (2009). Flexible modeling of conditional distributions using smooth mixtures of asymmetric student t densities. Journal of Statistical Planning & Inference, 140(12), 3638-3654.
}

\author{
Feng Li, Jiayue Zeng
}

\seealso{
\code{\link{split-t}()} for numerical characteristics of the Split student-t distribution.
}

\examples{
n <- 5
x <- c(0.25,0.5,0.75)
q <- c(0.25,0.5,0.75)
p <- c(0.25,0.5,0.75)
mu <- c(0,1,2)
df <- rep(10,3)
phi <- c(0.5,1,2)
lmd <- c(1,2,3)

dsplitt0 <- dsplitt(x, mu, df, phi, lmd, log0 = TRUE)
psplitt0 <- psplitt(q, mu, df, phi, lmd)
qsplitt0 <- qsplitt(p, mu, df, phi, lmd)
rsplitt0 <- rsplitt(n, mu, df, phi, lmd)
}

\keyword{distribution}
