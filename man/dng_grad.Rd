\name{dng_grad}
\alias{dng_grad}

\title{
Log likelihood for marginal densities
}
\description{
This function is used to compute the log gradient with respect to the parameters in the marginal part of copula model.
}
\usage{
dng_grad(y, par, type, parCaller, denscaller)
}

\arguments{
  \item{y}{
random variable vector.
}
  \item{par}{
a list contains location parameter \code{mu}, Degrees of freedom \code{df},Scaling Parameter \code{phi},Skewness Parameter \code{lmd}.
}
  \item{type}{
a character string giving the type of the marginal model. This must partially match one of "\code{gaussian}","\code{splitt}","\code{poisson}", with default"\code{gaussian}".
}
  \item{parCaller}{
a character string indicates which parameter is called. This must partially match one of "\code{mu}","\code{df}","\code{phi}","\code{lmd}".
}
  \item{denscaller}{
%%     ~~Describe \code{denscaller} here~~
}
}

\details{
This function is used to compute the log gradient with respect to the parameters in the marginal part of copula model.

The random variable y follows a split-t distribution with  \eqn{\nu}>0 degrees of freedom, y~t(\eqn{\mu}, \eqn{\phi}, \eqn{\lambda}, \eqn{\nu}), if its density function is of the form

\deqn{C K(\mu, \phi, \nu,)I(y\leq\mu) + C K(\mu, \lambda \phi, \nu)I(y>\mu),   }
where,
\deqn{K(\mu, \phi, \nu,) =[\nu/(\nu+(y-\mu)^2 /\phi ^2)]^{(\nu+1)/2}   }
is the kernel of a student \eqn{t} density with variance \eqn{\phi ^2\nu/(\nu-2)}
and
\deqn{c = 2[(1+\lambda)\phi (\sqrt \nu) Beta(\nu/2,1/2)]^{-1}  }is the normalization constant.

If y~t(\eqn{\mu}, \eqn{\phi}, \eqn{\lambda}, \eqn{\nu}) then,
\deqn{E(y) = \mu + h},
\deqn{V(y) = (1+\lambda ^3)/(1 + \lambda) \nu/(\nu-2)\phi^2 - h^2 },
\deqn{E[y - E(y)]^3 = 2h^3 + 2h\phi ^2(\lambda ^2+1)\nu/(\nu-3) - 3h\phi ^2(\lambda ^3 + 1)/(\lambda+1) \nu/(\nu-2) },
\deqn{E[y - E(y)]^4 =(3\nu ^2\phi ^4(1+\lambda ^5))/((1+\lambda)(\nu-2)(\nu-4)) - 4h^4 + (6h^(2)(1+\lambda ^3)\nu\phi ^2)/((1+\lambda)(\nu-2)) - (8h^2(\lambda ^2\nu\phi ^2))/(\nu-3). }

}
\value{
This gives log likelihood for marginal densities.
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y, par, type, parCaller, denscaller)
{
    .Call(dng_dng_grad, y, par, type, parCaller, denscaller)
  }
}

\keyword{gradient}
\keyword{distribution}
