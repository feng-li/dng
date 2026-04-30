#include <Rcpp.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <string>

using namespace Rcpp;

namespace {

std::string lower_string(std::string value)
{
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return value;
}

bool starts_with_flag(const GenericVector& values, const char flag)
{
  for (R_xlen_t i = 0; i < values.size(); ++i) {
    if (Rf_isNull(values[i])) {
      continue;
    }

    std::string value = lower_string(as<std::string>(values[i]));
    if (!value.empty() && value[0] == flag) {
      return true;
    }
  }

  return false;
}

bool invalid_parameters(const double sigma, const double lambda)
{
  return !std::isfinite(sigma) || !std::isfinite(lambda) ||
    sigma <= 0.0 || lambda <= 0.0;
}

} // namespace

//' @describeIn splitn Gradients for the split-normal CDF and log-density.
//' @param y vector of quantiles for gradient evaluation.
//' @param par list with \code{mu}, \code{sigma}, and \code{lmd} parameter vectors.
//' @param parCaller character scalar naming the parameter to differentiate:
//' \code{"mu"}, \code{"sigma"}, or \code{"lmd"}.
//' @param denscaller character vector selecting gradients to compute. Use
//' \code{"u"} for the CDF gradient and \code{"d"} for the log-density gradient.
//' @export
// [[Rcpp::export]]
List gsplitn(NumericVector y, List par, std::string parCaller, GenericVector denscaller)
{
  const int n = y.size();
  NumericVector mu = rep_len(as<NumericVector>(par["mu"]), n);
  NumericVector sigma = rep_len(as<NumericVector>(par["sigma"]), n);
  NumericVector lmd = rep_len(as<NumericVector>(par["lmd"]), n);

  NumericVector outu(n, R_NaN);
  NumericVector outd(n, R_NaN);

  const std::string parameter = lower_string(parCaller);
  const bool compute_u = denscaller.size() == 0 || starts_with_flag(denscaller, 'u');
  const bool compute_d = denscaller.size() == 0 || starts_with_flag(denscaller, 'd');

  if (parameter != "mu" && parameter != "sigma" && parameter != "lmd") {
    stop("No such parameter!");
  }

  for (int i = 0; i < n; ++i) {
    const double scale = sigma[i];
    const double lambda = lmd[i];

    if (invalid_parameters(scale, lambda)) {
      continue;
    }

    const double d = y[i] - mu[i];
    const bool left = y[i] <= mu[i];
    const double side = left ? 1.0 : lambda;
    const double side2 = side * side;
    const double z = d / (scale * side);
    const double phi_z = R::dnorm4(z, 0.0, 1.0, FALSE);
    const double one_plus_lambda = 1.0 + lambda;

    if (compute_u) {
      if (parameter == "mu") {
        outu[i] = -2.0 * phi_z / (one_plus_lambda * scale);
      } else if (parameter == "sigma") {
        outu[i] = -2.0 * d * phi_z /
          (one_plus_lambda * scale * scale);
      } else if (parameter == "lmd") {
        const double cdf_z = R::pnorm5(z, 0.0, 1.0, TRUE, FALSE);

        if (left) {
          outu[i] = -2.0 * cdf_z /
            (one_plus_lambda * one_plus_lambda);
        } else {
          outu[i] = 2.0 *
            (cdf_z - 1.0 - one_plus_lambda * z * phi_z) /
            (one_plus_lambda * one_plus_lambda);
        }
      }
    }

    if (compute_d) {
      if (parameter == "mu") {
        outd[i] = d / (scale * scale * side2);
      } else if (parameter == "sigma") {
        outd[i] = -1.0 / scale +
          d * d / (scale * scale * scale * side2);
      } else if (parameter == "lmd") {
        if (left) {
          outd[i] = -1.0 / one_plus_lambda;
        } else {
          outd[i] = -1.0 / one_plus_lambda +
            d * d / (scale * scale * lambda * lambda * lambda);
        }
      }
    }
  }

  return List::create(_["u"] = outu,
                      _["d"] = outd);
}
