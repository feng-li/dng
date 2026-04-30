#include <Rcpp.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

using namespace Rcpp;

namespace {

const double kTol = 1e-12;
const int kMaxTerms = 100000;

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

double hypergeo_pfq(const std::vector<double>& a,
                    const std::vector<double>& b,
                    const double z,
                    const int max_terms = kMaxTerms,
                    const double tol = kTol)
{
  if (!std::isfinite(z)) {
    return R_NaN;
  }

  double term = 1.0;
  double sum = 1.0;

  for (int n = 1; n <= max_terms; ++n) {
    double ratio = z / static_cast<double>(n);

    for (std::size_t i = 0; i < a.size(); ++i) {
      ratio *= a[i] + n - 1.0;
    }

    for (std::size_t i = 0; i < b.size(); ++i) {
      const double denominator = b[i] + n - 1.0;
      if (denominator == 0.0) {
        return R_NaN;
      }
      ratio /= denominator;
    }

    term *= ratio;
    sum += term;

    if (!std::isfinite(sum) || !std::isfinite(term)) {
      return R_NaN;
    }

    if (std::fabs(term) <= tol * std::max(1.0, std::fabs(sum))) {
      break;
    }
  }

  return sum;
}

double beta_density(const double x, const double a, const double b)
{
  if (x <= 0.0 || x >= 1.0) {
    return 0.0;
  }

  return std::exp((a - 1.0) * std::log(x) +
                  (b - 1.0) * std::log1p(-x) -
                  ::Rf_lbeta(a, b));
}

double regularized_beta_da(const double x, const double a, const double b)
{
  if (x <= 0.0 || x >= 1.0) {
    return 0.0;
  }

  const double ix = ::Rf_pbeta(x, a, b, TRUE, FALSE);
  const double h = hypergeo_pfq({a, a, 1.0 - b}, {a + 1.0, a + 1.0}, x);
  const double tail = std::exp(a * std::log(x) -
                               2.0 * std::log(a) -
                               ::Rf_lbeta(a, b)) * h;

  return ix * (std::log(x) - R::digamma(a) + R::digamma(a + b)) - tail;
}

double regularized_beta_dnu(const double x, const double nu)
{
  if (x <= 0.0 || x >= 1.0) {
    return 0.0;
  }

  const double a = 0.5 * nu;
  const double b = 0.5;
  const double d_da = regularized_beta_da(x, a, b);
  const double d_dx = beta_density(x, a, b);
  const double d_x_d_nu = x * (1.0 - x) / nu;

  return 0.5 * d_da + d_dx * d_x_d_nu;
}

bool invalid_parameters(const double nu, const double phi, const double lambda)
{
  return !std::isfinite(nu) || !std::isfinite(phi) || !std::isfinite(lambda) ||
    nu <= 0.0 || phi <= 0.0 || lambda <= 0.0;
}

} // namespace

//' Generalized hypergeometric function
//'
//' Evaluate generalized hypergeometric series used by the split-t gradient
//' calculations.
//'
//' @param a matrix of upper hypergeometric parameters.
//' @param b matrix of lower hypergeometric parameters.
//' @param z vector of hypergeometric function arguments.
//' @param k maximum number of hypergeometric series terms. Non-positive values
//' use the package default.
//' @return A one-column numeric matrix of generalized hypergeometric function
//' values. Rows of \code{a}, rows of \code{b}, and values of \code{z} are
//' recycled to the output length.
//' @seealso \code{\link{gsplitt}()}
//' @export
// [[Rcpp::export]]
NumericMatrix ghypergeo(NumericMatrix a, NumericMatrix b, NumericVector z, int k)
{
  const int n = std::max(static_cast<int>(z.size()),
                         std::max(a.nrow(), b.nrow()));
  const int max_terms = k > 0 ? k : kMaxTerms;
  NumericMatrix out(n, 1);

  for (int i = 0; i < n; ++i) {
    std::vector<double> upper(a.ncol());
    std::vector<double> lower(b.ncol());
    const int a_row = a.nrow() == 1 ? 0 : i % a.nrow();
    const int b_row = b.nrow() == 1 ? 0 : i % b.nrow();

    for (int j = 0; j < a.ncol(); ++j) {
      upper[j] = a(a_row, j);
    }

    for (int j = 0; j < b.ncol(); ++j) {
      lower[j] = b(b_row, j);
    }

    out(i, 0) = hypergeo_pfq(upper, lower, z[i % z.size()], max_terms);
  }

  return out;
}

//' @describeIn splitt Gradients for the split-t CDF and log-density.
//' @param y vector of quantiles for gradient evaluation.
//' @param par list with \code{mu}, \code{df}, \code{phi}, and \code{lmd} parameter vectors.
//' @param parCaller character scalar naming the parameter to differentiate:
//' \code{"mu"}, \code{"df"}, \code{"phi"}, or \code{"lmd"}.
//' @param denscaller character vector selecting gradients to compute. Use
//' \code{"u"} for the CDF gradient and \code{"d"} for the log-density gradient.
//' @export
// [[Rcpp::export]]
List gsplitt(NumericVector y, List par, std::string parCaller, GenericVector denscaller)
{
  const int n = y.size();
  NumericVector mu = rep_len(as<NumericVector>(par["mu"]), n);
  NumericVector df = rep_len(as<NumericVector>(par["df"]), n);
  NumericVector phi = rep_len(as<NumericVector>(par["phi"]), n);
  NumericVector lmd = rep_len(as<NumericVector>(par["lmd"]), n);

  NumericVector outu(n, R_NaN);
  NumericVector outd(n, R_NaN);

  const std::string parameter = lower_string(parCaller);
  const bool compute_u = denscaller.size() == 0 || starts_with_flag(denscaller, 'u');
  const bool compute_d = denscaller.size() == 0 || starts_with_flag(denscaller, 'd');

  if (parameter != "mu" && parameter != "df" &&
      parameter != "phi" && parameter != "lmd") {
    stop("No such parameter!");
  }

  for (int i = 0; i < n; ++i) {
    const double nu = df[i];
    const double scale = phi[i];
    const double lambda = lmd[i];

    if (invalid_parameters(nu, scale, lambda)) {
      continue;
    }

    const double d = y[i] - mu[i];
    const bool left = y[i] <= mu[i];
    const double side = left ? 1.0 : lambda;
    const double side2 = side * side;
    const double one_plus_lambda = 1.0 + lambda;
    const double beta0 = R::beta(0.5 * nu, 0.5);
    const double denominator = d * d + side2 * nu * scale * scale;
    const double a = side2 * nu * scale * scale / denominator;
    const double a_power = std::pow(a, 0.5 * nu);
    const double sqrt_denominator = std::sqrt(denominator);

    if (compute_u) {
      if (parameter == "mu") {
        outu[i] = -2.0 * side * a_power /
          (one_plus_lambda * beta0 * sqrt_denominator);
      } else if (parameter == "df") {
        const double multiplier = left ?
          1.0 / one_plus_lambda :
          -lambda / one_plus_lambda;
        outu[i] = multiplier * regularized_beta_dnu(a, nu);
      } else if (parameter == "phi") {
        outu[i] = -2.0 * side * d * a_power /
          (one_plus_lambda * scale * beta0 * sqrt_denominator);
      } else if (parameter == "lmd") {
        const double ibeta = ::Rf_pbeta(a, 0.5 * nu, 0.5, TRUE, FALSE);

        if (left) {
          outu[i] = -ibeta / (one_plus_lambda * one_plus_lambda);
        } else {
          outu[i] = -ibeta / (one_plus_lambda * one_plus_lambda) -
            2.0 * d * a_power /
              (one_plus_lambda * beta0 * sqrt_denominator);
        }
      }
    }

    if (compute_d) {
      if (parameter == "mu") {
        outd[i] = (1.0 + nu) * d / denominator;
      } else if (parameter == "df") {
        const double c0 = d * d - scale * scale * side2;
        const double c3 = std::log(nu / (nu + d * d / (scale * scale * side2)));
        const double digamma_diff = R::digamma(0.5 * nu) -
          R::digamma(0.5 * (nu + 1.0));

        outd[i] = 0.5 * (c0 / denominator + c3 - digamma_diff);
      } else if (parameter == "phi") {
        const double c0 = d * d - scale * scale * side2;
        outd[i] = nu * c0 / (scale * denominator);
      } else if (parameter == "lmd") {
        if (left) {
          outd[i] = -1.0 / one_plus_lambda;
        } else {
          outd[i] = -1.0 / one_plus_lambda +
            (1.0 + nu) * d * d / (lambda * denominator);
        }
      }
    }
  }

  return List::create(_["u"] = outu,
                      _["d"] = outd);
}
