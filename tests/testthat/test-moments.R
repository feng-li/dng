context("distribution moments")

test_that("symmetric split-normal moments reduce to normal moments", {
    mu <- c(-2, 0, 3)
    sigma <- c(0.5, 1.25, 2)
    lmd <- rep(1, length(mu))

    expect_equal(splitn_mean(mu, sigma, lmd), mu, tolerance = 1e-12)
    expect_equal(splitn_var(sigma, lmd), sigma^2, tolerance = 1e-12)
    expect_equal(splitn_skewness(sigma, lmd), rep(0, length(mu)), tolerance = 1e-12)
    expect_equal(splitn_kurtosis(lmd), rep(0, length(mu)), tolerance = 1e-12)
})

test_that("symmetric split-t moments reduce to Student t moments", {
    mu <- c(-2, 0, 3)
    df <- c(6, 8, 10)
    phi <- c(0.5, 1.25, 2)
    lmd <- rep(1, length(mu))

    expect_equal(splitt_mean(mu, df, phi, lmd), mu, tolerance = 1e-12)
    expect_equal(
        splitt_var(df, phi, lmd),
        df / (df - 2) * phi^2,
        tolerance = 1e-12
    )
    expect_equal(splitt_skewness(df, phi, lmd), rep(0, length(mu)), tolerance = 1e-12)
    expect_equal(splitt_kurtosis(df, phi, lmd), 6 / (df - 4), tolerance = 1e-12)
})
