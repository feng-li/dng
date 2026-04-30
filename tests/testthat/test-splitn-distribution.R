context("split-normal distribution")

test_that("symmetric split-normal matches base normal functions", {
    x <- c(-2.5, -0.25, 0, 0.75, 2.25)
    mu <- 0.4
    sigma <- 1.7

    expect_equal(
        dsplitn(x, mu, sigma, lmd = 1, logarithm = FALSE),
        dnorm(x, mu, sigma),
        tolerance = 1e-12
    )
    expect_equal(
        dsplitn(x, mu, sigma, lmd = 1, logarithm = TRUE),
        dnorm(x, mu, sigma, log = TRUE),
        tolerance = 1e-12
    )
    expect_equal(
        psplitn(x, mu, sigma, lmd = 1),
        pnorm(x, mu, sigma),
        tolerance = 1e-12
    )

    p <- c(0.05, 0.2, 0.5, 0.8, 0.95)
    expect_equal(
        qsplitn(p, mu, sigma, lmd = 1),
        qnorm(p, mu, sigma),
        tolerance = 1e-12
    )
})

test_that("split-normal CDF and quantile are inverse operations", {
    p <- c(0.01, 0.15, 0.33, 0.7, 0.95, 0.99)
    mu <- c(-1, 0.5)
    sigma <- c(0.75, 1.25, 2)
    lmd <- c(0.6, 1.5, 3)

    q <- qsplitn(p, mu, sigma, lmd)

    expect_equal(
        psplitn(q, mu, sigma, lmd),
        p,
        tolerance = 1e-12
    )
})

test_that("split-normal random generation uses inverse-CDF sampling", {
    n <- 8
    mu <- c(-1, 0.5)
    sigma <- c(0.75, 1.25, 2)
    lmd <- c(0.6, 1.5, 3, 1)

    set.seed(101)
    u <- runif(n)
    expected <- qsplitn(u, rep_len(mu, n), rep_len(sigma, n), rep_len(lmd, n))

    set.seed(101)
    actual <- rsplitn(n, mu, sigma, lmd)

    expect_equal(actual, expected, tolerance = 1e-12)
})
