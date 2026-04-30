context("split-t distribution")

test_that("symmetric split-t matches base Student t functions", {
    x <- c(-3, -0.5, 0, 1.25, 4)
    mu <- -0.25
    df <- 9
    phi <- 1.4

    z <- (x - mu) / phi

    expect_equal(
        dsplitt(x, mu, df, phi, lmd = 1, logarithm = FALSE),
        dt(z, df) / phi,
        tolerance = 1e-12
    )
    expect_equal(
        dsplitt(x, mu, df, phi, lmd = 1, logarithm = TRUE),
        dt(z, df, log = TRUE) - log(phi),
        tolerance = 1e-12
    )
    expect_equal(
        psplitt(x, mu, df, phi, lmd = 1),
        pt(z, df),
        tolerance = 1e-12
    )

    p <- c(0.05, 0.2, 0.5, 0.8, 0.95)
    expect_equal(
        qsplitt(p, mu, df, phi, lmd = 1),
        mu + phi * qt(p, df),
        tolerance = 1e-12
    )
})

test_that("split-t CDF and quantile are inverse operations", {
    p <- c(0.01, 0.15, 0.33, 0.7, 0.95, 0.99)
    mu <- c(-1, 0.5)
    df <- c(5, 9, 15)
    phi <- c(0.75, 1.25)
    lmd <- c(0.6, 1.5, 3)

    q <- qsplitt(p, mu, df, phi, lmd)

    expect_equal(
        psplitt(q, mu, df, phi, lmd),
        p,
        tolerance = 1e-12
    )
})

test_that("split-t random generation uses inverse-CDF sampling", {
    n <- 8
    mu <- c(-1, 0.5)
    df <- c(5, 9, 15)
    phi <- c(0.75, 1.25)
    lmd <- c(0.6, 1.5, 3, 1)

    set.seed(202)
    u <- runif(n)
    expected <- qsplitt(
        u,
        rep_len(mu, n),
        rep_len(df, n),
        rep_len(phi, n),
        rep_len(lmd, n)
    )

    set.seed(202)
    actual <- rsplitt(n, mu, df, phi, lmd)

    expect_equal(actual, expected, tolerance = 1e-12)
})
