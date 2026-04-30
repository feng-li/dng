context("split-t gradients")

finite_difference <- function(fun, at, step = 1e-6) {
    (fun(at + step) - fun(at - step)) / (2 * step)
}

test_that("gsplitt is exported", {
    expect_true("gsplitt" %in% getNamespaceExports("dng"))
})

test_that("ghypergeo evaluates a simple hypergeometric identity", {
    value <- ghypergeo(
        matrix(2, nrow = 1),
        matrix(numeric(), nrow = 1, ncol = 0),
        z = 0.25,
        k = 200
    )

    expect_equal(drop(value), (1 - 0.25)^(-2), tolerance = 1e-12)
})

test_that("gsplitt gradients match finite differences", {
    y <- c(-1.4, -0.7, 1.2, 3.0)
    params <- list(mu = 0.2, df = 7.5, phi = 1.3, lmd = 2.1)
    par_names <- c("mu", "df", "phi", "lmd")

    for (par_name in par_names) {
        analytic <- gsplitt(
            y,
            list(
                mu = rep(params$mu, length(y)),
                df = rep(params$df, length(y)),
                phi = rep(params$phi, length(y)),
                lmd = rep(params$lmd, length(y))
            ),
            par_name,
            c("u", "d")
        )

        step <- if (par_name == "df") 1e-5 else 1e-6
        numeric_u <- vapply(y, function(yy) {
            finite_difference(function(value) {
                shifted <- params
                shifted[[par_name]] <- value
                psplitt(yy, shifted$mu, shifted$df, shifted$phi, shifted$lmd)
            }, params[[par_name]], step)
        }, numeric(1))

        numeric_d <- vapply(y, function(yy) {
            finite_difference(function(value) {
                shifted <- params
                shifted[[par_name]] <- value
                dsplitt(yy, shifted$mu, shifted$df, shifted$phi, shifted$lmd, logarithm = TRUE)
            }, params[[par_name]], step)
        }, numeric(1))

        expect_lt(max(abs(analytic$u - numeric_u)), 1e-6)
        expect_lt(max(abs(analytic$d - numeric_d)), 1e-6)
    }
})
