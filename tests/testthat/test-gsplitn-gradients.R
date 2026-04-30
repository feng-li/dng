context("split-normal gradients")

finite_difference_splitn <- function(fun, at, step = 1e-6) {
    (fun(at + step) - fun(at - step)) / (2 * step)
}

test_that("gsplitn is exported", {
    expect_true("gsplitn" %in% getNamespaceExports("dng"))
})

test_that("gsplitn gradients match finite differences", {
    y <- c(-1.4, -0.7, 1.2, 3.0)
    params <- list(mu = 0.2, sigma = 1.3, lmd = 2.1)
    par_names <- c("mu", "sigma", "lmd")

    for (par_name in par_names) {
        analytic <- gsplitn(
            y,
            list(
                mu = rep(params$mu, length(y)),
                sigma = rep(params$sigma, length(y)),
                lmd = rep(params$lmd, length(y))
            ),
            par_name,
            c("u", "d")
        )

        numeric_u <- vapply(y, function(yy) {
            finite_difference_splitn(function(value) {
                shifted <- params
                shifted[[par_name]] <- value
                psplitn(yy, shifted$mu, shifted$sigma, shifted$lmd)
            }, params[[par_name]])
        }, numeric(1))

        numeric_d <- vapply(y, function(yy) {
            finite_difference_splitn(function(value) {
                shifted <- params
                shifted[[par_name]] <- value
                dsplitn(yy, shifted$mu, shifted$sigma, shifted$lmd, logarithm = TRUE)
            }, params[[par_name]])
        }, numeric(1))

        expect_lt(max(abs(analytic$u - numeric_u)), 1e-6)
        expect_lt(max(abs(analytic$d - numeric_d)), 1e-6)
    }
})
