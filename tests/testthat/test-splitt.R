context("splitt")

test_that("Splitt Densities",
{
    expect_equal(dsplit(-3, mu = 0, df = 10, phi = 1, lmd = 1, logarithm = FALSE)
                 1-dsplit(3, mu = 0, df = 10, phi = 1, lmd = 1, logarithm = FALSE))


})
