test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("demeaning works", {
  set.seed(999)
  x <- rnorm(10000, mean = 3, sd = 2)
  y <- detrendcl(x = x, tt = "constant")
  expect_equal(mean(y), 0)
})

test_that("demeaning works with NAs", {
  set.seed(999)
  x <- rnorm(10000, mean = 3, sd = 2)
  x[1] <- NA
  x[length(x)] <- NA
  y <- detrendcl(x = x, tt = "constant")
  expect_equal(mean(y, na.rm = TRUE), 0)
})

test_that(desc = "taking out a linear trend works", {
  set.seed(99)
  x <- -10:90
  y <- 2 + 3 * x + rnorm(length(x), mean = 0, sd = 3)
  z <- detrendcl(y, tt = "linear")

  m <- lm(y ~ x)

  expect_equal(z, as.vector(m$residuals))
})

test_that(desc = "taking out a linear trend with NAs works", {
  set.seed(99)
  x <- -10:90
  x[1] <- NA
  x[length(x)] <- NA
  y <- 2 + 3 * x + rnorm(length(x), mean = 0, sd = 3)
  z <- detrendcl(y, tt = "linear")

  expect_equal(is.na(z[1]), TRUE)
  expect_equal(is.na(z[length(x)]), TRUE)

  m <- lm(y ~ x)

  z_nona <- z[which(!is.na(z))]
  expect_equal(z_nona, as.vector(m$residuals))
})
