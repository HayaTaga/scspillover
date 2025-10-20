test_that("scspill object contract is satisfied", {
  fit <- sc_spillover(dat, ...) # スモールケース
  expect_s3_class(fit, "scspill")

  # 必須キー
  expect_true(is.vector(fit$inputs$times_post))
  expect_equal(length(fit$inputs$times_post), nrow(fit$inputs$Yc_post))

  expect_true(is.character(fit$inputs$units$control))
  expect_equal(length(fit$inputs$units$control), ncol(fit$inputs$Yc_post))

  # 型
  expect_type(fit$inputs$Y0_pre, "double")
  expect_s3_class(fit$inputs$Yc_pre, "matrix")
})
