test_that("the package loads and reports a version", {
  expect_s3_class(packageVersion("disturbanceMgmt"), "package_version")
})
