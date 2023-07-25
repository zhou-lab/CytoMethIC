context("data")

test_that("test='HM27.address' gives correct data", {
  dt <- sesameData::sesameDataGet("HM27.address")
  expect_is(dt, "list")
})
