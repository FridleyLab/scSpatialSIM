test_that("SpatSimObj Created", {
  dat = CreateSimulationObject()
  expect_equal(class(dat)[1], "SpatSimObj")
})
