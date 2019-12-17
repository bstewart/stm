context("Visualize results")

test_that("plot.STM doesn't throw error ", {
  load(url("http://goo.gl/VPdxlS"))
  tryCatch(plot.STM(poliblogContent, type = "perspectives", topics = 11),
           error = function(q) q) -> tca
  expect_true(!inherits(class(tca),"error"))
})