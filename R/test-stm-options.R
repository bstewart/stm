context("stm function options")

test_that("Test that CTM stm works", {
  test <- stm(poliblog5k.docs, poliblog5k.voc, K=20, max.em.its = 5)
  expect_identical(class(test), "STM")
})

test_that("Test that spectral init works", {
  test <- stm(poliblog5k.docs, poliblog5k.voc, K=20, max.em.its = 5, init.type = "Spectral")
  expect_identical(class(test), "STM")
})

test_that("Test that Prevalence works", {
  test <- stm(poliblog5k.docs, poliblog5k.voc, 
              prevalence=~s(day)*blog,K=20, max.em.its = 5, init.type = "Spectral",
              data=poliblog5k.meta)
  expect_identical(class(test), "STM")
})

test_that("Test that Content works", {
  test <- stm(poliblog5k.docs, poliblog5k.voc, 
              content=~blog,K=20, max.em.its = 5, init.type = "Spectral",
              data=poliblog5k.meta)
  expect_identical(class(test), "STM")
})

