library(stm)
context("stm function options")

test_that("Test that CTM stm works", {
  test <- stm(poliblog5k.docs, poliblog5k.voc, K=3, max.em.its = 1, init.type = "Random")
  expect_true(inherits(test,"STM"))
})

test_that("Test that spectral init works", {
  test <- stm(poliblog5k.docs, poliblog5k.voc, K=3, max.em.its = 1, init.type = "Spectral")
  expect_true(inherits(test, "STM"))
})

test_that("Test that Prevalence works", {
  test <- stm(poliblog5k.docs, poliblog5k.voc, 
              prevalence=~s(day)*blog,K=3, max.em.its = 2, init.type = "Spectral",
              data=poliblog5k.meta)
  expect_true(inherits(test, "STM"))
})

test_that("Test that Content works", {
  test <- stm(poliblog5k.docs, poliblog5k.voc, 
              content=~blog,K=3, max.em.its = 2, init.type = "Spectral",
              data=poliblog5k.meta)
  expect_true(inherits(test,"STM"))
})

test_that("Test that ngroups works", {
  test <- stm(poliblog5k.docs, poliblog5k.voc, K=3, init.type="Random",
              ngroups=2, max.em.its=2)
  expect_true(inherits(test,"STM"))
})

#test_that("Test that multi-core searchK works", {
#  test <- searchK(poliblog5k.docs, poliblog5k.voc, K=c(3,4), 
#                     init.type="Random", max.em.its=1, cores = 2)
#  expect_identical(class(test), "searchK")
#})
