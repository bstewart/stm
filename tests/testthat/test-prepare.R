context("Prepare text")

test_that("textProcessor and prepDocuments correctly subset metadata ", {
  data("gadarian")

  txtOut <- textProcessor(documents = gadarian$open.ended.response,
                          metadata = data.frame(MetaID = gadarian$MetaID),
                          sparselevel = .8)
  expect_equal(nrow(txtOut$meta), length(txtOut$documents))

  prepped <- prepDocuments(txtOut$documents, txtOut$vocab,
                           txtOut$meta, upper.thresh = 100)
  expect_equal(nrow(prepped$meta), length(prepped$documents))
})
