# test that stm works with quanteda

require(quanteda)

test_that("Test that stm works on a quanteda dfm", {
  require(quanteda)
  if(utils::compareVersion(as.character(utils::packageVersion("quanteda")), "0.9.9-31") >= 0) {
    gadarian_corpus <- corpus(gadarian, text_field = "open.ended.response")
    gadarian_dfm <- dfm(tokens(gadarian_corpus))
    set.seed(10012) # NYU :-)
    stm_from_dfm <- stm(gadarian_dfm,
                        K = 3,
                        prevalence = ~treatment + s(pid_rep),
                        data = docvars(gadarian_corpus),
                        max.em.its=2)
    expect_true(inherits(stm_from_dfm, "STM"))
  } else {
    #basically if the version is old, just skip this test for now.
    expect_identical("STM", "STM")
  }
})

if(requireNamespace("tm",quietly=TRUE) & utils::packageVersion("tm")>="0.6") {
  test_that("Test that stm works on a classic stm object structure", {
    temp <- textProcessor(documents = gadarian$open.ended.response,
                          metadata = gadarian)
    meta <- temp$meta
    vocab <- temp$vocab
    docs <- temp$documents
    out <- prepDocuments(docs, vocab, meta)
    docs <- out$documents
    vocab <- out$vocab
    meta <- out$meta
    set.seed(10012)
    stm_from_stmclassic <- 
      stm(docs, vocab, 3, prevalence = ~treatment + s(pid_rep), data = meta,
          max.em.its = 2)
    expect_true(inherits(stm_from_stmclassic, "STM"))
  })
}
## could add an additional test to compare the two outputs
## could differ based on tokenizer, although
## tm::stopwords("english") is the same set as quanteda::stopwords("english")