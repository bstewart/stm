library(readtext)
library(quanteda)
library(stm)

# DATA_DIR <- system.file("extdata/txt/UDHR", package = "readtext")
# DATA_FILES <- "D:\\git_checkouts\\GrimmerSenatePressReleases\\raw\\*"

# d <- readtext(DATA_FILES, docvarsfrom = "filepaths", dvsep = "[\\\\,/]", docvarnames = c("_", "_", "_", "_", "senator", "_"))
# d <- d[c("doc_id", "senator", "text")]
# c <- corpus(d, docid_field = "doc_id", text_field = "text")
# spr_dfm <- dfm(c, tolower = TRUE, stem = TRUE, remove = stopwords())
spr_dfm <- readRDS("d:\\git_checkouts\\GrimmerSenatePressReleases\\spr_dfm.rds")

K <- 30
ier <- 4
esr <- 5
ebr <- 5

print('----------')
print(K)
print('----------')


seed <- sample(1:1000000, 1)
print(paste0('seed = ', seed))
print(paste0('TRYING = ', ier, ',', esr, ',', ebr))
m1 <- stm(spr_dfm, K=K, seed=seed, cores=1, control=list(init.exp.round=ier, estep.sigma.round=esr, estep.beta.round=ebr))
print(paste0('m1 time ', m1$time))
m2 <- stm(spr_dfm, K=K, seed=seed, cores=7, control=list(init.exp.round=ier, estep.sigma.round=esr, estep.beta.round=ebr))
print(paste0('m2 time ', m2$time))
if (isTRUE(all.equal(m1$mu, m2$mu)) & isTRUE(all.equal(m1$sigma, m2$sigma)) & isTRUE(all.equal(m1$beta, m2$beta)) & isTRUE(all.equal(m1$theta, m2$theta)) & isTRUE(all.equal(m1$eta, m2$eta)) & isTRUE(all.equal(m1$invsigma, m2$invsigma))) {
  print(paste0('SUCCESS = ', ier, ',', esr, ',', ebr))
} else {
  print(all.equal(m1, m2))
}
