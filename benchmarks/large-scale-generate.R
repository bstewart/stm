#!/usr/bin/env Rscript

source(file.path("benchmarks", "large-scale-common.R"))

options <- parse_benchmark_args(list(
  output="large-scale-anchored-50k.rds",
  documents=50000L,
  topics=20L,
  vocabulary=5000L,
  mean_length=80L,
  anchors_per_topic=10L,
  anchor_mass=0.45,
  initialization_noise=0.55,
  seed=20260805L
))

N <- options$documents
K <- options$topics
V <- options$vocabulary
A <- options$anchors_per_topic
if(V <= K * A) stop("vocabulary must exceed topics * anchors_per_topic")

set.seed(options$seed)
background.size <- V - K * A
background.probability <- seq_len(background.size)^(-0.8)
background.probability <- background.probability/
  sum(background.probability)
anchors <- lapply(seq_len(K), function(topic) {
  (topic - 1L) * A + seq_len(A)
})

beta <- matrix(0, nrow=K, ncol=V)
beta[, K * A + seq_len(background.size)] <-
  (1 - options$anchor_mass) *
  matrix(background.probability, nrow=K, ncol=background.size, byrow=TRUE)
for(topic in seq_len(K)) {
  beta[topic, anchors[[topic]]] <- options$anchor_mass/A
}

eta <- matrix(rnorm(N * (K - 1L), sd=0.85), nrow=N, ncol=K - 1L)
theta <- softmax_rows(eta)
lengths <- pmax(20L, rpois(N, lambda=options$mean_length - 20L) + 20L)
total.tokens <- sum(lengths)

# Draw common token-level randomness in contiguous pools. Only the
# document-specific topic draw remains inside the loop.
anchor.draw <- runif(total.tokens) < options$anchor_mass
anchor.slot <- sample.int(A, total.tokens, replace=TRUE)
background.draw <- sample.int(
  background.size, total.tokens, replace=TRUE,
  prob=background.probability
)

documents <- vector("list", N)
offset <- 0L
progress.every <- max(1L, floor(N/20L))
for(document in seq_len(N)) {
  positions <- offset + seq_len(lengths[document])
  topics <- sample.int(K, lengths[document], replace=TRUE,
                       prob=theta[document,])
  words <- K * A + background.draw[positions]
  use.anchor <- anchor.draw[positions]
  words[use.anchor] <- (topics[use.anchor] - 1L) * A +
    anchor.slot[positions[use.anchor]]
  compressed <- rle(sort.int(words, method="quick"))
  documents[[document]] <- rbind(
    as.integer(compressed$values),
    as.integer(compressed$lengths)
  )
  offset <- offset + lengths[document]
  if(document %% progress.every == 0L) {
    cat(sprintf("Generated %d / %d documents\n", document, N))
  }
}

# STM's native corpus contract requires that every vocabulary index appear at
# least once. Tiny smoke cases may miss rare tail words, so trim and renumber
# them; at the intended 50K scale this is normally a no-op.
observed <- sort(unique(unlist(lapply(documents, function(document) {
  document[1,]
}), use.names=FALSE)))
original.vocab <- sprintf("term_%05d", seq_len(V))
if(length(observed) < V) {
  lookup <- integer(V)
  lookup[observed] <- seq_along(observed)
  documents <- lapply(documents, function(document) {
    document[1,] <- lookup[document[1,]]
    document
  })
  beta <- beta[,observed,drop=FALSE]
  beta <- beta/rowSums(beta)
  anchors <- lapply(anchors, function(words) {
    mapped <- lookup[words]
    unname(mapped[mapped > 0])
  })
  original.vocab <- original.vocab[observed]
  V <- length(observed)
}

random.beta <- matrix(rgamma(K * V, shape=0.2), nrow=K, ncol=V)
random.beta <- random.beta/rowSums(random.beta)
initial.beta <- (1 - options$initialization_noise) * beta +
  options$initialization_noise * random.beta
initial.beta <- initial.beta/rowSums(initial.beta)

corpus <- list(
  documents=documents,
  vocab=original.vocab,
  initial.logbeta=list(log(initial.beta)),
  truth=list(beta=beta, theta=theta, anchors=anchors),
  settings=options,
  tokens=total.tokens
)
saveRDS(corpus, options$output, compress=FALSE)
cat(sprintf(
  "Saved %s: N=%d, K=%d, V=%d, tokens=%d (%.1f tokens/document)\n",
  options$output, N, K, V, total.tokens, total.tokens/N
))
