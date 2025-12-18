# ============================================================
# UCI: "News Popularity in Multiple Social Media Platforms"
# ~93k documents with Title + Headline + metadata
#
# Requires: stm, data.table, lubridate, stringi
# install.packages(c("stm", "data.table", "lubridate", "stringi"))
# ============================================================

suppressPackageStartupMessages({
  library(stm)
  library(data.table)
  library(lubridate)
  library(stringi)
})

# ---- 1) Download + unzip ----------------------------------------------------
DATA_DIR <- file.path(tempdir(), "uci_news_popularity_432")
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

zip_url  <- "https://archive.ics.uci.edu/static/public/432/news%2Bpopularity%2Bin%2Bmultiple%2Bsocial%2Bmedia%2Bplatforms.zip"
zip_path <- file.path(DATA_DIR, "uci_432_news_popularity.zip")

if (!file.exists(zip_path)) {
  message("Downloading ZIP to: ", zip_path)
  download.file(zip_url, destfile = zip_path, mode = "wb", quiet = FALSE)
}

message("Unzipping into: ", DATA_DIR)
unzip(zip_path, exdir = DATA_DIR)

news_csv <- file.path(DATA_DIR, "Data", "News_Final.csv")
if (!file.exists(news_csv)) {
  stop("Could not find News_Final.csv after unzip. Look inside: ", DATA_DIR)
}

# ---- 2) Load data -----------------------------------------------------------
dt <- fread(news_csv, encoding = "UTF-8", na.strings = c("", "NA"))
message("Rows loaded: ", nrow(dt))

# The dataset provides: Title, Headline, Source, Topic, PublishDate, etc.
# We'll combine Title + Headline as the document text.

dt[, Title    := fifelse(is.na(Title),    "", Title)]
dt[, Headline := fifelse(is.na(Headline), "", Headline)]
dt[, text := paste(Title, Headline, sep = " . ")]

# Light normalization (optional but helpful for messy encodings)
dt[, text := stri_trans_general(text, "Latin-ASCII")]
dt[, text := gsub("\\s+", " ", text)]
dt <- dt[nchar(text) > 0]

# ---- 3) Build metadata ------------------------------------------------------
dt[, topic  := factor(Topic)]
dt[, source := factor(Source)]

# Parse PublishDate robustly (format can vary across datasets)
dt[, publish_time := suppressWarnings(parse_date_time(
  PublishDate,
  orders = c(
    "Ymd HMS", "Y-m-d H:M:S",
    "Ymd HM",  "Y-m-d H:M",
    "mdY HMS", "m/d/Y H:M:S",
    "mdY HM",  "m/d/Y H:M"
  ),
  tz = "UTC"
))]

# Fallback: try base parsing if lubridate didn't catch it
dt[is.na(publish_time), publish_time := as.POSIXct(PublishDate, tz = "UTC")]

# Drop rows where date parsing failed (should be rare)
dt <- dt[!is.na(publish_time)]

# Numeric day index is convenient for spline terms
dt[, day := as.numeric(as.Date(publish_time))]

meta <- dt[, .(
  topic,
  source,
  day,
  publish_time,
  Facebook,
  GooglePlus,
  LinkedIn
)]

# OPTIONAL: if you want a quick smoke test, uncomment a subsample:
# set.seed(1)
# dt   <- dt[sample(.N, 20000)]
# meta <- meta[dt$.I]  # keep aligned

# ---- 4) Preprocess for STM --------------------------------------------------
meta <- as.data.frame(meta)
processed <- textProcessor(
  documents = dt$text,
  metadata  = meta,
  lowercase        = TRUE,
  removestopwords  = TRUE,
  removenumbers    = TRUE,
  removepunctuation= TRUE,
  stem             = TRUE
)

# Drop rare terms and empty docs
out <- prepDocuments(
  processed$documents,
  processed$vocab,
  processed$meta,
  lower.thresh = 10   # increase to 20/50 for speed + smaller vocab
)

docs  <- out$documents
vocab <- out$vocab
meta  <- out$meta

message("Docs after prep: ", length(docs))
message("Vocab size: ", length(vocab))

# ---- 5) Fit STM -------------------------------------------------------------
# prevalence uses topic + a smooth over time (day)
set.seed(123)
K <- 20

fit <- stm_svi(
  documents  = docs,
  vocab      = vocab,
  K          = K,
  patience = 100
)

#trajectory of per word bound
fit$convergence$bound/sum(fit$settings$dim$wcounts$x)
