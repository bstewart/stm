// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace {

inline int document_ncol(SEXP document) {
  SEXP dimensions = Rf_getAttrib(document, R_DimSymbol);
  if (Rf_isNull(dimensions) || Rf_length(dimensions) != 2 ||
      INTEGER(dimensions)[0] != 2) {
    Rcpp::stop("Each document must be a two-row matrix");
  }
  return INTEGER(dimensions)[1];
}

inline double document_value(SEXP document, int index) {
  if (TYPEOF(document) == INTSXP) {
    return static_cast<double>(INTEGER(document)[index]);
  }
  if (TYPEOF(document) == REALSXP) {
    return REAL(document)[index];
  }
  Rcpp::stop("Document matrices must contain integer or numeric values");
  return 0.0;
}

inline std::uint64_t splitmix64(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31);
}

// Euclidean projection of z onto the probability simplex (Duchi, Shalev-Shwartz,
// Singer and Chandra 2008).  Writes the projection into out.
void project_simplex(const std::vector<double>& z, std::vector<double>& sorted,
                     std::vector<double>& out) {
  const int n = z.size();
  sorted.assign(z.begin(), z.end());
  std::sort(sorted.begin(), sorted.end(), std::greater<double>());
  double cumulative = 0.0;
  double theta = 0.0;
  for (int j = 0; j < n; ++j) {
    cumulative += sorted[j];
    const double candidate = (cumulative - 1.0) / static_cast<double>(j + 1);
    if (sorted[j] - candidate > 0.0) theta = candidate;
  }
  for (int j = 0; j < n; ++j) out[j] = std::max(z[j] - theta, 0.0);
}

// Exact solution of the RecoverL2 problem for one word,
//   minimize 0.5 * alpha' G alpha - b' alpha  subject to alpha in the simplex,
// by accelerated projected gradient (FISTA with adaptive restart) using the
// step 1/lambda_max(G).  G = anchor_gram and b = the word's anchor kernel.
// Used only when the legacy exponentiated-gradient iteration fails to meet its
// tolerance, which happens when the anchor rows are peaked enough that the
// fixed exponentiated-gradient step of 50 overshoots.  Returns true when the
// iterates stabilized to within the requested parameter tolerance.
bool solve_simplex_qp(const arma::mat& gram, const double* linear,
                      double inverse_lipschitz, std::vector<double>& alpha,
                      int maximum_iterations, int& iterations_used,
                      std::vector<double>& previous, std::vector<double>& extrapolated,
                      std::vector<double>& gradient, std::vector<double>& scratch,
                      std::vector<double>& candidate) {
  const int topics = alpha.size();
  const double parameter_tolerance = 1e-10;
  previous.assign(alpha.begin(), alpha.end());
  double momentum_old = 1.0;
  bool stabilized = false;
  int iteration = 0;
  while (iteration < maximum_iterations) {
    ++iteration;
    const double momentum_new =
        0.5 * (1.0 + std::sqrt(1.0 + 4.0 * momentum_old * momentum_old));
    const double weight = (momentum_old - 1.0) / momentum_new;
    for (int k = 0; k < topics; ++k) {
      extrapolated[k] = alpha[k] + weight * (alpha[k] - previous[k]);
    }
    for (int k = 0; k < topics; ++k) {
      double value = -linear[k];
      for (int j = 0; j < topics; ++j) value += gram(j, k) * extrapolated[j];
      gradient[k] = value;
    }
    for (int k = 0; k < topics; ++k) {
      scratch[k] = extrapolated[k] - inverse_lipschitz * gradient[k];
    }
    project_simplex(scratch, gradient, candidate);
    // Adaptive restart (O'Donoghue and Candes 2015): drop momentum when the
    // step turns against the previous direction of travel.
    double alignment = 0.0;
    double movement = 0.0;
    for (int k = 0; k < topics; ++k) {
      const double step = candidate[k] - alpha[k];
      alignment += (extrapolated[k] - candidate[k]) * step;
      movement = std::max(movement, std::abs(step));
    }
    previous.assign(alpha.begin(), alpha.end());
    alpha.assign(candidate.begin(), candidate.end());
    if (alignment > 0.0) {
      momentum_old = 1.0;
      previous.assign(alpha.begin(), alpha.end());
    } else {
      momentum_old = momentum_new;
    }
    if (movement < parameter_tolerance) {
      stabilized = true;
      break;
    }
  }
  iterations_used = iteration;
  return stabilized;
}

Rcpp::List recover_from_kernels(
    const arma::mat& word_anchor,
    const arma::mat& anchor_gram,
    const std::vector<int>& anchor_topic,
    Rcpp::NumericVector word_probability,
    const std::vector<unsigned char>& valid_word,
    double tolerance,
    int maximum_iterations,
    bool exact_only) {
  const int vocabulary_size = word_anchor.n_rows;
  const int topics = word_anchor.n_cols;
  Rcpp::NumericMatrix weights(vocabulary_size, topics);
  Rcpp::IntegerVector iterations(vocabulary_size);
  Rcpp::LogicalVector converged(vocabulary_size);
  Rcpp::LogicalVector fallback(vocabulary_size);
  std::vector<double> alpha(topics);
  std::vector<double> gradient(topics);
  std::vector<double> topic_probability(topics, 0.0);
  std::vector<double> qp_previous(topics), qp_extrapolated(topics),
      qp_gradient(topics), qp_scratch(topics), qp_candidate(topics),
      qp_linear(topics);
  const double uniform_weight = 1.0 / static_cast<double>(topics);
  const double step_size = 50.0;
  const double machine_epsilon = std::numeric_limits<double>::epsilon();
  // Lipschitz constant of the quadratic objective for the exact fallback.
  double inverse_lipschitz = 0.0;
  {
    arma::vec eigenvalues;
    if (arma::eig_sym(eigenvalues, anchor_gram) && eigenvalues.n_elem > 0 &&
        eigenvalues(eigenvalues.n_elem - 1) > 0.0) {
      inverse_lipschitz = 1.0 / eigenvalues(eigenvalues.n_elem - 1);
    }
  }
  const int fallback_iterations = 20000;

  for (int word = 0; word < vocabulary_size; ++word) {
    if ((word & 255) == 0) Rcpp::checkUserInterrupt();
    if (!valid_word[word]) continue;

    if (anchor_topic[word] >= 0) {
      weights(word, anchor_topic[word]) = 1.0;
      iterations[word] = 0;
      converged[word] = true;
    } else {
      std::fill(alpha.begin(), alpha.end(), uniform_weight);
      double previous_sse_root = std::numeric_limits<double>::infinity();
      bool word_converged = false;
      int iteration = exact_only ? 0 : 1;

      while (!exact_only) {
        double squared_error = 0.0;
        for (int topic = 0; topic < topics; ++topic) {
          double fitted = 0.0;
          for (int other = 0; other < topics; ++other) {
            fitted += alpha[other] * anchor_gram(other, topic);
          }
          gradient[topic] = word_anchor(word, topic) - fitted;
          squared_error += gradient[topic] * gradient[topic];
        }

        double maximum_derivative = -std::numeric_limits<double>::infinity();
        for (int topic = 0; topic < topics; ++topic) {
          gradient[topic] *= 2.0 * step_size;
          maximum_derivative =
              std::max(maximum_derivative, gradient[topic]);
        }

        double normalizer = 0.0;
        for (int topic = 0; topic < topics; ++topic) {
          alpha[topic] *= std::exp(gradient[topic] - maximum_derivative);
          normalizer += alpha[topic];
        }
        for (int topic = 0; topic < topics; ++topic) {
          alpha[topic] /= normalizer;
        }

        const double sse_root = std::sqrt(squared_error);
        word_converged =
            std::abs(previous_sse_root - sse_root) < tolerance;
        if (word_converged || iteration == maximum_iterations) break;
        previous_sse_root = sse_root;
        ++iteration;
      }

      if (!word_converged && inverse_lipschitz > 0.0) {
        // Either the exact solver was requested for every word, or the
        // fixed-step exponentiated gradient did not settle (typically a limit
        // cycle when the anchor rows are peaked).  Solve the same convex
        // problem exactly from the uniform start instead of returning the
        // arbitrary final iterate.  Words that did converge above keep the
        // legacy solution untouched.
        std::fill(alpha.begin(), alpha.end(), uniform_weight);
        for (int topic = 0; topic < topics; ++topic) {
          qp_linear[topic] = word_anchor(word, topic);
        }
        int qp_iterations = 0;
        word_converged = solve_simplex_qp(
            anchor_gram, qp_linear.data(), inverse_lipschitz, alpha,
            fallback_iterations, qp_iterations, qp_previous, qp_extrapolated,
            qp_gradient, qp_scratch, qp_candidate);
        iteration += qp_iterations;
        fallback[word] = true;
      }

      for (int topic = 0; topic < topics; ++topic) {
        weights(word, topic) =
            alpha[topic] > 0.0 ? alpha[topic] : machine_epsilon;
      }
      iterations[word] = iteration;
      converged[word] = word_converged;
    }

    const double probability = word_probability[word];
    for (int topic = 0; topic < topics; ++topic) {
      topic_probability[topic] += weights(word, topic) * probability;
    }
  }

  Rcpp::NumericMatrix beta(topics, vocabulary_size);
  for (int word = 0; word < vocabulary_size; ++word) {
    if (!valid_word[word]) continue;
    const double probability = word_probability[word];
    for (int topic = 0; topic < topics; ++topic) {
      beta(topic, word) =
          weights(word, topic) * probability / topic_probability[topic];
    }
  }

  return Rcpp::List::create(Rcpp::Named("A") = beta,
                            Rcpp::Named("weights") = weights,
                            Rcpp::Named("iterations") = iterations,
                            Rcpp::Named("converged") = converged,
                            Rcpp::Named("fallback") = fallback);
}

}  // namespace

// Directly construct a row-normalized random projection of the second-order
// moment.  The exact V-by-V moment is never materialized.  A fixed hash-based
// sparse sign projection makes the result reproducible without touching R's
// random-number stream.
//
// [[Rcpp::export]]
Rcpp::List spectralProjectionCpp(Rcpp::List documents,
                                 int vocabulary_size,
                                 int projection_dimension,
                                 double projection_sparsity = 0.05,
                                 int seed = 1) {
  if (vocabulary_size < 1 || projection_dimension < 1) {
    Rcpp::stop("Vocabulary and projection dimensions must be positive");
  }
  if (!(projection_sparsity > 0.0 && projection_sparsity <= 1.0)) {
    Rcpp::stop("projection_sparsity must be in (0, 1]");
  }

  const int number_documents = documents.size();
  const int nonzeros = std::min(
      projection_dimension,
      std::max(1, static_cast<int>(
                      std::ceil(projection_sparsity * projection_dimension))));
  const double projection_scale = 1.0 / std::sqrt(nonzeros);

  bool zero_indexed = false;
  for (int document_index = 0; document_index < number_documents;
       ++document_index) {
    SEXP document = documents[document_index];
    const int unique_terms = document_ncol(document);
    for (int term = 0; term < unique_terms; ++term) {
      if (document_value(document, 2 * term) == 0.0) {
        zero_indexed = true;
        break;
      }
    }
    if (zero_indexed) break;
  }
  const int index_offset = zero_indexed ? 0 : 1;

  std::vector<int> projection_index(
      static_cast<std::size_t>(vocabulary_size) * nonzeros);
  std::vector<double> projection_value(
      static_cast<std::size_t>(vocabulary_size) * nonzeros);
  const std::uint64_t unsigned_seed = static_cast<std::uint32_t>(seed);
  for (int word = 0; word < vocabulary_size; ++word) {
    for (int item = 0; item < nonzeros; ++item) {
      std::uint64_t state =
          unsigned_seed ^
          (static_cast<std::uint64_t>(word + 1) * 0xd2b74407b1ce6e93ULL) ^
          (static_cast<std::uint64_t>(item + 1) * 0xca5a826395121157ULL);
      std::uint64_t hashed = splitmix64(state);
      int coordinate = static_cast<int>(hashed % projection_dimension);
      bool duplicate = true;
      while (duplicate) {
        duplicate = false;
        for (int previous = 0; previous < item; ++previous) {
          if (projection_index[static_cast<std::size_t>(word) * nonzeros +
                               previous] == coordinate) {
            duplicate = true;
            hashed = splitmix64(hashed);
            coordinate = static_cast<int>(hashed % projection_dimension);
            break;
          }
        }
      }
      const std::size_t location =
          static_cast<std::size_t>(word) * nonzeros + item;
      projection_index[location] = coordinate;
      projection_value[location] =
          ((splitmix64(hashed) & 1ULL) ? 1.0 : -1.0) * projection_scale;
    }
  }

  arma::mat projected(projection_dimension, vocabulary_size, arma::fill::zeros);
  arma::vec sketched_counts(projection_dimension, arma::fill::zeros);
  Rcpp::NumericVector row_normalizer(vocabulary_size);
  Rcpp::NumericVector word_probability(vocabulary_size);
  double total_tokens = 0.0;

  for (int document_index = 0; document_index < number_documents;
       ++document_index) {
    if ((document_index & 63) == 0) Rcpp::checkUserInterrupt();
    SEXP document = documents[document_index];
    const int unique_terms = document_ncol(document);
    double document_length = 0.0;

    for (int term = 0; term < unique_terms; ++term) {
      const int word = static_cast<int>(document_value(document, 2 * term)) -
                       index_offset;
      const double count = document_value(document, 2 * term + 1);
      if (word < 0 || word >= vocabulary_size) {
        Rcpp::stop("Document contains a word index outside the vocabulary");
      }
      if (!std::isfinite(count) || count < 0.0) {
        Rcpp::stop("Document counts must be finite and non-negative");
      }
      word_probability[word] += count;
      document_length += count;
    }
    total_tokens += document_length;
    if (document_length < 2.0) continue;

    sketched_counts.zeros();
    for (int term = 0; term < unique_terms; ++term) {
      const int word = static_cast<int>(document_value(document, 2 * term)) -
                       index_offset;
      const double count = document_value(document, 2 * term + 1);
      const std::size_t projection_start =
          static_cast<std::size_t>(word) * nonzeros;
      for (int item = 0; item < nonzeros; ++item) {
        const std::size_t location = projection_start + item;
        sketched_counts[projection_index[location]] +=
            count * projection_value[location];
      }
    }

    const double inverse_divisor =
        1.0 / (document_length * (document_length - 1.0));
    for (int term = 0; term < unique_terms; ++term) {
      const int word = static_cast<int>(document_value(document, 2 * term)) -
                       index_offset;
      const double count = document_value(document, 2 * term + 1);
      const double multiplier = count * inverse_divisor;
      double* output = projected.colptr(word);
      for (int coordinate = 0; coordinate < projection_dimension;
           ++coordinate) {
        output[coordinate] += multiplier * sketched_counts[coordinate];
      }
      const std::size_t projection_start =
          static_cast<std::size_t>(word) * nonzeros;
      for (int item = 0; item < nonzeros; ++item) {
        const std::size_t location = projection_start + item;
        output[projection_index[location]] -=
            multiplier * projection_value[location];
      }
      row_normalizer[word] +=
          count * (document_length - 1.0) * inverse_divisor;
    }
  }

  if (!(total_tokens > 0.0)) Rcpp::stop("The corpus contains no tokens");
  for (int word = 0; word < vocabulary_size; ++word) {
    word_probability[word] /= total_tokens;
    if (row_normalizer[word] > 0.0) {
      projected.col(word) /= row_normalizer[word];
    }
  }

  return Rcpp::List::create(
      Rcpp::Named("Qbar") = projected.t(),
      Rcpp::Named("wprob") = word_probability,
      Rcpp::Named("row_normalizer") = row_normalizer,
      Rcpp::Named("dimension") = projection_dimension,
      Rcpp::Named("nonzeros") = nonzeros,
      Rcpp::Named("seed") = seed);
}

// Construct the same second-order moment as gram(), directly from stm's
// document representation.  This avoids building and transforming several
// sparse Matrix objects before the final dense result is needed.
//
// When keep is supplied (1-based vocabulary indices, in the order the caller
// wants the moment's rows and columns), the moment is restricted to those
// words exactly as the legacy code restricted the document-term matrix before
// calling gram(): document lengths count kept words only and documents with
// fewer than two kept tokens are skipped.  Word probabilities are always
// returned for the full vocabulary, computed over all tokens.
//
// [[Rcpp::export]]
Rcpp::List spectralGramCpp(Rcpp::List documents, int vocabulary_size,
                           Rcpp::Nullable<Rcpp::IntegerVector> keep = R_NilValue) {
  if (vocabulary_size < 1) {
    Rcpp::stop("vocabulary_size must be positive");
  }

  // Map each vocabulary word to its row in the moment (-1 = not kept).
  std::vector<int> position(vocabulary_size, -1);
  int moment_size = vocabulary_size;
  if (keep.isNotNull()) {
    Rcpp::IntegerVector kept(keep);
    moment_size = kept.size();
    if (moment_size < 1) Rcpp::stop("keep must contain at least one word");
    for (int j = 0; j < moment_size; ++j) {
      const int word = kept[j] - 1;
      if (word < 0 || word >= vocabulary_size) {
        Rcpp::stop("keep contains an index outside the vocabulary");
      }
      if (position[word] != -1) Rcpp::stop("keep must not contain duplicates");
      position[word] = j;
    }
  } else {
    for (int word = 0; word < vocabulary_size; ++word) position[word] = word;
  }

  const int number_documents = documents.size();
  bool zero_indexed = false;
  for (int document_index = 0; document_index < number_documents;
       ++document_index) {
    SEXP document = documents[document_index];
    const int unique_terms = document_ncol(document);
    for (int term = 0; term < unique_terms; ++term) {
      if (document_value(document, 2 * term) == 0.0) {
        zero_indexed = true;
        break;
      }
    }
    if (zero_indexed) break;
  }

  const int index_offset = zero_indexed ? 0 : 1;
  Rcpp::NumericMatrix moment(moment_size, moment_size);
  Rcpp::NumericVector word_probability(vocabulary_size);
  double total_tokens = 0.0;
  std::vector<int> rows;
  std::vector<double> counts;

  for (int document_index = 0; document_index < number_documents;
       ++document_index) {
    if ((document_index & 63) == 0) Rcpp::checkUserInterrupt();

    SEXP document = documents[document_index];
    const int unique_terms = document_ncol(document);
    double document_length = 0.0;
    rows.clear();
    counts.clear();

    for (int term = 0; term < unique_terms; ++term) {
      const int word = static_cast<int>(document_value(document, 2 * term)) -
                       index_offset;
      const double count = document_value(document, 2 * term + 1);
      if (word < 0 || word >= vocabulary_size) {
        Rcpp::stop("Document contains a word index outside the vocabulary");
      }
      if (!std::isfinite(count) || count < 0.0) {
        Rcpp::stop("Document counts must be finite and non-negative");
      }
      word_probability[word] += count;
      total_tokens += count;
      if (position[word] >= 0) {
        rows.push_back(position[word]);
        counts.push_back(count);
        document_length += count;
      }
    }

    if (document_length < 2.0) continue;
    const double inverse_divisor =
        1.0 / (document_length * (document_length - 1.0));
    const int kept_terms = rows.size();

    for (int first = 0; first < kept_terms; ++first) {
      const int first_row = rows[first];
      const double first_count = counts[first];
      moment(first_row, first_row) +=
          first_count * (first_count - 1.0) * inverse_divisor;

      for (int second = first + 1; second < kept_terms; ++second) {
        const int second_row = rows[second];
        const double contribution =
            first_count * counts[second] * inverse_divisor;
        moment(first_row, second_row) += contribution;
        moment(second_row, first_row) += contribution;
      }
    }
  }

  if (!(total_tokens > 0.0)) Rcpp::stop("The corpus contains no tokens");
  for (int word = 0; word < vocabulary_size; ++word) {
    word_probability[word] /= total_tokens;
  }

  return Rcpp::List::create(Rcpp::Named("Q") = moment,
                            Rcpp::Named("wprob") = word_probability);
}

// Find anchors with the same stabilized modified Gram-Schmidt procedure as
// fastAnchor().  Vocabulary points are stored in contiguous columns and the
// projection is applied in place, avoiding a V-by-q outer-product temporary.
//
// [[Rcpp::export]]
Rcpp::IntegerVector fastAnchorCpp(SEXP qbar_sexp, int topics) {
  Rcpp::NumericMatrix qbar_matrix(qbar_sexp);
  const int vocabulary_size = qbar_matrix.nrow();
  const int dimension = qbar_matrix.ncol();
  if (topics < 1 || topics > vocabulary_size) {
    Rcpp::stop("topics must be between one and the number of words");
  }

  arma::mat qbar(qbar_matrix.begin(), vocabulary_size, dimension, false, true);
  arma::mat points = qbar.t();
  std::vector<double> squared_norm(vocabulary_size, 0.0);
  std::vector<unsigned char> selected(vocabulary_size, 0);
  Rcpp::IntegerVector anchors(topics);

  for (int word = 0; word < vocabulary_size; ++word) {
    const double* point = points.colptr(word);
    double norm = 0.0;
    for (int coordinate = 0; coordinate < dimension; ++coordinate) {
      norm += point[coordinate] * point[coordinate];
    }
    squared_norm[word] = norm;
  }

  for (int topic = 0; topic < topics; ++topic) {
    if ((topic & 3) == 0) Rcpp::checkUserInterrupt();

    int anchor = 0;
    double maximum_norm = squared_norm[0];
    for (int word = 1; word < vocabulary_size; ++word) {
      if (squared_norm[word] > maximum_norm) {
        maximum_norm = squared_norm[word];
        anchor = word;
      }
    }
    if (!(maximum_norm > 0.0) || !std::isfinite(maximum_norm)) {
      Rcpp::stop("Anchor search encountered a degenerate residual matrix");
    }

    anchors[topic] = anchor + 1;
    selected[anchor] = 1;
    squared_norm[anchor] = 0.0;

    double* basis = points.colptr(anchor);
    const double inverse_norm = 1.0 / std::sqrt(maximum_norm);
    for (int coordinate = 0; coordinate < dimension; ++coordinate) {
      basis[coordinate] *= inverse_norm;
    }

    for (int word = 0; word < vocabulary_size; ++word) {
      if (selected[word]) continue;
      double* point = points.colptr(word);
      double inner_product = 0.0;
      for (int coordinate = 0; coordinate < dimension; ++coordinate) {
        inner_product += point[coordinate] * basis[coordinate];
      }
      double norm = 0.0;
      for (int coordinate = 0; coordinate < dimension; ++coordinate) {
        point[coordinate] -= inner_product * basis[coordinate];
        norm += point[coordinate] * point[coordinate];
      }
      squared_norm[word] = norm;
    }
  }

  return anchors;
}

// Recover all word-topic weights at once for a fixed representation and anchor
// set.  Qbar %*% t(X), the vocabulary-dependent part of every recovery
// problem, is shared across words; only the small K-by-K exponentiated-gradient
// iterations remain word-specific.
//
// [[Rcpp::export]]
Rcpp::List recoverL2Cpp(SEXP qbar_sexp,
                        Rcpp::IntegerVector anchors,
                        Rcpp::NumericVector word_probability,
                        double tolerance = 1e-7,
                        int maximum_iterations = 500,
                        bool exact_only = false) {
  Rcpp::NumericMatrix qbar_matrix(qbar_sexp);
  const int vocabulary_size = qbar_matrix.nrow();
  const int dimension = qbar_matrix.ncol();
  const int topics = anchors.size();
  if (topics < 1 || topics > vocabulary_size) {
    Rcpp::stop("anchors must contain between one and V entries");
  }
  if (word_probability.size() != vocabulary_size) {
    Rcpp::stop("word_probability must have one value per word");
  }
  if (!(tolerance >= 0.0) || maximum_iterations < 1) {
    Rcpp::stop("Invalid recovery convergence controls");
  }

  arma::mat qbar(qbar_matrix.begin(), vocabulary_size, dimension, false, true);
  arma::uvec anchor_indices(topics);
  std::vector<int> anchor_topic(vocabulary_size, -1);
  for (int topic = 0; topic < topics; ++topic) {
    const int anchor = anchors[topic] - 1;
    if (anchor < 0 || anchor >= vocabulary_size) {
      Rcpp::stop("Anchor index is outside the vocabulary");
    }
    if (anchor_topic[anchor] != -1) {
      Rcpp::stop("Anchor indices must be unique");
    }
    anchor_indices[topic] = anchor;
    anchor_topic[anchor] = topic;
  }

  const arma::mat anchor_rows = qbar.rows(anchor_indices);
  const arma::mat anchor_gram = anchor_rows * anchor_rows.t();
  const arma::mat word_anchor = qbar * anchor_rows.t();
  const std::vector<unsigned char> valid_word(vocabulary_size, 1);
  return recover_from_kernels(word_anchor, anchor_gram, anchor_topic,
                              word_probability, valid_word, tolerance,
                              maximum_iterations, exact_only);
}

// Use the projected representation only to select anchors, then solve the
// original RecoverL2 problem exactly for those anchors.  The necessary
// word-anchor kernels are accumulated in two sparse corpus passes, so the
// V-by-V moment is still avoided.
//
// [[Rcpp::export]]
Rcpp::List recoverL2OriginalCpp(Rcpp::List documents,
                                int vocabulary_size,
                                Rcpp::IntegerVector anchors,
                                Rcpp::NumericVector row_normalizer,
                                Rcpp::NumericVector word_probability,
                                double tolerance = 1e-7,
                                int maximum_iterations = 500,
                                bool exact_only = false) {
  const int number_documents = documents.size();
  const int topics = anchors.size();
  if (vocabulary_size < 1 || topics < 1 || topics > vocabulary_size) {
    Rcpp::stop("Invalid vocabulary size or anchor count");
  }
  if (row_normalizer.size() != vocabulary_size ||
      word_probability.size() != vocabulary_size) {
    Rcpp::stop("Moment normalizers and probabilities must have length V");
  }
  if (!(tolerance >= 0.0) || maximum_iterations < 1) {
    Rcpp::stop("Invalid recovery convergence controls");
  }

  bool zero_indexed = false;
  for (int document_index = 0; document_index < number_documents;
       ++document_index) {
    SEXP document = documents[document_index];
    const int unique_terms = document_ncol(document);
    for (int term = 0; term < unique_terms; ++term) {
      if (document_value(document, 2 * term) == 0.0) {
        zero_indexed = true;
        break;
      }
    }
    if (zero_indexed) break;
  }
  const int index_offset = zero_indexed ? 0 : 1;

  std::vector<int> anchor_topic(vocabulary_size, -1);
  for (int topic = 0; topic < topics; ++topic) {
    const int anchor = anchors[topic] - 1;
    if (anchor < 0 || anchor >= vocabulary_size ||
        !(row_normalizer[anchor] > 0.0)) {
      Rcpp::stop("Anchor is outside the valid moment rows");
    }
    if (anchor_topic[anchor] != -1) {
      Rcpp::stop("Anchor indices must be unique");
    }
    anchor_topic[anchor] = topic;
  }

  arma::mat anchor_rows(topics, vocabulary_size, arma::fill::zeros);
  for (int document_index = 0; document_index < number_documents;
       ++document_index) {
    if ((document_index & 127) == 0) Rcpp::checkUserInterrupt();
    SEXP document = documents[document_index];
    const int unique_terms = document_ncol(document);
    double document_length = 0.0;
    for (int term = 0; term < unique_terms; ++term) {
      document_length += document_value(document, 2 * term + 1);
    }
    if (document_length < 2.0) continue;
    const double inverse_divisor =
        1.0 / (document_length * (document_length - 1.0));

    for (int anchor_term = 0; anchor_term < unique_terms; ++anchor_term) {
      const int anchor_word =
          static_cast<int>(document_value(document, 2 * anchor_term)) -
          index_offset;
      if (anchor_word < 0 || anchor_word >= vocabulary_size) {
        Rcpp::stop("Document contains a word index outside the vocabulary");
      }
      const int topic = anchor_topic[anchor_word];
      if (topic < 0) continue;
      const double anchor_count =
          document_value(document, 2 * anchor_term + 1);
      const double multiplier = anchor_count * inverse_divisor;
      for (int term = 0; term < unique_terms; ++term) {
        const int word =
            static_cast<int>(document_value(document, 2 * term)) - index_offset;
        if (word < 0 || word >= vocabulary_size) {
          Rcpp::stop("Document contains a word index outside the vocabulary");
        }
        anchor_rows(topic, word) +=
            multiplier * document_value(document, 2 * term + 1);
      }
      anchor_rows(topic, anchor_word) -= multiplier;
    }
  }
  for (int topic = 0; topic < topics; ++topic) {
    anchor_rows.row(topic) /= row_normalizer[anchors[topic] - 1];
  }

  const arma::mat anchor_gram = anchor_rows * anchor_rows.t();
  arma::mat word_anchor(vocabulary_size, topics, arma::fill::zeros);
  std::vector<double> document_anchor(topics, 0.0);
  for (int document_index = 0; document_index < number_documents;
       ++document_index) {
    if ((document_index & 127) == 0) Rcpp::checkUserInterrupt();
    SEXP document = documents[document_index];
    const int unique_terms = document_ncol(document);
    double document_length = 0.0;
    for (int term = 0; term < unique_terms; ++term) {
      document_length += document_value(document, 2 * term + 1);
    }
    if (document_length < 2.0) continue;

    std::fill(document_anchor.begin(), document_anchor.end(), 0.0);
    for (int term = 0; term < unique_terms; ++term) {
      const int word =
          static_cast<int>(document_value(document, 2 * term)) - index_offset;
      if (word < 0 || word >= vocabulary_size) {
        Rcpp::stop("Document contains a word index outside the vocabulary");
      }
      const double count = document_value(document, 2 * term + 1);
      for (int topic = 0; topic < topics; ++topic) {
        document_anchor[topic] += count * anchor_rows(topic, word);
      }
    }

    const double inverse_divisor =
        1.0 / (document_length * (document_length - 1.0));
    for (int term = 0; term < unique_terms; ++term) {
      const int word =
          static_cast<int>(document_value(document, 2 * term)) - index_offset;
      if (word < 0 || word >= vocabulary_size) {
        Rcpp::stop("Document contains a word index outside the vocabulary");
      }
      const double multiplier =
          document_value(document, 2 * term + 1) * inverse_divisor;
      for (int topic = 0; topic < topics; ++topic) {
        word_anchor(word, topic) +=
            multiplier *
            (document_anchor[topic] - anchor_rows(topic, word));
      }
    }
  }

  std::vector<unsigned char> valid_word(vocabulary_size, 0);
  for (int word = 0; word < vocabulary_size; ++word) {
    if (row_normalizer[word] > 0.0) {
      word_anchor.row(word) /= row_normalizer[word];
      valid_word[word] = 1;
    }
  }

  return recover_from_kernels(word_anchor, anchor_gram, anchor_topic,
                              word_probability, valid_word, tolerance,
                              maximum_iterations, exact_only);
}
