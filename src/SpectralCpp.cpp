// [[Rcpp::depends(RcppArmadillo)]]

#include "RcppArmadillo.h"

// Fast Anchor Word Selection (C++ implementation)
//
// Implements the Gram-Schmidt-like anchor word selection algorithm
// from Arora et al. (2013) with optimized C++ operations.
//
// @param Qbar Row-normalized gram matrix (V x V), modified in-place
// @param K Number of topics (anchors to find)
// @param tol Numerical tolerance (unused, kept for compatibility)
// @param verbose Print progress dots
// @return List with 'basis' (anchor indices, 1-indexed) and 'Qbar' (modified matrix)
// [[Rcpp::export]]
Rcpp::List fastAnchorCpp(arma::mat Qbar, int K, double tol = 1e-3, bool verbose = true) {
  int V = Qbar.n_rows;
  arma::uvec basis(K);
  arma::vec rowSquaredSums = arma::sum(arma::square(Qbar), 1);

  for(int i = 0; i < K; i++) {
    // Find row with maximum squared sum
    arma::uword max_idx = rowSquaredSums.index_max();
    basis(i) = max_idx;

    // Normalize the selected row
    double maxval = rowSquaredSums(max_idx);
    double normalizer = 1.0 / std::sqrt(maxval);
    Qbar.row(max_idx) *= normalizer;

    // Compute inner products: Qbar %*% Qbar[basis[i],]
    arma::vec innerproducts = Qbar * Qbar.row(max_idx).t();

    // Project: innerproducts %o% Qbar[basis[i],]
    // This is an outer product creating a V x V matrix
    arma::mat project = innerproducts * Qbar.row(max_idx);

    // Zero out rows corresponding to already-selected basis vectors
    for(int j = 0; j <= i; j++) {
      project.row(basis(j)).zeros();
    }

    // Subtract projection
    Qbar -= project;

    // Update squared sums
    rowSquaredSums = arma::sum(arma::square(Qbar), 1);

    // Zero out basis entries
    for(int j = 0; j <= i; j++) {
      rowSquaredSums(basis(j)) = 0.0;
    }

    if(verbose) Rcpp::Rcout << ".";
  }

  if(verbose) Rcpp::Rcout << std::endl;

  // Convert to 1-indexed for R
  arma::uvec basis_r = basis + 1;

  return Rcpp::List::create(
    Rcpp::Named("basis") = basis_r,
    Rcpp::Named("Qbar") = Qbar
  );
}

// Exponentiated Gradient Descent (C++ implementation)
//
// Solves for optimal convex combination of anchor words to recover
// a non-anchor word using exponentiated gradient descent.
//
// @param X K x V matrix of anchor word rows from Qbar
// @param y V x 1 target vector (word to recover)
// @param XtX K x K precomputed X * X^T
// @param eta Learning rate (fixed)
// @param maxiter Maximum iterations
// @param rtol Convergence tolerance
// @return K x 1 vector of mixing weights (sums to 1)
// [[Rcpp::export]]
arma::vec expgradCpp(const arma::mat& X, const arma::vec& y, const arma::mat& XtX,
                     double eta = 50.0, int maxiter = 500, double rtol = 1e-7) {
  int K = X.n_rows;

  // Initialize alpha uniformly
  arma::rowvec alpha(K);
  alpha.fill(1.0 / K);

  // Precompute y^T X
  arma::rowvec ytX = y.t() * X.t();

  bool converged = false;
  double sse_old = arma::datum::inf;
  int its = 0;

  while(!converged && its < maxiter) {
    // Compute gradient: y^T X - alpha X^T X
    arma::rowvec grad = ytX - alpha * XtX;

    // Sum of squared errors (for convergence check)
    double sse = arma::accu(arma::square(grad));

    // Scale gradient
    grad *= 2.0 * eta;

    // Find max for numerical stability
    double maxderiv = grad.max();

    // Exponentiated gradient update
    alpha = alpha % arma::exp(grad - maxderiv);

    // Project back to simplex
    alpha /= arma::accu(alpha);

    // Check convergence
    converged = std::abs(std::sqrt(sse_old) - std::sqrt(sse)) < rtol;
    sse_old = sse;
    its++;
  }

  return alpha.t();
}

// Recover Topic-Word Distributions (C++ implementation)
//
// For each word in vocabulary, find optimal convex combination of anchor
// words that reconstructs it using exponentiated gradient descent.
//
// @param Qbar V x V row-normalized gram matrix
// @param anchors K x 1 vector of anchor indices (0-indexed)
// @param p_w V x 1 vector of word probabilities
// @param eta Learning rate for expgrad
// @param maxiter Maximum iterations for expgrad
// @param rtol Convergence tolerance for expgrad
// @param verbose Print progress
// @return K x V matrix of topic-word weights (matching R version)
// [[Rcpp::export]]
arma::mat recoverL2Cpp(const arma::mat& Qbar, const arma::uvec& anchors,
                       const arma::vec& p_w, double eta = 50.0,
                       int maxiter = 500, double rtol = 1e-7,
                       bool verbose = true) {
  int V = Qbar.n_rows;
  int K = anchors.n_elem;

  // Extract anchor rows: X = Qbar[anchors,]
  arma::mat X(K, V);
  for(int k = 0; k < K; k++) {
    X.row(k) = Qbar.row(anchors(k));
  }

  // Precompute X * X^T
  arma::mat XtX = X * X.t();

  // Result matrix: each ROW is mixing weights for a word (V x K, matching R)
  arma::mat A(V, K);

  // Progress tracking
  int progress_interval = V / 10;
  if(progress_interval == 0) progress_interval = 1;

  // For each word
  for(int i = 0; i < V; i++) {
    // Check if this word is an anchor
    bool is_anchor = false;
    int anchor_idx = -1;
    for(int k = 0; k < K; k++) {
      if(anchors(k) == static_cast<arma::uword>(i)) {
        is_anchor = true;
        anchor_idx = k;
        break;
      }
    }

    if(is_anchor) {
      // Anchor words get a one-hot vector
      A.row(i).zeros();
      A(i, anchor_idx) = 1.0;
    } else {
      // Non-anchor words: solve via expgrad
      arma::vec y = Qbar.row(i).t();
      arma::vec solution = expgradCpp(X, y, XtX, eta, maxiter, rtol);
      A.row(i) = solution.t();
    }

    // Print progress
    if(verbose && (i + 1) % progress_interval == 0) {
      Rcpp::Rcout << ".";
    }
  }

  if(verbose) Rcpp::Rcout << std::endl;

  // Apply final transformation to match R version (lines 190-192 in spectral.R):
  // weights <- do.call(rbind, condprob)  # V x K
  // A <- weights*wprob                    # V x K, multiply each row by wprob[i]
  // A <- t(A)/colSums(A)                  # K x V, transpose and normalize rows

  // Multiply each row by corresponding word probability
  for(int i = 0; i < V; i++) {
    A.row(i) *= p_w(i);
  }

  // Compute column sums before transpose (K values)
  arma::rowvec col_sums = arma::sum(A, 0);  // Sum down columns of V×K matrix

  // Transpose to K x V
  A = A.t();

  // Normalize each row by corresponding column sum
  // In R: t(A)/colSums(A) divides each row of t(A) by the K-length colSums vector
  for(int k = 0; k < K; k++) {
    A.row(k) /= col_sums(k);
  }

  // Return K x V matrix with row sums = 1
  return A;
}
