// [[Rcpp::depends(RcppArmadillo)]]

#include "RcppArmadillo.h"

#include <cfloat>
#include <cmath>

namespace {

class DocumentObjective {
 public:
  DocumentObjective(const arma::mat& beta,
                    const arma::vec& doc_ct,
                    const arma::vec& mu,
                    const arma::mat& siginv)
      : beta_(beta),
        doc_ct_(doc_ct),
        mu_(mu),
        siginv_(siginv),
        ndoc_(arma::sum(doc_ct)),
        expeta_(mu.n_elem + 1),
        mixture_(doc_ct.n_elem),
        diff_(mu.n_elem),
        prior_gradient_(mu.n_elem) {}

  double evaluate(const arma::vec& eta) {
    const arma::uword neta = eta.n_elem;
    expeta_.head(neta) = arma::exp(eta);
    expeta_(neta) = 1.0;
    mixture_ = beta_.t() * expeta_;
    expeta_sum_ = arma::sum(expeta_);
    diff_ = eta - mu_;
    prior_gradient_ = siginv_ * diff_;

    if (!expeta_.is_finite() || !mixture_.is_finite() ||
        !std::isfinite(expeta_sum_) || expeta_sum_ <= 0.0 ||
        arma::any(mixture_ <= 0.0)) {
      value_ = R_PosInf;
      return value_;
    }

    value_ = 0.5 * arma::dot(diff_, prior_gradient_) -
             arma::dot(arma::log(mixture_), doc_ct_) +
             ndoc_ * std::log(expeta_sum_);
    if (!std::isfinite(value_)) value_ = R_PosInf;
    return value_;
  }

  arma::vec gradient() const {
    arma::vec likelihood_gradient =
        expeta_ % (beta_ * (doc_ct_ / mixture_)) -
        (ndoc_ / expeta_sum_) * expeta_;
    return prior_gradient_ - likelihood_gradient.head(mu_.n_elem);
  }

 private:
  const arma::mat& beta_;
  const arma::vec& doc_ct_;
  const arma::vec& mu_;
  const arma::mat& siginv_;
  const double ndoc_;
  arma::vec expeta_;
  arma::vec mixture_;
  arma::vec diff_;
  arma::vec prior_gradient_;
  double expeta_sum_;
  double value_;
};

}  // namespace

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double lhoodcpp(SEXP eta,
                   SEXP beta,
                   SEXP doc_ct,
                   SEXP mu,
                   SEXP siginv){
   
   Rcpp::NumericVector etav(eta); 
   arma::vec etas(etav.begin(), etav.size(), false);
   Rcpp::NumericMatrix betam(beta);
   arma::mat betas(betam.begin(), betam.nrow(), betam.ncol(), false);
   Rcpp::NumericVector doc_ctv(doc_ct);
   arma::vec doc_cts(doc_ctv.begin(), doc_ctv.size(), false);
   Rcpp::NumericVector muv(mu);
   arma::vec mus(muv.begin(), muv.size(), false);
   Rcpp::NumericMatrix siginvm(siginv);
   arma::mat siginvs(siginvm.begin(), siginvm.nrow(), siginvm.ncol(), false);
   
   arma::rowvec expeta(etas.size()+1); 
   expeta.fill(1);
   int neta = etav.size(); 
   for(int j=0; j <neta;  j++){
     expeta(j) = exp(etas(j));
   }
   double ndoc = sum(doc_cts);
   double part1 = arma::as_scalar(log(expeta*betas)*doc_cts - ndoc*log(sum(expeta)));
   arma::vec diff = etas - mus;
   double part2 = .5*arma::as_scalar(diff.t()*siginvs*diff);
   double out = part2 - part1;
   return out;
}

// [[Rcpp::export]]
arma::vec gradcpp(SEXP eta,
                   SEXP beta,
                   SEXP doc_ct,
                   SEXP mu,
                   SEXP siginv){
   
   Rcpp::NumericVector etav(eta); 
   arma::vec etas(etav.begin(), etav.size(), false);
   Rcpp::NumericMatrix betam(beta);
   arma::mat betas(betam.begin(), betam.nrow(), betam.ncol());
   Rcpp::NumericVector doc_ctv(doc_ct);
   arma::vec doc_cts(doc_ctv.begin(), doc_ctv.size(), false);
   Rcpp::NumericVector muv(mu);
   arma::vec mus(muv.begin(), muv.size(), false);
   Rcpp::NumericMatrix siginvm(siginv);
   arma::mat siginvs(siginvm.begin(), siginvm.nrow(), siginvm.ncol(), false);
   
    arma::colvec expeta(etas.size()+1); 
    expeta.fill(1);
    int neta = etas.size(); 
    for(int j=0; j <neta;  j++){
       expeta(j) = exp(etas(j));
    }
    betas.each_col() %= expeta;
    arma::vec part1 = betas*(doc_cts/arma::trans(sum(betas,0))) - (sum(doc_cts)/sum(expeta))*expeta;
    arma::vec part2 = siginvs*(etas - mus);
    part1.shed_row(neta);
    return part2-part1;
}

// Algebraically equivalent gradient that avoids copying and rescaling the
// complete K by U document-specific beta matrix on every evaluation.
// [[Rcpp::export]]
arma::vec gradcpp_nocopy(SEXP eta,
                         SEXP beta,
                         SEXP doc_ct,
                         SEXP mu,
                         SEXP siginv) {
  Rcpp::NumericVector etav(eta);
  arma::vec etas(etav.begin(), etav.size(), false);
  Rcpp::NumericMatrix betam(beta);
  arma::mat betas(betam.begin(), betam.nrow(), betam.ncol(), false);
  Rcpp::NumericVector doc_ctv(doc_ct);
  arma::vec doc_cts(doc_ctv.begin(), doc_ctv.size(), false);
  Rcpp::NumericVector muv(mu);
  arma::vec mus(muv.begin(), muv.size(), false);
  Rcpp::NumericMatrix siginvm(siginv);
  arma::mat siginvs(siginvm.begin(), siginvm.nrow(), siginvm.ncol(), false);

  DocumentObjective objective(betas, doc_cts, mus, siginvs);
  objective.evaluate(etas);
  return objective.gradient();
}

// Compiled BFGS specialized for the document objective. The last objective
// evaluation in each accepted line search is cached and reused by the
// gradient, avoiding duplicate exponentiation and topic-word mixtures.
// [[Rcpp::export]]
Rcpp::List optimizecpp(SEXP eta,
                       SEXP beta,
                       SEXP doc_ct,
                       SEXP mu,
                       SEXP siginv,
                       Rcpp::List control) {
  Rcpp::NumericVector etav(eta);
  arma::vec current(etav.begin(), etav.size(), true);
  Rcpp::NumericMatrix betam(beta);
  arma::mat betas(betam.begin(), betam.nrow(), betam.ncol(), false);
  Rcpp::NumericVector doc_ctv(doc_ct);
  arma::vec doc_cts(doc_ctv.begin(), doc_ctv.size(), false);
  Rcpp::NumericVector muv(mu);
  arma::vec mus(muv.begin(), muv.size(), false);
  Rcpp::NumericMatrix siginvm(siginv);
  arma::mat siginvs(siginvm.begin(), siginvm.nrow(), siginvm.ncol(), false);

  const int maxit = control.containsElementNamed("maxit")
                        ? Rcpp::as<int>(control["maxit"])
                        : 100;
  const double reltol = control.containsElementNamed("reltol")
                            ? Rcpp::as<double>(control["reltol"])
                            : std::sqrt(DBL_EPSILON);
  const double abstol = control.containsElementNamed("abstol")
                            ? Rcpp::as<double>(control["abstol"])
                            : R_NegInf;

  DocumentObjective objective(betas, doc_cts, mus, siginvs);
  double value = objective.evaluate(current);
  if (!std::isfinite(value)) {
    Rcpp::stop("initial value in document BFGS is not finite");
  }
  if (maxit <= 0) {
    return Rcpp::List::create(
        Rcpp::Named("par") = Rcpp::NumericVector(current.begin(),
                                                 current.end()),
        Rcpp::Named("value") = value,
        Rcpp::Named("counts") = Rcpp::NumericVector::create(0, 0),
        Rcpp::Named("convergence") = 0);
  }

  arma::vec gradient = objective.gradient();
  int function_evaluations = 1;
  int gradient_evaluations = 1;
  int iteration = 1;
  int last_restart = gradient_evaluations;
  int unchanged_coordinates = 0;
  arma::mat inverse_hessian(current.n_elem, current.n_elem,
                            arma::fill::zeros);
  arma::vec base(current.n_elem);
  arma::vec old_gradient(current.n_elem);
  arma::vec direction(current.n_elem);
  arma::vec candidate(current.n_elem);
  const double armijo = 1e-4;
  const double step_reduction = 0.2;
  const double relative_step_test = 10.0;

  do {
    if ((iteration & 15) == 0) Rcpp::checkUserInterrupt();
    if (last_restart == gradient_evaluations) inverse_hessian.eye();
    base = current;
    old_gradient = gradient;

    double directional_derivative = 0.0;
    for (arma::uword row = 0; row < current.n_elem; ++row) {
      double component = 0.0;
      for (arma::uword column = 0; column <= row; ++column) {
        component -= inverse_hessian(row, column) * gradient(column);
      }
      for (arma::uword column = row + 1; column < current.n_elem; ++column) {
        component -= inverse_hessian(column, row) * gradient(column);
      }
      direction(row) = component;
      directional_derivative += component * gradient(row);
    }

    if (directional_derivative < 0.0) {
      double step = 1.0;
      double candidate_value = value;
      bool accepted = false;
      do {
        unchanged_coordinates = 0;
        for (arma::uword coordinate = 0;
             coordinate < current.n_elem; ++coordinate) {
          candidate(coordinate) =
              base(coordinate) + step * direction(coordinate);
          if (relative_step_test + base(coordinate) ==
              relative_step_test + candidate(coordinate)) {
            ++unchanged_coordinates;
          }
        }
        current = candidate;
        if (unchanged_coordinates < static_cast<int>(current.n_elem)) {
          candidate_value = objective.evaluate(current);
          ++function_evaluations;
          accepted = std::isfinite(candidate_value) &&
                     candidate_value <=
                         value + directional_derivative * step * armijo;
          if (!accepted) step *= step_reduction;
        }
      } while (unchanged_coordinates < static_cast<int>(current.n_elem) &&
               !accepted);

      const bool enough =
          candidate_value > abstol &&
          std::abs(candidate_value - value) >
              reltol * (std::abs(value) + reltol);
      if (!enough) {
        unchanged_coordinates = static_cast<int>(current.n_elem);
        value = candidate_value;
      }

      if (unchanged_coordinates < static_cast<int>(current.n_elem)) {
        value = candidate_value;
        gradient = objective.gradient();
        ++gradient_evaluations;
        ++iteration;

        double curvature = 0.0;
        for (arma::uword coordinate = 0;
             coordinate < current.n_elem; ++coordinate) {
          direction(coordinate) *= step;
          old_gradient(coordinate) =
              gradient(coordinate) - old_gradient(coordinate);
          curvature += direction(coordinate) * old_gradient(coordinate);
        }

        if (curvature > 0.0) {
          arma::vec hessian_gradient_change(current.n_elem);
          double quadratic = 0.0;
          for (arma::uword row = 0; row < current.n_elem; ++row) {
            double component = 0.0;
            for (arma::uword column = 0; column <= row; ++column) {
              component += inverse_hessian(row, column) *
                           old_gradient(column);
            }
            for (arma::uword column = row + 1;
                 column < current.n_elem; ++column) {
              component += inverse_hessian(column, row) *
                           old_gradient(column);
            }
            hessian_gradient_change(row) = component;
            quadratic += component * old_gradient(row);
          }
          const double update_scale = 1.0 + quadratic / curvature;
          for (arma::uword row = 0; row < current.n_elem; ++row) {
            for (arma::uword column = 0; column <= row; ++column) {
              inverse_hessian(row, column) +=
                  (update_scale * direction(row) * direction(column) -
                   hessian_gradient_change(row) * direction(column) -
                   direction(row) * hessian_gradient_change(column)) /
                  curvature;
            }
          }
        } else {
          last_restart = gradient_evaluations;
        }
      } else if (last_restart < gradient_evaluations) {
        unchanged_coordinates = 0;
        last_restart = gradient_evaluations;
      }
    } else {
      unchanged_coordinates = 0;
      if (last_restart == gradient_evaluations) {
        unchanged_coordinates = static_cast<int>(current.n_elem);
      } else {
        last_restart = gradient_evaluations;
      }
    }

    if (iteration >= maxit) break;
    if (gradient_evaluations - last_restart >
        2 * static_cast<int>(current.n_elem)) {
      last_restart = gradient_evaluations;
    }
  } while (unchanged_coordinates != static_cast<int>(current.n_elem) ||
           last_restart != gradient_evaluations);

  const int convergence = iteration < maxit ? 0 : 1;

  return Rcpp::List::create(
      Rcpp::Named("par") = Rcpp::NumericVector(current.begin(), current.end()),
      Rcpp::Named("value") = value,
      Rcpp::Named("counts") =
          Rcpp::NumericVector::create(function_evaluations,
                                      gradient_evaluations),
      Rcpp::Named("convergence") = convergence);
}

// [[Rcpp::export]]
SEXP hpbcpp(SEXP eta,
            SEXP beta,
            SEXP doc_ct,
            SEXP mu,
            SEXP siginv,
            SEXP sigmaentropy){
 
   Rcpp::NumericVector etav(eta); 
   arma::vec etas(etav.begin(), etav.size(), false);
   Rcpp::NumericMatrix betam(beta);
   arma::mat betas(betam.begin(), betam.nrow(), betam.ncol());
   Rcpp::NumericVector doc_ctv(doc_ct);
   arma::vec doc_cts(doc_ctv.begin(), doc_ctv.size(), false);
   Rcpp::NumericVector muv(mu);
   arma::vec mus(muv.begin(), muv.size(), false);
   Rcpp::NumericMatrix siginvm(siginv);
   arma::mat siginvs(siginvm.begin(), siginvm.nrow(), siginvm.ncol(), false);
   Rcpp::NumericVector sigmaentropym(sigmaentropy);
   arma::vec entropy(sigmaentropym);

   //Performance Nots from 3/6/2015
   //  I tried a few different variants and benchmarked this one as roughly twice as
   //  fast as the R code for a K=100 problem.  Key to performance was not creating
   //  too many objects and being selective in how things were flagged as triangular.
   //  Some additional notes in the code below.
   //
   //  Some things this doesn't have or I haven't tried
   //  - I didn't tweak the arguments much.  sigmaentropy is a double, and I'm still
   //    passing beta in the same way.  I tried doing a ", false" for beta but it didn't
   //    change much so I left it the same as in gradient.  
   //  - I tried treating the factors for doc_cts and colSums(EB) as a diagonal matrix- much slower.
   
   //  Haven't Tried/Done
   //  - each_row() might be much slower (not sure but arma is column order).  Maybe transpose in place?
   //  - depending on costs there are some really minor calculations that could be precomputed: 
   //     - sum(doc_ct)
   //     - sqrt(doc_ct)
   
   //  More on passing by reference here:
   //  - Hypothetically we could alter beta (because hessian is last thing we do) however down
   //    the road we may want to explore treating nonPD hessians by optimization at which point
   //    we would need it again.
   
   arma::colvec expeta(etas.size()+1); 
   expeta.fill(1);
   int neta = etas.size(); 
   for(int j=0; j <neta;  j++){
     expeta(j) = exp(etas(j));
   }
   arma::vec theta = expeta/sum(expeta);

   //create a new version of the matrix so we can mess with it
   arma::mat EB(betam.begin(), betam.nrow(), betam.ncol());
   //multiply each column by expeta
   EB.each_col() %= expeta; //this should be fastest as its column-major ordering
  
   //divide out by the column sums
   EB.each_row() %= arma::trans(sqrt(doc_cts))/sum(EB,0);
    
   //Combine the pieces of the Hessian which are matrices
   arma::mat hess = EB*EB.t() - sum(doc_cts)*(theta*theta.t());
  
   //we don't need EB any more so we turn it into phi
   EB.each_row() %= arma::trans(sqrt(doc_cts));
   
   //Now alter just the diagonal of the Hessian
   hess.diag() -= sum(EB,1) - sum(doc_cts)*theta;
   //Drop the last row and column
   hess.shed_row(neta);
   hess.shed_col(neta);
   //Now we can add in siginv
   hess = hess + siginvs;
   //At this point the Hessian is complete.
   
   //This next bit of code is from http://arma.sourceforge.net/docs.html#logging
   //It basically keeps arma from printing errors from chol to the console.
   std::ostream nullstream(0);
   //arma::set_stream_err2(nullstream);
   //arma::arma_cerr_stream<char>(&nullstream);
   
   ////
   //Invert via cholesky decomposition
   ////
   //Start by initializing an object
   arma::mat nu = arma::mat(hess.n_rows, hess.n_rows);
   //This version of chol generates a boolean which tells us if it failed.
   bool worked = arma::chol(nu,hess);
   if(!worked) {
     //It failed!  Oh Nos.
     // So the matrix wasn't positive definite.  In practice this means that it hasn't
     // converged probably along some minor aspect of the dimension.
     
     //Here we make it positive definite through diagonal dominance
     arma::vec dvec = hess.diag();
     //find the magnitude of the diagonal 
     arma::vec magnitudes = sum(abs(hess), 1) - abs(dvec);
     //iterate over each row and set the minimum value of the diagonal to be the magnitude of the other terms
     int Km1 = dvec.size();
     for(int j=0; j < Km1;  j++){
       if(arma::as_scalar(dvec(j)) < arma::as_scalar(magnitudes(j))) dvec(j) = magnitudes(j); //enforce diagonal dominance 
     }
     //overwrite the diagonal of the hessian with our new object
     hess.diag() = dvec;
     //that was sufficient to ensure positive definiteness so we now do cholesky
     nu = arma::chol(hess);
   }
   //compute 1/2 the determinant from the cholesky decomposition
   double detTerm = -sum(log(nu.diag()));
   
   //Now finish constructing nu
   nu = arma::inv(arma::trimatu(nu));
   nu = nu * nu.t(); //trimatu doesn't do anything for multiplication so it would just be timesink to signal here.
   
   //Precompute the difference since we use it twice
   arma::vec diff = etas - mus;
   //Now generate the bound and make it a scalar
   double bound = arma::as_scalar(log(arma::trans(theta)*betas)*doc_cts + detTerm - .5*diff.t()*siginvs*diff - entropy); 
   
   // Generate a return list that mimics the R output
   return Rcpp::List::create(
        Rcpp::Named("phis") = EB,
        Rcpp::Named("eta") = Rcpp::List::create(Rcpp::Named("lambda")=etas, Rcpp::Named("nu")=nu),
        Rcpp::Named("bound") = bound
        );
}

// Posterior calculation with a zero-copy view of beta. hpbcpp() first copies
// beta and then copies it again into EB; this version retains only the EB copy
// that is mutated into the assignment matrix.
// [[Rcpp::export]]
SEXP hpbcpp_nocopy(SEXP eta,
                   SEXP beta,
                   SEXP doc_ct,
                   SEXP mu,
                   SEXP siginv,
                   SEXP sigmaentropy) {
  Rcpp::NumericVector etav(eta);
  arma::vec etas(etav.begin(), etav.size(), false);
  Rcpp::NumericMatrix betam(beta);
  arma::mat betas(betam.begin(), betam.nrow(), betam.ncol(), false);
  Rcpp::NumericVector doc_ctv(doc_ct);
  arma::vec doc_cts(doc_ctv.begin(), doc_ctv.size(), false);
  Rcpp::NumericVector muv(mu);
  arma::vec mus(muv.begin(), muv.size(), false);
  Rcpp::NumericMatrix siginvm(siginv);
  arma::mat siginvs(siginvm.begin(), siginvm.nrow(), siginvm.ncol(), false);
  const double entropy = Rcpp::as<double>(sigmaentropy);

  const arma::uword neta = etas.n_elem;
  arma::vec expeta(neta + 1, arma::fill::ones);
  expeta.head(neta) = arma::exp(etas);
  const double expeta_sum = arma::sum(expeta);
  const arma::vec theta = expeta / expeta_sum;

  arma::mat EB = betas;
  EB.each_col() %= expeta;
  const arma::rowvec mixture = arma::sum(EB, 0);
  const arma::rowvec sqrt_counts = arma::trans(arma::sqrt(doc_cts));
  EB.each_row() %= sqrt_counts / mixture;

  arma::mat hess = EB * EB.t() -
                   arma::sum(doc_cts) * (theta * theta.t());
  EB.each_row() %= sqrt_counts;
  hess.diag() -= arma::sum(EB, 1) - arma::sum(doc_cts) * theta;
  hess.shed_row(neta);
  hess.shed_col(neta);
  hess += siginvs;

  arma::mat nu(hess.n_rows, hess.n_cols);
  bool worked = arma::chol(nu, hess);
  if (!worked) {
    arma::vec diagonal = hess.diag();
    const arma::vec magnitudes = arma::sum(arma::abs(hess), 1) -
                                 arma::abs(diagonal);
    for (arma::uword j = 0; j < diagonal.n_elem; ++j) {
      if (diagonal(j) < magnitudes(j)) diagonal(j) = magnitudes(j);
    }
    hess.diag() = diagonal;
    nu = arma::chol(hess);
  }

  const double determinant_term = -arma::sum(arma::log(nu.diag()));
  nu = arma::inv(arma::trimatu(nu));
  nu = nu * nu.t();

  const arma::vec diff = etas - mus;
  const double bound =
      arma::dot(arma::log(mixture / expeta_sum), doc_cts) +
      determinant_term -
      0.5 * arma::dot(diff, siginvs * diff) - entropy;

  return Rcpp::List::create(
      Rcpp::Named("phis") = EB,
      Rcpp::Named("eta") = Rcpp::List::create(
          Rcpp::Named("lambda") = etas,
          Rcpp::Named("nu") = nu),
      Rcpp::Named("bound") = bound);
}
