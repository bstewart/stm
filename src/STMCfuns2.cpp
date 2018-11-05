// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>


using namespace Eigen;

// [[Rcpp::export]]
double lhoodcpp2(Eigen::ArrayXd eta, Eigen::ArrayXXd beta, Eigen::ArrayXd doc_ct, Eigen::ArrayXd mu, Eigen::ArrayXXd siginv) {
  int nEta = eta.size();
  ArrayXd expeta = ArrayXd(nEta + 1);
  expeta << eta.exp(), 1;
  
  double doc_ct_sum = doc_ct.sum();
  double part1 = (((expeta.matrix().transpose() * beta.matrix()).array().log()).matrix() * doc_ct.matrix())(0, 0) - doc_ct_sum * log(expeta.sum());
  MatrixXd diff = (eta - mu).matrix();
  double part2 = 0.5 * (diff.transpose() * siginv.matrix() * diff)(0, 0);
  
  return part2-part1;
}

// [[Rcpp::export]]
Eigen::ArrayXXd gradcpp2(Eigen::ArrayXd eta, Eigen::ArrayXXd beta, Eigen::ArrayXd doc_ct, Eigen::ArrayXd mu, Eigen::ArrayXXd siginv) {
  int nEta = eta.size();
  ArrayXd expeta = ArrayXd(nEta + 1);
  expeta << eta.exp(), 1;
  
  double doc_ct_sum = doc_ct.sum();
  double expeta_sum = expeta.sum();
  
  beta.colwise() *= expeta;
  ArrayXd part1 = (beta.matrix() * (doc_ct / beta.colwise().sum().transpose()).matrix()).array() - (doc_ct_sum/expeta_sum)*expeta;
  ArrayXd part2 = siginv.matrix() * (eta - mu).matrix();
  part1 = part1.head(nEta).eval();
  
  return part2-part1;
}

// [[Rcpp::export]]
Rcpp::List hpbcpp2(Eigen::ArrayXd eta, Eigen::ArrayXXd beta, Eigen::ArrayXd doc_ct, Eigen::ArrayXd mu, Eigen::ArrayXXd siginv, double sigmaentropy) {
  
  int nEta = eta.size();
  ArrayXd expeta = ArrayXd(nEta + 1);
  expeta << eta.exp(), 1;
  ArrayXd theta = expeta/expeta.sum();
  
  //create a new version of the matrix so we can mess with it
  //multiply each column by expeta
  //this should be fastest as its column-major ordering
  ArrayXXd EB = beta.colwise() * expeta;
  
  //divide out by the column sums
  EB.rowwise() *= doc_ct.sqrt().transpose() / EB.colwise().sum();
  
  double doc_ct_sum = doc_ct.sum();
  //Combine the pieces of the Hessian which are matrices
  MatrixXd hess = (EB.matrix() * EB.matrix().transpose()) - doc_ct_sum*(theta.matrix() * theta.matrix().transpose());
  
  //we don't need EB any more so we turn it into phi
  EB.rowwise() *= doc_ct.sqrt().transpose();
  
  //Now alter just the diagonal of the Hessian
  ArrayXd xx = EB.rowwise().sum() - doc_ct_sum*theta;
  MatrixXd yy = xx.matrix().asDiagonal();
  hess.array() -= yy.array();
  
  //Drop the last row and column - notice the eval() here to avoid aliasing issues!
  hess = hess.topLeftCorner(nEta, nEta).eval();
  //Now we can add in siginv
  hess.array() += siginv;
  //At this point the Hessian is complete.
  
  ////
  //Invert via cholesky decomposition
  ////
  
  MatrixXd nu(nEta, nEta);
  LLT<MatrixXd, Upper> chol(hess);
  
  if (chol.info()!=Success) {
    //It failed!  Oh Nos.
    // So the matrix wasn't positive definite.  In practice this means that it hasn't
    // converged probably along some minor aspect of the dimension.
    
    //Here we make it positive definite through diagonal dominance
    ArrayXd diag = hess.matrix().diagonal();
    ArrayXd magnitudes = hess.array().abs().rowwise().sum() - diag.abs();
    hess.matrix().diagonal() = diag.max(magnitudes);
    chol = LLT<MatrixXd, Upper>(hess);
  }
  nu = chol.matrixU();
  
  //compute 1/2 the determinant from the cholesky decomposition
  double detTerm = -nu.diagonal().array().log().sum();
  // Though verbose, this is slightly faster than doing:
  // double detTerm = -log(nu.determinant());
  
  //Now finish constructing nu
  // TODO: There's got to be a way to get the inverse directly from chol!
  nu = chol.solve(MatrixXd::Identity(nEta, nEta));
  
  //Precompute the difference since we use it twice
  MatrixXd diff = (eta - mu).matrix();
  
  double t = ((theta.matrix().transpose() * beta.matrix()).array().log().matrix() * doc_ct.matrix())(0, 0);
  float distance = 0.5 * (diff.transpose() * siginv.matrix() * diff)(0, 0);
  double bound = t + detTerm - distance - sigmaentropy;
  
  return Rcpp::List::create(
    Rcpp::Named("phis") = EB,
    Rcpp::Named("eta") = Rcpp::List::create(
      Rcpp::Named("lambda") = MatrixXd(eta),
      Rcpp::Named("nu") = nu
    ),
    Rcpp::Named("bound")=bound
  );
}
