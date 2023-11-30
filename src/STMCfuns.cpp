// [[Rcpp::depends(RcppArmadillo)]]

#include "RcppArmadillo.h"

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
