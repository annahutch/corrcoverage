// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
using namespace Rcpp;


// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

//' Obtain pp from a matrix of Zj and ERR
//'
//' @param Zj vector of ...
//' @return pp Matrix of posterior probabilities (one row for each simulation)
//' @export
// [[Rcpp::export]]
NumericMatrix zj_pp_c(const NumericVector & zm,
		      int nrep,
		      const NumericMatrix & ERR,
		      const NumericVector & r){
  NumericMatrix out(nrep,zm.size());
  // process each rep in a loop
  for(int i=0; i < nrep; i++) {
    NumericVector zstar(zm.size());
    for(int j=0; j<zm.size(); j++)
      zstar(j) = zm(j) + ERR(i,j);
    // arma::vec bf(zstar.size());
    double maxbf=-1.0; // stores max log bf, to remove when calculating sum(exp(bf))
    for(int j=0; j<zstar.size(); j++) {
      out(i,j) = (log(1-r(j)) + r(j) * pow(zstar(j),2.0))/2.0;
      if(out(i,j) > maxbf)
	maxbf=out(i,j);
    }
    // sumbf 
    double sumbf=0.0;
    for(int j=0; j<zstar.size(); j++) {
      sumbf += exp(out(i,j)-maxbf);
    }
    sumbf=maxbf + log(sumbf);
    for(int j=0; j<zstar.size(); j++)
      out(i,j) = exp(out(i,j) - sumbf);
  }
  return out;
}

//' Obtain pp from a matrix of Zj and ERR
//'
//' @param Zj vector of ...
//' @return pp Matrix of posterior probabilities (one row for each simulation)
//' @export
// [[Rcpp::export]]
arma::mat zj_pp_arma(const arma::vec & Zj,
		   const arma::mat & sigma,
		   int nrep,
		   const arma::mat & ERR,
		   const arma::vec & r){
  arma::mat out(nrep,Zj.size());
  arma::mat exp_zm = Zj.t() * sigma; // 1 x Zj.size()
  // process each rep in a loop
  for(int i=0; i < nrep; i++) {
    arma::mat zstar = exp_zm + ERR.row(i);
    // arma::vec bf(zstar.size());
    double maxbf=-1.0; // stores max log bf, to remove when calculating sum(exp(bf))
    for(int j=0; j<zstar.size(); j++) {
      out(i,j) = (log(1-r(j)) + r(j) * pow(zstar(j),2.0))/2.0;
      if(out(i,j) > maxbf)
	maxbf=out(i,j);
    }
    // sumbf 
    double sumbf=0.0;
    for(int j=0; j<zstar.size(); j++) {
      sumbf += exp(out(i,j)-maxbf);
    }
    sumbf=maxbf + log(sumbf);
    for(int j=0; j<zstar.size(); j++)
      out(i,j) = exp(out(i,j) - sumbf);
  }
  return out;
}

//' Obtain pp from a matrix of Zj and ERR
//'
//' @param Zj vector of ...
//' @return pp Matrix of posterior probabilities (one row for each simulation)
//' @export
// [[Rcpp::export]]
arma::mat zj_pp_arma2(const arma::vec & zm,
		   int nrep,
		   const arma::mat & ERR,
		   const arma::vec & r){
  arma::mat out(nrep,zm.size());
  // process each rep in a loop
  for(int i=0; i < nrep; i++) {
    arma::mat zstar = zm.t() + ERR.row(i);
    // arma::vec bf(zstar.size());
    double maxbf=-1.0; // stores max log bf, to remove when calculating sum(exp(bf))
    for(int j=0; j<zstar.size(); j++) {
      out(i,j) = (log(1-r(j)) + r(j) * pow(zstar(j),2.0))/2.0;
      if(out(i,j) > maxbf)
	maxbf=out(i,j);
    }
    // sumbf 
    double sumbf=0.0;
    for(int j=0; j<zstar.size(); j++) {
      sumbf += exp(out(i,j)-maxbf);
    }
    sumbf=maxbf + log(sumbf);
    for(int j=0; j<zstar.size(); j++)
      out(i,j) = exp(out(i,j) - sumbf);
  }
  return out;
}
