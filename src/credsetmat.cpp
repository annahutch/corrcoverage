#include <Rcpp.h>
using namespace Rcpp;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

//' Obtain credible sets from a matrix of posterior probabilities
//'
//' @param pp Matrix of posterior probabilities (one row for each simulation)
//' @param iCV A vector of the indices of the CV
//' @param threshold The threshold to use to generate the credible set
//' @export
// [[Rcpp::export]]
List credsetmat(NumericMatrix pp, NumericVector iCV, double threshold) {
  int n=pp.nrow();
  int m=pp.ncol();
  NumericVector nvar(n); // number of variants in set
  NumericVector setsize(n); // size
  NumericVector covered(n); // start: iCV not in empty set
  // std::vector<int> V(n);
  for(int rep=0; rep<n; rep++) {
    int iiCV=iCV(rep)-1; // 0-based
    double csum=0.0;
    int i=0;
    NumericVector ppv = pp( rep, _);
    IntegerVector V=seq(0, m-1);
    // std::iota(V.begin(),V.end(),i++); //Initializing
    std::sort( V.begin(),V.end(), [&](int i,int j){
      return ppv(i)>ppv(j);
    }
    );
    // return V;
    for(i=0; i<m; i++) {
      csum+= ppv(V(i));
      if(iiCV==V(i))
        covered(rep)=1;
      if (csum > threshold) {
        nvar(rep)=i+1;
        setsize(rep)=csum;
        break;
      }
    }
  }
  return List::create(
    _["nvar"] = nvar,
    _["claimed"] = setsize,
    _["covered"] = covered
  );
}
