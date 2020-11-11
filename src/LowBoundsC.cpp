

#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


//[[Rcpp::export]]
double alphaArma(const arma::mat& X) {
  double k = X.n_cols;
    double out = k / (k - 1.0) * (1.0 - arma::trace(X) / arma::accu(X));
    return out;
}

//[[Rcpp::export]]
double l2Arma(const arma::mat& X) {
    double k  = X.n_cols,
           s  = 0.0,
           s2 = 0.0;
    for (unsigned int i = 0; i < X.n_cols - 1; i++)
        for (unsigned int j = i + 1; j < X.n_rows; j++)
        {
            s  += 2.0 * X(i, j);
            s2 += 2.0 * X(i, j) * X(i, j);
        }
    double s3 = s + accu(X.diag());
    double out = (s + std::sqrt(k / (k - 1) * s2)) / s3;
    return out;
}

//[[Rcpp::export]]
double l6Arma(const arma::mat& X) {
  // correlation matrix from covariance matrix:
  vec sds = 1/arma::sqrt(X.diag());
  mat Xcor = arma::diagmat(sds) * X * arma::diagmat(sds);
  Xcor.diag().ones();
  mat XCorInv = arma::inv_sympd(Xcor);
  vec smc = 1 - 1 / XCorInv.diag();
  double out = 1 - arma::accu(1 - smc) / arma::accu(Xcor);
  return out;
}


////[[Rcpp::export]]
//double alphaArmaOld(arma::mat X) {
//    double k = X.n_cols;
//    double tr = arma::trace(X);
//    double out = k/(k-1) * (1 - (tr/(arma::accu(X))));
//    return out;
//}

////[[Rcpp::export]]
//double l2ArmaOld(arma::mat X) {
//  double k = X.n_cols;
//  mat X0 = X;
//  X0.diag().zeros();
//  double out = (accu(X0) + sqrt(k/(k-1) * accu(square(X0)))) / accu(X);
//  return out;
//}

////[[Rcpp::export]]
//double l6ArmaOld(arma::mat X) {
//  // correlation matrix from covariance matrix:
//  vec sds = 1/sqrt(X.diag());
//  mat Xcor = diagmat(sds) * X * diagmat(sds);
//  Xcor.diag().ones();
//  mat XCorInv = inv_sympd(Xcor);
//  vec smc = 1 - 1 / XCorInv.diag();
//  double out = 1 - accu(1 - smc) / accu(Xcor);
//  return out;
//}

