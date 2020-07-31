#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
double sum_by_col(const arma::sp_mat& x) {
  double result = 0;
  for (size_t i = 0; i < x.n_cols; i++) {
    arma::sp_mat col(x.col(i));
    for (arma::sp_mat::iterator j = col.begin(); j != col.end(); ++j) {
      result += *j * (i + 1);
    }
  }
  return result;
}

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
double sum_by_row(const arma::sp_mat& x) {
  double result = 0;
  for (size_t i = 0; i < x.n_rows; i++) {
    arma::sp_mat row(x.row(i));
    for (arma::sp_mat::iterator j = row.begin(); j != row.end(); ++j) {
      result += *j * (i + 1);
    }
  }
  return result;
}


// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
