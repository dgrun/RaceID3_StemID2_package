#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
double myVar(double x, NumericVector v){
  return pow(2.0, v[0] + v[1]*log(x)/log(2.0) + v[2]*pow(log(x)/log(2.0),2.0)) ;
}

// [[Rcpp::export]]
double mySize(double x, NumericVector v){
  NumericVector w(2);
  w[0] = x + 1e-6;
  w[1] = myVar(x,v);
  return pow(x,2.0)/(max(w) - x);
}

// [[Rcpp::export]]
NumericVector ProbFun(NumericVector x,NumericVector v,NumericVector w){
  double mu = 0;
  for(int i=0;i<x.size();i++){
    mu += w[i]*x[i];
  }
  double size = mySize(mu,v);
  NumericVector y = round(x,0);
  NumericVector prob = pnbinom_mu( y, size, mu);
  for(int i=0;i<y.size();i++){
    if ( y[i] > mu ) prob[i] = 1 - prob[i];
  }
  return prob;
}

// [[Rcpp::export]]
NumericMatrix applyProb(NumericMatrix x,NumericVector v,NumericVector w){
  NumericMatrix output(x.nrow(),x.ncol());
   for(int i=0;i<x.nrow();i++){
    output(i,_)=ProbFun(x(i,_),v,w);
  } 
  return output;
}

// [[Rcpp::export]]
NumericVector applyNoise(IntegerMatrix x,NumericVector z,NumericVector co,double pv, NumericMatrix pvM){
  NumericVector output(x.ncol());
  for(int i=0;i<x.ncol();i++){
    IntegerVector ind = x(_,i) - 1;
    NumericVector pval = pvM(_,i);
    pval.push_front( 1 );
    ind = ind[pval > pv];

    NumericVector k = z[ind];
    double m = mean(k);
    double v = var(k);
    output(i) = v/myVar(m,co);
  }
  return output;
}
     
// [[Rcpp::export]]
NumericVector applyNoiseReg(IntegerMatrix x,NumericVector z,NumericVector co,double pv, NumericMatrix pvM){
  NumericVector output(x.ncol());
  for(int i=0;i<x.ncol();i++){
    IntegerVector ind = x(_,i) - 1;
    NumericVector pval = pvM(_,i);
    pval.push_front( 1 );
    ind = ind[pval > pv];

    NumericVector k = z[ind];
    double v = var(k);
    output(i) = v;
  }
  return output;
}
