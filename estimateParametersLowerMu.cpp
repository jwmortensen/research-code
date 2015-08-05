#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]] 
List amcmcUpdate(const arma::colvec& draw, const arma::colvec& curMn, const arma::mat& curVar, int curIt) {
  arma::colvec mn;
  arma::mat var;
  if (curIt > 0) {
    mn = ((curIt - 1) * curMn + draw) / curIt;
    if (curIt == 1) {
      var = arma::mat(draw.size(), draw.size(), arma::fill::zeros);
    } else {
      var = (curIt - 2) * curVar + (curIt - 1) * (curMn * curMn.t()) + draw * draw.t();
      var = (var - curIt*(mn * mn.t())) / (curIt - 1);
    }
  } else {
    mn = arma::colvec(curMn.size(), arma::fill::zeros);
    var = arma::mat(draw.size(), draw.size(), arma::fill::zeros);
  }
  return List::create(Named("mn") = mn, Named("var") = var);
}

// [[Rcpp::export]]
arma::colvec CalcLogLambda(const arma::colvec& lambdaStar, 
                           const arma::colvec& E) {
  return (log(E) + lambdaStar) - log(sum(E % exp(lambdaStar)));
}

// [[Rcpp::export]]
double LogLike(const arma::colvec& lambdaStar, 
               const arma::colvec& Nk,
               const arma::colvec& E) {
  arma::colvec logLambdas = CalcLogLambda(lambdaStar, E);
  return sum(Nk % logLambdas);
}

// [[Rcpp::export]]
double LogLambdaPrior(const arma::colvec& lambdaStar, 
                      const arma::colvec& lstarMu, 
                      const double& sig2, 
                      const arma::mat& lambdaInverseMatern) {
  return as_scalar(-0.5 * ((lambdaStar - lstarMu).t() * lambdaInverseMatern * (lambdaStar - lstarMu)) / sig2);
}

// [[Rcpp::export]]
double LogLambdaMuPrior(const arma::colvec& lambdaStar, 
                        const double& sig2, 
                        const arma::mat& lambdaInverseMatern,
  const arma::mat& intercept, const arma::colvec& beta) {
    return as_scalar(-0.5 * ((lambdaStar - intercept * beta).t() * lambdaInverseMatern * (lambdaStar - intercept * beta)) / sig2);
}


// [[Rcpp::export]]
NumericVector mvrnormC(int n, const arma::vec& mu, const arma::mat& sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return wrap(arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma));
}

// [[Rcpp::export]]
NumericVector mvrnormNoChol(int n, const arma::vec& mu, const arma::mat& cholSigma) {
  int ncols = cholSigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return wrap(arma::repmat(mu, 1, n).t() + Y * cholSigma);
}


// [[Rcpp::export]]
void proposeLambda(arma::colvec& propLambda, const arma::colvec& currLambda, const IntegerVector& ptsIndex, const arma::mat& propVar) {
  NumericVector nvLam = wrap(currLambda);
  nvLam[ptsIndex-1] = mvrnormC(1, as<arma::vec>(nvLam[ptsIndex-1]), propVar);
  propLambda = as<arma::colvec>(nvLam);
}


//List MHGibbsC(const int ndraws, const int thinFactor, const List& TempData, const NumericVector& initLambda, const NumericVector& initBeta,
//              const arma::mat& lambdaInverseMatern, const arma::mat& betaInverseMatern, const NumericMatrix& dateLocIndex,
//              const List& closePtsIndex) {
//  // Vars to use throughout
//  int nPredLocs = initLambda.size();
//  int nBlocks = closePtsIndex.size();
//  int nLags = 4;
//  
//  // Initialize storage containers
//  NumericMatrix lambdaStar(ndraws, initLambda.size());
//  NumericMatrix beta(ndraws, initBeta.size());
//  lambdaStar(0, _) = initLambda;
//  beta(0, _) = initBeta;
//  
//  double logMH = 0;
//  
//  // Lambda stuff
//  arma::colvec tempLambda = initLambda;
//  arma::colvec propLambda;
//  double tempLvar = 0.01;
//  double lvarA = 0.01;
//  double lvarB = 0.01;
//  double lamA = 0;
//  double lamB = 0;
//  double lambdaVarConst = 0.01;
//  arma::mat lambdaDiag = arma::mat(nPredLocs / nBlocks, nPredLocs / nBlocks, arma::fill::eye);
//  arma::mat lambdaPropVar = lambdaVarConst * lambdaDiag;
////  int lAMCMCIt = 500;
//  
//  // Beta stuff
//  arma::colvec tempBeta = initBeta;
//  arma::colvec propBeta;
//  double tempBvar = 0.1;
//  double bvarA = 0.01;
//  double bvarB = 0.01;
//  double betaA = 0;
//  double betaB = 0;
//  double betaVarConst = 0.00001;
//  arma::mat betaPropVar = betaVarConst * arma::mat(nLags, nLags, arma::fill::eye);  
//  
//  int n = 0;
//  bool update = false;  
//  for (int i=1; i < ndraws*thinFactor; i++) {
//    update = (i % thinFactor == 0) ? true : false;
//    if (update) n += 1;
//    
//    lamA = lvarA + nPredLocs/2;
//    lamB = lvarB + as_scalar(0.5*tempLambda.t()*lambdaInverseMatern*tempLambda);
//    tempLvar = 1 / R::rgamma(lamA, lamB);
//
//    for (int j=0; j < nBlocks; j++) {
////      if (i > lAMCMCIt) lambdaPropVar = (2.4 * 2.4 / (nPredLocs / nBlocks)) 
//      proposeLambda(propLambda, tempLambda, closePtsIndex[j], lambdaPropVar);
//      logMH = LogLike(TempData, propLambda, tempBeta, dateLocIndex) - LogLike(TempData, tempLambda, tempBeta, dateLocIndex);
//      logMH += LogLambdaPrior(propLambda, tempLvar, lambdaInverseMatern) - LogLambdaPrior(tempLambda, tempLvar, lambdaInverseMatern);
//      
//      tempLambda = (log(runif(1))[0] < logMH) ? propLambda : tempLambda;
//    }
//    if (update) lambdaStar(n, _) = as<NumericVector>(wrap(tempLambda));
//    
//    betaA = bvarA + nLags / 2;
//    betaB = bvarB + as_scalar(0.5*tempBeta.t()*betaInverseMatern*tempBeta);
//    tempBvar = 1/R::rgamma(betaA, betaB);
//    
//    propBeta = as<arma::colvec>(mvrnormC(1, tempBeta, betaPropVar));
//    logMH = LogLike(TempData, tempLambda, propBeta, dateLocIndex) - LogLike(TempData, tempLambda, tempBeta, dateLocIndex);
//    logMH += LogBetaPrior(propBeta, tempBvar, betaInverseMatern) - LogBetaPrior(tempBeta, tempBvar, betaInverseMatern);
//    
//    tempBeta = (log(runif(1))[0] < logMH) ? propBeta : tempBeta;
//    if (update) beta(n, _) = as<NumericVector>(wrap(tempBeta));
//  }
//  return List::create(Named("lambda.star") = lambdaStar, Named("beta") = beta);
//}
