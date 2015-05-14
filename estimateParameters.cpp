#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
NumericMatrix CalcLogLambda(arma::mat H, arma::colvec lambdaStar, arma::colvec beta) {
  return wrap(lambdaStar + H * beta - log(sum(exp(lambdaStar + H*beta))));
}

// [[Rcpp::export]]
double LogLike(List H, arma::colvec lambdaStar, arma::colvec beta, NumericMatrix dateLocIndex) {
  NumericMatrix logLambdas(lambdaStar.size(), H.size());
  for (int i = 0; i < H.size(); i++) {
    logLambdas(_, i) = CalcLogLambda(as<arma::mat>(H[i]), lambdaStar, beta);
  }
  double out = 0;
  for (int i = 0; i < dateLocIndex.nrow(); i++) {
    out += logLambdas(dateLocIndex(i, 1) - 1, dateLocIndex(i, 0) - 1);
  }
  return out;
}

// [[Rcpp::export]]
double LogLambdaPrior(arma::colvec lambdaStar, double sig2, arma::mat lambdaInverseMatern) {
  return as_scalar(-0.5*(lambdaStar.t() * lambdaInverseMatern * lambdaStar) / sig2);
}

// [[Rcpp::export]]
double LogBetaPrior(arma::colvec beta, double sig2, arma::mat betaInverseMatern) {
  return as_scalar(-0.5*(beta.t() * betaInverseMatern * beta) / sig2);
}

// [[Rcpp::export]]
NumericVector mvrnormC(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return wrap(arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma));
}

// [[Rcpp::export]]
arma::colvec proposeLambda(arma::colvec currLambda, IntegerVector ptsIndex, arma::mat propVar) {
  NumericVector nvLam = wrap(currLambda);
  nvLam[ptsIndex-1] = mvrnormC(1, as<arma::vec>(nvLam[ptsIndex-1]), propVar);
  return as<arma::colvec>(nvLam);
}

// [[Rcpp::export]]
List MHGibbsC(const int ndraws, const int thinFactor, const List TempData, const NumericVector initLambda, const NumericVector initBeta,
              const arma::mat lambdaInverseMatern, const arma::mat betaInverseMatern, const NumericMatrix dateLocIndex,
              const List closePtsIndex) {
  // Vars to use throughout
  int nPredLocs = initLambda.size();
  int nBlocks = closePtsIndex.size();
  int nLags = 4;
  
  // Initialize storage containers
  NumericMatrix lambdaStar(ndraws, initLambda.size());
  NumericMatrix beta(ndraws, initBeta.size());
  lambdaStar(0, _) = initLambda;
  beta(0, _) = initBeta;
  
  double logMH = 0;
  
  arma::colvec tempLambda = initLambda;
  double tempLvar = 0.01;
  double lvarA = 0.01;
  double lvarB = 0.01;
  double lamA = 0;
  double lamB = 0;
  double lambdaVarConst = 0.00000001;
  NumericMatrix lambdaPropTemp = NumericMatrix(nPredLocs/nBlocks, nPredLocs/nBlocks);
  lambdaPropTemp.fill_diag(lambdaVarConst);
  arma::mat lambdaPropVar = as<arma::mat>(lambdaPropTemp);
  
  arma::colvec tempBeta = initBeta;
  double tempBvar = 0.1;
  double bvarA = 0.01;
  double bvarB = 0.01;
  double betaA = 0;
  double betaB = 0;
  double betaVarConst = 0.000005;
  NumericMatrix betaPropTemp = NumericMatrix(nLags, nLags);
  betaPropTemp.fill_diag(betaVarConst);
  arma::mat betaPropVar = as<arma::mat>(betaPropTemp);  
  
  
  int n = 0;
  bool update = false;  
  for (int i=1; i < ndraws*thinFactor; i++) {
    update = (i % thinFactor == 0) ? true : false;
    if (update) n += 1;
    
    lamA = lvarA + nPredLocs/2;
    lamB = lvarB + as_scalar(0.5*tempLambda.t()*lambdaInverseMatern*tempLambda);
    tempLvar = 1 / R::rgamma(lamA, lamB);

    for (int j=0; j < nBlocks; j++) {
      arma::colvec propLambda = proposeLambda(tempLambda, closePtsIndex[j], lambdaPropVar);
      logMH = LogLike(TempData, propLambda, tempBeta, dateLocIndex) - LogLike(TempData, tempLambda, tempBeta, dateLocIndex);
      logMH += LogLambdaPrior(propLambda, tempLvar, lambdaInverseMatern) - LogLambdaPrior(tempLambda, tempLvar, lambdaInverseMatern);
      
      tempLambda = (log(runif(1))[0] < logMH) ? propLambda : tempLambda;
    }
    if (update) lambdaStar(n, _) = as<NumericVector>(wrap(tempLambda));
    
    betaA = bvarA + nLags / 2;
    betaB = bvarB + as_scalar(0.5*tempBeta.t()*betaInverseMatern*tempBeta);
    tempBvar = 1/R::rgamma(betaA, betaB);
    
    arma::colvec propBeta = as<arma::colvec>(mvrnormC(1, tempBeta, betaPropVar));
    logMH = LogLike(TempData, tempLambda, propBeta, dateLocIndex) - LogLike(TempData, tempLambda, tempBeta, dateLocIndex);
    logMH += LogBetaPrior(propBeta, tempBvar, betaInverseMatern) - LogBetaPrior(tempBeta, tempBvar, betaInverseMatern);
    
    tempBeta = (log(runif(1))[0] < logMH) ? propBeta : tempBeta;
    if (update) beta(n, _) = as<NumericVector>(wrap(tempBeta));
  }
  return List::create(Named("lambda.star") = lambdaStar, Named("beta") = beta);
}
