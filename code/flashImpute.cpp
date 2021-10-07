#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List SOLQL(
    arma::mat Y, // n by p
    arma::mat Af, // p by K
    arma::mat Cf, // p by K
    arma::umat Omega, // n by p
    arma::vec taul, // K
    double tau
  ){
  int n = size(Y)[0];
  //int p = size(Y)[1];
  int K = size(Af)[1];

  arma::mat Al(n,K); // n by K
  arma::mat Cl(n,K); // n by K

  for (int idxit = 0; idxit < n; idxit++){

    arma::uvec Omegait = Omega.row(idxit).t();
    arma::uvec idxobs = find(Omegait ==1);
    arma::mat Afi = Af.rows(idxobs);
    arma::mat Cfi = Cf.rows(idxobs);
    arma::mat Bfi = Afi.t()*Afi;
    Bfi.diag() = sum(Cfi,0);

    arma::mat Yi = Y.row(idxit);
    Yi = Yi.cols(idxobs);

    arma::mat Ali = Yi * Afi * inv(Bfi + diagmat(taul/tau));
    arma::mat Cli = pow(Ali, 2.0) + 1/(taul + tau * Bfi.diag()).t();

    Al.row(idxit) = Ali;
    Cl.row(idxit) = Cli;
  }

  return Rcpp::List::create(Rcpp::Named("A.l") = Al,
                            Rcpp::Named("C.l") = Cl
                              );
}


// [[Rcpp::export]]
double SOLTAU(
    arma::mat Y, // n by p
    arma::umat Omega, // n by p
    arma::mat Al, // n by K
    arma::mat Af, // p by K
    arma::mat Cl, // n by K
    arma::mat Cf // p by K
){
  arma::umat Omegaidx = find(Omega==1);

  double temp1= mean(square(arma::mat(Y - Al*Af.t())(Omegaidx)));
  double temp2 = mean(arma::mat(Cl*Cf.t())(Omegaidx));
  double temp3 = mean(arma::mat(square(Al)*square(Af).t())(Omegaidx));

  return 1/(temp1+temp2-temp3);
}



/*** R
# microbenchmark(SOLQL.R(Y, A.f, C.f, Omega, tau.l, tau),
#                SOLQL(  Y, A.f, C.f, Omega, tau.l, tau)
# )

# microbenchmark(SOL.TAU(Y, Omega, A.l, A.f, C.l, C.f),
#                SOLTAU(Y, Omega, A.l, A.f, C.l, C.f)
# )
*/

