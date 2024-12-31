#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

List createE(int N,int D,NumericVector qd,List T,NumericVector site,List A,List b) {	

	int id; List E(D);

	for(int d=0; d<D; d++){
		NumericMatrix Tmat = T(d);
		NumericMatrix Amat = A(d);
		NumericMatrix bmat = b(d);
		NumericMatrix E1(N,qd(d));
		for(int i=0; i<N; i++){
			id = site(i) - 1;	// site id
			E1(i,0) += Tmat(i,0) * bmat(id,0);
			for(int j=1; j<qd(d); j++){
				E1(i,j) += Tmat(i,j) * bmat(id,j);
				for(int k=0; k<j; k++){
					E1(i,j) += Tmat(i,j) * bmat(id,k) * Amat(j,k);
				}
			}
		}
		E[d] = E1;
	}
	return E;
}