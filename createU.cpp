#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

List createU(int N,int D,NumericVector qd,List T,NumericVector site,List V,List b) {	

	int id, index; List U(D);

	for(int d=0; d<D; d++){
		NumericMatrix Tmat = T(d);
		NumericMatrix Vmat = V(d);
		NumericMatrix bmat = b(d);
		NumericMatrix U1(N,qd(d)*(qd(d)-1)/2);
		for(int i=0; i<N; i++){
			id = site(i) - 1;	// site id
			index = 0;
			for(int j=0; j<(qd(d)-1); j++){
				for(int m=(j+1); m<qd(d); m++){
					U1(i,index) = bmat(id,j) * Tmat(i,m) * Vmat(m,m);
					index += 1;
				}
			}
		}
		U[d] = U1;
	}
	return U;
}