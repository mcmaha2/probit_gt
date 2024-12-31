#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericMatrix createTVAb(int N,int D,NumericVector qd,NumericMatrix TVAb,List TVA,List b,NumericVector site) {

	int end = 0; int id = 0;

	for(int i=0; i<N; i++){
		id = site(i) - 1;	// site id
		for(int d=0; d<D; d++){
			end += qd(d);
			NumericMatrix x1 = TVA(d);
			NumericMatrix x2 = b(d);
			end = qd(d);
			for(int j=0; j<end; j++){
				TVAb(i,d) += x1(i,j) * x2(id,j);
			}
		}
	}
	return TVAb;
}