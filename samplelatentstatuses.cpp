#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericVector samplelatents(int N,NumericVector pone,NumericVector pzero,NumericMatrix YY,NumericMatrix G,NumericVector W,NumericVector U,NumericVector se,NumericVector sp) {
	int np, tmp, Gk, ck, as, ysum, id;
	float sek, spk, pi1, pi2, pistar;

	for(int i=0; i<N; i++){
		pi1 = pone(i);
		pi2 = pzero(i);
		np = YY(i,1);	// number of pools this individual was assigned to
		for(int j=0; j<np; j++){
			tmp = YY(i,(j+2+1)) - 1;	// obtain a pool this individual was assigned to
			Gk = G(tmp,0);		// observation for that pool
			ck = G(tmp,1);		// number in that pool
			as = G(tmp,2) - 1;		// assay used
			sek = se(as);		// sensitivity for that pool
			spk = sp(as);		// specificity for that pool
			ysum = 0;
			YY(i,0) = 0;			// coding trick to 'remove' individual i
			for(int k=0; k<ck; k++){
				id = G(tmp,(k+3))-1;		// indices of individuals in pool
				ysum = ysum + YY(id,0);	// to decide if an individual in pool is positive
			}
			pi1 = pi1*(sek*Gk + (1-sek)*(1-Gk));
			if(ysum > 0)
				pi2 = pi2*(sek*Gk + (1-sek)*(1-Gk));
			else
				pi2 = pi2*((1-spk)*Gk + spk*(1-Gk));
		}
		pistar = (pi1 / (pi1 + pi2));
		if(U(i) < pistar)
			YY(i,0) = 1;
		else
			YY(i,0) = 0;
		W(i) = YY(i,0);
	}
	return W;
}