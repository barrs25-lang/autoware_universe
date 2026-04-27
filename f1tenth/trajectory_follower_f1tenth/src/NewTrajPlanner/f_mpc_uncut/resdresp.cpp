#include "structures.h"

/*
* resdresp(...)
*	Input
*		bag_of_many_things: structure containing system and problem parameters; see <structures.h> for details
*	Return
*		float containing sum of squares of norms of residuals
*/
float resdresp(struct Matrix_Set BoMT){
	
	float res;

	//Calculate THE SQUARE OF the Frobenius norm of the residual matrices, which is functionally equivalent to the Euclidean norm of the respective vector representation
	BoMT.res_p = BoMT.rp.squaredNorm();
	// cout << "BoMT.res_p:" << BoMT.res_p << endl;
	BoMT.res_dx = BoMT.rd_tilde_x.squaredNorm();
	// cout << "BoMT.res_dx:" << BoMT.res_dx << endl;
	BoMT.res_du = BoMT.rd_tilde_u.squaredNorm();

	res =BoMT.res_p + BoMT.res_dx + BoMT.res_du;

	if (isnan(res) || isinf(res))
	{
		cout << "BoMT.res_du:" << BoMT.res_du << " BoMT.res_dx:" << BoMT.res_dx << " BoMT.res_p:" << BoMT.res_p << endl;
	}

	// Return the overall residual
	return BoMT.res_p + BoMT.res_dx + BoMT.res_du;
}
