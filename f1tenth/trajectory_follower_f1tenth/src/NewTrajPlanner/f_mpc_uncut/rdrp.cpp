#include "structures.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

/*
* rd_tilde_rp(...)
*	Input
*		quadrotor, mpc_rules, bag_of_many_things: structures containing system and problem parameters; see <structures.h> for details
*		X, U: matrices containing the state and control inputs across the time horizon, respectively
*	Output
*		rd_tilde: soft-constrained dual residual, refer to equation 27 of Richards and equation 7 of Boyd
*		rp: primal residual, refer to equation 7 of Boyd
*/
void rd_tilde_rp(struct System* quadrotor, struct MPC_Params* mpc_rules, Eigen::MatrixXf* X, Eigen::MatrixXf* U, struct Matrix_Set* BoMT){


	// MULTIPLY MATRIX 'C^T' BY VECTOR 'nu'; STORE INTO VECTOR 'Ctnu'
	(*BoMT).Ctnu_x = (*BoMT).nu;

	// Iterate over the number of time steps
	for(unsigned short int i = 0; i < (*mpc_rules).T-1; ++i)
	{
		
		(*BoMT).Ctnu_x.col(i) -= (*quadrotor).eig_A.transpose() * (*BoMT).nu.col(i+1);

	} // for(int i = 0; i < (*mpc_rules).T-1; ++i)

	(*BoMT).Ctnu_u = (*quadrotor).eig_B.transpose() * (*BoMT).nu;

	// SUBTRACT VECTOR 'b' FROM VECTOR 'Cz'; STORE INTO VECTOR 'rp'
	(*BoMT).rp.noalias() = *X - (*quadrotor).eig_B * (*U);
	(*BoMT).rp.col(0).noalias() -= (*quadrotor).eig_A * (*quadrotor).X0.col(0);

	for(int i = 1; i < (*mpc_rules).T; ++i)
	{
		(*BoMT).rp.col(i).noalias() -= (*quadrotor).eig_A * (*X).col(i-1);
	}

	// MULTIPLY MATRIX 'P^T-tilde' BY VECTOR 'd-tilde'; STORE INTO VECTOR 'PtTdt'
	(*BoMT).PtTdt_x.noalias() = (*quadrotor).eig_Fx_soft.transpose() * (*BoMT).d_tilde_x;
	(*BoMT).PtTdt_u.noalias() = (*quadrotor).eig_Fu_soft.transpose() * (*BoMT).d_tilde_u;

	// STORE 'gf' + 'g' + 'mu_12' * 'gp' + 'Ctnu' + 'PtTdt' INTO 'rd_tilde'
	// Compute rd_tilde = 2Hz + g + mu_12P^Td + C^Tnu + (P_tilde^T)(d_tilde)
	for (int i = 0; i < (*mpc_rules).T - 1; ++i)
	{

		//r~_d_x_i = 2QX + (k * F_x * d_x_i) + (F~_x_i * d~_x_i) - (A^T * v_(i+1)) + v_i + q
		// (*BoMT).rd_tilde_x.col(i) = (*BoMT).TwoQX.col(i) + (*mpc_rules).mu_12*(*quadrotor).eig_Fx_hard.transpose()*(*BoMT).dx.col(i) + (*mpc_rules).mu_13*(*quadrotor).eig_Fx_soft.transpose()*(*BoMT).d_tilde_x.col(i) - (*quadrotor).eig_A.transpose() * (*BoMT).nu.col(i+1) + (*BoMT).nu.col(i) + (*mpc_rules).eig_q[i];
		(*BoMT).rd_tilde_x.col(i) = (*BoMT).TwoQX.col(i) + (*mpc_rules).mu_12*(*quadrotor).eig_Fx_hard.transpose()*(*BoMT).dx.col(i) + (*mpc_rules).mu_13*(*quadrotor).eig_Fx_soft.transpose()*(*BoMT).d_tilde_x.col(i) - (*quadrotor).eig_A.transpose() * (*BoMT).nu.col(i) + (*BoMT).nu.col(i+1) + (*mpc_rules).eig_q[i];

		//r~_d_u_i = 2RU + (k * F_u * d_u_i) + (F~_u_i * d~_u_i) - (B^T * v_i)
		// (*BoMT).rd_tilde_u.col(i) = (*BoMT).TwoRU.col(i) + (*mpc_rules).mu_12*(*quadrotor).eig_Fu_hard.transpose()*(*BoMT).du.col(i) + (*mpc_rules).mu_13*(*quadrotor).eig_Fu_soft.transpose()*(*BoMT).d_tilde_u.col(i) - (*quadrotor).eig_B.transpose() * (*BoMT).nu.col(i) + (*mpc_rules).eig_q_lambda_k[i];
		(*BoMT).rd_tilde_u.col(i) = (*BoMT).TwoRU.col(i) + (*mpc_rules).mu_12*(*quadrotor).eig_Fu_hard.transpose()*(*BoMT).du.col(i) + (*mpc_rules).mu_13*(*quadrotor).eig_Fu_soft.transpose()*(*BoMT).d_tilde_u.col(i) - (*quadrotor).eig_B.transpose() * (*BoMT).nu.col(i) + (*mpc_rules).eig_q_lambda_k[i];

	}

	//r_d~_x_i = 2QX + (k * F_x * d_x_i) + (F~_x_i * d~_x_i) - (A^T * v_(i+1)) + v_i + q 
	// (*BoMT).rd_tilde_x.col((*mpc_rules).T-1) = (*BoMT).TwoQXf + (*mpc_rules).mu_12*(*quadrotor).eig_Fx_hard.transpose()*(*BoMT).dx.col((*mpc_rules).T-1) + (*mpc_rules).mu_13*(*quadrotor).eig_Fx_soft.transpose()*(*BoMT).d_tilde_x.col((*mpc_rules).T-1) + (*BoMT).nu.col((*mpc_rules).T-1) + (*mpc_rules).eig_q[(*mpc_rules).T-1];
	(*BoMT).rd_tilde_x.col((*mpc_rules).T-1) = (*BoMT).TwoQXf + (*mpc_rules).mu_12*(*quadrotor).eig_Fx_hard.transpose()*(*BoMT).dx.col((*mpc_rules).T-1) + (*mpc_rules).mu_13*(*quadrotor).eig_Fx_soft.transpose()*(*BoMT).d_tilde_x.col((*mpc_rules).T-2) + (*BoMT).nu.col((*mpc_rules).T-1) + (*mpc_rules).eig_q[(*mpc_rules).T-1];

	//r~_d_u_i = 2RU + (k * F_u * d_u_i) + (F~_u_i * d~_u_i) - (B^T * v_i)
	// (*BoMT).rd_tilde_u.col((*mpc_rules).T-1) = (*BoMT).TwoRU.col((*mpc_rules).T-1) + (*mpc_rules).mu_12*(*quadrotor).eig_Fu_hard.transpose()*(*BoMT).du.col((*mpc_rules).T-1) + (*mpc_rules).mu_13*(*quadrotor).eig_Fu_soft.transpose()*(*BoMT).d_tilde_u.col((*mpc_rules).T-1) - (*quadrotor).eig_B.transpose() * (*BoMT).nu.col((*mpc_rules).T-1) + (*mpc_rules).eig_q_lambda_k[(*mpc_rules).T-1];
	(*BoMT).rd_tilde_u.col((*mpc_rules).T-1) = (*BoMT).TwoRU.col((*mpc_rules).T-1) + (*mpc_rules).mu_12*(*quadrotor).eig_Fu_hard.transpose()*(*BoMT).du.col((*mpc_rules).T-1) + (*mpc_rules).mu_13*(*quadrotor).eig_Fu_soft.transpose()*(*BoMT).d_tilde_u.col((*mpc_rules).T-1) - (*quadrotor).eig_B.transpose() * (*BoMT).nu.col((*mpc_rules).T-1) + (*mpc_rules).eig_q_lambda_k[(*mpc_rules).T-1];

	// Note: negating rd here saves a little time down the road
	(*BoMT).rd_tilde_x *= -1;
	(*BoMT).rd_tilde_u *= -1;

	return;
}
