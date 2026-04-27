
// List include files
#include "structures.h"

// List namespaces
using namespace std;

/*
* dnudz(...)
*	Input
*		quadrotor, mpc_rules, bag_of_many_things: structures containing system and problem parameters; see <structures.h> for details
*	Output
*		rd_tilde, rp
*/
void dnudz(struct System* quadrotor, struct MPC_Params* mpc_rules, struct Matrix_Set* BoMT)
{

	// Calculate kappa*PTdT and mu_13*PtTdhPt
	for ( unsigned short int i = 0; i < (*mpc_rules).T; ++i )
	{

		// Calculate diag(d_x)*F_x
		(*BoMT).set_PTdP_x[i].noalias() = (*BoMT).dx.col(i).asDiagonal() * (*quadrotor).eig_Fx_hard;

		// Calculate k * (diag(d_x)*F_x)^T * (diag(d_x)*F_x)
		(*BoMT).set_PTdP_x[i] = (*mpc_rules).mu_12 * (*BoMT).set_PTdP_x[i].transpose() * (*BoMT).set_PTdP_x[i];

		// Calculate diag(d_u)*F_u
		(*BoMT).set_PTdP_u[i].noalias() = (*BoMT).du.col(i).asDiagonal() * (*quadrotor).eig_Fu_hard;

		// Calculate k * (diag(d_u)*F_u)^T * (diag(d_u)*F_u)
		(*BoMT).set_PTdP_u[i] = (*mpc_rules).mu_12 * (*BoMT).set_PTdP_u[i].transpose() * (*BoMT).set_PTdP_u[i];

		// Calculate mu_13 * F~_x^T * diag(d_hat_x) * F~_x
		(*BoMT).set_PtTdhPt_x[i].noalias() = (*quadrotor).eig_Fx_soft.transpose() *  (*BoMT).rho_i_x_vec.asDiagonal() * (*BoMT).rho_i_x_vec.asDiagonal() * (*BoMT).d_hat_x.col(i).asDiagonal() * (*quadrotor).eig_Fx_soft;

		// Calculate mu_13 * F~_u^T * diag(d_hat_u) * F~_u
		(*BoMT).set_PtTdhPt_u[i].noalias() = (*mpc_rules).mu_13 * (*quadrotor).eig_Fu_soft.transpose() * (*BoMT).d_hat_u.col(i).asDiagonal() * (*BoMT).d_hat_u.col(i).asDiagonal() * (*quadrotor).eig_Fu_soft;

	} // for ( unsigned short int i = 0; i < (*mpc_rules).T; ++i )


	// Iterate over the number of time steps
	for ( unsigned short int i = 0; i < (*mpc_rules).T; ++i )
	{
		// Calculate Phi_Q
		(*BoMT).set_Phi_x[i].noalias() = 2*(*mpc_rules).eig_R_rk[i] + (*BoMT).set_PTdP_x[i] +  (*BoMT).set_PtTdhPt_x[i]; // this is the Hessian of the cost function which includes the log-barrier function and KS function (only the blocks corresponding to the state)

		// Calculate (Phi_Q)^-1
		(*BoMT).set_QInv[i].noalias() = (*BoMT).set_Phi_x[i].inverse();

		// Calculate A * (Phi_Q)^-1
		(*BoMT).set_AQInv[i].noalias() = (*quadrotor).eig_A * (*BoMT).set_QInv[i];

		// Calculate Phi_R
		(*BoMT).set_Phi_u[i].noalias() = 2*(*mpc_rules).eig_R_lambda + (*BoMT).set_PTdP_u[i] +  (*BoMT).set_PtTdhPt_u[i];

		// Calculate (Phi_R)^-1
		(*BoMT).set_RInv[i].noalias() = (*BoMT).set_Phi_u[i].inverse();

		// Calculate B * (Phi_R)^-1
		(*BoMT).set_BRInv[i].noalias() = (*quadrotor).eig_B * (*BoMT).set_RInv[i];

	} // for ( unsigned short int i = 0; i < (*mpc_rules).T; ++i )

	//Y = CPhiInv * C^T
	(*BoMT).set_Y_onDiag[0].noalias() = (*BoMT).set_BRInv[0] * (*quadrotor).eig_B.transpose() + (*BoMT).set_QInv[0]; // why is the last term here?
	
	// Iterate over the number of time steps
	for ( unsigned short int i = 1; i < (*mpc_rules).T; ++i )
	{
	
		(*BoMT).set_Y_onDiag[i].noalias() = (*BoMT).set_AQInv[i-1] * (*quadrotor).eig_A.transpose() + (*BoMT).set_BRInv[i-1] * (*quadrotor).eig_B.transpose() + (*BoMT).set_QInv[i];
	
	} // for ( unsigned short int i = 1; i < (*mpc_rules).T; ++i )

	//CPhiInvrd = C * Phi^-1 * rd (store in beta to save space)
	(*BoMT).beta.col(0).noalias() = (*BoMT).set_QInv[0]*(*BoMT).rd_tilde_x.col(0) - (*BoMT).set_BRInv[0]*(*BoMT).rd_tilde_u.col(0);
	
	// Iterate over the number of time steps
	for( unsigned short int i = 1; i < (*mpc_rules).T; ++i )
	{
	
		(*BoMT).beta.col(i).noalias() = -(*BoMT).set_BRInv[i]*(*BoMT).rd_tilde_u.col(i) + (*BoMT).set_QInv[i]*(*BoMT).rd_tilde_x.col(i) - (*BoMT).set_AQInv[i-1]*(*BoMT).rd_tilde_x.col(i-1);
	
	} // for( unsigned short int i = 1; i < (*mpc_rules).T; ++i )

	// Calculate (negative) Beta
	(*BoMT).beta += (*BoMT).rp;

	//solve L L^T dnu = -Beta (Part 1: Solving for L
	(*BoMT).cholYblox_onDiag[0] = (*BoMT).set_Y_onDiag[0].llt().matrixL();
	
	// Iterate over the number of time steps
	for( unsigned short int i = 0; i < (*mpc_rules).T-1; )
	{
		
		//solve for off diagonal (it's in set_AQInv now and also negative)
		(*BoMT).cholYblox_onDiag[i].transpose().triangularView<Eigen::Upper>().solveInPlace<Eigen::OnTheRight>((*BoMT).set_AQInv[i]);
		++i;
		//solve for on-diagonal
		(*BoMT).set_Y_onDiag[i] -= (*BoMT).set_AQInv[i-1]*(*BoMT).set_AQInv[i-1].transpose();
		(*BoMT).cholYblox_onDiag[i] = (*BoMT).set_Y_onDiag[i].llt().matrixL();

	} // for( unsigned short int i = 0; i < (*mpc_rules).T-1; )

	//solve L L^T dnu = -Beta (Part 2: Solving L Theta = -Beta
	(*BoMT).cholYblox_onDiag[0].triangularView<Eigen::Lower>().solveInPlace((*BoMT).beta.col(0)); // is this n^3*(numblocks) flops? Boyd's paper describes a particular method to find the Cholesky factorization. Are these methods the same?
	
	// Iterate over the number of time steps
	for( unsigned short int i = 1; i < (*mpc_rules).T; ++i )
	{

		(*BoMT).beta.col(i) += (*BoMT).set_AQInv[i-1] * (*BoMT).beta.col(i-1);
		(*BoMT).cholYblox_onDiag[i].triangularView<Eigen::Lower>().solveInPlace((*BoMT).beta.col(i));

	} // for( unsigned short int i = 1; i < (*mpc_rules).T; ++i )

	//solve L L^T dnu = -Beta (Part 3: Solving L^T dnu = Theta
	(*BoMT).cholYblox_onDiag[(*mpc_rules).T-1].transpose().triangularView<Eigen::Upper>().solveInPlace((*BoMT).beta.col((*mpc_rules).T-1));
	
	// Iterate over the number of time steps
	for(int i = (*mpc_rules).T-2; i >= 0; --i )
	{
		
		(*BoMT).beta.col(i) += (*BoMT).set_AQInv[i].transpose() * (*BoMT).beta.col(i+1);
		(*BoMT).cholYblox_onDiag[i].transpose().triangularView<Eigen::Upper>().solveInPlace((*BoMT).beta.col(i));

	} // for( unsigned short int i = (*mpc_rules).T-2; i >= 0; --i )

	// ************ //
	// solve for dz //
	// ************ //

	// nu = nu + s*dnu
	(*BoMT).dnu = (*BoMT).beta; 
	(*BoMT).dz_u.noalias() = (*BoMT).rd_tilde_u + (*quadrotor).eig_B.transpose() * (*BoMT).dnu;
	(*BoMT).dz_x.noalias() = (*BoMT).rd_tilde_x - (*BoMT).dnu;

	// (*BoMT).dz_u.noalias() = (*BoMT).rd_tilde_u + (*quadrotor).eig_B.transpose() * (*BoMT).nu;
	// (*BoMT).dz_x.noalias() = (*BoMT).rd_tilde_x - (*BoMT).nu;

	// Iterate over the number of time steps
	for ( unsigned short int i = 0; i < (*mpc_rules).T-1; ++i )
	{
		
		(*BoMT).dz_u.col(i) = (*BoMT).set_RInv[i] * (*BoMT).dz_u.col(i);
		(*BoMT).dz_x.col(i) += (*quadrotor).eig_A.transpose() * (*BoMT).dnu.col(i+1);
		(*BoMT).dz_x.col(i) = (*BoMT).set_QInv[i] * (*BoMT).dz_x.col(i);

	} // for ( unsigned short int i = 0; i < (*mpc_rules).T-1; ++i )
	
	(*BoMT).dz_u.col((*mpc_rules).T-1) = (*BoMT).set_RInv[(*mpc_rules).T-1] * (*BoMT).dz_u.col((*mpc_rules).T-1);
	(*BoMT).dz_x.col((*mpc_rules).T-1) = (*BoMT).set_QInv[(*mpc_rules).T-1] * (*BoMT).dz_x.col((*mpc_rules).T-1);

	return;

} // void dnudz(struct System* quadrotor, struct MPC_Params* mpc_rules, struct Matrix_Set* BoMT)
