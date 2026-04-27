#include "structures.h"
#include <stdio.h>
#include <iostream>
#include <fstream>

// Local method
// gfgphp: computes various values needed for solving mpc
// Input: quadrotor, mpc_rules, bag_of_many_things: structures containing system and problem parameters; see <structures.h> for details
// X, U: matrices containing the state and control inputs across the time horizon, respectively
// Output: P^Td: gradient of log barrier, see Boyd chapter 3 section A
// 2Hz: component of dual residual calculation
void gfgphp(struct System* quadrotor, struct MPC_Params* mpc_rules, Eigen::MatrixXf* X, Eigen::MatrixXf* U, struct Matrix_Set* BoMT)
{

	// **************** //
	// Calculation of d //
	// **************** //

	// Compute Fx * X
	(*BoMT).FxX.noalias() = (*quadrotor).eig_Fx_hard * (*X);

	// Iterate over the number of time steps
	for (unsigned short int i = 0; i < (*mpc_rules).T; i++)
	{
		(*BoMT).dx.col(i) = ((-(*BoMT).FxX).col(i) + (*quadrotor).eig_x_hard_bounds.col(i)).cwiseInverse();
	}

	// Compute Fu * U
	(*BoMT).FuU.noalias() = (*quadrotor).eig_Fu_hard * (*U);

	// Compute du[i] = 1 / (fu[i%m] - FuU[i])
	(*BoMT).du.noalias() = ((-(*BoMT).FuU).colwise() + (*quadrotor).eig_u_hard_bounds).cwiseInverse();

	// ******************* //
	// Calcluation of P^Td //
	// ******************* //

	// Compute Ptdx = Fx^T dx 
	(*BoMT).FxTdx.noalias() = (*quadrotor).eig_Fx_hard.transpose() * (*BoMT).dx;

	// Compute Ptdu = Fu^T du
	(*BoMT).FuTdu.noalias() = (*quadrotor).eig_Fu_hard.transpose() * (*BoMT).du;
 
 	// ****************** //
	// Calculation of 2Hz //
 	// ****************** //

	// Compute X_mod and U_mod
	// Iterate over the number of time steps minus one
	for (int i = 0; i < (*mpc_rules).T-1; i++)
	{

		(*BoMT).X_mod.col(i) = (*mpc_rules).eig_R_rk[i]*( (*X).col(i) - (*BoMT).X_goal);
		(*BoMT).U_mod.col(i) = (*U).col(i);

	} // for (int i = 0; i < (*mpc_rules).T-1; i++)

	// Compute X_mod or U_mod for the terminal time step
	(*BoMT).X_mod.col((*mpc_rules).T-1) = (*mpc_rules).eig_R_rf*( (*X).col((*mpc_rules).T-1) - (*BoMT).X_goal);
	(*BoMT).U_mod.col((*mpc_rules).T-1) = (*U).col((*mpc_rules).T-1);

	// Compute 2*R*U and 2*Q*X (as separated H)
	(*BoMT).TwoRU.noalias() = 2 * (*mpc_rules).eig_R_lambda * (*BoMT).U_mod;
	(*BoMT).TwoQX.noalias() = 2 * (*BoMT).X_mod;
	(*BoMT).TwoQXf.noalias() = 2 * (*BoMT).X_mod.col((*mpc_rules).T-1);

	return;

} // void gfgphp(struct System* quadrotor, struct MPC_Params* mpc_rules, Eigen::MatrixXf* X, Eigen::MatrixXf* U, struct Matrix_Set* BoMT)


// Local method
// dtildhat: compute dtilde for the soft constraints
// Input: quadrotor, mpc_rules, bag_of_many_things: structures containing system and problem parameters; see <structures.h> for details
// X, U: matrices containing the state and control inputs across the time horizon, respectively
// Output: d_tilde: component of soft penalty calculation, see equation 22 of Richards
// d_hat: component of soft penalty calculation, see equation 24 of Richards
void dtildhat(struct System* quadrotor, struct MPC_Params* mpc_rules, Eigen::MatrixXf* X, Eigen::MatrixXf* U, struct Matrix_Set* BoMT)
{

	// Compute Fz
	(*BoMT).FxtildeX.noalias() = (*quadrotor).eig_Fx_soft * (*X);
	(*BoMT).FutildeU.noalias() = (*quadrotor).eig_Fu_soft * (*U);

	// Declare float to store computation of exponential
	float e_i;

	// Float arrays to store the d values
	float* dtilde = new float[(*quadrotor).num_soft_x* (*mpc_rules).T]; 
	float* dhat = new float[(*quadrotor).num_soft_x* (*mpc_rules).T]; 

	// Booleans to indicate the success of finding a vector in the nullspace of P
	bool trying = true, success = false;
	Eigen::MatrixXf Px_kernel;
	Eigen::MatrixXf d_0; 
	e_i = 0.0;

	// find rho_i_x
	if ((*mpc_rules).centering_soft_constraints)
	{
		if (((*BoMT).v > 1) && ((*BoMT).initial_run))
		{
			Px_kernel = (*quadrotor).eig_Fx_soft.transpose().fullPivLu().kernel();

			for (int i = 0; i < (*quadrotor).num_soft_x; i++)
			{
				for (int j = 0; j < 10; j++)
				{
					if (Px_kernel.col(i).lpNorm<Eigen::Infinity>() < 1)
					{
						d_0 = Eigen::MatrixXf::Zero((*quadrotor).num_soft_x,1);
						d_0.col(0) = Px_kernel.col(i);
						trying = false;
						for (int k = 0; k < (*quadrotor).num_soft_x; k++)
						{
							d_0(k,0) = d_0(k,0) + 0.00001;
							if (d_0(k,0) <= 0)
							{
								trying = true;
								j = 9;
								break;
							}
						}
					}
					else
					{
						Px_kernel.col(i) = Px_kernel.col(i) * 0.999;
					}
				}

				if (!trying)
				{
					success = true;
					break;
				}
			}

			if (success)
			{
				//richards eqn 35
				for (int i = 0; i < (*quadrotor).num_soft_x; i++)
				{
					// (*BoMT).rho_i_x(i,0) = ( 1 / (*quadrotor).eig_x_soft_bounds[i] ) * log((1/d_0(i,0)) - 1);
					(*BoMT).rho_i_x(i,0) = ( 1 / (*quadrotor).eig_x_soft_bounds(i,0) ) * log((1/d_0(i,0)) - 1);
					(*BoMT).rho_i_x_vec[i] = (*BoMT).rho_i_x(i,0);
				}
			}
			else
			{
				for (int i = 0; i < (*quadrotor).num_soft_x; i++)
				{
					(*BoMT).rho_i_x(i,0) = (*mpc_rules).mu_13_initial;
					(*BoMT).rho_i_x_vec[i] = (*BoMT).rho_i_x(i,0);
				}
			}

		}
		else if ((*BoMT).v == 0)
		{
			for (int i = 0; i < (*quadrotor).num_soft_x; i++)
			{
				(*BoMT).rho_i_x(i,0) = (*mpc_rules).mu_13_initial;
				(*BoMT).rho_i_x_vec[i] = (*BoMT).rho_i_x(i,0);
			}
		}
		else
		{
			for (int i = 0; i < (*quadrotor).num_soft_x; i++)
			{
				(*BoMT).rho_i_x(i,0) = (*BoMT).rho_i_x(i,0) * (*mpc_rules).gamma;
				(*BoMT).rho_i_x_vec[i] = (*BoMT).rho_i_x(i,0);
			}
		}
	}
	else
	{
		for (int i = 0; i < (*quadrotor).num_soft_x; i++)
		{
			(*BoMT).rho_i_x(i,0) = (*mpc_rules).mu_13_initial;
			(*BoMT).rho_i_x_vec[i] = (*BoMT).rho_i_x(i,0);
		}
	}


	// Iterate over the number of time steps, starting from one
	for(unsigned short int i = 1; i < (*mpc_rules).T+1; ++i) //check optimality w.r.t. row-major/col-major order
	{

		// Iterate over the number of soft constraints, starting from the ith set of constraints
		for(unsigned short int j = (i-1)*(*quadrotor).num_soft_x; j < (i)*(*quadrotor).num_soft_x; ++j)
		{

			// Compute the exponential
			e_i = exp((*mpc_rules).mu_13 * ((*BoMT).FxtildeX.data()[j] - (*quadrotor).eig_x_soft_bounds.data()[j%(*quadrotor).num_soft_x])); //corresponds to e_i+ in eq19 of Richards
			
			// If the exponential is finite
			if(!isinf(e_i))
			{
			
				// Compute dtilde and dhat
				dtilde[j] = ( e_i < 1 ) ? ( e_i/(1.0f + e_i) ) : ( 1.0f/ (1.0f + 1.0f/e_i) ); //evaluates eq22 of Richards using the smaller of e_i+ and e_i- for state segment
				dhat[j] = e_i / ((1.0f + e_i) * (1.0f + e_i)); //evaluates eq24 of Richards for state segment
			
			} // if(!isinf(e_i))
			else 
			{
			
				dtilde[j] = 1.0f;
				dhat[j] = 0.0f;
			
			} // if(!isinf(e_i))

		} // for(unsigned short int j = (i-1)*(*quadrotor).num_soft_x; j < (i)*(*quadrotor).num_soft_x; ++j)

	} // for(int i = 1; i < (*mpc_rules).T+1; ++i)

	// Set d_tilde_x and d_hat_x matrices to dtilde and dhat arrays
	(*BoMT).d_tilde_x = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(dtilde, (*quadrotor).num_soft_x, (*mpc_rules).T);
	(*BoMT).d_hat_x = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(dhat, (*quadrotor).num_soft_x, (*mpc_rules).T);

	// New float arrays to store d for control constraints
	dtilde = new float[(*quadrotor).num_soft_u* (*mpc_rules).T]; 
	dhat = new float[(*quadrotor).num_soft_u* (*mpc_rules).T]; 

	e_i = 0.0;

	// Iterate over the number of time steps
	for(unsigned short int i = 0; i < (*mpc_rules).T; ++i)
	{

		// Iterate over the number of control constraints
		for(int j = 0; j < (*quadrotor).num_soft_u; ++j)
		{

			// Compute the exponential
			e_i = exp( (*mpc_rules).mu_13 * ((*BoMT).FutildeU(j,i) - (*quadrotor).eig_u_soft_bounds(j) ) ); //corresponds to e_i+ in eq19 of Richards
			

			if(!isinf(e_i))
			{
				
				// Compute dtilde and dhat
				dtilde[j+i*(*quadrotor).num_soft_u] = (e_i<1)?(e_i/(1.0f + e_i)):(1.0f/(1.0f + 1.0f/e_i)); //evaluates eq22 of richards using the smaller of e_i+ and e_i- for control segment
				dhat[j+i*(*quadrotor).num_soft_u] = e_i / ((1.0f + e_i) * (1.0f + e_i)); //evaluates eq24 of Richards for control segment
			
			} // if(!isinf(e_i))
			else 
			{
				
				// Compute dtilde and dhat
				dtilde[j+i*(*quadrotor).num_soft_u] = 1.0f;
				dhat[j+i*(*quadrotor).num_soft_u] = 0.0f;
			
			} // if(!isinf(e_i))

		} // for(int j = 0; j < (*quadrotor).num_soft_u; ++j)

	} // for(unsigned short int i = 0; i < (*mpc_rules).T; ++i)

	// Set d_tilde_u and d_hat_u matrices to dtilde and dhat arrays
	(*BoMT).d_tilde_u = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(dtilde, (*quadrotor).num_soft_u, (*mpc_rules).T);
	(*BoMT).d_hat_u = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(dhat, (*quadrotor).num_soft_u, (*mpc_rules).T);

	return;

} // void dtildhat(struct System* quadrotor, struct MPC_Params* mpc_rules, Eigen::MatrixXf* X, Eigen::MatrixXf* U, struct Matrix_Set* BoMT)

