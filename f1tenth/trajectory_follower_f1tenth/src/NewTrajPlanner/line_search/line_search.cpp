#include <line_search.h>

float objective_function_value(F_MPC_UNCUT* fmpc, Eigen::MatrixXf* X, Eigen::MatrixXf* U, Eigen::MatrixXf* dX, Eigen::MatrixXf* dU, float s, int segment_number)
{

	Eigen::MatrixXf Xs = (*X) - s*(*dX);
	Eigen::MatrixXf Us = (*U) - s*(*dU);

	Eigen::MatrixXf* tildeR = new Eigen::MatrixXf[(*fmpc).mpc_params.T];
	Eigen::MatrixXf* tildeq = new Eigen::MatrixXf[(*fmpc).mpc_params.T];
	Eigen::MatrixXf* Zk = new Eigen::MatrixXf[(*fmpc).mpc_params.T];
	Eigen::MatrixXf goal_n = Eigen::MatrixXf::Zero((*fmpc).quadrotor.n,1);
	Eigen::MatrixXf ZkTRZk;
	Eigen::MatrixXf ZfkT_Rf_Zfk;
	Eigen::MatrixXf qTZk;
	Eigen::MatrixXf qfT_Zfk;
	Eigen::MatrixXf ZkTRZk_qTZk;
	Eigen::MatrixXf ZfkT_Rf_Zfk__qfT_Zfk;

	bool constraint_violation;

	Eigen::MatrixXf Obj_vec = Eigen::MatrixXf::Zero((*fmpc).mpc_params.T,1);
	Eigen::MatrixXf Obj = Eigen::MatrixXf::Zero(1,1);
	Eigen::MatrixXf PhiZ = Eigen::MatrixXf::Zero(1,1);
	Eigen::MatrixXf ThetaZ = Eigen::MatrixXf::Zero(1,1);
	Eigen::MatrixXf kappaPhiZ = Eigen::MatrixXf::Zero(1,1);
	Eigen::MatrixXf rhoInverseThetaZ = Eigen::MatrixXf::Zero(1,1);
	Eigen::MatrixXf temp;
	temp = Eigen::MatrixXf::Zero(1,1);

	goal_n(0,0) = (*fmpc).goal(0,0);
	goal_n(1,0) = (*fmpc).goal(1,0);
	goal_n(2,0) = (*fmpc).goal(2,0);
	goal_n(12,0) = (*fmpc).goal(3,0);

	for (int i = 0; i < (*fmpc).mpc_params.T-1; i++)
	{
		tildeR[i] = Eigen::MatrixXf::Zero((*fmpc).quadrotor.n+(*fmpc).quadrotor.m,(*fmpc).quadrotor.n+(*fmpc).quadrotor.m);
		tildeq[i] = Eigen::MatrixXf::Zero((*fmpc).quadrotor.n+(*fmpc).quadrotor.m,1);
		Zk[i] =  Eigen::MatrixXf::Zero((*fmpc).quadrotor.n+(*fmpc).quadrotor.m,1);
		
		tildeR[i].block(0,0,(*fmpc).quadrotor.n,(*fmpc).quadrotor.n) = (*fmpc).mpc_params.eig_R_rk[i]; 
		tildeR[i].block(0,(*fmpc).quadrotor.n,(*fmpc).quadrotor.n,(*fmpc).quadrotor.m) = (*fmpc).mpc_params.eig_R_r_lambda_k[i]; 
		tildeR[i].block((*fmpc).quadrotor.n,0,(*fmpc).quadrotor.m,(*fmpc).quadrotor.n) = (*fmpc).mpc_params.eig_R_r_lambda_k[i].transpose();
		tildeR[i].block((*fmpc).quadrotor.n,(*fmpc).quadrotor.n,(*fmpc).quadrotor.m,(*fmpc).quadrotor.m) = (*fmpc).mpc_params.eig_R_lambda; 

		// cout << "q_rk: " << (*fmpc).mpc_params.eig_q_rk[i].transpose() << endl;
		tildeq[i].block(0,0,(*fmpc).quadrotor.n-2,1) = (*fmpc).mpc_params.eig_q_rk[i];
		tildeq[i](12,0) = (*fmpc).mpc_params.q_psi;
		tildeq[i](13,0) = 0;
		tildeq[i].block((*fmpc).quadrotor.n,0,(*fmpc).quadrotor.m,1) = (*fmpc).mpc_params.eig_q_lambda_k[i];

		Zk[i].block(0,0,(*fmpc).quadrotor.n,1) = Xs.col(i) - goal_n;
		Zk[i].block((*fmpc).quadrotor.n,0,(*fmpc).quadrotor.m,1) = Us.col(i);

		for (int ii = 0; ii < (*fmpc).mpc_params.h[segment_number].rows(); ii++)
		{
			temp = (*fmpc).mpc_params.P[segment_number].row(ii)*(Zk[i]);
			if ((*fmpc).mpc_params.h[segment_number](ii,i) - temp(0,0) <= 0)
			{
				PhiZ(0,0) += 1000000;
			}
			else
			{
				PhiZ(0,0) += -1*log((*fmpc).mpc_params.h[segment_number](ii,i) - temp(0,0) );
			}
			ThetaZ(0,0) += log(1 + exp((*fmpc).mpc_params.mu_13 * ( -(*fmpc).mpc_params.h_tilde[segment_number](ii,i) + temp(0,0) ) ) );
		}
		
		ZkTRZk = Zk[i].transpose()*tildeR[i]*Zk[i];
		qTZk = tildeq[i].transpose()*Zk[i];
		ZkTRZk_qTZk = ZkTRZk + qTZk;
		Obj_vec.block(i,0,1,1) = ZkTRZk_qTZk;

	}
	
	kappaPhiZ(0,0) = (*fmpc).mpc_params.mu_12*PhiZ(0,0); 
	rhoInverseThetaZ(0,0) = (1/(*fmpc).mpc_params.mu_13)*ThetaZ(0,0);

	ZfkT_Rf_Zfk = (Xs.col((*fmpc).mpc_params.T-1) - goal_n).transpose() * (*fmpc).mpc_params.eig_R_rf * (Xs.col((*fmpc).mpc_params.T-1) - goal_n);

	// cout << "Obj_vec.sum: " << Obj_vec.sum() << " term: " << ZfkT_Rf_Zfk(0,0) << " kappaPhiZ: " <<kappaPhiZ(0,0) << " rhoInverseThetaZ: " << rhoInverseThetaZ(0,0) << endl;

	Obj(0,0) = Obj_vec.sum() + ZfkT_Rf_Zfk(0,0) + kappaPhiZ(0,0) + rhoInverseThetaZ(0,0);

	delete[] tildeR;
	delete[] tildeq;
	delete[] Zk;

	return Obj(0,0);

} 

// void objective_function_gradient(F_MPC_UNCUT* fmpc, Eigen::MatrixXf* X, Eigen::MatrixXf* U)
// {

// 	Eigen::MatrixXf* tildeR = new Eigen::MatrixXf[(*fmpc).mpc_params.T];
// 	Eigen::MatrixXf* tildeq = new Eigen::MatrixXf[(*fmpc).mpc_params.T];
// 	Eigen::MatrixXf* Zk = new Eigen::MatrixXf[(*fmpc).mpc_params.T];
// 	Eigen::MatrixXf* gradObj = new Eigen::MatrixXf[(*fmpc).mpc_params.T];
// 	Eigen::MatrixXf goal_n = Eigen::MatrixXf::Zero((*fmpc).quadrotor.n,1);

// 	goal_n(0,0) = (*fmpc).goal(0,0);
// 	goal_n(1,0) = (*fmpc).goal(1,0);
// 	goal_n(2,0) = (*fmpc).goal(2,0);
// 	goal_n(12,0) = (*fmpc).goal(3,0);

// 	for (int i = 0; i < (*fmpc).mpc_params.T-1; i++)
// 	{
// 		tildeR[i] = Eigen::MatrixXf::Zero((*fmpc).quadrotor.n+(*fmpc).quadrotor.m,(*fmpc).quadrotor.n+(*fmpc).quadrotor.m);
// 		tildeq[i] = Eigen::MatrixXf::Zero((*fmpc).quadrotor.n+(*fmpc).quadrotor.m,1);
// 		Zk[i] =  Eigen::MatrixXf::Zero((*fmpc).quadrotor.n+(*fmpc).quadrotor.m,1);
		
// 		tildeR[i].block(0,0,(*fmpc).quadrotor.n,(*fmpc).quadrotor.n) = (*fmpc).mpc_params.eig_R_rk[i]; 
// 		tildeR[i].block(0,(*fmpc).quadrotor.n,(*fmpc).quadrotor.n,(*fmpc).quadrotor.m) = (*fmpc).mpc_params.eig_R_r_lambda_k[i]; 
// 		tildeR[i].block((*fmpc).quadrotor.n,0,(*fmpc).quadrotor.m,(*fmpc).quadrotor.n) = (*fmpc).mpc_params.eig_R_r_lambda_k[i].transpose();
// 		tildeR[i].block((*fmpc).quadrotor.n,(*fmpc).quadrotor.n,(*fmpc).quadrotor.m,(*fmpc).quadrotor.m) = (*fmpc).mpc_params.eig_R_lambda; 

// 		// cout << "q_rk: " << (*fmpc).mpc_params.eig_q_rk[i].transpose() << endl;
// 		tildeq[i].block(0,0,(*fmpc).quadrotor.n-2,1) = (*fmpc).mpc_params.eig_q_rk[i];
// 		tildeq[i](12,0) = (*fmpc).mpc_params.q_psi;
// 		tildeq[i](13,0) = 0;
// 		tildeq[i].block((*fmpc).quadrotor.n,0,(*fmpc).quadrotor.m,1) = (*fmpc).mpc_params.eig_q_lambda_k[i];

// 		Zk[i].block(0,0,(*fmpc).quadrotor.n,1) = (*X).col(i) - goal_n;
// 		Zk[i].block((*fmpc).quadrotor.n,0,(*fmpc).quadrotor.m,1) = (*U).col(i);

// 		for (int j = 0; j < (*fmpc).bomt.dx.rows(); j++)
// 		{
// 			= (*fmpc).mpc_params.mu_12*(*fmpc).mpc_params.P[segment_number]
// 		}

// 		gradObj[i] = 2*tildeR[i]*Zk[i] + tildeq[i];

// 	}

// }

// void objective_function_hessian()
// {



// }

float rd_tilde_rp_norm(F_MPC_UNCUT* fmpc)
{

	Eigen::MatrixXf rd_tilde = Eigen::MatrixXf::Zero((*fmpc).bomt.rd_tilde_x.rows() + (*fmpc).bomt.rd_tilde_u.rows(), 1);
	Eigen::MatrixXf rd_tilde_rp = Eigen::MatrixXf::Zero(1,1);

	// Need to input X,U, etc, and call fmpc computations to get correct rd_tilde and rp.

	rd_tilde.block(0,0,(*fmpc).bomt.rd_tilde_x.rows(),1) = (*fmpc).bomt.rd_tilde_x;
	rd_tilde.block((*fmpc).bomt.rd_tilde_x.rows(),0,(*fmpc).bomt.rd_tilde_u.rows(),1) = (*fmpc).bomt.rd_tilde_u;
	rd_tilde_rp(0,0) = ( rd_tilde*(*fmpc).bomt.rp.transpose() ).norm();

	return rd_tilde_rp(0,0);  


} 

float dichotomous_line_search(F_MPC_UNCUT* fmpc, Eigen::MatrixXf* X, Eigen::MatrixXf* U, Eigen::MatrixXf* dX, Eigen::MatrixXf* dU, float a1, float b1, float l, float epsilon, int segment_number, float start_s)
{

	float interval_of_uncertainty = b1-a1;
	float prev_interval_of_uncertainty = b1-a1;

	int k = 1;

	float theta_lambda_k;
	float theta_mu_k;

	vector<float> lambda_k;
	vector<float> mu_k;
	vector<float> a_k;
	vector<float> b_k;

	a_k.push_back(a1 - start_s);
	b_k.push_back(b1 - start_s);

	// cout << "initial interval: " << a1 << ", " << b1 << endl;

	while(interval_of_uncertainty >= l && k < 20)
	{

		lambda_k.push_back( 0.5 * ( a_k[k-1] + b_k[k-1] ) - epsilon );
		mu_k.push_back( 0.5 * ( a_k[k-1] + b_k[k-1] ) + epsilon );

		theta_lambda_k = objective_function_value(fmpc, X, U, dX, dU, lambda_k[k-1], segment_number);
		theta_mu_k = objective_function_value(fmpc, X, U, dX, dU, mu_k[k-1], segment_number);

		// theta_lambda_k = rd_tilde_rp_norm(fmpc);
		// theta_mu_k = rd_tilde_rp_norm(fmpc);

// cout << "obj lam: " << theta_lambda_k << "obj mu: " << theta_mu_k << endl;

		if (theta_lambda_k < theta_mu_k)
		{
			a_k.push_back( a_k[k-1] );
			b_k.push_back( mu_k[k-1] );
		}
		else if (theta_lambda_k > theta_mu_k)
		{
			a_k.push_back( lambda_k[k-1] );
			b_k.push_back( b_k[k-1] );
		}
		else
		{
			a_k.push_back( lambda_k[k-1] );
			b_k.push_back( mu_k[k-1] );
		}

	// cout << "new interval: " << a_k[k] << ", " << b_k[k] << endl;
		interval_of_uncertainty = b_k[k] - a_k[k];
		k++;

	}

	// cout << "exiting line search: " << a_k[k-1] << ", " << b_k[k-1] << endl;
	return (lambda_k.end()[-1] + mu_k.end()[-1]) / 2;

}

float newton_method_line_search()
{

}