#ifndef STRUCTURES_H_
#define STRUCTURES_H_

#include <Eigen/Dense>
#include <chrono>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

struct System {

	// Dynamical system parameters
	int n; 									// Dimension of the state
	int m; 									// Dimension of the input
	Eigen::MatrixXf eig_A;					// State transition matrix in discrete time (n X n)
	Eigen::MatrixXf eig_B;					// Control effectiveness matrix in discrete time (n X m)
	Eigen::MatrixXcf eig_A_complex;			// Dummy matrix used in Cayley-Hamilton Theorem for computing matrix exponential (converting continuous time eig_A to discrete time eig_A)		
	int noise_source;						// Noise in dynamical model source				
	Eigen::MatrixXf noise_w; 				// Vector storing the noise (n X 1)

	// Physical properties of UAV
	float length; 							// Length of UAV [meters]
	float width;							// Width of UAV [meters]
	float height;							// height of UAV [meters]
	float half_camera_FOV;					// FOV of camera [radians] this is a conical approximation of what may be a rectangular pyramid
	float mass;								// Mass of UAV [kg]
	float moment_arm_l;						// Length of UAV's moment arm [m]
	float cT;								// Drag coefficient of UAV's propellers
	Eigen::MatrixXf M; 						// Mixer used to convert control inputs to individual thrust forces
	Eigen::MatrixXf inertia_matrix; 		// Matrix of inertia of UAV [kg*m^2] 
	Eigen::MatrixXf inertia_matrix_inv; 	// Inverse of matrix of inertia of UAV

	// Hard constraint matrices and parameters
	int num_hard_x = 32; 					// Number of hard state constraints
	Eigen::MatrixXf eig_Fx_hard;			// F matrix for hard constraints on position, velocity, acceleration, and jerk (dynamic x n-2) (think of each row in the matrix being the normal to a plane in the state-space. All states must not "cross" the plane)
	Eigen::MatrixXf eig_Fpsi_hard;			// F matrix for hard constraints on psi and dpsi (2 x 2)
	Eigen::MatrixXf eig_Fz_hard_ceiling;	// Additional F matrix for altitude constraint (2 x n)
	// Eigen::VectorXf eig_x_hard_bounds;		// f vector for hard constraints on position, velocity, acceleration, and jerk (think of this as the offset for the plane defined by F. Really F and f form a plane equation F{1,1}x + F_{1,2}y + F_{1,3}z + ... = f_1)
	Eigen::MatrixXf eig_x_hard_bounds;		// f vector for hard constraints on position, velocity, acceleration, and jerk (think of this as the offset for the plane defined by F. Really F and f form a plane equation F{1,1}x + F_{1,2}y + F_{1,3}z + ... = f_1)
	Eigen::Vector2f eig_psi_hard_bounds;	// f vector for hard constraints on psi and dpsi (2 x 1)
	Eigen::Vector2f eig_z_hard_bounds;		// f vector for hard constraints on altitude (2 x 1)
	int num_hard_u; 						// Number of hard control constraints 
	Eigen::MatrixXf eig_Fu_hard;			// F matrix for hard constraints on control input (2*m X m)
	Eigen::VectorXf eig_u_hard_bounds;		// f vector for hard constraints on control input (2*m X 1)

	Eigen::MatrixXf eig_Fx_hard_boundary_conditions; // F matrix for boundary conditions (2*n X n) The reason it is 2*n is that we have to use inequality constraints to represent equality constraints. Here's an example: here's a boundary condition r(0) = r0, here's the equivalent condition that uses inequality constraints r(0) <= r0, r(0) >= r0. If both constraints are to be satisfied, then the only option is that r(0) = r0.
	Eigen::MatrixXf eig_x_hard_boundary_conditions; // f vector for boundary conditions (2*n X n) 
	int num_obstacle_constraints = 0;
	float phi_max, theta_max, T_min, T_max;				// Maximum angle in radians for roll and pitch angles

	// Soft constraint matrices and parameters
	int num_soft_x = 32; 					// Number of soft state constraints
	Eigen::MatrixXf eig_Fx_soft; 			// F matrix for soft constraints on position, velocity, acceleration, and jerk (dynamic x n-2) (think of each row in the matrix being the normal to a plane in the state-space. All states must not "cross" the plane)
	Eigen::MatrixXf eig_Fpsi_soft;			// F matrix for soft constraints on psi and dpsi (2 x 2)
	Eigen::MatrixXf eig_Fz_soft_ceiling;	// Additional F matrix for altitude constraint (2 x n)
	Eigen::MatrixXf eig_x_soft_bounds;		// f vector for soft constraints on position, velocity, acceleration, and jerk (think of this as the offset for the plane defined by F. Really F and f form a plane equation F{1,1}x + F_{1,2}y + F_{1,3}z + ... = f_1)
	// Eigen::VectorXf eig_x_soft_bounds;		// f vector for soft constraints on position, velocity, acceleration, and jerk (think of this as the offset for the plane defined by F. Really F and f form a plane equation F{1,1}x + F_{1,2}y + F_{1,3}z + ... = f_1)
	Eigen::Vector2f eig_psi_soft_bounds;	// f vector for soft constraints on psi and dpsi (2 x 1)
	Eigen::Vector2f eig_z_soft_bounds;		// f vector for soft constraints on altitude (2 x 1)
	int num_soft_u; 						// Number of soft control constraints 
	Eigen::MatrixXf eig_Fu_soft;			// F matrix for soft constraints on control input (2*m X m)
	Eigen::VectorXf eig_u_soft_bounds;		// f vector for soft constraints on control input (2*m X 1)
	float* eig_Fx_soft_array;
	float* eig_x_soft_bounds_array;

	Eigen::MatrixXf eig_Fx_soft_boundary_conditions; // F matrix for boundary conditions (2*n X n) The reason it is 2*n is that we have to use inequality constraints to represent equality constraints. Here's an example: here's a boundary condition r(0) = r0, here's the equivalent condition that uses inequality constraints r(0) <= r0, r(0) >= r0. If both constraints are to be satisfied, then the only option is that r(0) = r0.
	Eigen::MatrixXf eig_x_soft_boundary_conditions; // f vector for boundary conditions (2*n X n) 

	Eigen::VectorXf X0;						// State vector (n X 1)
	Eigen::VectorXf prevX0;						// State vector (n X 1)

	float* x0;
	Eigen::VectorXf eig_goal_difference;

	Eigen::MatrixXf delta_k; 				// State of thrust system (2 x n_{\rm t})
	Eigen::MatrixXf u1_k; 					// Total thrust to real system (1 x n_{\rm t})
	Eigen::MatrixXf T_k; 					// Individual thrust forces (4 x n_{\rm t})
	Eigen::MatrixXf T_k_min; 					// Individual thrust forces (4 x n_{\rm t})
	Eigen::MatrixXf T_k_max; 					// Individual thrust forces (4 x n_{\rm t})
	Eigen::MatrixXf u1_dot_k; 				// Time derivative of total thrust to real system (1 x n_{\rm t})
	Eigen::MatrixXf eig_phi;
	Eigen::MatrixXf eig_theta;

	Eigen::MatrixXf X_goal;					// Goal vector (n X 1)

};

struct MPC_Params {
	//Note: We *probably* don't need double-precision, so I'm switching to floats to speed things up
	float* Q; //the state costs
	// Eigen::MatrixXf eig_Q;
	Eigen::MatrixXf* eig_R_rk = new Eigen::MatrixXf[1];
	// Eigen::MatrixXf eig_Q_tilde;
	Eigen::MatrixXf eig_R_r_tilde;
	// Eigen::MatrixXf eig_Q_ru_tilde;
	Eigen::MatrixXf eig_R_r_lambda_tilde;
	Eigen::MatrixXf* eig_R_r_lambda_k = new Eigen::MatrixXf[1];
	Eigen::MatrixXf* eig_q = new Eigen::MatrixXf[1];
	Eigen::MatrixXf* eig_q_rk = new Eigen::MatrixXf[1];
	Eigen::MatrixXf eig_q_r_tilde;
	float q_psi;
	// Eigen::MatrixXf eig_q_u_tilde;
	Eigen::MatrixXf eig_q_lambda_tilde;
	Eigen::MatrixXf* eig_q_lambda_k = new Eigen::MatrixXf[1];
	float* R; //the control costs
	Eigen::MatrixXf eig_R_lambda;
	Eigen::MatrixXf eig_R_lambda_init;
	Eigen::MatrixXf eig_R_lambda_tilde;

	Eigen::MatrixXf* P = new Eigen::MatrixXf[1];
	Eigen::MatrixXf* h = new Eigen::MatrixXf[1];
	Eigen::MatrixXf* h_tilde = new Eigen::MatrixXf[1];

	// Total number of constraints
	int l;

	int T; //the time horizon for the controller
	float delta_t; //the time elapsed between each step
	int num_iterations; //the number of newton steps to perform at each timestep
	float tol_sq;
	float alpha;
	float beta;
	float mu_10, mu_11;
	float mu_12, mu_12_initial; //weighting coefficient for hard constraints AKA kappa
	float mu_13, mu_13_initial; //weighting coefficient for soft constraints AKA rho
	float gamma; //scaling factor for mu_13
	float rObs_saturation;
	bool centering_soft_constraints;
	//experiment number
	int exp_count;
	float collisionAvoidanceSoftConstraintOffset;

	float bc_epsilon; 

	int nu_X;

	Eigen::MatrixXf A_u1, B_u1;

	// velocity policy
	float mu_14, mu_15, mu_16, mu_17, mu_18, mu_19;
	float mu_20;

	// heading angle soft constraint offset
	float mu_24;
	
	// discount factor
	float beta_pl = 1.0;
	int discount_form; // 0: linear discount, 1: exponential discount
	float mu_21; // discount factor rate of change (linear)
	float mu_22; // discount factor (exponential)
	float mu_23; // discount factor (hyperbolic)

	float* Qf; //the state costs
	// Eigen::MatrixXf eig_Qf;
	Eigen::MatrixXf eig_R_rf_og;
	Eigen::MatrixXf eig_R_rf;
	Eigen::MatrixXf eig_qf;
	Eigen::MatrixXf eig_qf_tilde;

	bool preset_path = 0; // switch between A* online path planner (0) or offline preset path (1)
	int preset_path_length = 0; // length of preset path text file
	string preset_path_filename;

	float goal_tolerance;

	int line_search_style;

	int nu_sk = 0;

};

// struct Goal_Params {
// 	Goal_Params(){}
// 	//This constructor is specific to the quadrotor and uses \sout{10} 12 states
// 	Goal_Params(struct System system, bool isQuadrotor){
// 		if(isQuadrotor){
// 			float CGR[196] = { 0 };
// 			//eig_CGR = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> (CGR, 14, 14);
// 			//eig_AGR = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>(new float[392], 28, 14);
// 			// eig_wr = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, 1, Eigen::ColMajor>>(new float[28](), 28);
// 			// eig_xs = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, 1, Eigen::ColMajor>>(new float[14](), 14);
// 			// for(int i = 0; i < 14; ++i){
// 			// 	eig_AGR.row(i) = -system.eig_A.row(i);
// 			// 	eig_AGR(i,i) += 1;
// 			// 	eig_AGR.row(i+14) = eig_CGR.row(i);
// 			// }
// 			// eig_pinv = eig_AGR.completeOrthogonalDecomposition().pseudoInverse();
// 		}
// 		else throw;
// 	}

// 	// Eigen::MatrixXf eig_CGR;
// 	// Eigen::MatrixXf eig_AGR;
// 	// Eigen::VectorXf eig_wr;
// 	// Eigen::VectorXf eig_xs;
// 	// Eigen::MatrixXf eig_pinv;
// };

// A.K.A bag_of_many_things
struct Matrix_Set { //a collection of pre-allocated matrices to speed things up
	Matrix_Set(){}
	Matrix_Set(struct System system, struct MPC_Params rules){

		Tmn = rules.T * (system.m + system.n);
		nT = system.n * rules.T;

		_T = rules.T;
		_n = system.n;
		_m = system.m;

		//Note: Call for making a dynamic eigen matrix from an array ARR of N rows and M columns with elements of type TYPE

		n_by_n = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.n*system.n], system.n, system.n);
		m_by_m = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.m*system.m], system.m, system.m);

		b = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.n * rules.T](), system.n, rules.T);

		w = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, 1>>(new float[system.n], system.n);

		dnu = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.n * rules.T], system.n, rules.T);
		nu = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.n * rules.T](), system.n, rules.T);
		nu_prev = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.n * rules.T](), system.n, rules.T);
		rp = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.n * rules.T], system.n, rules.T);
		newnu = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.n * rules.T], system.n, rules.T);

		beta = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.n * rules.T], system.n, rules.T);

		rho_i_x_array = new float[system.num_soft_x];
		rho_i_x = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(rho_i_x_array, system.n, 1);
		rho_i_u_array = new float[system.num_soft_u];
		rho_i_u = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(rho_i_u_array, system.m, 1);

		FxX_array = new float[system.num_hard_x*rules.T];
		FxX = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(FxX_array, system.num_hard_x, rules.T);
		FuU = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.num_hard_u*rules.T], system.num_hard_u, rules.T);

		dx_array = new float[system.num_hard_x*rules.T];
		dx = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(dx_array, system.num_hard_x, rules.T);
		du = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.num_hard_u*rules.T], system.num_hard_u, rules.T);

		dz_x = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.n*rules.T], system.n, rules.T);
		dz_u = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.m*rules.T], system.m, rules.T);
		Ctnu_x = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.n*rules.T], system.n, rules.T);
		Ctnu_u = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.m*rules.T], system.m, rules.T);

		rd_tilde_x = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.n*rules.T], system.n, rules.T);
		rd_tilde_u = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.m*rules.T], system.m, rules.T);

		FxtildeX_array = new float[system.num_soft_x*rules.T];
		FxtildeX = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(FxtildeX_array, system.num_soft_x, rules.T);
		FutildeU = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.num_soft_u*rules.T], system.num_soft_u, rules.T);

		FxTdx = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.n*rules.T], system.n, rules.T);
		FuTdu = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.m*rules.T], system.m, rules.T);

		d_tilde_x_array = new float[system.num_soft_x*rules.T];
		d_tilde_x = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(d_tilde_x_array, system.num_soft_x, rules.T);
		d_tilde_u = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.num_soft_u*rules.T], system.num_soft_u, rules.T);

		d_hat_x_array = new float[system.num_soft_x*rules.T];
		d_hat_x = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(d_hat_x_array, system.num_soft_x, rules.T);
		d_hat_u = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.num_soft_u*rules.T], system.num_soft_u, rules.T);

		TwoQX = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.n*rules.T], system.n, rules.T);
		TwoRU = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.m*rules.T], system.m, rules.T);
		
		PtTdt_x = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.m*rules.T], system.n, rules.T);///
		PtTdt_u = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.m*rules.T], system.m, rules.T);///

		AX = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.n*rules.T], system.n, rules.T);
		BU = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.n*rules.T], system.n, rules.T);

		xdz_x = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.n*rules.T], system.n, rules.T);///
		udz_u = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.m*rules.T], system.m, rules.T);///

		newx = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.m*rules.T], system.n, rules.T);///
		newu = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.m*rules.T], system.m, rules.T);///

		set_QInv = new Eigen::MatrixXf[rules.T];
		set_RInv =  new Eigen::MatrixXf[rules.T];

		set_PTdP_x = new Eigen::MatrixXf[rules.T];
		set_PTdP_u = new Eigen::MatrixXf[rules.T];
		set_PtTdhPt_x = new Eigen::MatrixXf[rules.T];
		set_PtTdhPt_u = new Eigen::MatrixXf[rules.T];
		set_Phi_x = new Eigen::MatrixXf[rules.T];
		set_Phi_u = new Eigen::MatrixXf[rules.T];

		set_BRInv = new Eigen::MatrixXf[rules.T];
		set_AQInv = new Eigen::MatrixXf[rules.T];

		set_Y_onDiag = new Eigen::MatrixXf[rules.T];

		cholYblox_onDiag = new Eigen::MatrixXf[rules.T];

		X_goal = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.n*1], system.n, 1);
		
		X_mod = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.n*rules.T], system.n, rules.T);
		U_mod = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.m*rules.T], system.m, rules.T);

		TwoQXf = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>(new float[system.n], system.n, 1);

	}
	
	//TODO
	~Matrix_Set(){
		for(int i = 0; i < _T; i++){
			set_QInv[i].resize(0,0);
			set_RInv[i].resize(0,0);
			set_PTdP_x[i].resize(0,0);
			set_PTdP_u[i].resize(0,0);
			set_PtTdhPt_x[i].resize(0,0);
			set_PtTdhPt_u[i].resize(0,0);
			set_Phi_x[i].resize(0,0);
			set_Phi_u[i].resize(0,0);
			set_BRInv[i].resize(0,0);
			set_AQInv[i].resize(0,0);
			set_Y_onDiag[i].resize(0,0);
			cholYblox_onDiag[i].resize(0,0);
		}
	}


	Eigen::MatrixXf n_by_n;
	Eigen::MatrixXf m_by_m;

	Eigen::MatrixXf b;// = new double[T*n]();
	Eigen::MatrixXf w;// = new double[n];
	Eigen::MatrixXf dnu;// = new double[T*n];
	Eigen::MatrixXf dz_x;// = new double[nz]; //T*n
	Eigen::MatrixXf dz_u;// = new double[nz]; //T*m
	Eigen::MatrixXf nu;// = new double[T*n]();.
	Eigen::MatrixXf nu_prev;
	Eigen::MatrixXf Ctnu_x;// = new double[nz];
	Eigen::MatrixXf Ctnu_u;// = new double[nz];
	Eigen::MatrixXf rp;// = new double[T*n];
	Eigen::MatrixXf rd_tilde_x;// = new double[nz];
	Eigen::MatrixXf rd_tilde_u;// = new double[nz];
	Eigen::MatrixXf newnu;// = new double[T*n];
	Eigen::MatrixXf Pz;// = new double[Tl];
//	Eigen::MatrixXf g;// = params.g;
	Eigen::MatrixXf eX0;
	Eigen::MatrixXf beta;// = new double[n*T];

	//make some pointers so the corresponding matrices can be more easily resized
	float* FxtildeX_array;
	float* d_tilde_x_array;
	float* d_hat_x_array;
	float* FxX_array;
	float* rho_i_x_array;
	float* rho_i_u_array;
	float* dx_array;

	Eigen::MatrixXf FxX;
	Eigen::MatrixXf FuU;
	Eigen::MatrixXf FxtildeX;
	Eigen::MatrixXf FutildeU;
	Eigen::MatrixXf FxTdx;
	Eigen::MatrixXf FuTdu;

	Eigen::MatrixXf rho_i_x;
	Eigen::VectorXf rho_i_x_vec;
	Eigen::MatrixXf rho_i_u;

	Eigen::MatrixXf d_tilde_x;
	Eigen::MatrixXf d_tilde_u;
	Eigen::MatrixXf d_hat_x;
	Eigen::MatrixXf d_hat_u;
	Eigen::MatrixXf dx;
	Eigen::MatrixXf du;
	Eigen::MatrixXf TwoRU;
	Eigen::MatrixXf TwoQX;
	Eigen::MatrixXf PtTdt_x;
	Eigen::MatrixXf PtTdt_u;

	Eigen::MatrixXf AX;
	Eigen::MatrixXf BU;
//	Eigen::MatrixXf Xoff; //the X-state vector, but offset by one

	Eigen::MatrixXf xdz_x;
	Eigen::MatrixXf udz_u;

	Eigen::MatrixXf newx;
	Eigen::MatrixXf newu;

	Eigen::MatrixXf* set_PTdP_x; //Array containing the T instances of P^T d(t)^2 P for state
	Eigen::MatrixXf* set_PTdP_u; //Array containing the T instances of P^T d(t)^2 P for control
	Eigen::MatrixXf* set_PtTdhPt_x; //Array containing the T instances of P_tilde^T d_hat(t)^2 P_tilde for state
	Eigen::MatrixXf* set_PtTdhPt_u; //Array containing the T instances of P_tilde^T d_hat(t)^2 P_tilde for control

	Eigen::MatrixXf* set_Phi_x;
	Eigen::MatrixXf* set_Phi_u;

	Eigen::MatrixXf* set_QInv; //The set of matrices {Qinv1, Qinv2, Qinv3, ...}
	Eigen::MatrixXf* set_RInv; //The set of matrices {Rinv0, Rinv1, Rinv2, ...}

	Eigen::MatrixXf* set_BRInv;
	Eigen::MatrixXf* set_AQInv;

	Eigen::MatrixXf* set_BRInvBT;
	Eigen::MatrixXf* set_AQInvAT;

	Eigen::MatrixXf* set_Y_onDiag; //{Y_11, Y_22, ..., Y_TT}

	Eigen::MatrixXf* cholYblox_onDiag;

	Eigen::MatrixXf nu_proposed, x_proposed, u_proposed;

	int Tmn;
	int nT;
	float res_p, res_dx, res_du;

	int _T, _n, _m, v = 0;
	bool initial_run;

	Eigen::MatrixXf X_goal;
	Eigen::MatrixXf U_goal;
	Eigen::MatrixXf X_mod;
	Eigen::MatrixXf U_mod;

	Eigen::MatrixXf TwoQXf;
};

struct Line_Search_Params
{

	// Dichotomous line search parameters
	float dichoto_a1;		// left boundary of interval of uncertainty
	float dichoto_b1;		// right boundary of interval of uncertainty 
	float dichoto_l;		// minimum size of interval of uncertainty
	float dichoto_epsilon; 	// distinguishability constant

};

#endif

