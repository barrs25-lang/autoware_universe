// Simulation environment for the DARPA Tactical Mapping Project
// Author: Julius Allen Marshall
// Date Created: December 6th, 2021
// Last Modified: March 16th, 2022
// Contact: mjulius@vt.edu

// File Decscription ####################################################################################
// This source file is responsible for executing the model predictive control
// based trajectory planner.
// End File Decscription ################################################################################

#define DISABLE_COMMUNICATION
// List include files
#include <f_mpc_uncut.h>
#include <line_search.h>
#include "dnudz.cpp"
#include "gfgphp.cpp"
#include "rdrp.cpp"
#include "resdresp.cpp"

#ifndef DISABLE_COMMUNICATION
#include <my_client.h>
#endif
// Integers capturing thread status
int trajectory_status = 0, communication_status = 0;

// Boolean indicating if it is time to stop the program
bool time_to_exit = 0;

// LOGGER objects for logging messages, data, and runtimes.
LOGGER log_f_mpc, log_f_mpc_data, log_f_mpc_runtime, log_f_mpc_goal_data, log_f_mpc_policy_data;

// Constrant integers capturing the map's dimensions
const int grid_x = 100;
const int grid_y = 30;
const int grid_z = 100;

// Local method
// get_sign: this function gets the sign of the argument d and returns s, which is +/- 1 unless d is identically 0
// INPUTS: two template parameters, the address of the variable to store the sign of d, and d, whose parity is in question
template <class T> T get_sign(T& s, T d)
{

	if (d == 0)
	{
		
		//
		s = 0;

	} // if (d == 0)
	else
	{

		// Check the sign of d using signbit, which returns true if d is negative, 0 otherwise
		if (signbit(d))
		{

			s = -1;

		} // if (signbit(d))
		else
		{

			s = 1;

		} // if (signbit(d))

	} // if (d == 0)

	return s;

} // template <class T> T get_sign(T& s, T d)

uint64_t get_time_u()
{

	static struct timeval _time_stamp;
	gettimeofday(&_time_stamp, NULL);
	return _time_stamp.tv_sec*1000000 + _time_stamp.tv_usec;

} // uint64_t get_time_u()

// Constructor
// Reads parameter files, sets up matrices for solver, sets up log files
F_MPC_UNCUT::F_MPC_UNCUT()
{

	// Open log files for messages, planner data, and runtimes
	log_f_mpc.openLogFile("TRAJECTORY_PLANNER");
	log_f_mpc_data.openLogFile("TRAJECTORY_PLANNER_DATA");
	log_f_mpc_policy_data.openLogFile("TRAJECTORY_PLANNER_POLICY_DATA");
	log_f_mpc_goal_data.openLogFile("TRAJECTORY_PLANNER_GOAL_DATA");
	log_f_mpc_runtime.openLogFile("TRAJECTORY_PLANNER_RUNTIME");

	// Additional weighting on the control penalty
	float _weight = 8.5;

	// Pointers to floats, where arrays can be defined
	float *A_in, *B_in, *Fu_hard_in, *fu_hard_in, *Fu_soft_in, *fu_soft_in, heading_soft_offset_in;

	// Integers storing integer parameters
	int n_in, m_in, num_fu_hard_in, num_fu_soft_in, noise_source_in;

	// Floats storing float parameters: maximum roll angle, maximum pitch angle, minimum thrust, maximum thrust
	float phi_max_in, theta_max_in, T_min_in, T_max_in;

	// Floats capturing various parameters
	float noise_w_in, z_hard_ceiling_in, z_soft_ceiling_in, z_hard_floor_in, z_soft_floor_in, centering_soft_constraints_in;
	float mu_20_in, mu_24_in, bc_epsilon_in, goal_tolerance_in;

	// Eigen::MatrixXf capturing the state-transition matrix and control effectiveness matrix
	Eigen::MatrixXf eig_A_temp, eig_B_temp;

	// String capturing the parameter file's name	
	string param_filename;

	// Set the parameter file's name
	param_filename = "Parameter_Files/System_params_UGV.txt";

	try
	{
		//open a parameter file to read in
		ifstream file(param_filename);
		
		// String to store file lines
		string file_line;
		stringstream ss;
		
		//read system physical dimensions
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> quadrotor.width;
		ss >> quadrotor.length;
		ss >> quadrotor.height;

		// read quadrotor mass
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> quadrotor.mass;

		// read quadrotor length of moment arm
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> quadrotor.moment_arm_l;

		// read quadrotor coefficient of drag (props)
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> quadrotor.cT;

		// Define the size of the inertia matrix and its inverse
		quadrotor.inertia_matrix = Eigen::MatrixXf::Zero(3,3);
		quadrotor.inertia_matrix_inv = Eigen::MatrixXf::Zero(3,3);

		// read the inertia matrix
		// Iterate over the number of rows
		for(unsigned short int i = 0; i < 3; ++i)
		{
			//grab the next line of the file
			ss.clear(); getline(file, file_line); ss.str(file_line);
			
			//ignore comment lines
			if(file_line.at(0) == '/' && file_line.at(1) == '/')
			{
				
				--i;
				continue;

			} // if(file_line.at(0) == '/' && file_line.at(1) == '/')

			// Iterate over the number of columns
			for(unsigned short int j = 0; j < 3; ++j)
			{
				
				ss >> quadrotor.inertia_matrix(i,j);

			} // for(int j = 0; j < 3; ++j)

		} // for(int i = 0; i < 3; ++i)

		// Compute the inverse of the inertia matrix
		quadrotor.inertia_matrix_inv = quadrotor.inertia_matrix.inverse();

		// read the half-field-of-view angle (degrees)
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> quadrotor.half_camera_FOV;

		// convert the half-field-of-view angle (radians)
		quadrotor.half_camera_FOV = (quadrotor.half_camera_FOV*M_PI)/180;

		//read the number of states
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> n_in;

		//read the number of controls
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> m_in;

		// ************************************ //
		// read the A (state-transition) matrix //
		// ************************************ //

		// Define the size of the A matrix
		eig_A_temp = Eigen::MatrixXf::Zero(n_in,n_in);

		// Iterate over the number of states
		for(unsigned short int i = 0; i < n_in; ++i)
		{
			//grab the next line of the file
			ss.clear(); getline(file, file_line); ss.str(file_line);
			
			//ignore comment lines
			if(file_line.at(0) == '/' && file_line.at(1) == '/')
			{
				
				--i;
				continue;

			} // if(file_line.at(0) == '/' && file_line.at(1) == '/')

			// Iterate over the number of states
			for(int j = 0; j < n_in; ++j)
			{
				
				ss >> eig_A_temp(i,j);

			} // for(int j = 0; j < n_in; ++j)

		} // for(int i = 0; i < n_in; ++i)

		// **************************************** //
		// read the B (control effectivness) matrix //
		// **************************************** //
		
		// Define the size of the B matrix
		eig_B_temp = Eigen::MatrixXf::Zero(n_in,m_in);

		// Iterate over the number of states
		for(unsigned short int i = 0; i < n_in; ++i)
		{
			
			//grab the next line of the file
			ss.clear(); getline(file, file_line); ss.str(file_line);
			
			//ignore comment lines
			if(file_line.at(0) == '/' && file_line.at(1) == '/')
			{
				
				--i;
				continue;

			} // if(file_line.at(0) == '/' && file_line.at(1) == '/')

			// Iterate over the number of controls
			for(int j = 0; j < m_in; ++j)
			{
				
				ss >> eig_B_temp(i,j);

			} // for(int j = 0; j < m_in; ++j)

		} // for(int i = 0; i < n_in; ++i)

		//read noise element source
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> noise_source_in;

		//read noise element w_bar
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> noise_w_in;

		// Heading angle cone offset for heading soft constraints (in degrees): \mu_11 >= 0 (Equation 52)
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> mu_24_in; // heading_soft_offset_in
 	
 		// Convert to radians
 		mu_24_in = (mu_24_in*M_PI)/180;

		//read number of control constraints
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> num_fu_hard_in;

		// *************************************** //
		// read the hard control constraint matrix //
		// *************************************** //
		
		// Define a new float array
		Fu_hard_in = new float[num_fu_hard_in*m_in];

		// Iterate over the number of control constraints
		for(unsigned short int i = 0; i < num_fu_hard_in; ++i)
		{

			//grab the next line of the file
			ss.clear(); getline(file, file_line); ss.str(file_line);

			//ignore comment lines
			if(file_line.at(0) == '/' && file_line.at(1) == '/')
			{
				
				--i;
				continue;

			} // if(file_line.at(0) == '/' && file_line.at(1) == '/')

			// Iterate over the number of controls
			for(int j = 0; j < m_in; ++j)
			{
				
				// Read the elements of the constraint matrix
				ss >> Fu_hard_in[i*m_in+j];

			} // for(int j = 0; j < m_in; ++j)

		} // for(int i = 0; i < num_fu_hard_in; ++i)

		// ******************************* //
		// read hard control bounds vector //
		// ******************************* //

		// Define a new float array
		fu_hard_in = new float[num_fu_hard_in];

		// Iterate over the number of control constraints
		for(unsigned short int i = 0; i < num_fu_hard_in; ++i)
		{
			
			//grab the next line of the file
			ss.clear(); getline(file, file_line); ss.str(file_line);
			//ignore comment lines
			if(file_line.at(0) == '/' && file_line.at(1) == '/')
			{
				
				--i;
				continue;

			} // if(file_line.at(0) == '/' && file_line.at(1) == '/')
			
			// Read the rhs of the constraint
			ss >> fu_hard_in[i];

		} // for(int i = 0; i < num_fu_hard_in; ++i)


		// read number of soft control constraints
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> num_fu_soft_in;

		// **************************** //
		// read the soft control matrix //
		// **************************** //

		// Define a new float array 
		Fu_soft_in = new float[num_fu_soft_in*m_in];

		// Iterate over the number of soft control constraints
		for(unsigned short int i = 0; i < num_fu_soft_in; ++i)
		{
			//grab the next line of the file
			ss.clear(); getline(file, file_line); ss.str(file_line);
			//ignore comment lines
			if(file_line.at(0) == '/' && file_line.at(1) == '/'){
				--i;
				continue;
			}

			// Iterate over the number of controls
			for(int j = 0; j < m_in; ++j)
			{
				ss >> Fu_soft_in[i*m_in+j];
			}

		} // for(int i = 0; i < num_fu_soft_in; ++i)

		// ******************************* //
		// read soft control bounds vector //
		// ******************************* //
		
		// Define new float array
		fu_soft_in = new float[num_fu_soft_in];

		// Iterate over the number of soft control constraints
		for(unsigned short int i = 0; i < num_fu_soft_in; ++i)
		{
			//grab the next line of the file
			ss.clear(); getline(file, file_line); ss.str(file_line);
			
			//ignore comment lines
			if(file_line.at(0) == '/' && file_line.at(1) == '/')
			{

				--i;
				continue;

			} // if(file_line.at(0) == '/' && file_line.at(1) == '/')

			// Read the RHS of the soft control constraint
			ss >> fu_soft_in[i];

		} // for(int i = 0; i < num_fu_soft_in; ++i)

		// Read the altitude constraint
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> z_hard_ceiling_in;

		// Read the altitude soft constraint
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> z_soft_ceiling_in;

		// Read the floor constraint
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> z_hard_floor_in;

		// Read the floor soft constraint
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> z_soft_floor_in;		

		// Read the roll angle constraint
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> phi_max_in;		

		// Read the pitch angle constraint
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> theta_max_in;

		// Read the pitch angle constraint
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> T_min_in;

		// Read the pitch angle constraint
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> T_max_in;				

		// Read the boundary condition tolerance constraint
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> bc_epsilon_in;				

		// Read the flag indicating whether or not to center soft constraints
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> centering_soft_constraints_in;

		// Read the goal tolerance
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> goal_tolerance_in;

		// Close the parameter file
		file.close();

	} // try
	catch (...)
	{

		cout << "Error during System_params.txt reading" << endl;
		time_t datetime = time(0);
		tm* now = std::localtime(&datetime);
		string time = "(" + to_string(now->tm_year + 1900) + "-" + to_string(now->tm_mon + 1) + "-" + to_string(now->tm_mday) + " " + to_string(now->tm_hour) + ":" + to_string(now->tm_min) + ":" + to_string(now->tm_sec) + "): ";
		ofstream error_file;
		error_file.open("MPC_error_log.txt");
		error_file << '\n' << time << "An error has occurred while trying to parse " + param_filename;
		cin.get();
		throw;

	} // catch (...)

	// Save parameters to System (quadrotor) object
	quadrotor.n = n_in;
	quadrotor.m = m_in;
	quadrotor.num_hard_u = num_fu_hard_in;
	quadrotor.num_soft_u = num_fu_soft_in;
	quadrotor.eig_z_hard_bounds(0) = z_hard_ceiling_in;
	quadrotor.eig_z_hard_bounds(1) = -z_hard_floor_in;
	quadrotor.eig_z_soft_bounds(0) = z_soft_ceiling_in;
	quadrotor.eig_z_soft_bounds(1) = -z_soft_floor_in;
	quadrotor.phi_max = phi_max_in;
	quadrotor.theta_max = theta_max_in;
	quadrotor.T_min = T_min_in;
	quadrotor.T_max = T_max_in;

	// Define float arrays for storing input data
	float *tilde_R_r_in = new float[n_in*n_in], *R_r_f_in = new float[n_in*n_in], *R_lambda_in = new float[m_in*m_in], *tilde_R_r_lambda_in = new float[n_in*m_in]; 
	float *tilde_q_r_in = new float[n_in-2], *tilde_q_lambda_in = new float[m_in], *q_psi_in = new float[1];

	// Declare floats storing scalar input data
	float delta_T_in, tol_in, alpha_in, beta_in, kappa_in, rho_in, gamma_in, rObs_saturation_in;
	float mu_10_in, mu_11_in, mu_12_in, mu_13_in;
	float collisionAvoidanceSoftConstraintOffset_in;
	float mu_14_in, mu_15_in, mu_16_in, mu_17_in, mu_18_in, mu_19_in, mu_21_in, mu_22_in, mu_23_in;

	// Declare ints for storing scalar input data
	int T_in, num_iterations_in, discount_form_in, nu_X_in;

	// string count_string = to_string(mpc_params.exp_count);
	//param_filename = "Parameter_Files/MPC_params"+ count_string + ".txt";

	// Set the parameter file's name
	param_filename = "Parameter_Files/MPC_params.txt";

	// try opening and reading the file
	try{

		//open a parameter file to read in
		ifstream file(param_filename);

		// Declare strings to store lines of the file
		string file_line;
		stringstream ss;
		
		//read the tilde_R_r matrix
		for(unsigned short int i = 0; i < n_in; ++i)
		{

			//grab the next line of the file
			ss.clear(); getline(file, file_line); ss.str(file_line);
			
			//ignore comment lines
			if(file_line.at(0) == '/' && file_line.at(1) == '/')
			{
				
				--i;
				continue;

			} // if(file_line.at(0) == '/' && file_line.at(1) == '/')

			//parse the line for n-many floats for the Q matrix
			for(unsigned short int j = 0; j < n_in; ++j)
			{

				// Read in the \tilde_R_r matrix
				ss >> tilde_R_r_in[i*n_in+j];

			} // for(int j = 0; j < n_in; ++j)

		} // for(int i = 0; i < n_in; ++i){	

		//read the R_r,f matrix
		for(unsigned short int i = 0; i < n_in; ++i)
		{
			//grab the next line of the file
			ss.clear(); getline(file, file_line); ss.str(file_line);
			
			//ignore comment lines
			if(file_line.at(0) == '/' && file_line.at(1) == '/')
			{
				
				--i;
				continue;

			} // if(file_line.at(0) == '/' && file_line.at(1) == '/')
			
			//parse the line for n-many floats for the Q matrix
			for(unsigned short int j = 0; j < n_in; ++j)
			{
				
				// Read in the \_R_r_f matrix
				ss >> R_r_f_in[i*n_in+j];

			} // for(int j = 0; j < n_in; ++j)

		} // for(int i = 0; i < n_in; ++i)

		//read the R_lambda matrix
		for(unsigned short int i = 0; i < m_in; ++i)
		{

			//grab the next line of the file
			ss.clear(); getline(file, file_line); ss.str(file_line);

			//ignore comment lines
			if(file_line.at(0) == '/' && file_line.at(1) == '/')
			{
				--i;
				continue;
			} // if(file_line.at(0) == '/' && file_line.at(1) == '/')

			//parse the line for n-many floats for the Q matrix
			for(int j = 0; j < m_in; ++j)
			{
				
				ss >> R_lambda_in[i*m_in+j];

			} // for(int j = 0; j < m_in; ++j)

		} // for(int i = 0; i < m_in; ++i)

		//read the tilde_R_{r,\lambda} matrix
		for(unsigned short int i = 0; i < n_in; ++i)
		{

			//grab the next line of the file
			ss.clear(); getline(file, file_line); ss.str(file_line);

			//ignore comment lines
			if(file_line.at(0) == '/' && file_line.at(1) == '/')
			{
				
				--i;
				continue;

			} // if(file_line.at(0) == '/' && file_line.at(1) == '/')

			//parse the line for n-many floats for the Q matrix
			for(unsigned short int j = 0; j < m_in; ++j)
			{
				
				ss >> tilde_R_r_lambda_in[i*m_in+j];

			} // for(int j = 0; j < m_in; ++j)

		} // for(int i = 0; i < n_in; ++i)

		//read the tilde_q_r
		for(unsigned short int i = 0; i < n_in-2; ++i)
		{

			//grab the next line of the file
			ss.clear(); getline(file, file_line); ss.str(file_line);

			//ignore comment lines
			if(file_line.at(0) == '/' && file_line.at(1) == '/')
			{
				
				--i;
				continue;

			} // if(file_line.at(0) == '/' && file_line.at(1) == '/')

			ss >> tilde_q_r_in[i];

		} // for(int i = 0; i < n_in-2; ++i)

		// Read the weighting on heading angle in linear part of cost function
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> q_psi_in[0];		

		//read the tilde_q_lambda
		for(unsigned short int i = 0; i < m_in; ++i)
		{
			//grab the next line of the file
			ss.clear(); getline(file, file_line); ss.str(file_line);
			//ignore comment lines
			if(file_line.at(0) == '/' && file_line.at(1) == '/')
			{
				
				--i;
				continue;

			} // if(file_line.at(0) == '/' && file_line.at(1) == '/')

			ss >> tilde_q_lambda_in[i];

		} // for(int i = 0; i < m_in; ++i)

		//Read T, number of time steps
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> T_in;

		//Read delta_T, length of time step
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> delta_T_in;

		//Read rObs_saturation, distance threshold for obstacle influence
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> rObs_saturation_in;

		//Read num_interations, number of MPC iterations
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> num_iterations_in;

		//Read tol, residual tolerance
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> tol_in;

		// Read alpha, scaling parameter for MPC algorithm 
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> alpha_in;

		// Read beta, scaling parameter for \kappa (lagrange multiplier)
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> beta_in;

		// Read line_search_style, dictates which line search algorithm is used.
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> mpc_params.line_search_style;

		// read \mu_12, \kappa (lagrange multiplier for log barrier approximation of hard constraints)
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> mu_12_in;

		// read \mu_13, \rho (KS-approximation of soft constraints)
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> mu_13_in;

		// read \gamma, scaling parameter for \rho
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> gamma_in;

		// read \mu_10,\mu_11; \mu_10: weighing on goal or shelter priority in cost function; \mu_11: distance threshold for obstacle influence
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> mu_10_in;
		ss >> mu_11_in;

		// read \mu_14
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> mu_14_in;

		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');		
		ss >> mu_15_in;

		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');		
		ss >> mu_16_in;

		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');		
		ss >> mu_17_in;

		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');		
		ss >> mu_18_in;

		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');		
		ss >> mu_19_in;

		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');		
		ss >> discount_form_in;

		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');		
		ss >> mu_21_in;

		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');		
		ss >> mu_22_in;

		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');		
		ss >> mu_23_in;	

		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');		
		ss >> nu_X_in; // Path stride	

		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> collisionAvoidanceSoftConstraintOffset_in;	

		//close the parameter file
		file.close();
	}
	catch (...){
		time_t datetime = time(0);
		tm* now = std::localtime(&datetime);
		string time = "(" + to_string(now->tm_year + 1900) + "-" + to_string(now->tm_mon + 1) + "-" + to_string(now->tm_mday) + " " + to_string(now->tm_hour) + ":" + to_string(now->tm_min) + ":" + to_string(now->tm_sec) + "): ";
		ofstream error_file;
		error_file.open("MPC_error_log.txt");
		error_file << '\n' << time << "An error has occurred while trying to parse " + param_filename;
		throw;
	}


	// // Use Cayley-Hamilton Theorem to compute the matrix exponential
	quadrotor.eig_A = Eigen::MatrixXf::Zero(n_in,n_in);
	// quadrotor.eig_A_complex = Eigen::MatrixXcf::Zero(n_in,n_in);
	// Eigen::MatrixXf eig_A_temp_power;
	// // Step 1a: Find eigenvalues of the matrix
	// Eigen::VectorXcf eig_A_eigenvalues = eig_A_temp.eigenvalues();
	// Eigen::VectorXf eig_A_eigenvalues_real = Eigen::VectorXf::Zero(n_in);
	// Eigen::VectorXf eig_A_eigenvalues_imag = Eigen::VectorXf::Zero(n_in);
	// Eigen::VectorXi eig_A_eigenvalues_AM = Eigen::VectorXi::Ones(n_in);

	// Eigen::VectorXf temp_norm = Eigen::VectorXf::Zero(1); 

	// eig_A_eigenvalues_real = eig_A_eigenvalues.real();
	// eig_A_eigenvalues_imag = eig_A_eigenvalues.imag();

	// std::complex<float> c = 1i;

	// // Step 1b: Determine algebraic multiplicity
	// for (short int i = 0; i < n_in; i++)
	// {
	// 	for (short int j = 0; j < n_in; j++)
	// 	{
	// 		if (i != j)
	// 		{
	// 			temp_norm(0) = (eig_A_eigenvalues_real(j)-eig_A_eigenvalues_real(i))*(eig_A_eigenvalues_real(j)-eig_A_eigenvalues_real(i)) + (eig_A_eigenvalues_imag(j)-eig_A_eigenvalues_imag(i))*(eig_A_eigenvalues_imag(j)-eig_A_eigenvalues_imag(i));
	// 			if ( sqrt( temp_norm(0) ) < 1e-10)
	// 			{
	// 				eig_A_eigenvalues_AM(j)++;
	// 			}
	// 		}
	// 	}
	// }

	// // cout << "eig_A_eigenvalues:" << endl << eig_A_eigenvalues.transpose() << endl;
	// // cout << "Algebraic multiplicity: " << endl << eig_A_eigenvalues_AM.transpose() << endl;
	// // Step 2: Setup linear system to solve for coefficients of charateristic eqn
	// 		   // exp(\lambda_i*t) = sum_{k=0}^{n-1} \alpha_k * \lambda_i^k
	// 		   // where n is the dimension of A, t is a parameter (time for example) 
	// Eigen::VectorXcf lhs_coefficient_eqn = Eigen::VectorXcf::Zero(n_in);
	// for (short int i = 0; i < n_in; i++)
	// {
	// 	lhs_coefficient_eqn(i) = exp(eig_A_eigenvalues(i)*delta_T_in);
	// }

	// Eigen::MatrixXcf rhs_coefficient_eqn = Eigen::MatrixXcf::Zero(n_in,n_in);
	// Eigen::VectorXcf eig_A_product = Eigen::VectorXcf(1);

	// for (short int i = 0; i < n_in; i++)
	// {
	// 	eig_A_product(0) = 1.0;
	// 	for (short int j = 0; j < n_in; j++)
	// 	{
	// 		if (j > 0)
	// 		{
	// 			// if (eig_A_eigenvalues_imag(i) != 0)
	// 			// {
	// 			// 	eig_A_product(0) *= exp(eig_A_eigenvalues_real(i)*delta_T_in)*(cos(eig_A_eigenvalues_imag(i)*delta_T_in) + c*sin(eig_A_eigenvalues_imag(i)*delta_T_in));
	// 			// }
	// 			// else
	// 			// {
	// 				eig_A_product(0) *= eig_A_eigenvalues(i);
	// 			// }
	// 		}
			
	// 		rhs_coefficient_eqn(i,j) = eig_A_product(0);
	// 		// cout << "rhs_coefficient_eqn(" << i << "," << j <<"): " << endl << rhs_coefficient_eqn(i,j) << endl;
	// 	}
	// }


	// Step 3: Solve for the coefficients Lambda\alpha = exp(Lambda t)
	// Eigen::VectorXcf alpha_coefficients = rhs_coefficient_eqn.colPivHouseholderQr().solve(lhs_coefficient_eqn); 
	// assert(lhs_coefficient_eqn.isApprox(rhs_coefficient_eqn*alpha_coefficients));

	// Step 4: Compute the matrix exponential exp(A*t) = sum_{k=0}^{n-1} \alpha_k * A^k
	// eig_A_temp_power = Eigen::MatrixXf::Identity(n_in,n_in);
	// for (short int i = 0; i < n_in; i++)
	// {
	// 	if (i > 0)
	// 	{
	// 		eig_A_temp_power*=eig_A_temp;
	// 	}
	// 	// quadrotor.eig_A += (alpha_coefficients(i)*eig_A_temp_power);
	// 	quadrotor.eig_A_complex += (alpha_coefficients(i)*eig_A_temp_power);
	// }

	// quadrotor.eig_A = quadrotor.eig_A_complex.real();

	// cout << "eig_A_complex: " << endl << quadrotor.eig_A_complex << endl;
	// Eigen::VectorXcf discrete_eig_A_eigenvalues = quadrotor.eig_A.eigenvalues();

	// cout << "discrete_eig_A_eigenvalues: " << endl << discrete_eig_A_eigenvalues.transpose() << endl;

	quadrotor.eig_B = Eigen::MatrixXf::Zero(n_in,m_in);
	// quadrotor.eig_B = eig_A_temp.inverse()*(quadrotor.eig_A - Eigen::MatrixXf::Identity(n_in,n_in))*eig_B_temp / mpc_params.delta_t;

	quadrotor.eig_A = eig_A_temp;
	quadrotor.eig_B = eig_B_temp;

	// ****************** //
	// Set mpc parameters //
	// ****************** //

	mpc_params.T = T_in;
	mpc_params.num_iterations = num_iterations_in;
	mpc_params.tol_sq = tol_in*tol_in;
	mpc_params.alpha = alpha_in;
	mpc_params.beta = beta_in;
	mpc_params.mu_12 = mu_12_in;
	mpc_params.mu_12_initial = mu_12_in;
	mpc_params.mu_13 = mu_13_in;
	mpc_params.mu_13_initial = mu_13_in;
	mpc_params.gamma = gamma_in;
	mpc_params.delta_t = delta_T_in;
	mpc_params.rObs_saturation = rObs_saturation_in;
	mpc_params.mu_10 = mu_10_in;
	mpc_params.mu_11 = mu_11_in; 
	mpc_params.mu_14 = mu_14_in;
	mpc_params.mu_15 = mu_15_in;
	mpc_params.mu_16 = mu_16_in;
	mpc_params.mu_17 = mu_17_in;
	mpc_params.mu_18 = mu_18_in;
	mpc_params.mu_19 = mu_19_in;
	mpc_params.mu_20 = mu_20_in;
	mpc_params.discount_form = discount_form_in;
	mpc_params.mu_21 = mu_21_in;
	mpc_params.mu_22 = mu_22_in;
	mpc_params.mu_23 = mu_23_in;
	mpc_params.nu_X = nu_X_in;
	mpc_params.collisionAvoidanceSoftConstraintOffset = collisionAvoidanceSoftConstraintOffset_in;
	mpc_params.mu_24 = mu_24_in;
	mpc_params.bc_epsilon = bc_epsilon_in;
	mpc_params.centering_soft_constraints = centering_soft_constraints_in;
	mpc_params.goal_tolerance = goal_tolerance_in;

	quadrotor.noise_source = noise_source_in;
	quadrotor.noise_w = Eigen::MatrixXf::Ones(n_in, T_in);

	// Iterate over the number of states
	for (unsigned short int i = 0; i < n_in; i++)
	{
		
		// Iterate over the number of time steps
		for (unsigned short int j = 0; j < T_in; j++)
		{
			
			// Compute noise density
			quadrotor.noise_w(i,j) *= noise_w_in;	

		} // for (int j = 0; j < T_in; j++)

	} // for (int i = 0; i < n_in; i++)

	// Declare Eigen::MatrixXf to store Fpsi matrix (orientation of constraint on heading)
	Eigen::MatrixXf dummy_Fpsi;
	dummy_Fpsi = Eigen::MatrixXf::Zero(2,14);

	// Set the constraints to match Equation 34a and 34b of Marshall et. al. 2022
	dummy_Fpsi(0,12) = 1;
	dummy_Fpsi(1,12) = -1;

	// Declare Eigen::MatrixXf to store Fz matrix (orientation of constraint on altitude)
	Eigen::MatrixXf dummy_Fz_ceiling;
	dummy_Fz_ceiling = Eigen::MatrixXf::Zero(2,14);

	// Set the constraints to capture a floor and a ceiling
	dummy_Fz_ceiling(0,2) = 1;
	dummy_Fz_ceiling(1,2) = -1;

	// ******************************************************************************** //
	// The constraints below are used to capture the boundary conditions in Equation 12 // 
	// ******************************************************************************** //

	// Set the size of eig_Fx_hard_boundary_conditions 
	quadrotor.eig_Fx_hard_boundary_conditions = Eigen::MatrixXf::Zero(2*n_in,n_in);

	// Fill the matrix with 1s on diagonal 
	quadrotor.eig_Fx_hard_boundary_conditions.block(0,0,n_in,n_in) = Eigen::MatrixXf::Identity(n_in,n_in);  

	// Fill the matrix with -1s on diagonal 
	quadrotor.eig_Fx_hard_boundary_conditions.block(n_in,0,n_in,n_in) = -1*Eigen::MatrixXf::Identity(n_in,n_in);  
	
	// Set the size of eig_x_hard_boundary_conditions 
	quadrotor.eig_x_hard_boundary_conditions = Eigen::MatrixXf::Zero(2*n_in,mpc_params.T);

	// Set the size of eig_Fx_soft_boundary_conditions 
	quadrotor.eig_Fx_soft_boundary_conditions = Eigen::MatrixXf::Zero(2*n_in,n_in);

	// Fill the matrix with 1s on diagonal 
	quadrotor.eig_Fx_soft_boundary_conditions.block(0,0,n_in,n_in) = Eigen::MatrixXf::Identity(n_in,n_in);  

	// Fill the matrix with -1s on diagonal 
	quadrotor.eig_Fx_soft_boundary_conditions.block(n_in,0,n_in,n_in) = -1*Eigen::MatrixXf::Identity(n_in,n_in);  

	// Set the size of eig_x_hard_boundary_conditions 
	quadrotor.eig_x_soft_boundary_conditions = Eigen::MatrixXf::Zero(2*n_in,mpc_params.T);

	// Set the hard constraints on heading
	quadrotor.eig_Fpsi_hard = dummy_Fpsi;

	// Set the hard constraints on altitude
	quadrotor.eig_Fz_hard_ceiling = dummy_Fz_ceiling;

	// Set the hard constraints on the control input
	quadrotor.eig_Fu_hard = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (Fu_hard_in, quadrotor.num_hard_u, quadrotor.m);

	// Set the hard constraints RHS value
	quadrotor.eig_u_hard_bounds = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, 1>> (fu_hard_in, quadrotor.num_hard_u);

	// Set the hard constraints on heading
	quadrotor.eig_Fpsi_soft = quadrotor.eig_Fpsi_hard;
	
	// Set the hard constraints on altitude
	quadrotor.eig_Fz_soft_ceiling = quadrotor.eig_Fz_hard_ceiling;
	
	// Set the hard constraints on the control input
	quadrotor.eig_Fu_soft = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (Fu_soft_in, quadrotor.num_soft_u, quadrotor.m);
	
	// Set the soft constraints RHS value
	quadrotor.eig_u_soft_bounds = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, 1>> (fu_soft_in, quadrotor.num_soft_u);
	
	// Iterate over the size of the state^2
	for (unsigned short int i = 0; i < n_in*n_in; i++)
	{
		
		tilde_R_r_in[i] *= _weight;
		R_r_f_in[i] *= _weight;

	} // for (int i = 0; i < n_in*n_in; i++)

	// Set the R_r_tilde matrix, quadratic weighting on position error (set by user)
	mpc_params.eig_R_r_tilde = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (tilde_R_r_in, quadrotor.n, quadrotor.n);

	// Set the R_r_f matrix, quadratic weighting on final position error
	mpc_params.eig_R_rf = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(R_r_f_in, quadrotor.n, quadrotor.n);

	// Set the R_r_lambda matrix, quadratic weighting on coupling 
	mpc_params.eig_R_r_lambda_tilde =  Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (tilde_R_r_lambda_in, quadrotor.n, quadrotor.m);
	
	// Set the q_r_\tilde vector, linear weighting on position error
	mpc_params.eig_q_r_tilde = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (tilde_q_r_in, quadrotor.n-2, 1);

	// Set the q_psi parameter, linear weighting on heading error
	mpc_params.q_psi = (*q_psi_in);

	// Set the R_lambda_tilde matrix, quadratic weighting on control input (set by user)
	mpc_params.eig_R_lambda_tilde = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (R_lambda_in, quadrotor.m, quadrotor.m);

	// Set the R_lambda matrix, quadratic weighting on control input (varies with time)
	mpc_params.eig_R_lambda = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (R_lambda_in, quadrotor.m, quadrotor.m);
	mpc_params.eig_R_lambda_init = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (R_lambda_in, quadrotor.m, quadrotor.m);
	
	// Set the q_lambda vector, linear weighting on control input
	mpc_params.eig_q_lambda_tilde = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (tilde_q_lambda_in, quadrotor.m, 1);

	// Define matrix to store the RHS of the Schur complement (condition for positive-definiteness of block matrix tilde{eig_R})
	Eigen::MatrixXf RHS_LMI = Eigen::MatrixXf::Zero(n_in,n_in);

	// Compute the condition
	RHS_LMI = mpc_params.eig_R_r_tilde - 2*mpc_params.eig_R_r_lambda_tilde*mpc_params.eig_R_lambda_tilde.inverse()*mpc_params.eig_R_r_lambda_tilde.transpose();

	// Compute the eigenvalues
	Eigen::VectorXcf discrete_LMI_eigenvalues = RHS_LMI.eigenvalues();

	// Iterate over the number of states
	for (unsigned short int i = 0; i < n_in; i++)
	{
		
		if (discrete_LMI_eigenvalues(i).real() <= 0)
		{
			
			cout << "Schur complement condition on the positive-definiteness of block matrix tilde{eig_R} not satisfied!" << endl;
			exit(0);

		} // if (discrete_LMI_eigenvalues(i).real() <= 0)

	} // for (short int i = 0; i < n_in; i++)

	// Set a parameter file's name
	param_filename = "Parameter_Files/Line_Search_params.txt";

	try
	{
		//open a parameter file to read in
		ifstream file(param_filename);

		// Declare string to store lines of a file
		string file_line;
		stringstream ss;

		// Read the initial left boundary of the interval of uncertainty (for dichotomous line search)
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> lsp.dichoto_a1;

		// Read the initial right boundary of the interval of uncertainty (for dichotomous line search)
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> lsp.dichoto_b1;

		// Read the minimum length of the interval of uncertainty (for dichotomous line search)
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> lsp.dichoto_l;

		// Read the distuishability constant (for dichotomous line search)
		do{ss.clear(); getline(file, file_line); ss.str(file_line);}
		while(file_line.at(0) == '/' && file_line.at(1) == '/');
		ss >> lsp.dichoto_epsilon;			

		//close the parameter file
		file.close();
	}
	catch (...){
		time_t datetime = time(0);
		tm* now = std::localtime(&datetime);
		string time = "(" + to_string(now->tm_year + 1900) + "-" + to_string(now->tm_mon + 1) + "-" + to_string(now->tm_mday) + " " + to_string(now->tm_hour) + ":" + to_string(now->tm_min) + ":" + to_string(now->tm_sec) + "): ";
		ofstream error_file;
		error_file.open("MPC_error_log.txt");
		error_file << '\n' << time << "An error has occurred while trying to parse " + param_filename;
		throw;
	}

	// // Set the size of \Gamma, stores the 
	// Gamma = Eigen::MatrixXf::Zero(3,3);

	// // Set the size of \omega, stores the angular velocity of the body reference frame
	// omega = Eigen::MatrixXf::Zero(3,1);

	// // Set the size of dEuler, stores the time derivative of the Euler angles
	// dEuler = Eigen::MatrixXf::Zero(3,1);

	// Set the size of the quadrotor's state vector (its 6dof state not the ofl state)
	quadrotor.X0 = Eigen::VectorXf::Zero(n_in);

	// Set the size of the matrix storing the previous planned path
	prev_plannedPath = Eigen::MatrixXf::Zero(69,3);

	// Set the size of the vector storing the coordinates of the closest obstacle
	r_obs = Eigen::MatrixXf::Zero(3,1);

	// // Set the size of the state-transition matrix for the 
	// // higher order dynamics of u1
	// mpc_params.A_u1 = Eigen::MatrixXf::Zero(2,2);

	// // Set the elements (this is the discrete-time matrix)
	// mpc_params.A_u1(0,0) = 1.0;
	// mpc_params.A_u1(0,1) = 1.0;
	// mpc_params.A_u1(1,1) = 1.0;

	// // Set the size of the control effectiveness matrix for the
	// // higher dynamics of u1
	// mpc_params.B_u1 = Eigen::MatrixXf::Zero(2,1);

	// // Set the elements (this is the discrete-time matrix)
	// mpc_params.B_u1(0,0) = 0.5*mpc_params.delta_t*mpc_params.delta_t;
	// mpc_params.B_u1(1,0) = mpc_params.delta_t;

	// // Set the size of \delta_k, the vector storing u1 and u1dot
	// quadrotor.delta_k = Eigen::MatrixXf::Zero(2,mpc_params.T);

	// // Set the size of u1_k storing the value of u1 over the time horizon T
	// quadrotor.u1_k = Eigen::MatrixXf::Zero(mpc_params.T,1);

	// // Set the size of u1_dot_k storing the value of u1_dot over the time horizon T
	// quadrotor.u1_dot_k = Eigen::MatrixXf::Zero(mpc_params.T,1);

	// // Set the size of T_k storing the value of the thrust over the time horizon T
	// quadrotor.T_k = Eigen::MatrixXf::Zero(4,mpc_params.T);
	// quadrotor.T_k_min = Eigen::MatrixXf::Zero(4,mpc_params.T);
	// quadrotor.T_k_max = Eigen::MatrixXf::Zero(4,mpc_params.T);

	// Iterate over the number of time steps
	// for (unsigned short int i = 0; i < mpc_params.T; i++)
	// {
		
	// 	quadrotor.u1_k(i,0) = -quadrotor.mass*9.81;
		
	// 	// Iterate over the number of props
	// 	for (unsigned short int j = 0; j < 4; j++)
	// 	{
			
	// 		quadrotor.T_k(j,i) = quadrotor.mass*9.81*0.25;

	// 	} // for (int j = 0; j < 4; j++)

	// } // for (int i = 0; i < mpc_params.T; i++)

	// // Set the size of M, the mixer
	// quadrotor.M = Eigen::MatrixXf::Zero(4,4);

	// // Set mixer coefficients
	// quadrotor.M(0,0) = 0.25;
	// quadrotor.M(1,0) = 0.25;
	// quadrotor.M(2,0) = 0.25;
	// quadrotor.M(3,0) = 0.25;

	// // quadrotor.M(0,1) = 0;
	// // quadrotor.M(1,1) = -2*0.25/(quadrotor.moment_arm_l);
	// // quadrotor.M(2,1) = 0;
	// // quadrotor.M(3,1) = 2*0.25/(quadrotor.moment_arm_l);

	// // quadrotor.M(0,2) = 2*0.25/(quadrotor.moment_arm_l);
	// // quadrotor.M(1,2) = 0;
	// // quadrotor.M(2,2) = -2*0.25/(quadrotor.moment_arm_l);
	// // quadrotor.M(3,2) = 0;

	// quadrotor.M(0,1) = -1 / ( quadrotor.moment_arm_l * 4 );
	// quadrotor.M(1,1) = 1 / ( quadrotor.moment_arm_l * 4 );
	// quadrotor.M(2,1) = 1 / ( quadrotor.moment_arm_l * 4 );
	// quadrotor.M(3,1) = -1 / ( quadrotor.moment_arm_l * 4 );

	// quadrotor.M(0,2) = 1 / ( quadrotor.moment_arm_l * 4 );
	// quadrotor.M(1,2) = -1 / ( quadrotor.moment_arm_l * 4 );
	// quadrotor.M(2,2) = 1 / ( quadrotor.moment_arm_l * 4 );
	// quadrotor.M(3,2) = -1 / ( quadrotor.moment_arm_l * 4 );	

	// quadrotor.M(0,3) = -1 / ( 4 * quadrotor.cT );
	// quadrotor.M(1,3) = -1 / ( 4 * quadrotor.cT );
	// quadrotor.M(2,3) = 1 / ( 4 * quadrotor.cT );
	// quadrotor.M(3,3) = 1 / ( 4 * quadrotor.cT );

	// // Set the size of g_barrier, capturing constraints on phi,theta,thrust
	// g_barrier = Eigen::MatrixXf::Zero(mpc_params.T,1);

	// // Set pointers to new Eigen::MatrixXfs, one for each time step
	// mpc_params.eig_R_rk = new Eigen::MatrixXf[mpc_params.T];
	// mpc_params.eig_R_r_lambda_k = new Eigen::MatrixXf[mpc_params.T];
	// mpc_params.eig_q = new Eigen::MatrixXf[mpc_params.T];
	// mpc_params.eig_q_rk = new Eigen::MatrixXf[mpc_params.T];
	// mpc_params.eig_q_lambda_k = new Eigen::MatrixXf[mpc_params.T];

	// Set the size of the goal variable (xg,yg,zg,psi_g)
	goal = Eigen::MatrixXf::Zero(4,1);

	// // Set the matrices storing the roll and pitch angles resulting from a given control policy
	// quadrotor.eig_phi = Eigen::MatrixXf::Zero(mpc_params.T,1);
	// quadrotor.eig_theta = Eigen::MatrixXf::Zero(mpc_params.T,1);

	// Initialize an object which stores many variables needed to compute a solution to the MPC problem
	Matrix_Set temp_bomt(quadrotor, mpc_params);

	// Set the FMPC object's bomt object
	bomt = temp_bomt;

	// Set the b vector
	bomt.b = Eigen::MatrixXf::Zero(quadrotor.n,mpc_params.T);

	// // Allocate memory to Eigen::Matrices, one for each step in the path stride.
	// full_trajectory = new Eigen::MatrixXf[mpc_params.nu_X];
	// prevfull_trajectory = new Eigen::MatrixXf[mpc_params.nu_X];
	// full_trajectory_interf = new Eigen::MatrixXf[mpc_params.nu_X];
	// full_policy = new Eigen::MatrixXf[mpc_params.nu_X];
	// full_attitude = new Eigen::MatrixXf[mpc_params.nu_X];
	// prev_seg_eig_X = new Eigen::MatrixXf[mpc_params.nu_X];
	// prev_seg_eig_U = new Eigen::MatrixXf[mpc_params.nu_X];

	// Set the local path to zero
	localPath = Eigen::MatrixXf::Zero(1,3);

	// Set the parent ellipsoid to zero
	parentEllipsoid = Eigen::MatrixXf::Zero(2,2);

	// Set the size of the objective function value
	Obj = Eigen::MatrixXf::Zero(1,1);

	// u_k_nuX = new Eigen::MatrixXf[mpc_params.nu_X];
	// v_k_nuX = new Eigen::MatrixXf[mpc_params.nu_X];
	// du1_k_nuX = new Eigen::MatrixXf[mpc_params.nu_X];
	// ddu1_k_nuX = new Eigen::MatrixXf[mpc_params.nu_X];
	// lambda_k_nuX = new Eigen::MatrixXf[mpc_params.nu_X];
	// zeta_k_nuX = new Eigen::MatrixXf[mpc_params.nu_X];
	// T_k_nuX = new Eigen::MatrixXf[mpc_params.nu_X];

	// u_k = Eigen::MatrixXf::Zero(4,mpc_params.T);
	// v_k = Eigen::MatrixXf::Zero(4,mpc_params.T);
	// du1_k = Eigen::MatrixXf::Zero(1,mpc_params.T);
	// ddu1_k = Eigen::MatrixXf::Zero(1,mpc_params.T);
	// lambda_k = Eigen::MatrixXf::Zero(4,mpc_params.T);
	// zeta_k = Eigen::MatrixXf::Zero(4,mpc_params.T);

	// for (int i = 0; i < mpc_params.nu_X; i++)
	// {

	// 	u_k_nuX[i] = Eigen::MatrixXf::Zero(4,mpc_params.T);
	// 	v_k_nuX[i] = Eigen::MatrixXf::Zero(4,mpc_params.T);
	// 	du1_k_nuX[i] = Eigen::MatrixXf::Zero(1,mpc_params.T);
	// 	ddu1_k_nuX[i] = Eigen::MatrixXf::Zero(1,mpc_params.T);
	// 	lambda_k_nuX[i] = Eigen::MatrixXf::Zero(4,mpc_params.T);
	// 	zeta_k_nuX[i] = Eigen::MatrixXf::Zero(4,mpc_params.T);
	// 	T_k_nuX[i] = Eigen::MatrixXf::Zero(4,mpc_params.T);

	// }


} // F_MPC_UNCUT::F_MPC_UNCUT()


// Destructor
F_MPC_UNCUT::~F_MPC_UNCUT()
{
	
} // F_MPC_UNCUT::~F_MPC_UNCUT()


// Member function of F_MPC_UNCUT
// quit_handler: this function handles program behavior when certain signals are raised
// INPUTS: integer representing the raised signal
// OUTPUTS: none
void F_MPC_UNCUT::quit_handler( int sig )
{

	// If the user pressed ctrl+c
  	if (sig == SIGINT)
	{
		
		printf("\n");
		printf("<QUIT HANDLER> TERMINATING AT USER REQUEST\n");
		string message = "<QUIT HANDLER> TERMINATING AT USER REQUEST\n";
		int writeStatus = log_f_mpc.writeToLog(message);		
		printf("\n");

	} // if (sig == SIGINT)

	// If the pogram was aborted (from an error)
  	if (sig == SIGABRT)
	{
		
		printf("\n");
		printf("<QUIT HANDLER> TERMINATING AFTER ABORT RAISED\n");
		string message = "<QUIT HANDLER> TERMINATING AFTER ABORT RAISED\n";
		int writeStatus = log_f_mpc.writeToLog(message);		
		printf("\n");

	} // if (sig == SIGABRT)

	// If a segmentation fault was thrown
  	if (sig == SIGSEGV)
	{
		
		printf("\n");
		printf("<QUIT HANDLER> TERMINATING AFTER SEGMENTATION FAULT\n");
		printf("<QUIT HANDLER> PRESS ENTER TO CONTINUE\n");
		string message = "<QUIT HANDLER> AFTER SEGMENTATION FAULT\n";
		int writeStatus = log_f_mpc.writeToLog(message);		
		message = "<QUIT HANDLER> PRESS ENTER TO CONTINUE\n";
		writeStatus = log_f_mpc.writeToLog(message);		
		printf("\n");
		cin.ignore();

	} // if (sig == SIGSEGV)

	time_to_exit = 1;
	usleep(2000000);
	exit(0);

} // void F_MPC_UNCUT::quit_handler( int sig )


void F_MPC_UNCUT::verifyDimensions() const
{
    auto fail = [](const std::string& msg) {
        throw std::runtime_error("FMPC dimension error: " + msg);
    };

    if (quadrotor.eig_A.rows() != quadrotor.n ||
        quadrotor.eig_A.cols() != quadrotor.n)
        fail("A matrix must be n x n");

    if (quadrotor.eig_B.rows() != quadrotor.n ||
        quadrotor.eig_B.cols() != quadrotor.m)
        fail("B matrix must be n x m");

    if (quadrotor.num_hard_u <= 0)
        fail("No hard control constraints defined");

    if (mpc_params.T <= 0)
        fail("Time horizon T must be positive");
}


// Member function of F_MPC_UNCUT object
// compute_Gamma: this function computes the Gamma matrix used in rotational kinematic eqn
// INPUTS: none
// OUTPUTS: none (this function modifies member variables)
// void F_MPC_UNCUT::compute_Gamma()
// {
// 	float phi = pose[12], theta = pose[13];

// 	Gamma(0,0) = 1.0; 
// 	Gamma(1,0) = 0.0; 
// 	Gamma(2,0) = 0.0; 

// 	Gamma(0,1) = sin(phi)*tan(theta); 
// 	Gamma(1,1) = cos(phi); 
// 	Gamma(2,1) = sin(phi)*(1/cos(theta)); 

// 	Gamma(0,2) = cos(phi)*tan(theta); 
// 	Gamma(1,2) = -sin(phi); 
// 	Gamma(2,2) = cos(phi)*(1/cos(theta)); 	

// } // void F_MPC_UNCUT::compute_Gamma()


// Member function of F_MPC_UNCUT object
// compute_Gamma: this function computes the Gamma matrix used in rotational kinematic eqn
// INPUTS: float representing roll angle, float representing pitch angle, float representing heading angle
// pointer to a matrix representing Gamma
// OUTPUTS: none (this function modifies values pointed at by the last argument)
// void F_MPC_UNCUT::compute_Gamma(float phi, float theta, float psi, Eigen::MatrixXf* Gamma_k)
// {

// 	(*Gamma_k)(0,0) = 1.0; 
// 	(*Gamma_k)(1,0) = 0.0; 
// 	(*Gamma_k)(2,0) = 0.0; 

// 	(*Gamma_k)(0,1) = sin(phi)*tan(theta); 
// 	(*Gamma_k)(1,1) = cos(phi); 
// 	(*Gamma_k)(2,1) = sin(phi)*(1/cos(theta)); 

// 	(*Gamma_k)(0,2) = cos(phi)*tan(theta); 
// 	(*Gamma_k)(1,2) = -sin(phi); 
// 	(*Gamma_k)(2,2) = cos(phi)*(1/cos(theta)); 	

// } // void F_MPC_UNCUT::compute_Gamma(float phi, float theta, float psi, Eigen::MatrixXf* Gamma_k)


// Member function of F_MPC_UNCUT object
// get_Omega: gets the vector storing the angular velocity
// INPUTS: none
// OUTPUTS: none
// void F_MPC_UNCUT::get_Omega()
// {
	
// 	omega(0,0) = pose[15];
// 	omega(1,0) = pose[16];
// 	omega(2,0) = pose[17];

// } // void F_MPC_UNCUT::get_Omega()


// Member function of F_MPC_UNCUT object
// compute_dEuler: computes the rotational kinematic equation
// INPUTS: none
// OUTPUTS: none
// void F_MPC_UNCUT::compute_dEuler()
// {
	
// 	dEuler = Gamma*omega;

// } // void F_MPC_UNCUT::compute_dEuler()


// Member function of F_MPC_UNCUT object
// find_closest_obstacle: determines the closest obstacle to the UAVs position
// INPUTS: pointer to a matrix representing the goal point, integer representing the segment being planned
// OUTPUTS: none
void F_MPC_UNCUT::find_closest_obstacle(Eigen::MatrixXf* temp_goal_pos, int segment_number)
{

	// Declare a vector of pairs
	vector<pair<float,int>> d2;

	// Declare a float storing the distance
	float dist_norm;

	// Define matrices to store the orientation of the focal axis, the position, and a dot product
	Eigen::MatrixXf direction = Eigen::MatrixXf::Zero(3,1);
	Eigen::MatrixXf position = Eigen::MatrixXf::Zero(3,1);
	Eigen::MatrixXf dot_prod = Eigen::MatrixXf::Zero(1,1);

	direction(0,0) = cos(quadrotor.X0(12));
	direction(1,0) = sin(quadrotor.X0(12));

	position(0,0) = quadrotor.X0(0);
	position(1,0) = quadrotor.X0(1);

	// If there are obstacles
	if (local_obs.rows() > 0)
	{

		// Iterate over the number of obstacles
		for (int i = 0; i < local_obs.rows(); i++)
		{
			
			// Compute the norm of the distance
			dist_norm = ( local_obs.row(i) - (*temp_goal_pos) ).norm(); 

			// Compute the dot product
			dot_prod = ( (local_obs.row(i).transpose() - position).transpose()*direction).diagonal();
			
			// If the distance is less than some user-defined value and the obstacle is in the same direction the UAV is facing
			if ( ( dist_norm < 1/mpc_params.mu_11 ) && ( dot_prod(0,0) > 0 ) )
			{
				
				d2.push_back( make_pair( dist_norm , i ) );

			} // if ( ( dist_norm < 1/mpc_params.mu_11 ) && ( dot_prod(0,0) > 0 ) )

		} // for (int i = 0; i < local_obs.rows(); i++)

		// If there are obstacle in the same direction that the UAV is facing
		if (d2.size() > 0)
		{

			// Sort the vector storing the distance norms and labels
			sort(d2.begin(),d2.end());

			// Compute the nearest obstacle
			r_obs = local_obs.row(d2[0].second).transpose();

		} // if (d2.size() > 0)
		else if (d2.size() == 0 && localPath.rows() > 0 && segment_number < localPath.rows())
		{

			r_obs(0,0) = localPath(segment_number,0); 
			r_obs(1,0) = localPath(segment_number,1); 
			r_obs(2,0) = localPath(segment_number,2); 

		} // else if (d2.size() == 0 && localPath.rows() > 0 && traj_iterator < localPath.rows())
		else
		{

			r_obs(0,0) = quadrotor.X0(0); 
			r_obs(1,0) = quadrotor.X0(1); 
			r_obs(2,0) = quadrotor.X0(2); 

		} // if (d2.size() > 0)

	} // if (local_obs.rows() > 0)
	else
	{

		// If there are no obstacles, just use the planned path
		if (localPath.rows() > 0)
		{
			
			r_obs(0,0) = localPath(segment_number,0); 
			r_obs(1,0) = localPath(segment_number,1); 
			r_obs(2,0) = localPath(segment_number,2); 	

		} // if (localPath.rows() > 0)
		else
		{
			
			r_obs(0,0) = quadrotor.X0(0); 
			r_obs(1,0) = quadrotor.X0(1); 
			r_obs(2,0) = quadrotor.X0(2); 

		} // if (localPath.rows() > 0)

	} // if (local_obs.rows() > 0)

	cout << "<FMPC> find_closest_obstacle complete" << endl;

} // void F_MPC_UNCUT::find_closest_obstacle(Eigen::MatrixXf* temp_goal_pos, int segment_number)


// Member function of F_MPC_UNCUT object
// fsat: computes the fsat function, whose output is in the range [0,1]
// INPUTS: pointer to a matrix storing the value to check
// OUTPUTS: float representing the output of fsat function, [0,1]
float F_MPC_UNCUT::fsat(Eigen::MatrixXf* arg)
{

	// define float representing the input's norm (either Equi-induced Euclidean norm for matrix or Euclidean norm for vector)
	float arg_norm = (*arg).norm();

	// If the norm is less than 1, set the output to 1.0
	if (arg_norm < 1)
	{
		
		return 1.0f;
	
	} // if (arg_norm < 1)
	else
	{
		
		return 1.0f/arg_norm;

	} // if (arg_norm < 1)

	cout << "<FMPC> fsat complete" << endl;

} // float F_MPC_UNCUT::fsat(Eigen::MatrixXf* arg)


// Member function of F_MPC_UNCUT object
// hat_operator: computes a skew-symmetric matrix representing the
// hat operation of the input vectors
// INPUTS: pointer to a matrix storing the value to convert, pointer to a matrix
// which will store the result, integer representing the time step
// OUTPUTS: none
void hat_operator(Eigen::MatrixXf* omega, Eigen::MatrixXf* omega_cross, int i)
{

	(*omega_cross) = Eigen::MatrixXf::Zero(3,3);
	(*omega_cross)(0,1) = -(*omega)(2,i);
	(*omega_cross)(0,2) = (*omega)(1,i);
	(*omega_cross)(1,0) = (*omega)(2,i);
	(*omega_cross)(1,2) = -(*omega)(0,i);
	(*omega_cross)(2,0) = -(*omega)(1,i);
	(*omega_cross)(2,1) = (*omega)(0,i);

	// cout << "<FMPC> hat_operator complete" << endl;

} // void hat_operator(Eigen::MatrixXf* omega, Eigen::MatrixXf* omega_cross, int i)


// Member function of F_MPC_UNCUT object
// update_weighting_matrices: updates the tilde_R matrix in the cost function
// INPUTS: integer representing the segment number
// OUTPUTS: none
void F_MPC_UNCUT::update_weighting_matrices(int segment_number, bool success)
{

	// Declare float to store update weighting
	float recast_weight;

	// Matrix storing
	Eigen::MatrixXf priority_weight;
	Eigen::MatrixXf fsat_arg, temp;

	temp = Eigen::MatrixXf::Zero(3,1);
	fsat_arg = Eigen::MatrixXf::Zero(3,1);
	priority_weight = Eigen::MatrixXf::Zero(quadrotor.n,1);
	temp(0,0) = quadrotor.X0(0);
	temp(1,0) = quadrotor.X0(1);
	temp(2,0) = quadrotor.X0(2);

	// If the planned path exists, and we are not looking past the end of the path
	if (localPath.rows() > 0 && traj_iterator > 0 && traj_iterator < localPath.rows())
	{
	
		// Compute fsat's argument
		fsat_arg = mpc_params.mu_11*( localPath.row(traj_iterator).transpose() - r_obs );
	
	}
	else
	{
		
		fsat_arg = mpc_params.mu_11*( temp - r_obs ); 

	} // if (localPath.rows() > 0 && traj_iterator < localPath.rows())

	// Compute the weighting
	recast_weight = (1 + mpc_params.mu_10*(1 - fsat( &fsat_arg ) ) );

	cout << "pos: " << temp << endl;
	cout << "r_obs: " << r_obs << endl;

	// If the planned path exists
	if (localPath.rows() > 0 && traj_iterator > 0 && traj_iterator < localPath.rows())
	{ 

		priority_weight.block(0,0,3,1) = mpc_params.mu_10 * localPath.row(traj_iterator).transpose() + (1 - mpc_params.mu_10)*fsat( &fsat_arg ) * r_obs;
	
	}
	else
	{

		priority_weight.block(0,0,3,1) = mpc_params.mu_10 * temp + (1 - mpc_params.mu_10)*fsat( &fsat_arg ) * r_obs;
	
	} // if (localPath.rows() > 0 && traj_iterator > 0)

	// Reset the matrix weighting the control input (it does not have a time dependence unlike the others)
	mpc_params.eig_R_lambda = mpc_params.eig_R_lambda_tilde; // have to reset this matrix due to the discounting

	Eigen::MatrixXf RHS_LMI = Eigen::MatrixXf::Zero(quadrotor.n,quadrotor.n);

	Eigen::VectorXcf discrete_LMI_eigenvalues;

	string message;	

	// Iterate over the number of time steps
	for (unsigned short int i = 0; i < mpc_params.T; i++)
	{
		
		// Set the weighting matrices to zero
		mpc_params.eig_R_rk[i] = Eigen::MatrixXf::Zero(quadrotor.n,quadrotor.n);
		mpc_params.eig_R_r_lambda_k[i] = Eigen::MatrixXf::Zero(quadrotor.n,quadrotor.m);
		mpc_params.eig_q_rk[i] = Eigen::MatrixXf::Zero(quadrotor.n-2,1);
		mpc_params.eig_q_lambda_k[i] = Eigen::MatrixXf::Zero(quadrotor.m,1);
		mpc_params.eig_q[i] = Eigen::MatrixXf::Zero(quadrotor.n,1);

		// Compute the weighting matrices
		if (success || segment_number == 0)
		{

			// mpc_params.eig_R_rk[i] = recast_weight*(1/g_barrier(i,0))*mpc_params.eig_R_r_tilde;
			// mpc_params.eig_R_r_lambda_k[i] = recast_weight*(1/g_barrier(i,0))*mpc_params.eig_R_r_lambda_tilde;
			mpc_params.eig_q_rk[i] = recast_weight*( mpc_params.eig_q_r_tilde - 0 * mpc_params.eig_R_r_tilde.block(0,0,12,12) * priority_weight.block(0,0,quadrotor.n-2,1));
			mpc_params.eig_q_lambda_k[i] = mpc_params.eig_q_lambda_tilde - 0 * mpc_params.eig_R_r_lambda_tilde.transpose() * priority_weight;
			// mpc_params.eig_R_lambda = (1/g_barrier(i,0))*mpc_params.eig_R_lambda_init;
			(mpc_params.eig_q[i]).block(0,0,12,1) = mpc_params.eig_q_rk[i];
			mpc_params.eig_q[i](12,0) = mpc_params.q_psi;
			mpc_params.eig_q[i](13,0) = 0.0;

			// Compute the condition
			RHS_LMI = mpc_params.eig_R_rk[i] - 2*mpc_params.eig_R_r_lambda_k[i]*mpc_params.eig_R_lambda.inverse()*mpc_params.eig_R_r_lambda_k[i].transpose();

			// Compute the eigenvalues
			discrete_LMI_eigenvalues = RHS_LMI.eigenvalues();		

			// Iterate over the number of states
			for (unsigned short int j = 0; j < quadrotor.n; j++)
			{
				
				if (discrete_LMI_eigenvalues(j).real() <= 0.0)
				{
					cout << "discrete_LMI_eigenvalues(j).real(): " << discrete_LMI_eigenvalues(j).real() << endl;
					cout << "mpc_params.eig_R_rk[i]: " << endl << mpc_params.eig_R_rk[i] << endl;
					message = "Schur complement condition time step " + to_string(i) + " state " + to_string(j) + ", resetting " + to_string(i) + "th block.";
					cout << message << endl;
					log_f_mpc.writeToLog(message);

					mpc_params.eig_R_rk[i] = mpc_params.eig_R_r_tilde;
					mpc_params.eig_R_r_lambda_k[i] = mpc_params.eig_R_r_lambda_tilde;
					(mpc_params.eig_q[i]).block(0,0,12,1) = mpc_params.eig_q_r_tilde;
					mpc_params.eig_q[i](12,0) = mpc_params.q_psi;
					mpc_params.eig_q[i](13,0) = 0.0;

					break;

				} // if (discrete_LMI_eigenvalues(i).real() <= 0)

			} // for (short int i = 0; i < n_in; i++)

		}
		else
		{
			mpc_params.eig_R_rk[i] = mpc_params.eig_R_r_tilde;
			mpc_params.eig_R_r_lambda_k[i] = mpc_params.eig_R_r_lambda_tilde;
			mpc_params.eig_q_rk[i] = 0*recast_weight*( mpc_params.eig_q_r_tilde - 0 * mpc_params.eig_R_r_tilde.block(0,0,12,12) * priority_weight.block(0,0,quadrotor.n-2,1));
			mpc_params.eig_q_lambda_k[i] = 0*mpc_params.eig_q_lambda_tilde - 0 * mpc_params.eig_R_r_lambda_tilde.transpose() * priority_weight;
			mpc_params.eig_R_lambda = mpc_params.eig_R_lambda_init;
			(mpc_params.eig_q[i]).block(0,0,12,1) = mpc_params.eig_q_rk[i];
			mpc_params.eig_q[i](12,0) = mpc_params.q_psi;
			mpc_params.eig_q[i](13,0) = 0.0;
		}

	} // for (int i = 0; i < mpc_params.T; i++)
	
	cout << "<FMPC> update_weighting_matrices complete" << endl;

} // void F_MPC_UNCUT::update_weighting_matrices(int segment_number)


// Member function of F_MPC_UNCUT object
// update_and_discount_weighting_matrices: updates the tilde_R matrix in the cost function and discounts them
// INPUTS: integer representing the segment number
// OUTPUTS: none
void F_MPC_UNCUT::update_and_discount_weighting_matrices(int segment_number, bool success)
{

	// Update the \tilde{R} matrix
	update_weighting_matrices(segment_number, success);

	// Update the discount factor
	mpc_params.beta_pl = pow(mpc_params.mu_21,segment_number);

	// Iterate over the number of time steps 
	for (unsigned short int i = 0; i < mpc_params.T; i++)
	{
		
		mpc_params.eig_R_rk[i] = mpc_params.eig_R_rk[i]*mpc_params.beta_pl;
		mpc_params.eig_R_r_lambda_k[i] = mpc_params.eig_R_r_lambda_k[i]*mpc_params.beta_pl;
		mpc_params.eig_q_rk[i] = mpc_params.eig_q_rk[i]*mpc_params.beta_pl;
		mpc_params.eig_q_lambda_k[i] = mpc_params.eig_q_lambda_k[i]*mpc_params.beta_pl;
		mpc_params.eig_q[i] = mpc_params.eig_q[i]*mpc_params.beta_pl;
		mpc_params.eig_R_lambda = mpc_params.eig_R_lambda*mpc_params.beta_pl;

	} // for (int i = 0; i < mpc_params.T; i++)

	cout << "<FMPC> update_and_discount_weighting_matrices complete" << endl;
	
} // void F_MPC_UNCUT::update_and_discount_weighting_matrices(int segment_number)

bool F_MPC_UNCUT::validates_collision_constraints()
{

	Eigen::Vector3f n0, n1, n2, p_temp, proj_goal;
	Eigen::Vector3f v, proj_quad_on_plane_temp;
	Eigen::Vector3f	goal_temp = Eigen::Map<Eigen::Vector3f>(goal.data(),3);

	Eigen::MatrixXf	p = Eigen::MatrixXf::Zero(3,collisionConstraints.rows());
	Eigen::MatrixXf proj_quad_on_plane = Eigen::MatrixXf::Zero(3,collisionConstraints.rows());

	float d = 0.0, z1 = 0.0, plane_size = 2.0;
	float dist = 0.0;

	// Iterate over the number of constraints
	for (int j = 0; j < collisionConstraints.rows(); j++)
	{
		// find normal vector to the plane & offset
		n0(0) = collisionConstraints(j,0);
		n0(1) = collisionConstraints(j,1);
		n0(2) = collisionConstraints(j,2);
		d = collisionConstraints(j,3);

		// Find two vectors that form an orthonormal basis for the plane
		n1(0) = -n0(1);
		n1(1) = n0(0);
		n1(2) = z1;
		n1.normalize();

		n2 = n0.cross(n1); 

		// Find a point on the plane
		p_temp = quadrotor.X0.head(3) + n2*plane_size + n1*plane_size + n0*d;
		p.col(j) = Eigen::Map<Eigen::MatrixXf>(p_temp.data(),3,1);

	} // for (int j = 0; j < collisionConstraints.rows(); j++)

	// Iterate over the number of constraints
	for (int j = 0; j < collisionConstraints.rows(); j++)
	{
		// Project the quadrotor point onto each constraint
		v = quadrotor.X0.head(3) - Eigen::Map<Eigen::Vector3f>(p.col(j).data(),3);
		
		if (abs(collisionConstraints.block(j,0,1,3).norm()) < 10e-6)
		{
			continue;
		}

		n0(0) = collisionConstraints(j,0);
		n0(1) = collisionConstraints(j,1);
		n0(2) = collisionConstraints(j,2);
		dist = v.dot(-n0);

		proj_quad_on_plane_temp = quadrotor.X0.head(3) + dist*n0;
		proj_quad_on_plane.col(j) = Eigen::Map<Eigen::MatrixXf>(proj_quad_on_plane_temp.data(),3,1);

		// If the vectors do not point in opposite directions, 
		// the goal is not in the constraint set
		if (n0.dot(goal_temp) >= collisionConstraints(j,3))
		{
			// indicate that the goal is not in the constraint set
			return false;

		} // if (-n.dot(difference_temp) >= 0)

	} // for (int j = 0; j < collisionConstraints.rows(); j++)

	cout << "<FMPC> validates_collision_constraints complete" << endl;

	return true;

}

// Member function of F_MPC_UNCUT object
// get_Ginv_k: computes G_inv
// INPUTS: floats representing the roll, pitch, and heading angles,
// float representing the total thrust, pointer to a matrix representing 
// Ginv_k.
// OUTPUTS: none
void F_MPC_UNCUT::get_Ginv_k(float phi, float theta, float psi, float u1, Eigen::MatrixXf* G_inv)
{

	float cp = cos(phi);
	float ct = cos(theta);
	float ci = cos(psi);
	float sp = sin(phi);
	float st = sin(theta);
	float si = sin(psi);

	float I1 = quadrotor.inertia_matrix(0,0);
	float I2 = quadrotor.inertia_matrix(1,1);
	float I3 = quadrotor.inertia_matrix(2,2);

	(*G_inv)(0,0) = quadrotor.mass*(sp*si + cp*ci*st);  
	(*G_inv)(0,1) = quadrotor.mass*(cp*si*st - ci*sp);
	(*G_inv)(0,2) = quadrotor.mass*(cp*ct);
	(*G_inv)(0,3) = 0;

	(*G_inv)(1,0) = (quadrotor.mass*I1*(cp*si-ci*sp*st)) / (u1);
	(*G_inv)(1,1) = -(quadrotor.mass*I1*(cp*ci+si*sp*st)) / (u1);
	(*G_inv)(1,2) = -(quadrotor.mass*I1*(ct*sp)) / (u1);
	(*G_inv)(1,3) = 0;

	(*G_inv)(2,0) = (quadrotor.mass*I2*(ci*ct)) / (u1);
	(*G_inv)(2,1) = (quadrotor.mass*I2*(ct*si)) / (u1);
	(*G_inv)(2,2) = -(quadrotor.mass*I2*(st)) / (u1);
	(*G_inv)(2,3) = 0;

	(*G_inv)(3,0) = -(quadrotor.mass*I3*(ci*ct*sp)) / (u1*cp);
	(*G_inv)(3,1) = -(quadrotor.mass*I3*(ct*si*sp)) / (u1*cp);
	(*G_inv)(3,2) = (quadrotor.mass*I3*(st*sp)) / (u1*cp);
	(*G_inv)(3,3) = (quadrotor.mass*I3*(ct)) / cp;

	// cout << "<FMPC> get_Ginv_k complete" << endl;

} // void F_MPC_UNCUT::get_Ginv_k(float phi, float theta, float psi, float u1, Eigen::MatrixXf* G_inv)


// Member function of F_MPC_UNCUT object
// get_f_k: computes f
// INPUTS: floats representing the roll, pitch, and heading angles,
// pointer to a matrix the angular velocity, float representing the total thrust
// pointer to a matrix representing a float, integer representing the time step
// OUTPUTS: none
void F_MPC_UNCUT::get_f_k(float phi, float theta, float psi, Eigen::MatrixXf* omega, float u1, float u1_dot, Eigen::MatrixXf* f, int i)
{

	float cp = cos(phi);
	float ct = cos(theta);
	float ci = cos(psi);
	float sp = sin(phi);
	float st = sin(theta);
	float si = sin(psi);

	float I1 = quadrotor.inertia_matrix(0,0);
	float I2 = quadrotor.inertia_matrix(1,1);
	float I3 = quadrotor.inertia_matrix(2,2);

	float mass = quadrotor.mass;

	float omega1 = (*omega)(0,i);
	float omega2 = (*omega)(1,i);
	float omega3 = (*omega)(2,i);

	(*f)(0,0) = u1_dot*(((cp*si - ci*sp*st)*(omega1*ct + omega3*cp*st + omega2*sp*st))/(mass*ct) + ((ci*sp - cp*si*st)*(omega3*cp + omega2*sp))/(mass*ct) + (cp*ci*ct*(omega2*cp - omega3*sp))/mass) + ((omega3*cp + omega2*sp)*((u1_dot*(ci*sp - cp*si*st))/mass - (u1*(sp*si + cp*ci*st)*(omega3*cp + omega2*sp))/(mass*ct) + (u1*(cp*ci + sp*si*st)*(omega1*ct + omega3*cp*st + omega2*sp*st))/(mass*ct) - (u1*cp*ct*si*(omega2*cp - omega3*sp))/mass))/ct - (ci*(omega2*cp - omega3*sp)*(omega2*u1*st - u1_dot*cp*ct + omega1*u1*ct*sp))/mass - ((omega1*ct + omega3*cp*st + omega2*sp*st)*(omega1*u1*sp*si - u1_dot*cp*si + u1_dot*ci*sp*st + omega1*u1*cp*ci*st))/(mass*ct) + (omega2*omega3*u1*(I2 - I3)*(cp*si - ci*sp*st))/(I1*mass) - (omega1*omega3*u1*ci*ct*(I1 - I3))/(I2*mass);
	(*f)(1,0) = u1_dot*(((sp*si + cp*ci*st)*(omega3*cp + omega2*sp))/(mass*ct) - ((cp*ci + sp*si*st)*(omega1*ct + omega3*cp*st + omega2*sp*st))/(mass*ct) + (cp*ct*si*(omega2*cp - omega3*sp))/mass) + ((omega3*cp + omega2*sp)*((u1_dot*(sp*si + cp*ci*st))/mass + (u1*(ci*sp - cp*si*st)*(omega3*cp + omega2*sp))/(mass*ct) + (u1*(cp*si - ci*sp*st)*(omega1*ct + omega3*cp*st + omega2*sp*st))/(mass*ct) + (u1*cp*ci*ct*(omega2*cp - omega3*sp))/mass))/ct - ((omega1*ct + omega3*cp*st + omega2*sp*st)*(u1_dot*cp*ci + u1_dot*sp*si*st - omega1*u1*ci*sp + omega1*u1*cp*si*st))/(mass*ct) - (si*(omega2*cp - omega3*sp)*(omega2*u1*st - u1_dot*cp*ct + omega1*u1*ct*sp))/mass - (omega2*omega3*u1*(I2 - I3)*(cp*ci + sp*si*st))/(I1*mass) - (omega1*omega3*u1*ct*si*(I1 - I3))/(I2*mass);
	(*f)(2,0) = (omega1*omega3*u1*st*(I1 - I3))/(I2*mass) - ((u1_dot*sp + omega1*u1*cp)*(omega1*ct + omega3*cp*st + omega2*sp*st))/mass - ((omega2*cp - omega3*sp)*(omega2*u1*ct + u1_dot*cp*st - omega1*u1*sp*st))/mass - (u1_dot*(omega2*st + omega1*ct*sp))/mass - (omega2*omega3*u1*ct*sp*(I2 - I3))/(I1*mass);
	(*f)(3,0) = ((omega2*cp - omega3*sp)*(omega1*ct + omega3*cp*st + omega2*sp*st))/(ct*ct) + (st*(omega3*cp + omega2*sp)*(omega2*cp - omega3*sp))/(ct*ct) + (omega1*omega2*cp*(I1 - I2))/(I3*ct) - (omega1*omega3*sp*(I1 - I3))/(I2*ct);

	// cout << "<FMPC> get_f_k complete" << endl;

} // void F_MPC_UNCUT::get_f_k(float phi, float theta, float psi, Eigen::MatrixXf* omega, float u1, float u1_dot, Eigen::MatrixXf* f, int i)


rot_dyn_state sys_state;

// Member function of F_MPC_UNCUT object
// sys: copies ODEs over, ready for integration by runge_kutta4.
// INPUTS: rot_dyn_state variable representing the result of integration,
// rot_dyn_state variable representing the integrand, double representing time
// OUTPUTS: none
void F_MPC_UNCUT::sys(const rot_dyn_state &y, rot_dyn_state &ydot, double t)
{
	
	// Iterate over the number of states to integrate
	for (unsigned short int i = 0; i < rot_dyn_states; i++)
	{	

		ydot[i] = sys_state[i];

	} // for (short int i = 0; i < rot_dyn_states; i++)

} // void F_MPC_UNCUT::sys(const rot_dyn_state &y, rot_dyn_state &ydot, double t)


// Member function of F_MPC_UNCUT object
// compute_roll_pitch_thrust_for_lambda_k_and_g_barrier: computes the roll, pitch, thrust, and g_barrier for a given control policy
// INPUTS: pointer to a variable storing the ofl state, pointer to a matrix storing the control policy
// OUTPUTS: boolean indicating whether or not any element of g_barrier is <= 0
// bool F_MPC_UNCUT::compute_roll_pitch_thrust_for_lambda_k_and_g_barrier(Eigen::MatrixXf* chi, Eigen::MatrixXf* lambda)
// {

// 	// Double array storing ode rhs, and time
// 	double ode_rhs[8], t = 0.0;

// 	// Define eigen matrices to store various values necessary for computing roll, pitch, thrust, and g_barrier
// 	Eigen::MatrixXf temp_rhs = Eigen::MatrixXf::Zero(3,1);
// 	Eigen::MatrixXf chi_k = Eigen::MatrixXf::Zero(14,mpc_params.T);
// 	Eigen::MatrixXf chi_dot_k = Eigen::MatrixXf::Zero(14,mpc_params.T);
// 	Eigen::MatrixXf Gamma_k = Eigen::MatrixXf::Zero(3,3);
// 	Eigen::MatrixXf omega_k = Eigen::MatrixXf::Zero(3,mpc_params.T);
// 	Eigen::MatrixXf omega_k_cross = Eigen::MatrixXf::Zero(3,3);
// 	Eigen::MatrixXf G_inv_k = Eigen::MatrixXf::Zero(4,4);
// 	Eigen::MatrixXf f_k = Eigen::MatrixXf::Zero(4,1);
// 	Eigen::MatrixXf virtual_control_k = Eigen::MatrixXf::Zero(4,mpc_params.T);
// 	// Eigen::MatrixXf zeta_k = Eigen::MatrixXf::Zero(4,mpc_params.T);
// 	// Eigen::MatrixXf u_k = Eigen::MatrixXf::Zero(4,mpc_params.T);
// 	Eigen::MatrixXf local_u1_k = Eigen::MatrixXf::Zero(mpc_params.T,1);
// 	Eigen::MatrixXf local_du1_k = Eigen::MatrixXf::Zero(mpc_params.T,1);
// 	Eigen::MatrixXf phi_k = Eigen::MatrixXf::Zero(mpc_params.T,1);
// 	Eigen::MatrixXf theta_k = Eigen::MatrixXf::Zero(mpc_params.T,1);
// 	Eigen::MatrixXf psi_k = Eigen::MatrixXf::Zero(mpc_params.T,1);

// 	// Declare a runge_kutta object
	// runge_kutta4<rot_dyn_state> rk4;

// 	// Declare a rot_dyn_state data type to store output from integration
// 	rot_dyn_state y = {0.0,0.0,0.0,0.0,0.0,0.0,9.81,0.0};

// 	if (abs(y[6]) > 100)
// 	{
// 		y[6] = 9.81;
// 	}
// 	if (abs(y[7]) > 100)
// 	{
// 		y[7] = 0.0;
// 	}

// 	// Iterate over the number of time steps
// 	for (unsigned short int i = 0; i < mpc_params.T; i++)
// 	{

// 		// If this is the first step
// 		if (i == 0)
// 		{

// 			// Use information stored in the original state vector
// 			chi_k.block(0,0,12,1) = quadrotor.X0.block(0,0,12,1);
// 			phi_k(i,0) = pose[12];
// 			theta_k(i,0) = pose[13];
// 			psi_k(i,0) = pose[14];
// 			Gamma_k = Gamma;
// 			omega_k.col(i) = omega;
// 			local_u1_k(i,0) = quadrotor.u1_k(i,0);
// 			quadrotor.u1_k(i,0) = y[6];
// 			local_du1_k(i,0) = y[7];

// 			y[0] = phi_k(i,0);
// 			y[1] = theta_k(i,0);
// 			y[2] = psi_k(i,0);
// 			y[3] = omega_k(0,i);
// 			y[4] = omega_k(1,i);
// 			y[5] = omega_k(2,i);

// 		} // if (i == 0)
// 		else
// 		{

// 			// Otherwise use information stored in the integrator output
// 			chi_k.col(i) = (*chi).col(i);
// 			phi_k(i,0) = y[0];
// 			theta_k(i,0) = y[1];
// 			// psi_k(i,0) = chi_k(12,i);
// 			psi_k(i,0) = y[2];
// 			omega_k(0,i) = y[3]; // or y[2]?
// 			omega_k(1,i) = y[4];
// 			omega_k(2,i) = y[5];
// 			local_u1_k(i,0) = y[6];
// 			quadrotor.u1_k(i,0) = y[6];
// 			local_du1_k(i,0) = y[7];
// 			compute_Gamma(phi_k(i,0),theta_k(i,0),psi_k(i,0),&Gamma_k);

// 		} // if (i == 0)

// 		// Compute the lhs of the rotational kinematics
// 		temp_rhs = Gamma_k*omega_k.col(i);

// 		// Set the rotational kinematics
// 		for (unsigned short int j = 0; j < 3; j++)
// 		{
// 			ode_rhs[j] = temp_rhs(j,0);
// 		}

// 		// If this is the first step, set chi
// 		if (i == 0)
// 		{
// 			chi_k(12,0) = psi_k(i,0);
// 			chi_k(13,0) = temp_rhs(2,0);
// 		}

// 		// Call the hat operator to compute omega_cross
// 		hat_operator(&omega_k, &omega_k_cross, i);
		
// 		// Get the G_inv matrix
// 		get_Ginv_k(phi_k(i,0),theta_k(i,0),psi_k(i,0),quadrotor.u1_k(i,0),&G_inv_k);

// 		// Get the f vector
// 		get_f_k(phi_k(i,0),theta_k(i,0),psi_k(i,0),&omega_k,quadrotor.u1_k(i,0),quadrotor.u1_dot_k(i,0),&f_k,i);

// 		// Compute the RHS of the output feedback linearized dynamical system
// 		chi_dot_k.col(i) = quadrotor.eig_A*chi_k.col(i) + quadrotor.eig_B*(*lambda).col(i);

// 		// Get the virtual control, so it is easier to compute the actual control
// 		virtual_control_k(0,i) = chi_dot_k(9,i);
// 		virtual_control_k(1,i) = chi_dot_k(10,i);
// 		virtual_control_k(2,i) = chi_dot_k(11,i);
// 		virtual_control_k(3,i) = chi_dot_k(13,i);

// 		// Compute the actual control input
// 		zeta_k.col(i) = G_inv_k*(-f_k + virtual_control_k.col(i));

// 		// Compute the rhs of the rotational dynamic equations
// 		temp_rhs = quadrotor.inertia_matrix_inv*(zeta_k.block(1,i,3,1) - omega_k_cross*quadrotor.inertia_matrix*omega_k.col(i));
		
// 		for (unsigned short int j = 3; j < 6; j++)
// 		{
// 			ode_rhs[j] = temp_rhs(j-3,0);
// 		}

// 		// Two states of the integrand are reserved for the total thrust dynamics
// 		ode_rhs[6] = local_du1_k(i,0);
// 		ode_rhs[7] = zeta_k(0,i);
		
// 		for (unsigned short int j = 0; j < rot_dyn_states; j++)
// 		{
// 			sys_state[j] = ode_rhs[j];
// 		}

// 		// Compute the "time"
// 		t = i*mpc_params.delta_t;

// 		// Do a step
// 		rk4.do_step(sys, y, t, mpc_params.delta_t);

// 		// compute the control input
// 		// u_k(0,i) = quadrotor.u1_k(i,0);
// 		u_k(0,i) = y[6];
// 		u_k(1,i) = zeta_k(1,i);
// 		u_k(2,i) = zeta_k(2,i);
// 		u_k(3,i) = zeta_k(3,i);

// 		du1_k(0,i) = local_du1_k(i,0);
// 		ddu1_k(0,i) = (*lambda)(3,i);

// 		// Compute the total thrust
// 		quadrotor.T_k.col(i) = -quadrotor.M*u_k.col(i);
// 		quadrotor.T_k_min.col(i) = -quadrotor.M*u_k.col(i) - Eigen::MatrixXf::Ones(4,1)*quadrotor.T_min;
// 		quadrotor.T_k_max.col(i) = Eigen::MatrixXf::Ones(4,1)*quadrotor.T_max - quadrotor.M*u_k.col(i);

// 		g_barrier(i,0) = (quadrotor.phi_max*quadrotor.phi_max - phi_k(i,0)*phi_k(i,0)) * (quadrotor.theta_max*quadrotor.theta_max - theta_k(i,0)*theta_k(i,0) ) * quadrotor.T_k_min.col(i).prod()*quadrotor.T_k_max.col(i).prod(); 

// 		if (g_barrier(i,0) <= 0)
// 		{
// 			g_barrier = Eigen::MatrixXf::Ones(mpc_params.T,1);
// 			return true;
// 		}

// 	} // for (unsigned short int i = 0; i < mpc_params.T; i++)

// 	quadrotor.eig_phi = phi_k;	
// 	quadrotor.eig_theta = theta_k;
// 	lambda_k = virtual_control_k;

// 	// cout << "phi_k " << phi_k << endl;
// 	// cout << "theta_k " << theta_k << endl;
// 	// cout << "quadrotor.T_k: " << quadrotor.T_k << endl;

// 	cout << "<FMPC> compute_roll_pitch_thrust_for_lambda_k_and_g_barrier complete" << endl;

// 	return false;

// } // bool F_MPC_UNCUT::compute_roll_pitch_thrust_for_lambda_k_and_g_barrier(Eigen::MatrixXf* chi, Eigen::MatrixXf* lambda)



// Member function of F_MPC_UNCUT object
// compute_u1_u1dot: computes the first and second integral of u1ddot using a discrete time 
// model.
// INPUTS: pointer to a variable storing the virtual control input
// OUTPUTS: none
void F_MPC_UNCUT::compute_u1_u1dot(Eigen::MatrixXf* eig_U)
{

	// Compute \delta_k, see equation 6 of Ch3. Marshall et. al. 
	quadrotor.delta_k(0,0) = quadrotor.u1_k(0,0);
	quadrotor.delta_k(1,0) = quadrotor.u1_dot_k(0,0);

	// Iterate over the number of time steps (minus 1)
	for (int i = 0; i < mpc_params.T; i++)
	{
		
		quadrotor.delta_k.block(0,i,2,1) = mpc_params.A_u1*quadrotor.delta_k.block(0,i,2,1) + mpc_params.B_u1*(*eig_U)(0,i);
		quadrotor.u1_k(i,0) = quadrotor.delta_k(0,i); 
		quadrotor.u1_dot_k(i,0) = quadrotor.delta_k(1,i); 

	} // for (int i = 0; i < mpc_params.T-1; i++)

	cout << "<FMPC> compute_u1_u1dot complete" << endl;

} // void F_MPC_UNCUT::compute_u1_u1dot(Eigen::MatrixXf* eig_U)


// Member function of F_MPC_UNCUT object
// compute_box_constraints: computes box constraints used to compute
// later segements of the planned trajectory
// INPUTS: integer storing the current segment number, and an integer indicating
// the index of the goal position in the planned path
// OUTPUTS: none
void F_MPC_UNCUT::compute_box_constraints(int segment_number, int goalStride)
{

	// Define Eigen::MatrixXf to store the current and previous goal position 
	Eigen::MatrixXf rk1 = Eigen::MatrixXf::Zero(3,1);
	Eigen::MatrixXf rk2 = Eigen::MatrixXf::Zero(3,1);

	// Float storing the angle between the x-axis
	float dir;

	// Declare Eigen::Vector3f to store points defining box constraints
	Eigen::Vector3f n; 		// normal vector
	Eigen::Vector3f per;	// vector perpendicular to n
	Eigen::Vector3f p1;		// p1 is a point on the first box constraint
	Eigen::Vector3f p2;		// p2 is a point on the second box constraint
	Eigen::Vector3f p3;		// p3 is a point on the third box constraint
	Eigen::Vector3f p4;		// p4 is a point on the fourth box constraint

	// Floats capturing offsets for each box constraint 
	float dp1,dp2,dp3,dp4;

	// Define Eigen::MatrixXf to store the parameters for the box constraints (a,b,c,d as in a plane equation)
	Eigen::MatrixXf temp_Fx = Eigen::MatrixXf::Zero(8,quadrotor.n); // Stores ai, bi, ci in the first three columns, zeros in the remaining columns
	Eigen::MatrixXf temp_fx = Eigen::MatrixXf::Zero(8,1);			// Store di in the ith row

	// If the goal stride does not exceed the length of the path
	if (goalStride < localPath.rows() && goalStride-1 >= 0)
	{
		
		// Set rk1 and rk2 (previous waypoint and current waypoint for this segment)
		rk1 = localPath.row(goalStride-1).transpose();
		rk2 = localPath.row(goalStride).transpose();

	} // if (goalStride < localPath.rows())
	else
	{
		
		// Set rk1 and rk2 (previous waypoint and current waypoint for this segment)
		rk1 = localPath.row(localPath.rows()-2).transpose();
		rk2 = localPath.row(localPath.rows()-1).transpose();

	} // if (goalStride < localPath.rows())


	// If the points are not too close together (on voxel map with resolution = 0.2m, this means the points are at the same xy position, and different z-positions)
	if ( (rk1-rk2).norm() > 0.1 )
	{
		
		dir = atan2(rk2(1,0)-rk1(1,0),rk2(0,0)-rk1(0,0));
		// dir = atan2(rk2(0,0)-rk1(0,0),rk2(1,0)-rk1(1,0));

	} // if ( (rk1-rk2).norm() > 0.1 )
	else
	{
		
		dir = 0;
	
	} // if ( (rk1-rk2).norm() > 0.1 )

	// Set the coefficients of the normal vector
	n(0) = sin(dir);
	n(1) = cos(dir);
	n(2) = 0;

	// Set the coefficients of the perpendicular vector
	per(0) = -n(1);
	per(1) = n(0);
	per(2) = 0;

	float offset = 6.0;

	// Compute the points on the planes
	p1 = rk1 + offset * ( per - n );
	p2 = rk1 + offset * (- per - n );
	p3 = rk2 + offset * ( per + n );
	p4 = rk2 + offset * (- per + n );

	// Scalar representation of the points on the planes
	dp1 = n.dot(p1);
	dp2 = per.dot(p2);
	dp3 = -per.dot(p3);
	dp4 = -n.dot(p4);

	// Flip the sign for proper comparison
	dp1*=-1;dp2*=-1;dp3*=-1;dp4*=-1;

	// Set the matrix storing the plane coefficients
	temp_Fx.block(0,0,1,3) = -n.transpose();
	temp_Fx.block(1,0,1,3) = -per.transpose();
	temp_Fx.block(2,0,1,3) = per.transpose();
	temp_Fx.block(3,0,1,3) = n.transpose();
	temp_Fx.block(4,0,2,14) = quadrotor.eig_Fpsi_hard;
	temp_Fx.block(6,0,2,14) = quadrotor.eig_Fz_hard_ceiling;

	// Set the vector storing the offset coefficients
	temp_fx(0,0) = dp1;
	temp_fx(1,0) = dp2;
	temp_fx(2,0) = dp3;
	temp_fx(3,0) = dp4;
	// ... and the other constraints
	temp_fx.block(4,0,2,1) = quadrotor.eig_psi_hard_bounds;
	temp_fx.block(6,0,2,1) = quadrotor.eig_z_hard_bounds;

	// Set the hard constraint matrix
	quadrotor.eig_Fx_hard.block(0,0,8,quadrotor.n) = temp_Fx;

	// Store the hard constraint matrix in P
	mpc_params.P[segment_number].block(0,0,8,quadrotor.n) = temp_Fx;

	// Iterate over the number of time steps
	for (int i = 0; i < mpc_params.T; i++)
	{		

		quadrotor.eig_x_hard_bounds.block(0,i,8,1) = temp_fx;
		mpc_params.h[segment_number].block(0,i,8,1) = temp_fx;

	} // for (int i = 0; i < mpc_params.T; i++)

	cout << "<FMPC> compute_box_constraints complete" << endl;

} // void F_MPC_UNCUT::compute_box_constraints(int segment_number, int goalStride)


// Member function of F_MPC_UNCUT object
// fmpc_hard_constraints: computes box constraints used to compute
// later segements of the planned trajectory
// INPUTS: pointer to eigen matrix storing the state, integer representing the segment number, 
// integer representing the index of the goal in the planned path
// OUTPUTS: none
void F_MPC_UNCUT::fmpc_hard_constraints(Eigen::MatrixXf* eig_X, int segment_number, int goalStride, int numFailures)
{

	// Delete float arrays storing constraint data
	delete[] bomt.FxX_array;
	delete[] bomt.dx_array;
	
	// Floats storing the sign of state variables
	float sign_of_heading_constraint, state_temp = 0.0;

	// Define Eigen::VectorXf storing the orientation of the UAV and of the goal
	Eigen::VectorXf quad_dir = Eigen::VectorXf::Zero(3), goal_dir = Eigen::VectorXf::Zero(3);

	// Compute the hard constraints on the vehicle's heading

	// cout << "(*eig_X)(12,0): " << (*eig_X)(12,0) << endl;

	quadrotor.eig_psi_hard_bounds(0) = (float) quadrotor.half_camera_FOV + goal(3,0);
	quadrotor.eig_psi_hard_bounds(1) = (float) quadrotor.half_camera_FOV - goal(3,0);
	
	// get_sign(sign_of_heading_constraint,quadrotor.eig_psi_hard_bounds(0));
	// int turns = floor(quadrotor.eig_psi_hard_bounds(0)/((sign_of_heading_constraint)*M_PI));
	// if (turns > 0)
	// {
	// 	quadrotor.eig_psi_hard_bounds(0) = quadrotor.eig_psi_hard_bounds(0) - sign_of_heading_constraint*M_PI*turns;
	// }

	// get_sign(sign_of_heading_constraint,quadrotor.eig_psi_hard_bounds(1));
	// turns = floor(quadrotor.eig_psi_hard_bounds(1)/((sign_of_heading_constraint)*M_PI));
	// if (turns > 0)
	// {
	// 	quadrotor.eig_psi_hard_bounds(1) = quadrotor.eig_psi_hard_bounds(1) - sign_of_heading_constraint*M_PI*turns;
	// }

	// cout << "eig_psi_hard_bounds: " << quadrotor.eig_psi_hard_bounds << endl;
	// cout << "x0: " << quadrotor.X0.block(0,0,3,1) << endl;
	// cout << "goal" << goal.block(0,0,3,1) << endl;
	// cout << "goal heading: " << goal(3,0) << endl;


	// If the current segment is the first
	if (segment_number == 0)
	{
		// Losen the constraints if the path requires us to turn around
		// If the heading angle initially violates the hard constraints
		if ( quadrotor.X0(12) > quadrotor.eig_psi_hard_bounds(0) || -quadrotor.X0(12) > quadrotor.eig_psi_hard_bounds(1))
		{
			quadrotor.eig_psi_hard_bounds(0) = 10*3.1415;
			quadrotor.eig_psi_hard_bounds(1) = 10*3.1415;
		}

		// Floats capturing the distance between the UAV and the goal point in the horizontal plane
		float xdiff, ydiff;

		xdiff = goal(0,0) - quadrotor.X0(0);
		ydiff = goal(1,0) - quadrotor.X0(1);

		// If the distance between the UAV and the goal is small 
		if ( sqrt( xdiff*xdiff + ydiff*ydiff ) < 0.1 )
		{
			// Loosen the constraints
			quadrotor.eig_psi_hard_bounds(0) = 10*3.1415;
			quadrotor.eig_psi_hard_bounds(1) = 10*3.1415;
		}

		// Four for the heading and altitude constraints, 28 for the initial boundary conditions
		quadrotor.num_hard_x = collisionConstraints.rows() + 4 + 28;
		int num_obstacle_constraints = quadrotor.num_hard_x - 4 - 28;

		// Set new float arrays to store constraint data
		bomt.FxX_array = new float[quadrotor.num_hard_x*mpc_params.T];
		bomt.dx_array = new float[quadrotor.num_hard_x*mpc_params.T];

		// Set the values of the elements stored in the containers above to zero
		memset(bomt.FxX_array,0,quadrotor.num_hard_x*mpc_params.T*sizeof(*bomt.FxX_array));
		memset(bomt.dx_array,0,quadrotor.num_hard_x*mpc_params.T*sizeof(*bomt.dx_array));

		// Allocate memory for the Eigen::Matrix versions
		bomt.FxtildeX = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>> (bomt.FxX_array, quadrotor.num_hard_x, mpc_params.T);
		bomt.dx = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>> (bomt.dx_array, quadrotor.num_hard_x, mpc_params.T);	

		// Set the size of the constraint matrix Fx and fx
		quadrotor.eig_Fx_hard = Eigen::MatrixXf::Zero(quadrotor.num_hard_x,quadrotor.n);
		quadrotor.eig_x_hard_bounds = Eigen::MatrixXf::Zero(quadrotor.num_hard_x,mpc_params.T);
			
		// Total number of constraints for MPC problem
		mpc_params.l = quadrotor.eig_x_hard_bounds.rows() + quadrotor.eig_u_hard_bounds.size(); 

		// Iterate over the number of collision avoidance constraints
		for (unsigned short int i = 0; i < num_obstacle_constraints; i++)
		{
			// Set the coefficient of the collision avoidance constraints
			quadrotor.eig_Fx_hard(i,0) = collisionConstraints(i,0);
			quadrotor.eig_Fx_hard(i,1) = collisionConstraints(i,1);
			quadrotor.eig_Fx_hard(i,2) = collisionConstraints(i,2);
		
			// Iterate over the number of constraints
			for (unsigned short int j = 0; j < mpc_params.T; j++)
			{
				
				quadrotor.eig_x_hard_bounds(i,j) = collisionConstraints(i,3);

			} // for (unsigned short int j = 0; j < mpc_params.T; j++)

		} // for (unsigned short int i = 0; i < num_obstacle_constraints; i++)			

		// Set the other constraints
		quadrotor.eig_Fx_hard.block(num_obstacle_constraints,0,2,14) = quadrotor.eig_Fpsi_hard;												// Heading constraints
		quadrotor.eig_Fx_hard.block(num_obstacle_constraints + 2,0,2,14) = quadrotor.eig_Fz_hard_ceiling;		 							// Floor and ceiling constraints
		quadrotor.eig_Fx_hard.block(num_obstacle_constraints + 4,0,2*quadrotor.n,quadrotor.n) = quadrotor.eig_Fx_hard_boundary_conditions;	// Boundary conditions
		
		// Iterate over the number of time steps
		for (unsigned short int j = 0; j < mpc_params.T; j++)
		{

			// Iterate over the number of states
			for (unsigned short int i = 0; i < quadrotor.n; i++)
			{
				// if (j == 0 && i < 3)
				// {	
				// 	quadrotor.eig_x_hard_boundary_conditions(i,j) = (*eig_X)(i,0) + mpc_params.bc_epsilon*max(numFailures,1) + abs((*eig_X)(i+3))*mpc_params.delta_t;
				// 	quadrotor.eig_x_hard_boundary_conditions(i+quadrotor.n,j) = - (*eig_X)(i,0) + mpc_params.bc_epsilon*max(numFailures,1) + abs((*eig_X)(i+3))*mpc_params.delta_t;
				// }
				// else if (j == 0 && i == 12)
				// {	
				// 	quadrotor.eig_x_hard_boundary_conditions(i,j) = (*eig_X)(i,0) + mpc_params.bc_epsilon + abs((*eig_X)(i+1))*mpc_params.delta_t;
				// 	quadrotor.eig_x_hard_boundary_conditions(i+quadrotor.n,j) = - (*eig_X)(i,0) + mpc_params.bc_epsilon + abs((*eig_X)(i+1))*mpc_params.delta_t;
				// }		
				// else
				// {
					quadrotor.eig_x_hard_boundary_conditions(i,j) = 1000; // Aka big-M constraint
					quadrotor.eig_x_hard_boundary_conditions(i+quadrotor.n,j) = 1000; // Aka big-M constraint
				// }

			} // for (unsigned short int i = 0; i < quadrotor.n; i++)
			
			// Set the heading angle constraints
			// if (j == mpc_params.T-1)
			// {
				quadrotor.eig_x_hard_bounds.block(num_obstacle_constraints,j,2,1) = quadrotor.eig_psi_hard_bounds;
			// }
			// else
			// {
				// quadrotor.eig_x_hard_bounds.block(num_obstacle_constraints,j,2,1) = 3.14*Eigen::MatrixXf::Ones(2,1);
			// }

			quadrotor.eig_x_hard_bounds.block(num_obstacle_constraints + 2,j,2,1) = quadrotor.eig_z_hard_bounds;
			quadrotor.eig_x_hard_bounds.block(num_obstacle_constraints + 4,j,2*quadrotor.n,1) = quadrotor.eig_x_hard_boundary_conditions.col(j);
		}

		// Set the size of the P matrix and its elements to zero
		mpc_params.P[segment_number] = Eigen::MatrixXf::Zero(mpc_params.l, quadrotor.n+quadrotor.m);

		// Set the size of the h matrix and its elements to zero 
		mpc_params.h[segment_number] = Eigen::MatrixXf::Zero(mpc_params.l, mpc_params.T);

		// Fill the P and h matrices
		mpc_params.P[segment_number].block(0,0,quadrotor.num_hard_x,quadrotor.n) = quadrotor.eig_Fx_hard;
		mpc_params.P[segment_number].block(quadrotor.num_hard_x,quadrotor.n,quadrotor.num_hard_u,quadrotor.m) = quadrotor.eig_Fu_hard;
		mpc_params.h[segment_number].block(0,0,quadrotor.num_hard_x,mpc_params.T) = quadrotor.eig_x_hard_bounds;
		
		// Iterate over the number of time steps
		for (unsigned short int j = 0; j < mpc_params.T; j++)
		{
			
			mpc_params.h[segment_number].block(quadrotor.num_hard_x,j,quadrotor.num_hard_u,1) = quadrotor.eig_u_hard_bounds;

		} // for (unsigned short int j = 0; j < mpc_params.T; j++)

	} // if (segment_number == 0)
	else if (segment_number < mpc_params.nu_sk-1)
	{

		// Check if goal is in constraint set
		// if (validates_collision_constraints())
		// {
		// 	cout << "goal: " << goal << " validates constraint set" << endl;
		// 	quadrotor.num_hard_x = collisionConstraints.rows() + 4 + 28;
		// }
		// else
		// State constraints
		// 8 for box constraints
		quadrotor.num_hard_x = 8 + 28;

		// Set the size of the hard constraints and bounds. Their elements are set to zero
		quadrotor.eig_Fx_hard = Eigen::MatrixXf::Zero(quadrotor.num_hard_x,quadrotor.n);
		quadrotor.eig_x_hard_bounds = Eigen::MatrixXf::Zero(quadrotor.num_hard_x,mpc_params.T);

		// New float arrays to store fmpc data
		bomt.FxX_array = new float[quadrotor.num_hard_x*mpc_params.T];
		bomt.dx_array = new float[quadrotor.num_hard_x*mpc_params.T];

		// Set the elements of the containers above to zero
		memset(bomt.FxX_array, 0, quadrotor.num_hard_x*mpc_params.T*sizeof(*bomt.FxX_array));
		memset(bomt.dx_array, 0, quadrotor.num_hard_x*mpc_params.T*sizeof(*bomt.dx_array));

		// Map the containers to their eigen counterparts
		bomt.FxtildeX = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>> (bomt.FxX_array, quadrotor.num_hard_x, mpc_params.T);
		bomt.dx = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>> (bomt.dx_array, quadrotor.num_hard_x, mpc_params.T);

		// Set the size of P and h, as well as their elements to zero
		mpc_params.P[segment_number] = Eigen::MatrixXf::Zero(quadrotor.num_hard_x+quadrotor.num_hard_u, quadrotor.n+quadrotor.m);
		mpc_params.h[segment_number] = Eigen::MatrixXf::Zero(quadrotor.num_hard_x+quadrotor.num_hard_u, mpc_params.T);

		// ****************************************************** //
		// This block loosens constraints on heading angle if the //
		// goal position is too close 						      //
		// ****************************************************** //


		// If the current heading is outside of the constraint set (before planning), loosen the heading constraints
		if ( (*eig_X)(12,0) > quadrotor.eig_psi_hard_bounds(0) || -(*eig_X)(12,0) > quadrotor.eig_psi_hard_bounds(1))
		{
			
			quadrotor.eig_psi_hard_bounds(0) = 10*3.1415;
			quadrotor.eig_psi_hard_bounds(1) = 10*3.1415;

		} // if ( (*eig_X)(12,0) >= quadrotor.eig_psi_hard_bounds(0) || -(*eig_X)(12,0) >= quadrotor.eig_psi_hard_bounds(1))

		float xdiff, ydiff;

		xdiff = goal(0,0) - (*eig_X)(0,0);
		ydiff = goal(1,0) - (*eig_X)(1,0);

		if ( sqrt( xdiff*xdiff + ydiff*ydiff ) < 0.1 )
		{
			quadrotor.eig_psi_hard_bounds(0) = 10*3.1415;
			quadrotor.eig_psi_hard_bounds(1) = 10*3.1415;
		}

		// Compute box constraints for the segments occuring later in the planned trajectory
		compute_box_constraints(segment_number, goalStride);

		// Set the orientation of boundary conditions
		quadrotor.eig_Fx_hard.block(8,0,2*quadrotor.n,quadrotor.n) = quadrotor.eig_Fx_hard_boundary_conditions;	

		// Iterate over the number of time steps
		for (unsigned short int j = 0; j < mpc_params.T; j++)
		{

			// Iterate over the number of states
			for (unsigned short int i = 0; i < quadrotor.n; i++)
			{

				// else if (j == 0 && i == 12)
				// {	
				// 	quadrotor.eig_x_hard_boundary_conditions(i,j) = (*eig_X)(i,0) + mpc_params.bc_epsilon + abs((*eig_X)(i+1))*mpc_params.delta_t;
				// 	quadrotor.eig_x_hard_boundary_conditions(i+quadrotor.n,j) = - (*eig_X)(i,0) + mpc_params.bc_epsilon + abs((*eig_X)(i+1))*mpc_params.delta_t;
				// }	
				// if (j == 0)
				// {													// Position		// Offset which expands if planning fails  // Offset introduced as a function of velocity
				// 	quadrotor.eig_x_hard_boundary_conditions(i,j) = (*eig_X)(i,0) + mpc_params.bc_epsilon*max(numFailures,1) + abs((*eig_X)(i+3))*mpc_params.delta_t;
				// 	quadrotor.eig_x_hard_boundary_conditions(i+quadrotor.n,j) = - (*eig_X)(i,0) + mpc_params.bc_epsilon*max(numFailures,1) + abs((*eig_X)(i+3))*mpc_params.delta_t;
				// }
				// else 
				// {
					quadrotor.eig_x_hard_boundary_conditions(i,j) = 1000; // Aka big-M constraint
					quadrotor.eig_x_hard_boundary_conditions(i+quadrotor.n,j) = 1000; // Aka big-M constraint
				// }
				// else if (j == 0 && i > 5 && i < 9)
				// {
				// 	quadrotor.eig_x_hard_boundary_conditions(i,j) = (*eig_X)(i,0) + (mpc_params.bc_epsilon*max(numFailures,1) + abs((*eig_X)(i+3))*mpc_params.delta_t)*mpc_params.delta_t;
				// 	quadrotor.eig_x_hard_boundary_conditions(i+quadrotor.n,j) = - (*eig_X)(i,0) + (mpc_params.bc_epsilon*max(numFailures,1) + abs((*eig_X)(i+3))*mpc_params.delta_t)*mpc_params.delta_t;
				// }


			} // for (int i = 0; i < quadrotor.n; i++)
			
			// Set the boundary conditions
			quadrotor.eig_x_hard_bounds.block(8,j,2*quadrotor.n,1) = quadrotor.eig_x_hard_boundary_conditions.col(j);

		} // for (unsigned short int j = 0; j < mpc_params.T; j++)

		// Control constraints
		mpc_params.P[segment_number].block(8,0,2*quadrotor.n,quadrotor.n) = quadrotor.eig_Fx_hard_boundary_conditions;
		mpc_params.P[segment_number].block(quadrotor.num_hard_x,quadrotor.n,quadrotor.num_hard_u,quadrotor.m) = quadrotor.eig_Fu_hard;
		
		// Iterate over the number of time steps
		for (unsigned short int j = 0; j < mpc_params.T; j++)
		{
			
			mpc_params.h[segment_number].block(8,j,2*quadrotor.n,1) = quadrotor.eig_x_hard_boundary_conditions.col(j);
			mpc_params.h[segment_number].block(quadrotor.num_hard_x,j,quadrotor.num_hard_u,1) = quadrotor.eig_u_hard_bounds;

		} // for (unsigned short int j = 0; j < mpc_params.T; j++)

	} // else if (segment_number < mpc_params.nu_sk-1)
	else
	{
		// State constraints
		// 8 for box constraints
		quadrotor.num_hard_x = 8 + 28;

		// Set the size of the hard constraints and bounds. Their elements are set to zero
		quadrotor.eig_Fx_hard = Eigen::MatrixXf::Zero(quadrotor.num_hard_x,quadrotor.n);
		quadrotor.eig_x_hard_bounds = Eigen::MatrixXf::Zero(quadrotor.num_hard_x,mpc_params.T);

		// New float arrays to store fmpc data
		bomt.FxX_array = new float[quadrotor.num_hard_x*mpc_params.T];
		bomt.dx_array = new float[quadrotor.num_hard_x*mpc_params.T];

		// Set the elements of the containers above to zero
		memset(bomt.FxX_array, 0, quadrotor.num_hard_x*mpc_params.T*sizeof(*bomt.FxX_array));
		memset(bomt.dx_array, 0, quadrotor.num_hard_x*mpc_params.T*sizeof(*bomt.dx_array));

		// Map the containers to their eigen counterparts
		bomt.FxtildeX = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>> (bomt.FxX_array, quadrotor.num_hard_x, mpc_params.T);
		bomt.dx = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>> (bomt.dx_array, quadrotor.num_hard_x, mpc_params.T);

		// Set the size of P and h, as well as their elements to zero
		mpc_params.P[segment_number] = Eigen::MatrixXf::Zero(quadrotor.num_hard_x+quadrotor.num_hard_u, quadrotor.n+quadrotor.m);
		mpc_params.h[segment_number] = Eigen::MatrixXf::Zero(quadrotor.num_hard_x+quadrotor.num_hard_u, mpc_params.T);

		// If the path requires us to turn around, loosen the heading constraints
		if ( (*eig_X)(12,0) >= quadrotor.eig_psi_hard_bounds(0) || -(*eig_X)(12,0) >= quadrotor.eig_psi_hard_bounds(1))
		{

			quadrotor.eig_psi_hard_bounds(0) = 10*3.1415;
			quadrotor.eig_psi_hard_bounds(1) = 10*3.1415;

		} // if ( (*eig_X)(12,0) >= quadrotor.eig_psi_hard_bounds(0) || -(*eig_X)(12,0) >= quadrotor.eig_psi_hard_bounds(1))


		// ****************************************************** //
		// This block loosens constraints on heading angle if the //
		// goal position is too close 						      //
		// ****************************************************** //

		float xdiff, ydiff;

		xdiff = goal(0,0) - (*eig_X)(0,0);
		ydiff = goal(1,0) - (*eig_X)(1,0);

		if ( sqrt( xdiff*xdiff + ydiff*ydiff ) < 0.1 )
		{
			quadrotor.eig_psi_hard_bounds(0) = 10*3.1415;
			quadrotor.eig_psi_hard_bounds(1) = 10*3.1415;
		}
		// Compute box constraints for the segments occuring later in the planned trajectory
		compute_box_constraints(segment_number, goalStride);

		// Set the orientation of boundary conditions
		quadrotor.eig_Fx_hard.block(8,0,2*quadrotor.n,quadrotor.n) = quadrotor.eig_Fx_hard_boundary_conditions;	

		// Iterate over the number of time steps
		for (unsigned short int j = 0; j < mpc_params.T; j++)
		{

			// Iterate over the number of states 
			for (unsigned short int i = 0; i < quadrotor.n; i++)
			{

				// else if (j == 0 && i == 12)
				// {	
				// 	quadrotor.eig_x_hard_boundary_conditions(i,j) = (*eig_X)(i,0) + mpc_params.bc_epsilon + abs((*eig_X)(i+1))*mpc_params.delta_t;
				// 	quadrotor.eig_x_hard_boundary_conditions(i+quadrotor.n,j) = - (*eig_X)(i,0) + mpc_params.bc_epsilon + abs((*eig_X)(i+1))*mpc_params.delta_t;
				// }	
				// else 
				// if (j == mpc_params.T-1 && i < 3)
				// {
					// quadrotor.eig_x_hard_boundary_conditions(i,j) = localPath(goalStride-1,i) + mpc_params.bc_epsilon; // Aka big-M constraint
					// quadrotor.eig_x_hard_boundary_conditions(i+quadrotor.n,j) = -localPath(goalStride-1,i) + mpc_params.bc_epsilon; // Aka big-M constraint
				// }
				// if (j == 0)
				// {													// Position		// Offset which expands if planning fails  // Offset introduced as a function of velocity
				// 	quadrotor.eig_x_hard_boundary_conditions(i,j) = (*eig_X)(i,0) + mpc_params.bc_epsilon*max(numFailures,1) + abs((*eig_X)(i+3))*mpc_params.delta_t;
				// 	quadrotor.eig_x_hard_boundary_conditions(i+quadrotor.n,j) = - (*eig_X)(i,0) + mpc_params.bc_epsilon*max(numFailures,1) + abs((*eig_X)(i+3))*mpc_params.delta_t;
				// }
				// else
				// {
					quadrotor.eig_x_hard_boundary_conditions(i,j) = 1000; // Aka big-M constraint
					quadrotor.eig_x_hard_boundary_conditions(i+quadrotor.n,j) = 1000; // Aka big-M constraint
				// }
				// else if (j == 0 && i < 6)
				// {
				// 	quadrotor.eig_x_hard_boundary_conditions(i,j) = (*eig_X)(i,0) + (mpc_params.bc_epsilon*max(numFailures,1) + abs((*eig_X)(i+3))*mpc_params.delta_t)*mpc_params.delta_t;
				// 	quadrotor.eig_x_hard_boundary_conditions(i+quadrotor.n,j) = - (*eig_X)(i,0) + (mpc_params.bc_epsilon*max(numFailures,1) + abs((*eig_X)(i+3))*mpc_params.delta_t)*mpc_params.delta_t;
				// }
				// else if (j == 0 && i > 5 && i < 9)
				// {
				// 	quadrotor.eig_x_hard_boundary_conditions(i,j) = (*eig_X)(i,0) + (mpc_params.bc_epsilon*max(numFailures,1) + abs((*eig_X)(i+3))*mpc_params.delta_t)*mpc_params.delta_t;
				// 	quadrotor.eig_x_hard_boundary_conditions(i+quadrotor.n,j) = - (*eig_X)(i,0) + (mpc_params.bc_epsilon*max(numFailures,1) + abs((*eig_X)(i+3))*mpc_params.delta_t)*mpc_params.delta_t;
				// }

			} // for (unsigned short int i = 0; i < quadrotor.n; i++)
			
			quadrotor.eig_x_hard_bounds.block(8,j,2*quadrotor.n,1) = quadrotor.eig_x_hard_boundary_conditions.col(j);

		} // for (unsigned short int j = 0; j < mpc_params.T; j++)

		// Control constraints
		mpc_params.P[segment_number].block(8,0,2*quadrotor.n,quadrotor.n) = quadrotor.eig_Fx_hard_boundary_conditions;
		mpc_params.P[segment_number].block(quadrotor.num_hard_x,quadrotor.n,quadrotor.num_hard_u,quadrotor.m) = quadrotor.eig_Fu_hard;
		
		// Iterate over the number of time steps
		for (unsigned short int j = 0; j < mpc_params.T; j++)
		{
			
			mpc_params.h[segment_number].block(8,j,2*quadrotor.n,1) = quadrotor.eig_x_hard_boundary_conditions.col(j);
			mpc_params.h[segment_number].block(quadrotor.num_hard_x,j,quadrotor.num_hard_u,1) = quadrotor.eig_u_hard_bounds;

		} // for (unsigned short int j = 0; j < mpc_params.T; j++)

	} // if (segment_number == 0)

	cout << "<FMPC> fmpc_hard_constraints complete" << endl;

} // void F_MPC_UNCUT::fmpc_hard_constraints(Eigen::MatrixXf* eig_X, int segment_number, int goalStride, int numFailures)


// Member function of F_MPC_UNCUT object
// fmpc_soft_constraints: computes or organizes the constraints
// INPUTS: pointer to eigen matrix storing the state, integer representing the segment number, 
// integer representing the index of the goal in the planned path
// OUTPUTS: none
void F_MPC_UNCUT::fmpc_soft_constraints(Eigen::MatrixXf* eig_X, int segment_number, int goalStride)
{

	// Delete old storage containers
	delete[] bomt.FxtildeX_array;
	delete[] bomt.d_tilde_x_array;
	delete[] bomt.d_hat_x_array;
	delete[] bomt.rho_i_x_array;

	// Set the number of soft constraints on the state (the number of soft constraints on the control is constant)
	quadrotor.num_soft_x = quadrotor.num_hard_x;

	float sign_of_heading_constraint;

	// Set new storage contaniners for a few
	bomt.FxtildeX_array = new float[quadrotor.num_soft_x*mpc_params.T];
	bomt.d_tilde_x_array = new float[quadrotor.num_soft_x*mpc_params.T];
	bomt.d_hat_x_array = new float[quadrotor.num_soft_x*mpc_params.T];
	bomt.rho_i_x_array = new float[quadrotor.num_soft_x];

	// Set the elements of the new containers to zero
	memset(bomt.FxtildeX_array, 0, quadrotor.num_soft_x*mpc_params.T*sizeof(*bomt.FxtildeX_array));
	memset(bomt.d_tilde_x_array, 0, quadrotor.num_soft_x*mpc_params.T*sizeof(*bomt.d_tilde_x_array));
	memset(bomt.d_hat_x_array, 0, quadrotor.num_soft_x*mpc_params.T*sizeof(*bomt.d_hat_x_array));
	memset(bomt.rho_i_x_array, 0, quadrotor.num_soft_x*sizeof(*bomt.rho_i_x_array));

	// Map the containers to Eigen Matrices (easier to perform linear algebra using Eigen)
	bomt.FxtildeX = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>> (bomt.FxtildeX_array, quadrotor.num_hard_x, mpc_params.T);
	bomt.d_tilde_x = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>> (bomt.d_tilde_x_array, quadrotor.num_hard_x, mpc_params.T);
	bomt.d_hat_x = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>> (bomt.d_hat_x_array, quadrotor.num_hard_x, mpc_params.T);
	bomt.rho_i_x = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>> (bomt.rho_i_x_array, quadrotor.num_hard_x, 1);

	// Resize the vector storing the value of rho for each constraint
	bomt.rho_i_x_vec.resize(quadrotor.num_soft_x);

	// Set the soft constraint matrix and bound size and elements to zero
	quadrotor.eig_Fx_soft = Eigen::MatrixXf::Zero(quadrotor.num_soft_x,quadrotor.n);
	quadrotor.eig_x_soft_bounds = Eigen::MatrixXf::Zero(quadrotor.num_soft_x,mpc_params.T);

	// Set the heading angle soft constraints
	quadrotor.eig_psi_soft_bounds(0) = (float) quadrotor.half_camera_FOV*mpc_params.collisionAvoidanceSoftConstraintOffset + goal(3,0);
	quadrotor.eig_psi_soft_bounds(1) = (float) quadrotor.half_camera_FOV*mpc_params.collisionAvoidanceSoftConstraintOffset - goal(3,0);
	
	get_sign(sign_of_heading_constraint,quadrotor.eig_psi_soft_bounds(0));
	int turns = floor(quadrotor.eig_psi_soft_bounds(0)/((sign_of_heading_constraint)*M_PI));
	if (turns > 0)
	{
		quadrotor.eig_psi_soft_bounds(0) = quadrotor.eig_psi_soft_bounds(0) - sign_of_heading_constraint*M_PI*turns;
	}

	get_sign(sign_of_heading_constraint,quadrotor.eig_psi_soft_bounds(1));
	turns = floor(quadrotor.eig_psi_soft_bounds(1)/((sign_of_heading_constraint)*M_PI));
	if (turns > 0)
	{
		quadrotor.eig_psi_soft_bounds(1) = quadrotor.eig_psi_soft_bounds(1) - sign_of_heading_constraint*M_PI*turns;
	}

	if (segment_number == 0)
	{

		// ****************************************************** //
		// This block loosens constraints on heading angle if the //
		// goal position is too close 						      //
		// ****************************************************** //

		// Compute the number of collision avoidance constraints
		int num_obstacle_constraints = quadrotor.num_hard_x - 4 - 28;

		// Loosen the constraints if the path requires us to turn around
		if ( quadrotor.X0(12) >= quadrotor.eig_psi_soft_bounds(0) || -quadrotor.X0(12) >= quadrotor.eig_psi_soft_bounds(1))
		{
			quadrotor.eig_psi_soft_bounds(0) = 10*3.1415;
			quadrotor.eig_psi_soft_bounds(1) = 10*3.1415;
		}

		// ******************************************** //
		// Set the soft collision avoidance constraints //
		// ******************************************** //
		
		// Iterate over the number of collision avoidance constraints 
		for (int i = 0; i < num_obstacle_constraints; i++)
		{
		
			quadrotor.eig_Fx_soft(i,0) = collisionConstraints(i,0);
			quadrotor.eig_Fx_soft(i,1) = collisionConstraints(i,1);
			quadrotor.eig_Fx_soft(i,2) = collisionConstraints(i,2);

			// Iterate over the number of time steps
			for (int j = 0; j < mpc_params.T; j++)
			{
				
				// Set the bounds for the constraint. They are constant over the time horizon.
				quadrotor.eig_x_soft_bounds(i,j) = collisionConstraints(i,3)*mpc_params.collisionAvoidanceSoftConstraintOffset;

			} // for (int j = 0; j < mpc_params.T; j++)
		
		} // for (int i = 0; i < num_obstacle_constraints; i++)			

		// Set the remaining constraints on heading, altitude, and boundary conditions
		quadrotor.eig_Fx_soft.block(num_obstacle_constraints,0,2,14) = quadrotor.eig_Fpsi_hard;
		quadrotor.eig_Fx_soft.block(num_obstacle_constraints + 2,0,2,14) = quadrotor.eig_Fz_soft_ceiling;		
		quadrotor.eig_Fx_soft.block(num_obstacle_constraints + 4,0,2*quadrotor.n,quadrotor.n) = quadrotor.eig_Fx_hard_boundary_conditions;		

		// ******************************* //
		// Compute the boundary conditions //
		// ******************************* //

		// Iterate over the number of time steps
		for (int j = 0; j < mpc_params.T; j++)
		{

			// Iterate over the number of states
			for (int i = 0; i < quadrotor.n; i++)
			{
				
				// On first time step, for y and psi states						x,y,z and psi states
				if (j == 0 && (i == 1 || i == 3))
				{	
					
					quadrotor.eig_x_soft_boundary_conditions(i,j) = (*eig_X)(i,0) + mpc_params.bc_epsilon*mpc_params.collisionAvoidanceSoftConstraintOffset;
					quadrotor.eig_x_soft_boundary_conditions(i+quadrotor.n,j) = - (*eig_X)(i,0) + mpc_params.bc_epsilon*mpc_params.collisionAvoidanceSoftConstraintOffset;

				} // if (j == 0 && (i < 3 || i == 12))
				else
				{
					
					quadrotor.eig_x_soft_boundary_conditions(i,j) = 1000; // Aka big-M constraint
					quadrotor.eig_x_soft_boundary_conditions(i+quadrotor.n,j) = 1000; // Aka big-M constraint

				} // if (j == 0 && (i < 3 || i == 12))

			} // for (int i = 0; i < quadrotor.n; i++)
			
			// Set the bounds for the soft constraints on heading, altitude, and boundary conditions
			quadrotor.eig_x_soft_bounds.block(num_obstacle_constraints,j,2,1) = quadrotor.eig_psi_soft_bounds;
			quadrotor.eig_x_soft_bounds.block(num_obstacle_constraints + 2,j,2,1) = quadrotor.eig_z_soft_bounds;
			quadrotor.eig_x_soft_bounds.block(num_obstacle_constraints + 4,j,2*quadrotor.n,1) = quadrotor.eig_x_soft_boundary_conditions.col(j);

		} // for (int j = 0; j < mpc_params.T; j++)

		// ********************************************************** //
		// Compute the vector h, used in computation of the objective //
		// ********************************************************** //

		// Set the size of h and its elements to zero
		mpc_params.h_tilde[segment_number] = Eigen::MatrixXf::Zero(mpc_params.l, mpc_params.T);
		mpc_params.h_tilde[segment_number].block(0,0,quadrotor.num_soft_x,mpc_params.T) = quadrotor.eig_x_soft_bounds;
			
		// Iterate over the number of time steps
		for (unsigned short int j = 0; j < mpc_params.T; j++)
		{
			
			mpc_params.h_tilde[segment_number].block(quadrotor.num_soft_x,j,quadrotor.num_soft_u,1) = quadrotor.eig_u_soft_bounds;

		} // for (int j = 0; j < mpc_params.T; j++)

	}
	else if (segment_number < mpc_params.nu_sk-1)
	{

		// ****************************************************** //
		// This block loosens constraints on heading angle if the //
		// goal position is too close 						      //
		// ****************************************************** //

		// Loosen the constraints if the path requires us to turn around
		if ( (*eig_X)(12,0) >= quadrotor.eig_psi_soft_bounds(0) || -(*eig_X)(12,0) >= quadrotor.eig_psi_soft_bounds(1))
		{
			
			quadrotor.eig_psi_soft_bounds(0) = 10*3.1415;
			quadrotor.eig_psi_soft_bounds(1) = 10*3.1415;

		} //if ( (*eig_X)(12,0) >= quadrotor.eig_psi_soft_bounds(0) || -(*eig_X)(12,0) >= quadrotor.eig_psi_soft_bounds(1))

		// Set the soft constraints
		quadrotor.eig_Fx_soft = quadrotor.eig_Fx_hard;

		// ******************************* //
		// Compute the boundary conditions //
		// ******************************* //

		// Iterate over the number of time steps		
		for (unsigned short int j = 0; j < mpc_params.T; j++)
		{

			// Iterate over the number of states
			for (unsigned short int i = 0; i < quadrotor.n; i++)
			{
			
				if (j == 0)
				{	
					quadrotor.eig_x_soft_boundary_conditions(i,j) = (*eig_X)(i,0) + mpc_params.bc_epsilon*0.95;
					quadrotor.eig_x_soft_boundary_conditions(i+quadrotor.n,j) = - (*eig_X)(i,0) + mpc_params.bc_epsilon*0.95;
				}
				else
				{
					quadrotor.eig_x_soft_boundary_conditions(i,j) = 1000; // Aka big-M constraint
					quadrotor.eig_x_soft_boundary_conditions(i+quadrotor.n,j) = 1000; // Aka big-M constraint
				}
			
			} // for (unsigned short int i = 0; i < quadrotor.n; i++)
			
			// Set the soft bounds for the heading and altitude constraints
			quadrotor.eig_x_soft_bounds.block(0,j,8,1) = quadrotor.eig_x_hard_bounds.block(0,j,8,1);

		} // for (unsigned short int j = 0; j < mpc_params.T; j++)

		// Iterate over the number of time steps
		for (unsigned short int i = 0; i < mpc_params.T; i++)
		{
		
			// Set the soft bounds for the boundary conditions
			quadrotor.eig_x_soft_bounds.block(8,i,2*quadrotor.n,1) = quadrotor.eig_x_soft_boundary_conditions.col(i);
		
		}

		// ********************************************************** //
		// Compute the vector h, used in computation of the objective //
		// ********************************************************** //

		// Set the size of h and its elements to zero
		mpc_params.h_tilde[segment_number] = Eigen::MatrixXf::Zero(8+2*quadrotor.n+quadrotor.num_soft_u, mpc_params.T);
		mpc_params.h_tilde[segment_number].block(0,0,8+2*quadrotor.n,mpc_params.T) = quadrotor.eig_x_soft_bounds;
		
		// Iterate over the number of time steps
		for (unsigned short int j = 0; j < mpc_params.T; j++)
		{
		
			mpc_params.h_tilde[segment_number].block(8+2*quadrotor.n,j,quadrotor.num_soft_u,1) = quadrotor.eig_u_soft_bounds;
		
		} // for (unsigned short int j = 0; j < mpc_params.T; j++)

	} // else if (segment_number < mpc_params.nu_sk-1)
	else
	{

		// ****************************************************** //
		// This block loosens constraints on heading angle if the //
		// goal position is too close 						      //
		// ****************************************************** //

		// Loosen the constraints if the path requires us to turn around		
		if ( (*eig_X)(12,0) >= quadrotor.eig_psi_soft_bounds(0) || -(*eig_X)(12,0) >= quadrotor.eig_psi_soft_bounds(1))
		{

			quadrotor.eig_psi_soft_bounds(0) = 10*3.1415;
			quadrotor.eig_psi_soft_bounds(1) = 10*3.1415;

		} // if ( (*eig_X)(12,0) >= quadrotor.eig_psi_soft_bounds(0) || -(*eig_X)(12,0) >= quadrotor.eig_psi_soft_bounds(1))

		// Set the soft constraints
		quadrotor.eig_Fx_soft = quadrotor.eig_Fx_hard;

		// ******************************* //
		// Compute the boundary conditions //
		// ******************************* //

		// Iterate over the number of time steps
		for (int j = 0; j < mpc_params.T; j++)
		{

			// Iterate over the number of
			for (int i = 0; i < quadrotor.n; i++)
			{
				if (j == 0)
				{	
					quadrotor.eig_x_soft_boundary_conditions(i,j) = (*eig_X)(i,0) + mpc_params.bc_epsilon*0.95;
					quadrotor.eig_x_soft_boundary_conditions(i+quadrotor.n,j) = - (*eig_X)(i,0) + mpc_params.bc_epsilon*0.95;
				}
				else if (j == mpc_params.T-1 && i < 3)
				{
					
					quadrotor.eig_x_soft_boundary_conditions(i,j) = localPath(goalStride-1,i) + mpc_params.bc_epsilon*0.95; // Aka big-M constraint
					quadrotor.eig_x_soft_boundary_conditions(i+quadrotor.n,j) = -localPath(goalStride-1,i) + mpc_params.bc_epsilon*0.95; // Aka big-M constraint

				} // if (j == mpc_params.T-1 && i < 3)
				else
				{
					
					quadrotor.eig_x_soft_boundary_conditions(i,j) = 1000; // Aka big-M constraint
					quadrotor.eig_x_soft_boundary_conditions(i+quadrotor.n,j) = 1000; // Aka big-M constraint

				} // if (j == mpc_params.T-1 && i < 3)

			} // for (int i = 0; i < quadrotor.n; i++)

			// Set the bounds for the heading and altitude soft constraints
			quadrotor.eig_x_soft_bounds.block(0,j,8,1) = quadrotor.eig_x_hard_bounds.block(0,j,8,1);

		} // for (int j = 0; j < mpc_params.T; j++)

		// Set the bounds for soft constraints on the boundary conditions
		for (int i = 0; i < mpc_params.T; i++)
		{
			
			quadrotor.eig_x_soft_bounds.block(8,i,2*quadrotor.n,1) = quadrotor.eig_x_soft_boundary_conditions.col(i);

		} // for (int i = 0; i < mpc_params.T; i++)

		// ********************************************************** //
		// Compute the vector h, used in computation of the objective //
		// ********************************************************** //

		mpc_params.h_tilde[segment_number] = Eigen::MatrixXf::Zero(8+2*quadrotor.n+quadrotor.num_soft_u, mpc_params.T);
		mpc_params.h_tilde[segment_number].block(0,0,8+2*quadrotor.n,mpc_params.T) = quadrotor.eig_x_soft_bounds;
		for (int j = 0; j < mpc_params.T; j++)
		{
			
			mpc_params.h_tilde[segment_number].block(8+2*quadrotor.n,j,quadrotor.num_soft_u,1) = quadrotor.eig_u_soft_bounds;

		} // for (int j = 0; j < mpc_params.T; j++)

	} // if (segment_number == 0)

	cout << "<FMPC> fmpc_soft_constraints complete" << endl;

} // void F_MPC_UNCUT::fmpc_soft_constraints(Eigen::MatrixXf* eig_X, int segment_number, int goalStride)


// Member function of F_MPC_UNCUT
// fmpcsolve: Solves the model predictive control problem minimize (42) subject to (21), (23), and (34)
// input: quadrotor, mpc_params, bag_of_many_things: structures containing system and problem parameters; see <structures.h> for details
// output: X, U: matrices containing the state and control inputs across the time horizon, respectively. The ith column corresponds to the state or control input at time i*delta_t
bool F_MPC_UNCUT::fmpcsolve(Eigen::MatrixXf* X, Eigen::MatrixXf* U, int segment_number)
{

	// Declare floats to store the step length and residual values
	float s;
	float res, newres;

	//initialize / reset some variables
	mpc_params.mu_12 = mpc_params.mu_12_initial;
	mpc_params.mu_13 = mpc_params.mu_13_initial;

	// Iterate over the number of states
	for (int i = 0; i < quadrotor.n; i++)
	{

		// Iterate over the number of time steps
		for (int j = 0; j < mpc_params.T; j++)
		{
			
			bomt.nu(i,j) = 0.0;

		} // for (int j = 0; j < mpc_params.T; j++)

	} // for (int i = 0; i < quadrotor.n; i++)

	// Calculate Ax and store into b for the initial step
	bomt.b.col(0) = quadrotor.eig_A * (*X).col(0);

	// Indicate to underlying code that this is the first run
	bomt.initial_run = true;

	// Setup timing variables
	auto current_time = std::chrono::high_resolution_clock::now();
	auto end_time = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> timeElapsed = end_time - current_time;

	// Iterate over the number of iterations
	for ( unsigned short int i = 0; i < mpc_params.num_iterations; ++i )
	{		

		// These functions are used to determine feasible directions \deltaZ and \delta nu 
		gfgphp(&quadrotor, &mpc_params, X, U, &bomt);
		dtildhat(&quadrotor, &mpc_params, X, U, &bomt);
		rd_tilde_rp(&quadrotor, &mpc_params, X, U, &bomt);
		res = resdresp(bomt);

		// Indicate that the first run has begun
		bomt.initial_run = false;

		// If the residual is nan or inf, return
		if (isnan(res) || isinf(res))
		{
			log_f_mpc.writeToLog("<FMPC-SOLVE> Residual not finite.");
			cout << "<FMPC-SOLVE> Residual not finite." << endl;

			// Reset a few things
			bomt.res_p = 0.0;
			bomt.res_dx = 0.0;
			bomt.res_du = 0.0;
			bomt.rp = Eigen::MatrixXf::Zero(quadrotor.n,mpc_params.T);
			bomt.rd_tilde_x = Eigen::MatrixXf::Zero(quadrotor.n,mpc_params.T);
			bomt.Ctnu_x = Eigen::MatrixXf::Zero(quadrotor.n,mpc_params.T);
			bomt.dz_x = Eigen::MatrixXf::Zero(quadrotor.n,mpc_params.T);
			bomt.rd_tilde_u = Eigen::MatrixXf::Zero(quadrotor.m,mpc_params.T);
			bomt.Ctnu_u = Eigen::MatrixXf::Zero(quadrotor.m,mpc_params.T);
			bomt.dz_u = Eigen::MatrixXf::Zero(quadrotor.m,mpc_params.T);
			// bomt.TwoQX = Eigen::MatrixXf::Zero(quadrotor.n,mpc_params.T);
			// bomt.TwoRU = Eigen::MatrixXf::Zero(quadrotor.m,mpc_params.T);
			bomt.dx = Eigen::MatrixXf::Zero(quadrotor.num_hard_x,mpc_params.T);
			bomt.du = Eigen::MatrixXf::Zero(quadrotor.num_hard_u,mpc_params.T);
			bomt.d_tilde_x = Eigen::MatrixXf::Zero(quadrotor.num_soft_x,mpc_params.T);
			bomt.d_hat_x = Eigen::MatrixXf::Zero(quadrotor.num_soft_x,mpc_params.T);
			bomt.d_tilde_u = Eigen::MatrixXf::Zero(quadrotor.num_soft_u,mpc_params.T);
			bomt.d_hat_u = Eigen::MatrixXf::Zero(quadrotor.num_soft_u,mpc_params.T);
			bomt.nu = Eigen::MatrixXf::Zero(quadrotor.n,mpc_params.T);
			bomt.PtTdt_x = Eigen::MatrixXf::Zero(quadrotor.n,mpc_params.T);
			bomt.xdz_x = Eigen::MatrixXf::Zero(quadrotor.n,mpc_params.T);
			bomt.PtTdt_u = Eigen::MatrixXf::Zero(quadrotor.m,mpc_params.T);
			bomt.udz_u = Eigen::MatrixXf::Zero(quadrotor.m,mpc_params.T);
			bomt.newx = Eigen::MatrixXf::Zero(quadrotor.n,mpc_params.T);
			bomt.newu = Eigen::MatrixXf::Zero(quadrotor.m,mpc_params.T);

			return false;

		} // if (isnan(res) || isinf(res))

		// If the residual is small enough
		if (res < mpc_params.tol_sq)
			break;

		// Compute search directions
		dnudz(&quadrotor, &mpc_params, &bomt);
		s = 1.0f;

		// ****************************************************************** //
		// feasibility search (iterate until z satisfies all the constraints) //
		// ****************************************************************** //

		// Start a timer
		clock_t startTime = clock();
		double duration = 0.0;
		bool cont = true;
		int fx_count = 0;

		//this loop establishes a trajectory proposal from the start to the goal.
		while (cont)
		{
			// If, for whatever reason, a trajectory proposal cannot be found within a conservative amount of time, abort the process
			duration = (clock() - startTime) / (double)CLOCKS_PER_SEC;
			if(duration > mpc_params.delta_t)
			{

				// If the trajectory proposal after the time above is still infeasible (violates any constraints or boundary conditions)
				// then revert to the ray trajectory or the previous feasible trajectory.
				cout << "returning from duration" << endl;
				return false;

			}

			// Switch between different line search techniques
			switch(mpc_params.line_search_style)
			{
				case 0: 
					s *= mpc_params.beta;
					break;
				case 1: 
					// s = uniform_line_search();
					break;
				case 2: 
					s = dichotomous_line_search(this,X,U,&(bomt.dz_x),&(bomt.dz_u),lsp.dichoto_a1,lsp.dichoto_b1,lsp.dichoto_l,lsp.dichoto_epsilon, segment_number,s);
					break;
			}

			// Move z along the direction of dz by magnitude specified by s and store this temporary answer to zdz (p271)
			bomt.xdz_x = (*X) + s * bomt.dz_x;
			bomt.udz_u = (*U) + s * bomt.dz_u;

			// for (int k = 0; k < mpc_params.T; k++)
			// {
			// 	while (bomt.xdz_x(12,k) > M_PI){ bomt.xdz_x(12,k) -= 2*M_PI; }
			// 	while (bomt.xdz_x(12,k) < -M_PI){ bomt.xdz_x(12,k) += 2*M_PI; }
			// }			

			// Calculate Pz based on zdz calculated above
			bomt.FxX = quadrotor.eig_Fx_hard * bomt.xdz_x;
			bomt.FuU = quadrotor.eig_Fu_hard * bomt.udz_u;

			// Check for hard constraint violations. If any are present, revise the trajectory proposal.
			cont = false;

			// Iterate over the number of hard state constraints
			for(int i = 0; i < quadrotor.num_hard_x; ++i)
			{

				// Iterate over the number of time steps
				for(int j = 0; j <  mpc_params.T ; ++j)
				{

					// if (i != 1 && i != 2)
					// {

						// If the constraint is violated
						if(bomt.FxX(i,j) > quadrotor.eig_x_hard_bounds(i,j))
						{
						
							// if (segment_number > 0)
							// {
							// 	if (i < 8)
							// 	{
									// cout << "state: " << (*X).col(j).transpose() << endl;
									// cout << "goal: " << goal.transpose() << endl;
									// cout << "segment_nu1mber: " << segment_number << endl;
									// cout << "Fx rows" << i << ": " << quadrotor.eig_Fx_hard.row(i) << endl << quadrotor.eig_x_hard_bounds(i,j) << endl;		
									// std::cout << "(i, j)" << i << ", " << j << ": " << bomt.FxX(i,j) << ", " << quadrotor.eig_x_hard_bounds(i,j) << std::endl;
									// // fx_count = 0;
							// 	}
	
							// 	fx_count++;
							// }
							cont = true;
							break;

						} // if(bomt.FxX(i,j) > quadrotor.eig_x_hard_bounds(i,j))

					// } // if (i != 1 && i != 2)

				} // for(int j = 0; j <  mpc_params.T ; ++j)

				if (cont)
				{
					
					break;

				} // if (cont)

			} // for(int i = 0; i < quadrotor.num_hard_x; ++i)

			// Iterate over the number of hard control constraints
			for (int i = 0; i < quadrotor.num_hard_u; i++)
			{

				// Iterate over the number of time steps
				for(int j = 0; j < mpc_params.T; ++j)
				{

					// If the constraint is violated
					if(bomt.FuU(i,j) > quadrotor.eig_u_hard_bounds[i])
					{
						
						// if (fx_count == 5)
						// {
						// 	std::cout << "ctrl: (i, j)" << i << ", " << j << ": " << bomt.FuU(i,j) << " " << quadrotor.eig_u_hard_bounds[i] << std::endl;
							
						// 	fx_count=0;
						// }

						// fx_count++;
						cont = true;

					} // if(bomt.FuU(i,j) > quadrotor.eig_u_hard_bounds[i])

				} // for(int j = 0; j < mpc_params.T; ++j)
				
				if (cont)
				{
					
					break;

				} // if (cont)	

			} // for (int i = 0; i < quadrotor.num_hard_u; i++)

		} 

		if (mpc_params.line_search_style == 0)
		{
			s /= mpc_params.beta; //through the changes made i the preceding loop, s is multiplied by beta one more time on exit. This corrects that.
		}
		// s = golden_section_line_search();
		// s = fibonacci_line_search();

		// Move z along the direction of dz by magnitude specified by s and store this temporary answer to zdz (p271)
		bomt.xdz_x = (*X) + s * bomt.dz_x;
		bomt.udz_u = (*U) + s * bomt.dz_u;

		// for (int k = 0; k < mpc_params.T; k++)
		// {
		// 	while (bomt.xdz_x(12,k) > M_PI){ bomt.xdz_x(12,k) -= 2*M_PI; }
		// 	while (bomt.xdz_x(12,k) < -M_PI){ bomt.xdz_x(12,k) += 2*M_PI; }
		// }				

		// Calculate Pz based on zdz calculated above
		bomt.FxX = quadrotor.eig_Fx_hard * bomt.xdz_x;
		bomt.FuU = quadrotor.eig_Fu_hard * bomt.udz_u;

	 	// Store new nu and z as determined by dnu, dz, and the s found in the feasibility search
		bomt.nu_proposed.noalias() = bomt.nu;
		bomt.x_proposed.noalias() = (*X);
		bomt.u_proposed.noalias() = (*U);

		//nu := nu + s*dnu
		bomt.nu.noalias() += s*bomt.dnu;

		//z := z + s*dz
		(*X).noalias() += s*bomt.dz_x;
		(*U).noalias() += s*bomt.dz_u;

		// :Zero(quadrotor.n,mpc_params.T)

		// for (int k = 0; k < mpc_params.T; k++)
		// {
		// 	while ((*X)(12,k) > M_PI){ (*X)(12,k).noalias() -= 2*M_PI; }
		// 	while ((*X)(12,k) < -M_PI){ (*X)(12,k).noalias() += 2*M_PI; }
		// }

		// insert backtracking line search
		// while(1) vs. for loop has negligible effect on outcome of the path but while can get stuck in an infinite loop and runs much slower
		for (int j = 0; j < 10; j++)
		{
			gfgphp(&quadrotor, &mpc_params, X, U, &bomt);
			dtildhat(&quadrotor, &mpc_params, X, U, &bomt);
			rd_tilde_rp(&quadrotor, &mpc_params, X, U, &bomt);
			newres = resdresp(bomt);

			if (isnan(newres) || isinf(newres))
			{
				cout << "residual is inf" << endl;
				cout << "BoMT.res_p: " << bomt.res_p << " BoMT.res_dx: " << bomt.res_dx << " BoMT.res_du: " << bomt.res_du << endl;
				bomt.res_p = 0.0;
				bomt.res_dx = 0.0;
				bomt.res_du = 0.0;
				bomt.rp = Eigen::MatrixXf::Zero(quadrotor.n,mpc_params.T);
				bomt.rd_tilde_x = Eigen::MatrixXf::Zero(quadrotor.n,mpc_params.T);
				bomt.Ctnu_x = Eigen::MatrixXf::Zero(quadrotor.n,mpc_params.T);
				bomt.dz_x = Eigen::MatrixXf::Zero(quadrotor.n,mpc_params.T);
				bomt.rd_tilde_u = Eigen::MatrixXf::Zero(quadrotor.m,mpc_params.T);
				bomt.Ctnu_u = Eigen::MatrixXf::Zero(quadrotor.m,mpc_params.T);
				bomt.dz_u = Eigen::MatrixXf::Zero(quadrotor.m,mpc_params.T);
				// bomt.TwoQX = Eigen::MatrixXf::Zero(quadrotor.n,mpc_params.T);
				// bomt.TwoRU = Eigen::MatrixXf::Zero(quadrotor.m,mpc_params.T);
				bomt.dx = Eigen::MatrixXf::Zero(quadrotor.num_hard_x,mpc_params.T);
				bomt.du = Eigen::MatrixXf::Zero(quadrotor.num_hard_u,mpc_params.T);
				bomt.d_tilde_x = Eigen::MatrixXf::Zero(quadrotor.num_soft_x,mpc_params.T);
				bomt.d_hat_x = Eigen::MatrixXf::Zero(quadrotor.num_soft_x,mpc_params.T);
				bomt.d_tilde_u = Eigen::MatrixXf::Zero(quadrotor.num_soft_u,mpc_params.T);
				bomt.d_hat_u = Eigen::MatrixXf::Zero(quadrotor.num_soft_u,mpc_params.T);
				bomt.nu = Eigen::MatrixXf::Zero(quadrotor.n,mpc_params.T);
				bomt.PtTdt_x = Eigen::MatrixXf::Zero(quadrotor.n,mpc_params.T);
				bomt.xdz_x = Eigen::MatrixXf::Zero(quadrotor.n,mpc_params.T);
				bomt.PtTdt_u = Eigen::MatrixXf::Zero(quadrotor.m,mpc_params.T);
				bomt.udz_u = Eigen::MatrixXf::Zero(quadrotor.m,mpc_params.T);
				bomt.newx = Eigen::MatrixXf::Zero(quadrotor.n,mpc_params.T);
				bomt.newu = Eigen::MatrixXf::Zero(quadrotor.m,mpc_params.T);
				return false;
			}

			// search for better residual. If one is found, break
			if (newres <= (1 - mpc_params.alpha * s)*(1 - mpc_params.alpha * s)*res) //First factor is squared because we're using a squared norm now
				break;

			switch(mpc_params.line_search_style)
			{
				case 0: 
					s *= mpc_params.beta;
					break;
				case 1: 
					// s = uniform_line_search();
					break;
				case 2: 
					s = dichotomous_line_search(this,X,U,&(bomt.dz_x),&(bomt.dz_u),lsp.dichoto_a1,lsp.dichoto_b1,lsp.dichoto_l,lsp.dichoto_epsilon, segment_number,s);
					break;
			}			

			//nu := nu + s*dnu
			bomt.nu.noalias() = bomt.nu_proposed + s*bomt.dnu;

			//z := z + s*dz
			(*X).noalias() = bomt.x_proposed + s*bomt.dz_x;
			(*U).noalias() = bomt.u_proposed + s*bomt.dz_u;

			// for (int k = 0; k < mpc_params.T; k++)
			// {
			// 	while ((*X)(12,k) > M_PI){ (*X)(12,k).noalias() -= 2*M_PI; }
			// 	while ((*X)(12,k) < -M_PI){ (*X)(12,k).noalias() += 2*M_PI; }
			// }

		}
		// //move to next step with new rho value
		mpc_params.mu_12 *= mpc_params.beta;
		mpc_params.mu_13 *= mpc_params.gamma;

	} // for(int i = 0; i < mpc_params.num_iterations; ++i)

	 end_time = std::chrono::high_resolution_clock::now();
	 timeElapsed = end_time - current_time;
	 log_f_mpc_runtime.writeToLog( to_string( timeElapsed.count() ), 1);
	 // std::cout << ">> fMPC solve time: " << timeElapsed.count() << " sec" << std::endl;

	 cout << "<FMPC> fmpcsolve complete" << endl;
	return true;

} // bool F_MPC_UNCUT::fmpcsolve(Eigen::MatrixXf* X, Eigen::MatrixXf* U, int segment_number)

#ifndef DISABLE_COMMUNICATION
void* start_communication_thread(F_MPC_UNCUT* f_mpc_uncut)
{

	if( communication_status != 0 ){
		fprintf(stderr,"communication thread already running\n");
		return NULL;
	}
	else{
		f_mpc_uncut->communication_thread();
		return NULL;
	}
	return NULL;

}

void* start_comm_interface_thread(void *args)
{

	F_MPC_UNCUT* f_mpc_uncut = (F_MPC_UNCUT* )args;

	start_communication_thread(f_mpc_uncut);

	return NULL;

}
#endif

void* start_trajectory_planning_thread(F_MPC_UNCUT* f_mpc_uncut)
{

	if( trajectory_status != 0 ){
		fprintf(stderr,"trajectory planning thread already running\n");
		return NULL;
	}
	else{
		f_mpc_uncut->trajectory_planner_thread();
		return NULL;
	}
	return NULL;

}

void* start_traj_interface_thread(void *args)
{

	F_MPC_UNCUT* f_mpc_uncut = (F_MPC_UNCUT* )args;

	start_trajectory_planning_thread(f_mpc_uncut);

	return NULL;

}

void F_MPC_UNCUT::start()
{

	signal(SIGPIPE, quit_handler);
	signal(SIGINT, quit_handler);
	signal(SIGSEGV, quit_handler);
	signal(SIGABRT, quit_handler);

	pthread_mutex_init(&pose_lock, NULL);
	pthread_mutex_init(&path_lock, NULL);
	pthread_mutex_init(&map_lock, NULL);
	pthread_mutex_init(&constraint_lock, NULL);
	pthread_mutex_init(&trajectory_lock, NULL);
	pthread_mutex_init(&ellipsoid_lock, NULL);

	int	result_traj = pthread_create( &traj_tid, NULL, &start_traj_interface_thread, this);	
	if ( result_traj )
	{
		cout << "<TRAJ-PLANNER-FMPC> Unable to start trajectory planner thread." << endl;
		log_f_mpc.writeToLog( "<TRAJ-PLANNER-FMPC> Unable to start trajectory planner thread.");
		exit(0);
	}
	#ifndef DISABLE_COMMUNICATION
	int	result_comm = pthread_create( &comm_tid, NULL, &start_comm_interface_thread, this);	
	if ( result_comm )
	{
		cout << "<TRAJ-PLANNER-FMPC> Unable to start communication thread." << endl;
		log_f_mpc.writeToLog( "<TRAJ-PLANNER-FMPC> Unable to start communication thread.");
		exit(0);
	}

	while(!time_to_exit)
	{

	}
	#endif
}

void F_MPC_UNCUT::update_quadrotor_pose()
{

	quadrotor.X0 *= 0;
	quadrotor.X0(0) = pose[1]; // y
	quadrotor.X0(1) = -pose[0]; // dy
	quadrotor.X0(2) = -pose[2]; // psi
	quadrotor.X0(3) = pose[4]; // dpsi
	// quadrotor.X0(4) = -pose[3]; // dy
	// quadrotor.X0(5) = -pose[5]; // dz
	// quadrotor.X0(6) = pose[7]; // ddx
	// quadrotor.X0(7) = -pose[6]; // ddy
	// quadrotor.X0(8) = -pose[8]; // ddz
	// quadrotor.X0(9) = pose[10]; // dddx
	// quadrotor.X0(10) = -pose[9]; // dddy
	// quadrotor.X0(11) = -pose[11]; // dddz
	// quadrotor.X0(12) = pose[14]; // psi 	

	// get_Omega();
	// compute_Gamma();
	// compute_dEuler();
	// quadrotor.X0(13) = dEuler(2,0); // dpsi

	bomt.eX0 = Eigen::MatrixXf::Zero(quadrotor.n,1);
	for (int i = 0; i < quadrotor.n; i++)
	{
		bomt.eX0(i,0) = pose[i];
	}

	// cout << "<FMPC> update_quadrotor_pose complete" << endl;

}


void F_MPC_UNCUT::find_closest_waypoint(float* xdiff, float* ydiff, float* zdiff, int* closest_traj_iterator)
{
	float min_norm = 9999999;

	for(int i = 0; i < localPath.rows()-1; ++i)
	{

		(*xdiff) = quadrotor.X0(0) - localPath(i,0);
		(*ydiff) = quadrotor.X0(1) - localPath(i,1);
		(*zdiff) = quadrotor.X0(2) - localPath(i,2);

		if( (*xdiff)*(*xdiff) + (*ydiff)*(*ydiff) + (*zdiff)*(*zdiff) < min_norm )
		{

			min_norm = (*xdiff)*(*xdiff) + (*ydiff)*(*ydiff) + (*zdiff)*(*zdiff);
			traj_iterator = i;

		}
	}

	(*closest_traj_iterator) = traj_iterator;

	cout << "Closest point (supposedly): " << localPath.row(traj_iterator) << endl;
	cout << "<FMPC> find_closest_waypoint complete" << endl;

}


// Member function of F_MPC_UNCUT
// find_new_waypoint: this function is used to determine what the 
// next waypoint in the path should be
// INPUTS: pointer to an Eigen::MatrixXf capturing the goal,
// pointer to floats captuing the difference between the quadrotor's
// position and the goal
// OUTPUTS: none 
void F_MPC_UNCUT::find_new_waypoint(Eigen::MatrixXf* goal, float* xdiff, float* ydiff, float* zdiff)
{

	float yawdiff;
	int local_traj_iterator;
	local_traj_iterator = traj_iterator;
	Eigen::Vector3f quad_dir, yaw_des_dir;
	bool foundPoint = false;
	float norm_diff;

	// Integers capturing the voxel corresponding to the quadrotor's position
	int x,y,z,dx,dy,dz,nnwp, sx, sy, sz, exy, exz, ezy, ax, ay, az, bx, by, bz;

	// Compute a vector storing the orientation of the quadrotor's focal axis
	quad_dir(0) = cos(quadrotor.X0(12)); 
	quad_dir(1) = sin(quadrotor.X0(12)); 
	quad_dir(2) = 0.0; 

	(*xdiff) = quadrotor.X0(0) - localPath(local_traj_iterator,0);
	(*ydiff) = quadrotor.X0(1) - localPath(local_traj_iterator,1);
	(*zdiff) = quadrotor.X0(2) - localPath(local_traj_iterator,2);
	norm_diff = sqrt((*xdiff)*(*xdiff) + (*ydiff)*(*ydiff) + (*zdiff)*(*zdiff));

	// Search for the first point of the a* output that's 'goal_tolerance' away from the current position
	// Also, find the closest point of the A* output that is in the UAV's direction of motion
	// TODO:  We cannot assume that the path will never require the UAV to turn around
	// so make sure that if the closest waypoint is not in the UAV's direction of motion, that the next waypoint 
	// IN THE QUEUE is also not in the UAV's direction of motion AND not the previous waypoint
	for(; local_traj_iterator < localPath.rows() - 1; local_traj_iterator++)
	{

		(*xdiff) = quadrotor.X0(0) - localPath(local_traj_iterator,0);
		(*ydiff) = quadrotor.X0(1) - localPath(local_traj_iterator,1);
		(*zdiff) = quadrotor.X0(2) - localPath(local_traj_iterator,2);

		// yawdiff = atan2( plannedPath(local_traj_iterator,0) - quadrotor.X0(0), plannedPath(local_traj_iterator,1) - quadrotor.X0(1) );
		yawdiff = atan2( localPath(local_traj_iterator,0) - quadrotor.X0(0), localPath(local_traj_iterator,1) - quadrotor.X0(1) );

		yaw_des_dir(0) = cos(yawdiff);
		yaw_des_dir(1) = sin(yawdiff);
		yaw_des_dir(2) = 0.0;

		norm_diff = sqrt((*xdiff)*(*xdiff) + (*ydiff)*(*ydiff) + (*zdiff)*(*zdiff));

		// if(( norm_diff > mpc_params.goal_tolerance ) && (quad_dir.dot(yaw_des_dir)) >= 0.0 )
		if ( norm_diff > mpc_params.goal_tolerance )
		{

			x = floor(quadrotor.X0(0)*5.0 + 0.5);
			y = floor(quadrotor.X0(1)*5.0 + 0.5);
			z = floor(quadrotor.X0(2)*5.0 + 0.5);

			// Record the difference in position of the camera and sample point
			dx = (int) (floor(localPath(local_traj_iterator,0)*5.0 + 0.5) - x);
			dy = (int) (floor(localPath(local_traj_iterator,1)*5.0 + 0.5) - y);
			dz = (int) (floor(localPath(local_traj_iterator,2)*5.0 + 0.5) - z);			

			// Retrieve the sign of the change in each direction, store in sx/y/z
			get_sign(sx,dx); get_sign(sy,dy); get_sign(sz,dz);
			
			// magnitude of change in each direction
			ax = abs(dx); ay = abs(dy); az = abs(dz);
			
			bx = 2*ax; by = 2*ay; bz = 2*az;
			
			exy = ay-ax; exz = az-ax; ezy = ay-az; 
			
			nnwp = ax+ay+az;

			foundPoint = true;

			// If we have hit an obstacle, mark it explored, and exit the while loop
			while(pthread_mutex_trylock(&map_lock))
			{
				usleep(1);
			}
			pthread_mutex_unlock(&map_lock);
			pthread_mutex_lock(&map_lock);

				// Iterative voxel ray tracing See Ch V.3 of Graphics Gems 4
				while( nnwp-- && x < 100 && y < 30 && z < 100 && x >= 0 && y >= 0 && z >= 0)
				{
		
					if ( abs(x - floor(localPath(local_traj_iterator,0)*5+0.5)) <= 1 && abs(y - floor(localPath(local_traj_iterator,1)*5+0.5)) <= 1 && abs(z - floor(localPath(local_traj_iterator,2)*5+0.5) ) <= 1 )
					{
						cout << "ray trace got to goal!" << endl;
						break;
					}
					
					if (voxel_map[(y)][(z)][(x)] == 3)
					{

						cout << "ray trace did not get to goal!" << endl;
						local_traj_iterator-=1;
						local_traj_iterator = max(0,local_traj_iterator);
						foundPoint = true;
						break;

					} // if (voxel_map[(x)][(y)][(z)])
					// If ray tracing goes beyond the maps boundaries, stop
					else if (x + sx < 0 || y + sy < 0 || z + sz < 0 || x + sx >= 100 || y + sy >= 30 || z + sz >= 100)
					{
						break;

					} // if (x + sx < 0 || y + sy < 0 || z + sz < 0 || x + sx >= 100 || y + sy >= 30 || z + sz >= 100)
					// Otherwise, step along the ray
					else
					{

						// Determine which face of a voxel the 3D ray exits from (this allows to determine which voxel will be pierced by the ray next)
						if (exy < 0)
						{

							if (exz < 0)
							{

								x += sx;
								exy += by;
								exz += bz;

							} // if (exz < 0)
							else
							{

								z += sz;
								exz -= bx; 
								ezy += by;

							} // if (exz < 0)

						} // if (exy < 0)
						else
						{

							if (ezy < 0)
							{
								
								z += sz;
								exz -= bx; 
								ezy += by;	

							} // if (ezy < 0)
							else
							{
								
								y += sy;
								exy -= bx;
								ezy -= bz;

							} // if (ezy < 0)

						} // if (exy < 0)

					} // if (voxel_map[(x)][(y)][(z)])


				} // while( n-- && x < 100 && y < 30 && z < 100)

			pthread_mutex_unlock(&map_lock);
			
			if (foundPoint)
			{
				break;
			}

		}	

	}

	// If we did not find a point nearby that the UAV is facing, search again only for nearby points
	// if (!foundPoint)
	// {
	// 	cout << "Did not find a point that the UAV is facing" << endl;
	// 	local_traj_iterator = 1;
	// }

	// Set the goal position
	(*goal)(0,0) = localPath(local_traj_iterator,0);
	(*goal)(1,0) = localPath(local_traj_iterator,1);
	(*goal)(2,0) = localPath(local_traj_iterator,2);

	// ************************************************************ //
	// This block determines if the goal coincides with an obstacle //
	// ************************************************************ //

	// If there are some obstacle points
	if (local_obs.rows() > 0)
	{

		// Iterate over the number of obstacle points
		for (int i = 0; i < local_obs.rows(); i++)
		{

			// If the goal and the obstacle point are close enough
			if ( (local_obs.row(i)-(*goal).block(0,0,3,1).transpose()).norm() < 0.2)
			{
				cout << "goal coincides with obstacle" << endl;
				// Step back one point in the path
				local_traj_iterator--;

				// Update the goal point 
				if (local_traj_iterator >= 0 && local_traj_iterator < localPath.rows())
				{
				
					(*goal)(0,0) = localPath(local_traj_iterator,0);
					(*goal)(1,0) = localPath(local_traj_iterator,1);
					(*goal)(2,0) = localPath(local_traj_iterator,2);		
					break;

				} // if (local_traj_iterator >= 0)
				else
				{
					
					(*goal)(0,0) = localPath(0,0);
					(*goal)(1,0) = localPath(0,1);
					(*goal)(2,0) = localPath(0,2);
					break;

				} // if (local_traj_iterator >= 0)

			} // if ( (local_obs.row(i)-(*goal).block(0,0,3,1).transpose()).norm() < 0.2)

		} // for (int i = 0; i < local_obs.rows(); i++)

	} //if (local_obs.rows() > 0)

	// **************************************** //
	// Project the goal onto the constraint set //
	// **************************************** //
	while(pthread_mutex_trylock(&constraint_lock))
	{
		usleep(1);
	}
	pthread_mutex_unlock(&constraint_lock);
	pthread_mutex_lock(&constraint_lock);
						
		Eigen::Vector3f n0, n1, n2, p_temp, proj_goal;
		Eigen::MatrixXf	p = Eigen::MatrixXf::Zero(3,collisionConstraints.rows());

		float d = 0.0, z1 = 0.0, plane_size = 2.0;

		Eigen::Vector3f v, proj_quad_on_plane_temp, goal_temp, difference_temp;
		Eigen::MatrixXf proj_quad_on_plane = Eigen::MatrixXf::Zero(3,collisionConstraints.rows());
		float dist, lam = 1.0;

		bool early_break, not_in_constraint_set = false;

		// cout << "(*goal): " << (*goal) << endl;
		goal_temp = Eigen::Map<Eigen::Vector3f>((*goal).data(),3);

		// Iterate over the number of constraints
		for (int j = 0; j < collisionConstraints.rows(); j++)
		{
			// find normal vector to the plane & offset
			n0(0) = collisionConstraints(j,0);
			n0(1) = collisionConstraints(j,1);
			n0(2) = collisionConstraints(j,2);
			d = collisionConstraints(j,3);

			// Find two vectors that form an orthonormal basis for the plane
			n1(0) = -n0(1);
			n1(1) = n0(0);
			n1(2) = z1;
			n1.normalize();

			n2 = n0.cross(n1); 

			// Find a point on the plane
			p_temp = quadrotor.X0.head(3) + n2*plane_size + n1*plane_size + n0*d;
			p.col(j) = Eigen::Map<Eigen::MatrixXf>(p_temp.data(),3,1);

		} // for (int j = 0; j < collisionConstraints.rows(); j++)


		// Iterate over the number of constraints
		for (int j = 0; j < collisionConstraints.rows(); j++)
		{
			// Project the quadrotor point onto each constraint
			v = quadrotor.X0.head(3) - Eigen::Map<Eigen::Vector3f>(p.col(j).data(),3);
			
			if (abs(collisionConstraints.block(j,0,1,3).norm()) < 10e-6)
			{
				continue;
			}

			n0(0) = collisionConstraints(j,0);
			n0(1) = collisionConstraints(j,1);
			n0(2) = collisionConstraints(j,2);
			dist = v.dot(-n0);

			proj_quad_on_plane_temp = quadrotor.X0.head(3) + dist*n0;
			proj_quad_on_plane.col(j) = Eigen::Map<Eigen::MatrixXf>(proj_quad_on_plane_temp.data(),3,1);

			// Find the vector connecting the goal to the quadrotor projected onto the constraint
			difference_temp = goal_temp-proj_quad_on_plane_temp;
			difference_temp.normalize();

			// If the vectors do not point in opposite directions, 
			// the goal is not in the constraint set
			// if (-n0.dot(difference_temp) >= 0)
			if (n0.dot(goal_temp) >= collisionConstraints(j,3))
			{
				// indicate that the goal is not in the constraint set
				not_in_constraint_set = true;

				cout << "p: " << endl << p << endl;
				cout << "goal_temp: " << goal_temp << endl;
				cout << "n0: " << n0.transpose() << endl;
				cout << "dist: " << dist << endl;
				cout << "proj_quad_on_plane_temp: " << proj_quad_on_plane_temp.transpose() << endl;
				cout << "difference_temp: " << difference_temp.transpose() << endl;
				cout << "n.dot(difference_temp): " << n0.dot(difference_temp) << endl;
				// cin.ignore();

				// stop searching
				break;

			} // if (-n.dot(difference_temp) >= 0)

		} // for (int j = 0; j < collisionConstraints.rows(); j++)

		// While the goal is not in the constraint set
		while(not_in_constraint_set)
		{

			lam*=0.9;
			// cout << "lam: " << lam << endl;

			// Find a provisional goal point as the convex combination
			// of the goal and quadrotor point
			proj_goal = lam*goal_temp+(1-lam)*quadrotor.X0.head(3);

			// Indicate that we have not stopped checking the provisional
			// goal for constraint adherence
			early_break = false;
			
			// Iterate over the number of constraints
			for (int j = 0; j < collisionConstraints.rows(); j++)
			{
				
				// Find vector pointing from the provisional goal
				// to the quadrotor projected onto the plane
				difference_temp = proj_goal-proj_quad_on_plane.col(j);
				difference_temp.normalize();

				// If the vector points in the same direction as the constraint,
				// stop checking
				if (-n0.dot(difference_temp) >= 0)
				{
					
					// Indicate that we have stopped early
					early_break = true;

					// Break out of for loop
					break;

				} // if (-n.dot(difference_temp) >= 0)

			} // for (int j = 0; j < collisionConstraints.rows(); j++)

			// If we never broke the for loop above
			if (early_break)
			{
				
				// The provisional goal satisfies the constraints
				(*goal).block(0,0,3,1) = Eigen::Map<Eigen::MatrixXf>(proj_goal.data(),3,1);				

				// Indicate that the goal is in the constraint set
				not_in_constraint_set = false;
				cout << "proj_goal: " << (*goal) << endl;

			} // if (!early_break)

			if (abs(lam) < 10e-4)
			{
				break;
			}

		} // while(not_in_constraint_set)

	pthread_mutex_unlock(&constraint_lock);

	traj_iterator = local_traj_iterator;

	cout << "<FMPC> find_new_waypoint complete" << endl;

} // void F_MPC_UNCUT::find_new_waypoint(Eigen::MatrixXf* goal, float* xdiff, float* ydiff, float* zdiff)

bool closeEnough(const float& a, const float& b, const float& epsilon = std::numeric_limits<float>::epsilon()) 
{
    return (epsilon > std::abs(a - b));
}

// void F_MPC_UNCUT::compute_Euler_angles(Eigen::Matrix3f* R, Eigen::Vector3f* Eul) 
// {

//     //check for gimbal lock
//     if (closeEnough((*R)(0,2), -1.0f)) 
//     {

//         (*Eul)(0) = 0; //gimbal lock, value of x doesn't matter
//         (*Eul)(1) = M_PI / 2;
//         (*Eul)(2) = (*Eul)(0)+atan2((*R)(1,0), (*R)(2,0));
//         return;

//     } 
//     else if (closeEnough((*R)(0,2), 1.0f)) 
//     {

//         (*Eul)(0) = 0;
//         (*Eul)(1) = -M_PI / 2;
//         (*Eul)(2) = -(*Eul)(0)+atan2(-(*R)(1,0), -(*R)(1,0));
//         return;

//     } 
//     else 
//     { 

//     	//two solutions exist
//         float x1 = -asin((*R)(0,2));
//         float x2 = M_PI - x1;

//         float y1 = atan2((*R)(1,2) / cos(x1), (*R)(2,2) / cos(x1));
//         float y2 = atan2((*R)(1,2) / cos(x2), (*R)(2,2) / cos(x2));

//         float z1 = atan2((*R)(0,1) / cos(x1), (*R)(0,0) / cos(x1));
//         float z2 = atan2((*R)(0,1) / cos(x2), (*R)(0,0) / cos(x2));

//         //choose one solution to return
//         //for example the "shortest" rotation
//         if ((std::abs(x1) + std::abs(y1) + std::abs(z1)) <= (std::abs(x2) + std::abs(y2) + std::abs(z2))) 
//         {
//         	(*Eul)(0) = x1;
//         	(*Eul)(1) = y1;
//         	(*Eul)(2) = z1;
//             return;
//         } 
//         else 
//         {
//         	(*Eul)(0) = x2;
//         	(*Eul)(1) = y2;
//         	(*Eul)(2) = z2;
//             return;
//         }

//     }

//     cout << "<FMPC> compute_Euler_angles complete" << endl;

// } // void F_MPC_UNCUT::compute_Euler_angles(const Eigen::Matrix3f& R) 


void F_MPC_UNCUT::get_trajectory_goal()
{

	float xdiff, ydiff, zdiff, lambda;

	if (localPath.rows() >= 1)
	{

		for (int i = 0; i < mpc_params.nu_X; i++)
		{
			prev_seg_eig_X[i] = Eigen::MatrixXf::Zero(quadrotor.n,mpc_params.T);
			prev_seg_eig_U[i] = Eigen::MatrixXf::Zero(quadrotor.m,mpc_params.T);
		}

		// find new goal waypoint
		traj_iterator = 0;
		int closest_traj_iterator = 0;
		find_closest_waypoint(&xdiff, &ydiff, &zdiff, &closest_traj_iterator);
		find_new_waypoint(&goal, &xdiff, &ydiff, &zdiff);

		// Eigen::MatrixXf temp = Eigen::MatrixXf::Zero(1,3);
		// temp = (localPath.row(traj_iterator)-localPath.row(closest_traj_iterator)) / (localPath.row(traj_iterator)-localPath.row(closest_traj_iterator)).norm();
		
		// Eigen::Vector3f hat_n_xk = Eigen::Map<Eigen::Vector3f>(temp.data(),3);
		// Eigen::Vector3f hat_n_yk;
		// hat_n_yk(0) = 0.0; 
		// hat_n_yk(1) = 0.05; 
		// hat_n_yk(2) = 9.81;
		// hat_n_yk = (hat_n_xk.cross(hat_n_yk)) / ((hat_n_xk.cross(hat_n_yk)).norm());
		// Eigen::Vector3f hat_n_zk = hat_n_xk.cross(hat_n_yk);

		// Eigen::Matrix3f Rfk;
		// Rfk.col(0) = hat_n_xk;
		// Rfk.col(1) = hat_n_yk;
		// Rfk.col(2) = hat_n_zk;

		// if (Rfk.determinant() <= 1+0.05 && Rfk.determinant() >= 1-0.05)
		// {	
		// 	Eigen::Vector3f euler;
		// 	compute_Euler_angles(&Rfk,&euler);
		// 	// cout << "Heading from rotm: " << euler(2) << endl;
		// }

		goal(3,0) = atan2(goal(0,0)-quadrotor.X0(0), goal(1,0)-quadrotor.X0(1));
		
		prev_plannedPath = Eigen::MatrixXf::Zero(localPath.rows(),localPath.cols());
		prev_plannedPath = localPath;


		if (sqrt(xdiff*xdiff + ydiff*ydiff) < 0.313)
		{
			goal(3,0) = quadrotor.X0(12);
		}

		remaining_segments = localPath.rows() - traj_iterator;
		cout << "remaining_segments: " << remaining_segments << endl;

	}
	else
	{
		goal(0,0) = quadrotor.X0(0);
		goal(1,0) = quadrotor.X0(1);
		goal(2,0) = quadrotor.X0(2);
		goal(3,0) = quadrotor.X0(12);
		remaining_segments = 1;
	}

	cout << "<FMPC> get_trajectory_goal complete" << endl;

}

void F_MPC_UNCUT::trajectory_planner_thread()
{

	trajectory_status = 1;

	string message;
	message = "<TRAJ-PLANNER-FMPC> Starting trajectory planner thread.";
	cout << message << endl;
	log_f_mpc.writeToLog(message);

	Eigen::MatrixXf seg_state = Eigen::MatrixXf::Zero(quadrotor.n,1);
	Eigen::MatrixXf seg_control = Eigen::MatrixXf::Zero(quadrotor.m,1);
	Eigen::MatrixXf temp_goal_pos = Eigen::MatrixXf::Zero(1,3);
	Eigen::MatrixXf eig_X = Eigen::MatrixXf::Zero(quadrotor.n,mpc_params.T);
	Eigen::MatrixXf prev_eig_X = Eigen::MatrixXf::Zero(quadrotor.n,mpc_params.T);
	Eigen::MatrixXf eig_U = Eigen::MatrixXf::Zero(quadrotor.m,mpc_params.T);
	Eigen::MatrixXf prev_eig_U = Eigen::MatrixXf::Zero(quadrotor.m,mpc_params.T);
	Eigen::VectorXf temporary_pose = Eigen::VectorXf::Zero(quadrotor.n);
	Eigen::MatrixXf prev_goal = Eigen::MatrixXf::Zero(4,1);
	// Eigen::MatrixXf segmentGoals = Eigen::MatrixXf::Zero(0,4*mpc_params.nu_X);
	segmentGoals = Eigen::MatrixXf::Zero(0,4*mpc_params.nu_X);
	
	auto current_time = std::chrono::high_resolution_clock::now();
	auto current_time1 = std::chrono::high_resolution_clock::now();
	auto pauseTime_start = std::chrono::high_resolution_clock::now();
	auto end_time = std::chrono::high_resolution_clock::now();
	auto end_time1 = std::chrono::high_resolution_clock::now();
	auto pauseTime_end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> timeElapsed = end_time - current_time;
	std::chrono::duration<double> timeElapsed1 = end_time1 - current_time1;
	std::chrono::duration<double> pauseTime = pauseTime_end - pauseTime_start;
	double waitTime, sleepTime = 2500000;

	bool g_barrier_violation;
	int localPathSize, shifted_traj_iterator;
	bool badTraj = 0;
	int jjj = 0;

	for (int i = 0; i < mpc_params.nu_X; i++)
	{
		full_trajectory[i] = Eigen::MatrixXf::Zero(quadrotor.n,mpc_params.T);
		prevfull_trajectory[i] = Eigen::MatrixXf::Zero(quadrotor.n,mpc_params.T);
		full_trajectory_interf[i] = Eigen::MatrixXf::Zero(quadrotor.n,mpc_params.T);
		full_policy[i] = Eigen::MatrixXf::Zero(quadrotor.m,mpc_params.T);
		full_attitude[i] = Eigen::MatrixXf::Zero(2,mpc_params.T);
		prev_seg_eig_X[i] = Eigen::MatrixXf::Zero(quadrotor.n,mpc_params.T);
		prev_seg_eig_U[i] = Eigen::MatrixXf::Zero(quadrotor.m,mpc_params.T);
	}
	
	while(!firstPassComplete)
	{
		usleep(1000);
	}

	// usleep(100000000);

	message = "<TRAJ-PLANNER-FMPC> Entering trajectory planning loop.";
	cout << message << endl;
	log_f_mpc.writeToLog(message);	



	while(!time_to_exit)
	{

		if (firstPassComplete)
		{
			timeElapsed1 = end_time - current_time1;
			waitTime = sleepTime-(std::chrono::duration_cast<std::chrono::microseconds>(timeElapsed1)).count();
			if (waitTime > 0)
			{
				usleep(waitTime);
			}
		}
		else
		{
			usleep(sleepTime);
		}

		current_time1 = std::chrono::high_resolution_clock::now();

		end_time = std::chrono::high_resolution_clock::now();
		timeElapsed = end_time - current_time;

		while(pthread_mutex_trylock(&pose_lock))
		{
			usleep(1);
		}
		pthread_mutex_unlock(&pose_lock);
		pthread_mutex_lock(&pose_lock);


			// Store the pose
			quadrotor.prevX0 = quadrotor.X0;
			update_quadrotor_pose();

			if (firstPassComplete)
			{

				pauseTime_start = std::chrono::high_resolution_clock::now();

				while ( 1 )
				{

					pauseTime_end = std::chrono::high_resolution_clock::now();
					// check if the pose has changed enough
					if ( (quadrotor.prevX0.head(3) - quadrotor.X0.head(3)).norm() > 0.6 )
					{
						// If so, allow to continue
						break;
					}	
					
					update_quadrotor_pose();
					pauseTime = pauseTime_end - pauseTime_start; 

					// Check if enough time has passed
					if (pauseTime.count() > mpc_params.delta_t*mpc_params.T)
					{
						// If so, allow continue
						update_quadrotor_pose();
						break;
					}

				}
				// allow to continue 
			}

		pthread_mutex_unlock(&pose_lock);

		temporary_pose = quadrotor.X0;
		prev_goal = goal;

		while(pthread_mutex_trylock(&path_lock))
		{
			usleep(1);
		}
		pthread_mutex_unlock(&path_lock);
		pthread_mutex_lock(&path_lock);

			int pathSize = plannedPath.size()/3;

			if (pathSize > 0)
			{
				localPath = Eigen::MatrixXf::Zero(pathSize,3);
				localPath = plannedPath;
			}
			else
			{
				localPath = Eigen::MatrixXf::Zero(1,3);
				localPath(0,0) = quadrotor.X0(0);
				localPath(0,1) = quadrotor.X0(1);
				localPath(0,2) = quadrotor.X0(2);
			}

		pthread_mutex_unlock(&path_lock);

		while(pthread_mutex_trylock(&map_lock))
		{
			usleep(1);
		}
		pthread_mutex_unlock(&map_lock);
		pthread_mutex_lock(&map_lock);	

			// Copy the matrix of obstacle points to a "local"
			// variable. (its not local to the function, but
			// acts as if it were)
			local_obs = Eigen::MatrixXf::Zero(obs.rows(),obs.cols());
			local_obs = obs;

		pthread_mutex_unlock(&map_lock);
		
		get_trajectory_goal();

		segmentGoals.conservativeResize(segmentGoals.rows()+1,4*mpc_params.nu_X);
		segmentGoals(numTrajectoryPlans-1,0) = goal(0,0);
		segmentGoals(numTrajectoryPlans-1,1) = goal(1,0);
		segmentGoals(numTrajectoryPlans-1,2) = goal(2,0);
		segmentGoals(numTrajectoryPlans-1,3) = goal(3,0);

		mpc_params.nu_sk = max(min(remaining_segments,mpc_params.nu_X),1);
		cout << "remaining_segments: " << remaining_segments << endl;
		// cout << "mpc_params.nu_X: " << mpc_params.nu_X << endl;

		localPathSize = localPath.rows(); 

		prev_plannedPath.conservativeResizeLike(localPath);
		prev_plannedPath = localPath;

		bomt.X_goal(0,0) = goal(0,0); 
		bomt.X_goal(1,0) = goal(1,0); 
		bomt.X_goal(2,0) = goal(2,0); 
		bomt.X_goal(12,0) = goal(3,0); 	

		delete [] mpc_params.P;
		delete [] mpc_params.h;
		delete [] mpc_params.h_tilde;

		mpc_params.P = new Eigen::MatrixXf[mpc_params.nu_X];
		mpc_params.h = new Eigen::MatrixXf[mpc_params.nu_X];
		mpc_params.h_tilde = new Eigen::MatrixXf[mpc_params.nu_X];


		cout << "New planning episode!" << endl;

		int ii = 0;
		int iii = 0;

		// if ((Eigen::Map<Eigen::MatrixXf>(quadrotor.X0.head(3).data(),1,3) ).norm() < 0.6)
		// if ( (Eigen::Map<Eigen::MatrixXf>(quadrotor.X0.head(3).data(),1,3) - localPath.row(localPath.rows()-1)).norm() < 0.6)
		// {

		while(pthread_mutex_trylock(&trajectory_lock))
		{
			usleep(1);
		}
		pthread_mutex_unlock(&trajectory_lock);
		pthread_mutex_lock(&trajectory_lock);

			while ( ii < mpc_params.nu_sk )
			{

				cout << "ii: " << ii << endl;

				if (ii == 0)
				{

					Obj(0,0) = 0.0;

					temp_goal_pos(0,0) = goal(0,0);	
					temp_goal_pos(0,1) = goal(1,0);	
					temp_goal_pos(0,2) = goal(2,0);	

					cout << "goal: " << goal.transpose() << endl;

					eig_X = Eigen::MatrixXf::Zero(quadrotor.n,mpc_params.T);
					eig_U = Eigen::MatrixXf::Zero(quadrotor.m,mpc_params.T);

					for( int i = 0; i < mpc_params.T; i++ )
					{
						eig_X.col(i) = quadrotor.X0;
						eig_X.block(3,i,9,1)*=0;
						eig_X.block(13,i,1,1)*=0;
						eig_U.col(i) = Eigen::MatrixXf::Zero(4,1);
					}

					while(pthread_mutex_trylock(&map_lock))
					{
						usleep(1);
					}
					pthread_mutex_unlock(&map_lock);
					pthread_mutex_lock(&map_lock);
					
						// Find the closest obstacle to the previous goal
						find_closest_obstacle(&temp_goal_pos,ii);
					
					pthread_mutex_unlock(&map_lock);
					if (!firstSuccess)
					{

						for (int j = 0; j < mpc_params.T; j++)
						{
							g_barrier(j,0) = (quadrotor.phi_max*quadrotor.phi_max)*(quadrotor.theta_max*quadrotor.theta_max)*quadrotor.mass*9.81*0.25;
						}

						update_weighting_matrices(ii, 1);
					}

					while(pthread_mutex_trylock(&constraint_lock))
					{
						usleep(1);
					}
					pthread_mutex_unlock(&constraint_lock);
					pthread_mutex_lock(&constraint_lock);
					
						// Update the constraints
						fmpc_hard_constraints(&eig_X,ii,traj_iterator,numFailures);
						fmpc_soft_constraints(&eig_X,ii,traj_iterator);
					
					pthread_mutex_unlock(&constraint_lock);

					// If the goal for the first segment is far away, set the goal to the closest point in the path
					if (( quadrotor.X0(0) - goal(0,0) )*( quadrotor.X0(0) - goal(0,0) ) + ( quadrotor.X0(1) - goal(1,0) )*( quadrotor.X0(1) - goal(1,0) ) + ( quadrotor.X0(2) - goal(2,0) )*( quadrotor.X0(2) - goal(2,0) ) > 2.0)
					{				
					
						float xdifftemp0,ydifftemp0,zdifftemp0;
						float min_norm_temp = 99999;
						int temp_traj_iterator = 0;
						for(int i = 0; i < localPath.rows()-1; ++i)
						{

							xdifftemp0 = goal(0,0) - localPath(i,0);
							ydifftemp0 = goal(1,0) - localPath(i,1);
							zdifftemp0 = goal(2,0) - localPath(i,2);

							if( xdifftemp0*xdifftemp0 + ydifftemp0*ydifftemp0 + zdifftemp0*zdifftemp0 < min_norm_temp )
							{

								min_norm_temp = xdifftemp0*xdifftemp0 + ydifftemp0*ydifftemp0 + zdifftemp0*zdifftemp0;
								temp_traj_iterator = i;

							}
						}
						goal.block(0,0,3,1) = localPath.row(temp_traj_iterator).transpose();

					}

					// cout << "segment: " << ii << endl;
					// cout << "position: " << quadrotor.X0(0) << ", " << quadrotor.X0(1)  << ", " << quadrotor.X0(2) << " goal: " << goal.transpose() << endl;

					// cout << "eig_X: " << endl;
					// cout << eig_X << endl;
					// cout << "eig_U: " << endl;
					// cout << eig_U << endl;
					// cout << "traj_iterator: " << traj_iterator <<endl;
					// cout << "goal: " << goal << endl;
					// cout << "robs: " << r_obs.transpose() << endl;
					// cout << "eig_Fx_hard: " << endl << quadrotor.eig_Fx_hard << endl;
					// cout << "eig_x_hard_bounds: " << endl << quadrotor.eig_x_hard_bounds << endl;
					// cout << "eig_Fu_hard: " << endl << quadrotor.eig_Fu_hard << endl;
					// cout << "eig_u_hard_bounds: " << endl << quadrotor.eig_u_hard_bounds << endl;
					// cout << "g_barrier: " << endl << g_barrier << endl;

					// perform fmpc, check if it was successful
					intermediate_success = fmpcsolve(&eig_X, &eig_U, ii);
					// for (int ji = 0; ji < quadrotor.eig_Fx_hard.rows(); ji++)
					// {
					// 	for (int ij = 0; ij < mpc_params.T; ij++)
					// 	{
					// 		if (quadrotor.eig_Fx_hard.row(ji)*eig_X.col(ij) >= quadrotor.eig_x_hard_bounds(ji))
					// 		{
					// 			cout << "Constraint " << ji << " step " << ij << ": " <<  quadrotor.eig_Fx_hard.row(ji)*eig_X.col(ij) << ">=" << quadrotor.eig_x_hard_bounds(ji) << endl;
					// 			cout << "quadrotor.eig_Fx_hard: " << quadrotor.eig_Fx_hard.block(ji,0,1,quadrotor.eig_Fx_hard.cols()) << endl;
					// 		}
					// 	}
					// }

					if (intermediate_success)
					{
						
						cout << "fMPC success, checking g_barrier" << endl;
						compute_u1_u1dot(&eig_U);

						// posterior computation of roll, pitch angles, and thrust for a particular control policy
						// and evaluate g_barrier

						// g_barrier_violation = compute_roll_pitch_thrust_for_lambda_k_and_g_barrier(&eig_X, &eig_U);
						
						if ((eig_X.block(0,mpc_params.T-1,3,1) - goal.block(0,0,3,1)).norm() > 2.0)
						{
							cout << "fMPC failed: did not reach goal" << endl;
							// trajectory_plan_success = false;
						}
						
						// if (g_barrier_violation)
						// {
						// 	cout << " g_barrier violated" << endl;
						// 	// if violation, failure, adjust weighting matrices
						// 	update_weighting_matrices(ii, 1);
						// 	trajectory_plan_success = false;
						// } 
						// else
						// {
						// 	cout << " g_barrier ok" << endl;
						// 	// if no violation, success, continue to next segment
						// 	trajectory_plan_success = true;
						// 	u_k_nuX[ii] = u_k;
						// 	v_k_nuX[ii] = v_k; 
						// 	du1_k_nuX[ii] = du1_k;
						// 	ddu1_k_nuX[ii] = ddu1_k;
						// 	lambda_k_nuX[ii] = lambda_k;
						// 	zeta_k_nuX[ii] = zeta_k;
						// 	T_k_nuX[ii] = quadrotor.T_k;
						// }
					}
					else
					{
						cout << "fMPC failed" << endl;
						trajectory_plan_success = false;
					}

				}
				else
				{
					if (ii == 1)
					{
						firstSuccess = false;
						iii = 0;
					}

					if (trajectory_plan_success)
					{
						for (int j = 0; j < mpc_params.T; j++)
						{
							eig_X.block(0,j,14,1) = prev_eig_X.block(0,mpc_params.T-1,quadrotor.n,1);
							eig_U.col(j) = prev_eig_U.block(0,mpc_params.T-1,quadrotor.m,1);
							// eig_U.col(j) = Eigen::MatrixXf::Zero(4,1);
						}
						eig_X.block(3,0,9,mpc_params.T) *= 0;
						// eig_X.block(13,0,1,mpc_params.T) *= 0;
					}	

					if (trajectory_plan_success)
					{

						if (traj_iterator+ii+iii >= localPath.rows())
						{
							iii = localPath.rows()-traj_iterator-ii-1;
						}
						else
						{

							while ( !( ( localPath.block(traj_iterator+ii+iii,0,1,3).transpose() - eig_X.block(0,mpc_params.T-1,3,1) ).norm() > 0.4 ) )
							{
								iii++;
								// if (traj_iterator+ii+iii >= plannedPath.rows())
								if (traj_iterator+ii+iii >= localPath.rows())
								{
									iii--;
									break;
								}
							}

						}

						shifted_traj_iterator = traj_iterator+ii+iii;
						// cout << "shifted_traj_iterator: " << shifted_traj_iterator << endl;
						// cout << "ii: " << ii << endl;
						// cout << "iii: " << iii << endl;

						if (shifted_traj_iterator < localPath.rows())
						{
							goal.block(0,0,3,1) = localPath.row(shifted_traj_iterator).transpose();
						 	goal(3,0) = atan2(goal(0,0)- localPath(shifted_traj_iterator-1,0), goal(1,0) - localPath(shifted_traj_iterator-1,1)); 		
							
						}
						else
						{
							// Find waypoint closest to the goal of the previous segment.
							float xdifftemp,ydifftemp,zdifftemp;
							float min_norm_temp = 99999;
							for(int i = 0; i < localPath.rows()-1; ++i)
							{

								xdifftemp = goal(0,0) - localPath(i,0);
								ydifftemp = goal(1,0) - localPath(i,1);
								zdifftemp = goal(2,0) - localPath(i,2);

								if( xdifftemp*xdifftemp + ydifftemp*ydifftemp + zdifftemp*zdifftemp < min_norm_temp )
								{

									min_norm_temp = xdifftemp*xdifftemp + ydifftemp*ydifftemp + zdifftemp*zdifftemp;
									shifted_traj_iterator = i;

								}
							}
							goal.block(0,0,3,1) = localPath.row(shifted_traj_iterator).transpose();
						 	goal(3,0) = atan2(goal(0,0)- localPath(shifted_traj_iterator-1,0), goal(1,0) - localPath(shifted_traj_iterator-1,1)); 		
		
						}

						// If the goal point is an obstacle (perhaps the planned path is old) then pick one that isn't
						if (obs.rows() > 0)
						{

							for (int i = 0; i < obs.rows(); i++)
							{
								if ( (obs.row(i)-goal.block(0,0,3,1).transpose()).norm() <= 0.2)
								{
									shifted_traj_iterator--;
									if (shifted_traj_iterator >= 0)
									{			
										goal(0,0) = localPath(shifted_traj_iterator,0);
										goal(1,0) = localPath(shifted_traj_iterator,1);
										goal(2,0) = localPath(shifted_traj_iterator,2);			
										// goal(3,0) = atan2(goal(0,0)- quadrotor.X0(0), goal(1,0) - quadrotor.X0(1)); 	
										goal(3,0) = atan2( goal(1,0) - quadrotor.X0(1), goal(0,0)- quadrotor.X0(0)); 	
										break;
									}
									else
									{
										goal(0,0) = localPath(0,0);
										goal(1,0) = localPath(0,1);
										goal(2,0) = localPath(0,2);	
										// goal(3,0) = atan2(goal(0,0)- quadrotor.X0(0), goal(1,0) - quadrotor.X0(1)); 	
										goal(3,0) = atan2( goal(1,0) - quadrotor.X0(1), goal(0,0)- quadrotor.X0(0)); 	
										break;
									}
								}
							}

						}

						segmentGoals(numTrajectoryPlans-1,4*ii) = goal(0,0);
						segmentGoals(numTrajectoryPlans-1,4*ii+1) = goal(1,0);
						segmentGoals(numTrajectoryPlans-1,4*ii+2) = goal(2,0);
						segmentGoals(numTrajectoryPlans-1,4*ii+3) = goal(3,0);
					
						temp_goal_pos = goal.block(0,0,3,1).transpose();

						bomt.X_goal(0,0) = goal(0,0); 
						bomt.X_goal(1,0) = goal(1,0); 
						bomt.X_goal(2,0) = goal(2,0); 
						bomt.X_goal(12,0) = goal(3,0); 

						if (!firstSuccess)
						{
							for (int j = 0; j < mpc_params.T; j++)
							{
								g_barrier(j,0) = (quadrotor.phi_max*quadrotor.phi_max)*(quadrotor.theta_max*quadrotor.theta_max)*quadrotor.mass*9.81*0.25;
							}

							r_obs = temp_goal_pos.transpose();
						}


						while(pthread_mutex_trylock(&map_lock))
						{
							usleep(1);
						}
						pthread_mutex_unlock(&map_lock);
						pthread_mutex_lock(&map_lock);	

							// Find the closest obstacle to the previous goal
							find_closest_obstacle(&temp_goal_pos,shifted_traj_iterator);

						pthread_mutex_unlock(&map_lock);
						
						// Must update weighting matrices to account for new r_obs
						// and discount the cost by a factor depending on the 
						// segment number
						update_and_discount_weighting_matrices(ii, 1);

						while(pthread_mutex_trylock(&constraint_lock))
						{
							usleep(1);
						}
						pthread_mutex_unlock(&constraint_lock);
						pthread_mutex_lock(&constraint_lock);

							// Update the constraints
							fmpc_hard_constraints(&eig_X,ii,shifted_traj_iterator,numFailures);
							fmpc_soft_constraints(&eig_X,ii,shifted_traj_iterator);

							// cout << "collisionConstraints: " << endl << collisionConstraints << endl;

						pthread_mutex_unlock(&constraint_lock);
					
					}
					else
					{
						update_weighting_matrices(ii,0);
					}

					// cout << "eig_X: " << endl;
					// cout << eig_X << endl;
					// cout << "eig_U: " << endl;
					// cout << eig_U << endl;
					// cout << "traj_iterator: " << traj_iterator <<endl;
					// cout << "start: " << quadrotor.X0.head(3).transpose() << endl;					
					// cout << "goal: " << goal << endl;
					// cout << "robs: " << r_obs.transpose() << endl;
					// cout << "eig_Fx_hard: " << endl << quadrotor.eig_Fx_hard << endl;
					// cout << "eig_x_hard_bounds: " << endl << quadrotor.eig_x_hard_bounds << endl;
					// cout << "eig_Fu_hard: " << endl << quadrotor.eig_Fu_hard << endl;
					// cout << "eig_u_hard_bounds: " << endl << quadrotor.eig_u_hard_bounds << endl;
					// cout << "g_barrier: " << endl << g_barrier << endl;

					// perform fmpc, check if it was successful
					intermediate_success = fmpcsolve(&eig_X, &eig_U, ii);
					cout << "fmpc success: " << intermediate_success << endl;

					// cout << "Fx: " << endl << quadrotor.eig_Fx_hard.block(4,12,2,1) << endl;
					// cout << "eig_x_pos: " << endl << eig_X.block(12,0,1,mpc_params.T) << endl;
					// cout << "FxX: " << endl << quadrotor.eig_Fx_hard.block(4,12,2,1)*eig_X.block(12,0,1,mpc_params.T) << endl;
					// cout << "fx: " << endl << quadrotor.eig_x_hard_bounds.block(4,0,2,mpc_params.T) << endl;

					if (intermediate_success)
					{

						cout << "fMPC success, checking g barrier" << endl;
						// Use eig_U to determine u1 and u1_dot
						compute_u1_u1dot(&eig_U);

						// posterior computation of roll, pitch, and thrust for a particular control policy,
						// and evaluate g_barrier
						// g_barrier_violation = compute_roll_pitch_thrust_for_lambda_k_and_g_barrier(&eig_X, &eig_U);
						
						// if (g_barrier_violation)
						// {
						// 	// cout << "g_barrier < 0, external constraints violated." << endl;
						// 	// if violation, failure, adjust weighting matrices
						// 	cout << "g barrier violated" << endl;						
						// 	update_weighting_matrices(ii, 1);
						// 	trajectory_plan_success = false;
						// } 
						// else
						// {
						// 	// if no violation, success, continue to next segment
						// 	trajectory_plan_success = true;
						// 	u_k_nuX[ii] = u_k; // u_k is the same as zeta_k
						// 	v_k_nuX[ii] = v_k; 
						// 	du1_k_nuX[ii] = du1_k;
						// 	ddu1_k_nuX[ii] = ddu1_k;
						// 	lambda_k_nuX[ii] = lambda_k;
						// 	zeta_k_nuX[ii] = zeta_k;
						// 	T_k_nuX[ii] = quadrotor.T_k;
						// 	cout << "g barrier ok" << endl;
						// }

					}
					else
					{
						trajectory_plan_success = false;
					}

				}	

				if (trajectory_plan_success)
				{

					if ((eig_X.block(0,mpc_params.T-1,3,1) - goal.block(0,0,3,1)).norm() > 1.0)
					{
						message = "Plan #" + to_string(numTrajectoryPlans) + " suspicious.";
						cout << message << endl;
						log_f_mpc.writeToLog(message);
					}

					objective_function_value(&eig_X, &eig_U, ii);
					firstSuccess = true;
					prev_eig_X = eig_X;
					prev_eig_U = eig_U;
					prev_seg_eig_X[ii] = eig_X;
					prev_seg_eig_U[ii] = eig_U;
					full_trajectory[ii] = eig_X;
					full_policy[ii] = eig_U;
					ii++;
					numFailures = 0;						
				}
				else
				{
					numFailures++;
					if (numFailures == 10)
					{	

						message = "Plan #" + to_string(numTrajectoryPlans) + " failed at " + to_string(ii) + " segment.";
						cout << message << endl;
						log_f_mpc.writeToLog(message);
						prev_eig_X = eig_X;
						prev_eig_U = eig_U;

						if (ii == 0)
						{

							for (int j = 0; j < mpc_params.T; j++)
							{
								eig_X.col(j) = quadrotor.X0;
								eig_X.block(3,0,9,j)*=0;  
								eig_X.block(13,0,1,j)*=0;  
								eig_U.col(j) = Eigen::MatrixXf::Zero(4,1);
							}

							for (int jj = ii; jj < mpc_params.nu_X; jj++)
							{
								prev_seg_eig_X[jj] = eig_X;
								prev_seg_eig_U[jj] = eig_U;
								full_trajectory[jj] = eig_X;
								full_policy[jj] = eig_U;
								segmentGoals(numTrajectoryPlans-1,4*jj) = goal(0,0);
								segmentGoals(numTrajectoryPlans-1,4*jj+1) = goal(1,0);
								segmentGoals(numTrajectoryPlans-1,4*jj+2) = goal(2,0);
								segmentGoals(numTrajectoryPlans-1,4*jj+3) = goal(3,0);

							}
							// cout << "eig_X: " << eig_X << endl;
							// cout << "eig_U: " << eig_U << endl;
							firstSuccess = true;
							ii = mpc_params.nu_X;

						}
						else
						{

							for (int j = 0; j < mpc_params.T; j++)
							{
								eig_X.col(j) = prev_seg_eig_X[ii-1].col(mpc_params.T-1);
								// eig_U.col(j) = prev_seg_eig_U[ii-1].col(mpc_params.T-1);
								// eig_X.col(j) = quadrotor.X0;  
								eig_U.col(j) = Eigen::MatrixXf::Zero(4,1);
							}

							eig_X.block(3,0,9,mpc_params.T)*= 0; 
							eig_X.block(13,0,1,mpc_params.T)*=0;

							for (int jj = ii; jj < mpc_params.nu_X; jj++)
							{
								prev_seg_eig_X[jj] = eig_X;
								prev_seg_eig_U[jj] = eig_U;
								full_trajectory[jj] = eig_X;
								full_policy[jj] = eig_U;
								segmentGoals(numTrajectoryPlans-1,4*jj) = goal(0,0);
								segmentGoals(numTrajectoryPlans-1,4*jj+1) = goal(1,0);
								segmentGoals(numTrajectoryPlans-1,4*jj+2) = goal(2,0);
								segmentGoals(numTrajectoryPlans-1,4*jj+3) = goal(3,0);
							}

						}

						// eig_X.block(2,0,1,mpc_params.T) = Eigen::MatrixXf::Ones(1,mpc_params.T)*quadrotor.X0(2);
						badTraj = 1;
						numFailures = 0;
						firstSuccess = true;
						ii = mpc_params.nu_X;
						// ii++;
					}
					else
					{
						
						message = "Plan #" + to_string(numTrajectoryPlans) + " failed " + to_string(numFailures) + " times at " + to_string(ii) + " segment.";
						cout << message << endl;
						log_f_mpc.writeToLog(message);

						if (ii == 0)
						{

							for (int j = 0; j < mpc_params.T; j++)
							{
								eig_X.col(j) = quadrotor.X0;  
								eig_U.col(j) = Eigen::MatrixXf::Zero(4,1);
							}

						}
						else
						{
							for (int j = 0; j < mpc_params.T; j++)
							{
								// eig_X.col(j) = prev_seg_eig_X[ii-1].col(j);  
								eig_X.col(j) = prev_seg_eig_X[ii-1].col(mpc_params.T-1);  
								// eig_U.col(j) = prev_seg_eig_U[ii-1].col(mpc_params.T-1);
								eig_U.col(j) = Eigen::MatrixXf::Zero(4,1);
							}

						}

						eig_X.block(3,0,9,mpc_params.T)*= 0; 
						eig_X.block(13,0,1,mpc_params.T)*=0;

					}
					
				}

			}

			firstPlanningEpisode = true;
			if (mpc_params.nu_sk < mpc_params.nu_X)
			{
				for (int i = mpc_params.nu_sk; i < mpc_params.nu_X; i++)
				{
					for (unsigned int jj = 0; jj < mpc_params.T; jj++)
					{
		
						full_trajectory[i].block(0,jj,3,1) = full_trajectory[mpc_params.nu_sk-1].block(0,mpc_params.T-1,3,1);
						full_trajectory[i].block(12,jj,1,1) = full_trajectory[mpc_params.nu_sk-1].block(12,mpc_params.T-1,1,1);
						
					}
				
					full_trajectory[i].block(3,0,9,mpc_params.T) = Eigen::MatrixXf::Zero(9,mpc_params.T);
					full_trajectory[i].block(13,0,1,mpc_params.T) = Eigen::MatrixXf::Zero(1,mpc_params.T);
					full_policy[i].block(0,0,quadrotor.m,mpc_params.T) = Eigen::MatrixXf::Zero(quadrotor.m,mpc_params.T);
					u_k_nuX[i].block(0,0,1,mpc_params.T) = 9.81*Eigen::MatrixXf::Ones(1,mpc_params.T);
					u_k_nuX[i].block(1,0,quadrotor.m-1,mpc_params.T) = Eigen::MatrixXf::Zero(quadrotor.m-1,mpc_params.T);
					v_k_nuX[i].block(0,0,quadrotor.m,mpc_params.T) = Eigen::MatrixXf::Zero(quadrotor.m,mpc_params.T);
					du1_k_nuX[i].block(0,0,1,mpc_params.T) = Eigen::MatrixXf::Zero(1,mpc_params.T);
					ddu1_k_nuX[i].block(0,0,1,mpc_params.T) = Eigen::MatrixXf::Zero(1,mpc_params.T);
					lambda_k_nuX[i] = full_policy[i];
					zeta_k_nuX[i] = u_k_nuX[i];
					T_k_nuX[i] = 9.81*Eigen::MatrixXf::Ones(4,mpc_params.T)/4;
					// cout << "full_policy[i]: " << endl << full_policy[i] << endl; 
					// cout << "u_k_nuX[i]: " << endl << u_k_nuX[i] << endl; 
					// cout << "v_k_nuX[i]: " << endl << v_k_nuX[i] << endl; 
					// cout << "du1_k_nuX[i]: " << endl << du1_k_nuX[i] << endl; 
					// cout << "ddu1_k_nuX[i]: " << endl << ddu1_k_nuX[i] << endl; 
					// cout << "lambda_k_nuX[i]: " << endl << lambda_k_nuX[i] << endl; 
					// cout << "zeta_k_nuX[i]: " << endl << zeta_k_nuX[i] << endl; 
					// cout << "T_k_nuX[i]: " << endl << T_k_nuX[i] << endl; 
				}
	
				if (trajectory_plan_success)
				{

					for (int jj = mpc_params.nu_sk; jj < mpc_params.nu_X; jj++)
					{
						segmentGoals(numTrajectoryPlans-1,4*jj) = goal(0,0);
						segmentGoals(numTrajectoryPlans-1,4*jj+1) = goal(1,0);
						segmentGoals(numTrajectoryPlans-1,4*jj+2) = goal(2,0);
						segmentGoals(numTrajectoryPlans-1,4*jj+3) = goal(3,0);						
					}

				}

			}


			// cout << "Plan #: " << numTrajectoryPlans << " Pos: " << quadrotor.X0.head(3).transpose() << endl;
			// for (int i = 0; i < mpc_params.nu_X; i++)
			// {
			// 	cout << "Traj " << i << " - Goal: " << segmentGoals.block(numTrajectoryPlans-1,4*i,1,4) << endl << full_trajectory[i].block(0,0,3,mpc_params.T) << endl;
			// }
			// for (int i = 0; i < mpc_params.nu_X; i++)
			// {
			// 	cout << "Policy " << i << ": " << endl << full_policy[i].block(0,0,4,mpc_params.T) << endl;
			// }
			for (int i = 0; i < mpc_params.nu_X; i++)
			{

				if ((full_trajectory[i].block(0,0,3,mpc_params.T).array() < 0.0).any() || (full_trajectory[i].block(0,0,3,mpc_params.T).array() > 20.0).any())
				{
					badTraj = 1;
				}
			}

			// cout << "here" << endl;

			if (!badTraj)
			{
				prev_traj_iterator = traj_iterator;

				for (int i = 0; i < mpc_params.nu_X; i++)
				{
					prevfull_trajectory[i] = full_trajectory[i];
				}

				// Log the trajectory 
				log_f_mpc_data.writeToLog("Plan " + to_string(jjj));
				log_f_mpc_policy_data.writeToLog("Plan " + to_string(jjj));
				log_f_mpc_goal_data.writeToLog(to_string(get_time_u()) + ",",0);

				for (int ll = 0; ll < mpc_params.nu_X; ll++)
				{
				
					for (int kk = 0; kk < mpc_params.T; kk++)
					{

						for (int mm = 0; mm < quadrotor.n-1; mm++)
						{

							log_f_mpc_data.writeToLog(to_string(full_trajectory[ll](mm,kk)) + ",",0);

						}

						log_f_mpc_data.writeToLog(to_string(full_trajectory[ll](quadrotor.n-1,kk)) + ",",1);

						log_f_mpc_policy_data.writeToLog(to_string(u_k_nuX[ll](0,kk)) + ",",0);
						log_f_mpc_policy_data.writeToLog(to_string(du1_k_nuX[ll](0,kk)) + ",",0);
											
						for (int mm = 0; mm < quadrotor.m; mm++)
						{						
						
							log_f_mpc_policy_data.writeToLog(to_string(full_policy[ll](mm,kk)) + ",",0);
						
						}

						for (int mm = 0; mm < quadrotor.m; mm++)
						{
							log_f_mpc_policy_data.writeToLog(to_string(v_k_nuX[ll](mm,kk)) + ",",0);
						}

						for (int mm = 0; mm < quadrotor.m; mm++)
						{
							log_f_mpc_policy_data.writeToLog(to_string(zeta_k_nuX[ll](mm,kk)) + ",",0);
						}

						for (int mm = 0; mm < quadrotor.m-1; mm++)
						{
							log_f_mpc_policy_data.writeToLog(to_string(abs(T_k_nuX[ll](mm,kk))) + ",",0);
						}
						log_f_mpc_policy_data.writeToLog(to_string(abs(T_k_nuX[ll](quadrotor.m-1,kk))),1);


					}
					if (ll < mpc_params.nu_X - 1)
					{
						log_f_mpc_goal_data.writeToLog(to_string(segmentGoals(numTrajectoryPlans-1,4*ll)) + "," + to_string(segmentGoals(numTrajectoryPlans-1,4*ll+1)) + "," + to_string(segmentGoals(numTrajectoryPlans-1,4*ll+2)) + "," + to_string(segmentGoals(numTrajectoryPlans-1,4*ll+3)) + ",", 0);
					}
					else
					{
						log_f_mpc_goal_data.writeToLog(to_string(segmentGoals(numTrajectoryPlans-1,4*ll)) + "," + to_string(segmentGoals(numTrajectoryPlans-1,4*ll+1)) + "," + to_string(segmentGoals(numTrajectoryPlans-1,4*ll+2)) + "," + to_string(segmentGoals(numTrajectoryPlans-1,4*ll+3)), 1);
					}
				
				}
				
				jjj++;

			}
			else
			{
				for (int i = 0; i < mpc_params.nu_X; i++)
				{
					full_trajectory[i] = prevfull_trajectory[i];
				}
				badTraj = 0;
			}

			numTrajectoryPlans++;
			cout << "numTrajectoryPlans: " << numTrajectoryPlans << endl;

		pthread_mutex_unlock(&trajectory_lock);

		// cin.ignore();
		// exit(0);
	// }
	// else
	// {
	// 	cout << "Got to end of path!" << endl;
	// 	usleep(sleepTime);
	// }

	}


}


void F_MPC_UNCUT::objective_function_value(Eigen::MatrixXf* X, Eigen::MatrixXf* U, int segment_number)
{

	if (mpc_params.P[segment_number].cols() != quadrotor.n+quadrotor.m || mpc_params.P[segment_number].rows() > 1000)
	{
		cout << "Issue with the size of P at segment " << segment_number << ": " << mpc_params.P[segment_number].rows() << "x" << mpc_params.P[segment_number].cols() << endl;
		return;
	}

	Eigen::MatrixXf* tildeR = new Eigen::MatrixXf[mpc_params.T];
	Eigen::MatrixXf* tildeq = new Eigen::MatrixXf[mpc_params.T];
	Eigen::MatrixXf* Zk = new Eigen::MatrixXf[mpc_params.T];
	Eigen::MatrixXf goal_n = Eigen::MatrixXf::Zero(quadrotor.n,1);
	Eigen::MatrixXf ZkTRZk;
	Eigen::MatrixXf ZfkT_Rf_Zfk;
	Eigen::MatrixXf qTZk;
	Eigen::MatrixXf qfT_Zfk;
	Eigen::MatrixXf ZkTRZk_qTZk;
	Eigen::MatrixXf ZfkT_Rf_Zfk__qfT_Zfk;

	bool constraint_violation;

	Eigen::MatrixXf Obj_vec = Eigen::MatrixXf::Zero(mpc_params.T,1);
	// Eigen::MatrixXf Obj = Eigen::MatrixXf::Zero(1,1);
	Eigen::MatrixXf PhiZ = Eigen::MatrixXf::Zero(1,1);
	Eigen::MatrixXf ThetaZ = Eigen::MatrixXf::Zero(1,1);
	Eigen::MatrixXf kappaPhiZ = Eigen::MatrixXf::Zero(1,1);
	Eigen::MatrixXf rhoInverseThetaZ = Eigen::MatrixXf::Zero(1,1);
	Eigen::MatrixXf temp;
	temp = Eigen::MatrixXf::Zero(1,1);

	goal_n(0,0) = goal(0,0);
	goal_n(1,0) = goal(1,0);
	goal_n(2,0) = goal(2,0);
	goal_n(12,0) = goal(3,0);

	for (int i = 0; i < mpc_params.T-1; i++)
	{

		tildeR[i] = Eigen::MatrixXf::Zero(quadrotor.n+quadrotor.m,quadrotor.n+quadrotor.m);
		tildeq[i] = Eigen::MatrixXf::Zero(quadrotor.n+quadrotor.m,1);
		Zk[i] =  Eigen::MatrixXf::Zero(quadrotor.n+quadrotor.m,1);
		
		tildeR[i].block(0,0,quadrotor.n,quadrotor.n) = mpc_params.eig_R_rk[i]; 
		tildeR[i].block(0,quadrotor.n,quadrotor.n,quadrotor.m) = mpc_params.eig_R_r_lambda_k[i]; 
		tildeR[i].block(quadrotor.n,0,quadrotor.m,quadrotor.n) = mpc_params.eig_R_r_lambda_k[i].transpose();
		tildeR[i].block(quadrotor.n,quadrotor.n,quadrotor.m,quadrotor.m) = mpc_params.eig_R_lambda; 

		tildeq[i].block(0,0,quadrotor.n-2,1) = mpc_params.eig_q_rk[i];
		tildeq[i](12,0) = mpc_params.q_psi;
		tildeq[i](13,0) = 0;
		tildeq[i].block(quadrotor.n,0,quadrotor.m,1) = mpc_params.eig_q_lambda_k[i];

		Zk[i].block(0,0,quadrotor.n,1) = (*X).col(i) - goal_n;
		Zk[i].block(quadrotor.n,0,quadrotor.m,1) = (*U).col(i);

		for (int ii = 0; ii < mpc_params.h[segment_number].rows(); ii++)
		{
			temp = mpc_params.P[segment_number].row(ii)*(Zk[i]);
			if (mpc_params.h[segment_number](ii,i) - temp(0,0) <= 0)
			{
				PhiZ(0,0) += 1000000;
			}
			else
			{
				PhiZ(0,0) += -1*log(mpc_params.h[segment_number](ii,i) - temp(0,0) );
			}
			ThetaZ(0,0) += log(1 + exp(mpc_params.mu_13 * ( -mpc_params.h_tilde[segment_number](ii,i) + temp(0,0) ) ) );
		}
		
		ZkTRZk = Zk[i].transpose()*tildeR[i]*Zk[i];
		qTZk = tildeq[i].transpose()*Zk[i];
		ZkTRZk_qTZk = ZkTRZk + qTZk;
		Obj_vec.block(i,0,1,1) = ZkTRZk_qTZk;

	}
	
	kappaPhiZ(0,0) = mpc_params.mu_12*PhiZ(0,0); 
	rhoInverseThetaZ(0,0) = (1/mpc_params.mu_13)*ThetaZ(0,0);

	ZfkT_Rf_Zfk = ((*X).col(mpc_params.T-1) - goal_n).transpose() * mpc_params.eig_R_rf * ((*X).col(mpc_params.T-1) - goal_n);

	// cout << "Obj_vec.sum: " << Obj_vec.sum() << " term: " << ZfkT_Rf_Zfk(0,0) << " kappaPhiZ: " <<kappaPhiZ(0,0) << " rhoInverseThetaZ: " << rhoInverseThetaZ(0,0) << endl;

	Obj(0,0) += Obj_vec.sum() + ZfkT_Rf_Zfk(0,0) + kappaPhiZ(0,0) + rhoInverseThetaZ(0,0);

	// cout << "obj: " << Obj(0,0) << endl;

	delete[] tildeR;
	delete[] tildeq;
	delete[] Zk;

} 

#ifndef DISABLE_COMMUNICATION
void F_MPC_UNCUT::communication_thread()
{
	
	communication_status = 1;

    char comm_message[30] = "Hi, this is traj.\n";
	string message;
	message = "<TRAJ-PLANNER-FMPC> Starting communication thread.";
	cout << message << endl;
	log_f_mpc.writeToLog(message);

	CLIENT client("/home/julius/git/darpa_mapping_sim_server_testing/comm_server/build/socket_parameters3.txt", &time_to_exit);

	bool constraintPass = false, posePass = false, mapPass = false, pathPass = false;

	int prev_in_constraint_size = 0;
	int prev_in_path_size = 0;
	char map[300000];
	int turns;
	float sign_of_heading;
	Eigen::MatrixXf tempcollisionConstraints = Eigen::MatrixXf::Zero(1,3);

	vector<float> obs_x, obs_y, obs_z;

	float map_freq = 20; // Hz
	float pose_freq = 30; // Hz
	float goal_freq = 10; // Hz
	float path_freq = 2; // Hz
	float constraint_freq = 5; // Hz
	float trajectory_freq = 1; // Hz

	ifstream file("../../comm_server/build/communication_parameters.txt");
		
	// String to store file lines
	string file_line;
	stringstream ss;
	
	//read communication parameters
	do{ss.clear(); getline(file, file_line); ss.str(file_line);}
	while(file_line.at(0) == '/' && file_line.at(1) == '/');
	ss >> map_freq;	
	cout << "map_freq: " << map_freq << endl;

	//read communication parameters
	do{ss.clear(); getline(file, file_line); ss.str(file_line);}
	while(file_line.at(0) == '/' && file_line.at(1) == '/');
	ss >> pose_freq;	
	cout << "pose_freq: " << pose_freq << endl;

	do{ss.clear(); getline(file, file_line); ss.str(file_line);}
	while(file_line.at(0) == '/' && file_line.at(1) == '/');
	ss >> goal_freq;	
	cout << "goal_freq: " << goal_freq << endl;

	do{ss.clear(); getline(file, file_line); ss.str(file_line);}
	while(file_line.at(0) == '/' && file_line.at(1) == '/');
	ss >> path_freq;	
	cout << "path_freq: " << path_freq << endl;

	do{ss.clear(); getline(file, file_line); ss.str(file_line);}
	while(file_line.at(0) == '/' && file_line.at(1) == '/');
	ss >> constraint_freq;	
	cout << "constraint_freq: " << constraint_freq << endl;

	do{ss.clear(); getline(file, file_line); ss.str(file_line);}
	while(file_line.at(0) == '/' && file_line.at(1) == '/');
	ss >> trajectory_freq;	
	cout << "trajectory_freq: " << trajectory_freq << endl;	

	float map_time_us = 1000000/map_freq; // Microseconds	
	float pose_time_us = 1000000/(pose_freq+6); // Microseconds	
	float goal_time_us = 1000000/goal_freq; // Microseconds	
	float path_time_us = 1000000/(2+path_freq); // Microseconds	
	float constraint_time_us = 1000000/(4+constraint_freq); // Microseconds	
	float trajectory_time_us = 1000000/trajectory_freq; // Microseconds	

	auto current_time_map = std::chrono::high_resolution_clock::now();
	auto current_time_pose = std::chrono::high_resolution_clock::now();
	auto current_time_goal = std::chrono::high_resolution_clock::now();
	auto current_time_path = std::chrono::high_resolution_clock::now();
	auto current_time_constraints = std::chrono::high_resolution_clock::now();
	auto current_time_trajectory = std::chrono::high_resolution_clock::now();
	auto end_time = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> timeElapsed_map = end_time - current_time_map;	
	std::chrono::duration<double> timeElapsed_pose = end_time - current_time_pose;	
	std::chrono::duration<double> timeElapsed_goal = end_time - current_time_goal;	
	std::chrono::duration<double> timeElapsed_path = end_time - current_time_path;	
	std::chrono::duration<double> timeElapsed_constraints = end_time - current_time_constraints;	
	std::chrono::duration<double> timeElapsed_trajectory = end_time - current_time_trajectory;		

	// Pose communication vars
	std::vector<string> pose_recv;
	string pose_recv_string;
	char* pose_recv_buffer;
	pose_recv_buffer = new char[1000];

	// Path communication vars
	std::vector<string> path_recv;
	string path_recv_string;
	char* path_recv_buffer = new char[5000];	

	std::vector<string> constraints_recv;
	string constraints_recv_string;
	char* constraints_recv_buffer = new char[5000];	

	string trajectory_send;
	char* trajectory_send_buffer = new char[6000];
	char* cstr_trajectory = new char[6000];		

	int bytes_read_map = 0, bytes_read_pose = 0, bytes_read_goal = 0, bytes_read_path = 0, bytes_read_constraints = 0, bytes_read_trajectory = 0;
	int counter = 0, _interface_loops = 0;

	while(!time_to_exit)
	{

		// Get the "end" time
		end_time = std::chrono::high_resolution_clock::now();		

		/////////////////////////
		// POSE COMMUNICATIONS //
		/////////////////////////

		// if (std::chrono::duration_cast<std::chrono::microseconds>(timeElapsed_pose).count() >= 1.1*pose_time_us)
		// {
			// ioctl(client.sockets[2]->sockfd,TCFLSH,2);
			// current_time_pose = std::chrono::high_resolution_clock::now();
		// }
		// else if (std::chrono::duration_cast<std::chrono::microseconds>(timeElapsed_pose).count() >= pose_time_us)
		if (std::chrono::duration_cast<std::chrono::microseconds>(timeElapsed_pose).count() >= pose_time_us)
		{	

			// cout << "timeElapsed_pose: " << timeElapsed_pose.count() << endl;

			// RECEIVE THE POSE //
			if (client.socket_active[2])
			{

				// cout << "recv pose" << endl; 
				bytes_read_pose = client.sockets[2]->process_receiving(pose_recv_buffer,1000*sizeof(char),1);
				pose_recv_string.clear();
				// pose_recv_string.resize(1001);
				pose_recv_string = pose_recv_buffer;
				// cout << "pose_recv_string: " << pose_recv_string << endl;
				counter = 0;
				pose_recv.clear();

				if (pose_recv_string[0] == 'Q')
				{

					for (int i = 1; i < 1000; i++)
					{
						if (pose_recv_string[i] == '!')
						{
							break;
						}
						else if (pose_recv_string[i] == ',')
						{
							counter++;
							pose_recv.resize(counter);
						}
						else
						{
							pose_recv[counter-1] += pose_recv_string[i];
						}
					}
				
					while(pthread_mutex_trylock(&pose_lock))
					{
						usleep(1);
					}
					pthread_mutex_unlock(&pose_lock);
					pthread_mutex_lock(&pose_lock);

						if (counter > 0)
						{						

							for (int i = 0; i < 18; i++)
							{
								pose[i] = boost::lexical_cast<float>(pose_recv[i]);
							}

							get_sign(sign_of_heading,pose[14]);
							turns = floor(pose[14]/((sign_of_heading)*M_PI));
							if (turns > 0)
							{
								pose[14] = pose[14] - sign_of_heading*M_PI*turns;
							}

							if (!posePass)
							{
								for (int i = 0; i < 18; i++)
								{
									initialPose[i] = pose[i];
								}
							}
		
							posePass = true;
						}

					pthread_mutex_unlock(&pose_lock);
				}				

			}
			
			current_time_pose = std::chrono::high_resolution_clock::now();

		}
		else
		{
			ioctl(client.sockets[2]->sockfd,TCFLSH,2);
			// tcflush(client.sockets[2]->sockfd,TCIFLUSH);
		}

		////////////////////////
		// MAP COMMUNICATIONS //
		////////////////////////

		// If enough time has passed
		if (std::chrono::duration_cast<std::chrono::microseconds>(timeElapsed_map).count() >= map_time_us)
		{		
			// cout << "receiving map" << endl;
			// cout << "timeElapsed_map: " << timeElapsed_map.count() << endl;

			// RECEIVE THE MAP //
			if (client.socket_active[0])
			{

				counter = 0;

				bytes_read_map = client.sockets[0]->process_receiving(map,sizeof(map),0);
				// cout << "bytes_read_map: " << bytes_read_map << endl;
				// cout << "map[0]: " << +map[0] << endl;
				// if (bytes_read_map == sizeof(map))
				// {
					while(pthread_mutex_trylock(&map_lock))
					{
						usleep(1);
					}
					pthread_mutex_unlock(&map_lock);
					pthread_mutex_lock(&map_lock);
						
						obs_x.clear();
						obs_y.clear();
						obs_z.clear();
						for (int k = 0; k < grid_z; k++) 
						{
							for (int j = 0; j < grid_y; j++) 
							{
								for (int i = 0; i < grid_x; i++) 
								{
									voxel_map[i][j][k] = +map[counter];
									if (voxel_map[i][j][k] == 3)
									{
										obs_x.push_back(0.2*i);
										obs_y.push_back(0.2*k);
										obs_z.push_back(0.2*j);
									}
									counter++;

								}
							}
						}

						if (obs_x.size() > 0)
						{				
							obs = Eigen::MatrixXf::Zero(obs_x.size(),3);
							for (int i = 0; i < obs_x.size(); i++)
							{
								obs(i,0) = obs_x[i];
								obs(i,1) = obs_y[i];
								obs(i,2) = obs_z[i];
							}
						}
						
						mapPass = true;
						counter = 0;

					pthread_mutex_unlock(&map_lock);	
				// }

			} // if (client.socket_active[0])

			current_time_map = std::chrono::high_resolution_clock::now();

		} // if (std::chrono::duration_cast<std::chrono::microseconds>(timeElapsed_map).count() >= map_time_us)

		/////////////////////////
		// PATH COMMUNICATIONS //
		/////////////////////////

		if (std::chrono::duration_cast<std::chrono::microseconds>(timeElapsed_path).count() >= path_time_us)
		{

			// RECEIEVE THE PATH // 
			if (client.socket_active[6])
			{
				counter = 0;

				// cout << "Receiving the path" << endl;
				client.sockets[6]->process_receiving(path_recv_buffer,5000*sizeof(char),1);
				
				path_recv_string.clear();
				path_recv_string.resize(5000);
				path_recv_string = path_recv_buffer;

				path_recv.clear();

				if (path_recv_string[0] == 'P' && path_recv_string[2] != '!')
				{

					for (int i = 1; i < 5000; i++)
					{
						if (path_recv_string[i] == '!')
						{
							break;
						}
						else if (path_recv_string[i] == ',')
						{
							counter++;
							path_recv.resize(counter);
						}
						else
						{
							path_recv[counter-1] += path_recv_string[i];
						}
					}

				}

				// cout << "path: "; 
				// for (int i = 0; i < counter; i++)
				// {
				// 	cout << boost::lexical_cast<float>(path_recv[i]) << ",";
				// }
				// cout << endl;	

				while(pthread_mutex_trylock(&path_lock))
				{
					usleep(1);
				}
				pthread_mutex_unlock(&path_lock);
				pthread_mutex_lock(&path_lock);

					if (counter > 0)
					{

					plannedPath = Eigen::MatrixXf::Zero(path_recv.size()/3,3);
					for (int i = 0; i < path_recv.size()/3; i++)
					{
						plannedPath(i,0) = boost::lexical_cast<float>(path_recv[3*i]);
						plannedPath(i,1) = boost::lexical_cast<float>(path_recv[3*i+1]);
						plannedPath(i,2) = boost::lexical_cast<float>(path_recv[3*i+2]);
					}

						prev_in_path_size = path_recv.size()/3;

						pathPass = true;
						if (plannedPath.rows() > 0)
						{
							localPath = Eigen::MatrixXf::Zero(plannedPath.rows(),3);
							localPath = plannedPath;
						}
						else
						{
							if (!firstPassComplete)
							{
								localPath = Eigen::MatrixXf::Zero(1,3);
								localPath(0,0) = boost::lexical_cast<float>(pose_recv[0]); 
								localPath(0,1) = boost::lexical_cast<float>(pose_recv[1]); 
								localPath(0,2) = boost::lexical_cast<float>(pose_recv[2]); 
								pathPass = false;
							}
						}

					}

				pthread_mutex_unlock(&path_lock);

				counter = 0;

			}
			
			current_time_path = std::chrono::high_resolution_clock::now();

		}

		////////////////////////////////
		// CONSTRAINTS COMMUNICATIONS //
		////////////////////////////////	
		
		if (std::chrono::duration_cast<std::chrono::microseconds>(timeElapsed_constraints).count() >= constraint_time_us)
		{
			
			// RECEIEVE THE CONSTRAINTS // 
			if (client.socket_active[8])
			{
				counter = 0;
				constraints_recv.clear();

				// cout << "Receiving the constraints" << endl;
				client.sockets[8]->process_receiving(constraints_recv_buffer,5000*sizeof(char),1);
				constraints_recv_string = constraints_recv_buffer;
				// cout << "constraints_recv_string: " << endl << constraints_recv_string<<endl;

				if (constraints_recv_string[0] == 'C' && constraints_recv_string[2] != '!')
				{

					for (int i = 1; i < 5000; i++)
					{
						if (constraints_recv_string[i] == '!')
						{
							break;
						}
						else if (constraints_recv_string[i] == ',')
						{
							counter++;
							constraints_recv.resize(counter);
						}
						else
						{
							constraints_recv[counter-1] += constraints_recv_string[i];
						}
					}
	
					// cout << "constraints: "; 
					// for (int i = 0; i < counter; i++)
					// {
					// 	cout << boost::lexical_cast<float>(constraints_recv[i]) << ",";
					// }
					// cout << endl;	
			
					while(pthread_mutex_trylock(&constraint_lock))
					{
						usleep(1);
					} // while(pthread_mutex_trylock(&constraint_lock))
					pthread_mutex_unlock(&constraint_lock);
					pthread_mutex_lock(&constraint_lock);

						if (counter > 0)
						{

							Eigen::MatrixXf tempcollisionConstraints = Eigen::MatrixXf::Zero(constraints_recv.size()/4,4);
							collisionConstraints = Eigen::MatrixXf::Zero(constraints_recv.size()/4,4);
							// cout << "collisionConstraints.rows:"  << collisionConstraints.rows() << endl;
							// cout << "counter: " << counter << endl;
							// int j = -1;
							// Iterate over the number of constraints
							for (int i = 0; i < constraints_recv.size()/4; i++)
							{						
								// if (boost::lexical_cast<float>(constraints_recv[4*i+2]) < 0)
								// {
								// 	continue;
								// }
								// j++;
								// cout << "constraints_recv[" << 4*i << "]: " << constraints_recv[4*i] << endl;
								tempcollisionConstraints(i,0) = boost::lexical_cast<float>(constraints_recv[4*i]);
								// cout << "constraints_recv[" << 4*i+1 << "]: " << constraints_recv[4*i+1] << endl;
								tempcollisionConstraints(i,1) = boost::lexical_cast<float>(constraints_recv[4*i+1]);
								// cout << "constraints_recv[" << 4*i+2 << "]: " << constraints_recv[4*i+2] << endl;
								tempcollisionConstraints(i,2) = boost::lexical_cast<float>(constraints_recv[4*i+2]);
								tempcollisionConstraints(i,3) = 20+boost::lexical_cast<float>(constraints_recv[4*i+3]);
							}

							if (tempcollisionConstraints.norm() != 0)
							{
								collisionConstraints = tempcollisionConstraints;
							}


						}
						constraintPass = true;

					pthread_mutex_unlock(&constraint_lock);
			
					counter = 0;

				}


			}
			
			current_time_constraints = std::chrono::high_resolution_clock::now();

		}

		///////////////////////////////
		// TRAJECTORY COMMUNICATIONS //
		///////////////////////////////

		if (std::chrono::duration_cast<std::chrono::microseconds>(timeElapsed_trajectory).count() >= trajectory_time_us)
		{

			// SEND THE CONSTRAINTS // 
			if (client.socket_active[11])
			{
				
				// Make some test data and append to string buffer
				trajectory_send.clear();
				trajectory_send.append("T");
				trajectory_send.append(",");
				while(pthread_mutex_trylock(&trajectory_lock))
				{
					usleep(1);
				}
				pthread_mutex_unlock(&trajectory_lock);
				pthread_mutex_lock(&trajectory_lock);				
					for (int j = 0; j < mpc_params.nu_X; j++)
					{
						for (int i = 0; i < mpc_params.T; i++)
						{
							for (int k = 0; k < quadrotor.n; k++)
							{					
			
								float temp = full_trajectory[j](k,i);
								// full_trajectory_interf[j](k,i) = full_trajectory[j](k,i);
								trajectory_send.append( boost::lexical_cast<string>(temp) );
								if (j == mpc_params.nu_X-1 && i == mpc_params.T-1 && k == quadrotor.n-1)
								{
									trajectory_send.append("!");
								}
								else
								{
									trajectory_send.append(",");
								}
							}
						}
					}
				pthread_mutex_unlock(&trajectory_lock);

				// cout << "trajectory_send: " << trajectory_send << endl;

				// Copy the string to a char buffer
				// cstr_trajectory = new char[1000];
				strcpy(cstr_trajectory, trajectory_send.c_str());
				memcpy(trajectory_send_buffer, cstr_trajectory, strlen(cstr_trajectory)+1);

				// Send the path char buffer					
				client.sockets[11]->process_sending(trajectory_send_buffer,6000*sizeof(char));

			}
			
			current_time_trajectory = std::chrono::high_resolution_clock::now();

		}


		if (constraintPass && pathPass && mapPass && posePass)
		{
			firstPassComplete = true;
		}

		// Compute how much time has passed
		timeElapsed_map = end_time - current_time_map;
		timeElapsed_pose = end_time - current_time_pose;
		timeElapsed_path = end_time - current_time_path;
		timeElapsed_constraints = end_time - current_time_constraints;
		timeElapsed_trajectory = end_time - current_time_trajectory;

		// if (!(_interface_loops % 400))
		// {

		// 	cout << "\033[2J\033[1;1H";
		// 	cout << "Position: " << std::fixed << std::setprecision(5) << -pose[0] << ", " << -pose[1] << ", " << -pose[2] << " [m] Heading: " << pose[14] << " [rad]" << endl;
		// 	for (int i = 0; i < 6; i++)
		// 	{
		// 		// cout << "full_trajectory_interf " << i << " size: " << full_trajectory_interf[i].rows() << "x" << full_trajectory_interf[i].cols() << endl; 
		// 		if (segmentGoals.rows() > 0)
		// 		{
		// 			cout << "Traj " << std::fixed << std::setprecision(5) << i << " - Goal: " << segmentGoals.block(numTrajectoryPlans-1,4*i,1,4) << endl << full_trajectory_interf[i].block(0,0,3,5) << endl;
		// 		}
		// 		else
		// 		{
		// 			cout << "Traj " << std::fixed << std::setprecision(5) << i << endl << full_trajectory_interf[i].block(0,0,3,5) << endl;
		// 		}
		// 	}
		// 	cout << "Timers:   " << "  MAP  " << "  POS  " << "  PAT  " << "  CON  " << "  TRA  " << endl;
		// 	cout << "Time [s]: " << std::fixed << std::setprecision(5) << timeElapsed_map.count() << " " << timeElapsed_pose.count() << " " << timeElapsed_path.count() << " " << timeElapsed_constraints.count() << " " << timeElapsed_trajectory.count() << endl;

		// }


	}

}
#endif