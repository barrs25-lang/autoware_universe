#ifndef F_MPC_UNCUT_H
#define F_MPC_UNCUT_H


#include <read_files.h>
#include <logging.h>
#include <Eigen/Dense>
#include <structures.h>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#include <complex>
#include <cmath>
#include <algorithm>

#include <time.h>
#include <chrono>
#include <sys/time.h>
#include <ctime>

#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/interval.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <sys/ioctl.h>


using namespace std;
using namespace boost::numeric::odeint;
using namespace std::complex_literals;
    
const int rot_dyn_states = 8;
typedef boost::array< double, 8> rot_dyn_state; // phi,theta,psi,omega,u1,u1d

class F_MPC_UNCUT
{

public:

	F_MPC_UNCUT();
	~F_MPC_UNCUT();

	void start();

	void trajectory_planner_thread();
	void communication_thread();
	System quadrotor;
	MPC_Params mpc_params;
	Matrix_Set bomt;
	Line_Search_Params lsp;

	// Goal
	Eigen::MatrixXf goal;
	// Copy of the planned path
	Eigen::MatrixXf localPath;

	int numTrajectoryPlans = 1;

private:

	static void quit_handler(int sig);
	void verifyDimensions() const;
	void get_Omega();
	void compute_Gamma();
	void compute_Gamma(float phi, float theta, float psi, Eigen::MatrixXf* Gamma_k);
	void compute_dEuler();
	void find_closest_waypoint(float* xdiff, float* ydiff, float* zdiff, int* closest_traj_iterator);
	void find_new_waypoint(Eigen::MatrixXf* goal, float* xdiff, float* ydiff, float* zdiff);
	float fsat(Eigen::MatrixXf* arg);
	void find_closest_obstacle(Eigen::MatrixXf* temp_goal_pos, int segment_number);
	void compute_Euler_angles(Eigen::Matrix3f* R, Eigen::Vector3f* Eul);
	void update_weighting_matrices(int segment_number, bool success);
	void update_and_discount_weighting_matrices(int segment_number, bool success);
	bool compute_roll_pitch_thrust_for_lambda_k_and_g_barrier(Eigen::MatrixXf* chi, Eigen::MatrixXf* lambda);
	void get_Ginv_k(float phi, float theta, float psi, float u1, Eigen::MatrixXf* G_inv);	
	void get_f_k(float phi, float theta, float psi, Eigen::MatrixXf* omega, float u1, float u1_dot, Eigen::MatrixXf* f, int i);
	void compute_u1_u1dot(Eigen::MatrixXf* eig_U);
	void compute_box_constraints(int segment_number, int goalStride);
	bool fmpcsolve(Eigen::MatrixXf* X, Eigen::MatrixXf* U, int segment_number);
	void fmpc_hard_constraints(Eigen::MatrixXf* eig_X, int segment_number, int goalStride, int numFailures);
	void fmpc_soft_constraints(Eigen::MatrixXf* eig_X, int segment_number, int goalStride);
	static void sys(const rot_dyn_state &y, rot_dyn_state &ydot, double t);
	void update_quadrotor_pose();
	void get_trajectory_goal();
	void objective_function_value(Eigen::MatrixXf* X, Eigen::MatrixXf* U, int segment_number);
	bool validates_collision_constraints();

	pthread_t traj_tid;
	pthread_t comm_tid;

	pthread_mutex_t pose_lock;
	pthread_mutex_t path_lock;
	pthread_mutex_t map_lock;
	pthread_mutex_t constraint_lock;
	pthread_mutex_t trajectory_lock;
	pthread_mutex_t ellipsoid_lock;

	// Stores constraints sent from constraint generator
	Eigen::MatrixXf collisionConstraints;

	// Stores the parent ellipsoid from the collision avoidance algorithm
	Eigen::MatrixXf parentEllipsoid;

	// Stores the planned path from the path planner
	Eigen::MatrixXf plannedPath, prev_plannedPath;

	// Stores differences between consecutive planned paths
	Eigen::MatrixXf plannedPathDifference;

	// Stores the pose from the simulator
	float pose[18], initialPose[18]; // x,y,z,dx,dy,dz,ddx,ddy,ddz,dddx,dddy,dddz,phi,theta,psi,omega1,omega2,omega3

	// Gamma for rotational kinematic equation
	Eigen::MatrixXf Gamma;

	// Gamma for simulating the rotational kinematics using zeta_k
	Eigen::MatrixXf Gamma_k;

	// Angular velocity of body frame
	Eigen::MatrixXf omega;

	// Angular velocity of body frame simulated using zeta_k
	Eigen::MatrixXf omega_k;

	// Time derivative of euler angles
	Eigen::MatrixXf dEuler;

	// Obstacle coordinates
	Eigen::MatrixXf obs;
	Eigen::MatrixXf local_obs;

	// Map
	float voxel_map[100][30][100];

	// Nearest obstacle
	Eigen::MatrixXf r_obs;

	bool firstPassComplete = false;
	
	// Keeps track if the first planning episode was ever successful
	bool firstPlanningEpisode = false;

	// Keeps track of the success of planning segments
	bool trajectory_plan_success = false;

	// Keeps tracking of intermediate planning success
	bool intermediate_success = false;

	bool firstSuccess = false;

	// keeps track of nearest waypoint in the planned path
	int traj_iterator = 0;
	int prev_traj_iterator = -1;

	// Keeps track of how many times planning a segment failed
	int numFailures = 0;

	// Boolean indicating the unfortunate failure to plan
	bool failedPlan = false;

	// keeps track of the segments needed to be planned to interpolate the planned path 
	int remaining_segments = 0;

	// Barrier function used to ensure phi and theta remain in the intervals (-phi_max, phi_max) and (-theta_max, theta_max), resp.; and T_p > 0.
	Eigen::MatrixXf g_barrier;

	int newTraj = 0;

	// Matrix storing the full remaining trajectory which outlines the planned path
	Eigen::MatrixXf* full_trajectory;
	Eigen::MatrixXf* sendfull_trajectory;
	Eigen::MatrixXf* prevfull_trajectory;
	Eigen::MatrixXf* full_trajectory_interf;
	Eigen::MatrixXf segmentGoals;
	Eigen::MatrixXf* full_policy;
	Eigen::MatrixXf* full_attitude;
	Eigen::MatrixXf* prev_seg_eig_X;
	Eigen::MatrixXf* prev_seg_eig_U;

	Eigen::MatrixXf* u_k_nuX;
	Eigen::MatrixXf* v_k_nuX;
	Eigen::MatrixXf* du1_k_nuX;
	Eigen::MatrixXf* ddu1_k_nuX;
	Eigen::MatrixXf* lambda_k_nuX;
	Eigen::MatrixXf* zeta_k_nuX;
	Eigen::MatrixXf* T_k_nuX;

	Eigen::MatrixXf u_k;
	Eigen::MatrixXf v_k;
	Eigen::MatrixXf du1_k;
	Eigen::MatrixXf ddu1_k;
	Eigen::MatrixXf lambda_k;
	Eigen::MatrixXf zeta_k;

	// Matrix (will be sized as scalar) storing the objective function value
	Eigen::MatrixXf Obj;

};

#endif