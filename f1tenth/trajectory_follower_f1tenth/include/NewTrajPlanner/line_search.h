#ifndef LINE_SEARCH_H
#define LINE_SEARCH_H

#include <f_mpc_uncut.h>

float objective_function_value(F_MPC_UNCUT* fmpc, Eigen::MatrixXf* X, Eigen::MatrixXf* U, Eigen::MatrixXf* dX, Eigen::MatrixXf* dU, float s, int segment_number);
float rd_tilde_rp_norm(F_MPC_UNCUT* fmpc);
float dichotomous_line_search(F_MPC_UNCUT* fmpc, Eigen::MatrixXf* X, Eigen::MatrixXf* U, Eigen::MatrixXf* dX, Eigen::MatrixXf* dU, float a1, float b1, float l, float epsilon, int segment_number, float start_s);
float newton_method_line_search();

#endif