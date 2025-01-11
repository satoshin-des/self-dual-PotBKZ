#ifndef SDPB_HPP
#define SDPB_HPP

#include <iostream>
#include <eigen3/Eigen/Dense>

typedef Eigen::Matrix<long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXli;
typedef Eigen::Matrix<long, 1, Eigen::Dynamic> VectorXli;
typedef Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXld;
typedef Eigen::Matrix<long double, 1, Eigen::Dynamic> VectorXld;

void DualDeepInsertion(MatrixXli &basis, const int n, const int k, const int l);
MatrixXli Insert(const MatrixXli basis, const VectorXli x, const int n, const int m);
void GSO(const MatrixXli basis, VectorXld &squared_norm_of_gso_vec, MatrixXld &gso_coeff_mat, const int n, const int m);
void GSO(const MatrixXld basis, VectorXld &squared_norm_of_gso_vec, MatrixXld &gso_coeff_mat, const int n, const int m);
void GSO(const MatrixXli basis, VectorXld &squared_norm_of_gso_vec, VectorXld &log_squared_norm_of_gso_vec, MatrixXld &gso_coeff_mat, const int n, const int m);
void DualGSO(const VectorXld squared_norm_of_gso_vec, const VectorXld log_squared_norm_of_gso_vec, const MatrixXld gso_coeff_mat, VectorXld &dual_squared_norm_of_gso_vec, VectorXld &dual_log_squared_norm_of_gso_vec, MatrixXld &dual_gso_coeff_mat, const int n, const int m);
VectorXli ENUM(const MatrixXld gso_coeff_mat, const VectorXld squared_norm_of_gso_vec, VectorXld &rho, const int n, const double enumeration_upper_bound);
VectorXli enumerate(const MatrixXld gso_coeff_mat, const VectorXld squared_norm_of_gso_vec, VectorXld &rho, const int n);
VectorXli PotENUM(const MatrixXld gso_coeff_mat, const VectorXld squared_norm_of_gso_vec, const VectorXld log_squared_norm_of_gso_vec, const int n);
VectorXli DualPotENUM(const MatrixXld gso_coeff_mat, const VectorXld squared_norm_of_gso_vec, const VectorXld log_squared_norm_of_gso_vec, const int n);
void POT_LLL(MatrixXli &basis, const long double reduction_parameter, const int n, const int m);
void DUAL_POT_LLL(MatrixXli &basis, const double reduction_parameter, const int n, const int m);
void POT_BKZ(MatrixXli &basis, const int beta, const double reduction_parameter, const int n, const int m);
void DUAL_POT_BKZ(MatrixXli &basis, const int beta, const double delta, const int n, const int m);
void SELF_DUAL_POT_BKZ(MatrixXli &basis, const int beta, const double reduction_parameter, const int n, const int m);

extern "C" long **PotLLL(long **basis, const double reduction_parameter, const int n, const int m);
extern "C" long **DualPotLLL(long **basis, const double reduction_parameter, const int n, const int m);
extern "C" long **PotBKZ(long **basis, const int beta, const double reduction_parameter, const int n, const int m);
extern "C" long **DualPotBKZ(long **basis, const int beta, const double delta, const int n, const int m);
extern "C" long **SelfDualPotBKZ(long **basis, const int beta, const double reduction_parameter, const int n, const int m);

#endif // !SDPB_HPP
