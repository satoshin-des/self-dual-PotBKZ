#ifndef LATTICE
#define LATTICE

#include <iostream>
#include <eigen3/Eigen/Dense>

typedef Eigen::Matrix<long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXli;
typedef Eigen::Matrix<long, 1, Eigen::Dynamic> VectorXli;
typedef Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXld;
typedef Eigen::Matrix<long double, 1, Eigen::Dynamic> VectorXld;

void DualDeepInsertion(MatrixXli& b, const int n, const int k, const int l);
MatrixXli Insert(const MatrixXli b, const VectorXli x, const int n, const int m);
MatrixXld dual(const MatrixXli b, const int n, const int m);
void Randomization(MatrixXli& b, const int n, const int m);
void GSO(const MatrixXli b, VectorXld& B, MatrixXld& mu, const int n, const int m);
void GSO(const MatrixXld b, VectorXld& B, MatrixXld& mu, const int n, const int m);
void GSO(const MatrixXli b, VectorXld& B, VectorXld& logB, MatrixXld& mu, const int n, const int m);
void DualGSO(const VectorXld B, const MatrixXld mu, VectorXld& C, MatrixXld& hmu, const int n, const int m);
long double logPot(const VectorXld B, const int n);
long double rho(const VectorXld B, const int n, const int m);
VectorXli ENUM(const MatrixXld mu, const VectorXld B, VectorXld& rho, const int n, const double R);
VectorXli enumerate(const MatrixXld mu, const VectorXld B, VectorXld& rho, const int n);
VectorXli PotENUM(const MatrixXld mu, const VectorXld B, const VectorXld logB, const int n);
VectorXli DualPotENUM(const MatrixXld mu, const VectorXld B, const VectorXld logB, const int n);
void __PotLLL__(MatrixXli& b, const long double d, const int n, const int m);
void __DualPotLLL__(MatrixXli& b, const double d, const int n, const int m);
void __BKZ__(MatrixXli& b, const int beta, const double d, const int lp, const int n, const int m);
void __PotBKZ__(MatrixXli& b, const int beta, const double d, const int n, const int m);
void __DualPotBKZ__(MatrixXli& b, const int beta, const double delta, const int n, const int m);
void __SelfDualPotBKZ__(MatrixXli& b, const int beta, const double d, const int n, const int m);

extern "C" long **PotLLL(long **b, const double d, const int n, const int m);
extern "C" long **DualPotLLL(long **b, const double d, const int n, const int m);
extern "C" long **BKZ(long **b, const int beta, const double d, const int lp, const int n, const int m);
extern "C" long **PotBKZ(long **b, const int beta, const double d, const int n, const int m);
extern "C" long **DualPotBKZ(long **b, const int beta, const double delta, const int n, const int m);
extern "C" long **SelfDualPotBKZ(long **b, const int beta, const double d, const int n, const int m);

#endif // !LATTICE
