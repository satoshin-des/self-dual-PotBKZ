#include "Lattice.h"

#include <iostream>

#include <eigen3/Eigen/Dense>

void Lattice::GSO(VectorXld &B, MatrixXld &mu, const int n, const int m)
{
    MatrixXld gso_mat(n, m);

    for (int i = 0, j; i < n; ++i)
    {
        mu.coeffRef(i, i) = 1.0;
        gso_mat.row(i) = basis.row(i).cast<long double>();
        for (j = 0; j < i; ++j)
        {
            mu.coeffRef(i, j) = basis.row(i).cast<long double>().dot(gso_mat.row(j)) / gso_mat.row(j).dot(gso_mat.row(j));
            gso_mat.row(i) -= mu.coeff(i, j) * gso_mat.row(j);
        }
        B.coeffRef(i) = gso_mat.row(i).dot(gso_mat.row(i));
    }
}

void Lattice::GSO(VectorXld &B, VectorXld &logB, MatrixXld &mu, const int n, const int m)
{
    MatrixXld gso_mat(n, m);

    for (int i = 0, j; i < n; ++i)
    {
        mu.coeffRef(i, i) = 1.0;
        gso_mat.row(i) = basis.row(i).cast<long double>();
        for (j = 0; j < i; ++j)
        {
            mu.coeffRef(i, j) = basis.row(i).cast<long double>().dot(gso_mat.row(j)) / gso_mat.row(j).dot(gso_mat.row(j));
            gso_mat.row(i) -= mu.coeff(i, j) * gso_mat.row(j);
        }
        B.coeffRef(i) = gso_mat.row(i).dot(gso_mat.row(i));
        logB.coeffRef(i) = log(B.coeff(i));
    }
}

void Lattice::DualGSO(const VectorXld B, const VectorXld logB, const MatrixXld mu, VectorXld &C, VectorXld &logC, MatrixXld &hmu, const int n, const int m)
{
    for (int i = 0, j; i < n; ++i)
    {
        C.coeffRef(i) = 1.0 / B.coeff(i);
        logC.coeffRef(i) = -logB.coeff(i);
        for (j = i + 1; j < m; ++j)
        {
            hmu.coeffRef(i, j) = -mu.row(j).segment(i, j - i).dot(hmu.row(i).segment(i, j - i));
        }
    }
}
