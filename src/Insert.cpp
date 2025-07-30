#include "Lattice.h"

#include <iostream>
#include <cmath>

#include <eigen3/Eigen/Dense>

#include <NTL/ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>

MatrixXli Lattice::Insert(const MatrixXli b, const VectorXli x, const int n, const int m)
{
    int i, j;
    const long double beta = 1.35135135135135135135135135135135135135135135135135;
    long double gamma;
    MatrixXli U = MatrixXli::Zero(n, n);
    NTL::mat_ZZ temp_basis;
    temp_basis.SetDims(n, n + 1);

    /* Construction of gamma */
    m_temp = x.cast<double>().norm();
    m_temp *= pow(beta, (n - 2) >> 1);
    gamma = round(m_temp + m_temp);

    /* Construction of matrix */
    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < n; ++j)
        {
            temp_basis[i][j] = 0;
        }
        temp_basis[i][i] = 1;
        temp_basis[i][n] = gamma * x.coeff(i);
    }

    NTL::LLL_FP(temp_basis, 0.99);

    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < n; ++j)
        {
            U.coeffRef(i, j) = NTL::to_long(temp_basis[i][j]);
        }
    }

    return U * b;
}
