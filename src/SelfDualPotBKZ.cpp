/**
 * @file PotBKZ.cpp
 * @author Arata Sato (23lc002y@rikkyo.ac.jp)
 * @brief
 * @version 0.1
 * @date 2024-10-01
 *
 * @copyright Copyright (c) 2024
 *
 */

#include <iostream>
#include <random>
#include <algorithm>
#include <vector>
#include <cmath>

#include "Lattice.h"
#include "PotENUM.h"
#include "DualPotENUM.h"
#include "PotLLL.h"
#include "DualPotLLL.h"
#include "PotBKZ.h"
#include "DualPotBKZ.h"
#include "SelfDualPotBKZ.h"

extern "C" long **PotLLL(long **basis, const double reduction_parameter, const int n, const int m)
{
    Lattice B;
    B.setDims(n, m);
    int i, j;
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
        {
            B.basis.coeffRef(i, j) = basis[i][j];
        }

    B.POT_LLL(reduction_parameter);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            basis[i][j] = B.basis.coeff(i, j);
    return basis;
}

extern "C" long **DualPotLLL(long **basis, const double reduction_parameter, const int n, const int m)
{
    int i, j;
    Lattice B;
    B.setDims(n, m);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            B.basis.coeffRef(i, j) = basis[i][j];
    B.DUAL_POT_LLL(reduction_parameter);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            basis[i][j] = B.basis.coeff(i, j);
    return basis;
}

extern "C" long **PotBKZ(long **basis, const int beta, const double reduction_parameter, const int n, const int m)
{
    int i, j;
    Lattice B;
    B.setDims(n, m);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            B.basis.coeffRef(i, j) = basis[i][j];
    B.POT_BKZ(beta, reduction_parameter);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            basis[i][j] = B.basis.coeff(i, j);
    return basis;
}

extern "C" long **DualPotBKZ(long **basis, const int beta, const double delta, const int n, const int m)
{
    int i, j;
    Lattice B;
    B.setDims(n, m);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            B.basis.coeffRef(i, j) = basis[i][j];
    B.DUAL_POT_BKZ(beta, delta);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            basis[i][j] = B.basis.coeff(i, j);
    return basis;
}

extern "C" long **SelfDualPotBKZ(long **basis, const int beta, const double reduction_parameter, const int n, const int m)
{
    int i, j;
    Lattice B;
    B.setDims(n, m);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            B.basis.coeffRef(i, j) = basis[i][j];
    B.SELF_DUAL_POT_BKZ(beta, reduction_parameter);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            basis[i][j] = B.basis.coeff(i, j);
    return basis;
}
