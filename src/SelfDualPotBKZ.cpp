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

// External libraries
#include <eigen3/Eigen/Dense>
#include <NTL/RR.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>

#include "Lattice.h"

extern "C" long **PotLLL(long **basis, const double reduction_parameter, const int n, const int m)
{
    int i, j;
    Lattice B;
    B.setDims(n, m);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
        {
            B.basis.coeffRef(i, j) = basis[i][j];
        }

    B.PotLLL_(reduction_parameter, n, m);
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
    B.DualPotLLL_(reduction_parameter, n, m);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            basis[i][j] = B.basis.coeff(i, j);
    return basis;
}

extern "C" long **PotBKZ(long **basis, const int beta, const double reduction_parameter, const int n, const int m)
{
    FILE *_Pot_FILE_ = fopen(".data/potential_graph.csv", "wt");
    fprintf(_Pot_FILE_, "Potential\n");
    int i, j;
    Lattice B;
    B.setDims(n, m);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            B.basis.coeffRef(i, j) = basis[i][j];
    B.PotBKZ_(beta, reduction_parameter, n, m, _Pot_FILE_);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            basis[i][j] = B.basis.coeff(i, j);
    fclose(_Pot_FILE_);
    return basis;
}

extern "C" long **DualPotBKZ(long **basis, const int beta, const double delta, const int n, const int m)
{
    FILE *_Dual_Pot_FILE_ = fopen(".data/dual_potential_graph.csv", "wt");
    int i, j;
    Lattice B;
    B.setDims(n, m);
    fprintf(_Dual_Pot_FILE_, "Potential\n");

    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            B.basis.coeffRef(i, j) = basis[i][j];
    B.DualPotBKZ_(beta, delta, n, m, _Dual_Pot_FILE_);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            basis[i][j] = B.basis.coeff(i, j);
    fclose(_Dual_Pot_FILE_);
    return basis;
}

extern "C" long **SelfDualPotBKZ(long **basis, const int beta, const double reduction_parameter, const int n, const int m)
{
    FILE *_Dual_Primal_ = fopen(".data/potential_graph_DualENUM_and_PrimalLLL.csv", "wt");
    fprintf(_Dual_Primal_, "Potential\n");
    int i, j;
    Lattice B;
    B.setDims(n, m);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            B.basis.coeffRef(i, j) = basis[i][j];
    B.SelfDualPotBKZ_each_(beta, reduction_parameter, n, m, _Dual_Primal_);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            basis[i][j] = B.basis.coeff(i, j);
    fclose(_Dual_Primal_);
    return basis;
}
