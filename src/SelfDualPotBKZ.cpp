/**
 * @file SelfDualPotBKZ.cpp
 * @author 佐藤 新 (23lc002y@rikkyo.ac.jp)
 * @brief
 * @version 0.1
 * @date 2025-01-21
 *
 * @copyright Copyright (c) 2025
 *
 */

#include <iostream>
#include <random>
#include <algorithm>
#include <vector>
#include <cmath>

#include <eigen3/Eigen/Dense>
#include <NTL/RR.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>

#include "Lattice.h"
#include "BKZ.h"
#include "PotLLL.h"
#include "DualPotLLL.h"
#include "PotBKZ.h"
#include "DualPotBKZ.h"
#include "SelfDualPotBKZ.h"

extern "C" long **PotLLL(long **basis, const double reduction_parameter, const int n, const int m)
{
    int i, j;
    Lattice B;
    B.setDims(n, m);
    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < m; ++j)
        {
            B.basis.coeffRef(i, j) = basis[i][j];
        }
    }

    B.PotLLL_(reduction_parameter, n, m);

    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < m; ++j)
        {
            basis[i][j] = B.basis.coeff(i, j);
        }
    }
    return basis;
}

extern "C" long **BKZ(long **basis, const int block_size, const double reduction_parameter, const int max_loop, const int n, const int m)
{
    FILE *potential_of_bkz = fopen(".data/potential_of_BKZ.csv", "wt");
    fprintf(potential_of_bkz, "Potential\n");
    int i, j;
    Lattice B;
    B.setDims(n, m);
    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < m; ++j)
        {
            B.basis.coeffRef(i, j) = basis[i][j];
        }
    }

    B.BKZ_(block_size, reduction_parameter, max_loop, n, m, potential_of_bkz);

    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < m; ++j)
        {
            basis[i][j] = B.basis.coeff(i, j);
        }
    }
    fclose(potential_of_bkz);
    return basis;
}

extern "C" long **DualPotLLL(long **basis, const double reduction_parameter, const int n, const int m)
{
    int i, j;
    Lattice B;
    B.setDims(n, m);
    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < m; ++j)
        {
            B.basis.coeffRef(i, j) = basis[i][j];
        }
    }

    B.DualPotLLL_(reduction_parameter, n, m);

    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < m; ++j)
        {
            basis[i][j] = B.basis.coeff(i, j);
        }
    }
    return basis;
}

extern "C" long **PotBKZ(long **basis, const int beta, const double reduction_parameter, const int n, const int m)
{
    FILE *potential_of_pot_bkz = fopen(".data/potential_of_PotBKZ.csv", "wt");
    fprintf(potential_of_pot_bkz, "Potential\n");
    int i, j;
    Lattice B;
    B.setDims(n, m);

    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < m; ++j)
        {
            B.basis.coeffRef(i, j) = basis[i][j];
        }
    }

    B.PotBKZ_(beta, reduction_parameter, n, m, potential_of_pot_bkz);

    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < m; ++j)
        {
            basis[i][j] = B.basis.coeff(i, j);
        }
    }

    fclose(potential_of_pot_bkz);
    return basis;
}

extern "C" long **DualPotBKZ(long **basis, const int beta, const double delta, const int n, const int m)
{
    FILE *potential_of_dual_pot_bkz = fopen(".data/potential_of_DualPotBKZ.csv", "wt");
    int i, j;
    Lattice B;
    B.setDims(n, m);
    fprintf(potential_of_dual_pot_bkz, "Potential\n");

    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < m; ++j)
        {
            B.basis.coeffRef(i, j) = basis[i][j];
        }
    }

    B.DualPotBKZ_(beta, delta, n, m, potential_of_dual_pot_bkz);

    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < m; ++j)
        {
            basis[i][j] = B.basis.coeff(i, j);
        }
    }

    fclose(potential_of_dual_pot_bkz);
    return basis;
}

extern "C" long **SelfDualPotBKZ(long **basis, const int beta, const double reduction_parameter, const int n, const int m)
{
    FILE *potential_of_self_dual_pot_bkz = fopen(".data/potential_of_SelfDualPotBKZ.csv", "wt");
    fprintf(potential_of_self_dual_pot_bkz, "Potential\n");
    int i, j;
    Lattice B;
    B.setDims(n, m);
    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < m; ++j)
        {
            B.basis.coeffRef(i, j) = basis[i][j];
        }
    }

    B.SelfDualPotBKZ_(beta, reduction_parameter, n, m, potential_of_self_dual_pot_bkz);

    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < m; ++j)
        {
            basis[i][j] = B.basis.coeff(i, j);
        }
    }

    fclose(potential_of_self_dual_pot_bkz);
    return basis;
}
