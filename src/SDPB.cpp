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
#include <omp.h>

// External libraries
#include <eigen3/Eigen/Dense>
#include <NTL/RR.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>

#include "SDPB.hpp"

#define LOG099 -0.010050335853501441183548857558547706085515007674629873378

static int PotTour = 0, BKZTour = 0, Tr = 0;
static NTL::ZZ _;

/// @brief Dual version of deep-insertion
/// @param basis Lattice basis
/// @param n Rank of lattice
/// @param k
/// @param l
void DualDeepInsertion(MatrixXli &basis, const int n, const int k, const int l)
{
    const VectorXli t = basis.row(k);
    for (int j = k; j < l; ++j)
        basis.row(j) = basis.row(j + 1);
    basis.row(l) = t;
}

/// @brief P.239
/// @param NewBasis
/// @param Basis
/// @param x
MatrixXli Insert(const MatrixXli basis, const VectorXli x, const int n, const int m)
{
    int i, j;
    double beta = 1.35135135135135135135135135135135135135135135135135;
    double tmp, gamma;
    MatrixXli U(n, n), c;
    U.setZero();
    NTL::mat_ZZ tmp_base;
    tmp_base.SetDims(n, n + 1);

    /* Construction of gamma */
    tmp = x.cast<double>().norm();
    tmp *= pow(beta, (n - 2) * 0.5);
    gamma = round(tmp + tmp);

    /* Construction of matrix */
    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < n; ++j)
            tmp_base[i][j] = 0;
        tmp_base[i][i] = 1;
        tmp_base[i][n] = gamma * x.coeff(i);
    }
    NTL::LLL(_, tmp_base, 99, 100);

    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j)
            U.coeffRef(i, j) = NTL::to_long(tmp_base[i][j]);

    /* Construction of a new basis */
    return U * basis;
}

void GSO(const MatrixXli basis, VectorXld &squared_norm_of_gso_vec, MatrixXld &gso_coeff_mat, const int n, const int m)
{
    MatrixXld gso_vec_mat(n, m);

    for (int i = 0, j; i < n; ++i)
    {
        gso_coeff_mat.coeffRef(i, i) = 1.0;
        gso_vec_mat.row(i) = basis.row(i).cast<long double>();
        for (j = 0; j < i; ++j)
        {
            gso_coeff_mat.coeffRef(i, j) = basis.row(i).cast<long double>().dot(gso_vec_mat.row(j)) / gso_vec_mat.row(j).dot(gso_vec_mat.row(j));
            gso_vec_mat.row(i) -= gso_coeff_mat.coeff(i, j) * gso_vec_mat.row(j);
        }
        squared_norm_of_gso_vec.coeffRef(i) = gso_vec_mat.row(i).dot(gso_vec_mat.row(i));
    }
}

void GSO(const MatrixXld basis, VectorXld &squared_norm_of_gso_vec, MatrixXld &gso_coeff_mat, const int n, const int m)
{
    MatrixXld gso_vec_mat(n, m);

    for (int i = 0, j; i < n; ++i)
    {
        gso_coeff_mat.coeffRef(i, i) = 1.0;
        gso_vec_mat.row(i) = basis.row(i);
        for (j = 0; j < i; ++j)
        {
            gso_coeff_mat.coeffRef(i, j) = basis.row(i).dot(gso_vec_mat.row(j)) / gso_vec_mat.row(j).dot(gso_vec_mat.row(j));
            gso_vec_mat.row(i) -= gso_coeff_mat.coeff(i, j) * gso_vec_mat.row(j);
        }
        squared_norm_of_gso_vec.coeffRef(i) = gso_vec_mat.row(i).dot(gso_vec_mat.row(i));
    }
}

void GSO(const MatrixXli basis, VectorXld &squared_norm_of_gso_vec, VectorXld &log_squared_norm_of_gso_vec, MatrixXld &gso_coeff_mat, const int n, const int m)
{
    MatrixXld gso_vec_mat(n, m);

    for (int i = 0, j; i < n; ++i)
    {
        gso_coeff_mat.coeffRef(i, i) = 1.0;
        gso_vec_mat.row(i) = basis.row(i).cast<long double>();
        for (j = 0; j < i; ++j)
        {
            gso_coeff_mat.coeffRef(i, j) = basis.row(i).cast<long double>().dot(gso_vec_mat.row(j)) / gso_vec_mat.row(j).dot(gso_vec_mat.row(j));
            gso_vec_mat.row(i) -= gso_coeff_mat.coeff(i, j) * gso_vec_mat.row(j);
        }
        squared_norm_of_gso_vec.coeffRef(i) = gso_vec_mat.row(i).dot(gso_vec_mat.row(i));
        log_squared_norm_of_gso_vec.coeffRef(i) = log(squared_norm_of_gso_vec.coeff(i));
    }
}

/// @brief Computes GSO-informations of dual-lattice basis using GSO-information of lattice
/// @param squared_norm_of_gso_vec Squared norms of GSO-vectors
/// @param gso_coeff_mat GSO-coefficient mattix
/// @param dual_squared_norm_of_gso_vec Squared norms of dual-GSO-vectors
/// @param hmu dual-GSO-coefficient mattix
/// @param n
/// @param m
void DualGSO(const VectorXld squared_norm_of_gso_vec, const VectorXld log_squared_norm_of_gso_vec, const MatrixXld gso_coeff_mat, VectorXld &dual_squared_norm_of_gso_vec, VectorXld &logC, MatrixXld &hmu, const int n, const int m)
{
    for (int i = 0, j; i < n; ++i)
    {
        dual_squared_norm_of_gso_vec.coeffRef(i) = 1.0 / squared_norm_of_gso_vec.coeff(i);
        logC.coeffRef(i) = -log_squared_norm_of_gso_vec.coeff(i);
        for (j = i + 1; j < m; ++j)
        {
            hmu.coeffRef(i, j) = -gso_coeff_mat.row(j).segment(i, j - i).dot(hmu.row(i).segment(i, j - i));
        }
    }
}

/// @brief Enumerates a vector whose norm is shorter than R on the lattice
/// @param gso_coeff_mat GSO-coefficient matrix
/// @param squared_norm_of_gso_vec Squared-norms of row vectors of GSO-matrix
/// @param n Rank of lattice
/// @param R An upper bound of enumeration
/// @return
VectorXli ENUM(const MatrixXld gso_coeff_mat, const VectorXld squared_norm_of_gso_vec, VectorXld &rho, const int n, const double R)
{
    const int n1 = n + 1;
    int i, r[n1];
    double tmp;
    VectorXli w(n), v(n);
    w.setZero();
    v.setZero();
    v.coeffRef(0) = 1; // w: ジグザグに動く変分
    Eigen::VectorXd c(n);
    c.setZero();
    rho.setZero();
    Eigen::MatrixXd sigma(n1, n);
    sigma.setZero();
    for (i = 0; i < n; ++i)
    {
        r[i] = i;
    }

    for (int k = 0, last_nonzero = 0;;)
    {
        tmp = (double)v.coeff(k) - c.coeff(k);
        tmp *= tmp;
        rho.coeffRef(k) = rho.coeff(k + 1) + tmp * squared_norm_of_gso_vec.coeff(k); // rho[k]=∥πₖ(v)∥
        if (rho.coeff(k) <= R)
        {
            if (k == 0)
            {
                return v;
            }
            else
            {
                --k;
                r[k] = (r[k] > r[k + 1] ? r[k] : r[k + 1]);
                for (i = r[k]; i > k; --i)
                {
                    sigma.coeffRef(i, k) = sigma.coeff(i + 1, k) + gso_coeff_mat.coeff(i, k) * v.coeff(i);
                }
                c.coeffRef(k) = -sigma.coeff(k + 1, k);
                v.coeffRef(k) = round(c.coeff(k));
                w.coeffRef(k) = 1; // 小さいやつが見つかったら、変分を元に戻す
            }
        }
        else
        {
            ++k;
            if (k == n)
            {
                v.setZero();
                return v;
            }
            else
            {
                r[k] = k;
                if (k >= last_nonzero)
                {
                    last_nonzero = k;
                    ++v.coeffRef(k);
                }
                else
                {
                    v.coeff(k) > c.coeff(k) ? v.coeffRef(k) -= w.coeff(k) : v.coeffRef(k) += w.coeff(k);
                    ++w.coeffRef(k);
                }
            }
        }
    }
}

/// @brief Enumerates the shortest vector on the lattice
/// @param gso_coeff_mat GSO-coefficient matrix
/// @param squared_norm_of_gso_vec Squared-norms of row vectors of GSO-matrix
/// @param n Rank of lattice
/// @return VectorXli the shortest vector
VectorXli enumerate(const MatrixXld gso_coeff_mat, const VectorXld squared_norm_of_gso_vec, VectorXld &rho, const int n)
{
    VectorXli enumerated_vector(n), old_enumerated_vector(n);
    enumerated_vector.setZero();
    old_enumerated_vector.setZero();
    VectorXld pre_rho(n + 1);
    pre_rho.setZero();
    for (double R = squared_norm_of_gso_vec.coeff(0);;)
    {
        pre_rho = rho;
        old_enumerated_vector = enumerated_vector;
        enumerated_vector = ENUM(gso_coeff_mat, squared_norm_of_gso_vec, rho, n, R);
        if (enumerated_vector.isZero())
        {
            rho = pre_rho;
            return old_enumerated_vector;
        }
        R *= 0.99;
    }
}

/// @brief Enumerates a short and small potential lattice vector
/// @param gso_coeff_mat GSO-coefficient matrix.
/// @param squared_norm_of_gso_vec Squared norms of GSO vectors.
/// @param log_squared_norm_of_gso_vec Logarithm valude of squared norms of GSO-vectors.
/// @param n Rank of lattice.
/// @return VectorXli
VectorXli PotENUM(const MatrixXld gso_coeff_mat, const VectorXld squared_norm_of_gso_vec, const VectorXld log_squared_norm_of_gso_vec, const int n)
{
    int i, r[n + 1];
    double tmp, R = log_squared_norm_of_gso_vec.coeff(0), P = 0;
    VectorXli w(n), v(n);
    w.setZero();
    v.setZero();
    Eigen::VectorXd c(n), D(n + 1);
    c.setZero();
    D.setZero();
    Eigen::MatrixXd sigma(n + 1, n);
    sigma.setZero();

    v.coeffRef(0) = 1;

    for (i = 0; i <= n; ++i)
        r[i] = i;
    for (int k = 0, last_nonzero = 0;;)
    {
        tmp = (double)v.coeff(k) - c.coeff(k);
        tmp *= tmp;
        D.coeffRef(k) = D.coeff(k + 1) + tmp * squared_norm_of_gso_vec.coeff(k);

        if ((k + 1) * log(D.coeff(k)) + P < (k + 1) * LOG099 + R)
        {
            if (k == 0)
                return v;
            else
            {
                P += log(D.coeff(k));
                --k;
                r[k] = (r[k] > r[k + 1] ? r[k] : r[k + 1]);
                for (i = r[k]; i > k; --i)
                    sigma.coeffRef(i, k) = sigma.coeff(i + 1, k) + gso_coeff_mat.coeff(i, k) * v(i);
                c.coeffRef(k) = -sigma.coeff(k + 1, k);
                v.coeffRef(k) = round(c.coeff(k));
                w.coeffRef(k) = 1;
            }
        }
        else
        {
            ++k;
            if (k == n)
            {
                v.setZero();
                return v;
            }
            else
            {
                r[k - 1] = k;
                if (k >= last_nonzero)
                {
                    last_nonzero = k;
                    ++v.coeffRef(k);
                    if (v.coeff(last_nonzero) >= 2)
                    {
                        ++k;
                        if (k == n)
                        {
                            v.setZero();
                            return v;
                        }
                        else
                        {
                            r[k - 1] = k;
                            last_nonzero = k;
                            v.coeffRef(last_nonzero) = 1;
                        }
                    }
                    P = 0;
                    R = log_squared_norm_of_gso_vec.head(last_nonzero + 1).array().sum();
                }
                else
                {
                    v.coeff(k) > c.coeff(k) ? v.coeffRef(k) -= w.coeff(k) : v.coeffRef(k) += w.coeff(k);
                    ++w.coeffRef(k);
                    P -= log(D.coeff(k));
                }
            }
        }
    }
}

VectorXli DualPotENUM(const MatrixXld gso_coeff_mat, const VectorXld squared_norm_of_gso_vec, const VectorXld log_squared_norm_of_gso_vec, const int n)
{
    int i, r[n + 1];
    double tmp, R = log_squared_norm_of_gso_vec.coeff(0), P = 0;
    VectorXli w(n), v(n);
    w.setZero();
    v.setZero();
    Eigen::VectorXd c(n), D(n + 1);
    c.setZero();
    D.setZero();
    Eigen::MatrixXd sigma(n + 1, n);
    sigma.setZero();

    v.coeffRef(0) = 1;

    for (i = 0; i <= n; ++i)
    {
        r[i] = i;
    }
    for (int k = 0, last_nonzero = 0;;)
    {
        tmp = (double)v.coeff(k) - c.coeff(k);
        tmp *= tmp;
        D.coeffRef(k) = D.coeff(k + 1) + tmp * squared_norm_of_gso_vec.coeff(k);

        if (last_nonzero == 0)
        {
            P = log(D.coeff(0));
        }
        else if (last_nonzero == 1)
        {
            P = log(D.coeff(0));
            P += P;
        }
        else
        {
            P = 0;
            for (i = 0; i <= last_nonzero - 2; ++i)
            {
                P += log(D.coeff(i));
            }
            tmp = log(D.coeff(last_nonzero - 1));
            tmp += tmp;
            P += tmp;
        }

        R = log_squared_norm_of_gso_vec.head(last_nonzero + 1).array().sum();
        if (LOG099 + P > R)
        {
            if (k == 0)
                return v;
            else
            {
                P += log(D.coeff(k));
                --k;
                r[k] = (r[k] > r[k + 1] ? r[k] : r[k + 1]);
                for (i = r[k]; i > k; --i)
                    sigma.coeffRef(i, k) = sigma.coeff(i + 1, k) + gso_coeff_mat.coeff(i, k) * v(i);
                c.coeffRef(k) = -sigma.coeff(k + 1, k);
                v.coeffRef(k) = round(c.coeff(k));
                w.coeffRef(k) = 1;
            }
        }
        else
        {
            ++k;
            if (k == n)
            {
                v.setZero();
                return v;
            }
            else
            {
                r[k - 1] = k;
                if (k >= last_nonzero)
                {
                    last_nonzero = k;
                    ++v.coeffRef(k);
                    if (v.coeff(last_nonzero) >= 2)
                    {
                        ++k;
                        if (k == n)
                        {
                            v.setZero();
                            return v;
                        }
                        else
                        {
                            r[k - 1] = k;
                            last_nonzero = k;
                            v.coeffRef(last_nonzero) = 1;
                        }
                    }
                }
                else
                {
                    v.coeff(k) > c.coeff(k) ? v.coeffRef(k) -= w.coeff(k) : v.coeffRef(k) += w.coeff(k);
                    ++w.coeffRef(k);
                    P -= log(D.coeff(k));
                }
            }
        }
    }
}

/// @brief PotLLL
/// @param basis lattice basis
/// @param reduction_parameter Reduction parameter
/// @param n
/// @param m
void POT_LLL(MatrixXli &basis, const long double reduction_parameter, const int n, const int m)
{
    long double P, P_min, S;
    VectorXli t;
    VectorXld squared_norm_of_gso_vec(n);
    NTL::mat_ZZ c;
    c.SetDims(n, m);
    MatrixXld gso_coeff_mat(n, n);

    // LLL基底簡約
    for (int i = 0, j; i < n; ++i)
    {
        for (j = 0; j < m; ++j)
            c[i][j] = basis.coeff(i, j);
    }
    NTL::LLL(_, c, 99, 100);
    for (int i = 0, j; i < n; ++i)
    {
        for (j = 0; j < m; ++j)
            basis.coeffRef(i, j) = NTL::to_long(c[i][j]);
    }

    GSO(basis, squared_norm_of_gso_vec, gso_coeff_mat, n, m);

    for (int l = 0, l1, j, i, k, q; l < n;)
    {
        l1 = l - 1;
        // 部分サイズ基底簡約
        for (j = l1; j > -1; --j)
            if (gso_coeff_mat.coeff(l, j) > 0.5 || gso_coeff_mat.coeff(l, j) < -0.5)
            {
                q = round(gso_coeff_mat.coeff(l, j));
                basis.row(l) -= q * basis.row(j);
                gso_coeff_mat.row(l).head(j + 1) -= (long double)q * gso_coeff_mat.row(j).head(j + 1);
            }

        P = P_min = 1.0;
        k = 0;
        for (j = l1; j >= 0; --j)
        {
            S = (gso_coeff_mat.row(l).segment(j, l - j).array().square() * squared_norm_of_gso_vec.segment(j, l - j).array()).sum();
            P *= (squared_norm_of_gso_vec.coeff(l) + S) / squared_norm_of_gso_vec.coeff(j);
            if (P < P_min)
            {
                k = j;
                P_min = P;
            }
        }

        if (reduction_parameter > P_min)
        {
            // Deep insertion
            t = basis.row(l);
            for (j = l; j > k; --j)
                basis.row(j) = basis.row(j - 1);
            basis.row(k) = t;

            GSO(basis, squared_norm_of_gso_vec, gso_coeff_mat, n, m);
            l = k;
        }
        else
            ++l;
    }
}

/// @brief
/// @param basis
/// @param reduction_parameter
/// @param n
/// @param m
void DUAL_POT_LLL(MatrixXli &basis, const double reduction_parameter, const int n, const int m)
{
    double P, P_max, P_min, s;
    MatrixXld gso_coeff_mat(n, n), nu(n, n);
    gso_coeff_mat.setZero();
    nu.setZero();
    VectorXld squared_norm_of_gso_vec(n);
    squared_norm_of_gso_vec.setZero();
    NTL::mat_ZZ c;
    c.SetDims(n, m);

    // LLL基底簡約
    for (int i = 0, j; i < n; ++i)
    {
        for (j = 0; j < m; ++j)
            c[i][j] = basis.coeff(i, j);
    }
    NTL::LLL(_, c, 99, 100);
    for (int i = 0, j; i < n; ++i)
    {
        for (j = 0; j < m; ++j)
            basis.coeffRef(i, j) = NTL::to_long(c[i][j]);
    }

    GSO(basis, squared_norm_of_gso_vec, gso_coeff_mat, n, m);

    for (int k = n - 1, j, i, l, q; k >= 0;)
    {
        nu.coeffRef(k, k) = 1.0;

        // Dual size reduction
        for (j = k + 1; j < n; ++j)
        {
            nu.coeffRef(k, j) = 0;
            for (i = k; i < j; ++i)
                nu.coeffRef(k, j) -= gso_coeff_mat.coeff(j, i) * nu.coeff(k, i);

            if (nu.coeff(k, j) > 0.5 || nu.coeff(k, j) < -0.5)
            {
                q = round(nu.coeff(k, j));
                basis.row(j) += q * basis.row(k);
                nu.row(k).tail(n - j + 1) -= q * nu.row(j).tail(n - j + 1);
                gso_coeff_mat.row(j).head(k + 1) += q * gso_coeff_mat.row(k).head(k + 1);
            }
        }

        // Potential

        P = P_min = 1.0;
        l = n - 1;
        for (j = k + 1; j < n; ++j)
        {
            s = 0.0;
            for (i = k; i <= j; ++i)
                s += nu.coeff(k, i) * nu.coeff(k, i) / squared_norm_of_gso_vec.coeff(i);
            P *= squared_norm_of_gso_vec.coeff(j);
            P *= s;

            if (P < P_min)
            {
                l = j;
                P_min = P;
            }
        }

        if (reduction_parameter > P_min)
        {
            DualDeepInsertion(basis, m, k, l);
            GSO(basis, squared_norm_of_gso_vec, gso_coeff_mat, n, m);
            k = l;
        }
        else
            --k;
    }
}

/// @brief PotBKZ-reduces the lattice basis
/// @param basis lattice basis
/// @param beta Block size
/// @param reduction_parameter Reduction parameter
/// @param n Rank of lattice
/// @param m
/// @param fp Files of potentials of lattice
void POT_BKZ(MatrixXli &basis, const int beta, const double reduction_parameter, const int n, const int m)
{
    const int n1 = n - 1, n2 = n - 2;
    VectorXli v, w;
    MatrixXli tmp_b(n, n);
    VectorXld squared_norm_of_gso_vec(n), log_squared_norm_of_gso_vec(n);
    MatrixXld gso_coeff_mat(n, n);
    NTL::mat_ZZ cc;

    POT_LLL(basis, reduction_parameter, n, m);
    GSO(basis, squared_norm_of_gso_vec, log_squared_norm_of_gso_vec, gso_coeff_mat, n, m);

    for (int z = 0, j = 0, i, k, l, kj1; z < n - 1;)
    {
        //        printf("z = %reduction_parameter\n", z);
        ++Tr;

        if (j == n2)
        {
            j = 0;
            ++PotTour;
            ++Tr;
        }
        ++j;
        k = (j + beta - 1 < n1 ? j + beta - 1 : n1);
        kj1 = k - j + 1;
        v.resize(kj1);
        v.setZero();

        /* enumerate a shortest vector*/
        v = PotENUM(gso_coeff_mat.block(j, j, kj1, kj1), squared_norm_of_gso_vec.segment(j, kj1), log_squared_norm_of_gso_vec.segment(j, kj1), kj1);

        if (!v.isZero())
        {
            z = 0;

            w = v * basis.block(j, 0, kj1, m);
            cc.SetDims(n + 1, m);
            for (l = 0; l < m; ++l)
            {
                for (i = 0; i < j; ++i)
                    cc[i][l] = basis.coeffRef(i, l);
                cc[j][l] = w[l];
                for (i = j + 1; i < n + 1; ++i)
                    cc[i][l] = basis.coeffRef(i - 1, l);
            }
            NTL::LLL(_, cc, 99, 100);

            for (i = 0; i < n; ++i)
                for (l = 0; l < m; ++l)
                    basis.coeffRef(i, l) = NTL::to_long(cc[i + 1][l]);

            POT_LLL(basis, reduction_parameter, n, m);
            GSO(basis, squared_norm_of_gso_vec, log_squared_norm_of_gso_vec, gso_coeff_mat, n, m);
        }
        else
            ++z;
    }
}

/// @brief Dual version of PotBKZ
/// @param basis
/// @param beta
/// @param delta
/// @param n
/// @param m
/// @param fp
void DUAL_POT_BKZ(MatrixXli &basis, const int beta, const double delta, const int n, const int m)
{
    VectorXli x;
    MatrixXli tmp_b;
    VectorXld squared_norm_of_gso_vec(n), log_squared_norm_of_gso_vec(n), dual_squared_norm_of_gso_vec, logC;
    squared_norm_of_gso_vec.setZero();
    log_squared_norm_of_gso_vec.setZero();
    MatrixXld gso_coeff_mat(n, n), hmu, BB;
    gso_coeff_mat.setZero();

    DUAL_POT_LLL(basis, 0.99, n, m);

    GSO(basis, squared_norm_of_gso_vec, log_squared_norm_of_gso_vec, gso_coeff_mat, n, m);

    for (int z = n, i, j = n, k, dim_of_block_lattice; z > 2;)
    {
        if (j == 1)
        {
            ++PotTour;
            j = n;
        }
        --j;
        k = (j - beta + 1 > 0 ? j - beta + 1 : 0);
        dim_of_block_lattice = j - k + 1;

        //        printf("z = %dim_of_block_lattice\n", z);
        ++Tr;

        dual_squared_norm_of_gso_vec.resize(dim_of_block_lattice);
        logC.resize(dim_of_block_lattice);
        hmu.resize(dim_of_block_lattice, dim_of_block_lattice);
        DualGSO(squared_norm_of_gso_vec.segment(k, dim_of_block_lattice), log_squared_norm_of_gso_vec.segment(k, dim_of_block_lattice), gso_coeff_mat.block(k, k, dim_of_block_lattice, dim_of_block_lattice), dual_squared_norm_of_gso_vec, logC, hmu, dim_of_block_lattice, dim_of_block_lattice);

        // Dual Enumeration
        x = DualPotENUM(hmu, dual_squared_norm_of_gso_vec, logC, dim_of_block_lattice);

        if (x.isZero())
        {
            --z;
        }
        else
        {
            z = n;

            tmp_b = Insert(basis.block(k, 0, dim_of_block_lattice, m), x, dim_of_block_lattice, m);
            basis.block(k, 0, dim_of_block_lattice, m) = tmp_b.block(0, 0, dim_of_block_lattice, m);

            POT_LLL(basis, delta, n, m);
            GSO(basis, squared_norm_of_gso_vec, log_squared_norm_of_gso_vec, gso_coeff_mat, n, m);
        }
    }
}

/**
 * @brief
 *
 * @param basis
 * @param beta
 * @param reduction_parameter
 * @param n
 * @param m
 * @param fp
 */
void SELF_DUAL_POT_BKZ(MatrixXli &basis, const int beta, const double reduction_parameter, const int n, const int m)
{
    const int n1 = n - 1, n2 = n - 2;
    VectorXli v, w;
    MatrixXli tmp_b(n, n);
    VectorXld squared_norm_of_gso_vec(n), log_squared_norm_of_gso_vec(n);
    MatrixXld gso_coeff_mat(n, n);
    NTL::mat_ZZ cc;
    VectorXld dual_squared_norm_of_gso_vec, logC;
    MatrixXld hmu, BB;

    GSO(basis, squared_norm_of_gso_vec, log_squared_norm_of_gso_vec, gso_coeff_mat, n, m);

    DUAL_POT_LLL(basis, 0.99, n, m);

    for (int zp = 0, jp = 0, i, k, l, kj1, zd = n, jd = n, IsPrimal = 1; zp < n - 1 && zd > 1;)
    {
        /// ================================
        /// Primal part
        /// ================================
        //        printf("Primal = %3d, Dual = %3d\n", zp, zd);

        if (IsPrimal)
        {
            if (jp == n2)
            {
                jp = 0;
                ++PotTour;
                IsPrimal = 0;
            }
            ++jp;
            k = (jp + beta - 1 < n1 ? jp + beta - 1 : n1);
            kj1 = k - jp + 1;
            v.resize(kj1);
            v.setZero();

            /* enumerate a shortest vector*/
            v = PotENUM(gso_coeff_mat.block(jp, jp, kj1, kj1), squared_norm_of_gso_vec.segment(jp, kj1), log_squared_norm_of_gso_vec.segment(jp, kj1), kj1);

            if (!v.isZero())
            {
                zp = 0;

                w = v * basis.block(jp, 0, kj1, m);
                cc.SetDims(n + 1, m);
                for (l = 0; l < m; ++l)
                {
                    for (i = 0; i < jp; ++i)
                        cc[i][l] = basis.coeffRef(i, l);
                    cc[jp][l] = w[l];
                    for (i = jp + 1; i < n + 1; ++i)
                        cc[i][l] = basis.coeffRef(i - 1, l);
                }
                NTL::LLL(_, cc, 99, 100);

                for (i = 0; i < n; ++i)
                    for (l = 0; l < m; ++l)
                        basis.coeffRef(i, l) = NTL::to_long(cc[i + 1][l]);

                DUAL_POT_LLL(basis, reduction_parameter, n, m);
                GSO(basis, squared_norm_of_gso_vec, log_squared_norm_of_gso_vec, gso_coeff_mat, n, m);
            }
            else
                ++zp;
        }

        /// ================================
        /// Dual part
        /// ================================
        if ((!IsPrimal))
        {
            if (jd == 1)
            {
                ++PotTour;
                jd = n;
                IsPrimal = 1;
            }
            --jd;
            k = (jd - beta + 1 > 0 ? jd - beta + 1 : 0);
            kj1 = jd - k + 1;

            dual_squared_norm_of_gso_vec.resize(kj1);
            logC.resize(kj1);
            hmu.resize(kj1, kj1);
            DualGSO(squared_norm_of_gso_vec.segment(k, kj1), log_squared_norm_of_gso_vec.segment(k, kj1), gso_coeff_mat.block(k, k, kj1, kj1), dual_squared_norm_of_gso_vec, logC, hmu, kj1, kj1);

            // Dual Enumeration
            v = DualPotENUM(hmu, dual_squared_norm_of_gso_vec, logC, kj1);

            if (v.isZero())
            {
                --zd;
            }
            else
            {
                zd = n;

                tmp_b = Insert(basis.block(k, 0, kj1, m), v, kj1, m);
                basis.block(k, 0, kj1, m) = tmp_b.block(0, 0, kj1, m);

                POT_LLL(basis, reduction_parameter, n, m);
                GSO(basis, squared_norm_of_gso_vec, log_squared_norm_of_gso_vec, gso_coeff_mat, n, m);
            }
        }
    }
}

extern "C" long **PotLLL(long **basis, const double reduction_parameter, const int n, const int m)
{
    int i, j;
    MatrixXli B(n, m);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
        {
            B.coeffRef(i, j) = basis[i][j];
        }

    POT_LLL(B, reduction_parameter, n, m);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            basis[i][j] = B.coeff(i, j);
    return basis;
}

extern "C" long **DualPotLLL(long **basis, const double reduction_parameter, const int n, const int m)
{
    int i, j;
    MatrixXli B(n, m);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            B.coeffRef(i, j) = basis[i][j];
    DUAL_POT_LLL(B, reduction_parameter, n, m);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            basis[i][j] = B.coeff(i, j);
    return basis;
}

extern "C" long **PotBKZ(long **basis, const int beta, const double reduction_parameter, const int n, const int m)
{
    int i, j;
    MatrixXli B(n, m);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            B.coeffRef(i, j) = basis[i][j];
    POT_BKZ(B, beta, reduction_parameter, n, m);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            basis[i][j] = B.coeff(i, j);
    return basis;
}

extern "C" long **DualPotBKZ(long **basis, const int beta, const double delta, const int n, const int m)
{
    int i, j;
    MatrixXli B(n, m);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            B.coeffRef(i, j) = basis[i][j];
    DUAL_POT_BKZ(B, beta, delta, n, m);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            basis[i][j] = B.coeff(i, j);
    return basis;
}

extern "C" long **SelfDualPotBKZ(long **basis, const int beta, const double reduction_parameter, const int n, const int m)
{
    int i, j;
    MatrixXli B(n, m);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            B.coeffRef(i, j) = basis[i][j];
    SELF_DUAL_POT_BKZ(B, beta, reduction_parameter, n, m);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            basis[i][j] = B.coeff(i, j);
    return basis;
}
