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

static int PotTour = 0, BKZTour = 0, SelfPotTour = 0, Tr = 0;
static NTL::ZZ _;
static FILE *_TOUR_DTA_ = fopen("tour_data.csv", "wt");

#define LOG099 -0.010050335853501441183548857558547706085515007674629873378

/// @brief Dual version of deep-insertion
/// @param b Lattice basis
/// @param n Rank of lattice
/// @param k
/// @param l
void DualDeepInsertion(MatrixXli &b, const int n, const int k, const int l)
{
    const VectorXli t = b.row(k);
    for (int j = k; j < l; ++j)
        b.row(j) = b.row(j + 1);
    b.row(l) = t;
}

/// @brief P.239
/// @param NewBasis
/// @param Basis
/// @param x
MatrixXli Insert(const MatrixXli b, const VectorXli x, const int n, const int m)
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
    return U * b;
}

MatrixXld dual(const MatrixXli b, const int n, const int m)
{
    return b.cast<long double>().inverse().transpose();
}

void Randomization(MatrixXli &b, const int n, const int m)
{
    int i;
    MatrixXli c = b;
    std::vector<int> index(n);
    std::random_device seed_gen;
    std::mt19937 get_rand_mt(seed_gen());

    for (i = 0; i < n; ++i)
        index[i] = i;

    std::shuffle(index.begin(), index.end(), get_rand_mt);

    for (i = 0; i < n; ++i)
        b.row(i) = c.row(index[i]);
}

void GSO(const MatrixXli b, VectorXld &B, MatrixXld &mu, const int n, const int m)
{
    MatrixXld GSOb(n, m);

    for (int i = 0, j; i < n; ++i)
    {
        mu.coeffRef(i, i) = 1.0;
        GSOb.row(i) = b.row(i).cast<long double>();
        for (j = 0; j < i; ++j)
        {
            mu.coeffRef(i, j) = b.row(i).cast<long double>().dot(GSOb.row(j)) / GSOb.row(j).dot(GSOb.row(j));
            GSOb.row(i) -= mu.coeff(i, j) * GSOb.row(j);
        }
        B.coeffRef(i) = GSOb.row(i).dot(GSOb.row(i));
    }
}

void GSO(const MatrixXld b, VectorXld &B, MatrixXld &mu, const int n, const int m)
{
    MatrixXld GSOb(n, m);

    for (int i = 0, j; i < n; ++i)
    {
        mu.coeffRef(i, i) = 1.0;
        GSOb.row(i) = b.row(i);
        for (j = 0; j < i; ++j)
        {
            mu.coeffRef(i, j) = b.row(i).dot(GSOb.row(j)) / GSOb.row(j).dot(GSOb.row(j));
            GSOb.row(i) -= mu.coeff(i, j) * GSOb.row(j);
        }
        B.coeffRef(i) = GSOb.row(i).dot(GSOb.row(i));
    }
}

void GSO(const MatrixXli b, VectorXld &B, VectorXld &logB, MatrixXld &mu, const int n, const int m)
{
    MatrixXld GSOb(n, m);

    for (int i = 0, j; i < n; ++i)
    {
        mu.coeffRef(i, i) = 1.0;
        GSOb.row(i) = b.row(i).cast<long double>();
        for (j = 0; j < i; ++j)
        {
            mu.coeffRef(i, j) = b.row(i).cast<long double>().dot(GSOb.row(j)) / GSOb.row(j).dot(GSOb.row(j));
            GSOb.row(i) -= mu.coeff(i, j) * GSOb.row(j);
        }
        B.coeffRef(i) = GSOb.row(i).dot(GSOb.row(i));
        logB.coeffRef(i) = log(B.coeff(i));
    }
}

/// @brief Computes GSO-informations of dual-lattice basis using GSO-information of lattice
/// @param B Squared norms of GSO-vectors
/// @param mu GSO-coefficient mattix
/// @param C Squared norms of dual-GSO-vectors
/// @param hmu dual-GSO-coefficient mattix
/// @param n
/// @param m
void DualGSO(const VectorXld B, const VectorXld logB, const MatrixXld mu, VectorXld &C, VectorXld &logC, MatrixXld &hmu, const int n, const int m)
{
    for (int i = 0, j; i < n; ++i)
    {
        C.coeffRef(i) = 1.0 / B.coeff(i);
        logC.coeffRef(i) = -logB.coeff(i);
        for (j = i + 1; j < m; ++j)
            hmu.coeffRef(i, j) = -mu.row(j).segment(i, j - i).dot(hmu.row(i).segment(i, j - i));
    }
}

/// @brief Conputes log value of potential of basis
/// @param B Squared norms of GSO-vectors
/// @param n Rank of lattice
/// @return double
long double logPot(const VectorXld B, const int n)
{
    NTL::RR logp;
    logp = 0;
    for (int i = 0; i < n; ++i)
        logp += (n - i) * NTL::log(NTL::to_RR((double)B.coeff(i)));
    return NTL::to_double(logp);
}

long double rho(const VectorXld B, const int n, const int m)
{
    long double S = 0, T = 0;
    for (int i = 0; i < n; ++i)
    {
        S += (i + 1) * log(B.coeff(i));
        T += log(B.coeff(i));
    }
    return 12 * (S - 0.5 * (n + 1) * T) / (n * (n * n - 1));
}

/// @brief Enumerates a vector whose norm is shorter than R on the lattice
/// @param mu GSO-coefficient matrix
/// @param B Squared-norms of row vectors of GSO-matrix
/// @param n Rank of lattice
/// @param R An upper bound of enumeration
/// @return
VectorXli ENUM(const MatrixXld mu, const VectorXld B, VectorXld &rho, const int n, const double R)
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
        r[i] = i;

    for (int k = 0, last_nonzero = 0;;)
    {
        tmp = (double)v.coeff(k) - c.coeff(k);
        tmp *= tmp;
        rho.coeffRef(k) = rho.coeff(k + 1) + tmp * B.coeff(k); // rho[k]=∥πₖ(v)∥
        if (rho.coeff(k) <= R)
        {
            if (k == 0)
                return v;
            else
            {
                --k;
                r[k] = (r[k] > r[k + 1] ? r[k] : r[k + 1]);
                for (i = r[k]; i > k; --i)
                    sigma.coeffRef(i, k) = sigma.coeff(i + 1, k) + mu.coeff(i, k) * v.coeff(i);
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
/// @param mu GSO-coefficient matrix
/// @param B Squared-norms of row vectors of GSO-matrix
/// @param n Rank of lattice
/// @return VectorXli the shortest vector
VectorXli enumerate(const MatrixXld mu, const VectorXld B, VectorXld &rho, const int n)
{
    VectorXli enum_v(n), pre_enum_v(n);
    enum_v.setZero();
    pre_enum_v.setZero();
    VectorXld pre_rho(n + 1);
    pre_rho.setZero();
    for (double R = B.coeff(0);;)
    {
        pre_rho = rho;
        pre_enum_v = enum_v;
        enum_v = ENUM(mu, B, rho, n, R);
        if (enum_v.isZero())
        {
            rho = pre_rho;
            return pre_enum_v;
        }
        R *= 0.99;
    }
}

/// @brief Enumerates a short and small potential lattice vector
/// @param mu GSO-coefficient matrix.
/// @param B Squared norms of GSO vectors.
/// @param logB Logarithm valude of squared norms of GSO-vectors.
/// @param n Rank of lattice.
/// @return VectorXli
VectorXli PotENUM(const MatrixXld mu, const VectorXld B, const VectorXld logB, const int n)
{
    int i, r[n + 1];
    double tmp, R = logB.coeff(0), P = 0;
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
        D.coeffRef(k) = D.coeff(k + 1) + tmp * B.coeff(k);

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
                    sigma.coeffRef(i, k) = sigma.coeff(i + 1, k) + mu.coeff(i, k) * v(i);
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
#if 0
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
#endif
                    P = 0;
                    R = logB.head(last_nonzero + 1).array().sum();
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

VectorXli DualPotENUM(const MatrixXld mu, const VectorXld B, const VectorXld logB, const int n)
{
    int i, r[n + 1];
    double tmp, R = logB.coeff(0), P = 0;
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
        D.coeffRef(k) = D.coeff(k + 1) + tmp * B.coeff(k);

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

        R = logB.head(last_nonzero + 1).array().sum();
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
                    sigma.coeffRef(i, k) = sigma.coeff(i + 1, k) + mu.coeff(i, k) * v(i);
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
#if 0
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
#endif
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
/// @param b lattice basis
/// @param d Reduction parameter
/// @param n
/// @param m
void PotLLL_(MatrixXli &b, const long double d, const int n, const int m)
{
    long double P, P_min, S;
    VectorXli t;
    VectorXld B(n);
    NTL::mat_ZZ c;
    c.SetDims(n, m);
    MatrixXld mu(n, n);

    // LLL基底簡約
    for (int i = 0, j; i < n; ++i)
    {
        for (j = 0; j < m; ++j)
            c[i][j] = b.coeff(i, j);
    }
    NTL::LLL(_, c, 99, 100);
    for (int i = 0, j; i < n; ++i)
    {
        for (j = 0; j < m; ++j)
            b.coeffRef(i, j) = NTL::to_long(c[i][j]);
    }

    GSO(b, B, mu, n, m);

    for (int l = 0, l1, j, i, k, q; l < n;)
    {
        l1 = l - 1;
        // 部分サイズ基底簡約
        for (j = l1; j > -1; --j)
            if (mu.coeff(l, j) > 0.5 || mu.coeff(l, j) < -0.5)
            {
                q = round(mu.coeff(l, j));
                b.row(l) -= q * b.row(j);
                mu.row(l).head(j + 1) -= (long double)q * mu.row(j).head(j + 1);
            }

        P = P_min = 1.0;
        k = 0;
        for (j = l1; j >= 0; --j)
        {
            S = (mu.row(l).segment(j, l - j).array().square() * B.segment(j, l - j).array()).sum();
            P *= (B.coeff(l) + S) / B.coeff(j);
            if (P < P_min)
            {
                k = j;
                P_min = P;
            }
        }

        if (d > P_min)
        {
            // Deep insertion
            t = b.row(l);
            for (j = l; j > k; --j)
                b.row(j) = b.row(j - 1);
            b.row(k) = t;

            GSO(b, B, mu, n, m);
            l = k;
        }
        else
            ++l;
    }
}

/// @brief
/// @param b
/// @param d
/// @param n
/// @param m
void DualPotLLL_(MatrixXli &b, const double d, const int n, const int m)
{
    double P, P_max, P_min, s;
    MatrixXld mu(n, n), nu(n, n);
    mu.setZero();
    nu.setZero();
    VectorXld B(n);
    B.setZero();
    NTL::mat_ZZ c;
    c.SetDims(n, m);

    // LLL基底簡約
    for (int i = 0, j; i < n; ++i)
    {
        for (j = 0; j < m; ++j)
            c[i][j] = b.coeff(i, j);
    }
    NTL::LLL(_, c, 99, 100);
    for (int i = 0, j; i < n; ++i)
    {
        for (j = 0; j < m; ++j)
            b.coeffRef(i, j) = NTL::to_long(c[i][j]);
    }

    GSO(b, B, mu, n, m);

    for (int k = n - 1, j, i, l, q; k >= 0;)
    {
        nu.coeffRef(k, k) = 1.0;

        // Dual size reduction
        for (j = k + 1; j < n; ++j)
        {
            nu.coeffRef(k, j) = 0;
            for (i = k; i < j; ++i)
                nu.coeffRef(k, j) -= mu.coeff(j, i) * nu.coeff(k, i);

            if (nu.coeff(k, j) > 0.5 || nu.coeff(k, j) < -0.5)
            {
                q = round(nu.coeff(k, j));
                b.row(j) += q * b.row(k);
                nu.row(k).tail(n - j + 1) -= q * nu.row(j).tail(n - j + 1);
                mu.row(j).head(k + 1) += q * mu.row(k).head(k + 1);
            }
        }

        // Potential

        P = P_min = 1.0;
        l = n - 1;
        for (j = k + 1; j < n; ++j)
        {
            s = 0.0;
            for (i = k; i <= j; ++i)
                s += nu.coeff(k, i) * nu.coeff(k, i) / B.coeff(i);
            P *= B.coeff(j);
            P *= s;

            if (P < P_min)
            {
                l = j;
                P_min = P;
            }
        }

        if (d > P_min)
        {
            DualDeepInsertion(b, m, k, l);
            GSO(b, B, mu, n, m);
            k = l;
        }
        else
            --k;
    }
}

/// @brief BKZ-reduces the lattice basis
/// @param b Lattice basis matrix
/// @param beta Block size
/// @param d Reduction parameter
/// @param lp Tour limit
/// @param n Rank of lattice
/// @param m
void BKZ_(MatrixXli &b, const int beta, const double d, const int lp, const int n, const int m, FILE *fp)
{
    const int n1 = n - 1;
    VectorXli v, w;
    NTL::mat_ZZ c;
    VectorXld B(n), s;
    B.setZero();
    MatrixXld mu(n, n);
    mu.setIdentity();

    GSO(b, B, mu, n, m);

    for (int z = 0, j, t = 0, i, k = 0, h, lk1, l; z < n - 1;)
    {
        if (BKZTour >= lp)
            break;
        fprintf(fp, "%Lf\n", logPot(B, n));
        printf("z = %d\n", z);

        if (k == n1)
        {
            k = 0;
            ++BKZTour;
        }
        ++k;
        l = std::min(k + beta - 1, n);
        h = std::min(l + 1, n);
        lk1 = l - k + 1;

        s.resize(lk1 + 1);
        s.setZero();
        w = enumerate(mu.block(k - 1, k - 1, lk1, lk1), B.segment(k - 1, lk1), s, lk1);
        if (B.coeff(k - 1) > s.coeff(k - 1) && (!w.isZero()))
        {
            z = 0;

            v = w * b.block(k - 1, 0, lk1, m);

            // Inserts and removes linearly-dependentness by LLL-reducing
            c.SetDims(h + 1, m);
            for (j = 0; j < m; ++j)
            {
                for (i = 0; i < k - 1; ++i)
                    c[i][j] = b.coeff(i, j);
                c[k - 1][j] = v.coeff(j);
                for (i = k; i < h + 1; ++i)
                    c[i][j] = b.coeff(i - 1, j);
            }
            NTL::LLL(_, c, 99, 100);
            for (i = 1; i <= h; ++i)
                for (j = 0; j < m; ++j)
                    b.coeffRef(i - 1, j) = NTL::to_long(c[i][j]);
            GSO(b, B, mu, n, m);
        }
        else
        {
            ++z;

            c.SetDims(h, m);
            for (j = 0; j < m; ++j)
                for (i = 0; i < h; ++i)
                    c[i][j] = b.coeff(i, j);
            NTL::LLL(_, c, 99, 100);
            for (i = 0; i < h; ++i)
                for (j = 0; j < m; ++j)
                    b.coeffRef(i, j) = NTL::to_long(c[i][j]);
            GSO(b, B, mu, n, m);
        }
    }
}

/// @brief PotBKZ-reduces the lattice basis
/// @param b lattice basis
/// @param beta Block size
/// @param d Reduction parameter
/// @param n Rank of lattice
/// @param m
/// @param fp Files of potentials of lattice
void PotBKZ_(MatrixXli &b, const int beta, const double d, const int n, const int m, FILE *fp)
{
    const int n1 = n - 1, n2 = n - 2;
    VectorXli v, w;
    MatrixXli tmp_b(n, n);
    VectorXld B(n), logB(n);
    MatrixXld mu(n, n);
    NTL::mat_ZZ cc;

    GSO(b, B, logB, mu, n, m);
    fprintf(fp, "%Lf\n", logPot(B, n));
    for (int z = 0, j = 0, i, k, l, kj1; z < n - 1;)
    {
        printf("%d: z = %d\n", PotTour, z);
        ++Tr;
        fprintf(fp, "%Lf\n", logPot(B, n));

        if (j == n2)
        {
            j = 0;
            ++PotTour;
            ++Tr;
        }
        ++j;
        k = (j + beta - 1 < n1 ? j + beta - 1 : n1);
        // k = std::min(j + beta - 1, n1);
        kj1 = k - j + 1;
        v.resize(kj1);
        v.setZero();

        /* enumerate a shortest vector*/
        v = PotENUM(mu.block(j, j, kj1, kj1), B.segment(j, kj1), logB.segment(j, kj1), kj1);

        if (!v.isZero())
        {
            z = 0;

            w = v * b.block(j, 0, kj1, m);
            cc.SetDims(n + 1, m);
            for (l = 0; l < m; ++l)
            {
                for (i = 0; i < j; ++i)
                    cc[i][l] = b.coeffRef(i, l);
                cc[j][l] = w[l];
                for (i = j + 1; i < n + 1; ++i)
                    cc[i][l] = b.coeffRef(i - 1, l);
            }
            NTL::LLL(_, cc, 99, 100);

            for (i = 0; i < n; ++i)
                for (l = 0; l < m; ++l)
                    b.coeffRef(i, l) = NTL::to_long(cc[i + 1][l]);

            PotLLL_(b, d, n, m);
            GSO(b, B, logB, mu, n, m);
        }
        else
            ++z;
    }
    fprintf(_TOUR_DTA_, "%d\n", Tr);
    // printf("Primal: %d\n", Tr);
}


/// @brief Dual version of PotBKZ
/// @param b
/// @param beta
/// @param delta
/// @param n
/// @param m
/// @param fp
void DualPotBKZ_(MatrixXli &b, const int beta, const double delta, const int n, const int m, FILE *fp)
{
    VectorXli x;
    MatrixXli tmp_b;
    VectorXld B(n), logB(n), C, logC;
    B.setZero();
    logB.setZero();
    MatrixXld mu(n, n), hmu, BB;
    mu.setZero();

    GSO(b, B, logB, mu, n, m);
    fprintf(fp, "%Lf\n", logPot(B, n));

    // DualPotLLL(b, 0.99, n, m);
    for (int z = n, i, j = n, k, d; z > 1;)
    {
        if (j == 1)
        {
            ++PotTour;
            j = n;
        }
        --j;
        k = (j - beta + 1 > 0 ? j - beta + 1 : 0);
        // k = std::max(j - beta + 1, 0);
        d = j - k + 1;

        printf("z = %d\n", z);
        ++Tr;
        fprintf(fp, "%Lf\n", logPot(B, n));

        C.resize(d);
        logC.resize(d);
        hmu.resize(d, d);
        DualGSO(B.segment(k, d), logB.segment(k, d), mu.block(k, k, d, d), C, logC, hmu, d, d);

        // Dual Enumeration
        x = DualPotENUM(hmu, C, logC, d);

        if (x.isZero())
        {
            --z;
        }
        else
        {
            z = n;

            tmp_b = Insert(b.block(k, 0, d, m), x, d, m);
            b.block(k, 0, d, m) = tmp_b.block(0, 0, d, m);

            // DualPotLLL(b, delta, n, m);
            PotLLL_(b, delta, n, m);
            GSO(b, B, logB, mu, n, m);
        }
    }
    fprintf(_TOUR_DTA_, "%d\n", Tr);
    printf("Dual: %d\n", Tr);
}


/**
 * @brief
 *
 * @param b
 * @param beta
 * @param d
 * @param n
 * @param m
 * @param fp
 */
void SelfDualPotBKZ_each_(MatrixXli &b, const int beta, const double d, const int n, const int m, FILE *fp)
{
    const int n1 = n - 1, n2 = n - 2;
    VectorXli v, w;
    MatrixXli tmp_b(n, n);
    VectorXld B(n), logB(n);
    MatrixXld mu(n, n);
    NTL::mat_ZZ cc;
    VectorXld C, logC;
    MatrixXld hmu, BB;

    GSO(b, B, logB, mu, n, m);
    fprintf(fp, "%Lf\n", logPot(B, n));

    DualPotLLL_(b, 0.99, n, m);

    for (int zp = 0, jp = 0, i, k, l, kj1, zd = n, jd = n, IsPrimal = 1; zp < n - 1 && zd > 1;)
    {
        /// ================================
        /// Primal part
        /// ================================
        printf("Primal = %3d, Dual = %3d\n", zp, zd);

        if (IsPrimal)
        {
            if (jp == n2)
            {
                jp = 0;
                ++SelfPotTour;
                IsPrimal = 0;
            }
            ++jp;
            k = (jp + beta - 1 < n1 ? jp + beta - 1 : n1);
            kj1 = k - jp + 1;

            fprintf(fp, "%Lf\n", logPot(B, n));

            v.resize(kj1);
            v.setZero();

            /* enumerate a shortest vector*/
            v = PotENUM(mu.block(jp, jp, kj1, kj1), B.segment(jp, kj1), logB.segment(jp, kj1), kj1);

            if (!v.isZero())
            {
                zp = 0;

                w = v * b.block(jp, 0, kj1, m);
                cc.SetDims(n + 1, m);
                for (l = 0; l < m; ++l)
                {
                    for (i = 0; i < jp; ++i)
                        cc[i][l] = b.coeffRef(i, l);
                    cc[jp][l] = w[l];
                    for (i = jp + 1; i < n + 1; ++i)
                        cc[i][l] = b.coeffRef(i - 1, l);
                }
                NTL::LLL(_, cc, 99, 100);

                for (i = 0; i < n; ++i)
                    for (l = 0; l < m; ++l)
                        b.coeffRef(i, l) = NTL::to_long(cc[i + 1][l]);

                DualPotLLL_(b, d, n, m);
                GSO(b, B, logB, mu, n, m);
            }
            else
                ++zp;
        }

        /// ================================
        /// Dual part
        /// ================================
        if (!IsPrimal)
        {
            if (jd == 1)
            {
                ++SelfPotTour;
                jd = n;
                IsPrimal = 1;
            }
            --jd;
            k = (jd - beta + 1 > 0 ? jd - beta + 1 : 0);
            kj1 = jd - k + 1;

            fprintf(fp, "%Lf\n", logPot(B, n));

            C.resize(kj1);
            logC.resize(kj1);
            hmu.resize(kj1, kj1);
            DualGSO(B.segment(k, kj1), logB.segment(k, kj1), mu.block(k, k, kj1, kj1), C, logC, hmu, kj1, kj1);

            // Dual Enumeration
            v = DualPotENUM(hmu, C, logC, kj1);

            if (v.isZero())
            {
                --zd;
            }
            else
            {
                zd = n;

                tmp_b = Insert(b.block(k, 0, kj1, m), v, kj1, m);
                b.block(k, 0, kj1, m) = tmp_b.block(0, 0, kj1, m);

                PotLLL_(b, d, n, m);
                GSO(b, B, logB, mu, n, m);
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

    PotLLL_(B, reduction_parameter, n, m);
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
    DualPotLLL_(B, reduction_parameter, n, m);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            basis[i][j] = B.coeff(i, j);
    return basis;
}

extern "C" long **PotBKZ(long **basis, const int beta, const double reduction_parameter, const int n, const int m)
{
    FILE *_Pot_FILE_ = fopen("potential_graph.csv", "wt");
    fprintf(_Pot_FILE_, "Potential\n");
    int i, j;
    MatrixXli B(n, m);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            B.coeffRef(i, j) = basis[i][j];
    PotBKZ_(B, beta, reduction_parameter, n, m, _Pot_FILE_);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            basis[i][j] = B.coeff(i, j);
    fclose(_Pot_FILE_);
    return basis;
}

extern "C" long **DualPotBKZ(long **basis, const int beta, const double delta, const int n, const int m)
{
    FILE *_Dual_Pot_FILE_ = fopen("dual_potential_graph.csv", "wt");
    int i, j;
    MatrixXli B(n, m);
    fprintf(_Dual_Pot_FILE_, "Potential\n");

    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            B.coeffRef(i, j) = basis[i][j];
    DualPotBKZ_(B, beta, delta, n, m, _Dual_Pot_FILE_);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            basis[i][j] = B.coeff(i, j);
    fclose(_Dual_Pot_FILE_);
    return basis;
}

extern "C" long **SelfDualPotBKZ(long **basis, const int beta, const double reduction_parameter, const int n, const int m)
{
    FILE *_Dual_Primal_ = fopen("potential_graph_DualENUM_and_PrimalLLL.csv", "wt");
    fprintf(_Dual_Primal_, "Potential\n");
    int i, j;
    MatrixXli B(n, m);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            B.coeffRef(i, j) = basis[i][j];
    SelfDualPotBKZ_each_(B, beta, reduction_parameter, n, m, _Dual_Primal_);
    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            basis[i][j] = B.coeff(i, j);
    fclose(_Dual_Primal_);
    return basis;
}
