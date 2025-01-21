#ifndef SDPB_HPP
#define SDPB_HPP

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <NTL/RR.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>

typedef Eigen::Matrix<long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXli;
typedef Eigen::Matrix<long, 1, Eigen::Dynamic> VectorXli;
typedef Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXld;
typedef Eigen::Matrix<long double, 1, Eigen::Dynamic> VectorXld;

#define LOG099 -0.010050335853501441183548857558547706085515007674629873378

class Lattice
{
private:
    int PotTour = 0, BKZTour = 0, SelfPotTour = 0, Tr = 0;
    NTL::ZZ _;
    FILE *_TOUR_DTA_ = fopen("tour_data.csv", "wt");

public:
    MatrixXli basis;

    void setDims(const int n, const int m)
    {
        basis.resize(n, m);
    }

    /// @brief Dual version of deep-insertion
    /// @param basis Lattice basis
    /// @param n Rank of lattice
    /// @param k
    /// @param l
    void DualDeepInsertion(const int n, const int k, const int l)
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

    void GSO(VectorXld &B, MatrixXld &mu, const int n, const int m)
    {
        MatrixXld GSOb(n, m);

        for (int i = 0, j; i < n; ++i)
        {
            mu.coeffRef(i, i) = 1.0;
            GSOb.row(i) = basis.row(i).cast<long double>();
            for (j = 0; j < i; ++j)
            {
                mu.coeffRef(i, j) = basis.row(i).cast<long double>().dot(GSOb.row(j)) / GSOb.row(j).dot(GSOb.row(j));
                GSOb.row(i) -= mu.coeff(i, j) * GSOb.row(j);
            }
            B.coeffRef(i) = GSOb.row(i).dot(GSOb.row(i));
        }
    }

    void GSO(VectorXld &B, VectorXld &logB, MatrixXld &mu, const int n, const int m)
    {
        MatrixXld GSOb(n, m);

        for (int i = 0, j; i < n; ++i)
        {
            mu.coeffRef(i, i) = 1.0;
            GSOb.row(i) = basis.row(i).cast<long double>();
            for (j = 0; j < i; ++j)
            {
                mu.coeffRef(i, j) = basis.row(i).cast<long double>().dot(GSOb.row(j)) / GSOb.row(j).dot(GSOb.row(j));
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

    VectorXli PotENUM(const MatrixXld mu, const VectorXld B, const VectorXld logB, const int n);

    VectorXli DualPotENUM(const MatrixXld mu, const VectorXld B, const VectorXld logB, const int n);
    
    void PotLLL_(const long double d, const int n, const int m);

    void DualPotLLL_(const double d, const int n, const int m);

    void BKZ_(const int beta, const double d, const int lp, const int n, const int m, FILE *fp)
    {
        const int n1 = n - 1;
        VectorXli v, w;
        NTL::mat_ZZ c;
        VectorXld B(n), s;
        B.setZero();
        MatrixXld mu(n, n);
        mu.setIdentity();

        GSO(B, mu, n, m);

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

                v = w * basis.block(k - 1, 0, lk1, m);

                // Inserts and removes linearly-dependentness by LLL-reducing
                c.SetDims(h + 1, m);
                for (j = 0; j < m; ++j)
                {
                    for (i = 0; i < k - 1; ++i)
                        c[i][j] = basis.coeff(i, j);
                    c[k - 1][j] = v.coeff(j);
                    for (i = k; i < h + 1; ++i)
                        c[i][j] = basis.coeff(i - 1, j);
                }
                NTL::LLL(_, c, 99, 100);
                for (i = 1; i <= h; ++i)
                    for (j = 0; j < m; ++j)
                        basis.coeffRef(i - 1, j) = NTL::to_long(c[i][j]);
                GSO(B, mu, n, m);
            }
            else
            {
                ++z;

                c.SetDims(h, m);
                for (j = 0; j < m; ++j)
                    for (i = 0; i < h; ++i)
                        c[i][j] = basis.coeff(i, j);
                NTL::LLL(_, c, 99, 100);
                for (i = 0; i < h; ++i)
                    for (j = 0; j < m; ++j)
                        basis.coeffRef(i, j) = NTL::to_long(c[i][j]);
                GSO(B, mu, n, m);
            }
        }
    }

    void PotBKZ_(const int beta, const double d, const int n, const int m, FILE *fp);

    void DualPotBKZ_(const int beta, const double delta, const int n, const int m, FILE *fp);

    void SelfDualPotBKZ_each_(const int beta, const double d, const int n, const int m, FILE *fp);
};

#endif // !SDPB_HPP