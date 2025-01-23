#ifndef LATTICE_H
#define LATTICE_H

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

/**
 * @brief class for lattice basis reduction
 *
 */
class Lattice
{
private:
    int m_tours_of_bkz = 0;
    NTL::ZZ _;

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
        {
            basis.row(j) = basis.row(j + 1);
        }
        basis.row(l) = t;
    }

    /**
     * @brief insert a vector into dual lattice basis
     *
     * @param b basis
     * @param x coefficient vector to insert
     * @param n number of rows of lattice basis matrix
     * @param m number of columns of lattice basis matrix
     * @return MatrixXli new basis whose dual basis was inserted vector
     */
    MatrixXli Insert(const MatrixXli b, const VectorXli x, const int n, const int m)
    {
        int i, j;
        const double beta = 1.35135135135135135135135135135135135135135135135135;
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
            {
                tmp_base[i][j] = 0;
            }
            tmp_base[i][i] = 1;
            tmp_base[i][n] = gamma * x.coeff(i);
        }

        NTL::LLL(_, tmp_base, 99, 100);

        for (i = 0; i < n; ++i)
        {
            for (j = 0; j < n; ++j)
            {
                U.coeffRef(i, j) = NTL::to_long(tmp_base[i][j]);
            }
        }

        return U * b;
    }

    /**
     * @brief computes GSO-informations
     *
     * @param B vector of squared norm of GSO-vectors
     * @param mu GSO-coefficient-matrix
     * @param n number of rows of lattice basis matrix
     * @param m number of columns of lattice basis matrix
     */
    void GSO(VectorXld &B, MatrixXld &mu, const int n, const int m)
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

    /**
     * @brief computes GSO-informations
     *
     * @param B vector of squared norm of GSO-vectors
     * @param logB vector of logarithm value of squared norm of GSO-vectors
     * @param mu GSO-coefficient-matrix
     * @param n number of rows of lattice basis matrix
     * @param m number of columns of lattice basis matrix
     */
    void GSO(VectorXld &B, VectorXld &logB, MatrixXld &mu, const int n, const int m)
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

    /// @brief Computes GSO-informations of dual-lattice basis using GSO-information of lattice
    /// @param B Squared norms of GSO-vectors
    /// @param mu GSO-coefficient mattix
    /// @param C Squared norms of dual-GSO-vectors
    /// @param hmu dual-GSO-coefficient mattix
    /// @param n number of rows of lattice basis matrix
    /// @param m number of columns of lattice basis matrix
    void DualGSO(const VectorXld B, const VectorXld logB, const MatrixXld mu, VectorXld &C, VectorXld &logC, MatrixXld &hmu, const int n, const int m)
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

    /// @brief Conputes log value of potential of basis
    /// @param B Squared norms of GSO-vectors
    /// @param n number of rows of lattice basis matrix
    /// @return double logarithm value of potential of basis
    long double logPot(const VectorXld B, const int n)
    {
        NTL::RR logp;
        logp = 0;
        for (int i = 0; i < n; ++i)
        {
            logp += (n - i) * NTL::log(NTL::to_RR((double)B.coeff(i)));
        }
        return NTL::to_double(logp);
    }

    /**
     * @brief computes slope of GSA-slope
     *
     * @param B vector of squared norm of GSO-vectors
     * @param n number of rows of lattice basis matrix
     * @param m number of columns of lattice basis matrix
     * @return long double slope of GSA-slope
     */
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

    /**
     * @brief Enumerates a vector whose norm is shorter than R on the lattice
     * 
     * @param mu GSO-coefficient-matrix
     * @param B vector of squared norm of GSO-vectors
     * @param rho projected vector
     * @param n number of rows of lattice basis matrix
     * @param R upper bounds of enumeration
     * @return VectorXli vector whose norm is shorter that given upper bounds
     */
    VectorXli ENUM(const MatrixXld mu, const VectorXld B, VectorXld &rho, const int n, const double R);

    /// @brief Enumerates the shortest vector on the lattice
    /// @param mu GSO-coefficient matrix
    /// @param B Squared-norms of row vectors of GSO-matrix
    /// @param n Rank of lattice
    /// @return VectorXli the shortest vector
    VectorXli enumerate(const MatrixXld mu, const VectorXld B, VectorXld &rho, const int n);

    VectorXli PotENUM(const MatrixXld mu, const VectorXld B, const VectorXld logB, const int n);

    VectorXli DualPotENUM(const MatrixXld mu, const VectorXld B, const VectorXld logB, const int n);

    void PotLLL_(const long double d, const int n, const int m);

    void DualPotLLL_(const double d, const int n, const int m);

    void BKZ_(const int block_size, const double reduction_parameter, const int max_loop, const int n, const int m, FILE *potential_file);

    void PotBKZ_(const int beta, const double d, const int n, const int m, FILE *fp);

    void DualPotBKZ_(const int beta, const double delta, const int n, const int m, FILE *fp);

    void SelfDualPotBKZ_(const int beta, const double d, const int n, const int m, FILE *fp);
};

#endif // !LATTICE_H