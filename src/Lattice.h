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
    double m_temp;
    NTL::ZZ _;

    /**
     * @brief Updates GSO-informations by deep-insertion without computing GSO directly
     *
     * @param i
     * @param k
     * @param B vector of squared norm of GSO-vectors
     * @param mu GSO-coefficient-matrix
     * @param n number of rows of lattice basis matrix
     */
    void updateDeepInsGSO(const int i, const int k, VectorXld &B, MatrixXld &mu, const int n)
    {
        int j, l;
        double T, eps;
        VectorXld P(n), D(n), S(n);
        P.setZero();
        D.setZero();
        S.setZero();

        P.coeffRef(k) = D.coeffRef(k) = B.coeff(k);
        for (j = k - 1; j >= i; --j)
        {
            P.coeffRef(j) = mu.coeff(k, j) * B.coeff(j);
            D.coeffRef(j) = D.coeff(j + 1) + mu.coeff(k, j) * P.coeff(j);
        }

        for (j = k; j > i; --j)
        {
            T = mu.coeff(k, j - 1) / D.coeff(j);
            for (l = n - 1; l > k; --l)
            {
                S.coeffRef(l) += mu.coeff(l, j) * P.coeff(j);
                mu.coeffRef(l, j) = mu.coeff(l, j - 1) - T * S.coeff(l);
            }
            for (l = k; l > j; --l)
            {
                S.coeffRef(l) += mu.coeff(l - 1, j) * P.coeff(j);
                mu.coeffRef(l, j) = mu.coeff(l - 1, j - 1) - T * S.coeff(l);
            }
        }

        T = 1.0 / D.coeff(i);

        for (l = n - 1; l > k; --l)
        {
            mu.coeffRef(l, i) = T * (S.coeff(l) + mu.coeff(l, i) * P.coeff(i));
        }
        for (l = k; l >= i + 2; --l)
        {
            mu.coeffRef(l, i) = T * (S.coeff(l) + mu.coeff(l - 1, i) * P.coeff(i));
        }

        mu.coeffRef(i + 1, i) = T * P.coeff(i);
        for (j = 0; j < i; ++j)
        {
            eps = mu.coeff(k, j);
            for (l = k; l > i; --l)
            {
                mu.coeffRef(l, j) = mu.coeff(l - 1, j);
            }
            mu.coeffRef(i, j) = eps;
        }

        for (j = k; j > i; --j)
        {
            B.coeffRef(j) = D.coeff(j) * B.coeff(j - 1) / D.coeff(j - 1);
        }
        B.coeffRef(i) = D.coeff(i);
    }

    void updateDualDeepInsGSO(const int k, const int l, VectorXld &B, MatrixXld &mu, MatrixXld &dual_mu, const VectorXld dual_D, const int n)
    {
        int i, j, h;
        double sum;
        MatrixXld xi = mu;

        for (i = l + 1; i < n; ++i)
        {
            sum = 0;
            for (h = k; h <= l; ++h)
            {
                sum += dual_mu.coeff(k, h) * mu.coeff(i, h);
            }
            xi.coeffRef(i, l) = sum;
        }

        for (j = k; j < l; ++j)
        {
            for (i = j + 1; i < l; ++i)
            {
                sum = 0;
                for (h = k; h <= j; ++h)
                {
                    sum += dual_mu.coeff(k, h) * mu.coeff(i + 1, h);
                }

                xi.coeffRef(i, j) = mu.coeff(i + 1, j + 1) * dual_D.coeff(j) / dual_D.coeff(j + 1) - dual_mu.coeff(k, j + 1) / (dual_D.coeff(j + 1) * B.coeff(j + 1)) * sum;
            }
            xi.coeffRef(l, j) = -dual_mu.coeff(k, j + 1) / (dual_D.coeff(j + 1) * B.coeff(j + 1));
            for (i = l + 1; i < n; ++i)
            {
                sum = 0;
                for (h = k; h <= j; ++h)
                {
                    sum += dual_mu.coeff(k, h) * mu.coeff(i, h);
                }

                xi.coeffRef(i, j) = mu.coeff(i, j + 1) * dual_D.coeff(j) / dual_D.coeff(j + 1) - dual_mu.coeff(k, j + 1) / (dual_D.coeff(j + 1) * B.coeff(j + 1)) * sum;
            }
        }

        for (j = 0; j < k; ++j)
        {
            for (i = k; i < l; ++i)
            {
                xi.coeffRef(i, j) = mu.coeff(i + 1, j);
            }
            xi.coeffRef(l, j) = mu.coeff(k, j);
        }

        mu = xi;
        for (j = k; j < l; ++j)
        {
            B.coeffRef(j) = dual_D.coeff(j + 1) * B.coeff(j + 1) / dual_D.coeff(j);
        }
        B.coeffRef(l) = 1 / dual_D.coeff(l);
    }

public:
    MatrixXli basis;

    /**
     * @brief set numbers of rows and columns
     *
     * @param n number of rows to set
     * @param m number of columns to set
     */
    void setDims(const int n, const int m)
    {
        basis.resize(n, m);
    }

    /**
     * @brief dual version of deep-insertion
     *
     * @param n number of rows of lattice basis matrix
     * @param k
     * @param l
     */
    void dualDeepInsertion(const int n, const int k, const int l)
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
        double gamma;
        MatrixXli U(n, n);
        U.setZero();
        NTL::mat_ZZ temp_basis;
        temp_basis.SetDims(n, n + 1);

        /* Construction of gamma */
        m_temp = x.cast<double>().norm();
        m_temp *= pow(beta, (n - 2) * 0.5);
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

        NTL::LLL(_, temp_basis, 99, 100);

        for (i = 0; i < n; ++i)
        {
            for (j = 0; j < n; ++j)
            {
                U.coeffRef(i, j) = NTL::to_long(temp_basis[i][j]);
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
    /// @param n number of rows of lattice basis matrix
    /// @return VectorXli the shortest vector
    VectorXli enumerate(const MatrixXld mu, const VectorXld B, VectorXld &rho, const int n);

    /**
     * @brief PotENUM algorithm
     *
     * @param mu GSO-coefficient-matrix
     * @param B vector of squared norm of GSO-vectors
     * @param logB vector of logarithm value of squared norm of GSO-vectors
     * @param n number of rows of lattice basis matrix
     * @return VectorXli delta-anomalous vector
     */
    VectorXli PotENUM(const MatrixXld mu, const VectorXld B, const VectorXld logB, const int n);

    /**
     * @brief DualPotENUM algorithm
     *
     * @param mu GSO-coefficient-matrix
     * @param B vector of squared norm of GSO-vectors
     * @param logB vector of logarithm value of squared norm of GSO-vectors
     * @param n number of rows of lattice basis matrix
     * @return VectorXli vector that decreases potential of lattice inserting to dual basis
     */
    VectorXli DualPotENUM(const MatrixXld mu, const VectorXld B, const VectorXld logB, const int n);

    /**
     * @brief PotLLL reduction
     *
     * @param reduction_parameter reduction parameter
     * @param n number of rows of lattice basis matrix
     * @param m number of columns of lattice basis matrix
     */
    void PotLLL_(const double reduction_parameter, const int n, const int m);

    /**
     * @brief dual version of PotLLL reduction
     *
     * @param reduction_parameter reduction parameter
     * @param n number of rows of lattice basis matrix
     * @param m number of columns of lattice basis matrix
     */
    void DualPotLLL_(const double reduction_parameter, const int n, const int m);

    /**
     * @brief BKZ reduction
     *
     * @param block_size block size
     * @param reduction_parameter reduction parameter
     * @param max_loop limit number of tours
     * @param n number of rows of lattice basis matrix
     * @param m number of columns of lattice basis matrix
     * @param potential_file file to output logarithm values of potential
     */
    void BKZ_(const int block_size, const double reduction_parameter, const int max_loop, const int n, const int m, FILE *potential_file);

    /**
     * @brief PotBKKZ reduction
     *
     * @param block_size block size
     * @param reduction_parameter reduction parameter
     * @param n number of rows of lattice basis matrix
     * @param m number of columns of lattice basis matrix
     * @param potential_file file to output logarithm values of potential
     */
    void PotBKZ_(const int block_size, const double reduction_parameter, const int n, const int m, FILE *potential_file);

    /**
     * @brief dual version of PotBKZ algorithm
     *
     * @param block_size block size
     * @param reduction_parameter reduction parameter
     * @param n number of rows of lattice basis matrix
     * @param m number of columns of lattice basis matrix
     * @param potential_file file to output logarithm values of potential
     */
    void DualPotBKZ_(const int block_size, const double reduction_parameter, const int n, const int m, FILE *potential_file);

    /**
     * @brief self-dual version of PotBKZ algorithm
     *
     * @param block_size block size
     * @param reduction_parameter reduction parameterr
     * @param n number of rows of lattice basis matrix
     * @param m number of columns of lattice basis matrix
     * @param potential_file file to output logarithm values of potential
     */
    void SelfDualPotBKZ_(const int block_size, const double reduction_parameter, const int n, const int m, FILE *potential_file);
};

#endif // !LATTICE_H