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
    long double m_temp;
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
    void updateDeepInsGSO(const int i, const int k, VectorXld &B, MatrixXld &mu, const int n);

    void updateDualDeepInsGSO(const int k, const int l, VectorXld &B, MatrixXld &mu, MatrixXld &dual_mu, const VectorXld dual_D, const int n);

public:
    MatrixXli basis;

    /**
     * @brief set numbers of rows and columns
     *
     * @param n number of rows to set
     * @param m number of columns to set
     */
    void setDims(const int n, const int m);

    /**
     * @brief dual version of deep-insertion
     *
     * @param n number of rows of lattice basis matrix
     * @param k
     * @param l
     */
    void dualDeepInsertion(const int n, const int k, const int l);

    /**
     * @brief insert a vector into dual lattice basis
     *
     * @param b basis
     * @param x coefficient vector to insert
     * @param n number of rows of lattice basis matrix
     * @param m number of columns of lattice basis matrix
     * @return MatrixXli new basis whose dual basis was inserted vector
     */
    MatrixXli Insert(const MatrixXli b, const VectorXli x, const int n, const int m);

    /**
     * @brief computes GSO-informations
     *
     * @param B vector of squared norm of GSO-vectors
     * @param mu GSO-coefficient-matrix
     * @param n number of rows of lattice basis matrix
     * @param m number of columns of lattice basis matrix
     */
    void GSO(VectorXld &B, MatrixXld &mu, const int n, const int m);

    /**
     * @brief computes GSO-informations
     *
     * @param B vector of squared norm of GSO-vectors
     * @param logB vector of logarithm value of squared norm of GSO-vectors
     * @param mu GSO-coefficient-matrix
     * @param n number of rows of lattice basis matrix
     * @param m number of columns of lattice basis matrix
     */
    void GSO(VectorXld &B, VectorXld &logB, MatrixXld &mu, const int n, const int m);

    /// @brief Computes GSO-informations of dual-lattice basis using GSO-information of lattice
    /// @param B Squared norms of GSO-vectors
    /// @param mu GSO-coefficient mattix
    /// @param C Squared norms of dual-GSO-vectors
    /// @param hmu dual-GSO-coefficient mattix
    /// @param n number of rows of lattice basis matrix
    /// @param m number of columns of lattice basis matrix
    void DualGSO(const VectorXld B, const VectorXld logB, const MatrixXld mu, VectorXld &C, VectorXld &logC, MatrixXld &hmu, const int n, const int m);

    /// @brief Conputes log value of potential of basis
    /// @param B Squared norms of GSO-vectors
    /// @param n number of rows of lattice basis matrix
    /// @return double logarithm value of potential of basis
    long double logPot(const VectorXld B, const int n);

    /**
     * @brief computes slope of GSA-slope
     *
     * @param B vector of squared norm of GSO-vectors
     * @param n number of rows of lattice basis matrix
     * @param m number of columns of lattice basis matrix
     * @return long double slope of GSA-slope
     */
    long double rho(const VectorXld B, const int n, const int m);

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
    bool ENUM(VectorXli &u, const MatrixXld mu, const VectorXld B, VectorXld &rho, const int n, double R);

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

extern "C" long **generator(long** basis, const int n, const int seed);

#endif // !LATTICE_H
