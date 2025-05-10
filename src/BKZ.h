#ifndef BKZ_H
#define BKZ_H

#include "Lattice.h"

inline VectorXli Lattice::ENUM(const MatrixXld mu, const VectorXld B, VectorXld &rho, const int n, const double R)
{
    int i, r[n + 1];
    int last_nonzero = 0;                              // index of last non-zero elements
    VectorXli weight = VectorXli::Zero(n);             // weight of zigzag searching
    VectorXli coeff_vector = VectorXli::Zero(n);       // coefficient vector to putput
    Eigen::VectorXd center = Eigen::VectorXd::Zero(n); // center of zigzag serching
    Eigen::MatrixXd sigma = Eigen::MatrixXd::Zero(n + 1, n);
    coeff_vector.coeffRef(0) = 1;
    rho.setZero();
    for (i = 0; i < n; ++i)
    {
        r[i] = i;
    }

    for (int k = 0;;)
    {
        m_temp = static_cast<double>(coeff_vector.coeff(k)) - center.coeff(k);
        m_temp *= m_temp;
        rho.coeffRef(k) = rho.coeff(k + 1) + m_temp * B.coeff(k); // rho[k]=∥πₖ(shortest_vec)∥
        if (rho.coeff(k) <= R)
        {
            if (k == 0)
            {
                return coeff_vector;
            }
            else
            {
                --k;
                r[k] = (r[k] > r[k + 1] ? r[k] : r[k + 1]);
                for (i = r[k]; i > k; --i)
                {
                    sigma.coeffRef(i, k) = sigma.coeff(i + 1, k) + mu.coeff(i, k) * coeff_vector.coeff(i);
                }
                center.coeffRef(k) = -sigma.coeff(k + 1, k);
                coeff_vector.coeffRef(k) = round(center.coeff(k));
                weight.coeffRef(k) = 1; // 小さいやつが見つかったら、変分を元に戻す
            }
        }
        else
        {
            ++k;
            if (k == n)
            { // no solution
                return VectorXli::Zero(n);
            }
            else
            {
                r[k] = k;
                if (k >= last_nonzero)
                {
                    last_nonzero = k;
                    ++coeff_vector.coeffRef(k);
                }
                else
                {
                    coeff_vector.coeff(k) > center.coeff(k) ? coeff_vector.coeffRef(k) -= weight.coeff(k) : coeff_vector.coeffRef(k) += weight.coeff(k);
                    ++weight.coeffRef(k);
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
inline VectorXli Lattice::enumerate(const MatrixXld mu, const VectorXld B, VectorXld &rho, const int n)
{
    VectorXli enum_vector = VectorXli::Zero(n);
    VectorXli old_enum_vector = VectorXli::Zero(n);
    VectorXld old_rho = VectorXld::Zero(n + 1);

    for (double R = B.coeff(0);;)
    {
        old_rho = rho;
        old_enum_vector = enum_vector;
        enum_vector = ENUM(mu, B, rho, n, R);
        if (enum_vector.isZero())
        {
            rho = old_rho;
            return old_enum_vector;
        }
        R *= 0.99;
    }
}

inline void Lattice::BKZ_(const int block_size, const double reduction_parameter, const int max_loop, const int n, const int m, FILE *potential_file)
{
    int consecutive_non_insertion_count = 0; // number of consecutive not inserting vector to lattice basis
    int dim_of_local_block_lattice;          // dimesion of local projected block lattice
    VectorXli shortest_vec;                  // shortest vector on local projected block lattice
    VectorXli coeff_vec;                     // coefficients of shortest vector on local projected block lattice
    NTL::mat_ZZ inserted_vecs;               // vectors of lattice basis vector and inserted shortest vector
    VectorXld B = VectorXld::Zero(n), s;
    MatrixXld mu = MatrixXld::Identity(n, n);

    GSO(B, mu, n, m);

    for (int j, i, k = 0, h, l; consecutive_non_insertion_count < n - 1;)
    {
        if (m_tours_of_bkz >= max_loop)
        {
            break;
        }

        fprintf(potential_file, "%Lf\n", logPot(B, n));

        if (k == n - 1)
        {
            k = 0;
            ++m_tours_of_bkz;
        }
        ++k;
        l = std::min(k + block_size - 1, n);
        h = std::min(l + 1, n);
        dim_of_local_block_lattice = l - k + 1;

        s.resize(dim_of_local_block_lattice + 1);
        s.setZero();
        coeff_vec = enumerate(mu.block(k - 1, k - 1, dim_of_local_block_lattice, dim_of_local_block_lattice), B.segment(k - 1, dim_of_local_block_lattice), s, dim_of_local_block_lattice);
        if (B.coeff(k - 1) > s.coeff(k - 1) && (!coeff_vec.isZero()))
        {
            consecutive_non_insertion_count = 0;

            shortest_vec = coeff_vec * basis.block(k - 1, 0, dim_of_local_block_lattice, m);

            // Inserts and removes linearly-dependentness by LLL-reducing
            inserted_vecs.SetDims(h + 1, m);
            for (j = 0; j < m; ++j)
            {
                for (i = 0; i < k - 1; ++i)
                {
                    inserted_vecs[i][j] = basis.coeff(i, j);
                }
                inserted_vecs[k - 1][j] = shortest_vec.coeff(j);
                for (i = k; i < h + 1; ++i)
                {
                    inserted_vecs[i][j] = basis.coeff(i - 1, j);
                }
            }

            NTL::LLL(_, inserted_vecs, 99, 100);

            for (i = 1; i <= h; ++i)
            {
                for (j = 0; j < m; ++j)
                {
                    basis.coeffRef(i - 1, j) = NTL::to_long(inserted_vecs[i][j]);
                }
            }
            GSO(B, mu, n, m);
        }
        else
        {
            ++consecutive_non_insertion_count;

            inserted_vecs.SetDims(h, m);
            for (j = 0; j < m; ++j)
            {
                for (i = 0; i < h; ++i)
                {
                    inserted_vecs[i][j] = basis.coeff(i, j);
                }
            }

            NTL::LLL(_, inserted_vecs, 99, 100);

            for (i = 0; i < h; ++i)
            {
                for (j = 0; j < m; ++j)
                {
                    basis.coeffRef(i, j) = NTL::to_long(inserted_vecs[i][j]);
                }
            }

            GSO(B, mu, n, m);
        }
    }
}

#endif // !BKZ_H