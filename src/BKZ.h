#ifndef BKZ_H
#define BKZ_H

#include "Lattice.h"

inline bool Lattice::ENUM(VectorXli &u, const MatrixXld mu, const VectorXld B, VectorXld &rho, const int n, double R)
{
    bool has_solution = false;
    int i, r[n + 1];
    int last_nonzero = 0;                              // index of last non-zero elements
    VectorXli weight = VectorXli::Zero(n);             // weight of zigzag searching
    VectorXli coeff_vector = VectorXli::Zero(n);       // coefficient vector to putput
    Eigen::VectorXd center = Eigen::VectorXd::Zero(n); // center of zigzag serching
    Eigen::MatrixXd sigma = Eigen::MatrixXd::Zero(n + 1, n);
    u = VectorXli::Zero(n);
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
                has_solution = true;
                u = coeff_vector;
                R = fmin(R, 0.99 * rho.coeff(0));
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
                return has_solution;
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

inline void Lattice::BKZ_(const int block_size, const double reduction_parameter, const int max_loop, const int n, const int m, FILE *potential_file)
{
    int consecutive_non_insertion_count = 0; // number of consecutive not inserting vector to lattice basis
    int dim_of_local_block_lattice;          // dimesion of local projected block lattice
    VectorXli shortest_vec;                  // shortest vector on local projected block lattice
    VectorXli coeff_vec;                     // coefficients of shortest vector on local projected block lattice
    NTL::mat_ZZ inserted_vecs;               // vectors of lattice basis vector and inserted shortest vector
    VectorXld B = VectorXld::Zero(n), s;
    MatrixXld mu = MatrixXld::Identity(n, n);
    double b1_norm = basis.row(0).cast<double>().norm();

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

        if(basis.row(0).cast<double>().norm() < b1_norm)
        {
            b1_norm = basis.row(0).cast<double>().norm();
            printf("%d tours: A shorter vector found: %lf\n", m_tours_of_bkz, b1_norm);
        }

        //coeff_vec = enumerate(mu.block(k - 1, k - 1, dim_of_local_block_lattice, dim_of_local_block_lattice), B.segment(k - 1, dim_of_local_block_lattice), s, dim_of_local_block_lattice);
        if (ENUM(coeff_vec, mu.block(k - 1, k - 1, dim_of_local_block_lattice, dim_of_local_block_lattice), B.segment(k - 1, dim_of_local_block_lattice), s, dim_of_local_block_lattice, 0.99 * B.coeff(k - 1)))
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

            NTL::LLL_FP(inserted_vecs, reduction_parameter);

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

            NTL::LLL_FP(inserted_vecs, reduction_parameter);
            //NTL::LLL(_, inserted_vecs, 99, 100);

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