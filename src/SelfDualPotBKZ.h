#ifndef SELF_DUAL_POT_BKZ_H
#define SELF_DUAL_POT_BKZ_H

#include "Lattice.h"
#include "PotENUM.h"
#include "DualPotENUM.h"
#include "PotLLL.h"
#include "DualPotLLL.h"

inline void Lattice::SelfDualPotBKZ_(const int block_size, const double reduction_parameter, const int n, const int m, FILE *fp)
{
    bool is_primal = true;                     // current tour is primal part or not
    int primal_consecutive_solution_count = 0; // consecutive numbers of that PotENUM has solution
    int dual_consecutive_solution_count = 0;   // consecutive numbers of that DualPotENUM has solution
    int dim_of_local_block_lattice;            // dimension of local projected block lattice
    VectorXli v;                               // enumerated vector by PotENUM or DualPotENUM
    VectorXli w;                               // enumerated coefficient vector by PotENUM or DualPotENUM
    MatrixXli temp_basis(n, n);
    VectorXld B(n), logB(n);
    MatrixXld mu(n, n);
    NTL::mat_ZZ inserted_vecs;
    VectorXld C, logC;
    MatrixXld hmu, BB;

    GSO(B, logB, mu, n, m);
    fprintf(fp, "%Lf\n", logPot(B, n));

    DualPotLLL_(0.99, n, m);

    for (int primal_j = 0, i, k, l, dual_j = n; std::max(primal_consecutive_solution_count, dual_consecutive_solution_count) < n - 1;)
    {
        if (is_primal)
        { // primal part
            if (primal_j == n - 2)
            {
                primal_j = 0;
                is_primal = false;
            }
            ++primal_j;
            k = (primal_j + block_size - 1 < n - 1 ? primal_j + block_size - 1 : n - 1);
            dim_of_local_block_lattice = k - primal_j + 1;

            fprintf(fp, "%Lf\n", logPot(B, n));

            v.resize(dim_of_local_block_lattice);
            v.setZero();

            /* enumerate a shortest vector*/
            v = PotENUM(mu.block(primal_j, primal_j, dim_of_local_block_lattice, dim_of_local_block_lattice), B.segment(primal_j, dim_of_local_block_lattice), logB.segment(primal_j, dim_of_local_block_lattice), dim_of_local_block_lattice);

            if (!v.isZero())
            {
                primal_consecutive_solution_count = 0;

                w = v * basis.block(primal_j, 0, dim_of_local_block_lattice, m);
                inserted_vecs.SetDims(n + 1, m);
                for (l = 0; l < m; ++l)
                {
                    for (i = 0; i < primal_j; ++i)
                    {
                        inserted_vecs[i][l] = basis.coeffRef(i, l);
                    }
                    inserted_vecs[primal_j][l] = w[l];
                    for (i = primal_j + 1; i < n + 1; ++i)
                    {
                        inserted_vecs[i][l] = basis.coeffRef(i - 1, l);
                    }
                }

                NTL::LLL(_, inserted_vecs, 99, 100);

                for (i = 0; i < n; ++i)
                {
                    for (l = 0; l < m; ++l)
                    {
                        basis.coeffRef(i, l) = NTL::to_long(inserted_vecs[i + 1][l]);
                    }
                }

                DualPotLLL_(reduction_parameter, n, m);
                GSO(B, logB, mu, n, m);
            }
            else
            {
                ++primal_consecutive_solution_count;
            }
        }
        else
        { // dual part
            if (dual_j == 1)
            {
                dual_j = n;
                is_primal = true;
            }
            --dual_j;
            k = (dual_j - block_size + 1 > 0 ? dual_j - block_size + 1 : 0);
            dim_of_local_block_lattice = dual_j - k + 1;

            fprintf(fp, "%Lf\n", logPot(B, n));

            C.resize(dim_of_local_block_lattice);
            logC.resize(dim_of_local_block_lattice);
            hmu.resize(dim_of_local_block_lattice, dim_of_local_block_lattice);
            DualGSO(B.segment(k, dim_of_local_block_lattice),
                    logB.segment(k, dim_of_local_block_lattice),
                    mu.block(k, k, dim_of_local_block_lattice, dim_of_local_block_lattice),
                    C, logC, hmu,
                    dim_of_local_block_lattice, dim_of_local_block_lattice);

            // Dual Enumeration
            v = DualPotENUM(hmu, C, logC, dim_of_local_block_lattice);

            if (v.isZero())
            {
                ++dual_consecutive_solution_count;
            }
            else
            {
                dual_consecutive_solution_count = 0;

                temp_basis = Insert(basis.block(k, 0, dim_of_local_block_lattice, m), v, dim_of_local_block_lattice, m);
                basis.block(k, 0, dim_of_local_block_lattice, m) = temp_basis.block(0, 0, dim_of_local_block_lattice, m);

                PotLLL_(reduction_parameter, n, m);
                GSO(B, logB, mu, n, m);
            }
        }
    }
}

#endif // !SELF_DUAL_POT_BKZ_H