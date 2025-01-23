#ifndef DUAL_POT_BKZ_H
#define DUAL_POT_BKZ_H

#include "Lattice.h"
#include "DualPotENUM.h"
#include "DualPotLLL.h"
#include "PotLLL.h"

inline void Lattice::DualPotBKZ_(const int beta, const double delta, const int n, const int m, FILE *fp)
{
    int consecutive_solution_count = 0; // consecutive numbers of that DualPotENUM has solution
    int dim_of_local_block_lattice;     // dimension of local projected block lattice
    VectorXli enum_coeff_vector;        // coefficient vector enumerated by DualPotENUM
    MatrixXli tmp_b;
    VectorXld B(n), logB(n), C, logC;
    MatrixXld mu(n, n), hmu, BB;
    B.setZero();
    logB.setZero();
    mu.setZero();

    GSO(B, logB, mu, n, m);
    fprintf(fp, "%Lf\n", logPot(B, n));

    for (int j = n, k; consecutive_solution_count < n - 1;)
    {
        if (j == 1)
        {
            j = n;
        }
        --j;
        k = (j - beta + 1 > 0 ? j - beta + 1 : 0);
        dim_of_local_block_lattice = j - k + 1;

        fprintf(fp, "%Lf\n", logPot(B, n));

        C.resize(dim_of_local_block_lattice);
        logC.resize(dim_of_local_block_lattice);
        hmu.resize(dim_of_local_block_lattice, dim_of_local_block_lattice);
        DualGSO(B.segment(k, dim_of_local_block_lattice),
                logB.segment(k, dim_of_local_block_lattice),
                mu.block(k, k, dim_of_local_block_lattice, dim_of_local_block_lattice),
                C, logC, hmu,
                dim_of_local_block_lattice,
                dim_of_local_block_lattice);

        enum_coeff_vector = DualPotENUM(hmu, C, logC, dim_of_local_block_lattice);

        if (enum_coeff_vector.isZero())
        {
            ++consecutive_solution_count;
        }
        else
        {
            consecutive_solution_count = 0;

            tmp_b = Insert(basis.block(k, 0, dim_of_local_block_lattice, m), enum_coeff_vector, dim_of_local_block_lattice, m);
            basis.block(k, 0, dim_of_local_block_lattice, m) = tmp_b.block(0, 0, dim_of_local_block_lattice, m);

            DualPotLLL_(delta, n, m);
            GSO(B, logB, mu, n, m);
        }
    }
}

#endif // !DUAL_POT_BKZ_H