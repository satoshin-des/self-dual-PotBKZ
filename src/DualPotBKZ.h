#ifndef DUAL_POT_BKZ_H
#define DUAL_POT_BKZ_H

#include "Lattice.h"
#include "DualPotENUM.h"
#include "DualPotLLL.h"
#include "PotLLL.h"

inline void Lattice::DualPotBKZ_(const int block_size, const double reduction_parameter, const int n, const int m, FILE *potential_file)
{
    int consecutive_solution_count = 0; // consecutive numbers of that DualPotENUM has solution
    int dim_of_local_block_lattice;     // dimension of local projected block lattice
    VectorXli enum_coeff_vector;        // coefficient vector enumerated by DualPotENUM
    VectorXld B(n), logB(n), C, logC;
    MatrixXld mu(n, n), hmu, BB;
    MatrixXli temp_basis;
    B.setZero();
    logB.setZero();
    mu.setZero();

    GSO(B, logB, mu, n, m);
    fprintf(potential_file, "%Lf\n", logPot(B, n));

    for (int j = n, k; consecutive_solution_count < n - 1;)
    {
        if (j == 1)
        {
            j = n;
        }
        --j;
        k = (j - block_size + 1 > 0 ? j - block_size + 1 : 0);
        dim_of_local_block_lattice = j - k + 1;

        fprintf(potential_file, "%Lf\n", logPot(B, n));

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

            temp_basis = Insert(basis.block(k, 0, dim_of_local_block_lattice, m), enum_coeff_vector, dim_of_local_block_lattice, m);
            basis.block(k, 0, dim_of_local_block_lattice, m) = temp_basis.block(0, 0, dim_of_local_block_lattice, m);

            DualPotLLL_(reduction_parameter, n, m);
            GSO(B, logB, mu, n, m);
        }
    }
}

#endif // !DUAL_POT_BKZ_H