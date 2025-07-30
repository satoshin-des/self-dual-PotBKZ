#include "Lattice.h"

#include <iostream>
#include <cmath>

#include <eigen3/Eigen/Dense>

#include <NTL/ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>

void Lattice::PotBKZ_(const int block_size, const double reduction_parameter, const int n, const int m, FILE *potential_file)
{
    int n_tour = 0;
    const int n1 = n - 1, n2 = n - 2;
    VectorXli v, w;
    MatrixXli tmp_b(n, n);
    VectorXld B(n), logB(n);
    MatrixXld mu(n, n);
    NTL::mat_ZZ cc;
    double b1_norm = basis.row(0).cast<double>().norm();

    GSO(B, logB, mu, n, m);
    fprintf(potential_file, "%Lf\n", logPot(B, n));
    for (int z = 0, j = 0, i, k, l, kj1; z < n - 1;)
    {
        fprintf(potential_file, "%Lf\n", logPot(B, n));

        if (j == n2)
        {
            j = 0;
            ++n_tour;
        }
        ++j;
        k = (j + block_size - 1 < n1 ? j + block_size - 1 : n1);
        kj1 = k - j + 1;
        v.resize(kj1);
        v.setZero();

        /* enumerate a shortest vector*/
        v = PotENUM(mu.block(j, j, kj1, kj1), B.segment(j, kj1), logB.segment(j, kj1), kj1);

        if(basis.row(0).cast<double>().norm() < b1_norm)
        {
            b1_norm = basis.row(0).cast<double>().norm();
            printf("%d tours: A shorter vector found: %lf\n", n_tour, b1_norm);
        }

        if (!v.isZero())
        {
            z = 0;

            w = v * basis.block(j, 0, kj1, m);
            cc.SetDims(n + 1, m);
            for (l = 0; l < m; ++l)
            {
                for (i = 0; i < j; ++i)
                {
                    cc[i][l] = basis.coeffRef(i, l);
                }
                cc[j][l] = w[l];
                for (i = j + 1; i < n + 1; ++i)
                {
                    cc[i][l] = basis.coeffRef(i - 1, l);
                }
            }

            NTL::LLL_FP(cc, 0.99);

            for (i = 0; i < n; ++i)
            {
                for (l = 0; l < m; ++l)
                {
                    basis.coeffRef(i, l) = NTL::to_long(cc[i + 1][l]);
                }
            }

            PotLLL_(reduction_parameter, n, m);
            GSO(B, logB, mu, n, m);
        }
        else
        {
            ++z;
        }
    }
}
