#ifndef POT_BKZ_H
#define POT_BKZ_H

#include "Lattice.h"
#include "PotENUM.h"
#include "PotLLL.h"

inline void Lattice::PotBKZ_(const int block_size, const double reduction_parameter, const int n, const int m, FILE *fp)
{
    const int n1 = n - 1, n2 = n - 2;
    VectorXli v, w;
    MatrixXli tmp_b(n, n);
    VectorXld B(n), logB(n);
    MatrixXld mu(n, n);
    NTL::mat_ZZ cc;

    GSO(B, logB, mu, n, m);
    fprintf(fp, "%Lf\n", logPot(B, n));
    for (int z = 0, j = 0, i, k, l, kj1; z < n - 1;)
    {
        fprintf(fp, "%Lf\n", logPot(B, n));

        if (j == n2)
        {
            j = 0;
        }
        ++j;
        k = (j + block_size - 1 < n1 ? j + block_size - 1 : n1);
        kj1 = k - j + 1;
        v.resize(kj1);
        v.setZero();

        /* enumerate a shortest vector*/
        v = PotENUM(mu.block(j, j, kj1, kj1), B.segment(j, kj1), logB.segment(j, kj1), kj1);

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

            NTL::LLL(_, cc, 99, 100);

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

#endif // !POT_BKZ_H