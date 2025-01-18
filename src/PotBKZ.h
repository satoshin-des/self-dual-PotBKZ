#ifndef POT_BKZ_H
#define POT_BKZ_H

#include "Lattice.h"

inline void Lattice::POT_BKZ(const int beta, const double reduction_parameter)
{
    const int n1 = _n - 1, n2 = _n - 2;
    VectorXli v, w;
    MatrixXli tmp_b(_n, _n);
    NTL::mat_ZZ cc;

    POT_LLL(reduction_parameter);
    GSO();

    for (int z = 0, j = 0, i, k, l, kj1; z < _n - 1;)
    {
        if (j == n2)
        {
            j = 0;
        }
        ++j;
        k = (j + beta - 1 < n1 ? j + beta - 1 : n1);
        kj1 = k - j + 1;
        v.resize(kj1);
        v.setZero();

        /* enumerate a shortest vector*/
        v = PotENUM(_gso_coeff_mat.block(j, j, kj1, kj1), _squared_norm_of_gso_vec.segment(j, kj1), _log_squared_norm_of_gso_vec.segment(j, kj1), kj1);

        if (!v.isZero())
        {
            z = 0;

            w = v * basis.block(j, 0, kj1, _m);
            cc.SetDims(_n + 1, _m);
            for (l = 0; l < _m; ++l)
            {
                for (i = 0; i < j; ++i)
                    cc[i][l] = basis.coeffRef(i, l);
                cc[j][l] = w[l];
                for (i = j + 1; i < _n + 1; ++i)
                    cc[i][l] = basis.coeffRef(i - 1, l);
            }
            NTL::LLL(_, cc, 99, 100);

            for (i = 0; i < _n; ++i)
                for (l = 0; l < _m; ++l)
                    basis.coeffRef(i, l) = NTL::to_long(cc[i + 1][l]);

            POT_LLL(reduction_parameter);
            GSO();
        }
        else
        {
            ++z;
        }
    }
}

#endif // !POT_BKZ_H