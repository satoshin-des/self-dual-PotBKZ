#include "Lattice.h"

inline void Lattice::POT_LLL(const double reduction_parameter)
{
    long double P, P_min, S;
    VectorXli t;
    NTL::mat_ZZ c;
    c.SetDims(_n, _m);

    // LLL基底簡約
    for (int i = 0, j; i < _n; ++i)
    {
        for (j = 0; j < _m; ++j)
        {
            c[i][j] = basis.coeff(i, j);
        }
    }
    NTL::LLL(_, c, 99, 100);
    for (int i = 0, j; i < _n; ++i)
    {
        for (j = 0; j < _m; ++j)
        {
            basis.coeffRef(i, j) = NTL::to_long(c[i][j]);
        }
    }

    GSO();

    for (int l = 0, l1, j, i, k, q; l < _n;)
    {
        l1 = l - 1;
        // 部分サイズ基底簡約
        for (j = l1; j > -1; --j)
            if (_gso_coeff_mat.coeff(l, j) > 0.5 || _gso_coeff_mat.coeff(l, j) < -0.5)
            {
                q = round(_gso_coeff_mat.coeff(l, j));
                basis.row(l) -= q * basis.row(j);
                _gso_coeff_mat.row(l).head(j + 1) -= static_cast<long double>(q) * _gso_coeff_mat.row(j).head(j + 1);
            }

        P = P_min = 1.0;
        k = 0;
        for (j = l1; j >= 0; --j)
        {
            S = (_gso_coeff_mat.row(l).segment(j, l - j).array().square() * _squared_norm_of_gso_vec.segment(j, l - j).array()).sum();
            P *= (_squared_norm_of_gso_vec.coeff(l) + S) / _squared_norm_of_gso_vec.coeff(j);
            if (P < P_min)
            {
                k = j;
                P_min = P;
            }
        }

        if (reduction_parameter > P_min)
        {
            // Deep insertion
            t = basis.row(l);
            for (j = l; j > k; --j)
            {
                basis.row(j) = basis.row(j - 1);
            }
            basis.row(k) = t;

            GSO();
            l = k;
        }
        else
        {
            ++l;
        }
    }
}