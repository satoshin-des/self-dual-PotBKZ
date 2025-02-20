#ifndef POT_LLL_H
#define POT_LLL_H

#include "Lattice.h"

inline void Lattice::PotLLL_(const double reduction_parameter, const int n, const int m)
{
    long double P, P_min, S;
    VectorXli t;
    VectorXld B(n);
    NTL::mat_ZZ c;
    c.SetDims(n, m);
    MatrixXld mu(n, n);

    // LLL基底簡約
    for (int i = 0, j; i < n; ++i)
    {
        for (j = 0; j < m; ++j)
        {
            c[i][j] = basis.coeff(i, j);
        }
    }

    NTL::LLL(_, c, 99, 100);

    for (int i = 0, j; i < n; ++i)
    {
        for (j = 0; j < m; ++j)
        {
            basis.coeffRef(i, j) = NTL::to_long(c[i][j]);
        }
    }

    GSO(B, mu, n, m);

    for (int l = 0, j, i, k, q; l < n;)
    {
        // 部分サイズ基底簡約
        for (j = l - 1; j > -1; --j)
        {
            if (mu.coeff(l, j) > 0.5 || mu.coeff(l, j) < -0.5)
            {
                q = round(mu.coeff(l, j));
                basis.row(l) -= q * basis.row(j);
                mu.row(l).head(j + 1) -= (long double)q * mu.row(j).head(j + 1);
            }
        }

        P = P_min = 1.0;
        k = 0;
        for (j = l - 1; j >= 0; --j)
        {
            S = (mu.row(l).segment(j, l - j).array().square() * B.segment(j, l - j).array()).sum();
            P *= (B.coeff(l) + S) / B.coeff(j);
            if (P < P_min)
            {
                k = j;
                P_min = P;
            }
        }

        if (reduction_parameter > P_min)
        {
            // deep insertion
            t = basis.row(l);
            for (j = l; j > k; --j)
            {
                basis.row(j) = basis.row(j - 1);
            }
            basis.row(k) = t;

            updateDeepInsGSO(k, l, B, mu, n);
            // GSO(B, mu, n, m);
            l = k;
        }
        else
        {
            ++l;
        }
    }
}

#endif // !POT_LLL_H
