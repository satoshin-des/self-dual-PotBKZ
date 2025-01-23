#ifndef DUAL_POT_LLL_H
#define DUAL_POT_LLL_H

#include "Lattice.h"

inline void Lattice::DualPotLLL_(const double reduction_parameter, const int n, const int m)
{
    double potential;
    double minimal_potential;
    double s;
    MatrixXld mu(n, n), nu(n, n);
    mu.setZero();
    nu.setZero();
    VectorXld B(n);
    B.setZero();
    NTL::mat_ZZ c;
    c.SetDims(n, m);

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

    for (int k = n - 1, j, i, l, q; k >= 0;)
    {
        nu.coeffRef(k, k) = 1.0;

        // Dual size reduction
        for (j = k + 1; j < n; ++j)
        {
            nu.coeffRef(k, j) = 0;
            for (i = k; i < j; ++i)
            {
                nu.coeffRef(k, j) -= mu.coeff(j, i) * nu.coeff(k, i);
            }

            if (nu.coeff(k, j) > 0.5 || nu.coeff(k, j) < -0.5)
            {
                q = round(nu.coeff(k, j));
                basis.row(j) += q * basis.row(k);
                nu.row(k).tail(n - j + 1) -= q * nu.row(j).tail(n - j + 1);
                mu.row(j).head(k + 1) += q * mu.row(k).head(k + 1);
            }
        }

        potential = 1.0;
        minimal_potential = 1.0;
        l = n - 1;
        for (j = k + 1; j < n; ++j)
        {
            s = 0.0;
            for (i = k; i <= j; ++i)
            {
                s += nu.coeff(k, i) * nu.coeff(k, i) / B.coeff(i);
            }
            potential *= B.coeff(j);
            potential *= s;

            if (potential < minimal_potential)
            {
                l = j;
                minimal_potential = potential;
            }
        }

        if (reduction_parameter > minimal_potential)
        {
            dualDeepInsertion(m, k, l);
            GSO(B, mu, n, m);
            k = l;
        }
        else
        {
            --k;
        }
    }
}

#endif // !DUAL_POT_LLL_H