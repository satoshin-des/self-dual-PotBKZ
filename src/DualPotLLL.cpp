#include "Lattice.h"

#include <iostream>
#include <cmath>

#include <eigen3/Eigen/Dense>

#include <NTL/ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>

void Lattice::DualPotLLL_(const double reduction_parameter, const int n, const int m)
{
    long double potential;
    long double minimal_potential;
    long double s;
    long double D;
    MatrixXld mu = MatrixXld::Zero(n, n), nu = MatrixXld::Zero(n, n);
    VectorXld dual_D(n);
    VectorXld B = VectorXld::Zero(n);
    NTL::mat_ZZ c;
    c.SetDims(n, m);

    for (int i = 0, j; i < n; ++i)
    {
        for (j = 0; j < m; ++j)
        {
            c[i][j] = basis.coeff(i, j);
        }
    }
    NTL::LLL_FP(c, 0.99);
    for (int i = 0, j; i < n; ++i)
    {
        for (j = 0; j < m; ++j)
        {
            basis.coeffRef(i, j) = NTL::to_long(c[i][j]);
        }
    }

    GSO(B, mu, n, m);

    for (int k = n - 1, j, i, l, q, h; k >= 0;)
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
            D = 1.0 / B.coeff(k);
            dual_D.setZero();
            dual_D.coeffRef(k) = D;
            for (h = k + 1; h < n; ++h)
            {
                D += nu.coeff(k, h) * nu.coeff(k, h) / B.coeff(h);
                dual_D.coeffRef(h) = D;
            }

            dualDeepInsertion(m, k, l);
            updateDualDeepInsGSO(k, l, B, mu, nu, dual_D, n);

            k = l;
        }
        else
        {
            --k;
        }
    }
}
