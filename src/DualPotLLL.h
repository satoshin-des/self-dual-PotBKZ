#ifndef DUAL_POT_LLL_H
#define DUAL_POT_LLL_H

#include "Lattice.h"

inline void Lattice::DUAL_POT_LLL(const double reduction_parameter)
{
    double P, P_max, P_min, s;
    NTL::mat_ZZ c;
    c.SetDims(_n, _m);

    _dual_gso_coeff_mat.resize(_n, _n);
    _dual_gso_coeff_mat.setZero();
    _dual_squared_norm_of_gso_vec.resize(_n);
    _dual_squared_norm_of_gso_vec.setZero();
    _gso_coeff_mat.setZero();
    _gso_vec_mat.setZero();
    _squared_norm_of_gso_vec.setZero();

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

    for (int k = _n - 1, j, i, l, q; k >= 0;)
    {
        _dual_gso_coeff_mat.coeffRef(k, k) = 1.0;

        // Dual size reduction
        for (j = k + 1; j < _n; ++j)
        {
            _dual_gso_coeff_mat.coeffRef(k, j) = 0;
            for (i = k; i < j; ++i)
            {
                _dual_gso_coeff_mat.coeffRef(k, j) -= _gso_coeff_mat.coeff(j, i) * _dual_gso_coeff_mat.coeff(k, i);
            }

            if (_dual_gso_coeff_mat.coeff(k, j) > 0.5 || _dual_gso_coeff_mat.coeff(k, j) < -0.5)
            {
                q = round(_dual_gso_coeff_mat.coeff(k, j));
                basis.row(j) += q * basis.row(k);
                _dual_gso_coeff_mat.row(k).tail(_n - j + 1) -= static_cast<long double>(q) * _dual_gso_coeff_mat.row(j).tail(_n - j + 1);
                _gso_coeff_mat.row(j).head(k + 1) += q * _gso_coeff_mat.row(k).head(k + 1);
            }
        }

        P = P_min = 1.0;
        l = _n - 1;
        for (j = k + 1; j < _n; ++j)
        {
            s = 0.0;
            for (i = k; i <= j; ++i)
            {
                s += _dual_gso_coeff_mat.coeff(k, i) * _dual_gso_coeff_mat.coeff(k, i) / _squared_norm_of_gso_vec.coeff(i);
            }
            P *= _squared_norm_of_gso_vec.coeff(j);
            P *= s;

            if (P < P_min)
            {
                l = j;
                P_min = P;
            }
        }

        if (reduction_parameter > P_min)
        {
            DualDeepInsertion(k, l);
            GSO();
            k = l;
        }
        else
        {
            --k;
        }
    }
}

#endif // !DUAL_POT_LLL_H