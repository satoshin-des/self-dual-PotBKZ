#ifndef SELF_DUAL_POT_BKZ_H
#define SELF_DUAL_POT_BKZ_H

#include "Lattice.h"

inline void Lattice::SELF_DUAL_POT_BKZ(const int beta, const double reduction_parameter)
{
    const int n1 = _n - 1, n2 = _n - 2;
    VectorXli v, w;
    MatrixXli tmp_b(_n, _n);
    NTL::mat_ZZ cc;

    GSO();

    DUAL_POT_LLL(0.99);

    for (int zp = 0, jp = 0, i, k, l, kj1, zd = _n, jd = _n, is_primal_part = 1; zp < _n - 1 && zd > 1;)
    {
        /// ================================
        /// Primal part
        /// ================================
        if (is_primal_part)
        {
            if (jp == n2)
            {
                jp = 0;
                is_primal_part = 0;
            }
            ++jp;
            k = (jp + beta - 1 < n1 ? jp + beta - 1 : n1);
            kj1 = k - jp + 1;
            v.resize(kj1);
            v.setZero();
            
            /* enumerate a shortest vector*/
            v = PotENUM(_gso_coeff_mat.block(jp, jp, kj1, kj1), _squared_norm_of_gso_vec.segment(jp, kj1), _log_squared_norm_of_gso_vec.segment(jp, kj1), kj1);

            if (!v.isZero())
            {
                zp = 0;

                w = v * basis.block(jp, 0, kj1, _m);
                cc.SetDims(_n + 1, _m);

                for (l = 0; l < _m; ++l)
                {
                    for (i = 0; i < jp; ++i)
                    {
                        cc[i][l] = basis.coeffRef(i, l);
                    }
                    cc[jp][l] = w[l];
                    for (i = jp + 1; i < _n + 1; ++i)
                    {
                        cc[i][l] = basis.coeffRef(i - 1, l);
                    }
                }

                NTL::LLL(_, cc, 99, 100);

                for (i = 0; i < _n; ++i)
                {
                    for (l = 0; l < _m; ++l)
                    {
                        basis.coeffRef(i, l) = NTL::to_long(cc[i + 1][l]);
                    }
                }

                DUAL_POT_LLL(reduction_parameter);

                GSO();
            }
            else
            {
                ++zp;
            }
        }

        /// ================================
        /// Dual part
        /// ================================
        if ((!is_primal_part))
        {
            if (jd == 1)
            {
                jd = _n;
                is_primal_part = 1;
            }
            --jd;
            k = (jd - beta + 1 > 0 ? jd - beta + 1 : 0);
            kj1 = jd - k + 1;

            _dual_squared_norm_of_gso_vec.resize(kj1);
            _dual_log_squared_norm_of_gso_vec.resize(kj1);
            _dual_gso_coeff_mat.resize(kj1, kj1);
            DualGSO(_squared_norm_of_gso_vec.segment(k, kj1),
                    _log_squared_norm_of_gso_vec.segment(k, kj1),
                    _gso_coeff_mat.block(k, k, kj1, kj1),
                    _dual_squared_norm_of_gso_vec,
                    _dual_log_squared_norm_of_gso_vec,
                    _dual_gso_coeff_mat,
                    kj1, kj1);

            // Dual Enumeration
            v = DualPotENUM(_dual_gso_coeff_mat, _dual_squared_norm_of_gso_vec, _dual_log_squared_norm_of_gso_vec, kj1);

            if (v.isZero())
            {
                --zd;
            }
            else
            {
                zd = _n;

                tmp_b = Insert(basis.block(k, 0, kj1, _m), v, kj1, _m);
                basis.block(k, 0, kj1, _m) = tmp_b.block(0, 0, kj1, _m);

                POT_LLL(reduction_parameter);
                GSO();
            }
        }
    }
}

#endif // !SELF_DUAL_POT_BKZ_H