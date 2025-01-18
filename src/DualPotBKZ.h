#ifndef DUAL_POT_BKZ_H
#define DUAL_POT_BKZ_H

#include "Lattice.h"

inline void Lattice::DUAL_POT_BKZ(const int beta, const double delta)
{
    VectorXli x;
    MatrixXli tmp_b;
    _log_squared_norm_of_gso_vec.setZero();
    MatrixXld BB;

    DUAL_POT_LLL(0.99);

    GSO();

    for (int z = _n, i, j = _n, k, dim_of_block_lattice; z > 2;)
    {
        if (j == 1)
        {
            j = _n;
        }
        --j;
        k = std::max(j - beta + 1, 0);
        dim_of_block_lattice = j - k + 1;

        _dual_squared_norm_of_gso_vec.resize(dim_of_block_lattice);
        _dual_log_squared_norm_of_gso_vec.resize(dim_of_block_lattice);
        _dual_gso_coeff_mat.resize(dim_of_block_lattice, dim_of_block_lattice);
        DualGSO(_squared_norm_of_gso_vec.segment(k, dim_of_block_lattice), _log_squared_norm_of_gso_vec.segment(k, dim_of_block_lattice), _gso_coeff_mat.block(k, k, dim_of_block_lattice, dim_of_block_lattice), _dual_squared_norm_of_gso_vec, _dual_log_squared_norm_of_gso_vec, _dual_gso_coeff_mat, dim_of_block_lattice, dim_of_block_lattice);

        // Dual Enumeration
        x = DualPotENUM(_dual_gso_coeff_mat, _dual_squared_norm_of_gso_vec, _dual_log_squared_norm_of_gso_vec, dim_of_block_lattice);

        if (x.isZero())
        {
            --z;
        }
        else
        {
            z = _n;

            tmp_b = Insert(basis.block(k, 0, dim_of_block_lattice, _m), x, dim_of_block_lattice, _m);
            basis.block(k, 0, dim_of_block_lattice, _m) = tmp_b.block(0, 0, dim_of_block_lattice, _m);

            POT_LLL(delta);
            GSO();
        }
    }
}

#endif // !DUAL_POT_BKZ_H