#ifndef SELF_DUAL_POT_BKZ_H
#define SELF_DUAL_POT_BKZ_H

#include "Lattice.h"
#include "PotENUM.h"
#include "DualPotENUM.h"
#include "PotLLL.h"
#include "DualPotLLL.h"

inline void Lattice::SelfDualPotBKZ_(const int beta, const double d, const int n, const int m, FILE *fp)
{
    VectorXli v, w;
    MatrixXli tmp_b(n, n);
    VectorXld B(n), logB(n);
    MatrixXld mu(n, n);
    NTL::mat_ZZ cc;
    VectorXld C, logC;
    MatrixXld hmu, BB;

    GSO(B, logB, mu, n, m);
    fprintf(fp, "%Lf\n", logPot(B, n));

    DualPotLLL_(0.99, n, m);

    for (int primal_z = 0, jp = 0, i, k, l, kj1, dual_z = n, jd = n, is_primal = 1; primal_z < n - 1 && dual_z > 1;)
    {
        /// ================================
        /// Primal part
        /// ================================
        if (is_primal)
        {
            if (jp == n - 2)
            {
                jp = 0;
                is_primal = 0;
            }
            ++jp;
            k = (jp + beta - 1 < n - 1 ? jp + beta - 1 : n - 1);
            kj1 = k - jp + 1;

            fprintf(fp, "%Lf\n", logPot(B, n));

            v.resize(kj1);
            v.setZero();

            /* enumerate a shortest vector*/
            v = PotENUM(mu.block(jp, jp, kj1, kj1), B.segment(jp, kj1), logB.segment(jp, kj1), kj1);

            if (!v.isZero())
            {
                primal_z = 0;

                w = v * basis.block(jp, 0, kj1, m);
                cc.SetDims(n + 1, m);
                for (l = 0; l < m; ++l)
                {
                    for (i = 0; i < jp; ++i)
                    {
                        cc[i][l] = basis.coeffRef(i, l);
                    }
                    cc[jp][l] = w[l];
                    for (i = jp + 1; i < n + 1; ++i)
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

                DualPotLLL_(d, n, m);
                GSO(B, logB, mu, n, m);
            }
            else
            {
                ++primal_z;
            }
        }

        /// ================================
        /// Dual part
        /// ================================
        if (!is_primal)
        {
            if (jd == 1)
            {
                jd = n;
                is_primal = 1;
            }
            --jd;
            k = (jd - beta + 1 > 0 ? jd - beta + 1 : 0);
            kj1 = jd - k + 1;

            fprintf(fp, "%Lf\n", logPot(B, n));

            C.resize(kj1);
            logC.resize(kj1);
            hmu.resize(kj1, kj1);
            DualGSO(B.segment(k, kj1), logB.segment(k, kj1), mu.block(k, k, kj1, kj1), C, logC, hmu, kj1, kj1);

            // Dual Enumeration
            v = DualPotENUM(hmu, C, logC, kj1);

            if (v.isZero())
            {
                --dual_z;
            }
            else
            {
                dual_z = n;

                tmp_b = Insert(basis.block(k, 0, kj1, m), v, kj1, m);
                basis.block(k, 0, kj1, m) = tmp_b.block(0, 0, kj1, m);

                PotLLL_(d, n, m);
                GSO(B, logB, mu, n, m);
            }
        }
    }
}

#endif // !SELF_DUAL_POT_BKZ_H