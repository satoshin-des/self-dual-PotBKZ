#ifndef SELF_DUAL_POT_BKZ_H
#define SELF_DUAL_POT_BKZ_H

#include "Lattice.h"
#include "PotENUM.h"
#include "DualPotENUM.h"
#include "PotLLL.h"
#include "DualPotLLL.h"

inline void Lattice::SelfDualPotBKZ_each_(const int beta, const double d, const int n, const int m, FILE *fp)
{
    const int n1 = n - 1, n2 = n - 2;
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

    for (int zp = 0, jp = 0, i, k, l, kj1, zd = n, jd = n, IsPrimal = 1; zp < n - 1 && zd > 1;)
    {
        /// ================================
        /// Primal part
        /// ================================
        // printf("Primal = %3d, Dual = %3d\n", zp, zd);

        if (IsPrimal)
        {
            if (jp == n2)
            {
                jp = 0;
                ++SelfPotTour;
                IsPrimal = 0;
            }
            ++jp;
            k = (jp + beta - 1 < n1 ? jp + beta - 1 : n1);
            kj1 = k - jp + 1;

            fprintf(fp, "%Lf\n", logPot(B, n));

            v.resize(kj1);
            v.setZero();

            /* enumerate a shortest vector*/
            v = PotENUM(mu.block(jp, jp, kj1, kj1), B.segment(jp, kj1), logB.segment(jp, kj1), kj1);

            if (!v.isZero())
            {
                zp = 0;

                w = v * basis.block(jp, 0, kj1, m);
                cc.SetDims(n + 1, m);
                for (l = 0; l < m; ++l)
                {
                    for (i = 0; i < jp; ++i)
                        cc[i][l] = basis.coeffRef(i, l);
                    cc[jp][l] = w[l];
                    for (i = jp + 1; i < n + 1; ++i)
                        cc[i][l] = basis.coeffRef(i - 1, l);
                }
                NTL::LLL(_, cc, 99, 100);

                for (i = 0; i < n; ++i)
                    for (l = 0; l < m; ++l)
                        basis.coeffRef(i, l) = NTL::to_long(cc[i + 1][l]);

                DualPotLLL_(d, n, m);
                GSO(B, logB, mu, n, m);
            }
            else
                ++zp;
        }

        /// ================================
        /// Dual part
        /// ================================
        if (!IsPrimal)
        {
            if (jd == 1)
            {
                ++SelfPotTour;
                jd = n;
                IsPrimal = 1;
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
                --zd;
            }
            else
            {
                zd = n;

                tmp_b = Insert(basis.block(k, 0, kj1, m), v, kj1, m);
                basis.block(k, 0, kj1, m) = tmp_b.block(0, 0, kj1, m);

                PotLLL_(d, n, m);
                GSO(B, logB, mu, n, m);
            }
        }
    }
}

#endif // !SELF_DUAL_POT_BKZ_H