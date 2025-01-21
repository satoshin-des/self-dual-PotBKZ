#ifndef DUAL_POT_BKZ_H
#define DUAL_POT_BKZ_H

#include "Lattice.h"
#include "DualPotENUM.h"
#include "DualPotLLL.h"
#include "PotLLL.h"

inline void Lattice::DualPotBKZ_(const int beta, const double delta, const int n, const int m, FILE *fp)
{
    VectorXli x;
    MatrixXli tmp_b;
    VectorXld B(n), logB(n), C, logC;
    B.setZero();
    logB.setZero();
    MatrixXld mu(n, n), hmu, BB;
    mu.setZero();

    GSO(B, logB, mu, n, m);
    fprintf(fp, "%Lf\n", logPot(B, n));

    // DualPotLLL(basis, 0.99, n, m);
    for (int z = n, i, j = n, k, d; z > 1;)
    {
        if (j == 1)
        {
            ++PotTour;
            j = n;
        }
        --j;
        k = (j - beta + 1 > 0 ? j - beta + 1 : 0);
        d = j - k + 1;

        printf("z = %d\n", z);
        ++Tr;
        fprintf(fp, "%Lf\n", logPot(B, n));

        C.resize(d);
        logC.resize(d);
        hmu.resize(d, d);
        DualGSO(B.segment(k, d), logB.segment(k, d), mu.block(k, k, d, d), C, logC, hmu, d, d);

        // Dual Enumeration
        x = DualPotENUM(hmu, C, logC, d);

        if (x.isZero())
        {
            --z;
        }
        else
        {
            z = n;

            tmp_b = Insert(basis.block(k, 0, d, m), x, d, m);
            basis.block(k, 0, d, m) = tmp_b.block(0, 0, d, m);

            DualPotLLL_(delta, n, m);
            // PotLLL_(delta, n, m);
            GSO(B, logB, mu, n, m);
        }
    }
    fprintf(_TOUR_DTA_, "%d\n", Tr);
    printf("Dual: %d\n", Tr);
}

#endif // !DUAL_POT_BKZ_H