#include <sstream>

#include <NTL/LLL.h>

#include <tools.h>

#include "Lattice.h"

extern "C" long **generator(long **basis, const int n, const int seed)
{
    long bit = 10;

    NTL::vec_ZZ v;
    generate_random_HNF(v, n, bit, NTL::to_ZZ(seed));
    NTL::mat_ZZ B;
    B.SetDims(n, n);
    NTL::clear(B);
    B(1, 1) = v(1);
    for (int i = 2; i <= n; i++)
    {
        B(i, 1) = v(i);
        B(i, i) = 1;
    }
    NTL::LLL_XD(B);
    for (int i = 0, j; i < n; ++i)
    {
        for (j = 0; j < n; ++j)
        {
            basis[i][j] = NTL::to_long(B(i + 1, j + 1));
        }
    }

    return basis;
}
