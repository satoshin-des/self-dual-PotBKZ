#ifndef SD_POT_BKZ_H_
#define SD_POT_BKZ_H_

#include <iostream>
#include <string>

extern "C" long **BKZ(long **basis, const int block_size, const double reduction_parameter, const int max_loop, const int n, const int m, const int seed);

extern "C" long **PotBKZ(long **basis, const int beta, const double reduction_parameter, const int n, const int m, const int seed);

extern "C" long **DualPotBKZ(long **basis, const int beta, const double delta, const int n, const int m, const int seed);

extern "C" long **SelfDualPotBKZ(long **basis, const int beta, const double reduction_parameter, const int n, const int m, const int seed);

#endif // !SD_POT_BKZ_H_
