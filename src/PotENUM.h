#ifndef POT_ENUM_H
#define POT_ENUM_H

#include <iostream>
#include <random>
#include <algorithm>
#include <vector>
#include <cmath>
#include <eigen3/Eigen/Dense>

#include "Lattice.h"

/// @brief Enumerates a short and small potential lattice vector
/// @param mu GSO-coefficient matrix.
/// @param B Squared norms of GSO vectors.
/// @param logB Logarithm valude of squared norms of GSO-vectors.
/// @param n Rank of lattice.
/// @return VectorXli
inline VectorXli Lattice::PotENUM(const MatrixXld mu, const VectorXld B, const VectorXld logB, const int n)
{
    int i, r[n + 1];
    double R = logB.coeff(0), P = 0;
    VectorXli w(n), v(n);
    w.setZero();
    v.setZero();
    Eigen::VectorXd c(n), D(n + 1);
    c.setZero();
    D.setZero();
    Eigen::MatrixXd sigma(n + 1, n);
    sigma.setZero();

    v.coeffRef(0) = 1;

    for (i = 0; i <= n; ++i)
    {
        r[i] = i;
    }

    for (int k = 0, last_nonzero = 0;;)
    {
        m_temp = static_cast<double>(v.coeff(k)) - c.coeff(k);
        m_temp *= m_temp;
        D.coeffRef(k) = D.coeff(k + 1) + m_temp * B.coeff(k);

        if ((k + 1) * log(D.coeff(k)) + P < (k + 1) * LOG099 + R)
        {
            if (k == 0)
                return v;
            else
            {
                P += log(D.coeff(k));
                --k;
                r[k] = (r[k] > r[k + 1] ? r[k] : r[k + 1]);
                for (i = r[k]; i > k; --i)
                    sigma.coeffRef(i, k) = sigma.coeff(i + 1, k) + mu.coeff(i, k) * v(i);
                c.coeffRef(k) = -sigma.coeff(k + 1, k);
                v.coeffRef(k) = round(c.coeff(k));
                w.coeffRef(k) = 1;
            }
        }
        else
        {
            ++k;
            if (k == n)
            {
                v.setZero();
                return v;
            }
            else
            {
                r[k - 1] = k;
                if (k >= last_nonzero)
                {
                    last_nonzero = k;
                    ++v.coeffRef(k);
#if 0
                    if (v.coeff(last_nonzero) >= 2)
                    {
                        ++k;
                        if (k == n)
                        {
                            v.setZero();
                            return v;
                        }
                        else
                        {
                            r[k - 1] = k;
                            last_nonzero = k;
                            v.coeffRef(last_nonzero) = 1;
                        }
                    }
#endif
                    P = 0;
                    R = logB.head(last_nonzero + 1).array().sum();
                }
                else
                {
                    v.coeff(k) > c.coeff(k) ? v.coeffRef(k) -= w.coeff(k) : v.coeffRef(k) += w.coeff(k);
                    ++w.coeffRef(k);
                    P -= log(D.coeff(k));
                }
            }
        }
    }
}

#endif // !POT_ENUM_H