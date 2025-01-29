#ifndef DUAL_POT_ENUM
#define DUAL_POT_ENUM

#include "Lattice.h"

inline VectorXli Lattice::DualPotENUM(const MatrixXld mu, const VectorXld B, const VectorXld logB, const int n)
{
    int i, r[n + 1];
    int last_nonzero = 0;      // index of last non-zero elements
    double potential = 0;      // potential-like information of lattice
    double R = logB.coeff(0);  // upper bounds
    VectorXli weight(n);       // weight of zigzag searching
    Eigen::VectorXd center(n); // center of zegzag searching
    VectorXli v(n);            // coefficient vector to output
    Eigen::VectorXd D(n + 1);
    Eigen::MatrixXd sigma(n + 1, n);
    weight.setZero();
    v.setZero();
    center.setZero();
    D.setZero();
    sigma.setZero();

    v.coeffRef(0) = 1;

    for (i = 0; i <= n; ++i)
    {
        r[i] = i;
    }

    for (int k = 0;;)
    {
        m_temp = static_cast<double>(v.coeff(k)) - center.coeff(k);
        m_temp *= m_temp;
        D.coeffRef(k) = D.coeff(k + 1) + m_temp * B.coeff(k);

        if (last_nonzero == 0)
        {
            potential = log(D.coeff(0));
        }
        else if (last_nonzero == 1)
        {
            potential = log(D.coeff(0));
            potential += potential;
        }
        else
        {
            potential = 0;
            for (i = 0; i <= last_nonzero - 2; ++i)
            {
                potential += log(D.coeff(i));
            }
            m_temp = log(D.coeff(last_nonzero - 1));
            m_temp += m_temp;
            potential += m_temp;
        }

        R = logB.head(last_nonzero + 1).array().sum();
        if (LOG099 + potential > R)
        {
            if (k == 0)
            {
                return v;
            }
            else
            {
                potential += log(D.coeff(k));
                --k;
                r[k] = (r[k] > r[k + 1] ? r[k] : r[k + 1]);
                for (i = r[k]; i > k; --i)
                {
                    sigma.coeffRef(i, k) = sigma.coeff(i + 1, k) + mu.coeff(i, k) * v(i);
                }
                center.coeffRef(k) = -sigma.coeff(k + 1, k);
                v.coeffRef(k) = round(center.coeff(k));
                weight.coeffRef(k) = 1;
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
#if 1
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
                }
                else
                {
                    v.coeff(k) > center.coeff(k) ? v.coeffRef(k) -= weight.coeff(k) : v.coeffRef(k) += weight.coeff(k);
                    ++weight.coeffRef(k);
                    potential -= log(D.coeff(k));
                }
            }
        }
    }
}

#endif // !DUAL_POT_ENUM