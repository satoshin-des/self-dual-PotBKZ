#ifndef DUAL_POT_ENUM
#define DUAL_POT_ENUM

#include "Lattice.h"

inline VectorXli Lattice::DualPotENUM(const MatrixXld mu, const VectorXld B, const VectorXld logB, const int n)
{
    int i, r[n + 1];
    double tmp, R = logB.coeff(0), P = 0;
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
        r[i] = i;
    for (int k = 0, last_nonzero = 0;;)
    {
        tmp = (double)v.coeff(k) - c.coeff(k);
        tmp *= tmp;
        D.coeffRef(k) = D.coeff(k + 1) + tmp * B.coeff(k);

        if (last_nonzero == 0)
        {
            P = log(D.coeff(0));
        }
        else if (last_nonzero == 1)
        {
            P = log(D.coeff(0));
            P += P;
        }
        else
        {
            P = 0;
            for (i = 0; i <= last_nonzero - 2; ++i)
            {
                P += log(D.coeff(i));
            }
            tmp = log(D.coeff(last_nonzero - 1));
            tmp += tmp;
            P += tmp;
        }

        R = logB.head(last_nonzero + 1).array().sum();
        if (LOG099 + P > R)
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

#endif // !DUAL_POT_ENUM