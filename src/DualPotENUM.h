#ifndef DUAL_POT_ENUM
#define DUAL_POT_ENUM

#include "Lattice.h"

inline VectorXli Lattice::DualPotENUM(const MatrixXld _gso_coeff_mat, const VectorXld _squared_norm_of_gso_vec, const VectorXld _log_squared_norm_of_gso_vec, const int n)
{
    int i, r[n + 1];
    double tmp, enumeration_upper_bound = _log_squared_norm_of_gso_vec.coeff(0), P = 0;
    VectorXli weight(n), v(n);
    weight.setZero();
    v.setZero();
    Eigen::VectorXd center(n), D(n + 1);
    center.setZero();
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
        tmp = static_cast<double>(v.coeff(k)) - center.coeff(k);
        tmp *= tmp;
        D.coeffRef(k) = D.coeff(k + 1) + tmp * _squared_norm_of_gso_vec.coeff(k);

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

        enumeration_upper_bound = _log_squared_norm_of_gso_vec.head(last_nonzero + 1).array().sum();
        if (LOG099 + P > enumeration_upper_bound)
        {
            if (k == 0)
            {
                return v;
            }
            else
            {
                P += log(D.coeff(k));
                --k;
                r[k] = (r[k] > r[k + 1] ? r[k] : r[k + 1]);
                for (i = r[k]; i > k; --i)
                {
                    sigma.coeffRef(i, k) = sigma.coeff(i + 1, k) + _gso_coeff_mat.coeff(i, k) * v(i);
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
                }
                else
                {
                    v.coeff(k) > center.coeff(k) ? v.coeffRef(k) -= weight.coeff(k) : v.coeffRef(k) += weight.coeff(k);
                    ++weight.coeffRef(k);
                    P -= log(D.coeff(k));
                }
            }
        }
    }
}

#endif // !DUAL_POT_ENUM