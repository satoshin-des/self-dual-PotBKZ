#include "Lattice.h"

#include <iostream>
#include <cmath>

#include <eigen3/Eigen/Dense>

bool Lattice::ENUM(VectorXli &u, const MatrixXld mu, const VectorXld B, VectorXld &rho, const int n, double R)
{
    bool has_solution = false;
    int i, r[n + 1];
    int last_nonzero = 0;                              // index of last non-zero elements
    VectorXli weight = VectorXli::Zero(n);             // weight of zigzag searching
    VectorXli coeff_vector = VectorXli::Zero(n);       // coefficient vector to putput
    Eigen::VectorXd center = Eigen::VectorXd::Zero(n); // center of zigzag serching
    Eigen::MatrixXd sigma = Eigen::MatrixXd::Zero(n + 1, n);
    u = VectorXli::Zero(n);
    coeff_vector.coeffRef(0) = 1;
    rho.setZero();
    for (i = 0; i < n; ++i)
    {
        r[i] = i;
    }

    for (int k = 0;;)
    {
        m_temp = static_cast<long double>(coeff_vector.coeff(k)) - center.coeff(k);
        m_temp *= m_temp;
        rho.coeffRef(k) = rho.coeff(k + 1) + m_temp * B.coeff(k); // rho[k]=∥πₖ(shortest_vec)∥
        if (rho.coeff(k) <= R)
        {
            if (k == 0)
            {
                has_solution = true;
                u = coeff_vector;
                R = fmin(R, 0.99 * rho.coeff(0));
            }
            else
            {
                --k;
                r[k] = (r[k] > r[k + 1] ? r[k] : r[k + 1]);
                for (i = r[k]; i > k; --i)
                {
                    sigma.coeffRef(i, k) = sigma.coeff(i + 1, k) + mu.coeff(i, k) * coeff_vector.coeff(i);
                }
                center.coeffRef(k) = -sigma.coeff(k + 1, k);
                coeff_vector.coeffRef(k) = round(center.coeff(k));
                weight.coeffRef(k) = 1; // 小さいやつが見つかったら、変分を元に戻す
            }
        }
        else
        {
            ++k;
            if (k == n)
            { // no solution
                return has_solution;
            }
            else
            {
                r[k] = k;
                if (k >= last_nonzero)
                {
                    last_nonzero = k;
                    ++coeff_vector.coeffRef(k);
                }
                else
                {
                    coeff_vector.coeff(k) > center.coeff(k) ? coeff_vector.coeffRef(k) -= weight.coeff(k) : coeff_vector.coeffRef(k) += weight.coeff(k);
                    ++weight.coeffRef(k);
                }
            }
        }
    }
}
