#include "Lattice.h"

#include <iostream>

#include <eigen3/Eigen/Dense>

#include <NTL/RR.h>

void Lattice::setDims(const int n, const int m)
{
    basis.resize(n, m);
}

void Lattice::dualDeepInsertion(const int n, const int k, const int l)
{
    const VectorXli t = basis.row(k);
    for (int j = k; j < l; ++j)
    {
        basis.row(j) = basis.row(j + 1);
    }
    basis.row(l) = t;
}

void Lattice::updateDeepInsGSO(const int i, const int k, VectorXld &B, MatrixXld &mu, const int n)
{
    int j, l;
    long double T, eps;
    VectorXld P = VectorXld::Zero(n), D = VectorXld::Zero(n), S = VectorXld::Zero(n);

    P.coeffRef(k) = D.coeffRef(k) = B.coeff(k);
    for (j = k - 1; j >= i; --j)
    {
        P.coeffRef(j) = mu.coeff(k, j) * B.coeff(j);
        D.coeffRef(j) = D.coeff(j + 1) + mu.coeff(k, j) * P.coeff(j);
    }

    for (j = k; j > i; --j)
    {
        T = mu.coeff(k, j - 1) / D.coeff(j);
        for (l = n - 1; l > k; --l)
        {
            S.coeffRef(l) += mu.coeff(l, j) * P.coeff(j);
            mu.coeffRef(l, j) = mu.coeff(l, j - 1) - T * S.coeff(l);
        }
        for (l = k; l > j; --l)
        {
            S.coeffRef(l) += mu.coeff(l - 1, j) * P.coeff(j);
            mu.coeffRef(l, j) = mu.coeff(l - 1, j - 1) - T * S.coeff(l);
        }
    }

    T = 1.0 / D.coeff(i);

    for (l = n - 1; l > k; --l)
    {
        mu.coeffRef(l, i) = T * (S.coeff(l) + mu.coeff(l, i) * P.coeff(i));
    }
    for (l = k; l >= i + 2; --l)
    {
        mu.coeffRef(l, i) = T * (S.coeff(l) + mu.coeff(l - 1, i) * P.coeff(i));
    }

    mu.coeffRef(i + 1, i) = T * P.coeff(i);
    for (j = 0; j < i; ++j)
    {
        eps = mu.coeff(k, j);
        for (l = k; l > i; --l)
        {
            mu.coeffRef(l, j) = mu.coeff(l - 1, j);
        }
        mu.coeffRef(i, j) = eps;
    }

    for (j = k; j > i; --j)
    {
        B.coeffRef(j) = D.coeff(j) * B.coeff(j - 1) / D.coeff(j - 1);
    }
    B.coeffRef(i) = D.coeff(i);
}

void Lattice::updateDualDeepInsGSO(const int k, const int l, VectorXld &B, MatrixXld &mu, MatrixXld &dual_mu, const VectorXld dual_D, const int n)
{
    int i, j, h;
    long double sum;
    MatrixXld xi = mu;

    for (i = l + 1; i < n; ++i)
    {
        sum = 0;
        for (h = k; h <= l; ++h)
        {
            sum += dual_mu.coeff(k, h) * mu.coeff(i, h);
        }
        xi.coeffRef(i, l) = sum;
    }

    for (j = k; j < l; ++j)
    {
        for (i = j + 1; i < l; ++i)
        {
            sum = 0;
            for (h = k; h <= j; ++h)
            {
                sum += dual_mu.coeff(k, h) * mu.coeff(i + 1, h);
            }

            xi.coeffRef(i, j) = mu.coeff(i + 1, j + 1) * dual_D.coeff(j) / dual_D.coeff(j + 1) - dual_mu.coeff(k, j + 1) / (dual_D.coeff(j + 1) * B.coeff(j + 1)) * sum;
        }
        xi.coeffRef(l, j) = -dual_mu.coeff(k, j + 1) / (dual_D.coeff(j + 1) * B.coeff(j + 1));
        for (i = l + 1; i < n; ++i)
        {
            sum = 0;
            for (h = k; h <= j; ++h)
            {
                sum += dual_mu.coeff(k, h) * mu.coeff(i, h);
            }

            xi.coeffRef(i, j) = mu.coeff(i, j + 1) * dual_D.coeff(j) / dual_D.coeff(j + 1) - dual_mu.coeff(k, j + 1) / (dual_D.coeff(j + 1) * B.coeff(j + 1)) * sum;
        }
    }

    for (j = 0; j < k; ++j)
    {
        for (i = k; i < l; ++i)
        {
            xi.coeffRef(i, j) = mu.coeff(i + 1, j);
        }
        xi.coeffRef(l, j) = mu.coeff(k, j);
    }

    mu = xi;
    for (j = k; j < l; ++j)
    {
        B.coeffRef(j) = dual_D.coeff(j + 1) * B.coeff(j + 1) / dual_D.coeff(j);
    }
    B.coeffRef(l) = 1 / dual_D.coeff(l);
}

long double Lattice::logPot(const VectorXld B, const int n)
{
    NTL::RR logp;
    logp = 0;
    for (int i = 0; i < n; ++i)
    {
        logp += (n - i) * NTL::log(NTL::to_RR((double)B.coeff(i)));
    }
    return NTL::to_double(logp);
}

long double Lattice::rho(const VectorXld B, const int n, const int m)
{
    long double S = 0, T = 0;
    for (int i = 0; i < n; ++i)
    {
        S += (i + 1) * log(B.coeff(i));
        T += log(B.coeff(i));
    }
    return 12 * (S - (n + 1) * T * 0.5) / (n * (n * n - 1));
}
