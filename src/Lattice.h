/**
 * @file Lattice.h
 * @author 佐藤 新 (Arata Sato) (transcendentalliouville@gmail.com)
 * @brief　PotBKZ，及び自己双対型PotBKZ，及びそれを開発するに当たり必要となったアルゴリズムが詰めてあるもの
 * @version 0.1
 * @date 2025-01-18
 *
 * @copyright Copyright (c) 2025
 *
 */

#ifndef LATTICE_H
#define LATTICE_H

#define LOG099 -0.010050335853501441183548857558547706085515007674629873378

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <NTL/RR.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>

typedef Eigen::Matrix<long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXli;
typedef Eigen::Matrix<long, 1, Eigen::Dynamic> VectorXli;
typedef Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXld;
typedef Eigen::Matrix<long double, 1, Eigen::Dynamic> VectorXld;

/**
 * @brief 格子に関わるclass
 *
 */
class Lattice
{
private:
    NTL::ZZ _;
    int _n;                                      // 基底行列の行数
    int _m;                                      // 基底行列の列数
    VectorXld _squared_norm_of_gso_vec;          // GSOベクトルの二乗ノルム
    VectorXld _log_squared_norm_of_gso_vec;      // GSOベクトルの二乗ノルムの対数値
    MatrixXld _gso_coeff_mat;                    // GSO係数行列
    MatrixXld _gso_vec_mat;                      // GSO行列
    VectorXld _dual_squared_norm_of_gso_vec;     // 双対基底のGSOベクトルの二乗ノルム
    VectorXld _dual_log_squared_norm_of_gso_vec; // 双対基底のGSOベクトルの二乗ノルムの対数値
    MatrixXld _dual_gso_coeff_mat;               // 双対基底のGSO係数行列

    /**
     * @brief 双対型deep-insetionを行う関数
     *
     */
    void DualDeepInsertion(const int k, const int l)
    {
        const VectorXli t = basis.row(k);
        for (int j = k; j < l; ++j)
        {
            basis.row(j) = basis.row(j + 1);
        }
        basis.row(l) = t;
    }

    /**
     * @brief 双対基底への格子ベクトルの挿入を行う関数
     *
     * @param _basis 基底
     * @param x 格子ベクトル
     * @param n 格子の次元
     * @param m 基底行列の列
     * @return MatrixXli 挿入された基底
     */
    MatrixXli Insert(const MatrixXli _basis, const VectorXli x, const int n, const int m)
    {
        int i, j;
        double beta = 1.35135135135135135135135135135135135135135135135135;
        double tmp, gamma;
        MatrixXli U(n, n), c;
        U.setZero();
        NTL::mat_ZZ tmp_base;
        tmp_base.SetDims(n, n + 1);

        /* Construction of gamma */
        tmp = x.cast<double>().norm();
        tmp *= pow(beta, (n - 2) * 0.5);
        gamma = round(tmp + tmp);

        /* Construction of matrix */
        for (i = 0; i < n; ++i)
        {
            for (j = 0; j < n; ++j)
            {
                tmp_base[i][j] = 0;
            }
            tmp_base[i][i] = 1;
            tmp_base[i][n] = gamma * x.coeff(i);
        }
        NTL::LLL(_, tmp_base, 99, 100);

        for (i = 0; i < n; ++i)
        {
            for (j = 0; j < n; ++j)
            {
                U.coeffRef(i, j) = NTL::to_long(tmp_base[i][j]);
            }
        }

        /* Construction of a new basis */
        return U * _basis;
    }

    void GSO()
    {
        _gso_vec_mat.resize(_n, _m);

        for (int i = 0, j; i < _n; ++i)
        {
            _gso_coeff_mat.coeffRef(i, i) = 1.0;
            _gso_vec_mat.row(i) = basis.row(i).cast<long double>();
            for (j = 0; j < i; ++j)
            {
                _gso_coeff_mat.coeffRef(i, j) = basis.row(i).cast<long double>().dot(_gso_vec_mat.row(j)) / _gso_vec_mat.row(j).dot(_gso_vec_mat.row(j));
                _gso_vec_mat.row(i) -= _gso_coeff_mat.coeff(i, j) * _gso_vec_mat.row(j);
            }
            _squared_norm_of_gso_vec.coeffRef(i) = _gso_vec_mat.row(i).dot(_gso_vec_mat.row(i));
            _log_squared_norm_of_gso_vec.coeffRef(i) = log(_squared_norm_of_gso_vec.coeff(i));
        }
    }

    /**
     * @brief primalな基底のGSO情報から双対基底のGSO情報を計算する関数
     *
     * @param squared_norm_of_gso_vec GSOベクトルの二乗ノルム
     * @param log_squared_norm_of_gso_vec GSOベクトルの二乗ノルムの対数値
     * @param gso_coeff_mat GSO係数行列
     * @param dual_squared_norm_of_gso_vec 双対基底のGSOベクトルの二乗ノルム
     * @param dual_log_squared_norm_of_gso_vec 双対基底のGSOベクトルの二乗ノルムの対数値
     * @param dual_gso_coeff_mat 双対基底のGSO係数行列
     * @param n 格子の次元
     * @param m 基底行列の列数
     */
    void DualGSO(const VectorXld squared_norm_of_gso_vec,
                 const VectorXld log_squared_norm_of_gso_vec,
                 const MatrixXld gso_coeff_mat,
                 VectorXld &dual_squared_norm_of_gso_vec,
                 VectorXld &dual_log_squared_norm_of_gso_vec,
                 MatrixXld &dual_gso_coeff_mat,
                 const int n, const int m)
    {
        for (int i = 0, j; i < n; ++i)
        {
            dual_squared_norm_of_gso_vec.coeffRef(i) = 1.0 / squared_norm_of_gso_vec.coeff(i);
            dual_log_squared_norm_of_gso_vec.coeffRef(i) = -log_squared_norm_of_gso_vec.coeff(i);
            for (j = i + 1; j < m; ++j)
            {
                dual_gso_coeff_mat.coeffRef(i, j) = -gso_coeff_mat.row(j).segment(i, j - i).dot(dual_gso_coeff_mat.row(i).segment(i, j - i));
            }
        }
    }

    /**
     * @brief PotENUMアルゴリズム
     *
     * @param _gso_coeff_mat GSO係数行列
     * @param _squared_norm_of_gso_vec GSOベクトルの二乗ノルム
     * @param _log_squared_norm_of_gso_vec GSOベクトルの二乗ノルムの対数値
     * @param n 格子の次元
     * @return VectorXli $delta$-特異な格子ベクトル
     */
    VectorXli PotENUM(const MatrixXld _gso_coeff_mat, const VectorXld _squared_norm_of_gso_vec, const VectorXld _log_squared_norm_of_gso_vec, const int n);

    /**
     * @brief 双対型PotENUM
     *
     * @param _gso_coeff_mat GSO係数行列
     * @param _squared_norm_of_gso_vec GSOベクトルの二乗ノルム
     * @param _log_squared_norm_of_gso_vec GSOベクトルの二乗ノルムの対数値
     * @param n 格子次元
     * @return VectorXli
     */
    VectorXli DualPotENUM(const MatrixXld _gso_coeff_mat, const VectorXld _squared_norm_of_gso_vec, const VectorXld _log_squared_norm_of_gso_vec, const int n);

public:
    MatrixXli basis; // 基底行列

    /**
     * @brief 基底行列の行数と列数を設定する関数
     *
     * @param n 基底行列の行数
     * @param m 基底行列の列数
     */
    void setDims(const int n, const int m)
    {
        _n = m;                                      // 基底行列の行数
        _m = m;                                      // 基底行列の列数
        basis.resize(n, m);                          // 基底行列
        _squared_norm_of_gso_vec.resize(n);          // GSOベクトルの二乗ノルム
        _log_squared_norm_of_gso_vec.resize(n);      // GSOベクトルの二乗ノルムの対数値
        _gso_coeff_mat.resize(n, n);                 // GSO係数行列
        _gso_vec_mat.resize(n, m);                   // GSO行列
        _dual_squared_norm_of_gso_vec.resize(n);     // 双対基底のGSOベクトルの二乗ノルム
        _dual_log_squared_norm_of_gso_vec.resize(n); // 双対基底のGSOベクトルの二乗ノルムの対数値
        _dual_gso_coeff_mat.resize(n, n);            // 双対基底のGSO係数行列
    }

    /**
     * @brief PotLLLアルゴリズム
     *
     * @param reduction_parameter 簡約変数
     */
    void POT_LLL(const double reduction_parameter);

    /**
     * @brief 双対型PotLLLアルゴリズム
     *
     * @param reduction_parameter 簡約変数
     */
    void DUAL_POT_LLL(const double reduction_parameter);

    /**
     * @brief PotBKZアルゴリズム
     *
     * @param beta ブロックサイズ
     * @param reduction_parameter 簡約変数
     */
    void POT_BKZ(const int beta, const double reduction_parameter);

    /**
     * @brief 双対型PotBKZアルゴリズム
     *
     * @param beta ブロックサイズ
     * @param delta 簡約変数
     */
    void DUAL_POT_BKZ(const int beta, const double delta);

    /**
     * @brief 自己双対型PotBKZアルゴリズム
     *
     * @param beta ブロックサイズ
     * @param reduction_parameter 簡約変数
     */
    void SELF_DUAL_POT_BKZ(const int beta, const double reduction_parameter);
};

#endif // !LATTICE_H
