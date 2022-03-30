#include "ImplicitIntegrator.h"
#include "Deformation.hpp"
#include <Eigen/Dense>
#include <iostream>

using Eigen::MatrixXd;

namespace
{
    using FloatT = float;
    using Mat3 = Eigen::Matrix<FloatT, 3, 3>;
    using Mat9x12 = Eigen::Matrix<FloatT, 9, 12>;
    using Mat9 = Eigen::Matrix<FloatT, 9, 9>;
    using Mat12 = Eigen::Matrix<FloatT, 12, 12>;

    Deformation::TetraGroup* group;
    using namespace Deformation;

    Mat3 Calc_Ds(const Deformation::Tetrahedra& tet, const std::vector<Vertex> & vertices)
    {
        Mat3 ds;
        for (int i = 1; i < 4; i++)
            for (int j = 0; j < 3; j++)
                ds(i - 1, j) = vertices[tet.mIndices[i]].mPosition[j] - vertices[tet.mIndices[0]].mPosition[j];
        return ds;
    }

    Mat3 Calc_DmInv(const Deformation::Tetrahedra& tet, const std::vector<Vertex>& vertices)
    {
        return Calc_Ds(tet, vertices).inverse();
    }

    Mat3 Calc_F(const Mat3& ds, const Mat3& dmInv)
    {
        return ds * dmInv;
    }

    Mat3 Calc_E(const Mat3& f)
    {
        return 0.5f * (f * f.transpose() - Mat3::Identity());
    }

    Mat3 Calc_dPsiDf(const Mat3& f, const Mat3& e)
    {
        FloatT youngsModulus = 8 * 10e8f;
        FloatT poissonRatio = 0.45f;
        FloatT lamba = youngsModulus * poissonRatio / ((1.0f + poissonRatio) * (1.0f - 2.0f * poissonRatio));
        FloatT mu = youngsModulus / (2.0f + 2.0f * poissonRatio);
        return mu * f * e + lamba * e.trace() * f;
    }

    Mat9x12 Calc_dFdx(const Mat3& dmInv)
    {
        auto m = dmInv(1, 1);
        auto n = dmInv(1, 2);
        auto o = dmInv(1, 3);
        auto p = dmInv(2, 1);
        auto q = dmInv(2, 2);
        auto r = dmInv(2, 3);
        auto s = dmInv(3, 1);
        auto t = dmInv(3, 2);
        auto u = dmInv(3, 3);

        auto t1 = -m - p - s;
        auto t2 = -n - q - t;
        auto t3 = -o - r - u;

        Mat9x12 dFdx = Mat9x12::Zero();

        dFdx(1, 1) = t1;
        dFdx(1, 4) = m;
        dFdx(1, 7) = p;
        dFdx(1, 10) = s;
        dFdx(2, 2) = t1;
        dFdx(2, 5) = m;
        dFdx(2, 8) = p;
        dFdx(2, 11) = s;
        dFdx(3, 3) = t1;
        dFdx(3, 6) = m;
        dFdx(3, 9) = p;
        dFdx(3, 12) = s;
        dFdx(4, 1) = t2;
        dFdx(4, 4) = n;
        dFdx(4, 7) = q;
        dFdx(4, 10) = t;
        dFdx(5, 2) = t2;
        dFdx(5, 5) = n;
        dFdx(5, 8) = q;
        dFdx(5, 11) = t;
        dFdx(6, 3) = t2;
        dFdx(6, 6) = n;
        dFdx(6, 9) = q;
        dFdx(6, 12) = t;
        dFdx(7, 1) = t3;
        dFdx(7, 4) = o;
        dFdx(7, 7) = r;
        dFdx(7, 10) = u;
        dFdx(8, 2) = t3;
        dFdx(8, 5) = o;
        dFdx(8, 8) = r;
        dFdx(8, 11) = u;
        dFdx(9, 3) = t3;
        dFdx(9, 6) = o;
        dFdx(9, 9) = r;
        dFdx(9, 12) = u;

        return dFdx;
    }

    Eigen::Matrix<FloatT, 9, 1> Calc_g_i(const Mat3& f)
    {
        Eigen::Matrix<FloatT, 9, 1> g_i;
        g_i(0, 0) = 2 * f(0, 0);
        g_i(1, 0) = 2 * f(1, 0);
        g_i(2, 0) = 2 * f(2, 0);
                    
        g_i(3, 0) = 2 * f(0, 1);
        g_i(4, 0) = 2 * f(1, 1);
        g_i(5, 0) = 2 * f(2, 1);
                    
        g_i(6, 0) = 2 * f(0, 2);
        g_i(7, 0) = 2 * f(1, 2);
        g_i(8, 0) = 2 * f(2, 2);
        return g_i;
    }

    FloatT Calc_i_c(const Mat3& f)
    {
        auto frob = f.norm();
        return frob * frob;
    }

    Mat9 Calc_h_i()
    {
        return 2 * Mat9::Identity();
    }

    Mat9 Calc_D(const Mat3& f)
    {
        Mat9 d = Mat9::Zero();
        const auto& f0 = f.col(0);
        const auto& f1 = f.col(1);
        const auto& f2 = f.col(2);

        d.block(0, 0, 3, 3) = f0 * f0.transpose();
        d.block(3, 0, 3, 3) = f0 * f1.transpose();
        d.block(6, 0, 3, 3) = f0 * f2.transpose();

        d.block(0, 3, 3, 3) = f1 * f0.transpose();
        d.block(3, 3, 3, 3) = f1 * f1.transpose();
        d.block(6, 3, 3, 3) = f1 * f2.transpose();

        d.block(0, 6, 3, 3) = f2 * f0.transpose();
        d.block(3, 6, 3, 3) = f2 * f1.transpose();
        d.block(6, 6, 3, 3) = f2 * f2.transpose();

        return d;
    }

    Mat9 KroneckerProduct(const Mat3& a, const Mat3& b)
    {
        Mat9 product = Mat9::Zero();
        product.block(0, 0, 3, 3) = a(0, 0) * b;
        product.block(3, 0, 3, 3) = a(1, 0) * b;
        product.block(6, 0, 3, 3) = a(2, 0) * b;

        product.block(0, 3, 3, 3) = a(0, 1) * b;
        product.block(3, 3, 3, 3) = a(1, 1) * b;
        product.block(6, 3, 3, 3) = a(2, 1) * b;

        product.block(0, 6, 3, 3) = a(0, 2) * b;
        product.block(3, 6, 3, 3) = a(1, 2) * b;
        product.block(6, 6, 3, 3) = a(2, 2) * b;

        return product;
    }

    Mat9 Calc_h_2(const Mat3& f, const Mat9& d)
    {
        return 4 * (KroneckerProduct(Mat3::Identity(), f * f.transpose()) + KroneckerProduct(f * f.transpose(), Mat3::Identity()) + d);
    }

    Mat9 Calc_vec_dPsi2_dF2(const Mat9& g_i, FloatT i_c, const Mat9& h_i, const Mat9& h_2, FloatT mu, FloatT lambda)
    {
        return lambda / 4.0f * (g_i * g_i.transpose()) + (-mu / 2.0f + lambda / 4.0f * (i_c - 3.0f)) * h_i + mu / 4.0f * h_2;
    }

    Mat12 Calc_dfdx(const Mat9x12& dFdx, const Mat9 vec_dPsi2_dF2, FloatT a)
    {
        return -a * dFdx.transpose() * vec_dPsi2_dF2 * dFdx;
    }
}

////////////////////////////////////////////////////////////////////////////////