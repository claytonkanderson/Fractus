#include "ImplicitIntegrator.h"
#include "Deformation.hpp"
#include "ProtoConverter.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SVD>
#include <Eigen/IterativeLinearSolvers>
#include <Mathematics/SymmetricEigensolver3x3.h>

#include <iostream>
#include <glm/gtc/constants.hpp>

using Eigen::MatrixXd;

namespace Deformation
{
    using FloatT = double;
    using Mat3 = Eigen::Matrix<FloatT, 3, 3>;
    using Mat3x4 = Eigen::Matrix<FloatT, 3, 4>;
    using Mat9x12 = Eigen::Matrix<FloatT, 9, 12>;
    using Mat9 = Eigen::Matrix<FloatT, 9, 9>;
    using Mat12 = Eigen::Matrix<FloatT, 12, 12>;
    using Vec3 = Eigen::Matrix<FloatT, 3, 1>;
    using Vec9 = Eigen::Matrix<FloatT, 9, 1>;
    using Vec12 = Eigen::Matrix<FloatT, 12, 1>;

    Mat3 Calc_Ds(const Deformation::Tetrahedra& tet, const std::vector<Vertex>& vertices, bool restState)
    {
        Mat3 ds = Mat3::Zero();

        for (int row = 0; row < 3; row++)
        {
            for (int col = 1; col < 4; col++)
            {
                if (restState)
					ds(row, col - 1) = vertices[tet.mIndices[col]].mMaterialCoordinates[row] - vertices[tet.mIndices[0]].mMaterialCoordinates[row];
                else
                    ds(row, col - 1) = vertices[tet.mIndices[col]].mPosition[row] - vertices[tet.mIndices[0]].mPosition[row];
            }
        }

        return ds;
    }

    Mat3 Calc_Ds(const std::array<glm::dvec3,4> & positions)
    {
        Mat3 ds = Mat3::Zero();

        for (int row = 0; row < 3; row++)
        {
            for (int col = 1; col < 4; col++)
            {
				ds(row, col - 1) = positions[col][row] - positions[0][row];
            }
        }

        return ds;
    }

    Mat3 Calc_DmInv(const Deformation::Tetrahedra& tet, const std::vector<Vertex>& vertices)
    {
        auto mat = Calc_Ds(tet, vertices, true);
        if (std::abs(mat.determinant()) < 1e-10f)
            throw std::exception("Small");
        return mat.inverse();
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

    Vec9 Reshape3x3(const Mat3& m)
    {
        Vec9 vec = Vec9::Zero();
        vec(0, 0) = m(0, 0);
        vec(1, 0) = m(1, 0);
        vec(2, 0) = m(2, 0);

        vec(3, 0) = m(0, 1);
        vec(4, 0) = m(1, 1);
        vec(5, 0) = m(2, 1);

        vec(6, 0) = m(0, 2);
        vec(7, 0) = m(1, 2);
        vec(8, 0) = m(2, 2);

        return vec;
    }

    Vec12 Reshape3x4(const Mat3x4& m)
    {
        Vec12 vec = Vec12::Zero();
        vec(0, 0) = m(0, 0);
        vec(1, 0) = m(1, 0);
        vec(2, 0) = m(2, 0);

        vec(3, 0) = m(0, 1);
        vec(4, 0) = m(1, 1);
        vec(5, 0) = m(2, 1);

        vec(6, 0) = m(0, 2);
        vec(7, 0) = m(1, 2);
        vec(8, 0) = m(2, 2);

        vec(9, 0) = m(0, 3);
        vec(10, 0) = m(1, 3);
        vec(11, 0) = m(2, 3);
        return vec;
    }

    Mat9x12 Calc_dFdx(const Mat3& dmInv)
    {
        auto m = dmInv(0, 0);
        auto n = dmInv(0, 1);
        auto o = dmInv(0, 2);
        auto p = dmInv(1, 0);
        auto q = dmInv(1, 1);
        auto r = dmInv(1, 2);
        auto s = dmInv(2, 0);
        auto t = dmInv(2, 1);
        auto u = dmInv(2, 2);

        auto t1 = -m - p - s;
        auto t2 = -n - q - t;
        auto t3 = -o - r - u;

        Mat9x12 dFdx = Mat9x12::Zero();

        dFdx(0, 0) = t1;
        dFdx(0, 3) = m;
        dFdx(0, 6) = p;
        dFdx(0, 9) = s;
        dFdx(1, 1) = t1;
        dFdx(1, 4) = m;
        dFdx(1, 7) = p;
        dFdx(1, 10) = s;

        dFdx(2, 2) = t1;
        dFdx(2, 5) = m;
        dFdx(2, 8) = p;
        dFdx(2, 11) = s;
        dFdx(3, 0) = t2;
        dFdx(3, 3) = n;
        dFdx(3, 6) = q;
        dFdx(3, 9) = t;
        dFdx(4, 1) = t2;
        dFdx(4, 4) = n;

        dFdx(4, 7) = q;
        dFdx(4, 10) = t;
        dFdx(5, 2) = t2;
        dFdx(5, 5) = n;
        dFdx(5, 8) = q;
        dFdx(5, 11) = t;
        dFdx(6, 0) = t3;
        dFdx(6, 3) = o;
        dFdx(6, 6) = r;
        dFdx(6, 9) = u;
        dFdx(7, 1) = t3;
        dFdx(7, 4) = o;

        dFdx(7, 7) = r;
        dFdx(7, 10) = u;
        dFdx(8, 2) = t3;
        dFdx(8, 5) = o;
        dFdx(8, 8) = r;
        dFdx(8, 11) = u;

        return dFdx;
    }

    Vec9 Calc_g_i(const Mat3& f)
    {
        Vec9 g_i = Vec9::Zero();
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

    Mat3 Calc_dj_df(const Mat3& f)
    {
        Mat3 dj_df = Mat3::Zero();
        dj_df.col(0) = f.col(1).cross(f.col(2));
        dj_df.col(1) = f.col(2).cross(f.col(0));
        dj_df.col(2) = f.col(0).cross(f.col(1));

        //dj_df.block(0, 0, 3, 1) = f.col(1).cross(f.col(2));
        //dj_df.block(0, 1, 3, 1) = f.col(2).cross(f.col(0));
        //dj_df.block(0, 2, 3, 1) = f.col(0).cross(f.col(1));
        return dj_df;
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

    Mat3 Calc_f_hat(const Vec3& f_i)
    {
        Mat3 f_hat = Mat3::Zero();
        // col 0
        f_hat(1, 0) = f_i(2);
        f_hat(2, 0) = -f_i(1);
        
        // col 1
        f_hat(0, 1) = -f_i(2);
        f_hat(2, 1) = f_i(0);

        // col 2
        f_hat(0, 2) = f_i(1);
        f_hat(1, 2) = -f_i(0);

        return f_hat;
    }

    Mat9 Calc_h_j(const Mat3& f)
    {
        Mat9 h_j = Mat9::Zero();
        auto f0_hat = Calc_f_hat(f.col(0));
        auto f1_hat = Calc_f_hat(f.col(1));
        auto f2_hat = Calc_f_hat(f.col(2));

        // col 0
        h_j.block(3, 0, 3, 3) = f2_hat;
        h_j.block(6, 0, 3, 3) = -f1_hat;

        // col 1
        h_j.block(0, 3, 3, 3) = -f2_hat;
        h_j.block(6, 3, 3, 3) = f0_hat;

        // col 2
        h_j.block(0, 6, 3, 3) = f1_hat;
        h_j.block(3, 6, 3, 3) = -f0_hat;
        
        return h_j;
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

    Mat9 Calc_vec_dPsi2_dF2(const Vec9& g_i, FloatT i_c, const Mat9& h_i, const Mat9& h_2, FloatT mu, FloatT lambda)
    {
        return lambda / 4.0f * (g_i * g_i.transpose()) + (-mu / 2.0f + lambda / 4.0f * (i_c - 3.0f)) * h_i + mu / 4.0f * h_2;
    }

    Mat12 Calc_dfdx(const Mat9x12& dFdx, const Mat9 & vec_dPsi2_dF2, FloatT restVol)
    {
        return -restVol * dFdx.transpose() * vec_dPsi2_dF2 * dFdx;
    }

    Mat3 Get_Q(int i, const std::array<FloatT, 3>& epsilon, const std::array<FloatT, 3>& sigma, const Mat3& u, const Mat3& v_transpose)
    {
        Mat3 d = Mat3::Zero();
        d(0, 0) = sigma[0] * sigma[2] + sigma[1] * epsilon[i];
        d(1, 1) = sigma[1] * sigma[2] + sigma[0] * epsilon[i];
		d(2, 2) = epsilon[i] * epsilon[i] - sigma[2] * sigma[2];
        return 1.0f / d.norm() * u * d * v_transpose;
    }

    Mat3 Get_Q(const Mat3 & U, const Mat3 & V, int i)
    {
        //Mat3 q = Mat3::Zero();

        //if (i == 3)
        //{
        //    q(0, 1) = -1;
        //    q(1, 0) = 1;
        //}
        //else if (i == 4)
        //{
        //    q(1, 2) = 1;
        //    q(2, 1) = -1;
        //}
        //else if (i == 5)
        //{
        //    q(0, 2) = 1;
        //    q(2, 0) = -1;
        //}
        //else if (i == 6)
        //{
        //    q(0, 1) = 1;
        //    q(1, 0) = 1;
        //}
        //else if (i == 7)
        //{
        //    q(1, 2) = 1;
        //    q(2, 1) = 1;
        //}
        //else if (i == 8)
        //{
        //    q(0, 2) = 1;
        //    q(2, 0) = 1;
        //}

        //return 1 / sqrt(2.0f) * u * q * v_transpose;

        static const FloatT scale = 1.0 / std::sqrt(2.0);
        const Mat3 sV = scale * V;

        using M3 = Eigen::Matrix<FloatT, 3, 3, Eigen::ColMajor>;

        M3 A;
        A << sV(0, 2) * U(0, 1), sV(1, 2)* U(0, 1), sV(2, 2)* U(0, 1),
            sV(0, 2)* U(1, 1), sV(1, 2)* U(1, 1), sV(2, 2)* U(1, 1),
            sV(0, 2)* U(2, 1), sV(1, 2)* U(2, 1), sV(2, 2)* U(2, 1);

        M3 B;
        B << sV(0, 1) * U(0, 2), sV(1, 1)* U(0, 2), sV(2, 1)* U(0, 2),
            sV(0, 1)* U(1, 2), sV(1, 1)* U(1, 2), sV(2, 1)* U(1, 2),
            sV(0, 1)* U(2, 2), sV(1, 1)* U(2, 2), sV(2, 1)* U(2, 2);

        M3 C;
        C << sV(0, 2) * U(0, 0), sV(1, 2)* U(0, 0), sV(2, 2)* U(0, 0),
            sV(0, 2)* U(1, 0), sV(1, 2)* U(1, 0), sV(2, 2)* U(1, 0),
            sV(0, 2)* U(2, 0), sV(1, 2)* U(2, 0), sV(2, 2)* U(2, 0);

        M3 D;
        D << sV(0, 0) * U(0, 2), sV(1, 0)* U(0, 2), sV(2, 0)* U(0, 2),
            sV(0, 0)* U(1, 2), sV(1, 0)* U(1, 2), sV(2, 0)* U(1, 2),
            sV(0, 0)* U(2, 2), sV(1, 0)* U(2, 2), sV(2, 0)* U(2, 2);

        M3 E;
        E << sV(0, 1) * U(0, 0), sV(1, 1)* U(0, 0), sV(2, 1)* U(0, 0),
            sV(0, 1)* U(1, 0), sV(1, 1)* U(1, 0), sV(2, 1)* U(1, 0),
            sV(0, 1)* U(2, 0), sV(1, 1)* U(2, 0), sV(2, 1)* U(2, 0);

        M3 F;
        F << sV(0, 0) * U(0, 1), sV(1, 0)* U(0, 1), sV(2, 0)* U(0, 1),
            sV(0, 0)* U(1, 1), sV(1, 0)* U(1, 1), sV(2, 0)* U(1, 1),
            sV(0, 0)* U(2, 1), sV(1, 0)* U(2, 1), sV(2, 0)* U(2, 1);

        if (i == 3)
            return B - A;
        else if (i == 4)
            return D - C;
        else if (i == 5)
            return F - E;
        else if (i == 6)
            return A + B;
        else if (i == 7)
            return C + D;
        else if (i == 8)
            return E + F;
    }

    FloatT Calc_Epsilon(int i, FloatT I2, FloatT I3)
    {
        return 2 * sqrtf(I2 / 3.0f) * cosf(1.0f / 3.0f * (acosf(3*I3/I2*sqrtf(3/I2)) + 2*glm::pi<FloatT>()*(i-1)));
    }

    void BuildTwistAndFlipEigenvectors(const Mat3& U, const Mat3& V, Mat9& Q)
    {
        static const FloatT scale = 1.0 / std::sqrt(2.0);
        const Mat3 sV = scale * V;

        using M3 = Eigen::Matrix<FloatT, 3, 3, Eigen::ColMajor>;

        M3 A;
        A << sV(0, 2) * U(0, 1), sV(1, 2)* U(0, 1), sV(2, 2)* U(0, 1),
            sV(0, 2)* U(1, 1), sV(1, 2)* U(1, 1), sV(2, 2)* U(1, 1),
            sV(0, 2)* U(2, 1), sV(1, 2)* U(2, 1), sV(2, 2)* U(2, 1);

        M3 B;
        B << sV(0, 1) * U(0, 2), sV(1, 1)* U(0, 2), sV(2, 1)* U(0, 2),
            sV(0, 1)* U(1, 2), sV(1, 1)* U(1, 2), sV(2, 1)* U(1, 2),
            sV(0, 1)* U(2, 2), sV(1, 1)* U(2, 2), sV(2, 1)* U(2, 2);

        M3 C;
        C << sV(0, 2) * U(0, 0), sV(1, 2)* U(0, 0), sV(2, 2)* U(0, 0),
            sV(0, 2)* U(1, 0), sV(1, 2)* U(1, 0), sV(2, 2)* U(1, 0),
            sV(0, 2)* U(2, 0), sV(1, 2)* U(2, 0), sV(2, 2)* U(2, 0);

        M3 D;
        D << sV(0, 0) * U(0, 2), sV(1, 0)* U(0, 2), sV(2, 0)* U(0, 2),
            sV(0, 0)* U(1, 2), sV(1, 0)* U(1, 2), sV(2, 0)* U(1, 2),
            sV(0, 0)* U(2, 2), sV(1, 0)* U(2, 2), sV(2, 0)* U(2, 2);

        M3 E;
        E << sV(0, 1) * U(0, 0), sV(1, 1)* U(0, 0), sV(2, 1)* U(0, 0),
            sV(0, 1)* U(1, 0), sV(1, 1)* U(1, 0), sV(2, 1)* U(1, 0),
            sV(0, 1)* U(2, 0), sV(1, 1)* U(2, 0), sV(2, 1)* U(2, 0);

        M3 F;
        F << sV(0, 0) * U(0, 1), sV(1, 0)* U(0, 1), sV(2, 0)* U(0, 1),
            sV(0, 0)* U(1, 1), sV(1, 0)* U(1, 1), sV(2, 0)* U(1, 1),
            sV(0, 0)* U(2, 1), sV(1, 0)* U(2, 1), sV(2, 0)* U(2, 1);

        // Twist eigenvectors
        //std::cout << "eigenvectors before mapping" << std::endl;
        //std::cout << Q << std::endl;
        //std::cout << "B-A" << std::endl;
        //std::cout << B - A << std::endl;
        Eigen::Map<M3>(Q.data()) = B - A;
        //std::cout << "eigenvectors after single mapping" << std::endl;
        //std::cout << Q << std::endl;
        Eigen::Map<M3>(Q.data() + 9) = D - C;
        Eigen::Map<M3>(Q.data() + 18) = F - E;

        // Flip eigenvectors
        Eigen::Map<M3>(Q.data() + 27) = A + B;
        Eigen::Map<M3>(Q.data() + 36) = C + D;
        Eigen::Map<M3>(Q.data() + 45) = E + F;
    }

    Mat9 CalcProjectedHessian(FloatT mu, FloatT lambda, const Mat3& F, const Mat3& U, const Mat3& V, const Vec3& S)
    {
        Vec9 eigenvalues;
        Mat9 eigenvectors;

        const FloatT J = F.determinant();

        // Compute the twist and flip eigenvalues
        {
            // Twist eigenvalues
            eigenvalues.segment<3>(0) = S;
            // Flip eigenvalues
            eigenvalues.segment<3>(3) = -S;
            const FloatT evScale = lambda * (J - 1.0) - mu;
            eigenvalues.segment<6>(0) *= evScale;
            eigenvalues.segment<6>(0).array() += mu;
        }

        // Compute the twist and flip eigenvectors
        BuildTwistAndFlipEigenvectors(U, V, eigenvectors);

        // Compute the remaining three eigenvalues and eigenvectors
        {
            Mat3 A;
            const FloatT s0s0 = S(0) * S(0);
            const FloatT s1s1 = S(1) * S(1);
            const FloatT s2s2 = S(2) * S(2);
            A(0, 0) = mu + lambda * s1s1 * s2s2;
            A(1, 1) = mu + lambda * s0s0 * s2s2;
            A(2, 2) = mu + lambda * s0s0 * s1s1;
            const FloatT evScale = lambda * (2.0 * J - 1.0) - mu;
            A(0, 1) = evScale * S(2);
            A(1, 0) = A(0, 1);
            A(0, 2) = evScale * S(1);
            A(2, 0) = A(0, 2);
            A(1, 2) = evScale * S(0);
            A(2, 1) = A(1, 2);

            //std::cout << "A" << std::endl;
            //std::cout << A << std::endl;

            const Eigen::SelfAdjointEigenSolver<Mat3> Aeigs(A);
            eigenvalues.segment<3>(6) = Aeigs.eigenvalues();

   //         std::cout << "u" << std::endl;
   //         std::cout << U << std::endl;

   //         std::cout << "v" << std::endl;
   //         std::cout << V << std::endl;

   //         std::cout << "A eigs" << std::endl;
   //         std::cout << Aeigs.eigenvalues() << std::endl;
   //         std::cout << "Unrotated eig vectors" << std::endl;
			//std::cout << Aeigs.eigenvectors() << std::endl;

   //         std::cout << "Rotated eigen 0 matrix" << std::endl;
   //         std::cout << U * Aeigs.eigenvectors().col(0).asDiagonal() * V.transpose() << std::endl;

   //         std::cout << "eigenvectors before mapping" << std::endl;
   //         std::cout << eigenvectors << std::endl;
            Eigen::Map<Mat3>(eigenvectors.data() + 54) = U * Aeigs.eigenvectors().col(0).asDiagonal() * V.transpose();
            //std::cout << "eigenvectors after single mapping" << std::endl;
            //std::cout << eigenvectors << std::endl;
            Eigen::Map<Mat3>(eigenvectors.data() + 63) = U * Aeigs.eigenvectors().col(1).asDiagonal() * V.transpose();
            Eigen::Map<Mat3>(eigenvectors.data() + 72) = U * Aeigs.eigenvectors().col(2).asDiagonal() * V.transpose();
        }

        // Clamp the eigenvalues
        for (int i = 0; i < 9; i++)
        {
            //std::cout << "lambda_" << i << " " << eigenvalues(i) << std::endl;

            if (eigenvalues(i) < 0.0)
            {
                eigenvalues(i) = 0.0;
                //eigenvalues(i) = -eigenvalues(i);

            }
        }
        //std::cout << "eigenvalues" << std::endl;
        //std::cout << eigenvalues << std::endl;

        //std::cout << "eigenvectors" << std::endl;
        //std::cout << eigenvectors << std::endl;

        return eigenvectors * eigenvalues.asDiagonal() * eigenvectors.transpose();
    }

    Mat3x4 CalcBm(const std::array<size_t, 4>& vertIndices, const std::vector<Vertex>& vertices)
    {
        std::vector<Vec3> vertices_converted;
        glm::ivec4 tetIndices({ vertIndices[0], vertIndices[1], vertIndices[2], vertIndices[3] });

        for (int i = 0; i < 4; i++)
        {
            const auto& pos = vertices[tetIndices[i]].mMaterialCoordinates;
            vertices_converted.push_back(Vec3(pos.x, pos.y, pos.z));
        }

        // TODO: Eliminate this class, it is just used to get the area normal
        class Triangle final
        {
        public:
            Triangle(const Vec3& v0, const Vec3& v1, const Vec3& v2)
                : _x0(v0)
                , _x1(v1)
                , _x2(v2)
            {}

            // TODO: Make this a constructor
            static Triangle GenerateFace(const int faceNum, const glm::ivec4& tet, const std::vector<Vec3>& verts)
            {
                assert(faceNum >= 0); assert(faceNum < 4);
                if (faceNum == 0)
                {
                    return Triangle(verts[static_cast<unsigned long>(tet[0])], verts[static_cast<unsigned long>(tet[1])], verts[static_cast<unsigned long>(tet[3])]);
                }
                else if (faceNum == 1)
                {
                    return Triangle(verts[static_cast<unsigned long>(tet[0])], verts[static_cast<unsigned long>(tet[2])], verts[static_cast<unsigned long>(tet[1])]);
                }
                else if (faceNum == 2)
                {
                    return Triangle(verts[static_cast<unsigned long>(tet[3])], verts[static_cast<unsigned long>(tet[2])], verts[static_cast<unsigned long>(tet[0])]);
                }
                else if (faceNum == 3)
                {
                    return Triangle(verts[static_cast<unsigned long>(tet[1])], verts[static_cast<unsigned long>(tet[2])], verts[static_cast<unsigned long>(tet[3])]);
                }
                else
                {
                    std::cerr << "Error, impossible code path hit." << std::endl;
                    std::exit(EXIT_FAILURE);
                }
            }

            // TODO: Roll these into one computation
            Vec3 normal() const
            {
                return ((_x1 - _x0).cross(_x2 - _x0)).normalized();
            }
            FloatT area() const
            {
                return 0.5 * ((_x1 - _x0).cross(_x2 - _x0)).norm();
            }
        private:
            const Vec3& _x0;
            const Vec3& _x1;
            const Vec3& _x2;
        };

        Mat3x4 Bm;
		const Triangle tri0 = Triangle::GenerateFace(0, tetIndices, vertices_converted);
		const Triangle tri1 = Triangle::GenerateFace(1, tetIndices, vertices_converted);
		const Triangle tri2 = Triangle::GenerateFace(2, tetIndices, vertices_converted);
		const Triangle tri3 = Triangle::GenerateFace(3, tetIndices, vertices_converted);

		// Calculate the area vectors
		// v0 is incident on faces (0,1,2)
		Bm.col(0) = tri0.normal() * tri0.area() + tri1.normal() * tri1.area() + tri2.normal() * tri2.area();
		// v1 is incident on faces (0,1,3)
		Bm.col(1) = tri0.normal() * tri0.area() + tri1.normal() * tri1.area() + tri3.normal() * tri3.area();
		// v2 is incident on faces (1,2,3)
		Bm.col(2) = tri1.normal() * tri1.area() + tri2.normal() * tri2.area() + tri3.normal() * tri3.area();
		// v3 is incident on faces (0,2,3)
		Bm.col(3) = tri0.normal() * tri0.area() + tri2.normal() * tri2.area() + tri3.normal() * tri3.area();
		Bm /= -3.0;

        return Bm;
    }

    Eigen::VectorXf ComputePerturbedForce(const TetraGroup& group, const Eigen::VectorXf& deltaX)
    {
        size_t numVertices = group.mVertices.size();
        Eigen::VectorXf globalForce = Eigen::VectorXf::Zero(3 * numVertices);

        std::array<glm::dvec3, 4> perturbedPositions;

        for (const auto& pair : group.mIdToTetrahedra)
        {
            const auto& tet = pair.second;

            for (int i = 0; i < 4; i++)
            {
                auto id = tet.mIndices[i];
                perturbedPositions[i] = group.mVertices[id].mPosition + glm::dvec3(deltaX[3*id], deltaX[3 * id + 1], deltaX[3 * id + 2]);
            }

            auto ds = Calc_Ds(perturbedPositions);
            auto dmInv = Calc_DmInv(tet, group.mVertices);
            auto f = Calc_F(ds, dmInv);

            auto j = f.determinant();
            auto dj_df = Calc_dj_df(f);
            auto g_j = Reshape3x3(dj_df);
            Mat3 dPsi_dF = group.mMu * f + (group.mLambda * (j - 1.0f) - group.mMu) * dj_df;

            auto dF_dx = Calc_dFdx(dmInv);
            auto dPsi_dx = dF_dx.transpose() * Reshape3x3(dPsi_dF);
            Vec12 force = -tet.mRestVolume * dPsi_dx; // 12x1 atm

            Eigen::JacobiSVD<Mat3, Eigen::NoQRPreconditioner> svd_f(f, Eigen::ComputeFullU | Eigen::ComputeFullV);
            auto sigmas = svd_f.singularValues();
            auto u = svd_f.matrixU();
            auto v = svd_f.matrixV();

            if (u.determinant() < 0.0)
            {
                u.col(0) *= -1.0;
                sigmas(0) *= -1.0;
            }
            if (v.determinant() < 0.0)
            {
                v.col(0) *= -1.0;
                sigmas(0) *= -1.0;
            }

            // Contribution to globalB
            for (size_t i = 0; i < 4; i++)
            {
                auto vert_i = tet.mIndices[i];
                globalForce(3 * vert_i + 0) += force(3 * i + 0);
                globalForce(3 * vert_i + 1) += force(3 * i + 1);
                globalForce(3 * vert_i + 2) += force(3 * i + 2);
            }
        }

        return globalForce;
    }

    bool ImplicitUpdate(TetraGroup& group, float timestep, bool saveFrame, IronGames::SimulationFrame* frame, float deltaVThreshold)
    {
        size_t numVertices = group.mVertices.size();
        Eigen::SparseMatrix<float> globalA(3* numVertices, 3* numVertices);
        Eigen::VectorXf globalB = Eigen::VectorXf::Zero(3 * numVertices);
		Eigen::VectorXf globalX = Eigen::VectorXf::Zero(3 * numVertices);
        Vec12 nodeVelocities = Vec12::Zero();
        Vec12 localB = Vec12::Zero();

        Eigen::VectorXf globalForce = Eigen::VectorXf::Zero(3 * numVertices);
        Eigen::VectorXf globalVelocities = Eigen::VectorXf::Zero(3 * numVertices);
        Eigen::SparseMatrix<float> globalDfDx(3 * numVertices, 3 * numVertices);
        Eigen::SparseMatrix<float> massMatrix(3 * numVertices, 3 * numVertices);

        //std::cout << "timestep " << timestep << std::endl;

        for (size_t i = 0; i < group.mVertices.size(); i++)
        {
            globalA.coeffRef(3 * i, 3 * i) = group.mVertices[i].mMass;
            globalA.coeffRef(3 * i + 1, 3 * i + 1) = group.mVertices[i].mMass;
			globalA.coeffRef(3 * i + 2, 3 * i + 2) = group.mVertices[i].mMass;
        }

        for (const auto& pair : group.mIdToTetrahedra)
        {
            const auto& tet = pair.second;

            auto ds = Calc_Ds(tet, group.mVertices, false);
            //std::cout << "ds" << std::endl;
            //std::cout << ds << std::endl;
            auto dmInv = Calc_DmInv(tet, group.mVertices);
            //std::cout << "dmInv" << std::endl;
            //std::cout << dmInv << std::endl;
            auto f = Calc_F(ds, dmInv);
            //std::cout << "f" << std::endl;
            //std::cout << f << std::endl;

            // Stable Neo-hookean potential
			auto j = f.determinant();
			auto dj_df = Calc_dj_df(f); // SEEMS OK
            //std::cout << "dj_df" << std::endl;
            //std::cout << dj_df << std::endl;
			auto g_j = Reshape3x3(dj_df);
            //auto dPsi_dF = group.mMu * Reshape3x3(f) + group.mLambda * (j - 1.0f - group.mMu/group.mLambda) * g_j;
            Mat3 dPsi_dF = group.mMu * f + (group.mLambda * (j - 1.0f) - group.mMu) * dj_df;
            //auto elementForce = dPsi_dF * CalcBm(tet.mIndices, group.mVertices);

			auto dF_dx = Calc_dFdx(dmInv);
            //std::cout << "dfdx" << std::endl;
            //std::cout << dF_dx << std::endl;
			auto dPsi_dx = dF_dx.transpose() * Reshape3x3(dPsi_dF);
			Vec12 force = -tet.mRestVolume * dPsi_dx; // 12x1 atm

            //std::cout << "Force" << std::endl;
            //std::cout << force << std::endl;
            //force = -Reshape3x4(elementForce); // 12x1 atm

			auto h_j = Calc_h_j(f);
            //std::cout << "h_j" << std::endl;
            //std::cout << h_j << std::endl;
            //std::cout << "g_j * g_j^T" << std::endl;
            //std::cout << (g_j * g_j.transpose()) << std::endl;
			//Mat9 dPsi2_dF2 = group.mMu * Mat9::Identity() + group.mLambda * (j - 1.0f - group.mMu / group.mLambda) * h_j + group.mLambda * g_j * g_j.transpose();
            
            Eigen::JacobiSVD<Mat3, Eigen::NoQRPreconditioner> svd_f(f, Eigen::ComputeFullU | Eigen::ComputeFullV);
            auto sigmas = svd_f.singularValues();
            auto u = svd_f.matrixU();
            auto v = svd_f.matrixV();

            if (u.determinant() < 0.0)
            {
                u.col(0) *= -1.0;
                sigmas(0) *= -1.0;
            }
            if (v.determinant() < 0.0)
            {
                v.col(0) *= -1.0;
                sigmas(0) *= -1.0;
            }

            auto dPsi2_dF2 = CalcProjectedHessian(group.mMu, group.mLambda, f, u, v, sigmas);

			auto dfdx = Calc_dfdx(dF_dx, dPsi2_dF2, tet.mRestVolume);

            //std::cout << "dfdx" << std::endl;
            //std::cout << dfdx << std::endl;

            // Contribution to globalA
            // 12 x 12 matrix df/dx
            // the force on one of the 4 vertices is a 3-vector
            // the derivative of that with respect to the spatial coordinates of all four node spatial coordinates should give me a 3x12 matrix (?)
            // then we have that for each of the nodes
            // and we get a 12x12
            // so a single element tells me how the force in coordinate (x,y,z) on node i changes with respect to the (x,y,z) position of node j, makes sense
            // and a 3x3 block tells me how the force on node i changes with respect to the position of node j
            // so it's probably easier to consider the global matrix as n x n where elements are 3x3 rather than 3n x 3n with scalar elements.
            for (size_t i = 0; i < 4; i++)
            {
                auto vert_i = tet.mIndices[i];

                for (size_t j = 0; j < 4; j++)
                {
                    auto vert_j = tet.mIndices[j];

                    // we want to set a 3x3 block in globalA here from (3*vertId to 3*vertId + 2, 3*otherVertId to 3*otherVertId + 2)
                    // and we should draw from (3*i to 3*i + 2, 3*j to 3*j + 2)
                    for (size_t sub_i = 0; sub_i < 3; sub_i++)
                    {
                        for (size_t sub_j = 0; sub_j < 3; sub_j++)
                        {
                            globalA.coeffRef(3 * vert_i + sub_i, 3 * vert_j + sub_j) -= timestep * timestep * dfdx(3 * i + sub_i, 3 * j + sub_j);
                            globalDfDx.coeffRef(3 * vert_i + sub_i, 3 * vert_j + sub_j) += dfdx(3 * i + sub_i, 3 * j + sub_j);
                        }
                    }
                }
            }
            
            for (size_t i = 0; i < 4; i++)
            {
                nodeVelocities(3 * i + 0) = group.mVertices[tet.mIndices[i]].mVelocity.x;
                nodeVelocities(3 * i + 1) = group.mVertices[tet.mIndices[i]].mVelocity.y;
                nodeVelocities(3 * i + 2) = group.mVertices[tet.mIndices[i]].mVelocity.z;
            }

            localB = timestep * force + timestep * timestep * dfdx * nodeVelocities;

            if (isnan(localB.norm()))
                throw std::exception("nan detected in local forces");

            // Contribution to globalB
            for (size_t i = 0; i < 4; i++)
            {
                auto vert_i = tet.mIndices[i];
                globalB(3 * vert_i + 0) += localB(3 * i + 0);
                globalB(3 * vert_i + 1) += localB(3 * i + 1);
                globalB(3 * vert_i + 2) += localB(3 * i + 2);

                globalForce(3 * vert_i + 0) += force(3 * i + 0);
                globalForce(3 * vert_i + 1) += force(3 * i + 1);
                globalForce(3 * vert_i + 2) += force(3 * i + 2);
            }
        }

        // Apply additional per-node forces
        for (int i = 0; i < group.mVertices.size(); i++)
        {
            // Gravity
            //globalB(3 * i + 1) += -9.8 * timestep * group.mVertices[i].mMass;

            // Collision forces
            //globalB(3 * i + 0) += group.mVertices[i].mForce.x;
            //globalB(3 * i + 1) += group.mVertices[i].mForce.y;
            //globalB(3 * i + 2) += group.mVertices[i].mForce.z;
        }

        if (isnan(globalA.norm()))
            throw std::exception("Nan detected.");

        if (isnan(globalB.norm()))
            throw std::exception("Nan detected.");

        //{
        //    auto denseGlobalA = globalA.toDense();

        //    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> solver(denseGlobalA);
        //    auto eigs = solver.eigenvalues();
        //    auto minVal = 1000.0f;
        //    for (int i = 0; i < eigs.size(); i++)
        //        minVal = std::min(minVal, eigs[i]);

        //    if (minVal < 0)
        //    {
        //        std::cout << "Timestep : " << timestep << std::endl;
        //        std::cout << "Negative eig" << std::endl;
        //    }
        //}

        //std::cout << "globalA" << std::endl;
        //std::cout << globalA << std::endl;

        //std::cout << "globalB" << std::endl;
        //std::cout << globalB << std::endl;

        // PCG Solver
        {
            Eigen::ConjugateGradient<Eigen::SparseMatrix<float>, Eigen::Lower | Eigen::Upper> solver;
            solver.compute(globalA);

            globalX = solver.solve(globalB);
            if (isnan(solver.error()))
                throw std::exception("Solver failed.");

            //std::cout << "globalX" << std::endl;
            //std::cout << globalX << std::endl;
            //std::cout << "globalX inf norm: " << globalX.lpNorm<Eigen::Infinity>() << std::endl;
            //if (globalX.lpNorm<Eigen::Infinity>() > deltaVThreshold)
            //    return false;
            // lets see, make the threshold value here 1

            std::string solverInfo;
            if (solver.info() == 0)
                solverInfo = std::string("Success");
            else if (solver.info() == 1)
                solverInfo = std::string("Numerical issue");
            else if (solver.info() == 2)
                solverInfo = std::string("No convergence");
            else if (solver.info() == 3)
                solverInfo = std::string("Invalid input");

            if (solver.info() != 0)
				std::cout << "Solver info : " << solverInfo << std::endl;
        }

        // Check error against nonlinear equation
        {
            for (int i = 0; i < numVertices; i++)
            {
                globalVelocities(3 * i + 0) = group.mVertices[i].mVelocity[0];
                globalVelocities(3 * i + 1) = group.mVertices[i].mVelocity[1];
                globalVelocities(3 * i + 2) = group.mVertices[i].mVelocity[2];
                massMatrix.coeffRef(3 * i, 3 * i) = 1.0/group.mVertices[i].mMass;
                massMatrix.coeffRef(3 * i + 1, 3 * i + 1) = 1.0/group.mVertices[i].mMass;
                massMatrix.coeffRef(3 * i + 2, 3 * i + 2) = 1.0/group.mVertices[i].mMass;
            }

            // dv = h M^-1 (f0 + df/dx * dx)
            auto dx = timestep * (globalVelocities + globalX);
            auto perturbedForce = ComputePerturbedForce(group, dx);
            auto residual = timestep * massMatrix * perturbedForce - globalX;

            std::cout << "timestep: " << timestep << std::endl;
            std::cout << "infinity error: " << residual.lpNorm<Eigen::Infinity>() << std::endl;
            std::cout << "l2 error: " << residual.norm() << std::endl;
            if (residual.lpNorm<Eigen::Infinity>() > 10.0)
            {
                return false;
            }
        }

        // SparseLU Solver
        //{
        //    Eigen::SparseLU<Eigen::SparseMatrix<float>, Eigen::COLAMDOrdering<int> >   solver;
        //    solver.analyzePattern(globalA);
        //    solver.factorize(globalA);
        //    globalX = solver.solve(globalB);

        //    if (isnan(globalX.norm()))
        //        throw std::exception("Solver failed.");
        //}

        // Integrate
        for (int i = 0; i < group.mVertices.size(); i++)
        {
            auto& vertex = group.mVertices[i];
            auto deltaV = glm::dvec3(globalX(3 * i + 0), globalX(3 * i + 1), globalX(3 * i + 2));

            if (saveFrame)
            {
                auto force = glm::dvec3(globalForce(3 * i + 0), globalForce(3 * i + 1), globalForce(3 * i + 2));
                auto& vert = *frame->mutable_vertices(i);
                *vert.mutable_force() = ProtoConverter::Convert(force);
            }

            vertex.mVelocity += deltaV;
            vertex.mPosition += (double)timestep * vertex.mVelocity;
        }

        return true;

        // things to try:
        // - 'solveWithGuess' and providing the previous timestep's change in velocities as the guess?
        // - other preconditioners, solvers
        // - the doc's have globalX = solver.solve(globalB) twice, so maybe it tries for a specified number of
        //   iterations each time and refines the solution.. ?
    }
}

////////////////////////////////////////////////////////////////////////////////