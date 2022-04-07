#include "ImplicitIntegrator.h"
#include "Deformation.hpp"
#include "ProtoConverter.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <iostream>

using Eigen::MatrixXd;

namespace Deformation
{
    using FloatT = float;
    using Mat3 = Eigen::Matrix<FloatT, 3, 3>;
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
        dj_df.block(0, 0, 3, 1) = f.col(1).cross(f.col(2));
        dj_df.block(0, 1, 3, 1) = f.col(2).cross(f.col(0));
        dj_df.block(0, 2, 3, 1) = f.col(0).cross(f.col(1));
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
        f_hat(1, 0) = f_i(2);
        f_hat(2, 0) = -f_i(1);
        
        f_hat(0, 1) = -f_i(2);
        f_hat(2, 1) = -f_i(0);

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

        h_j.block(3, 0, 3, 3) = f2_hat;
        h_j.block(6, 0, 3, 3) = -f1_hat;

        h_j.block(0, 3, 3, 3) = -f2_hat;
        h_j.block(6, 3, 3, 3) = f0_hat;

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

    Mat12 Calc_dfdx(const Mat9x12& dFdx, const Mat9 vec_dPsi2_dF2, FloatT restVol)
    {
        return -restVol * dFdx.transpose() * vec_dPsi2_dF2 * dFdx;
    }

    void ImplicitUpdate(TetraGroup& group, float timestep, bool saveFrame, IronGames::SimulationFrame* frame)
    {
        size_t numVertices = group.mVertices.size();
        Eigen::SparseMatrix<float> globalA(3* numVertices, 3* numVertices);
        Eigen::VectorXf globalB = Eigen::VectorXf::Zero(3 * numVertices);
		Eigen::VectorXf globalX = Eigen::VectorXf::Zero(3 * numVertices);
        Vec12 nodeVelocities = Vec12::Zero();
        Vec12 localB = Vec12::Zero();

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
            auto dmInv = Calc_DmInv(tet, group.mVertices);
            auto f = Calc_F(ds, dmInv);

            // Stable Neo-hookean potential
			auto j = f.determinant();
			auto dj_df = Calc_dj_df(f);
			auto g_j = Reshape3x3(dj_df);
			auto dPsi_dF = group.mMu * Reshape3x3(f) + (group.mLambda * (j - 1.0f) - group.mMu) * g_j;
			auto dF_dx = Calc_dFdx(dmInv);
			auto dPsi_dx = dF_dx.transpose() * dPsi_dF;
			Vec12 force = -tet.mRestVolume * dPsi_dx; // 12x1 atm

			auto h_j = Calc_h_j(f);
			auto dPsi2_dF2 = group.mMu * Mat9::Identity() + (-group.mMu + group.mLambda * (j - 1.0f)) * h_j + group.mLambda * g_j * g_j.transpose();
			auto dfdx = Calc_dfdx(dF_dx, dPsi2_dF2, tet.mRestVolume);

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
            }
        }

        // Apply additional per-node forces
        for (int i = 0; i < group.mVertices.size(); i++)
        {
            // Gravity
            globalB(3 * i + 1) += -9.8 * group.mVertices[i].mMass;

            // Collision forces
            //globalB(3 * i + 0) += group.mVertices[i].mForce.x;
            //globalB(3 * i + 1) += group.mVertices[i].mForce.y;
            //globalB(3 * i + 2) += group.mVertices[i].mForce.z;
        }

        //std::cout << "globalA norm: " << globalA.norm() << std::endl;
        if (isnan(globalA.norm()))
            throw std::exception("Nan detected.");

        //std::cout << "globalB norm: " << globalB.norm() << std::endl;
        if (isnan(globalB.norm()))
            throw std::exception("Nan detected.");

        Eigen::ConjugateGradient<Eigen::SparseMatrix<float>, Eigen::Lower | Eigen::Upper> solver;
        solver.compute(globalA);
        globalX = solver.solve(globalB);

        //std::cout << "Solver error : " << solver.error() << std::endl;

        if (isnan(solver.error()))
            throw std::exception("Solver failed.");

        //std::cout << "GlobalX" << std::endl;
        //std::cout << globalX << std::endl;

        // Integrate
        for (int i = 0; i < group.mVertices.size(); i++)
        {
            auto& vertex = group.mVertices[i];
            auto deltaV = glm::dvec3(globalX(3 * i + 0), globalX(3 * i + 1), globalX(3 * i + 2));

            if (saveFrame)
            {
                auto& vert = *frame->mutable_vertices(i);
                *vert.mutable_force() = ProtoConverter::Convert(deltaV);
            }

            vertex.mVelocity += deltaV;
            vertex.mPosition += (double)timestep * vertex.mVelocity;
        }

        //std::cout << "solution" << std::endl;
        //std::cout << globalX << std::endl;

        // things to try:
        // - 'solveWithGuess' and providing the previous timestep's change in velocities as the guess?
        // - other preconditioners, solvers
        // - the doc's have globalX = solver.solve(globalB) twice, so maybe it tries for a specified number of
        //   iterations each time and refines the solution.. ?
    }
}

////////////////////////////////////////////////////////////////////////////////