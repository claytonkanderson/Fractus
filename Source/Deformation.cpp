#include "Deformation.hpp"
#include "FractureUtil.h"
#include "FractureContext.h"
#include "ConvexIntersection.h"
#include "ProtoConverter.hpp"
#include "ImplicitIntegrator.h"

#include <Mathematics/Delaunay3.h>
#include <Mathematics/SymmetricEigensolver3x3.h>
#include <Mathematics/DistPointHyperplane.h>
#include <Mathematics/IntrLine3Plane3.h>

#include <iostream>
#include <fstream>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/mat4x3.hpp>
#include <glm/gtx/hash.hpp>

using namespace glm;

namespace Deformation
{
    ////////////////////////////////////////////////////////////////////////////////

    void TetraGroup::Initialize(double lambda, double psi, double phi, double mu, double density, double toughness)
    {
        mLambda = lambda;
        mPsi = psi;
        mPhi = phi;
        mMu = mu;
        mDensity = density;
        mToughness = toughness;

        mPoissonRatio = mLambda / (2 * (mLambda + mMu));

        mSummary.set_lambda(mLambda);
        mSummary.set_psi(mPsi);
        mSummary.set_phi(mPhi);
        mSummary.set_mu(mMu);
        mSummary.set_density(mDensity);
        mSummary.set_toughness(mToughness);
    }

    ////////////////////////////////////////////////////////////////////////////////

    void TetraGroup::Update(double timestep)
	{
        auto frame = mSummary.add_frames();

        bool saveFrame = false;
        if (saveFrame)
            frame->set_time(mSimulationTime);

        ClearState(saveFrame, frame);
        // Apply collision forces from the current state
        ConvexIntersection::ResolveCollisions(mVertices, mIdToTetrahedra);
        // Calculate new positions and velocities based on the deformed position
        // - also includes gravity
        // - also includes collision forces
        const float deltaVThreshold = 70;
        const int maxNumTries = 20;

        for (int i = 0; i < maxNumTries; i++)
        {
            if (Deformation::ImplicitUpdate(*this, timestep, saveFrame, frame, deltaVThreshold))
                break;

            timestep /= 2.0f;

            if (i == maxNumTries - 1)
                throw std::exception("Simulation unstable despite small timestep.");
        }

        // Calculate separation tensor based on new deformed state
        CalculateSeparationTensor(saveFrame, frame);
        // Apply fracture based on separation tensor
        //Fracture();
        // Apply ground response
        ApplyGroundCollision();

        mStepNum++;
        mSimulationTime += timestep;

        if (saveFrame)
            mStepsSinceLastSave = 0;
        else
            mStepsSinceLastSave++;
	}

    ////////////////////////////////////////////////////////////////////////////////

    void TetraGroup::ComputeDerivedQuantities()
    {
        for (auto & vert : mVertices)
        {
            vert.mMass = 0.0f;
            vert.mInvMass = 0.0f;
        }

        for (auto & pair : mIdToTetrahedra)
        {
            auto& tetrahedra = pair.second;
            auto& m0 = mVertices[tetrahedra.mIndices[0]].mMaterialCoordinates;
            auto& m1 = mVertices[tetrahedra.mIndices[1]].mMaterialCoordinates;
            auto& m2 = mVertices[tetrahedra.mIndices[2]].mMaterialCoordinates;
            auto& m3 = mVertices[tetrahedra.mIndices[3]].mMaterialCoordinates;

            glm::dmat4 m = glm::dmat4(
                glm::dvec4(m0, 1.f),
                glm::dvec4(m1, 1.f),
                glm::dvec4(m2, 1.f),
                glm::dvec4(m3, 1.f)
            );

            tetrahedra.mBeta = glm::inverse(m);

            tetrahedra.mVolume = 1.0f / 6.0f * fabs(
                glm::dot(
                    glm::cross(
                        m1 - m0,
                        m2 - m0
                    ),
                    m3 - m0));

            if (tetrahedra.mVolume <= 0.0f)
            {
                throw std::exception("Zero volume tetrahedra encountered.");
                std::cout << "Tetrahedra volume " << tetrahedra.mVolume << std::endl;
            }

            if (tetrahedra.mRestVolume == -1)
                tetrahedra.mRestVolume = tetrahedra.mVolume;

            tetrahedra.mMass = mDensity * tetrahedra.mVolume;

            for (int j = 0; j < 4; j++)
            {
                mVertices[tetrahedra.mIndices[j]].mMass += 0.25f * tetrahedra.mMass;
            }
        }

        for (auto & vert : mVertices)
        {
            vert.mInvMass = 1.0f / vert.mMass;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////

    void TetraGroup::FractureNode(size_t fractureNodeIdx, const glm::dvec3& fracturePlaneNormal)
    {
        std::cout << "Fracture Node " << fractureNodeIdx << std::endl;
        FractureContext context(fracturePlaneNormal, fractureNodeIdx, mIdToTetrahedra, mVertices, mTetIdCounter);
        context.Fracture();
    }
    
    ////////////////////////////////////////////////////////////////////////////////

    size_t TetraGroup::SplitEdge(const glm::ivec2 & edgeIdx, size_t fractureNodeIdx, const glm::dvec3& planeNormal)
    {
        const auto& planePosition = mVertices[fractureNodeIdx].mPosition;
        const auto& edgePos0 = mVertices[edgeIdx.x].mPosition;
        const auto& edgePos1 = mVertices[edgeIdx.y].mPosition;

        glm::dvec3 intersectionPos;
        double d;
        if (!PlaneIntersectEdge(planePosition, planeNormal, edgePos0, edgePos1, d, &intersectionPos))
            std::cout << "Expected edge to be split by plane but the two don't intersect." << std::endl;

        Vertex vertex;
        vertex.mPosition = intersectionPos;

        const auto& m0 = mVertices[edgeIdx.x].mMaterialCoordinates;
        const auto& m1 = mVertices[edgeIdx.y].mMaterialCoordinates;
        vertex.mMaterialCoordinates = m0 + (m1 - m0) * d;

        const auto& v0 = mVertices[edgeIdx.x].mVelocity;
        const auto& v1 = mVertices[edgeIdx.y].mVelocity;
        // Not sure this is the right way to do this but it seems rather reasonable
        vertex.mVelocity = v0 + (v1 - v0) * d;

        mVertices.push_back(vertex);

        return mVertices.size() - 1;
    }

    ////////////////////////////////////////////////////////////////////////////////

    void TetraGroup::OutputSaveFile()
    {
        std::ofstream ofs("simulation.summary", std::ios_base::out | std::ios_base::binary);
        IronGames::SimulationSummaries summaries;
        *summaries.add_summaries() = mSummary;
        summaries.SerializeToOstream(&ofs);
    }

    ////////////////////////////////////////////////////////////////////////////////

    void TetraGroup::ClearState(bool saveFrame, IronGames::SimulationFrame* frame)
    {
        for (auto& vertex : mVertices)
        {
            vertex.mCompressiveForces.clear();
            vertex.mTensileForces.clear();
            vertex.mCollisionForces.clear();
            vertex.mForce = vec3(0.0);
            vertex.mLargestEigenvalue = 0.0f;
            vertex.mPrincipalEigenVector = vec3(0.0f);

            if (saveFrame)
            {
                auto vert = frame->add_vertices();
                *vert->mutable_position() = ProtoConverter::Convert(vertex.mPosition);
                *vert->mutable_material_coordinates() = ProtoConverter::Convert(vertex.mMaterialCoordinates);
                *vert->mutable_velocity() = ProtoConverter::Convert(vertex.mVelocity);
                vert->set_mass(vertex.mMass);
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////////////

    void TetraGroup::CalculateSeparationTensor(bool saveFrame, IronGames::SimulationFrame* frame)
    {
        dmat4x3 pMat;
        dmat4x3 vMat;
        dmat3 sigmaPlus = dmat3(0);
        dmat3 sigmaMinus = dmat3(0);
        dmat3 mMat;
        dmat3 sigma;
        std::array<double, 3> eigenvalues;
        std::array<std::array<double, 3>, 3> eigenvectors;
        double separation[3][3];
        const int32_t sortType = -1; // -1 is decreasing order, so the first is the largest
        const bool aggressive = false;
        gte::SymmetricEigensolver3x3<double> solver;

        for (auto& pair : mIdToTetrahedra)
        {
            size_t tetIdx = pair.first;
            auto& tetrahedra = pair.second;
            auto& v0 = mVertices[tetrahedra.mIndices[0]];
            auto& v1 = mVertices[tetrahedra.mIndices[1]];
            auto& v2 = mVertices[tetrahedra.mIndices[2]];
            auto& v3 = mVertices[tetrahedra.mIndices[3]];

            pMat = dmat4x3(v0.mPosition,
                v1.mPosition,
                v2.mPosition,
                v3.mPosition);

            vMat = dmat4x3(v0.mVelocity,
                v1.mVelocity,
                v2.mVelocity,
                v3.mVelocity);

            dmat3x4 pBeta = pMat * tetrahedra.mBeta;

            dvec3 dx_u1 = pBeta * dvec4(1, 0, 0, 0);
            dvec3 dx_u2 = pBeta * dvec4(0, 1, 0, 0);
            dvec3 dx_u3 = pBeta * dvec4(0, 0, 1, 0);

            std::array<dvec3, 3> dx{ dx_u1, dx_u2, dx_u3 };

            dmat3x4 vBeta = vMat * tetrahedra.mBeta;
            dvec3 dxd_u1 = vBeta * dvec4(1, 0, 0, 0);
            dvec3 dxd_u2 = vBeta * dvec4(0, 1, 0, 0);
            dvec3 dxd_u3 = vBeta * dvec4(0, 0, 1, 0);

            std::array<dvec3, 3 > dxd{ dxd_u1, dxd_u2, dxd_u3 };

            dmat3 strainTensor;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    strainTensor[i][j] = dot(dx[i], dx[j]) - (i == j ? 1 : 0);

            dmat3 rateOfStrainTensor;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    rateOfStrainTensor[i][j] = dot(dx[i], dxd[j]) + dot(dxd[i], dx[j]);

            dmat3 elasticStress(0.0);
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                        elasticStress[i][j] += mLambda * strainTensor[k][k] * (i == j ? 1 : 0) + 2 * mMu * strainTensor[i][j];

            dmat3 viscousStress(0.0);
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                        viscousStress[i][j] += mPhi * rateOfStrainTensor[k][k] * (i == j ? 1 : 0) + 2 * mPsi * rateOfStrainTensor[i][j];

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    sigma[i][j] = (double)elasticStress[i][j] + (double)viscousStress[i][j];

                    if (isnan(sigma[i][j]))
                        throw std::exception("Nan value in sigma detected.");
                }
            }

            for (int i = 0; i < 4; i++)
            {
                dvec3 forceOnNode = dvec3(0);
                for (int j = 0; j < 4; j++)
                {
                    double innerProduct = 0;
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            // Reversed indices from paper because glm is column, row indexing
                            // and the conventional matrix indexing is (row, column)
                            // - should not matter in a lot of places above because the matrices are symmetric
                            innerProduct += tetrahedra.mBeta[l][j] * tetrahedra.mBeta[k][i] * sigma[l][k];
                        }
                    }
                    forceOnNode += pMat[j] * innerProduct;
                }

                forceOnNode *= -tetrahedra.mVolume * 0.5f;

                if (isnan(forceOnNode.x) || isnan(forceOnNode.y) || isnan(forceOnNode.z))
                {
                    throw std::exception("Nan value in forceOnNode detected.");
                }

                mVertices[tetrahedra.mIndices[i]].mForce += forceOnNode;
            }

            solver(sigma[0][0], sigma[0][1], sigma[0][2],
                sigma[1][1], sigma[1][2], sigma[2][2],
                aggressive, sortType, eigenvalues, eigenvectors);

            sigmaPlus = dmat3(0.0);
            sigmaMinus = dmat3(0.0);

            for (int i = 0; i < 3; i++)
            {
                dvec3 eigenVec(eigenvectors[i][0], eigenvectors[i][1], eigenvectors[i][2]);
                ComputeMMat(eigenVec, mMat);

                sigmaPlus += fmax(0.0f, eigenvalues[i]) * mMat;
                sigmaMinus += fmin(0.0f, eigenvalues[i]) * mMat;
            }

            std::array<dvec3, 4> fPlus = { dvec3(0.0),dvec3(0.0),dvec3(0.0),dvec3(0.0) };
            std::array<dvec3, 4> fMinus = { dvec3(0.0),dvec3(0.0),dvec3(0.0),dvec3(0.0) };;

            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    double innerProductPlus = 0;
                    double innerProductMinus = 0;
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            innerProductPlus += tetrahedra.mBeta[k][i] * tetrahedra.mBeta[l][j] * sigmaPlus[l][k];
                            innerProductMinus += tetrahedra.mBeta[k][i] * tetrahedra.mBeta[l][j] * sigmaMinus[l][k];
                        }
                    }
                    fPlus[i] += pMat[j] * innerProductPlus;
                    fMinus[i] += pMat[j] * innerProductMinus;
                }

                fPlus[i] *= -tetrahedra.mVolume * 0.5f;
                fMinus[i] *= -tetrahedra.mVolume * 0.5f;
                //fPlus[i] *= tetrahedra.mVolume * 0.5f;
                //fMinus[i] *= tetrahedra.mVolume * 0.5f;
            }

            for (int i = 0; i < 4; i++)
            {
                mVertices[tetrahedra.mIndices[i]].mCompressiveForces.push_back(fMinus[i]);
                mVertices[tetrahedra.mIndices[i]].mTensileForces.push_back(fPlus[i]);
                //mVertices[tetrahedra.mIndices[i]].mCompressiveForces.push_back(fPlus[i]);
                //mVertices[tetrahedra.mIndices[i]].mTensileForces.push_back(fMinus[i]);
            }

            if (saveFrame)
            {
                auto& tet = *frame->add_tetrahedra();
                tet.set_mass(tetrahedra.mMass);
                tet.set_volume(tetrahedra.mVolume);
                *tet.mutable_strain_tensor() = ProtoConverter::Convert(strainTensor);
                *tet.mutable_stress_tensor() = ProtoConverter::Convert(sigma);
                for (auto& idx : tetrahedra.mIndices)
                    tet.add_indices(idx);
            }
        }

        dvec3 a;
        dmat3 mSeparation;
        dmat3 m_Mat;
        dvec3 compressiveForceSum;
        dvec3 tensileForceSum;

        for (size_t vertexIdx = 0; vertexIdx < mVertices.size(); vertexIdx++)
        {
            auto& vertex = mVertices[vertexIdx];
            //vertex.mForce += glm::dvec3(0, -9.8, 0) * vertex.mMass;

            //a = vertex.mInvMass * vertex.mForce;

            //if (isnan(a[0]) || isnan(a[1]) || isnan(a[2]))
            //    throw std::exception("Nan found in acceleration.");

            //vertex.mVelocity += a * timestep;
            //vertex.mPosition += vertex.mVelocity * timestep;

            compressiveForceSum = vec3(0.0);
            tensileForceSum = vec3(0.0);

            // Tensile is f+
            // Compressive is f-
            // Formula is -m(f+) + m(f-) + sum(m(f+)) - sum(m(f-))

            for (const auto& compressiveForce : vertex.mCompressiveForces)
                compressiveForceSum += compressiveForce;
            for (const auto& tensileForce : vertex.mTensileForces)
                tensileForceSum += tensileForce;

            mSeparation = dmat3(0.0);

            ComputeMMat(tensileForceSum, m_Mat);
            mSeparation -= m_Mat;
            ComputeMMat(compressiveForceSum, m_Mat);
            mSeparation += m_Mat;

            for (const auto& compressiveForce : vertex.mCompressiveForces)
            {
                ComputeMMat(compressiveForce, m_Mat);
                mSeparation -= m_Mat;
            }
            for (const auto& tensileForce : vertex.mTensileForces)
            {
                ComputeMMat(tensileForce, m_Mat);
                mSeparation += m_Mat;
            }

            mSeparation *= 0.5f;

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    separation[i][j] = mSeparation[i][j];
                }
            }

            solver(separation[0][0], separation[0][1], separation[0][2],
                separation[1][1], separation[1][2], separation[2][2],
                aggressive, sortType, eigenvalues, eigenvectors);

            // Look for largest eigenvalue and corresponding eigen vector
            double largestEigenvalue = eigenvalues[0];
            dvec3 principalEigenVector = dvec3(eigenvectors[0][0], eigenvectors[0][1], eigenvectors[0][2]);

            mVertices[vertexIdx].mPrincipalEigenVector = principalEigenVector;
            mVertices[vertexIdx].mLargestEigenvalue = largestEigenvalue;

            if (saveFrame)
            {
                auto& vert = *frame->mutable_vertices(vertexIdx);
                //*vert.mutable_force() = ProtoConverter::Convert(mVertices[vertexIdx].mForce);
                vert.set_largest_eigenvalue(mVertices[vertexIdx].mLargestEigenvalue);
                *vert.mutable_principal_eigenvector() = ProtoConverter::Convert(mVertices[vertexIdx].mPrincipalEigenVector);

                for (const auto& compressiveForce : mVertices[vertexIdx].mCompressiveForces)
                    *vert.add_compressive_forces() = ProtoConverter::Convert(compressiveForce);
                for (const auto& tensileForce : mVertices[vertexIdx].mTensileForces)
                    *vert.add_tensile_forces() = ProtoConverter::Convert(tensileForce);
                for (const auto& collisionForce : mVertices[vertexIdx].mCollisionForces)
                    *vert.add_collision_forces() = ProtoConverter::Convert(collisionForce);

                *vert.mutable_separation_tensor() = ProtoConverter::Convert(mSeparation);
            }
        }
    }

    void TetraGroup::Fracture()
    {
        // Casually prevent fracturing newly created vertices
        size_t numVertices = mVertices.size();
        if (numVertices < mMaxNumVertices)
        {
            //for (size_t idx = 0; idx < numVertices; idx++)
            //{
            //    if (mVertices[idx].mLargestEigenvalue > mToughness)
            //        FractureNode(idx, mVertices[idx].mPrincipalEigenVector);

            //    if (mVertices.size() >= mMaxNumVertices)
            //        break;
            //}

            // let's try fracturing only the highest eigenvalue and only allowing one fracture per frame
            struct FractureAttempt
            {
                size_t mVertexId = 0;
                double mLargestEigenValue = 0;
            };

            std::vector<FractureAttempt> fractureAttempts;

            for (size_t idx = 0; idx < numVertices; idx++)
            {
                if (mVertices[idx].mLargestEigenvalue > mToughness)
                {
                    fractureAttempts.push_back({ idx, mVertices[idx].mLargestEigenvalue });
                }
            }

            std::sort(fractureAttempts.begin(), fractureAttempts.end(), [](const FractureAttempt& a, const FractureAttempt& b)
                {
                    return a.mLargestEigenValue > b.mLargestEigenValue;
                });

            bool fractured = fractureAttempts.empty();
            for (const auto& attempt : fractureAttempts)
            {
                //std::cout << "Fracture Node " << attempt.mVertexId << std::endl;
                FractureContext context(mVertices[attempt.mVertexId].mPrincipalEigenVector, attempt.mVertexId, mIdToTetrahedra, mVertices, mTetIdCounter);
                if (context.Fracture())
                {
                    fractured = true;
                    break;
                }
            }

            if (!fractured)
                std::cout << "All fracture attempts failed so no fracture occurred." << std::endl;

            ComputeDerivedQuantities();
        }
    }

    ////////////////////////////////////////////////////////////////////////////////

    void TetraGroup::ApplyGroundCollision()
    {
        for (auto& vertex : mVertices)
        {
            if (vertex.mPosition.y < 0)
            {
                vertex.mPosition.y = 0.0;
                double elasticity = 0.4f;
                double friction = 0.1f;
                const auto& velocity = vertex.mVelocity;
                vertex.mVelocity = (glm::dvec3((1 - friction) * velocity.x, -elasticity * velocity.y, (1 - friction) * velocity.z));
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////////////

    std::vector<size_t> TetraGroup::GetTetrahedrasFromNode(size_t nodeIdx) const
    {
        std::vector<size_t> indices;

        for (const auto & pair : mIdToTetrahedra)
        {
            size_t tetIdx = pair.first;
            const auto& tet = pair.second;
            if (tet.ContainsVertexIndex(nodeIdx))
                indices.push_back(tetIdx);
        }

        return indices;
    }

    ////////////////////////////////////////////////////////////////////////////////

    bool TetraGroup::EdgeIntersectsPlane(const glm::ivec2& edgeIdx, int fractureNodeIdx, const glm::dvec3& planeNormal) const
    {
        const auto& planePosition = mVertices[fractureNodeIdx].mPosition;
        const auto& edgePos0 = mVertices[edgeIdx.x].mPosition;
        const auto& edgePos1 = mVertices[edgeIdx.y].mPosition;
        double d;
        return PlaneIntersectEdge(planePosition, planeNormal, edgePos0, edgePos1, d);
    }

    ////////////////////////////////////////////////////////////////////////////////

    bool TetraGroup::PlaneIntersectEdge(const glm::dvec3& planePos, const glm::dvec3& planeNormal, const glm::dvec3& edgePos0, const glm::dvec3& edgePos1, double& d, glm::dvec3* intersectionPos) const
    {
        // plane equation is (p - p0) * n = 0, where n is the normal vector
        // line equation is p = l0 + v * t, where v is the direction vector of the line
        // compute d = (p0 - l0) * n / (l * n)
        // intersection is at p = l0 + v * d
        // if d == 0, plane contains line
        auto edgeDirVec = glm::normalize(edgePos1 - edgePos0);

        d = glm::dot((planePos - edgePos0), planeNormal) / glm::dot(edgeDirVec, planeNormal);
        
        // Outside the line segment
        if (d <= 0 || d > 1)
            return false;

        if (intersectionPos)
            *intersectionPos = edgePos0 + d * edgeDirVec;

        return true;
    }

    ////////////////////////////////////////////////////////////////////////////////

    double TetraGroup::GetSignedDistanceToPlane(int nodeIdx, int fractureNodeIdx, const glm::dvec3& planeNormal) const
    {
        const auto& planePosition = mVertices[fractureNodeIdx].mPosition;
        const auto& nodePosition = mVertices[nodeIdx].mPosition;
        return glm::dot(planeNormal, nodePosition - planePosition);
    }
    
    ////////////////////////////////////////////////////////////////////////////////

    std::vector<size_t> TetraGroup::GetTetrahedraNeighbors(size_t tetrahedraIdx) const
    {
        std::vector<size_t> neighborIds;
        const auto& tet = mIdToTetrahedra.at(tetrahedraIdx);

        for (const auto & pair : mIdToTetrahedra)
        {
            size_t tetIdx = pair.first;
            if (tetIdx == tetrahedraIdx)
                continue;

            const auto& neighborCandidate = pair.second;
            for (const auto& nodeIdx : neighborCandidate.mIndices)
            {
                if (tet.ContainsVertexIndex(nodeIdx))
                    neighborIds.push_back(tetIdx);
            }
        }

        return neighborIds;
    }

    ////////////////////////////////////////////////////////////////////////////////

    double TetraGroup::GetMinVertexMass() const
    {
        double val = DBL_MAX;
        for (const auto& vert : mVertices)
            val = std::min(vert.mMass, val);
        return val;
    }
    
    ////////////////////////////////////////////////////////////////////////////////

    double TetraGroup::GetMaxVertexMass() const
    {
        double val = DBL_MIN;
        for (const auto& vert : mVertices)
            val = std::max(vert.mMass, val);
        return val;
    }

    ////////////////////////////////////////////////////////////////////////////////

    double TetraGroup::GetTotalMass() const
    {
        return mDensity * GetTotalVolume();
    }

    ////////////////////////////////////////////////////////////////////////////////

    double TetraGroup::GetTotalVolume() const
    {
        double vol = 0;
        for (const auto& pair : mIdToTetrahedra)
        {
            vol += pair.second.mVolume;
        }

        return vol;
    }

    ////////////////////////////////////////////////////////////////////////////////

    double TetraGroup::GetAverageMaxEigenvalue() const
    {
        double avg = 0;
        for (const auto& vert : mVertices)
            avg += vert.mLargestEigenvalue;

        avg /= mVertices.size();
        return avg;
    }

    ////////////////////////////////////////////////////////////////////////////////

    double TetraGroup::GetMaxEigenvalue() const
    {
        double max = 0;
        for (const auto& vert : mVertices)
            max = std::max(max, vert.mLargestEigenvalue);

        return max;
    }

    ////////////////////////////////////////////////////////////////////////////////

    std::array<size_t, 3> Tetrahedra::GetOtherVertices(size_t vertexId) const
    {
        if (!ContainsVertexIndex(vertexId))
            throw std::exception("Tetrahedron for Get Other Vertices did not contain the specified vertex.");

        if (mIndices[0] == vertexId)
            return { mIndices[1], mIndices[2], mIndices[3] };

        if (mIndices[1] == vertexId)
            return { mIndices[0], mIndices[2], mIndices[3] };

        if (mIndices[2] == vertexId)
            return { mIndices[1], mIndices[0], mIndices[3] };

        if (mIndices[3] == vertexId)
            return { mIndices[0], mIndices[1], mIndices[2] };

        throw std::exception("Tetrahedron for Get Other Vertices encountered an unexpected error.");
        return std::array<size_t, 3>();
    }

    ////////////////////////////////////////////////////////////////////////////////

    std::array<size_t, 2> Tetrahedra::GetOtherVertices(size_t vertexId1, size_t vertexId2) const
    {
        std::array<size_t, 2> otherVertices{};
        size_t counter = 0;
        for (size_t i = 0; i < 4; i++)
        {
            if (mIndices[i] == vertexId1 || mIndices[i] == vertexId2)
                continue;

            otherVertices[counter] = mIndices[i];
            counter++;
        }

        if (counter != 2)
            throw std::exception("Tetrahedron for Get Other Vertices (2) encountered an unexpected error.");

        return otherVertices;
    }

    ////////////////////////////////////////////////////////////////////////////////

    size_t Tetrahedra::GetOtherVertex(const std::array<size_t, 3>& vertices) const
    {
        for (size_t i = 0; i < 4; i++)
        {
            if (mIndices[i] == vertices[0] || mIndices[i] == vertices[1] || mIndices[i] == vertices[2])
                continue;

            return mIndices[i];
        }

        throw std::exception("Tetrahedron for Get Other Vertex encountered an unexpected error.");
        return size_t();
    }
    
    ////////////////////////////////////////////////////////////////////////////////

    bool Tetrahedra::ContainsVertexIndex(size_t idx) const
    {
        for (const auto& index : mIndices)
        {
            if (idx == index)
                return true;
        }

        return false;
    }

    ////////////////////////////////////////////////////////////////////////////////

    bool Tetrahedra::ContainsEdgeIndex(const glm::ivec2& edgeId) const
    {
        const auto& edges = GetEdges();
        for (const auto& myEdge : edges)
        {
            if (myEdge == edgeId)
                return true;
        }

        return false;
    }

    ////////////////////////////////////////////////////////////////////////////////

    bool Tetrahedra::ContainsFaceIndex(const glm::ivec3& faceId) const
    {
        // every vertex is in 3 triangles
        // there are four total triangles
        // v0, v1, v2
        // v0, v1, v3
        // v0, v2, v3
        // v1, v2, v3

        if (faceId == GetFaceId(mIndices[0], mIndices[1], mIndices[2]))
            return true;
        if (faceId == GetFaceId(mIndices[0], mIndices[1], mIndices[3]))
            return true;
        if (faceId == GetFaceId(mIndices[0], mIndices[2], mIndices[3]))
            return true;
        if (faceId == GetFaceId(mIndices[1], mIndices[2], mIndices[3]))
            return true;

        return false;
    }

    ////////////////////////////////////////////////////////////////////////////////

    std::array<glm::ivec2, 6> Tetrahedra::GetEdges() const
    {
        std::array<glm::ivec2, 6> edges;
        edges[0] = GetEdgeId(mIndices[0], mIndices[1]);
        edges[1] = GetEdgeId(mIndices[0], mIndices[2]);
        edges[2] = GetEdgeId(mIndices[0], mIndices[3]);
        edges[3] = GetEdgeId(mIndices[1], mIndices[2]);
        edges[4] = GetEdgeId(mIndices[1], mIndices[3]);
        edges[5] = GetEdgeId(mIndices[2], mIndices[3]);

        return edges;
    }

    ////////////////////////////////////////////////////////////////////////////////
    
    void Tetrahedra::ReplaceVertex(size_t oldVertexId, size_t newVertexId)
    {
        bool setVertex = false;
        for (int i = 0; i < 4; i++)
        {
            if (mIndices[i] == oldVertexId)
            {
                mIndices[i] = newVertexId;
                setVertex = true;
                break;
            }
        }

        if (!setVertex)
            std::cout << "Attempted to replace a vertex but the specified was not in the tetrahedron." << std::endl;
    }

    ////////////////////////////////////////////////////////////////////////////////

    glm::dvec3 Tetrahedra::GetCentroid(const std::vector<Vertex>& vertices) const
    {
        glm::dvec3 center(0);
        for (auto id : mIndices)
            center += 0.25 * vertices[id].mPosition;

        return center;
    }

    ////////////////////////////////////////////////////////////////////////////////

    double Tetrahedra::GetMinDihedralAngle(const std::vector<Vertex>& vertices) const
    {
        // Angle between every pair of edges that shares a node
        // All Edges:
        // 01 = 0
        // 02 = 1
        // 03 = 2
        // 12 = 3
        // 13 = 4
        // 23 = 5
        //         edges[0] = GetEdgeId(mIndices[0], mIndices[1]);
        //edges[1] = GetEdgeId(mIndices[0], mIndices[2]);
        //edges[2] = GetEdgeId(mIndices[0], mIndices[3]);
        //edges[3] = GetEdgeId(mIndices[1], mIndices[2]);
        //edges[4] = GetEdgeId(mIndices[1], mIndices[3]);
        //edges[5] = GetEdgeId(mIndices[2], mIndices[3]);
        // 
        // Edges that share a node:
        // 01, 02 | 0, 1
        // 01, 03 | 0, 2
        // 01, 12 | 0, 3
        // 01, 13 | 0, 4
        // 02, 12 | 1, 3
        // 02, 23 | 1, 5
        // 03, 13 | 2, 4
        // 03, 23 | 2, 5
        // 12, 23 | 3, 5
        // 13, 23 | 4, 5

        auto angleFunc = [&](const auto & e1, const auto & e2) 
        {
            auto dir1 = vertices[e1.x].mPosition - vertices[e1.y].mPosition;
            auto dir2 = vertices[e2.x].mPosition - vertices[e2.y].mPosition;
            auto angle = glm::acos(glm::dot(dir1, dir2) / glm::length(dir1) / glm::length(dir2));

            return glm::degrees(std::min(angle, glm::pi<float>() - angle));
        };

        std::vector<float> angles;

        const auto & edges = GetEdges();
        angles.push_back(angleFunc(edges[0], edges[1]));
        angles.push_back(angleFunc(edges[0], edges[2]));
        angles.push_back(angleFunc(edges[0], edges[3]));
        angles.push_back(angleFunc(edges[0], edges[4]));
        angles.push_back(angleFunc(edges[1], edges[3]));
        angles.push_back(angleFunc(edges[1], edges[5]));
        angles.push_back(angleFunc(edges[2], edges[4]));
        angles.push_back(angleFunc(edges[2], edges[5]));
        angles.push_back(angleFunc(edges[3], edges[5]));
        angles.push_back(angleFunc(edges[4], edges[5]));
        return *std::min_element(angles.begin(), angles.end());
    }

    ////////////////////////////////////////////////////////////////////////////////
}