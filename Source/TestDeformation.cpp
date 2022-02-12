#include "TestDeformation.hpp"
#include "ProtoConverter.hpp"

#include <Mathematics/Delaunay3.h>
#include <Mathematics/SymmetricEigensolver3x3.h>
#include <Mathematics/DistPointHyperplane.h>

#include <iostream>
#include <fstream>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/mat4x3.hpp>
#include <glm/gtx/hash.hpp>

using namespace glm;

namespace TestDeformation
{
    ////////////////////////////////////////////////////////////////////////////////
    
    namespace
    {
        glm::ivec2 GetEdgeId(size_t idx1, size_t idx2)
        {
            if (idx1 == idx2)
            {
                std::cout << "Attempting to get edge for the same vertice twice." << std::endl;
                return glm::ivec2();
            }

            if (idx1 < idx2)
                return glm::ivec2(idx1, idx2);
            return glm::ivec2(idx2, idx1);
        }

        glm::ivec3 GetFaceId(size_t idx1, size_t idx2, size_t idx3)
        {
            if (idx1 == idx2 || idx1 == idx3 || idx2 == idx3)
                throw std::exception("Attempted to get a face with duplicate vertex ids.");

            std::array<int, 3> verts{ idx1, idx2, idx3 };
            std::sort(verts.begin(), verts.end());

            return glm::ivec3(verts[0], verts[1], verts[2]);
        }
    }

    void ComputeMMat(const dvec3& eigenVector, dmat3& outMat)
    {
        if (glm::length(eigenVector) < 1e-6)
        {
            outMat = dmat3(0.0);
            return;
        }

        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                outMat[j][k] = eigenVector[j] * eigenVector[k];

        outMat /= glm::length(eigenVector);
    }

    ////////////////////////////////////////////////////////////////////////////////

    void TetraGroup::Initialize(double lambda, double psi, double phi, double mu, double density, double toughness)
    {
        mLambda = lambda;
        mPsi = psi;
        mPhi = phi;
        mMu = mu;
        mDensity = density;
        mToughness = toughness;

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

        //auto* frame = (mSummary ? mSummary->add_frames() : nullptr);
        auto frame = mSummary.add_frames();

        //bool saveFrame = (mStepsSinceLastSave >= mSaveEveryXSteps && mSummary);
        bool saveFrame = true;

        if (saveFrame)
            frame->set_time(mSimulationTime);

        //
        gte::SymmetricEigensolver3x3<double> solver;

        for (auto& vertex : mVertices)
        {
            vertex.mCompressiveForces.clear();
            vertex.mTensileForces.clear();
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

		for (size_t tetIdx = 0; tetIdx < mTetrahedra.size(); tetIdx++)
		{
            auto& tetrahedra = mTetrahedra[tetIdx];
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
                    for(int k = 0; k < 3; k++)
                        elasticStress[i][j] += mLambda * strainTensor[k][k] * (i == j ? 1 : 0) + 2 * mMu * strainTensor[i][j];

            dmat3 viscousStress(0.0);
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    for(int k = 0; k < 3; k++)
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

                //fPlus[i] *= -tetrahedra.mVolume * 0.5f;
                //fMinus[i] *= -tetrahedra.mVolume * 0.5f;
                fPlus[i] *= tetrahedra.mVolume * 0.5f;
                fMinus[i] *= tetrahedra.mVolume * 0.5f;
            }

            for (int i = 0; i < 4; i++)
            {
                //mVertices[tetrahedra.mIndices[i]].mCompressiveForces.push_back(fMinus[i]);
                //mVertices[tetrahedra.mIndices[i]].mTensileForces.push_back(fPlus[i]);
                mVertices[tetrahedra.mIndices[i]].mCompressiveForces.push_back(fPlus[i]);
                mVertices[tetrahedra.mIndices[i]].mTensileForces.push_back(fMinus[i]);
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
            vertex.mForce += glm::dvec3(0, -9.8, 0) * vertex.mMass;

            a = vertex.mInvMass * vertex.mForce;
            vertex.mVelocity += a * timestep;
            vertex.mPosition += vertex.mVelocity * timestep;

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

            //if (largestEigenvalue > mToughness)
            //    std::cout << " Large eigen " << std::endl;

            mVertices[vertexIdx].mPrincipalEigenVector = principalEigenVector;
            mVertices[vertexIdx].mLargestEigenvalue = largestEigenvalue;

            if (saveFrame)
            {
                auto& vert = *frame->mutable_vertices(vertexIdx);
                *vert.mutable_force() = ProtoConverter::Convert(mVertices[vertexIdx].mForce);
                vert.set_largest_eigenvalue(mVertices[vertexIdx].mLargestEigenvalue);
                *vert.mutable_principal_eigenvector() = ProtoConverter::Convert(mVertices[vertexIdx].mPrincipalEigenVector);

                for (const auto& compressiveForce : mVertices[vertexIdx].mCompressiveForces)
                    *vert.add_compressive_forces() = ProtoConverter::Convert(compressiveForce);
                for (const auto& tensileForce : mVertices[vertexIdx].mTensileForces)
                    *vert.add_tensile_forces() = ProtoConverter::Convert(tensileForce);

                *vert.mutable_separation_tensor() = ProtoConverter::Convert(mSeparation);
            }
        }

        // Casually prevent fracturing newly created vertices
        size_t numVertices = mVertices.size();
        if (numVertices < mMaxNumVertices)
        {
            for (size_t idx = 0; idx < numVertices; idx++)
            {
                if (mVertices[idx].mLargestEigenvalue > mToughness)
                    FractureNode(idx, mVertices[idx].mPrincipalEigenVector);

                if (mVertices.size() >= mMaxNumVertices)
                    break;
            }

            /*
            
                        // let's try fracturing only the highest eigenvalue and only allowing one fracture per frame

            float maxEigenValue = -1;
            size_t maxEigenValueNodeId = -1;

            for (size_t idx = 0; idx < numVertices; idx++)
            {
                if (mVertices[idx].mLargestEigenvalue > mToughness)
                {
                    if (mVertices[idx].mLargestEigenvalue > maxEigenValue)
                    {
                        maxEigenValue = mVertices[idx].mLargestEigenvalue;
                        maxEigenValueNodeId = idx;
                    }
                }

                if (mVertices.size() >= mMaxNumVertices)
                    break;
            }

            if (maxEigenValueNodeId != -1)
                FractureNode(maxEigenValueNodeId, mVertices[maxEigenValueNodeId].mPrincipalEigenVector);
            
            */

            ComputeDerivedQuantities();
        }

        for (auto& vertex : mVertices)
        {
            if (vertex.mPosition.y < 0)
            {
                vertex.mPosition.y = -vertex.mPosition.y;
                double elasticity = 0.9f;
                double friction = 0.1f;
                const auto & velocity = vertex.mVelocity;
                vertex.mVelocity = (glm::dvec3((1 - friction)* velocity.x, -elasticity * velocity.y, (1 - friction)* velocity.z));
            }
        }

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

        for (auto & tetrahedra : mTetrahedra)
        {
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
				std::cout << "Tetrahedra volume " << tetrahedra.mVolume << std::endl;

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
        FractureContext context(fracturePlaneNormal, fractureNodeIdx, mTetrahedra, mVertices);
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

    std::vector<size_t> TetraGroup::GetTetrahedrasFromNode(size_t nodeIdx) const
    {
        std::vector<size_t> indices;

        for (size_t tetIdx = 0; tetIdx < mTetrahedra.size(); tetIdx++)
        {
            if (mTetrahedra[tetIdx].ContainsVertexIndex(nodeIdx))
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
        const auto& tet = mTetrahedra[tetrahedraIdx];

        for (int i = 0; i < mTetrahedra.size(); i++)
        {
            if (i == tetrahedraIdx)
                continue;

            const auto& neighborCandidate = mTetrahedra[i];
            for (const auto& nodeIdx : neighborCandidate.mIndices)
            {
                if (tet.ContainsVertexIndex(nodeIdx))
                    neighborIds.push_back(i);
            }
        }

        return neighborIds;
    }

    ////////////////////////////////////////////////////////////////////////////////

    namespace
    {
        size_t GetCommonVertexFromEdges(const glm::ivec2& edge0, const glm::ivec2& edge1)
        {
            if (edge0.x == edge1.x || edge0.x == edge1.y)
                return edge0.x;
            if (edge0.y == edge1.x || edge0.y == edge1.y)
                return edge0.y;

            throw std::exception("Failed to find common vertices for provided edges.");
            return 0;
        }
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

    void FractureContext::Fracture()
    {
        // create new node (the 'negative' one)
        mNegativeFractureNodeIdx = CloneVertex(mFractureNodeIdx);

        bool fractured = true;

        // We loop through each tetrahedra that owns the fracture node, splitting one and all their neighbors
        // per iteration of the loop.
        while (fractured)
        {
            fractured = false;

            const auto& neighbors = GetTetrahedraNeighbors(mFractureNodeIdx);

            for (size_t tetIdx = 0; tetIdx < neighbors.size(); tetIdx++)
            {
                if (mTetrahedra[neighbors[tetIdx]].mFracturedThisFrame)
                    continue;

                fractured = true;

                FractureTetrahedra(neighbors[tetIdx]);

                // deletion in descending order so we don't invalidate the indices
                std::sort(mTetrahedraToDelete.begin(), mTetrahedraToDelete.end(), std::greater<size_t>());
                for (size_t idx = 0; idx < mTetrahedraToDelete.size(); idx++)
                    mTetrahedra.erase(mTetrahedra.begin() + mTetrahedraToDelete[idx]);

                for (auto& tet : mNewTetrahedra)
                {
                    tet.mFracturedThisFrame = true;
                    mTetrahedra.push_back(tet);
                }

                mTetrahedraToDelete.clear();
                mNewTetrahedra.clear();
            }
        }

        for (auto& tet : mTetrahedra)
            tet.mFracturedThisFrame = false;
    }

    ////////////////////////////////////////////////////////////////////////////////

    namespace
    {
        double DistPointPlane(
            const glm::vec3& point,
            const glm::vec3& planeNormal,
            const glm::vec3& planePosition
        )
        {
            gte::Vector3<float> normal
            {
                planeNormal.x,
                planeNormal.y,
                planeNormal.z
            };
            gte::Vector3<float> origin
            {
               planePosition.x,
               planePosition.y,
               planePosition.z
            };
            gte::Plane3<float> plane(normal, origin);

            gte::Vector3<float> position
            {
               point.x,
               point.y,
               point.z
            };

            // Get signed distance of all vertices to plane
            gte::DCPQuery<float, gte::Vector3<float>, gte::Plane3<float>> distanceQuery;
            auto results = distanceQuery(position, plane);
            return results.signedDistance;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////

    void FractureContext::DetermineSnapping()
    {
        mSnapToFace = false;
        mSnapToEdge = false;

        const double cDistanceTolerance = 0.05;
        const double cAngularSeparation = 0.1; // radians? double check what values paper uses
        const auto& tetrahedron = mTetrahedra[mTetrahedraIdx];

        size_t numDistanceTriggered = 0;
        std::vector<size_t> closeVertexIndices;

        // first check how close the other vertices are to the plane
        // and snap to the close ones

        for (size_t i = 0; i < 4; i++)
        {
            auto otherVertexId = tetrahedron.mIndices[i];
            if (otherVertexId == mFractureNodeIdx)
                continue;

            double dist = DistPointPlane(mVertices[otherVertexId].mPosition, mFracturePlaneNormal, mFractureNodePosition);
            bool triggeredTolerance = (std::abs(dist) <= cDistanceTolerance);

            if (triggeredTolerance)
            {
                numDistanceTriggered++;
                closeVertexIndices.push_back(otherVertexId);
            }
        }

        if (numDistanceTriggered != 0)
        {
			if (numDistanceTriggered == 1)
			{
				mSnapToEdge = true;
                mSnappingEdgeId = GetEdgeId(closeVertexIndices[0], mFractureNodeIdx);
				return;
			}
			else if (numDistanceTriggered == 2)
			{
				mSnapToFace = true;
                mSnappingFaceId = GetFaceId(closeVertexIndices[0], closeVertexIndices[1], mFractureNodeIdx);
				return;
			}
			else
			{
				// warn small tet
				return;
			}
        }

        // no nodes were close so we go to angular tolerance between the fracture plane and 
        // edges of the tetrahedron

        // angle between line and plane
        // project line onto plane
        // get angle between projection and line

        size_t numAngularTriggered = 0;

        // need to get indices of the other three vertices then get the edges
        // check angular distance from edge to plane

        const auto& otherVertices = tetrahedron.GetOtherVertices(mFractureNodeIdx);
        std::vector<size_t> snappedVertexIds;

        for (int i = 0; i < 3; i++)
        {
            auto otherVertexId = otherVertices[i];
            auto edgeDir = glm::normalize(mVertices[mFractureNodeIdx].mPosition - mVertices[otherVertexId].mPosition);
            auto d = glm::dot(mFracturePlaneNormal, edgeDir);
            // if we dot with the normal again
            // dot(normal, edgeDir) - dot(normal, edgeDir)*dot(normal, normal) = 0
            auto projectedEdgeDir = glm::normalize(edgeDir - d * mFracturePlaneNormal);
            // angle is between 0 and pi
            auto angle = glm::acos(glm::dot(projectedEdgeDir, edgeDir));

            if (angle < cAngularSeparation || angle >(glm::pi<float>() - cAngularSeparation))
            {
                numAngularTriggered++;
                snappedVertexIds.push_back(otherVertexId);
            }
        }

        if (numAngularTriggered == 0)
        {
            // this is the no snap, regular fracture case
            return;
        }
        else if (numAngularTriggered == 1)
        {
            mSnapToEdge = true;
            mSnappingEdgeId = GetEdgeId(mFractureNodeIdx, snappedVertexIds[0]);
        }
        
        // warn poorly conditioned tet, snap to one of the closest edges
        mSnapToEdge = true;
        mSnappingEdgeId = GetEdgeId(mFractureNodeIdx, snappedVertexIds[0]);

    }

    ////////////////////////////////////////////////////////////////////////////////

    void FractureContext::FaceSnappedFracture()
    {
        AssignTetToSide(mSnappedFaceNormal);
    }

    ////////////////////////////////////////////////////////////////////////////////

    void FractureContext::EdgeSnappedFracture()
    {
        // get single intersection between not the snapped node and not the fracture not

        auto& tetrahedron = mTetrahedra[mTetrahedraIdx];
        const auto& otherNodes = tetrahedron.GetOtherVertices(mSnappingEdgeId[0], mSnappingEdgeId[1]);
        auto splitEdgeId = GetEdgeId(otherNodes[0], otherNodes[1]);
        const auto& edgePos0 = mVertices[otherNodes[0]].mPosition;
        const auto& edgePos1 = mVertices[otherNodes[1]].mPosition;

        double d;
        if (!PlaneIntersectEdge(mFractureNodePosition, mSnappedEdgeNormal, edgePos0, edgePos1, d, nullptr))
            throw std::exception("Failure during edge snapped fracture, plane did not intersect remaining edge.");

        mTetrahedraToDelete.push_back(mTetrahedraIdx);

        // need to create a new node on the edge
        size_t newEdgeNode = CreateEdgeVertex(mSnappingEdgeId, d);
        size_t snappedNodeId = (mSnappingEdgeId[0] == mFractureNodeIdx ? mSnappingEdgeId[1] : mSnappingEdgeId[0]);

        const auto& otherVertices = tetrahedron.GetOtherVertices(mFractureNodeIdx, snappedNodeId);
        size_t positiveOtherNode;
        size_t negativeOtherNode;
        SeparateNodes(positiveOtherNode, negativeOtherNode, otherVertices, mFractureNodePosition, mSnappedEdgeNormal);

        if (IsIsolatedEdge(splitEdgeId))
        {
            size_t positiveEdgeNode = newEdgeNode;
            size_t negativeEdgeNode = CloneVertex(positiveEdgeNode);

            mNewTetrahedra.push_back(Tetrahedra(mFractureNodeIdx, snappedNodeId, positiveEdgeNode, positiveOtherNode));
            mNewTetrahedra.push_back(Tetrahedra(mFractureNodeIdx, snappedNodeId, negativeEdgeNode, negativeOtherNode));

            // no neighbor fracturing necessary because this edge was isolated, which also means there's no face neighbor
            return;
        }

        mNewTetrahedra.push_back(Tetrahedra(mFractureNodeIdx, snappedNodeId, newEdgeNode, positiveOtherNode));
        mNewTetrahedra.push_back(Tetrahedra(mFractureNodeIdx, snappedNodeId, newEdgeNode, negativeOtherNode));

        const auto& fracturedFaceIds = tetrahedron.GetOtherVertices(mFractureNodeIdx);
        NeighborFaceFracture(GetFaceId(fracturedFaceIds[0], fracturedFaceIds[1], fracturedFaceIds[2]), { newEdgeNode }, { splitEdgeId }, mSnappedEdgeNormal);
        NeighborEdgeFracture(splitEdgeId, newEdgeNode, mSnappedEdgeNormal);
    }

    ////////////////////////////////////////////////////////////////////////////////

    void FractureContext::FractureTetrahedra(size_t tetIdx)
    {
        if (tetIdx >= mTetrahedra.size())
            throw std::exception("Fracture tet index is invalid.");

        mTetrahedraIdx = tetIdx;

        std::cout << "Fracturing Node " << mFractureNodeIdx << " and Tetrahedra " << mTetrahedraIdx << std::endl;

        // done, but would be good to test
        DetermineSnapping();

        if (mSnapToFace)
        {
            CalculateSnappedFaceNormal();
            return FaceSnappedFracture();
        }

        if (mSnapToEdge)
        { 
            CalculateSnappedEdgeNormal();
            return EdgeSnappedFracture();
        }

        // todo NeighborFaceFracture still needs work

        std::vector<glm::ivec2> splitEdges;
        std::vector<double> parametricDistances;
        PlaneIntersectTetrahedraEdges(splitEdges, parametricDistances);
        
        // done
        if (splitEdges.size() == 0)
            return AssignTetToSide(mFracturePlaneNormal);

        if (splitEdges.size() == 1)
            throw std::exception("Edge tolerance failed to detect single split edge case.");

        // still need impl of edge case 2
        // tetrahedra creation code from intersection context
        if (splitEdges.size() == 2)
            return RegularFracture({ splitEdges[0], splitEdges[1] }, { parametricDistances[0], parametricDistances[1] });

        throw std::exception("More than 2 split edges, serious edge case detected.");
    }

    ////////////////////////////////////////////////////////////////////////////////

    void FractureContext::AssignTetToSide(const glm::dvec3& normal)
    {
        size_t numPosition = 0;
        size_t numNegative = 0;

        auto& tetrahedra = mTetrahedra[mTetrahedraIdx];

        for (size_t i = 0; i < 4; i++)
        {
            double dist = DistPointPlane(mVertices[tetrahedra.mIndices[i]].mPosition, normal, mFractureNodePosition);

            if (dist >= -1e-7)
                numPosition++;
            if (dist <= 1e-7)
                numNegative++;
        }

        tetrahedra.mFracturedThisFrame = true;

        if (numPosition == 4)
        {
            // assign tet to + side
            // nothing to do here (it already has the + one)
            return;
        }

        if (numNegative == 4)
        {
            // assign tet to - side
            tetrahedra.ReplaceVertex(mFractureNodeIdx, mNegativeFractureNodeIdx);
            return;
        }
        
        throw std::exception("Attempted to assign tetrahedra to plane side but vertices were on both sides.");
    }

    ////////////////////////////////////////////////////////////////////////////////

    void FractureContext::RegularFracture(
        const std::array<glm::ivec2, 2>& edgeIds,
        const std::array<double, 2>& parametricDistance)
    {
        const auto& tet = mTetrahedra[mTetrahedraIdx];
        const auto& otherVertices = tet.GetOtherVertices(mFractureNodeIdx);

        // Create two new nodes
        auto edgeVertex0 = CreateEdgeVertex(edgeIds[0], parametricDistance[0]);
        auto edgeVertex1 = CreateEdgeVertex(edgeIds[1], parametricDistance[1]);
        
        std::vector<size_t> positiveVertices;
        std::vector<size_t> negativeVertices;
        SeparateNodes(positiveVertices, negativeVertices, { otherVertices[0], otherVertices[1], otherVertices[2] }, mFractureNodePosition, mFracturePlaneNormal);

        if (!
            ((positiveVertices.size() == 1 && negativeVertices.size() == 2)
            || 
            (positiveVertices.size() == 2 && negativeVertices.size() == 1))
            )
        {
            throw std::exception("Expected 2 positive, 1 negative or 1 positive, 2 negative vertices during regular fracture.");
        }

        auto dp = mFractureNodeIdx;
        auto dm = mNegativeFractureNodeIdx;
        auto ep = edgeVertex0;
        auto em = edgeVertex0;
        auto fp = edgeVertex1;
        auto fm = edgeVertex1;

        if (IsIsolatedEdge(edgeIds[0]))
            em = CloneVertex(ep);
        
        if (IsIsolatedEdge(edgeIds[1]))
            fm = CloneVertex(fp);

        if (positiveVertices.size() == 1 && negativeVertices.size() == 2)
        {
            auto a = positiveVertices[0];
            auto b = negativeVertices[0];
            auto c = negativeVertices[1];
            
            mNewTetrahedra.push_back(Tetrahedra(a, dp, ep, fp));
            mNewTetrahedra.push_back(Tetrahedra(b, fm, em, dm));
            mNewTetrahedra.push_back(Tetrahedra(b, c, fm, dm));
        }
        else if (positiveVertices.size() == 2 && negativeVertices.size() == 1)
        {
            auto a = negativeVertices[0];
            auto b = positiveVertices[0];
            auto c = positiveVertices[1];
            
            mNewTetrahedra.push_back(Tetrahedra(a, dm, em, fm));
            mNewTetrahedra.push_back(Tetrahedra(b, fp, ep, dp));
            mNewTetrahedra.push_back(Tetrahedra(b, c, fp, dp));
        }

        // Mark this one for deletion
        mTetrahedraToDelete.push_back(mTetrahedraIdx);

        NeighborFaceFracture(GetFaceId(otherVertices[0], otherVertices[1], otherVertices[2]), { edgeVertex0, edgeVertex1 }, { edgeIds[0], edgeIds[1] }, mFracturePlaneNormal);
        NeighborEdgeFracture(edgeIds[0], edgeVertex0, mFracturePlaneNormal);
        NeighborEdgeFracture(edgeIds[1], edgeVertex1, mFracturePlaneNormal);
    }

    ////////////////////////////////////////////////////////////////////////////////

    void FractureContext::NeighborFaceFracture(const glm::ivec3& faceId, const std::vector<size_t>& newNodeIds, const std::vector<glm::ivec2>& splitEdges, const glm::dvec3& planeNormal)
    {
        // Find face neighbor if it exists
        size_t neighborTetId = -1;
        for (size_t tetId = 0; tetId < mTetrahedra.size(); tetId++)
        {
            if (tetId == mTetrahedraIdx)
                continue;

            if (mTetrahedra[tetId].ContainsFaceIndex(faceId))
            {
                neighborTetId = tetId;
                break;
            }
        }

        // no neighbor exists on the fractured face
        if (neighborTetId == -1)
            return;

        mFracturedNeighborFaceTetId = neighborTetId;
        mTetrahedraToDelete.push_back(mFracturedNeighborFaceTetId);
        const auto& neighborTet = mTetrahedra[mFracturedNeighborFaceTetId];

        if (newNodeIds.size() == 1)
        {
            if (!mSnapToEdge)
                throw std::exception("Neighbor face fracture expected snapped to edge but it was not.");

            auto snappedNodeId = GetNonFractureNode(mSnappingEdgeId);
            if (!neighborTet.ContainsVertexIndex(snappedNodeId))
                throw std::exception("Single new node neighbor face fracture, expected the neighbor tet to contain the snapping node but it did not.");

            // the snapped vertex should be in the face that is being fracture
            std::vector<size_t> vertices;
            if (faceId[0] != snappedNodeId)
                vertices.push_back(faceId[0]);
            if (faceId[1] != snappedNodeId)
                vertices.push_back(faceId[1]);
            if (faceId[2] != snappedNodeId)
                vertices.push_back(faceId[2]);

            size_t positiveVertex;
            size_t negativeVertex;
            SeparateNodes(positiveVertex, negativeVertex, { vertices[0], vertices[1] }, mFractureNodePosition, mSnappedEdgeNormal);

            auto otherNode = neighborTet.GetOtherVertex({vertices[0], vertices[1], snappedNodeId});


            mNewTetrahedra.push_back(Tetrahedra(positiveVertex, snappedNodeId, newNodeIds[0], otherNode));
            mNewTetrahedra.push_back(Tetrahedra(negativeVertex, snappedNodeId, newNodeIds[0], otherNode));

            return;
        }
        
        if (newNodeIds.size() == 2)
        {
            // these three vertices form a triangle and should form a tetrahedra with the vertex not present in the fractured face
            auto a = GetCommonVertexFromEdges(splitEdges[0], splitEdges[1]);
            auto e = newNodeIds[0];
            auto f = newNodeIds[1];

            auto o = neighborTet.GetOtherVertex({ (size_t)faceId[0], (size_t)faceId[1], (size_t)faceId[2] });

            // let c be on the edge that e split,
            // let d be on the edge that f split.
            // a is on both these edges
            // newNodesIds[0] is on the 0th split edge, so splitEdges[0] 
            auto c = (splitEdges[0][0] == a ? splitEdges[0][1] : splitEdges[0][0]);
            auto d = (splitEdges[1][0] == a ? splitEdges[1][1] : splitEdges[1][0]);

            mNewTetrahedra.push_back(Tetrahedra(a, e, f, o));
            mNewTetrahedra.push_back(Tetrahedra(e, f, d, o));
            mNewTetrahedra.push_back(Tetrahedra(e, d, c, o));

            return;
        }

        throw std::exception("Attempted neighbor face fracture with invalid number of new node ids.");
    }

    ////////////////////////////////////////////////////////////////////////////////

    void FractureContext::NeighborEdgeFracture(const glm::ivec2& edgeId, size_t newNodeId, const glm::dvec3& planeNormal)
    {
        // need to get all tetrahedra that contain this edge and are not the original tet or the face tet
        std::vector<size_t> neighborTetrahedraIds;
        for (size_t neighborTetId = 0; neighborTetId < mTetrahedra.size(); neighborTetId++)
        {
            if (neighborTetId == mTetrahedraIdx || (mFracturedNeighborFaceTet && neighborTetId == mFracturedNeighborFaceTetId))
                continue;
            
            const auto& neighborTet = mTetrahedra[neighborTetId];
            if (!neighborTet.ContainsEdgeIndex(edgeId))
                continue;

            mTetrahedraToDelete.push_back(neighborTetId);
            
            auto otherNodes = neighborTet.GetOtherVertices(edgeId[0], edgeId[1]);
            size_t posEdgeNode;
            size_t negEdgeNode;
            SeparateNodes(posEdgeNode, negEdgeNode, { (size_t)edgeId[0], (size_t)edgeId[1] }, mFractureNodePosition, planeNormal);

            mNewTetrahedra.push_back(Tetrahedra(newNodeId, otherNodes[0], otherNodes[1], posEdgeNode));
            mNewTetrahedra.push_back(Tetrahedra(newNodeId, otherNodes[0], otherNodes[1], negEdgeNode));
        }
    }

    ////////////////////////////////////////////////////////////////////////////////

    void FractureContext::CalculateSnappedFaceNormal()
    {
        const auto& p0 = mVertices[mSnappingFaceId.x].mPosition;
        const auto& p1 = mVertices[mSnappingFaceId.y].mPosition;
        const auto& p2 = mVertices[mSnappingFaceId.z].mPosition;

        mSnappedFaceNormal = glm::normalize(glm::cross(p2 - p0, p1 - p0));

        if (glm::dot(mSnappedFaceNormal, mFracturePlaneNormal) < 0)
            mSnappedFaceNormal *= -1.0;
    }

    ////////////////////////////////////////////////////////////////////////////////

    void FractureContext::CalculateSnappedEdgeNormal()
    {
        // we want the minimal change that makes the edge lay in the plane
        // subtract from fracture plane normal the portion that dots with edge
        // so snappedFracturePlane = fracturePlane - dot(fracturePlane, edge) * edge
        // so that dot(snappedFracturePlane, edge) = dot(fracturePlane, edge) - dot(fracturePlane, edge) * dot(edge, edge) = 0
        auto edge = glm::normalize(mVertices[mSnappingEdgeId.x].mPosition - mVertices[mSnappingEdgeId.y].mPosition);
        mSnappedEdgeNormal = glm::normalize(mFracturePlaneNormal - glm::dot(mFracturePlaneNormal, edge) * edge);
    }

    ////////////////////////////////////////////////////////////////////////////////

    size_t FractureContext::CloneVertex(size_t vertexId)
    {
        Vertex vertex;
        vertex.mPosition = mVertices[vertexId].mPosition;
        vertex.mMaterialCoordinates = mVertices[vertexId].mMaterialCoordinates;
        vertex.mVelocity = mVertices[vertexId].mVelocity;

        mVertices.push_back(vertex);

        return (mVertices.size() - 1);
    }

    ////////////////////////////////////////////////////////////////////////////////

    size_t FractureContext::CreateEdgeVertex(const glm::ivec2& edgeId, double parametricDistance)
    {
        Vertex vertex;
        const auto& p0 = mVertices[edgeId.x].mPosition;
        const auto& p1 = mVertices[edgeId.y].mPosition;

        vertex.mPosition = p0 + (p1 - p0) * parametricDistance;

        const auto& m0 = mVertices[edgeId.x].mMaterialCoordinates;
        const auto& m1 = mVertices[edgeId.y].mMaterialCoordinates;
        vertex.mMaterialCoordinates = m0 + (m1 - m0) * parametricDistance;

        const auto& v0 = mVertices[edgeId.x].mVelocity;
        const auto& v1 = mVertices[edgeId.y].mVelocity;
        // Not sure this is the right way to do this but it seems rather reasonable
        vertex.mVelocity = v0 + (v1 - v0) * parametricDistance;

        mVertices.push_back(vertex);

        return (mVertices.size() - 1);
    }

    ////////////////////////////////////////////////////////////////////////////////

    void FractureContext::SeparateNodes(size_t& outPositiveNode, size_t& outNegativeNode, const std::array<size_t, 2>& nodes, const glm::dvec3& planePos, const glm::dvec3& planeNormal) const
    {
        auto d0 = DistPointPlane(mVertices[nodes[0]].mPosition, planeNormal, planePos);
        auto d1 = DistPointPlane(mVertices[nodes[1]].mPosition, planeNormal, planePos);

        size_t numPos = 0;
        size_t numNeg = 0;

        if (d0 < 0)
            numNeg++;
        if (d1 < 0)
            numNeg++;
        if (d0 > 0)
            numPos++;
        if (d1 > 0)
            numPos++;

        if (!(numPos == 1 && numNeg == 1))
            throw std::exception("Expected one positive and one negative vertex in isolated edge fracture.");

        outPositiveNode = (d0 > 0 ? nodes[0] : nodes[1]);
        outNegativeNode = (d0 < 0 ? nodes[0] : nodes[1]);
    }

    ////////////////////////////////////////////////////////////////////////////////

    void FractureContext::SeparateNodes(std::vector<size_t>& outPositiveNodes, std::vector<size_t>& outNegativeNodes, const std::vector<size_t>& nodes, const glm::dvec3& planePos, const glm::dvec3& planeNormal) const
    {
        for (size_t idx : nodes)
        {
            auto d = DistPointPlane(mVertices[idx].mPosition, planeNormal, planePos);
            
            if (d > 0)
                outPositiveNodes.push_back(idx);
            else if (d < 0)
                outNegativeNodes.push_back(idx);
        }
    }

    ////////////////////////////////////////////////////////////////////////////////

    bool FractureContext::PlaneIntersectEdge(const glm::dvec3& planePos, const glm::dvec3& planeNormal, const glm::dvec3& edgePos0, const glm::dvec3& edgePos1, double& d, glm::dvec3* intersectionPos) const
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

    bool FractureContext::IsIsolatedEdge(const glm::ivec2& edgeId) const
    {
        size_t counter = 0;

        for (const auto& tet : mTetrahedra)
        {
            if (tet.ContainsEdgeIndex(edgeId))
                counter++;
        }

        if (counter == 0)
            throw std::exception("Failed to find isolated edge id in any tetrahedra.");

        return (counter == 1);
    }

    ////////////////////////////////////////////////////////////////////////////////

    std::vector<size_t> FractureContext::GetTetrahedraNeighbors(size_t nodeIdx) const
    {
        std::vector<size_t> neighborIds;

        for (size_t tetId = 0; tetId < mTetrahedra.size(); tetId++)
        {
            if (mTetrahedra[tetId].ContainsVertexIndex(nodeIdx))
                neighborIds.push_back(tetId);
        }

        return neighborIds;
    }

    ////////////////////////////////////////////////////////////////////////////////

    size_t FractureContext::GetNonFractureNode(const glm::ivec2& edge) const
    {
        if (!(edge[0] == mFractureNodeIdx || edge[1] == mFractureNodeIdx))
            throw std::exception("Expected fracture node to be present in GetNonFractureNode.");

        return (edge[0] == mFractureNodeIdx ? edge[1] : edge[0]);
    }

    ////////////////////////////////////////////////////////////////////////////////

    void FractureContext::PlaneIntersectTetrahedraEdges(std::vector<glm::ivec2>& outVertexIds, std::vector<double>& parametricDistance) const
    {
        const auto& tet = mTetrahedra[mTetrahedraIdx];
        const auto& otherVertices = tet.GetOtherVertices(mFractureNodeIdx);

        // the three possible edges are those that do not include the fracture node
        auto edge0 = GetEdgeId(otherVertices[0], otherVertices[1]);
        auto edge1 = GetEdgeId(otherVertices[0], otherVertices[2]);
        auto edge2 = GetEdgeId(otherVertices[1], otherVertices[2]);

        double d;

        if (PlaneIntersectEdge(mFractureNodePosition, mFracturePlaneNormal, mVertices[edge0[0]].mPosition, mVertices[edge0[1]].mPosition, d, nullptr))
        {
            outVertexIds.push_back(edge0);
            parametricDistance.push_back(d);
        }

        if (PlaneIntersectEdge(mFractureNodePosition, mFracturePlaneNormal, mVertices[edge1[0]].mPosition, mVertices[edge1[1]].mPosition, d, nullptr))
        {
            outVertexIds.push_back(edge1);
            parametricDistance.push_back(d);
        }

        if (PlaneIntersectEdge(mFractureNodePosition, mFracturePlaneNormal, mVertices[edge2[0]].mPosition, mVertices[edge2[1]].mPosition, d, nullptr))
        {
            outVertexIds.push_back(edge2);
            parametricDistance.push_back(d);
        }
    }

    ////////////////////////////////////////////////////////////////////////////////
}