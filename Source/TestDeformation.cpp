#include "TestDeformation.hpp"
#include "ProtoConverter.hpp"

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

		for (auto & pair : mIdToTetrahedra)
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
            //for (size_t idx = 0; idx < numVertices; idx++)
            //{
            //    if (mVertices[idx].mLargestEigenvalue > mToughness)
            //        FractureNode(idx, mVertices[idx].mPrincipalEigenVector);

            //    if (mVertices.size() >= mMaxNumVertices)
            //        break;
            //}
            
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
			}

            if (maxEigenValueNodeId != -1)
            {
                std::cout << "Max eigen value : " << maxEigenValue << std::endl;
                FractureNode(maxEigenValueNodeId, mVertices[maxEigenValueNodeId].mPrincipalEigenVector);
            }

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

    namespace
    {
        struct EdgeIntersection
        {
            // Interpreted as from mEdgeId[0] to mEdgeId[1]
            // so that a value of 0 corresponds to mEdgeId[0]
            double mParam = 0;
            glm::ivec2 mEdgeId;
        };

        // unsigned angle, our edges don't have a direction because we sort based on vertex index
        double AngleEdgePlane(const glm::dvec3& edgeDir, const glm::dvec3& planeNormal)
        {
            auto d = glm::dot(planeNormal, edgeDir);
            auto projectedEdgeDir = glm::normalize(edgeDir - d * planeNormal);
            return std::abs(glm::acos(glm::dot(projectedEdgeDir, edgeDir)));
        }

        size_t FindFaceNeighbor(
            const std::unordered_map<size_t, Tetrahedra>& tetrahedra,
            size_t mainTetId,
            const glm::ivec3 & faceId
        )
        {
            for (const auto& pair : tetrahedra)
            {
                auto tetId = pair.first;
                if (tetId == mainTetId)
                    continue;

                const auto& tet = pair.second;
                if (tet.ContainsFaceIndex(faceId))
                    return tetId;
            }

            return -1;
        }

        std::vector<size_t> FindEdgeNeighbors(
            const std::unordered_map<size_t, Tetrahedra>& tetrahedra,
            const glm::ivec2& edgeId
        )
        {
            std::vector<size_t> neighborIds;

            for (const auto& pair : tetrahedra)
            {
                auto tetId = pair.first;
                const auto& tet = pair.second;
                if (tet.ContainsEdgeIndex(edgeId))
                    neighborIds.push_back(tetId);
            }

            return neighborIds;
        }

        enum class FractureOptions
        {
            ePositiveSide,
            eNegativeSide,
            eFracture
        };

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

        // If this returns false we will create a new node
        // instead of snapping to the existing one
        // This should only return false if there's another tet that will be assigned ot the + side or
        // that the fracture plane will split the tet.
        bool NoPositiveNodeEdgeCase(
            const std::unordered_set<size_t>& neighborTets,
            const std::unordered_map<size_t, Tetrahedra>& idToTetrahedra,
            const std::vector<Vertex> & vertices,
            const glm::vec3& fracturePlanePosition,
            const glm::vec3& fracturePlaneNormal)
        {
            for (const auto& neighborId : neighborTets)
            {
                const auto& tet = idToTetrahedra.at(neighborId);

                for (const auto& nodeId : tet.mIndices)
                {
                    auto d = DistPointPlane(vertices[nodeId].mPosition, fracturePlaneNormal, fracturePlanePosition);
                    if (d > 0.2)
                    {
                        std::cout << "Returning false for the no positive node edge case because d = " << d << "." << std::endl;
                        return false;
                    }
                }
            }

            std::cout << "Returning true for the no positive node edge case." << std::endl;
            return true;
        }
    }

    void FractureContext::Fracture()
    {
		const auto& neighbors = GetTetrahedraNeighbors(mFractureNodeIdx);
        std::unordered_set<size_t> neighborSet(neighbors.begin(), neighbors.end());

		const double cDistanceTolerance = 0.2;
		const double cAngularSeparation = 0.34; // radians? double check what values paper uses

		// it would be nice to be able to set these tolerances to large numbers to force always snapping
		// i don't think we have this capabilitiy atm

		std::array<double, 4> nodeDistToPlane;
		std::array<double, 3> edgeAngleToPlane; // really only care about the 3 edges connected to the fracture node
		std::vector<EdgeIntersection> edgeIntersections;

		bool snapToEdge = false;
		glm::ivec2 snapEdgeId;
		bool snapToFace = false;
		glm::ivec3 snapFaceId;

		size_t negativeFractureNodeId = -1;

        bool noPositiveEdgeCase = NoPositiveNodeEdgeCase(neighborSet, mIdToTetrahedra, mVertices, mFractureNodePosition, mFracturePlaneNormal);

		auto inputStateFunctor = [&](const Tetrahedra& tet, size_t fractureNodeId)
		{
			const auto& edges = tet.GetEdges();

			for (size_t i = 0; i < 4; i++)
				nodeDistToPlane[i] = DistPointPlane(mVertices[tet.mIndices[i]].mPosition, mFracturePlaneNormal, mFractureNodePosition);

			const auto& otherNodes = tet.GetOtherVertices(fractureNodeId);

			for (size_t i = 0; i < 3; i++)
			{
				auto edgeId = GetEdgeId(fractureNodeId, otherNodes[i]);
				edgeAngleToPlane[i] = AngleEdgePlane(glm::normalize(mVertices[edgeId[1]].mPosition - mVertices[edgeId[0]].mPosition), mFracturePlaneNormal);
			}
		};

		// we won't always fracture, so the negative node is computed lazily
		auto getNegativeFractureNodeId = [&]()
		{
            if (negativeFractureNodeId == -1)
            {
                if (noPositiveEdgeCase)
                    negativeFractureNodeId = mFractureNodeIdx;
                else
					negativeFractureNodeId = CloneVertex(mFractureNodeIdx);
            }

			return negativeFractureNodeId;
		};

		auto checkTolerances = [&](const auto& tet, size_t fractureNodeId)
		{
			size_t numDistanceTriggered = 0;
			std::vector<size_t> closeVertexIndices;

			for (size_t i = 0; i < nodeDistToPlane.size(); i++)
			{
				if (tet.mIndices[i] == fractureNodeId)
					continue;

				if (std::abs(nodeDistToPlane[i]) <= cDistanceTolerance)
				{
					numDistanceTriggered++;
					closeVertexIndices.push_back(tet.mIndices[i]);
				}
			}

			if (numDistanceTriggered != 0)
			{
				if (numDistanceTriggered == 1)
				{
					snapToEdge = true;
					snapEdgeId = GetEdgeId(fractureNodeId, closeVertexIndices[0]);
					return;
				}
				else if (numDistanceTriggered == 2)
				{
					snapToFace = true;
					snapFaceId = GetFaceId(fractureNodeId, closeVertexIndices[0], closeVertexIndices[1]);
					return;
				}
				else
				{
					std::cout << "Small tet!" << std::endl;
					return;
				}
			}

			const auto& otherNodes = tet.GetOtherVertices(fractureNodeId);
			std::vector<size_t> snappedVertexIds;

			for (size_t i = 0; i < 3; i++)
			{
				auto angle = edgeAngleToPlane[i];
				if (angle < cAngularSeparation || angle >(glm::pi<float>() - cAngularSeparation))
					snappedVertexIds.push_back(otherNodes[i]);
			}

			if (snappedVertexIds.empty())
			{
				// this is the no snap, regular fracture case
				return;
			}
			else if (snappedVertexIds.size() == 1)
			{
				snapToEdge = true;
				snapEdgeId = GetEdgeId(mFractureNodeIdx, snappedVertexIds[0]);
				return;
			}
			else if (snappedVertexIds.size() == 2)
			{
				snapToFace = true;
				snapFaceId = GetFaceId(fractureNodeId, snappedVertexIds[0], snappedVertexIds[1]);
				return;
			}
			else
			{
				// poorly conditioned tet
				std::cout << "Very poorly conditioned tet as determined by angular threshold." << std::endl;
			}
		};

		auto determineFractureOption = [&]()
		{
			size_t numPos = 0;
			size_t numNeg = 0;

			for (const auto& dist : nodeDistToPlane)
			{
				if (dist <= 1e-5)
					numNeg++;
				if (dist >= -1e-5)
					numPos++;
			}

			if (numPos == 4)
				return FractureOptions::ePositiveSide;
			if (numNeg == 4)
				return FractureOptions::eNegativeSide;

			return FractureOptions::eFracture;
		};

        while (!neighborSet.empty())
		{
            std::cout << "Num neighbors : " << neighborSet.size() << std::endl;

            auto fracturingTetId = *neighborSet.begin();

			// Clear shared state
			snapToEdge = false;
			snapToFace = false;
			edgeIntersections.clear();

			// Start fracturing
			const auto& tet = mIdToTetrahedra[fracturingTetId];
			inputStateFunctor(tet, mFractureNodeIdx);

			checkTolerances(tet, mFractureNodeIdx);
			if (snapToEdge)
			{
				auto edge = glm::normalize(mVertices[snapEdgeId[0]].mPosition - mVertices[snapEdgeId[1]].mPosition);
				mFracturePlaneNormal = glm::normalize(mFracturePlaneNormal - glm::dot(mFracturePlaneNormal, edge) * edge);

				inputStateFunctor(tet, mFractureNodeIdx);
			}
			else if (snapToFace)
			{
				const auto& p0 = mVertices[snapFaceId.x].mPosition;
				const auto& p1 = mVertices[snapFaceId.y].mPosition;
				const auto& p2 = mVertices[snapFaceId.z].mPosition;

				auto faceNormal = glm::normalize(glm::cross(p2 - p0, p1 - p0));

				if (glm::dot(faceNormal, mFracturePlaneNormal) < 0)
					faceNormal *= -1.0;

				mFracturePlaneNormal = faceNormal;
				inputStateFunctor(tet, mFractureNodeIdx);
			}

            if (snapToFace)
                std::cout << "Snapping to Face" << std::endl;
            else if (snapToEdge)
                std::cout << "Snapping to Edge" << std::endl;
            else
                std::cout << "No snapping" << std::endl;

			// Check whether we want to fracture this tet
			auto option = determineFractureOption();
			if (option == FractureOptions::ePositiveSide)
			{
                std::cout << "Assigning to positive side" << std::endl;
				// Assign Tet to Side
				mNewTetrahedra.push_back(tet);
				mTetrahedraIdsToDelete.insert(fracturingTetId);
			}
			else if (option == FractureOptions::eNegativeSide)
			{
                std::cout << "Assigning to negative side" << std::endl;

				auto negativeFractureNodeIdx = getNegativeFractureNodeId();

				// Assign Tet to Side
				auto copy = tet;
				copy.ReplaceVertex(mFractureNodeIdx, negativeFractureNodeIdx);
				mNewTetrahedra.push_back(copy);
				mTetrahedraIdsToDelete.insert(fracturingTetId);
			}
            else
            {
                // if snap to edge, only intersect the one other edge

                if (snapToEdge)
                {
                    const auto& otherNodes = tet.GetOtherVertices(mFractureNodeIdx, snapEdgeId[0] == mFractureNodeIdx ? snapEdgeId[1] : snapEdgeId[0]);
                    auto edgeId = GetEdgeId(otherNodes[0], otherNodes[1]);
                    double d;
                    if (PlaneIntersectEdge(mFractureNodePosition, mFracturePlaneNormal, mVertices[edgeId[0]].mPosition, mVertices[edgeId[1]].mPosition, d, nullptr))
                        edgeIntersections.push_back({ d, edgeId });
                }
                else
                {
                    // Intersect edges opposite to fracture node
                    const auto& otherNodes = tet.GetOtherVertices(mFractureNodeIdx);

                    {
                        auto edgeId = GetEdgeId(otherNodes[0], otherNodes[1]);
                        double d;
                        if (PlaneIntersectEdge(mFractureNodePosition, mFracturePlaneNormal, mVertices[edgeId[0]].mPosition, mVertices[edgeId[1]].mPosition, d, nullptr))
                            edgeIntersections.push_back({ d, edgeId });
                    }
                    {
                        auto edgeId = GetEdgeId(otherNodes[0], otherNodes[2]);
                        double d;
                        if (PlaneIntersectEdge(mFractureNodePosition, mFracturePlaneNormal, mVertices[edgeId[0]].mPosition, mVertices[edgeId[1]].mPosition, d, nullptr))
                            edgeIntersections.push_back({ d, edgeId });
                    }
                    {
                        auto edgeId = GetEdgeId(otherNodes[1], otherNodes[2]);
                        double d;
                        if (PlaneIntersectEdge(mFractureNodePosition, mFracturePlaneNormal, mVertices[edgeId[0]].mPosition, mVertices[edgeId[1]].mPosition, d, nullptr))
                            edgeIntersections.push_back({ d, edgeId });
                    }
                }

                std::cout << "Num edge intersections : " << edgeIntersections.size() << std::endl;

                if (edgeIntersections.size() == 0)
                {
                    // failed to assign tet to side
                    throw std::exception("Failed to assign tet to side despite plane not intersecting the tet.");
                }
                else if (edgeIntersections.size() == 1)
                {
                    if (!snapToEdge)
                        throw std::exception("One split edge, expected edge snap but it was false.");

                    auto negativeFractureNodeIdx = getNegativeFractureNodeId();

                    // edge snapped fracture
                    const auto& intersection = edgeIntersections.front();
                    const auto& splitEdgeId = intersection.mEdgeId;
                    size_t positiveSplitEdgeNode, negativeSplitEdgeNode;
                    SeparateNodes(positiveSplitEdgeNode, negativeSplitEdgeNode, { (size_t)splitEdgeId[0], (size_t)splitEdgeId[1] }, mFractureNodePosition, mFracturePlaneNormal);

                    size_t newEdgeNode = CreateEdgeVertex(splitEdgeId, intersection.mParam);
                    size_t snappedNodeId = (snapEdgeId[0] == mFractureNodeIdx ? snapEdgeId[1] : snapEdgeId[0]);

                    mTetrahedraIdsToDelete.insert(fracturingTetId);

                    if (IsIsolatedEdge(splitEdgeId))
                    {
                        size_t positiveNewEdgeNode = newEdgeNode;
                        size_t negativeNewEdgeNode = CloneVertex(positiveNewEdgeNode);

                        mNewTetrahedra.push_back(Tetrahedra(mFractureNodeIdx, snappedNodeId, positiveNewEdgeNode, positiveSplitEdgeNode));
                        mNewTetrahedra.push_back(Tetrahedra(negativeFractureNodeIdx, snappedNodeId, negativeNewEdgeNode, negativeSplitEdgeNode));
                    }
                    else
                    {
                        mNewTetrahedra.push_back(Tetrahedra(mFractureNodeIdx, snappedNodeId, newEdgeNode, positiveSplitEdgeNode));
                        mNewTetrahedra.push_back(Tetrahedra(negativeFractureNodeIdx, snappedNodeId, newEdgeNode, negativeSplitEdgeNode));

                        // case 1 - the face neighbor does not contain the snapped node, but does contain the fracture node
                        // this face contains the fracture node, and both nodes from the split edge
                        // for this case, tet0 contains (mFractureNodeIdx, the new edge node, the positive edge node, and the other node in the neighbor tet)
                        //              , tet1 contains (negative fracture node, the new edge node, the negative edge node, and the other node in the neighbor tet)
                        auto faceId0 = GetFaceId(mFractureNodeIdx, splitEdgeId[0], splitEdgeId[1]);
                        auto faceNeighborId0 = FindFaceNeighbor(mIdToTetrahedra, fracturingTetId, faceId0);
                        if (faceNeighborId0 != -1)
                        {
                            const auto& neighborTet = mIdToTetrahedra.at(faceNeighborId0);
                            mTetrahedraIdsToDelete.insert(faceNeighborId0);
                            mNewTetrahedra.push_back(Tetrahedra(mFractureNodeIdx, newEdgeNode, positiveSplitEdgeNode, neighborTet.GetOtherVertex({ (size_t)faceId0[0], (size_t)faceId0[1], (size_t)faceId0[2] })));
                            mNewTetrahedra.push_back(Tetrahedra(negativeFractureNodeIdx, newEdgeNode, negativeSplitEdgeNode, neighborTet.GetOtherVertex({ (size_t)faceId0[0], (size_t)faceId0[1], (size_t)faceId0[2] })));
                        }

                        // case 2 - the face neighbor does not contain the fracture node, but does contain the snapped node
                        // this face contains the snapped node, and both nodes from the split edge
                        // for this case, tet0 contains (snapped edge node, the new edge node, the positive edge node, and the other node in the neighbor tet)
                        //              , tet1 contains (snapped edge node, the new edge node, the negative edge node, and the other node in the neighbor tet)
                        auto faceId1 = GetFaceId(snappedNodeId, splitEdgeId[0], splitEdgeId[1]);
                        auto faceNeighborId1 = FindFaceNeighbor(mIdToTetrahedra, fracturingTetId, faceId1);
                        if (faceNeighborId1 != -1)
                        {
                            const auto& neighborTet = mIdToTetrahedra.at(faceNeighborId1);
                            mTetrahedraIdsToDelete.insert(faceNeighborId1);
                            mNewTetrahedra.push_back(Tetrahedra(snappedNodeId, newEdgeNode, positiveSplitEdgeNode, neighborTet.GetOtherVertex({ (size_t)faceId1[0], (size_t)faceId1[1], (size_t)faceId1[2] })));
                            mNewTetrahedra.push_back(Tetrahedra(snappedNodeId, newEdgeNode, negativeSplitEdgeNode, neighborTet.GetOtherVertex({ (size_t)faceId1[0], (size_t)faceId1[1], (size_t)faceId1[2] })));
                        }

                        const auto& edgeNeighbors = FindEdgeNeighbors(mIdToTetrahedra, splitEdgeId);

                        for (size_t neighborId : edgeNeighbors)
                        {
                            if (neighborId == fracturingTetId || neighborId == faceNeighborId0 || neighborId == faceNeighborId1)
                                continue;

                            const auto& edgeNeighorTet = mIdToTetrahedra[neighborId];
                            const auto& remainingNodes = edgeNeighorTet.GetOtherVertices(splitEdgeId[0], splitEdgeId[1]);
                            size_t posEdgeNode;
                            size_t negEdgeNode;
                            SeparateNodes(posEdgeNode, negEdgeNode, { (size_t)splitEdgeId[0], (size_t)splitEdgeId[1] }, mFractureNodePosition, mFracturePlaneNormal);

                            mTetrahedraIdsToDelete.insert(neighborId);

                            mNewTetrahedra.push_back(Tetrahedra(newEdgeNode, remainingNodes[0], remainingNodes[1], posEdgeNode));
                            mNewTetrahedra.push_back(Tetrahedra(newEdgeNode, remainingNodes[0], remainingNodes[1], negEdgeNode));
                        }
                    }
                }
                else if (edgeIntersections.size() == 2)
                {
                    const auto& otherNodes = tet.GetOtherVertices(mFractureNodeIdx);

                    // regular fracture
                    auto negativeFractureNodeIdx = getNegativeFractureNodeId();

                    const auto& edgeId0 = edgeIntersections[0].mEdgeId;
                    const auto& edgeId1 = edgeIntersections[1].mEdgeId;

                    // Create two new nodes
                    auto edgeVertex0 = CreateEdgeVertex(edgeId0, edgeIntersections[0].mParam);
                    auto edgeVertex1 = CreateEdgeVertex(edgeId1, edgeIntersections[1].mParam);

                    std::unordered_set<size_t> nonEdgeNeighbors{ fracturingTetId };

                    // This remains just for validation, can remove later
                    std::vector<size_t> positiveVertices;
                    std::vector<size_t> negativeVertices;
                    SeparateNodes(positiveVertices, negativeVertices, { otherNodes[0], otherNodes[1], otherNodes[2] }, mFractureNodePosition, mFracturePlaneNormal);

                    if (!
                        ((positiveVertices.size() == 1 && negativeVertices.size() == 2)
                            ||
                            (positiveVertices.size() == 2 && negativeVertices.size() == 1))
                        )
                    {
                        throw std::exception("Expected 2 positive, 1 negative or 1 positive, 2 negative vertices during regular fracture.");
                    }

                    // 7 unique nodes for the regular fracture case
                    // + fracture node
                    // - fracture node
                    // the node common to the two split edges
                    // new edge node 0
                    // the tetrahedral node that is on the same edge as the new edge node 0 and is not the common node
                    // new edge node 1
                    // the tetrahedral node ... new edge node 1 ... not common node
                    //
                    // if the common node is on the + side of the fracture plane, then we're in 'config 1' if not, we're in 'config 2'
                    // config 1 = 1 positive node , 2 negative
                    // config 2 = 2 positive nodes, 1 negative

					auto posFractureNode = mFractureNodeIdx;
					auto negFractureNode = getNegativeFractureNodeId();
					auto commonNode = GetCommonVertexFromEdges(edgeId0, edgeId1);
					auto newEdgeNode0 = edgeVertex0;
					auto edge0Node = (commonNode == edgeId0[0] ? edgeId0[1] : edgeId0[0]);
					auto newEdgeNode1 = edgeVertex1;
					auto edge1Node = (commonNode == edgeId1[0] ? edgeId1[1] : edgeId1[0]);

					// Only for isolated edges, doesn't affect neighbors
					size_t negNewEdgeNode0 = newEdgeNode0;
					size_t negNewEdgeNode1 = newEdgeNode1;

					if (IsIsolatedEdge(edgeId0))
						negNewEdgeNode0 = CloneVertex(newEdgeNode0);

					if (IsIsolatedEdge(edgeId1))
						negNewEdgeNode1 = CloneVertex(newEdgeNode1);

					auto d = DistPointPlane(mVertices[commonNode].mPosition, mFracturePlaneNormal, mFractureNodePosition);
					if (d < 0)
					{
						if (!(negativeVertices.size() == 1 && positiveVertices.size() == 2))
							throw std::exception("Expected 2 positive, 1 negative vertices during regular fracture.");

						std::swap(posFractureNode, negFractureNode);
						std::swap(newEdgeNode0, negNewEdgeNode0);
						std::swap(newEdgeNode1, negNewEdgeNode1);
					}

					// We choose the (newEdgeNode0, edge1Node) to be the diagonal edge
                    mTetrahedraIdsToDelete.insert(fracturingTetId);
					mNewTetrahedra.push_back(Tetrahedra(posFractureNode, commonNode, newEdgeNode0, newEdgeNode1));
					mNewTetrahedra.push_back(Tetrahedra(negFractureNode, negNewEdgeNode0, edge0Node, edge1Node));
					mNewTetrahedra.push_back(Tetrahedra(negFractureNode, negNewEdgeNode0, negNewEdgeNode1, edge1Node));

					// the first face neighbor, containing the fracture node, and both nodes from a single split edge
					{
						auto fracturedFaceId = GetFaceId(mFractureNodeIdx, edgeId0[0], edgeId0[1]);
						auto faceNeighborId = FindFaceNeighbor(mIdToTetrahedra, fracturingTetId, fracturedFaceId);
						if (faceNeighborId != -1)
						{
							const auto& faceNeighborTet = mIdToTetrahedra[faceNeighborId];
							mNewTetrahedra.push_back(Tetrahedra(posFractureNode, newEdgeNode0, commonNode, faceNeighborTet.GetOtherVertex({ (size_t)mFractureNodeIdx, (size_t)edgeId0[0], (size_t)edgeId0[1] })));
							mNewTetrahedra.push_back(Tetrahedra(negFractureNode, newEdgeNode0, edge0Node, faceNeighborTet.GetOtherVertex( { (size_t)mFractureNodeIdx, (size_t)edgeId0[0], (size_t)edgeId0[1] })));
							mTetrahedraIdsToDelete.insert(faceNeighborId);
                            nonEdgeNeighbors.insert(faceNeighborId);
						}
					}
					// the second face neighbor, containing the fracture node, and both nodes from the other split edge
					{
						auto fracturedFaceId = GetFaceId(mFractureNodeIdx, edgeId1[0], edgeId1[1]);
						auto faceNeighborId = FindFaceNeighbor(mIdToTetrahedra, fracturingTetId, fracturedFaceId);
						if (faceNeighborId != -1)
						{
							const auto& faceNeighborTet = mIdToTetrahedra[faceNeighborId];
							mNewTetrahedra.push_back(Tetrahedra(posFractureNode, newEdgeNode1, commonNode, faceNeighborTet.GetOtherVertex({ (size_t)mFractureNodeIdx, (size_t)edgeId1[0], (size_t)edgeId1[1] })));
							mNewTetrahedra.push_back(Tetrahedra(negFractureNode, newEdgeNode1, edge1Node, faceNeighborTet.GetOtherVertex( { (size_t)mFractureNodeIdx, (size_t)edgeId1[0], (size_t)edgeId1[1] })));
							mTetrahedraIdsToDelete.insert(faceNeighborId);
                            nonEdgeNeighbors.insert(faceNeighborId);
						}
					}
					// the third face neighbor, containing the nodes that are not the fracture node
					{
						auto fracturedFaceId = GetFaceId(commonNode, edge0Node, edge1Node);
						auto faceNeighborId = FindFaceNeighbor(mIdToTetrahedra, fracturingTetId, fracturedFaceId);
						if (faceNeighborId != -1)
						{
							const auto& faceNeighborTet = mIdToTetrahedra[faceNeighborId];
							auto otherNode = faceNeighborTet.GetOtherVertex({ (size_t)commonNode, (size_t)edge0Node, (size_t)edge1Node });
							mNewTetrahedra.push_back(Tetrahedra(commonNode, newEdgeNode0, newEdgeNode1, otherNode));
							// newEdgeNode0 and edge1Node form the diagonal
							mNewTetrahedra.push_back(Tetrahedra(edge0Node, newEdgeNode0, edge1Node, otherNode));
							mNewTetrahedra.push_back(Tetrahedra(newEdgeNode1, newEdgeNode0, edge1Node, otherNode));
							mTetrahedraIdsToDelete.insert(faceNeighborId);
                            nonEdgeNeighbors.insert(faceNeighborId);
						}
					}

                    // TODO think about whether we can remove the separate nodes below in favor of using the previously defined nodes above
                    {
                        const auto& edgeNeighbors = FindEdgeNeighbors(mIdToTetrahedra, edgeId0);

                        for (size_t neighborId : edgeNeighbors)
                        {
                            if (nonEdgeNeighbors.find(neighborId) != nonEdgeNeighbors.end())
                                continue;

                            const auto& edgeNeighorTet = mIdToTetrahedra[neighborId];
                            const auto& remainingNodes = edgeNeighorTet.GetOtherVertices(edgeId0[0], edgeId0[1]);
                            size_t posEdgeNode;
                            size_t negEdgeNode;
                            SeparateNodes(posEdgeNode, negEdgeNode, { (size_t)edgeId0[0], (size_t)edgeId0[1] }, mFractureNodePosition, mFracturePlaneNormal);

                            mTetrahedraIdsToDelete.insert(neighborId);

                            mNewTetrahedra.push_back(Tetrahedra(edgeVertex0, remainingNodes[0], remainingNodes[1], posEdgeNode));
                            mNewTetrahedra.push_back(Tetrahedra(edgeVertex0, remainingNodes[0], remainingNodes[1], negEdgeNode));
                        }
                    }
                    {
                        const auto& edgeNeighbors = FindEdgeNeighbors(mIdToTetrahedra, edgeId1);

                        for (size_t neighborId : edgeNeighbors)
                        {
                            if (nonEdgeNeighbors.find(neighborId) != nonEdgeNeighbors.end())
                                continue;

                            const auto& edgeNeighorTet = mIdToTetrahedra[neighborId];
                            const auto& remainingNodes = edgeNeighorTet.GetOtherVertices(edgeId1[0], edgeId1[1]);
                            size_t posEdgeNode;
                            size_t negEdgeNode;
                            SeparateNodes(posEdgeNode, negEdgeNode, { (size_t)edgeId1[0], (size_t)edgeId1[1] }, mFractureNodePosition, mFracturePlaneNormal);

                            mTetrahedraIdsToDelete.insert(neighborId);

                            mNewTetrahedra.push_back(Tetrahedra(edgeVertex1, remainingNodes[0], remainingNodes[1], posEdgeNode));
                            mNewTetrahedra.push_back(Tetrahedra(edgeVertex1, remainingNodes[0], remainingNodes[1], negEdgeNode));
                        }
                    }
                }
                else
                {
                    // tetrahedra is small
                }
            }

            std::cout << "Num tet to delete : " << mTetrahedraIdsToDelete.size() << std::endl;
            std::cout << "Num tet to add : " << mNewTetrahedra.size() << std::endl;

            for (auto tetId : mTetrahedraIdsToDelete)
            {
                mIdToTetrahedra.erase(tetId);
                neighborSet.erase(tetId);
            }

			for (auto& tet : mNewTetrahedra)
				mIdToTetrahedra[mTetIdCounter++] = tet;

            // verify that all nodes are referenced by at least one tet
            {
                std::unordered_set<size_t> vertexIds;

                for (const auto& pair : mIdToTetrahedra)
                {
                    for (const auto& vertexId : pair.second.mIndices)
                        vertexIds.insert(vertexId);
                }

                for (size_t id = 0; id < mVertices.size(); id++)
                {
                    if (vertexIds.find(id) == vertexIds.end())
                    {
						std::cout << "Vertex id : " << id << std::endl;
						throw std::exception("Isolated vertex detected.");
                    }
                }
            }

            mTetrahedraIdsToDelete.clear();
            mNewTetrahedra.clear();
		}
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
        {
            std::cout << "Num Pos : " << numPos << " Num Neg : " << numNeg << std::endl;
            std::cout << "d0 : " << d0 << " d1 : " << d1 << std::endl;
            throw std::exception("Expected one positive and one negative vertex in isolated edge fracture.");
        }

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
		gte::Vector3<double> normal
		{
			planeNormal.x,
			planeNormal.y,
			planeNormal.z
		};
		gte::Vector3<double> planeOrigin
		{
		   planePos.x,
		   planePos.y,
		   planePos.z
		};
		gte::Plane3<double> plane(normal, planeOrigin);

		auto dir = glm::normalize(edgePos1 - edgePos0);

		gte::Vector3<double> lineDir
		{
			dir.x, dir.y, dir.z
		};

		gte::Vector3<double> lineOrigin
		{
			edgePos0.x, edgePos0.y, edgePos0.z
		};

		gte::Line3<double> line(lineOrigin, lineDir);

		// Get signed distance of all vertices to plane
		gte::FIQuery<double, gte::Line3<double>, gte::Plane3<double>> findIntersectionQuery;
		auto results = findIntersectionQuery(line, plane);

		if (!results.intersect)
			return false;

        if (isnan(results.parameter))
            throw std::exception("Nan detected in plane-edge intersection.");

        std::cout << "Result parameter : " << results.parameter << std::endl;

        double param = results.parameter / (glm::distance(edgePos0, edgePos1));

        std::cout << "Scaled Result parameter : " << param << std::endl;

        if (param <= 0 || param >= 1.0)
            return false;

		d = param;

		if (intersectionPos)
		{
			intersectionPos->x = results.point[0];
			intersectionPos->y = results.point[1];
			intersectionPos->z = results.point[2];
		}

		return true;
    }

    ////////////////////////////////////////////////////////////////////////////////

    bool FractureContext::IsIsolatedEdge(const glm::ivec2& edgeId) const
    {
        size_t counter = 0;

        for (const auto& pair : mIdToTetrahedra)
        {
            const auto& tet = pair.second;
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

        for (const auto & pair : mIdToTetrahedra)
        {
            size_t tetId = pair.first;
            const auto& tetrahedra = pair.second;
            if (tetrahedra.ContainsVertexIndex(nodeIdx))
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
}