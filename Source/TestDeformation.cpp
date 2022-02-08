#include "TestDeformation.hpp"
#include "dsyevq3.h"
#include "dsyevh3.h"

#include <Mathematics/Delaunay3.h>
#include <Mathematics/SymmetricEigensolver3x3.h>
#include <Mathematics/DistPointHyperplane.h>

#include <iostream>
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
            return glm::ivec3();
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

    namespace 
    {
        IronGames::Vector3 Convert(const glm::dvec3& vec)
        {
            IronGames::Vector3 v;
            v.set_x(vec.x);
            v.set_y(vec.y);
            v.set_z(vec.z);
            return v;
        }

        IronGames::Matrix3 Convert(const glm::dmat3& mat)
        {
            IronGames::Matrix3 m;
            *m.mutable_col0() = Convert(mat[0]);
            *m.mutable_col1() = Convert(mat[1]);
            *m.mutable_col2() = Convert(mat[2]);
            return m;
        }
    }

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

        auto* frame = (mSummary ? mSummary->add_frames() : nullptr);

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
                *vert->mutable_position() = Convert(vertex.mPosition);
                *vert->mutable_material_coordinates() = Convert(vertex.mMaterialCoordinates);
                *vert->mutable_velocity() = Convert(vertex.mVelocity);
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
                *tet.mutable_strain_tensor() = Convert(strainTensor);
                *tet.mutable_stress_tensor() = Convert(sigma);
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

            if (largestEigenvalue > mToughness)
                std::cout << " Large eigen " << std::endl;

            mVertices[vertexIdx].mPrincipalEigenVector = principalEigenVector;
            mVertices[vertexIdx].mLargestEigenvalue = largestEigenvalue;

            if (saveFrame)
            {
                auto& vert = *frame->mutable_vertices(vertexIdx);
                *vert.mutable_force() = Convert(mVertices[vertexIdx].mForce);
                vert.set_largest_eigenvalue(mVertices[vertexIdx].mLargestEigenvalue);
                *vert.mutable_principal_eigenvector() = Convert(mVertices[vertexIdx].mPrincipalEigenVector);

                for (const auto& compressiveForce : mVertices[vertexIdx].mCompressiveForces)
                    *vert.add_compressive_forces() = Convert(compressiveForce);
                for (const auto& tensileForce : mVertices[vertexIdx].mTensileForces)
                    *vert.add_tensile_forces() = Convert(tensileForce);

                *vert.mutable_separation_tensor() = Convert(mSeparation);
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

        Vertex negativeCopiedVertex;
        negativeCopiedVertex.mPosition = mVertices[fractureNodeIdx].mPosition;
        negativeCopiedVertex.mMaterialCoordinates = mVertices[fractureNodeIdx].mMaterialCoordinates;
        negativeCopiedVertex.mVelocity = mVertices[fractureNodeIdx].mVelocity;
        mVertices.push_back(negativeCopiedVertex);

        size_t negativeCopiedVertexIdx = mVertices.size() - 1;

        // so, n+ is fractureNodeIndex and n- is negativeCopiedVertexIdx

        std::unordered_map<glm::ivec2, int> splitEdgeToNewVertices;
        std::vector<int> elementVertexIndices;
        std::vector<int> topVertexIndices;
        std::vector<int> botVertexIndices;

        std::vector<int> tetrahedraToRemove;
        std::vector<int> neighborTetrahedraToRemove;
        std::vector<Tetrahedra> tetrahedraToAdd;

        for (size_t elementIdx : GetTetrahedrasFromNode(fractureNodeIdx))
        {
            auto& element = mTetrahedra[elementIdx];

            std::vector<glm::ivec2> splitEdges;

            for (const auto& edge : element.GetEdges())
            {
                auto iter = splitEdgeToNewVertices.find(edge);
                if (iter != splitEdgeToNewVertices.end())
                {
                    splitEdges.push_back(edge);
                }
                else if (EdgeIntersectsPlane(edge, fractureNodeIdx, fracturePlaneNormal))
                {
                    auto newVertexId = SplitEdge(edge, fractureNodeIdx, fracturePlaneNormal);
                    splitEdgeToNewVertices[edge] = newVertexId;
                    splitEdges.push_back(edge);
                }
            }

            if (splitEdges.size() == 0)
            {
                // we need to assign this tet to n+ or n-
                auto d0 = GetSignedDistanceToPlane(element.mIndices[0], fractureNodeIdx, fracturePlaneNormal);
                auto d1 = GetSignedDistanceToPlane(element.mIndices[1], fractureNodeIdx, fracturePlaneNormal);
                auto d2 = GetSignedDistanceToPlane(element.mIndices[2], fractureNodeIdx, fracturePlaneNormal);
                auto d3 = GetSignedDistanceToPlane(element.mIndices[3], fractureNodeIdx, fracturePlaneNormal);

                if (d0 >= 0 && d1 >= 0 && d2 >= 0 && d3 >= 0)
                {
                    // It gets the positive one, which it already has
                    continue;
                }
                else if (d0 <= 0 && d1 <= 0 && d2 <= 0 && d3 <= 0)
                {
                    element.ReplaceVertex(fractureNodeIdx, negativeCopiedVertexIdx);
                    continue;
                }
                else
                    std::cout << "No edges were split but the tetrahedra has vertices on both sides of the fracture plane." << std::endl;
            }
            else if (splitEdges.size() == 1)
            {
                // TODO edge case == 1 
                // This can only happen for neighbors, or the edge case where the plane is along an edge?
            }

            if (splitEdges.size() != 2)
            {
                std::cout << "Expected 2 split edges, but found " << splitEdges.size() << "." << std::endl;
                continue;
            }

            elementVertexIndices.clear();
            topVertexIndices.clear();
            botVertexIndices.clear();

            tetrahedraToRemove.push_back(elementIdx);
            for (const auto& idx : element.mIndices)
            {
                if (idx != fractureNodeIdx)
					elementVertexIndices.push_back(idx);
            }
            
            for (const auto& vertexIdx : elementVertexIndices)
            {
                if (GetSignedDistanceToPlane(vertexIdx, fractureNodeIdx, fracturePlaneNormal) < 0.0f)
                    botVertexIndices.push_back(vertexIdx);
                else
                    topVertexIndices.push_back(vertexIdx);
            }

            if (topVertexIndices.size() == 1 && botVertexIndices.size() == 2)
            {
                auto a = topVertexIndices[0];
                auto b = botVertexIndices[0];
                auto c = botVertexIndices[1];
                auto dp = fractureNodeIdx;
                auto dm = negativeCopiedVertexIdx;
                auto e = splitEdgeToNewVertices[splitEdges[0]];
                auto f = splitEdgeToNewVertices[splitEdges[1]];

                tetrahedraToAdd.push_back(Tetrahedra(a, dp, e, f));
                tetrahedraToAdd.push_back(Tetrahedra(b, f, e, dm));
                tetrahedraToAdd.push_back(Tetrahedra(b, c, f, dm));
            }
            else if (topVertexIndices.size() == 2 && botVertexIndices.size() == 1)
            {
                auto a = botVertexIndices[0];
                auto b = topVertexIndices[0];
                auto c = topVertexIndices[1];
                auto dp = fractureNodeIdx;
                auto dm = negativeCopiedVertexIdx;
                auto e = splitEdgeToNewVertices[splitEdges[0]];
                auto f = splitEdgeToNewVertices[splitEdges[1]];

                tetrahedraToAdd.push_back(Tetrahedra(a, dm, e, f));
                tetrahedraToAdd.push_back(Tetrahedra(b, f, e, dp));
                tetrahedraToAdd.push_back(Tetrahedra(b, c, f, dp));
            }
            else
            {
                std::cout << "Unexpected number of vertices above (" << topVertexIndices.size() << "), and below (" << botVertexIndices.size() << ")." << std::endl;
            }
        }

        std::vector<int> splitNeighbors;

        for (auto splitTetIdx : tetrahedraToRemove)
        {
            for (auto neighborTetIdx : GetTetrahedraNeighbors(splitTetIdx))
            {
                const auto& neighborTet = mTetrahedra[neighborTetIdx];
                if (neighborTet.ContainsVertexIndex(fractureNodeIdx))
                    continue;

                // Also skip neighbors we have already split
                auto iter = std::find(splitNeighbors.begin(), splitNeighbors.end(), neighborTetIdx);
                if (iter != splitNeighbors.end())
                    continue;

                const auto& edges = neighborTet.GetEdges();

                std::vector<glm::ivec2> splitNeighborEdges;
                for (const auto& edge : edges)
                {
                    if (splitEdgeToNewVertices.find(edge) != splitEdgeToNewVertices.end())
                        splitNeighborEdges.push_back(edge);
                }

                if (splitNeighborEdges.size() == 0)
                {
                    // Nothing to do here
                }
                else if (splitNeighborEdges.size() == 1)
                {
                    const auto& edge = splitNeighborEdges[0];
                    const auto& newNode = splitEdgeToNewVertices[edge];

                    std::vector<int> sharedNodes;
                    for (auto nodeIdx : neighborTet.mIndices)
                    {
                        if (nodeIdx == edge.x || nodeIdx == edge.y)
                            continue;
                        
                        sharedNodes.push_back(nodeIdx);
                    }

                    // actually, this neighbor may be a face neighbor
                    // so we are definitely missing that case here
                    if (sharedNodes.size() != 2)
                        std::cout << "Expected 2 shared nodes during neighbor splitting, found " << sharedNodes.size() << "." << std::endl;

                    auto a = sharedNodes[0];
                    auto b = sharedNodes[1];
                    auto c = edge.x;
                    auto d = edge.y;
                    auto e = newNode;

                    tetrahedraToAdd.push_back(Tetrahedra(a, c, b, e));
                    tetrahedraToAdd.push_back(Tetrahedra(a, b, d, e));

                    neighborTetrahedraToRemove.push_back(neighborTetIdx);
                }
                // this is a face neighbor
                else if (splitNeighborEdges.size() == 2)
                {
                    auto a = GetCommonVertexFromEdges(splitNeighborEdges[0], splitNeighborEdges[1]);
                    auto e = splitEdgeToNewVertices[splitNeighborEdges[0]];
                    auto f = splitEdgeToNewVertices[splitNeighborEdges[1]];

                    std::vector<int> remainingVertices;

                    for (auto vert : neighborTet.mIndices)
                    {
                        if (vert == a)
                            continue;
                        if (vert == splitNeighborEdges[0].x || vert == splitNeighborEdges[0].y)
                            continue;
                        if (vert == splitNeighborEdges[1].x || vert == splitNeighborEdges[1].y)
                            continue;

                        remainingVertices.push_back(vert);
                    }

                    if (remainingVertices.size() != 1)
                        std::cout << "Expected one remaining vertex in neighbor split but found " << remainingVertices.size() << "." << std::endl;

                    auto b = remainingVertices[0];

                    std::vector<int> cdVertices;
                    for (auto vert : neighborTet.mIndices)
                    {
                        if (vert == a || vert == b)
                            continue;
                        cdVertices.push_back(vert);
                    }

                    if (cdVertices.size() == 2)
                        std::cout << " Expected two remaining vertices in neighbor split but found " << cdVertices.size() << "." << std::endl;

                    auto c = cdVertices[0];
                    auto d = cdVertices[1];

                    std::vector<glm::ivec2> diagonalEdges;

                    if (splitNeighborEdges[0].x == c || splitNeighborEdges[0].y == c)
                        diagonalEdges.push_back(GetEdgeId(c, f));
                    else
                        diagonalEdges.push_back(GetEdgeId(c, e));

                    if (splitNeighborEdges[0].x == d || splitNeighborEdges[0].y == d)
                        diagonalEdges.push_back(GetEdgeId(d, f));
                    else
                        diagonalEdges.push_back(GetEdgeId(d, e));

                    if (diagonalEdges.size() != 2)
                        std::cout << "Expected 2 diagonal edges during neighbor split, found " << diagonalEdges.size() << "." << std::endl;

                    int diagonalEdgeIdx = -1;

                    for (const auto& newTet : tetrahedraToAdd)
                    {
                        if (newTet.ContainsEdgeIndex(diagonalEdges[0]))
                        {
                            diagonalEdgeIdx = 0;
                            break;
                        }
                        else if (newTet.ContainsEdgeIndex(diagonalEdges[1]))
                        {
                            diagonalEdgeIdx = 1;
                            break;
                        }
                    }

                    auto diag = diagonalEdges[diagonalEdgeIdx];

                    std::vector<int> offDiagNodes;

                    if (c != diag.x && c != diag.y)
                        offDiagNodes.push_back(c);
                    if (d != diag.x && d != diag.y)
                        offDiagNodes.push_back(d);
                    if (e != diag.x && e != diag.y)
                        offDiagNodes.push_back(e);
                    if (f != diag.x && f != diag.y)
                        offDiagNodes.push_back(f);

                    if (offDiagNodes.size() != 2)
                        std::cout << "Incorrect number of off diagonal nodes, expected 2, found " << offDiagNodes.size() << "." << std::endl;

                    tetrahedraToAdd.push_back(Tetrahedra(a, e, b, f));
                    tetrahedraToAdd.push_back(Tetrahedra(b, diag.y, diag.x, offDiagNodes[0]));
                    tetrahedraToAdd.push_back(Tetrahedra(b, diag.x, diag.y, offDiagNodes[1]));

                    neighborTetrahedraToRemove.push_back(neighborTetIdx);
                }
                else
                {
                    std::cout << "Too many neighbor edges were split by the fracture plane, num split = " << splitNeighborEdges.size() << "." << std::endl;
                }

                splitNeighbors.push_back(neighborTetIdx);
            }
        }

        {
		    std::vector<int> tetIndicesToRemove;
            for (auto tetIdx : tetrahedraToRemove)
                tetIndicesToRemove.push_back(tetIdx);
            for (auto neighborTetIdx : neighborTetrahedraToRemove)
                tetIndicesToRemove.push_back(neighborTetIdx);

            std::sort(tetIndicesToRemove.begin(), tetIndicesToRemove.end(), std::greater<int>());

            for (auto idx : tetIndicesToRemove)
            {
                mTetrahedra.erase(mTetrahedra.begin() + idx);
            }
        }

        // TODO - optimize, recomputing quantities for all tet instead of just those that were modified.
        // - Don't currently have a vertex to tetrahedra correspondence so redoing the mass would require
        //   a linear search.
        {
            for (auto& tet : tetrahedraToAdd)
            {
                mTetrahedra.push_back(tet);
            }

        }
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

    size_t TetraGroup::GetCommonVertexFromEdges(const glm::ivec2& edge0, const glm::ivec2& edge1) const
    {
        if (edge0.x == edge1.x || edge0.x == edge1.y)
            return edge0.x;
        if (edge0.y == edge1.x || edge0.y == edge1.y)
            return edge0.y;

        std::cout << "Failed to find common vertices for provided edges." << std::endl;
        return 0;
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
            return;

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
                if (mTetrahedra[tetIdx].mCreatedThisFrame)
                    continue;

                fractured = true;

                FractureTetrahedra(tetIdx);

                // deletion in descending order so we don't invalidate the indices
                std::sort(mTetrahedraToDelete.begin(), mTetrahedraToDelete.end(), std::greater<size_t>());
                for (size_t idx = 0; idx < mTetrahedraToDelete.size(); idx++)
                    mTetrahedra.erase(mTetrahedra.begin() + mTetrahedraToDelete[idx]);

                for (auto& tet : mNewTetrahedra)
                    mTetrahedra.push_back(tet);

                mTetrahedraToDelete.clear();
                mNewTetrahedra.clear();
            }
        }

        for (auto& tet : mTetrahedra)
            tet.mCreatedThisFrame = false;
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
            return results.distance;
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

        if (IsIsolatedEdge(GetEdgeId(otherNodes[0], otherNodes[1])))
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
        NeighborFaceFracture(GetFaceId(fracturedFaceIds[0], fracturedFaceIds[1], fracturedFaceIds[2]), { newEdgeNode });
        NeighborEdgeFracture();
    }

    ////////////////////////////////////////////////////////////////////////////////

    void FractureContext::FractureTetrahedra(size_t tetIdx)
    {
        mTetrahedraIdx = tetIdx;

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

            if (dist >= 0.0)
                numPosition++;
            if (dist <= 0.0)
                numNegative++;
        }

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

        if (IsIsolatedEdge(GetEdgeId(mFractureNodeIdx, edgeVertex0)))
            em = CloneVertex(ep);
        
        if (IsIsolatedEdge(GetEdgeId(mFractureNodeIdx, edgeVertex1)))
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

        NeighborFaceFracture(GetFaceId(otherVertices[0], otherVertices[1], otherVertices[2]), { edgeVertex0, edgeVertex1 });
        NeighborEdgeFracture(edgeIds[0]);
        NeighborEdgeFracture(edgeIds[1]);
    }

    ////////////////////////////////////////////////////////////////////////////////

    void FractureContext::NeighborFaceFracture(const glm::ivec3& faceId, const std::vector<size_t>& newNodeIds)
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

        if (newNodeIds.size() == 1)
        {
            // split neighborTetId into two tetrahedra using it's 4 vertices and
            // the single newNodeId
            return;
        }
        
        if (newNodeIds.size() == 2)
        {
            // split neighborTetId into three tetrahedra using it's 4 vertices and
            // the two newNodeIds

            {
                /*
                
                                    auto a = GetCommonVertexFromEdges(splitNeighborEdges[0], splitNeighborEdges[1]);
                    auto e = splitEdgeToNewVertices[splitNeighborEdges[0]];
                    auto f = splitEdgeToNewVertices[splitNeighborEdges[1]];

                    std::vector<int> remainingVertices;

                    for (auto vert : neighborTet.mIndices)
                    {
                        if (vert == a)
                            continue;
                        if (vert == splitNeighborEdges[0].x || vert == splitNeighborEdges[0].y)
                            continue;
                        if (vert == splitNeighborEdges[1].x || vert == splitNeighborEdges[1].y)
                            continue;

                        remainingVertices.push_back(vert);
                    }

                    if (remainingVertices.size() != 1)
                        std::cout << "Expected one remaining vertex in neighbor split but found " << remainingVertices.size() << "." << std::endl;

                    auto b = remainingVertices[0];

                    std::vector<int> cdVertices;
                    for (auto vert : neighborTet.mIndices)
                    {
                        if (vert == a || vert == b)
                            continue;
                        cdVertices.push_back(vert);
                    }

                    if (cdVertices.size() == 2)
                        std::cout << " Expected two remaining vertices in neighbor split but found " << cdVertices.size() << "." << std::endl;

                    auto c = cdVertices[0];
                    auto d = cdVertices[1];

                    std::vector<glm::ivec2> diagonalEdges;

                    if (splitNeighborEdges[0].x == c || splitNeighborEdges[0].y == c)
                        diagonalEdges.push_back(GetEdgeId(c, f));
                    else
                        diagonalEdges.push_back(GetEdgeId(c, e));

                    if (splitNeighborEdges[0].x == d || splitNeighborEdges[0].y == d)
                        diagonalEdges.push_back(GetEdgeId(d, f));
                    else
                        diagonalEdges.push_back(GetEdgeId(d, e));

                    if (diagonalEdges.size() != 2)
                        std::cout << "Expected 2 diagonal edges during neighbor split, found " << diagonalEdges.size() << "." << std::endl;

                    int diagonalEdgeIdx = -1;

                    for (const auto& newTet : tetrahedraToAdd)
                    {
                        if (newTet.ContainsEdgeIndex(diagonalEdges[0]))
                        {
                            diagonalEdgeIdx = 0;
                            break;
                        }
                        else if (newTet.ContainsEdgeIndex(diagonalEdges[1]))
                        {
                            diagonalEdgeIdx = 1;
                            break;
                        }
                    }

                    auto diag = diagonalEdges[diagonalEdgeIdx];

                    std::vector<int> offDiagNodes;

                    if (c != diag.x && c != diag.y)
                        offDiagNodes.push_back(c);
                    if (d != diag.x && d != diag.y)
                        offDiagNodes.push_back(d);
                    if (e != diag.x && e != diag.y)
                        offDiagNodes.push_back(e);
                    if (f != diag.x && f != diag.y)
                        offDiagNodes.push_back(f);

                    if (offDiagNodes.size() != 2)
                        std::cout << "Incorrect number of off diagonal nodes, expected 2, found " << offDiagNodes.size() << "." << std::endl;

                    tetrahedraToAdd.push_back(Tetrahedra(a, e, b, f));
                    tetrahedraToAdd.push_back(Tetrahedra(b, diag.y, diag.x, offDiagNodes[0]));
                    tetrahedraToAdd.push_back(Tetrahedra(b, diag.x, diag.y, offDiagNodes[1]));

                    neighborTetrahedraToRemove.push_back(neighborTetIdx);
                */
            }

            // may need the id's of the split edges
            return;
        }

        throw std::exception("Attempted neighbor face fracture with invalid number of new node ids.");
    }

    ////////////////////////////////////////////////////////////////////////////////

    void FractureContext::NeighborEdgeFracture(const glm::ivec2& edgeId)
    {
        // need to get all tetrahedra that contain this edge and are not the original tet or the face tet
        std::vector<size_t> neighborTetrahedraIds;
        for (size_t neighborTetId = 0; neighborTetId < mTetrahedra.size(); neighborTetId++)
        {
            if (neighborTetId == mTetrahedraIdx || (mFracturedNeighborFaceTet && neighborTetId == mFracturedNeighborFaceTetId))
                continue;

            neighborTetrahedraIds.push_back(neighborTetId);
        }
        
        // need to think through the order of creation and deletion
        // - safe to create vertices whenever because it's additive
        // - safe to create tetrahedra whenever because it's additive
        // - not safe to delete tetrahedra, so the id's of tetrahedra that need to be deleted should be
        //   stored and deleted after... all elements are considered from a single fracture?
        // - let's say yes

        /*
                            const auto& edge = splitNeighborEdges[0];
                    const auto& newNode = splitEdgeToNewVertices[edge];

                    std::vector<int> sharedNodes;
                    for (auto nodeIdx : neighborTet.mIndices)
                    {
                        if (nodeIdx == edge.x || nodeIdx == edge.y)
                            continue;
                        
                        sharedNodes.push_back(nodeIdx);
                    }

                    // actually, this neighbor may be a face neighbor
                    // so we are definitely missing that case here
                    if (sharedNodes.size() != 2)
                        std::cout << "Expected 2 shared nodes during neighbor splitting, found " << sharedNodes.size() << "." << std::endl;

                    auto a = sharedNodes[0];
                    auto b = sharedNodes[1];
                    auto c = edge.x;
                    auto d = edge.y;
                    auto e = newNode;

                    tetrahedraToAdd.push_back(Tetrahedra(a, c, b, e));
                    tetrahedraToAdd.push_back(Tetrahedra(a, b, d, e));

                    neighborTetrahedraToRemove.push_back(neighborTetIdx);
        */

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

    std::vector<size_t> FractureContext::GetTetrahedraNeighbors(size_t tetrahedraIdx) const
    {
        return std::vector<size_t>();
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