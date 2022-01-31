#include "TestDeformation.hpp"
#include "dsyevq3.h"
#include "dsyevh3.h"

#include <Mathematics/Delaunay3.h>

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
    }

    void ComputeMMat(const vec3& eigenVector, mat3& outMat)
    {
        if (glm::length(eigenVector) < 1e-6)
        {
            outMat = mat3(0.0);
            return;
        }

        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                outMat[j][k] = eigenVector[j] * eigenVector[k];

        outMat /= glm::length(eigenVector);
    }

    ////////////////////////////////////////////////////////////////////////////////

	void TetraGroup::Update(float timestep)
	{
		mat4x3 pMat;
		mat4x3 vMat;
        mat3 sigmaPlus = mat3(0);
        mat3 sigmaMinus = mat3(0);
        mat3 mMat;
        double sigma[3][3];
        double eigenVectors[3][3];
        double eigenValues[3];
        double separation[3][3];

        for (auto& vertex : mVertices)
        {
            vertex.mCompressiveForces.clear();
            vertex.mTensileForces.clear();
            vertex.mForce = vec3(0.0);
            vertex.mLargestEigenvalue = 0.0f;
            vertex.mPrincipalEigenVector = vec3(0.0f);
        }

		for (auto& tetrahedra : mTetrahedra)
		{
            auto& v0 = mVertices[tetrahedra.mIndices[0]];
            auto& v1 = mVertices[tetrahedra.mIndices[1]];
            auto& v2 = mVertices[tetrahedra.mIndices[2]];
            auto& v3 = mVertices[tetrahedra.mIndices[3]];

			pMat = mat4x3(v0.mPosition,
				v1.mPosition,
				v2.mPosition,
				v3.mPosition);

			vMat = mat4x3(v0.mVelocity,
				v1.mVelocity,
				v2.mVelocity,
				v3.mVelocity);

			mat3x4 pBeta = pMat * tetrahedra.mBeta;

            vec3 dx_u1 = pBeta * vec4(1, 0, 0, 0);
            vec3 dx_u2 = pBeta * vec4(0, 1, 0, 0);
            vec3 dx_u3 = pBeta * vec4(0, 0, 1, 0);

            std::array<vec3, 3> dx{ dx_u1, dx_u2, dx_u3 };

            mat3x4 vBeta = vMat * tetrahedra.mBeta;
            vec3 dxd_u1 = vBeta * vec4(1, 0, 0, 0);
            vec3 dxd_u2 = vBeta * vec4(0, 1, 0, 0);
            vec3 dxd_u3 = vBeta * vec4(0, 0, 1, 0);

            std::array<vec3, 3 > dxd{ dxd_u1, dxd_u2, dxd_u3 };
            std::array<glm::vec3, 3> dots;

            mat3 strainTensor;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    strainTensor[i][j] = dot(dx[i], dx[j]) - (i == j ? 1 : 0);

            mat3 rateOfStrainTensor;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    rateOfStrainTensor[i][j] = dot(dx[i], dxd[j]) + dot(dxd[i], dx[j]);

            mat3 elasticStress;
            float strainTrace = strainTensor[0][0] + strainTensor[1][1] + strainTensor[2][2];
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    elasticStress[i][j] = mLambda * strainTrace * (i == j ? 1 : 0) + 2 * mMu * strainTensor[i][j];

            mat3 viscousStress;
            float rateTrace = rateOfStrainTensor[0][0] + rateOfStrainTensor[1][1] + rateOfStrainTensor[2][2];
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    viscousStress[i][j] = mPhi * rateTrace * (i == j ? 1 : 0) + 2 * mPsi * rateOfStrainTensor[i][j];

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
                vec3 forceOnNode = vec3(0);
                for (int j = 0; j < 4; j++)
                {
                    float innerProduct = 0;
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

                if (isnan(forceOnNode.x) || isnan(forceOnNode.x) || isnan(forceOnNode.x))
                {
                    throw std::exception("Nan value in forceOnNode detected.");
                }

                mVertices[tetrahedra.mIndices[i]].mForce += forceOnNode;
            }

            dsyevh3(sigma, eigenVectors, eigenValues);

            sigmaPlus = mat3(0.0);
            sigmaMinus = mat3(0.0);

            for (int i = 0; i < 3; i++)
            {
                // Not completely sure that this is right (might be [i][0], etc)
                dvec3 eigenVec(eigenVectors[0][i], eigenVectors[1][i], eigenVectors[2][i]);
                mat3 mMat;
                ComputeMMat(eigenVec, mMat);

                sigmaPlus += fmax(0.0f, (float)eigenValues[i]) * mMat;
                sigmaMinus += fmin(0.0f, (float)eigenValues[i]) * mMat;
            }

            std::array<vec3, 4> fPlus = { vec3(0.0),vec3(0.0),vec3(0.0),vec3(0.0) };
            std::array<vec3, 4> fMinus = { vec3(0.0),vec3(0.0),vec3(0.0),vec3(0.0) };;

            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    float innerProductPlus = 0;
                    float innerProductMinus = 0;
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

                fPlus[i] *= tetrahedra.mVolume * 0.5f;
                fMinus[i] *= tetrahedra.mVolume * 0.5f;
            }

            for (int i = 0; i < 4; i++)
            {
                mVertices[tetrahedra.mIndices[i]].mCompressiveForces.push_back(fMinus[i]);
                mVertices[tetrahedra.mIndices[i]].mTensileForces.push_back(fPlus[i]);
            }
		}

        vec3 a;
        mat3 mSeparation;
        mat3 m_Mat;
        vec3 fPlus;
        vec3 fMinus;
        
        for (size_t vertexIdx = 0; vertexIdx < mVertices.size(); vertexIdx++)
        {
            auto& vertex = mVertices[vertexIdx];
            vertex.mForce += glm::vec3(0, -9.8, 0) * vertex.mMass;

            a = vertex.mInvMass * vertex.mForce;
            vertex.mVelocity += a * timestep;
            vertex.mPosition += vertex.mVelocity * timestep;

            fMinus = vec3(0.0);
            fPlus = vec3(0.0);

            for (const auto& compressiveForce : vertex.mCompressiveForces)
                fMinus += compressiveForce;
            for (const auto& tensileForce : vertex.mTensileForces)
                fPlus += tensileForce;

            mSeparation = mat3(0.0);

            ComputeMMat(fMinus, m_Mat);
            mSeparation += m_Mat;
            ComputeMMat(fPlus, m_Mat);
            mSeparation -= m_Mat;

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

            dsyevh3(separation, eigenVectors, eigenValues);

            // Look for largest eigenvalue and corresponding eigen vector
            float largestEigenvalue = (float)eigenValues[0];
            vec3 principalEigenVector = vec3(eigenVectors[0][0], eigenVectors[0][1], eigenVectors[0][2]);
            if (largestEigenvalue < eigenValues[1]) {
                largestEigenvalue = (float)eigenValues[1];
                principalEigenVector = vec3(eigenVectors[1][0], eigenVectors[1][1], eigenVectors[1][2]);
            }
            if (largestEigenvalue < eigenValues[2]) {
                largestEigenvalue = (float)eigenValues[2];
                principalEigenVector = vec3(eigenVectors[2][0], eigenVectors[2][1], eigenVectors[2][2]);
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
        }

        for (auto& vertex : mVertices)
        {
            if (vertex.mPosition.y < 0)
            {
                vertex.mPosition.y = -vertex.mPosition.y;
                float elasticity = 0.9f;
                float friction = 0.1f;
                const auto & velocity = vertex.mVelocity;
                vertex.mVelocity = (glm::vec3((1 - friction)* velocity.x, -elasticity * velocity.y, (1 - friction)* velocity.z));
            }
        }
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

            glm::mat4 m = glm::mat4(
                glm::vec4(m0, 1.f),
                glm::vec4(m1, 1.f),
                glm::vec4(m2, 1.f),
                glm::vec4(m3, 1.f)
            );

            tetrahedra.mBeta = glm::inverse(m);

            tetrahedra.mVolume = 1.0f / 6.0f * fabs(
                glm::dot(
                    glm::cross(
                        m1 - m0,
                        m2 - m0
                    ),
                    m3 - m0));

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
        std::cout << "Fracture Node" << std::endl;

        Vertex negativeCopiedVertex;
        negativeCopiedVertex.mPosition = mVertices[fractureNodeIdx].mPosition;
        negativeCopiedVertex.mMaterialCoordinates = mVertices[fractureNodeIdx].mMaterialCoordinates;
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
                }
                else if (d0 < 0 && d1 < 0 && d2 < 0 && d3 < 0)
                {
                    element.ReplaceVertex(fractureNodeIdx, negativeCopiedVertexIdx);
                }
                else
                    std::cout << "Not edges were split but the tetrahedra has vertices on both sides of the fracture plane." << std::endl;
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

            ComputeDerivedQuantities();
        }
    }
    
    ////////////////////////////////////////////////////////////////////////////////

    size_t TetraGroup::SplitEdge(const glm::ivec2 & edgeIdx, size_t fractureNodeIdx, const glm::vec3& planeNormal)
    {
        const auto& planePosition = mVertices[fractureNodeIdx].mPosition;
        const auto& edgePos0 = mVertices[edgeIdx.x].mPosition;
        const auto& edgePos1 = mVertices[edgeIdx.y].mPosition;

        glm::vec3 intersectionPos;
        if (!PlaneIntersectEdge(planePosition, planeNormal, edgePos0, edgePos1, &intersectionPos))
            std::cout << "Expected edge to be split by plane but the two don't intersect." << std::endl;

        Vertex vertex;
        vertex.mPosition = intersectionPos;

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

    bool TetraGroup::EdgeIntersectsPlane(const glm::ivec2& edgeIdx, int fracutreNodeIdx, const glm::vec3& planeNormal) const
    {
        const auto& planePosition = mVertices[fracutreNodeIdx].mPosition;
        const auto& edgePos0 = mVertices[edgeIdx.x].mPosition;
        const auto& edgePos1 = mVertices[edgeIdx.y].mPosition;
        return PlaneIntersectEdge(planePosition, planeNormal, edgePos0, edgePos1);
    }

    ////////////////////////////////////////////////////////////////////////////////

    bool TetraGroup::PlaneIntersectEdge(const glm::vec3& planePos, const glm::vec3& planeNormal, const glm::vec3& edgePos0, const glm::vec3& edgePos1, glm::vec3* intersectionPos) const
    {
        // plane equation is (p - p0) * n = 0, where n is the normal vector
        // line equation is p = l0 + v * t, where v is the direction vector of the line
        // compute d = (p0 - l0) * n / (l * n)
        // intersection is at p = l0 + v * d
        // if d == 0, plane contains line
        auto edgeDirVec = glm::normalize(edgePos1 - edgePos0);

        float d = glm::dot((planePos - edgePos0), planeNormal) / glm::dot(edgeDirVec, planeNormal);
        
        // Outside the line segment
        if (d <= 0 || d > 1)
            return false;

        if (intersectionPos)
            *intersectionPos = edgePos0 + d * edgeDirVec;

        return true;
    }

    ////////////////////////////////////////////////////////////////////////////////

    float TetraGroup::GetSignedDistanceToPlane(int nodeIdx, int fractureNodeIdx, const glm::vec3& planeNormal) const
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
}