#include "TestDeformation.hpp"
#include "dsyevq3.h"
#include "dsyevh3.h"

#include <glm/mat4x3.hpp>
#include <iostream>

using namespace glm;

namespace TestDeformation
{
	////////////////////////////////////////////////////////////////////////////////

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
        
        for (auto& vertex : mVertices)
        {
            vertex.mForce += glm::vec3(0, -9.8, 0) * vertex.mMass;

            a = vertex.mInvMass * vertex.mForce;
            vertex.mVelocity += a * timestep;
            vertex.mPosition += vertex.mVelocity * timestep;

            //std::cout << "Position " << vertex.mPosition[0] << ", " << vertex.mPosition[1] << ", " << vertex.mPosition[2] << std::endl;

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
}