#include "TestDeformation.hpp"
#include <glm/vec4.hpp>
#include <glm/common.hpp>
#include <iostream>

////////////////////////////////////////////////////////////////////////////////
//
extern "C" void __declspec(dllexport) __stdcall Deform(
	void* vertexPositions,
	int numVertices,
	void* tetrahedraIndices,
	int numTetrahedra,
	float timestep,
	float lambda,
	float psi,
	float mu,
	float phi,
	float density)
{
	auto positions = reinterpret_cast<float*>(vertexPositions);
	auto indices = reinterpret_cast<int*>(tetrahedraIndices);

	TestDeformation::TetraGroup group;
	group.mLambda = lambda;
	group.mMu = mu;
	group.mPhi = phi;
	group.mPsi = psi;
	group.mTetrahedra.resize(numTetrahedra);
	group.mVertices.resize(numVertices);

	for (int i = 0; i < numVertices; i++)
	{
		group.mVertices[i].mPosition = glm::vec3(positions[3 * i], positions[3 * i + 1], positions[3 * i + 2]);
		group.mVertices[i].mMaterialCoordinates = group.mVertices[i].mPosition;
		group.mVertices[i].mMass = 0.0f;
		group.mVertices[i].mInvMass = 0.0f;
	}

	for (int i = 0; i < numTetrahedra; i++)
	{
		auto& tetrahedra = group.mTetrahedra[i];
		tetrahedra.mIndices[0] = indices[4 * i];
		tetrahedra.mIndices[1] = indices[4 * i + 1];
		tetrahedra.mIndices[2] = indices[4 * i + 2];
		tetrahedra.mIndices[3] = indices[4 * i + 3];

		auto& m0 = group.mVertices[tetrahedra.mIndices[0]].mMaterialCoordinates;
		auto& m1 = group.mVertices[tetrahedra.mIndices[1]].mMaterialCoordinates;
		auto& m2 = group.mVertices[tetrahedra.mIndices[2]].mMaterialCoordinates;
		auto& m3 = group.mVertices[tetrahedra.mIndices[3]].mMaterialCoordinates;

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

		tetrahedra.mMass = density * tetrahedra.mVolume;

		for (int j = 0; j < 4; j++)
		{
			group.mVertices[tetrahedra.mIndices[j]].mMass += 0.25f * tetrahedra.mMass;
		}
	}

	for (int i = 0; i < numVertices; i++)
	{
		std::cout << "Mass : " << group.mVertices[i].mMass << std::endl;
		group.mVertices[i].mInvMass = 1.0f / group.mVertices[i].mMass;
	}

	// End initialization

	group.Update(timestep);

	// End simulation

	for (int i = 0; i < group.mVertices.size(); i++)
	{
		positions[3 * i] = group.mVertices[i].mPosition.x;
		positions[3 * i + 1] = group.mVertices[i].mPosition.y;
		positions[3 * i + 2 ] = group.mVertices[i].mPosition.z;
	}
}

#include <iostream>

int main(int argc, const char* argv[]) {

	float lambda = 2.65e6f;
	float psi = 397.f;
	float phi = 264.f;
	float mu = 3.97e6f;
	float density = 1013.f;
	float timestep = 0.00001f;

	TestDeformation::TetraGroup group;

	std::vector<float> positions = { 0, 0, 0, 1, 0.5f, 0, 0, 1, 0, 0, 0.5f, 1 };
	std::vector<int> indices = { 0, 1, 2, 3 };

	Deform(positions.data(), 4, indices.data(), 1, timestep, lambda, psi, mu, phi, density);

	for (int i = 0; i < 12; i++)
		std::cout << positions[i] << std::endl;
}

////////////////////////////////////////////////////////////////////////////////