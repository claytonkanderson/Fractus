#include "TestDeformation.hpp"
#include <glm/vec4.hpp>
#include <glm/common.hpp>
#include <iostream>

static void* mData = nullptr;

////////////////////////////////////////////////////////////////////////////////

extern "C" void __declspec(dllexport) __stdcall Initialize(
	void* vertexPositions,
	int numVertices,
	void* tetrahedraIndices,
	int numTetrahedra,
	float lambda,
	float psi,
	float mu,
	float phi,
	float density)
{
	TestDeformation::TetraGroup* group = new TestDeformation::TetraGroup();
	mData = group;

	group->mLambda = lambda;
	group->mPsi = psi;
	group->mPhi = phi;
	group->mMu = mu;
	group->mDensity = density;
	group->mVertices.resize(numVertices);
	group->mTetrahedra.resize(numTetrahedra);

	auto positions = reinterpret_cast<float*>(vertexPositions);
	auto indices = reinterpret_cast<int*>(tetrahedraIndices);

	for (int i = 0; i < numVertices; i++)
	{
		group->mVertices[i].mPosition = glm::vec3(positions[3 * i], positions[3 * i + 1], positions[3 * i + 2]);
		group->mVertices[i].mMaterialCoordinates = group->mVertices[i].mPosition;
		group->mVertices[i].mMass = 0.0f;
		group->mVertices[i].mInvMass = 0.0f;
	}

	for (int i = 0; i < numTetrahedra; i++)
	{
		auto& tetrahedra = group->mTetrahedra[i];
		tetrahedra.mIndices[0] = indices[4 * i];
		tetrahedra.mIndices[1] = indices[4 * i + 1];
		tetrahedra.mIndices[2] = indices[4 * i + 2];
		tetrahedra.mIndices[3] = indices[4 * i + 3];

		auto& m0 = group->mVertices[tetrahedra.mIndices[0]].mMaterialCoordinates;
		auto& m1 = group->mVertices[tetrahedra.mIndices[1]].mMaterialCoordinates;
		auto& m2 = group->mVertices[tetrahedra.mIndices[2]].mMaterialCoordinates;
		auto& m3 = group->mVertices[tetrahedra.mIndices[3]].mMaterialCoordinates;

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
			group->mVertices[tetrahedra.mIndices[j]].mMass += 0.25f * tetrahedra.mMass;
		}
	}

	for (int i = 0; i < numVertices; i++)
	{
		group->mVertices[i].mInvMass = 1.0f / group->mVertices[i].mMass;
	}
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void __declspec(dllexport) __stdcall Destroy()
{
	TestDeformation::TetraGroup* group = (TestDeformation::TetraGroup*)mData;
	delete group;
}

////////////////////////////////////////////////////////////////////////////////
//
extern "C" void __declspec(dllexport) __stdcall Deform(
	void* vertexPositions,
	void* vertexVelocities,
	bool useVelocities,
	void* vertexForces,
	bool useForces,
	int numVertices,
	float timestep,
	int numSteps)
{
	auto positions = reinterpret_cast<float*>(vertexPositions);
	auto velocities = reinterpret_cast<float*>(vertexVelocities);
	auto forces = reinterpret_cast<float*>(vertexForces);

	TestDeformation::TetraGroup* group = (TestDeformation::TetraGroup*)mData;

	for (int i = 0; i < numSteps; i++)
	{
		try {
			group->Update(timestep);
		}
		catch (const std::exception& e)
		{
			std::cout << "Update failed due to " << e.what() << std::endl;
			break;
		}
	}

	// End simulation

	for (int i = 0; i < group->mVertices.size(); i++)
	{
		auto& vertex = group->mVertices[i];

		positions[3 * i]      = vertex.mPosition.x;
		positions[3 * i + 1]  = vertex.mPosition.y;
		positions[3 * i + 2 ] = vertex.mPosition.z;

		if (useVelocities)
		{
			velocities[3 * i] = vertex.mVelocity.x;
			velocities[3 * i + 1] = vertex.mVelocity.y;
			velocities[3 * i + 2] = vertex.mVelocity.z;
		}

		if (useForces)
		{
			forces[3 * i] = vertex.mForce.x;
			forces[3 * i + 1] = vertex.mForce.y;
			forces[3 * i + 2] = vertex.mForce.z;
		}
	}
}

//#include <iostream>
//
//int main(int argc, const char* argv[]) {
//
//	float lambda = 2.65e6f;
//	float psi = 397.f;
//	float phi = 264.f;
//	float mu = 3.97e6f;
//	float density = 1013.f;
//	float timestep = 0.001f;
//
//	std::vector<float> positions = { 0, 0, 0, 1, 0.5f, 0, 0, 1, 0, 0, 0.5f, 1 };
//	std::vector<int> indices = { 0, 1, 2, 3 };
//
//	for (int i = 0; i < 4; i++)
//	{
//		positions[3 * i + 1] += 1;
//	}
//
//	Initialize(positions.data(), 4, indices.data(), 1, lambda, psi, mu, phi, density);
//	TestDeformation::TetraGroup* group = (TestDeformation::TetraGroup*)mData;
//
//
//	Deform(positions.data(), 4, timestep, 40);
//
//	for (const auto& vertex : group->mVertices)
//	{
//		std::cout << "Position : (" << vertex.mPosition.x << ", " << vertex.mPosition.y << ", " << vertex.mPosition.z << ")" << std::endl;
//		std::cout << "Velocity : (" << vertex.mVelocity.x << ", " << vertex.mVelocity.y << ", " << vertex.mVelocity.z << ")" << std::endl;
//		std::cout << std::endl;
//	}
//
//	Destroy();
//}

////////////////////////////////////////////////////////////////////////////////