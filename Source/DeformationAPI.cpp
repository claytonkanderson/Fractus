#include "TestDeformation.hpp"
#include <glm/vec4.hpp>
#include <glm/common.hpp>
#include <iostream>

#include <Mathematics/Delaunay3.h>
#include <Mathematics/IntrSegment3Plane3.h>
#include <Mathematics/AlignedBox.h>

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

////////////////////////////////////////////////////////////////////////////////
//
// cubeMinMax should be {minX, minY, minZ, maxX, maxY, maxZ}
// planeNormalOrigin should be {pX, pY, pZ, nX, nY, nZ}
// outVertices should be preallocated to the maximum number of vertices, which I think is 10
// outTetrahedraIndices should be preallocated to 4 * maximum number of tetrahedra, which I think is 3?
//
extern "C" void __declspec(dllexport) __stdcall TetrahedralizeCubeIntersection(
	float* cubeMinMax, // x, y, z, x, y, z
	float* planeNormalOrigin,
	float* outVertices,
	int* outNumVertices,
	int* outTetrahedraIndices,
	int* outNumTetrahedra,
	bool positivePlane
)
{
	glm::vec3 min(cubeMinMax[0], cubeMinMax[1], cubeMinMax[2]);
	glm::vec3 max(cubeMinMax[3], cubeMinMax[4], cubeMinMax[5]);
	gte::AlignedBox3<float> box3d({{min.x, min.y, min.z}, {max.x, max.y, max.z}});

	gte::Vector3<float> planeNormal
	{
		planeNormalOrigin[0],
		planeNormalOrigin[1],
		planeNormalOrigin[2]
	};
	gte::Vector3<float> planeOrigin
	{
		planeNormalOrigin[3],
		planeNormalOrigin[4],
		planeNormalOrigin[5]
	};
	gte::Plane3<float> plane(planeNormal, planeOrigin);

	std::vector<gte::Vector3<float>> intersections;
	float dX = max.x - min.x;
	float dY = max.y - min.y;
	float dZ = max.z - min.z;

	std::vector<gte::Segment3<float>> edges(12);
	// X edges
	edges[0] = { { min.x, min.y, min.z }, { min.x + dX, min.y, min.z } };
	edges[1] = { { min.x, min.y + dY, min.z }, { min.x + dX, min.y + dY, min.z } };
	edges[2] = { { min.x, min.y, min.z + dZ}, { min.x + dX, min.y, min.z + dZ} };
	edges[3] = { { min.x, min.y + dY, min.z + dZ}, { min.x + dX, min.y + dY, min.z + dZ} };
	// Y edges
	edges[4] = { { min.x, min.y, min.z }, { min.x, min.y + dY, min.z } };
	edges[5] = { { min.x + dX, min.y, min.z }, { min.x + dX, min.y + dY, min.z } };
	edges[6] = { { min.x, min.y, min.z + dZ}, { min.x, min.y + dY, min.z + dZ} };
	edges[7] = { { min.x + dX, min.y, min.z + dZ}, { min.x + dX, min.y + dY, min.z + dZ} };
	// Z edges
	edges[8 ] = { { min.x, min.y, min.z }, { min.x, min.y, min.z + dZ } };
	edges[9 ] = { { min.x + dX, min.y, min.z }, { min.x + dX, min.y, min.z + dZ} };
	edges[10] = { { min.x, min.y + dY, min.z }, { min.x, min.y + dY, min.z + dZ} };
	edges[11] = { { min.x + dX, min.y + dY, min.z }, { min.x + dX, min.y + dY, min.z + dZ } };

	gte::FIQuery<float, gte::Segment3<float>, gte::Plane3<float>> query;
	for (const auto& edge : edges)
	{
		auto result = query(edge, plane);
		if (result.intersect)
			intersections.push_back(result.point);
	}

	std::vector<gte::Vector3<float>>  positiveVertices(intersections);
	std::vector<gte::Vector3<float>>  negativeVertices(intersections);
	std::array<gte::Vector3<float>, 8> cubeVertices;
	box3d.GetVertices(cubeVertices);

	for (const auto& vertex : cubeVertices)
	{
		// vertex - origin is a vector from the plane to the vertex
		// if we dot with the normal we get a scalar distance
		auto displacement = vertex - planeOrigin;
		auto signedDistance = displacement[0] * planeNormal[0] + displacement[1] * planeNormal[1] + displacement[2] * planeNormal[2];
		if (signedDistance > 0)
			positiveVertices.push_back(vertex);
		else
			negativeVertices.push_back(vertex);
	}

	gte::Delaunay3<float> triangulator;

	if (positivePlane)
		triangulator(positiveVertices);
	else
		triangulator(negativeVertices);

	*outNumVertices = triangulator.GetNumVertices();
	*outNumTetrahedra = triangulator.GetNumTetrahedra();

	for (int i = 0; i < triangulator.GetNumVertices(); i++)
	{
		outVertices[3 * i] = triangulator.GetVertices()[i][0];
		outVertices[3 * i + 1] = triangulator.GetVertices()[i][1];
		outVertices[3 * i + 2] = triangulator.GetVertices()[i][2];
	}

	for (int i = 0; i < triangulator.GetIndices().size(); i++)
	{
		outTetrahedraIndices[i] = triangulator.GetIndices()[i];
	}
}

////////////////////////////////////////////////////////////////////////////////

#include <iostream>

//int main(int argc, const char* argv[]) {
//
//	{
//		//float lambda = 2.65e6f;
//		//float psi = 397.f;
//		//float phi = 264.f;
//		//float mu = 3.97e6f;
//		//float density = 1013.f;
//		//float timestep = 0.001f;
//
//		//std::vector<float> positions = { 0, 0, 0, 1, 0.5f, 0, 0, 1, 0, 0, 0.5f, 1 };
//		//std::vector<int> indices = { 0, 1, 2, 3 };
//
//		//for (int i = 0; i < 4; i++)
//		//{
//		//	positions[3 * i + 1] += 1;
//		//}
//
//		//Initialize(positions.data(), 4, indices.data(), 1, lambda, psi, mu, phi, density);
//		//TestDeformation::TetraGroup* group = (TestDeformation::TetraGroup*)mData;
//
//		//Deform(positions.data(), );
//
//		//for (const auto& vertex : group->mVertices)
//		//{
//		//	std::cout << "Position : (" << vertex.mPosition.x << ", " << vertex.mPosition.y << ", " << vertex.mPosition.z << ")" << std::endl;
//		//	std::cout << "Velocity : (" << vertex.mVelocity.x << ", " << vertex.mVelocity.y << ", " << vertex.mVelocity.z << ")" << std::endl;
//		//	std::cout << std::endl;
//		//}
//
//		//Destroy();
//	}
//	{
//		std::vector<float> minMax = { 0, 0, 0, 2, 2, 2 };
//		std::vector<float> normalOrigin = { 0, 0, 1, 0, 0, 1 };
//		std::vector<float> vertices(10);
//		std::vector<int> indices(12);
//		int numTetrahedra;
//		int numVertices;
//		bool positivePlane = true;
//		TetrahedralizeCubeIntersection(minMax.data(), normalOrigin.data(), vertices.data(), &numVertices, indices.data(), &numTetrahedra, positivePlane);
//
//		std::cout << "Num tetrahedra : " << numTetrahedra << std::endl;
//		std::cout << "Num vertices : " << numVertices << std::endl;
//
//		for (int i = 0; i < numVertices; i++)
//		{
//			
//			std::cout << "Position : (" << vertices[3 * i] << ", " << vertices[3 * i + 1] << ", " << vertices[3 * i + 2] << ")" << std::endl;
//		}
//	}
//
//}

////////////////////////////////////////////////////////////////////////////////