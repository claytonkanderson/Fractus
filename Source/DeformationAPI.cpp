#include "TestDeformation.hpp"
#include "DeformationAPI.h"
#include "ProtoConverter.hpp"
#include <glm/vec4.hpp>
#include <glm/common.hpp>
#include <iostream>
#include <fstream>

#include <Mathematics/Delaunay3.h>
#include <Mathematics/IntrSegment3Plane3.h>
#include <Mathematics/AlignedBox.h>
#include <Mathematics/DistPointHyperplane.h>

void* mData = nullptr;

////////////////////////////////////////////////////////////////////////////////

extern "C" int __declspec(dllexport) __stdcall Initialize(
	void* vertexPositions,
	int numVertices,
	int maxNumVertices,
	void* tetrahedraIndices,
	int numTetrahedra,
	int maxNumTetrahedra,
	float lambda,
	float psi,
	float mu,
	float phi,
	float toughness,
	float density)
{
	if (mData != nullptr)
		Destroy();

	TestDeformation::TetraGroup* group = new TestDeformation::TetraGroup();
	mData = group;

	if (density <= 1e-6)
		return -5;

	group->Initialize(lambda, psi, phi, mu, density, toughness);
	group->mMaxNumTetrahedra = maxNumTetrahedra;
	group->mMaxNumVertices = maxNumVertices;
	// Could reserve the max memory here
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
	}

	if (group->mDensity <= 1e-6)
		return -1;

	group->ComputeDerivedQuantities();

	for (const auto& tet : group->mTetrahedra)
	{
		if (tet.mVolume <= 1e-6)
			return 1;
		if (tet.mMass <= 1e-6)
			return -2;
	}

	for (const auto& vert : group->mVertices)
	{
		if (vert.mMass <= 1e-6)
			return 2;
	}

	return 0;
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void __declspec(dllexport) __stdcall Destroy()
{
	TestDeformation::TetraGroup* group = (TestDeformation::TetraGroup*)mData;
	group->OutputSaveFile();
	delete group;
	mData = nullptr;
}

////////////////////////////////////////////////////////////////////////////////

extern "C" int __declspec(dllexport) __stdcall Deform(
	void* vertexPositions,
	void* vertexVelocities,
	bool useVelocities,
	void* vertexForces,
	bool useForces,
	void* tetrahedraIndices,
	int* numVertices,
	int* numTetrahedra,
	float timestep,
	int numSteps)
{
	auto positions = reinterpret_cast<float*>(vertexPositions);
	auto velocities = reinterpret_cast<float*>(vertexVelocities);
	auto forces = reinterpret_cast<float*>(vertexForces);
	auto indices = reinterpret_cast<int*>(tetrahedraIndices);

	TestDeformation::TetraGroup* group = (TestDeformation::TetraGroup*)mData;

	for (const auto& tet : group->mTetrahedra)
	{
		if (tet.mVolume <= 1e-6f)
			return 1;
	}

	for (int i = 0; i < numSteps; i++)
	{
		try {
			group->Update(timestep);
		}
		catch (const std::exception& e)
		{
			std::cout << "Update failed due to " << e.what() << std::endl;
			return 4;
		}
	}

	for (const auto& tet : group->mTetrahedra)
	{
		if (tet.mVolume <= 1e-6f)
			return 2;
	}

	for (const auto& vertex : group->mVertices)
	{
		if (vertex.mMass <= 1e-6f)
			return 3;
	}

	// End simulation

	*numVertices = group->mVertices.size();
	*numTetrahedra = group->mTetrahedra.size();

	for (int i = 0; i < group->mVertices.size(); i++)
	{
		auto& vertex = group->mVertices[i];

		positions[3 * i + 0] = (float)vertex.mPosition.x;
		positions[3 * i + 1] = (float)vertex.mPosition.y;
		positions[3 * i + 2] = (float)vertex.mPosition.z;

		if (useVelocities)
		{
			velocities[3 * i + 0] = (float)vertex.mVelocity.x;
			velocities[3 * i + 1] = (float)vertex.mVelocity.y;
			velocities[3 * i + 2] = (float)vertex.mVelocity.z;
		}

		if (useForces)
		{
			forces[3 * i + 0] = (float)vertex.mForce.x;
			forces[3 * i + 1] = (float)vertex.mForce.y;
			forces[3 * i + 2] = (float)vertex.mForce.z;
		}
	}

	// Optimize - don't need to loop over all tet's, only the new ones
	for (int i = 0; i < group->mTetrahedra.size(); i++)
	{
		indices[4 * i + 0] = group->mTetrahedra[i].mIndices[0];
		indices[4 * i + 1] = group->mTetrahedra[i].mIndices[1];
		indices[4 * i + 2] = group->mTetrahedra[i].mIndices[2];
		indices[4 * i + 3] = group->mTetrahedra[i].mIndices[3];
	}

	return 0;
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void __declspec(dllexport) __stdcall Delaunay3D(
	float* inoutVertices,
	int* inoutNumVertices,
	int maxNumVertices,
	int* outTetrahedraIndices,
	int maxNumTetrahedra,
	int* outNumTetrahedra,
	bool* success
)
{
	gte::Delaunay3<float> triangulator;

	int numInitialVertices = (*inoutNumVertices);
	std::vector<gte::Vector3<float>> vertices(numInitialVertices);

	for (int i = 0; i < numInitialVertices; i++)
		vertices[i] = gte::Vector3<float>({inoutVertices[3*i], inoutVertices[3 * i + 1] , inoutVertices[3 * i + 2] });

	*success = triangulator(vertices);

	if (!*success)
		return;

	if (maxNumVertices < triangulator.GetNumVertices() || maxNumTetrahedra < triangulator.GetNumTetrahedra())
	{
		*success = false;
		return;
	}

	*inoutNumVertices = triangulator.GetNumVertices();
	*outNumTetrahedra = triangulator.GetNumTetrahedra();

	for (int i = 0; i < triangulator.GetNumVertices(); i++)
	{
		inoutVertices[3 * i] = triangulator.GetVertices()[i][0];
		inoutVertices[3 * i + 1] = triangulator.GetVertices()[i][1];
		inoutVertices[3 * i + 2] = triangulator.GetVertices()[i][2];
	}

	for (int i = 0; i < triangulator.GetIndices().size(); i++)
	{
		outTetrahedraIndices[i] = triangulator.GetIndices()[i];
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

	std::array<gte::Vector3<float>, 8> cubeVertices;
	box3d.GetVertices(cubeVertices);

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

	// Check that the box is actually intersected by the plane, otherwise return
	// We perform this check by comparing the signed distance of each cube vertex 
	{
		gte::DCPQuery<float, gte::Vector3<float>, gte::Plane3<float>> distanceQuery;
		bool first = true;
		bool positive = false;
		bool intersected = false;
		for (int i = 0; i < 8; i++)
		{
			const auto& vertex = cubeVertices[i];
			auto result = distanceQuery(vertex, plane);

			if (result.distance <= 1e-6f && result.distance >= -1e-6f)
				continue;

			if (first)
			{
				positive = (result.signedDistance > 0.0f);
				first = false;
				continue;
			}

			if (positive && (result.signedDistance < 0.0f))
			{
				intersected = true;
				break;
			}

			if (!positive && (result.signedDistance > 0.0f))
			{
				intersected = true;
				break;
			}
		}

		if (!intersected)
		{
			*outNumVertices = 0;
			*outNumTetrahedra = 0;
			return;
		}
	}

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
		// Skip cube edges that are contained in the plane. We check for these by 
		// checking if both edge points are contained in the plane.
		gte::DCPQuery<float, gte::Vector3<float>, gte::Plane3<float>> distanceQuery;
		auto p0Query = distanceQuery(edge.p[0], plane);
		auto p1Query = distanceQuery(edge.p[1], plane);

		if (p0Query.distance <= 1e-6f && p1Query.distance <= 1e-6f)
			continue;

		auto result = query(edge, plane);
		if (result.intersect)
			intersections.push_back(result.point);
	}

	std::vector<gte::Vector3<float>>  positiveVertices(intersections);
	std::vector<gte::Vector3<float>>  negativeVertices(intersections);

	for (const auto& vertex : cubeVertices)
	{
		// vertex - origin is a vector from the plane to the vertex
		// if we dot with the normal we get a scalar distance
		auto displacement = vertex - planeOrigin;
		auto signedDistance = displacement[0] * planeNormal[0] + displacement[1] * planeNormal[1] + displacement[2] * planeNormal[2];
		if (signedDistance >= 0)
			positiveVertices.push_back(vertex);
		if (signedDistance <= 0)
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