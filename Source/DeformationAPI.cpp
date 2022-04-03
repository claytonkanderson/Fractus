#include "Deformation.hpp"
#include "DeformationAPI.h"
#include "MeshVolume.h"
#include "FractureUtil.h"
#include "ProtoConverter.hpp"

#include <glm/vec4.hpp>
#include <glm/common.hpp>
#include <iostream>
#include <fstream>
#include <Mathematics/Delaunay3.h>
#include <Mathematics/IntrSegment3Plane3.h>
#include <Mathematics/AlignedBox.h>
#include <Mathematics/DistPointHyperplane.h>
#include <tetgen/tetgen.h>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/hash.hpp>
#include <glm/gtx/vector_angle.hpp>

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

	Deformation::TetraGroup* group = new Deformation::TetraGroup();
	mData = group;

	if (density <= 1e-6)
		return -5;

	group->Initialize(lambda, psi, phi, mu, density, toughness);
	group->mMaxNumTetrahedra = maxNumTetrahedra;
	group->mMaxNumVertices = maxNumVertices;
	// Could reserve the max memory here
	group->mVertices.resize(numVertices);

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
		Deformation::Tetrahedra tetrahedra;

		tetrahedra.mIndices[0] = indices[4 * i];
		tetrahedra.mIndices[1] = indices[4 * i + 1];
		tetrahedra.mIndices[2] = indices[4 * i + 2];
		tetrahedra.mIndices[3] = indices[4 * i + 3];

		group->mIdToTetrahedra[group->mTetIdCounter++] = tetrahedra;
	}

	if (group->mDensity <= 1e-6)
		return -1;

	group->ComputeDerivedQuantities();

	for (auto& pair : group->mIdToTetrahedra)
	{
		if (pair.second.mVolume <= 1e-6)
			return 1;
		if (pair.second.mMass <= 1e-6)
			return -2;

		// TODO - not a good place to put this!
		pair.second.mRestVolume = pair.second.mVolume;
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
	Deformation::TetraGroup* group = (Deformation::TetraGroup*)mData;
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

	Deformation::TetraGroup* group = (Deformation::TetraGroup*)mData;

	//for (const auto& tet : group->mTetrahedra)
	//{
	//	if (tet.mVolume <= 1e-6f)
	//		return 1;
	//}

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

	//for (const auto& tet : group->mTetrahedra)
	//{
	//	if (tet.mVolume <= 1e-6f)
	//		return 2;
	//}

	//for (const auto& vertex : group->mVertices)
	//{
	//	if (vertex.mMass <= 1e-6f)
	//		return 3;
	//}

	// End simulation

	*numVertices = group->mVertices.size();
	*numTetrahedra = group->mIdToTetrahedra.size();

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
	size_t counter = 0;
	for (const auto & pair : group->mIdToTetrahedra)
	{
		const auto& tetrahedra = pair.second;
		indices[4 * counter + 0] = tetrahedra.mIndices[0];
		indices[4 * counter + 1] = tetrahedra.mIndices[1];
		indices[4 * counter + 2] = tetrahedra.mIndices[2];
		indices[4 * counter + 3] = tetrahedra.mIndices[3];
		
		counter++;
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

extern "C" int __declspec(dllexport) __stdcall CreateSurfaceMesh(
	float* cubeCenters, int numCubes, float* tetrahedraPositions, int numTetrahedra,
	float* outVertexPositions, int* outNumVertices, int numMaxVertices,
	int* outTriangleIndices, int* outNumTriangles, int numMaxTriangles)
{
	std::vector<glm::vec3> vertices(8*numCubes + 4*numTetrahedra);
	std::vector<glm::ivec3> notOrientedFaceIds(12 * numCubes + 4 * numTetrahedra);
	std::unordered_map<glm::ivec3, glm::ivec3> faceIdToOrientedFace;

	int vertexCounter = 0;
	int triCounter = 0;
	for (int i = 0; i < numCubes; i++)
	{
		for (int j = 0; j < 12; j++)
		{
			auto orientedFace = glm::ivec3(
				vertexCounter + Deformation::CubeTriangleIndices[j][0],
				vertexCounter + Deformation::CubeTriangleIndices[j][1],
				vertexCounter + Deformation::CubeTriangleIndices[j][2]);
			auto id = Deformation::GetFaceId(orientedFace[0], orientedFace[1], orientedFace[2]);
			notOrientedFaceIds[triCounter++] = id;
			faceIdToOrientedFace[id] = orientedFace;
		}

		glm::vec3 cubeCenter(cubeCenters[3 * i], cubeCenters[3 * i + 1], cubeCenters[3 * i + 2]);
		for (int j = 0; j < 8; j++)
		{
			vertices[vertexCounter++] = cubeCenter + Deformation::CubeVertexPositions[j] - 0.5f*glm::vec3(1.0f);
		}
	}

	for (int i = 0; i < numTetrahedra; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			auto orientedFace = vertexCounter + Deformation::TetrahedraIndices[j];
			auto id = Deformation::GetFaceId(orientedFace[0], orientedFace[1], orientedFace[2]);
			notOrientedFaceIds[triCounter++] = id;
			faceIdToOrientedFace[id] = orientedFace;
		}

		vertices[vertexCounter++] = glm::vec3(tetrahedraPositions[i * 12 + 0], tetrahedraPositions[i * 12 + 1], tetrahedraPositions[i * 12 + 2]);
		vertices[vertexCounter++] = glm::vec3(tetrahedraPositions[i * 12 + 3], tetrahedraPositions[i * 12 + 4], tetrahedraPositions[i * 12 + 5]);
		vertices[vertexCounter++] = glm::vec3(tetrahedraPositions[i * 12 + 6], tetrahedraPositions[i * 12 + 7], tetrahedraPositions[i * 12 + 8]);
		vertices[vertexCounter++] = glm::vec3(tetrahedraPositions[i * 12 + 9], tetrahedraPositions[i * 12 + 10], tetrahedraPositions[i * 12 + 11]);
	}

	std::unordered_map<int,int> oldVertIdToNewVertId;
	std::vector<glm::vec3> uniqueVertices;

	const float spatialThreshold = 0.05f;

	for (int i = 0; i < vertices.size(); i++)
	{
		bool foundNearest = false;
		int nearestId = -1;

		for (int j = 0; j < uniqueVertices.size(); j++)
		{
			if (glm::distance(vertices[i], uniqueVertices[j]) < spatialThreshold)
			{
				foundNearest = true;
				nearestId = j;
				break;
			}
		}

		if (foundNearest)
		{
			oldVertIdToNewVertId[i] = nearestId;
			continue;
		}

		uniqueVertices.push_back(vertices[i]);
		oldVertIdToNewVertId[i] = uniqueVertices.size() - 1;
	}

	if (uniqueVertices.size() > numMaxVertices)
		return -1;

	std::unordered_map<glm::ivec3, glm::ivec3> reindexedFaceIdToOrientedFaceId;

	for (int i = 0; i < notOrientedFaceIds.size(); i++)
	{
		const auto& faceId = notOrientedFaceIds[i];
		const auto& orientedFaceId = faceIdToOrientedFace[faceId];

		// Can fail if any input triangle has vertices closer together than the spatial tolerance
		glm::ivec3 newFaceId;
		try
		{
			newFaceId = Deformation::GetFaceId(oldVertIdToNewVertId[faceId[0]], oldVertIdToNewVertId[faceId[1]], oldVertIdToNewVertId[faceId[2]]);
		}
		catch (const std::exception& e)
		{
			return -10;
		}

		notOrientedFaceIds[i] = newFaceId;
		reindexedFaceIdToOrientedFaceId[newFaceId] = glm::ivec3(oldVertIdToNewVertId[orientedFaceId[0]], oldVertIdToNewVertId[orientedFaceId[1]], oldVertIdToNewVertId[orientedFaceId[2]]);
	}

	std::unordered_set<glm::ivec3> notOrientedUniqueFaceIds;

	for (const auto& faceId : notOrientedFaceIds)
		notOrientedUniqueFaceIds.insert(faceId);

	// okay now we have all unique vertices and the correct triangle indices
	Deformation::MeshVolume meshVolume;
	meshVolume.mCubeCenters.resize(numCubes);
	for (int i = 0; i < numCubes; i++)
		meshVolume.mCubeCenters[i] = glm::vec3(cubeCenters[3 * i + 0], cubeCenters[3 * i + 1], cubeCenters[3 * i + 2]);

	meshVolume.mTetrahedra.resize(numTetrahedra);
	for (int i = 0; i < numTetrahedra; i++)
		for (int j = 0; j < 4; j++)
			meshVolume.mTetrahedra[i].mVertices[j] = glm::vec3(tetrahedraPositions[12 * i + 3*j + 0], tetrahedraPositions[12 * i + 3 * j+ 1], tetrahedraPositions[12 * i + 3 * j + 2]);

	// next step is to evaluate the center of each triangle for containment

	for (auto iter = notOrientedUniqueFaceIds.begin(); iter != notOrientedUniqueFaceIds.end(); )
	{
		// triangle center
		auto triangle = *iter;
		glm::vec3 center = (uniqueVertices[triangle[0]] + uniqueVertices[triangle[1]] + uniqueVertices[triangle[2]]) / 3.0f;
		if (meshVolume.Contains(center))
			iter = notOrientedUniqueFaceIds.erase(iter);
		else
			++iter;
	}

	// finally write out data

	if (notOrientedUniqueFaceIds.size() > numMaxTriangles)
		return -2;

	outNumVertices[0] = uniqueVertices.size();
	outNumTriangles[0] = notOrientedUniqueFaceIds.size();

	// could do this with direct assignment instead of loop?
	for (int i = 0; i < uniqueVertices.size(); i++)
	{
		outVertexPositions[3 * i + 0] = uniqueVertices[i].x;
		outVertexPositions[3 * i + 1] = uniqueVertices[i].y;
		outVertexPositions[3 * i + 2] = uniqueVertices[i].z;
	}

	int outTriCounter = 0;
	for (auto iter = notOrientedUniqueFaceIds.begin(); iter != notOrientedUniqueFaceIds.end(); ++iter)
	{
		const auto& faceId = (*iter);
		const auto& orientedFace = reindexedFaceIdToOrientedFaceId[faceId];
		outTriangleIndices[3 * outTriCounter + 0] = orientedFace[0];
		outTriangleIndices[3 * outTriCounter + 1] = orientedFace[1];
		outTriangleIndices[3 * outTriCounter + 2] = orientedFace[2];
		outTriCounter++;
	}

	return 0;
}

////////////////////////////////////////////////////////////////////////////////

extern "C" int __declspec(dllexport) __stdcall CreateTetrahedralMesh(
	float* vertexPositions, int numVertices, int* triangleIndices, int numTriangles,
	float* outVertexPositions, int* outNumVertices, int numMaxVertices,
	int* outTetrahedralIndices, int* outNumTetrahedra, int numMaxTetrahedra)
{
	// lets try assuming each triangle is a facet with one polygon

	tetgenio in, out;
	tetgenio::facet* f;
	tetgenio::polygon* p;
	in.firstnumber = 0;

	in.numberofpoints = numVertices;
	in.pointlist = new REAL[in.numberofpoints * 3];
	for (int i = 0; i < numVertices; i++)
	{
		in.pointlist[3 * i + 0] = vertexPositions[3 * i + 0];
		in.pointlist[3 * i + 1] = vertexPositions[3 * i + 1];
		in.pointlist[3 * i + 2] = vertexPositions[3 * i + 2];
	}

	in.numberoffacets = numTriangles;
	in.facetlist = new tetgenio::facet[in.numberoffacets];
	//in.facetmarkerlist = new int[in.numberoffacets];

	for (int i = 0; i < numTriangles; i++)
	{
		f = &in.facetlist[i];
		f->numberofpolygons = 1;
		f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
		f->numberofholes = 0;
		f->holelist = NULL;
		p = &f->polygonlist[0];
		p->numberofvertices = 3;
		p->vertexlist = new int[p->numberofvertices];
		p->vertexlist[0] = triangleIndices[3 * i + 0];
		p->vertexlist[1] = triangleIndices[3 * i + 2];
		p->vertexlist[2] = triangleIndices[3 * i + 1];

		//in.facetmarkerlist[i] = 0;
	}

	// Tetrahedralize the PLC. Switches are chosen to read a PLC (p),
	//   do quality mesh generation (q) with a specified quality bound
	//   (1.414), and apply a maximum volume constraint (a0.1).
	//const char* switches = "pq1.414a2.0";
	const char* switches = "p";
	char* nonConstSwitches = const_cast<char*>(switches);

	tetrahedralize(nonConstSwitches, &in, &out);

	if (out.numberofpoints > numMaxVertices)
		return -1;

	if (out.numberoftetrahedra > numMaxTetrahedra)
		return -2;

	outNumVertices[0] = out.numberofpoints;
	outNumTetrahedra[0] = out.numberoftetrahedra;

	auto tetCornerPtr = out.tetrahedronlist;
	auto vertexPtr = out.pointlist;

	for (int i = 0; i < outNumVertices[0]; i++)
	{
		outVertexPositions[3 * i + 0] = vertexPtr[3 * i + 0];
		outVertexPositions[3 * i + 1] = vertexPtr[3 * i + 1];
		outVertexPositions[3 * i + 2] = vertexPtr[3 * i + 2];
	}

	for (int i = 0; i < out.numberoftetrahedra; i++)
	{
		outTetrahedralIndices[4 * i + 0] = tetCornerPtr[4 * i + 0];
		outTetrahedralIndices[4 * i + 1] = tetCornerPtr[4 * i + 1];
		outTetrahedralIndices[4 * i + 2] = tetCornerPtr[4 * i + 2];
		outTetrahedralIndices[4 * i + 3] = tetCornerPtr[4 * i + 3];
	}

	return 0;
}

////////////////////////////////////////////////////////////////////////////////