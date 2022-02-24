#include "FractureUtil.h"

#include <Mathematics/DistPointHyperplane.h>
#include <Mathematics/IntrLine3Plane3.h>

#include <array>
#include <algorithm>
#include <iostream>

////////////////////////////////////////////////////////////////////////////////

glm::ivec2 Deformation::GetEdgeId(size_t idx1, size_t idx2)
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

////////////////////////////////////////////////////////////////////////////////

glm::ivec3 Deformation::GetFaceId(size_t idx1, size_t idx2, size_t idx3)
{
	if (idx1 == idx2 || idx1 == idx3 || idx2 == idx3)
		throw std::exception("Attempted to get a face with duplicate vertex ids.");

	std::array<int, 3> verts{ idx1, idx2, idx3 };
	std::sort(verts.begin(), verts.end());

	return glm::ivec3(verts[0], verts[1], verts[2]);
}

////////////////////////////////////////////////////////////////////////////////

size_t Deformation::GetCommonVertexFromEdges(const glm::ivec2& edge0, const glm::ivec2& edge1)
{
	if (edge0.x == edge1.x || edge0.x == edge1.y)
		return edge0.x;
	if (edge0.y == edge1.x || edge0.y == edge1.y)
		return edge0.y;

	throw std::exception("Failed to find common vertices for provided edges.");
	return 0;
}

////////////////////////////////////////////////////////////////////////////////

void Deformation::ComputeMMat(const glm::dvec3& eigenVector, glm::dmat3& outMat)
{
	if (glm::length(eigenVector) < 1e-6)
	{
		outMat = glm::dmat3(0.0);
		return;
	}

	for (int j = 0; j < 3; j++)
		for (int k = 0; k < 3; k++)
			outMat[j][k] = eigenVector[j] * eigenVector[k];

	outMat /= glm::length(eigenVector);
}

////////////////////////////////////////////////////////////////////////////////

double Deformation::DistPointPlane(const glm::vec3& point, const glm::vec3& planeNormal, const glm::vec3& planePosition)
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

////////////////////////////////////////////////////////////////////////////////
