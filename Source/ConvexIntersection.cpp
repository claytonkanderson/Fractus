#include "ConvexIntersection.h"

#include <ConvexPolyhedron.h>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/hash.hpp>
#include <glm/gtx/vector_angle.hpp>

////////////////////////////////////////////////////////////////////////////////

bool Deformation::ConvexIntersection::TetrahedraIntersection(
	const std::array<glm::vec3, 4>& tetVertices, 
	const std::array<glm::vec3, 4>& otherVertices,
	CollisionResults& results)
{
	std::vector<gte::Vector3<float>> convertedVerts
	{
		gte::Vector3<float>{tetVertices[0].x,tetVertices[0].y, tetVertices[0].z},
		gte::Vector3<float>{tetVertices[1].x,tetVertices[1].y, tetVertices[1].z},
		gte::Vector3<float>{tetVertices[2].x,tetVertices[2].y, tetVertices[2].z},
		gte::Vector3<float>{tetVertices[3].x,tetVertices[3].y, tetVertices[3].z}
	};

	std::vector<gte::Vector3<float>> convertedOtherVerts
	{
		gte::Vector3<float>{otherVertices[0].x, otherVertices[0].y, otherVertices[0].z},
		gte::Vector3<float>{otherVertices[1].x, otherVertices[1].y, otherVertices[1].z},
		gte::Vector3<float>{otherVertices[2].x, otherVertices[2].y, otherVertices[2].z},
		gte::Vector3<float>{otherVertices[3].x, otherVertices[3].y, otherVertices[3].z}
	};

	// Must match tetrahedra indices above
	static const std::vector<int> indices =
	{
		0, 2, 1,
		0, 1, 3,
		0, 3, 2,
		1, 2, 3
	};

	ConvexPolyhedron<float> poly(convertedVerts, indices);
	ConvexPolyhedron<float> other(convertedOtherVerts, indices);
	ConvexPolyhedron<float> intersection;

	if (!poly.FindIntersection(other, intersection))
		return false;

	results.mVolume = intersection.GetVolume();
	results.mForce = glm::vec3(0.0f);

	gte::Vector3<float> forceDirection{ 0,0,0 };
	// Tolerance for triangles to be considered the same
	const float cAngleTolerance = 0.01f;

	for (int i = 0; i < intersection.GetNumTriangles(); i++)
	{
		const auto& triangleNormal = intersection.GetPlane(i).normal;

		for (int j = 0; j < 4; j++)
		{
			const auto& v0 = convertedVerts[indices[3 * j + 0]];
			const auto& v1 = convertedVerts[indices[3 * j + 1]];
			const auto& v2 = convertedVerts[indices[3 * j + 2]];
			const auto& d0 = glm::normalize(glm::vec3(v1[0], v1[1], v1[2]) - glm::vec3(v0[0], v0[1], v0[2]));
			const auto& d1 = glm::normalize(glm::vec3(v2[0], v2[1], v2[2]) - glm::vec3(v0[0], v0[1], v0[2]));
			auto n = glm::cross(d0, d1);
			auto angle = glm::angle(n, glm::vec3(triangleNormal[0], triangleNormal[1], triangleNormal[2]));

			// This means that we think this facet belongs to the first tet
			if (std::abs(angle) < cAngleTolerance)
			{
				// calculate area of intersection triangle
				const auto& triangle = intersection.GetTriangle(i);
				const auto& p0 = intersection.GetPoint(triangle.GetVertex(0));
				const auto& p1 = intersection.GetPoint(triangle.GetVertex(1));
				const auto& p2 = intersection.GetPoint(triangle.GetVertex(2));
				auto area = 0.5f * gte::Length(gte::Cross(p1 - p0, p2 - p0));

				results.mForce += area * n;
			}
		}
	}

	results.mForce = glm::normalize(results.mForce);
}

////////////////////////////////////////////////////////////////////////////////
