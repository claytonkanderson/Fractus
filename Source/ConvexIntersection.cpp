#include "ConvexIntersection.h"
#include "FractureUtil.h"

#include <iostream>
#include <ConvexPolyhedron.h>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/hash.hpp>
#include <glm/gtx/vector_angle.hpp>

////////////////////////////////////////////////////////////////////////////////

// apply forces directly to vertices
// take tetrahedra and vertices as input
// no output

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

	// TODO - cache in static
	const std::vector<int> indices =
	{
		TetrahedraIndices[0][0],TetrahedraIndices[0][1],TetrahedraIndices[0][2],
		TetrahedraIndices[1][0],TetrahedraIndices[1][1],TetrahedraIndices[1][2],
		TetrahedraIndices[2][0],TetrahedraIndices[2][1],TetrahedraIndices[2][2],
		TetrahedraIndices[3][0],TetrahedraIndices[3][1],TetrahedraIndices[3][2],
	};

	ConvexPolyhedron<float> poly(convertedVerts, indices);
	ConvexPolyhedron<float> other(convertedOtherVerts, indices);

	for (const auto& gon : { poly, other })
	{
		if (!gon.ValidateHalfSpaceProperty(-1e-6))
		{
			for (const auto& pt : gon.GetPoints())
				std::cout << "pos : " << pt[0] << ", " << pt[1] << ", " << pt[2] << std::endl;

			std::cout << "vol " << gon.GetVolume() << std::endl;
			throw std::exception("Invalid polyhedron.");
		}
	}

	ConvexPolyhedron<float> intersection;

	if (!poly.FindIntersection(other, intersection))
		return false;

	results.mVolume = intersection.GetVolume();
	results.mForce = glm::vec3(0.0f);

	gte::Vector3<float> forceDirection{ 0,0,0 };
	// Tolerance for triangles to be considered the same
	const float cAngleTolerance = 0.05f;
	float bestAlignedPlane = 10.0f;
	bool matchedAtLeastOnePlane = false;

	for (int i = 0; i < intersection.GetNumTriangles(); i++)
	{
		const auto& triangleNormal = intersection.GetPlane(i).normal;
		
		// compare the four tet triangles to the planes of the intersection polygon
		for (int j = 0; j < 4; j++)
		{
			const auto& v0 = convertedVerts[indices[3 * j + 0]];
			const auto& v1 = convertedVerts[indices[3 * j + 1]];
			const auto& v2 = convertedVerts[indices[3 * j + 2]];
			const auto& d0 = glm::vec3(v1[0], v1[1], v1[2]) - glm::vec3(v0[0], v0[1], v0[2]);
			const auto& d1 = glm::vec3(v2[0], v2[1], v2[2]) - glm::vec3(v0[0], v0[1], v0[2]);
			auto n = glm::normalize(glm::cross(d0, d1));
			auto angle = glm::angle(n, glm::vec3(triangleNormal[0], triangleNormal[1], triangleNormal[2]));

			bool acceptPlane = (angle < cAngleTolerance || angle > glm::pi<float>() - cAngleTolerance);

			if (std::abs(angle) < bestAlignedPlane)
				bestAlignedPlane = std::abs(angle);

			// This means that we think this facet belongs to the first tet
			if (acceptPlane)
			{
				// calculate area of intersection triangle
				const auto& triangle = intersection.GetTriangle(i);
				const auto& p0 = intersection.GetPoint(triangle.GetVertex(0));
				const auto& p1 = intersection.GetPoint(triangle.GetVertex(1));
				const auto& p2 = intersection.GetPoint(triangle.GetVertex(2));
				auto area = 0.5f * gte::Length(gte::Cross(p1 - p0, p2 - p0));

				results.mForce += area * n;
				matchedAtLeastOnePlane = true;
			}
		}
	}

	if (results.mVolume < 1e-8f)
		return false;

	if (matchedAtLeastOnePlane && glm::length(results.mForce) < 1e-6f)
		return false;

	if (!matchedAtLeastOnePlane)
	{
		std::cout << "First Tet Normals" << std::endl;
		{
			for (int j = 0; j < 4; j++)
			{
				const auto& v0 = convertedVerts[indices[3 * j + 0]];
				const auto& v1 = convertedVerts[indices[3 * j + 1]];
				const auto& v2 = convertedVerts[indices[3 * j + 2]];
				const auto& d0 = glm::vec3(v1[0], v1[1], v1[2]) - glm::vec3(v0[0], v0[1], v0[2]);
				const auto& d1 = glm::vec3(v2[0], v2[1], v2[2]) - glm::vec3(v0[0], v0[1], v0[2]);
				auto n = glm::normalize(glm::cross(d0, d1));

				std::cout << n.x << ", " << n.y << ", " << n.z << std::endl;
			}
		}

		std::cout << "Second Tet Normals" << std::endl;
		{
			for (int j = 0; j < 4; j++)
			{
				const auto& v0 = convertedOtherVerts[indices[3 * j + 0]];
				const auto& v1 = convertedOtherVerts[indices[3 * j + 1]];
				const auto& v2 = convertedOtherVerts[indices[3 * j + 2]];
				const auto& d0 = glm::vec3(v1[0], v1[1], v1[2]) - glm::vec3(v0[0], v0[1], v0[2]);
				const auto& d1 = glm::vec3(v2[0], v2[1], v2[2]) - glm::vec3(v0[0], v0[1], v0[2]);
				auto n = glm::normalize(glm::cross(d0, d1));

				std::cout << n.x << ", " << n.y << ", " << n.z << std::endl;
			}
		}

		std::cout << "Intersection Normals" << std::endl;
		{
			for (int i = 0; i < intersection.GetNumTriangles(); i++)
			{
				const auto& n = intersection.GetPlane(i).normal;
				std::cout << n[0] << ", " << n[1] << ", " << n[2] << std::endl;
			}
		}

		std::cout << glm::angle(glm::vec3(0, 0, 1), glm::vec3(0, 0, -1)) << std::endl;

		std::cout << "Best Angle : " << bestAlignedPlane << std::endl;
		std::cout << "Collision volume : " << results.mVolume << std::endl;
		throw std::exception("No faces aligned during collision.");
	}

	results.mForce = glm::normalize(results.mForce);

	return true;
}

////////////////////////////////////////////////////////////////////////////////

void Deformation::ConvexIntersection::ResolveCollisions(std::vector<Vertex>& vertices, std::unordered_map<size_t, Tetrahedra>& tetrahedra)
{
	const float cCollisionPerUnitVolume = 1e8f;
	// Tolerance for triangles to be considered the same
	const float cAngleTolerance = 0.05f;

	const std::vector<int> indices =
	{
		TetrahedraIndices[0][0],TetrahedraIndices[0][1],TetrahedraIndices[0][2],
		TetrahedraIndices[1][0],TetrahedraIndices[1][1],TetrahedraIndices[1][2],
		TetrahedraIndices[2][0],TetrahedraIndices[2][1],TetrahedraIndices[2][2],
		TetrahedraIndices[3][0],TetrahedraIndices[3][1],TetrahedraIndices[3][2],
	};

	std::vector<gte::Vector3<float>> convertedVertices(4);
	std::vector<gte::Vector3<float>> convertedOtherVerts(4);

	std::unordered_map<size_t, ConvexPolyhedron<float>> tetIdToConvexRep;
	ConvexPolyhedron<float> intersection;

	for (auto& pair : tetrahedra)
	{
		auto& tet = pair.second;

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				convertedVertices[i][j] = vertices[tet.mIndices[i]].mPosition[j];
			}
		}

		tetIdToConvexRep.insert({ pair.first, ConvexPolyhedron<float>(convertedVertices, indices) });
	}

	ConvexPolyhedron<float>::Clipper clipper;
	auto findIntersection = [&](const ConvexPolyhedron<float> & tet, const ConvexPolyhedron<float>& otherTet, ConvexPolyhedron<float> & intersection)
	{
		clipper.Clear();
		clipper.Initialize(tet, 1e-6f);

		for (auto const& plane : otherTet.GetPlanes())
		{
			if (clipper.Clip(plane) < 0)
			{
				return false;
			}
		}

		clipper.Convert(intersection);
		return true;
	};
	
	std::vector<size_t> tetrahedraIds;
	tetrahedraIds.reserve(tetrahedra.size());
	for (const auto& pair : tetrahedra)
		tetrahedraIds.push_back(pair.first);

	size_t tetCounter = 0;
	for (size_t tetIdx0 = 0; tetIdx0 < tetrahedraIds.size(); tetIdx0++)
	{
		std::cout << tetCounter++ << std::endl;

		for (size_t tetIdx1 = tetIdx0 + 1; tetIdx1 < tetrahedraIds.size(); tetIdx1++)
		{
			auto& firstTet = tetrahedra[tetrahedraIds[tetIdx0]];
			auto& secondTet = tetrahedra[tetrahedraIds[tetIdx1]];

			// TODO remove this?
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					convertedVertices[i][j] = vertices[firstTet.mIndices[i]].mPosition[j];
					convertedOtherVerts[i][j] = vertices[secondTet.mIndices[i]].mPosition[j];
				}
			}

			const auto& poly = tetIdToConvexRep[tetrahedraIds[tetIdx0]];
			const auto& other = tetIdToConvexRep[tetrahedraIds[tetIdx1]];

			intersection.Reset();

			if (!findIntersection(poly, other, intersection))
				continue;

			auto intersectionVolume = intersection.GetVolume();

			if (intersectionVolume < 1e-8f)
				continue;

			auto intersectionSurfaceArea = intersection.GetSurfaceArea();

			for (auto& intersectionPlane : intersection.GetPlanes())
			{
				bool matchedAtLeastOnePlane = false;

				float bestAlignedPlane = 10.0f;
				const auto & planeNormal = intersectionPlane.normal;

				// Loop over the individual triangles of each tetrahedra that are colliding
				// to see if they're part of the intersection, and if they are,
				// apply force to the vertices in that triangle proportional to 
				// volume of intersection and surface area of triangle.

				for (auto& tet : { firstTet, secondTet })
				{
					for (int i = 0; i < 4; i++)
					{
						// These are the three vertices of the triangle
						auto & v0 = vertices[tet.mIndices[indices[3 * i + 0]]];
						auto & v1 = vertices[tet.mIndices[indices[3 * i + 1]]];
						auto & v2 = vertices[tet.mIndices[indices[3 * i + 2]]];
						const auto& p0 = v0.mPosition;
						const auto& p1 = v1.mPosition;
						const auto& p2 = v2.mPosition;

						const auto& d0 = glm::vec3(p1[0], p1[1], p1[2]) - glm::vec3(p0[0], p0[1], p0[2]);
						const auto& d1 = glm::vec3(p2[0], p2[1], p2[2]) - glm::vec3(p0[0], p0[1], p0[2]);
						auto n = glm::normalize(glm::cross(d0, d1));
						auto angle = glm::angle(n, glm::vec3(planeNormal[0], planeNormal[1], planeNormal[2]));

						bool acceptPlane = (angle < cAngleTolerance || angle > glm::pi<float>() - cAngleTolerance);

						if (!acceptPlane)
							continue;

						// calculate area of intersection triangle
						const auto& triangle = intersection.GetTriangle(i);
						const auto& triangleP0 = intersection.GetPoint(triangle.GetVertex(0));
						const auto& triangleP1 = intersection.GetPoint(triangle.GetVertex(1));
						const auto& triangleP2 = intersection.GetPoint(triangle.GetVertex(2));
						auto intersectionTriangleArea = 0.5f * gte::Length(gte::Cross(triangleP1 - triangleP0, triangleP2 - triangleP0));

						auto weightedNormal = intersectionTriangleArea * n;
						auto directedForce = weightedNormal / intersectionSurfaceArea * intersectionVolume * cCollisionPerUnitVolume;
						// apply this weightedNormal * collision strength to the nodes in this triangle
						v0.mForce += directedForce;
						v1.mForce += directedForce;
						v2.mForce += directedForce;

						matchedAtLeastOnePlane = true;
					}
				}

				if (!matchedAtLeastOnePlane)
				{
					std::cout << "Best Angle : " << bestAlignedPlane << std::endl;
					std::cout << "Collision volume : " << intersectionVolume << std::endl;
					throw std::exception("No faces aligned during collision.");
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
