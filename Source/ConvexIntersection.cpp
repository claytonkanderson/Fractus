#include "ConvexIntersection.h"
#include "FractureUtil.h"

#include <iostream>
#include <ConvexPolyhedron.h>
#include <Mathematics/IntrAlignedBox3AlignedBox3.h>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/hash.hpp>
#include <glm/gtx/vector_angle.hpp>

////////////////////////////////////////////////////////////////////////////////

void Deformation::ConvexIntersection::ResolveCollisions(std::vector<Vertex>& vertices, std::unordered_map<size_t, Tetrahedra>& tetrahedra)
{
	const float cCollisionPerUnitVolume = 1e8f;
	// Tolerance for triangles to be considered the same
	const float cAngleTolerance = 0.05f;
	const float cIntersectionVolumeThreshold = 1e-8f;
	const float cBoxIntersectionVolumeThreshold = 1e-4f;

	const std::vector<int> indices =
	{
		TetrahedraIndices[0][0],TetrahedraIndices[0][1],TetrahedraIndices[0][2],
		TetrahedraIndices[1][0],TetrahedraIndices[1][1],TetrahedraIndices[1][2],
		TetrahedraIndices[2][0],TetrahedraIndices[2][1],TetrahedraIndices[2][2],
		TetrahedraIndices[3][0],TetrahedraIndices[3][1],TetrahedraIndices[3][2],
	};

	std::vector<gte::Vector3<float>> convertedVertices(4);
	std::vector<gte::Vector3<float>> convertedOtherVerts(4);

	ConvexPolyhedron<float> intersection;

	ConvexPolyhedron<float>::Clipper clipper;
	auto findIntersection = [&](const ConvexPolyhedron<float> & tet, const ConvexPolyhedron<float>& otherTet, ConvexPolyhedron<float> & intersection)
	{
		intersection.Reset();
		clipper.Clear();
		clipper.Initialize(tet, 1e-6f);

		try
		{
			for (auto const& plane : otherTet.GetPlanes())
			{
				if (clipper.Clip(plane) < 0)
				{
					return false;
				}
			}
			clipper.Convert(intersection);
		}
		catch (const std::exception& e)
		{
			std::cout << e.what() << std::endl;
			std::cout << "Intersection has " << clipper.mVertices.size() << " points." << std::endl;
			std::cout << "Intersection has " << clipper.mFaces.size() << " planes." << std::endl;
			
			for (int i = 0; i < 4; i++)
				std::cout << "Tet Vert " << i << " is (" << tet.GetPoints()[i][0] << ", " << tet.GetPoints()[i][1] << ", " << tet.GetPoints()[i][2] << ")" << std::endl;
			for (int i = 0; i < 4; i++)
				std::cout << "Other Tet Vert " << i << " is (" << otherTet.GetPoints()[i][0] << ", " << otherTet.GetPoints()[i][1] << ", " << otherTet.GetPoints()[i][2] << ")" << std::endl;

			//throw std::exception(e);
			return false;
		}
		return true;
	};
	
	std::vector<ConvexPolyhedron<float>> convexPolyhedra;
	convexPolyhedra.reserve(tetrahedra.size());
	
	std::vector<gte::AlignedBox3<float>> axisAlignedBoundingBoxes;
	axisAlignedBoundingBoxes.reserve(tetrahedra.size());
	
	std::vector<size_t> tetrahedraIds;
	tetrahedraIds.reserve(tetrahedra.size());

	for (const auto& pair : tetrahedra)
	{
		tetrahedraIds.push_back(pair.first);

		auto& tet = pair.second;
		axisAlignedBoundingBoxes.emplace_back();
		auto& box = axisAlignedBoundingBoxes.back();

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				convertedVertices[i][j] = vertices[tet.mIndices[i]].mPosition[j];
			}

			if (i == 0)
			{
				box.min = convertedVertices[i];
				box.max = convertedVertices[i];
			}
			else
			{
				// Expand to include
				for (int j = 0; j < 3; j++)
				{
					box.min[j] = std::min(box.min[j], convertedVertices[i][j]);
					box.max[j] = std::max(box.max[j], convertedVertices[i][j]);
				}
			}
		}

		convexPolyhedra.push_back(ConvexPolyhedron<float>(convertedVertices, indices));
	}

	gte::FIQuery<float, gte::AlignedBox3<float>, gte::AlignedBox3<float>> boundingBoxQuery;

	size_t numValidCollisions = 0;
	size_t numCollisionsPastAABB = 0;
	for (size_t tetIdx0 = 0; tetIdx0 < tetrahedraIds.size(); tetIdx0++)
	{
		const auto& firstBox = axisAlignedBoundingBoxes[tetIdx0];

		for (size_t tetIdx1 = tetIdx0 + 1; tetIdx1 < tetrahedraIds.size(); tetIdx1++)
		{
			const auto& secondBox = axisAlignedBoundingBoxes[tetIdx1];
			auto result = boundingBoxQuery(firstBox, secondBox);
			if (!result.intersect)
				continue;

			const auto& intersectedBox = result.box;
			auto intersectedBoxVolume = (intersectedBox.max[0] - intersectedBox.min[0]) * (intersectedBox.max[1] - intersectedBox.min[1]) * (intersectedBox.max[2] - intersectedBox.min[2]);
			if (intersectedBoxVolume < cBoxIntersectionVolumeThreshold)
				continue;

			numCollisionsPastAABB++;
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

			const auto& poly = convexPolyhedra[tetIdx0];
			const auto& other = convexPolyhedra[tetIdx1];

			if (!findIntersection(poly, other, intersection))
				continue;

			auto intersectionVolume = intersection.GetVolume();

			if (intersectionVolume < cIntersectionVolumeThreshold)
				continue;

			numValidCollisions++;

			auto intersectionSurfaceArea = intersection.GetSurfaceArea();

			// we need to loop over the triangles of the intersection region
			// for each triangle we need the area and the normal
			for (int intersectionTriId = 0; intersectionTriId < intersection.GetNumTriangles(); intersectionTriId++)
			{
				auto& intersectionPlane = intersection.GetPlanes()[intersectionTriId];
				const auto& intersectionTri = intersection.GetTriangle(intersectionTriId);
				const auto& p0 = intersection.GetPoint(intersectionTri.GetVertex(0));
				const auto& p1 = intersection.GetPoint(intersectionTri.GetVertex(1));
				const auto& p2 = intersection.GetPoint(intersectionTri.GetVertex(2));

				auto intersectionTriangleArea = intersection.GetTriangleArea(intersectionPlane.normal,
					p0,
					p1,
					p2);

				//std::cout << "Triangle Area " << intersectionTriangleArea << std::endl;

				bool matchedAtLeastOnePlane = false;

				float bestAlignedPlane = 10.0f;
				const auto & planeNormal = intersectionPlane.normal;

				// Loop over the individual triangles of each tetrahedra that are colliding
				// to see if they're part of the intersection, and if they are,
				// apply force to the vertices in that triangle proportional to 
				// volume of intersection and surface area of triangle.

				for (auto& tet : { firstTet, secondTet })
				{
					const auto& tetCenter = tet.GetCentroid(vertices);

					for (int i = 0; i < 4; i++)
					{
						// These are the three vertices of the triangle
						auto & v0 = vertices[tet.mIndices[indices[3 * i + 0]]];
						auto & v1 = vertices[tet.mIndices[indices[3 * i + 1]]];
						auto & v2 = vertices[tet.mIndices[indices[3 * i + 2]]];
						const auto& p0 = v0.mPosition;
						const auto& p1 = v1.mPosition;
						const auto& p2 = v2.mPosition;

						// we need the outward normal here
						const auto& d0 = glm::vec3(p1[0], p1[1], p1[2]) - glm::vec3(p0[0], p0[1], p0[2]);
						const auto& d1 = glm::vec3(p2[0], p2[1], p2[2]) - glm::vec3(p0[0], p0[1], p0[2]);
						auto n = glm::normalize(glm::cross(d0, d1));

						// ensure outward normal
						auto cd = (glm::vec3)glm::normalize(v0.mPosition - tetCenter);
						if (glm::dot(n, cd) > 0)
							n = -n;

						auto angle = glm::angle(n, glm::vec3(planeNormal[0], planeNormal[1], planeNormal[2]));

						bool acceptPlane = (angle < cAngleTolerance || angle > glm::pi<float>() - cAngleTolerance);

						if (!acceptPlane)
							continue;


						//std::cout << "Angle " << angle << std::endl;
						//std::cout << "Triangle area " << intersectionTriangleArea << std::endl;

						auto weightedNormal = intersectionTriangleArea * n;
						auto directedForce = weightedNormal / intersectionSurfaceArea * intersectionVolume * cCollisionPerUnitVolume;

						// testing
						directedForce = n * intersectionVolume * cCollisionPerUnitVolume;

						// apply this weightedNormal * collision strength to the nodes in this triangle
						v0.mForce += directedForce;
						v1.mForce += directedForce;
						v2.mForce += directedForce;

						if (v0.mCollisionForces.empty())
							v0.mCollisionForces.push_back(directedForce);
						else
							v0.mCollisionForces.back() += directedForce;

						if (v1.mCollisionForces.empty())
							v1.mCollisionForces.push_back(directedForce);
						else
							v1.mCollisionForces.back() += directedForce;

						if (v2.mCollisionForces.empty())
							v2.mCollisionForces.push_back(directedForce);
						else
							v2.mCollisionForces.back() += directedForce;

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

	//std::cout << std::endl;

	//std::cout << "Num collisions past AABB : " << numCollisionsPastAABB << std::endl;
	//std::cout << "Num accepted collisions : " << numValidCollisions << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
