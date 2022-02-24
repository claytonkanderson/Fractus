#include "MeshVolume.h"
#include "FractureUtil.h"

#include <iostream>

////////////////////////////////////////////////////////////////////////////////

bool Deformation::MeshVolume::Tetrahedron::Contains(const glm::vec3& pos) const
{
	for (int i = 0; i < 4; i++)
	{
		const auto& triangle = TetrahedraIndices[i];
		// not sure this is how the normals should be calculated
		const auto& d0 = glm::normalize(mVertices[triangle[1]] - mVertices[triangle[0]]);
		const auto& d1 = glm::normalize(mVertices[triangle[2]] - mVertices[triangle[0]]);
		auto n = glm::cross(d0, d1);
		auto signedDist = DistPointPlane(pos, n, mVertices[triangle[0]]);

		if (signedDist >= cRadius())
			return false;
	}

	return true;
}

bool Deformation::MeshVolume::Contains(const glm::vec3& pos) const
{
	size_t numContainers = 0;

	for (const auto& cubeCenter : mCubeCenters)
	{
		if (CubeContains(pos, cubeCenter))
		{
			numContainers++;
		}
	}

	for (const auto& tet : mTetrahedra)
	{
		if (tet.Contains(pos))
		{
			numContainers++;
		}
	}

	if (numContainers > 1)
	{
		std::cout << " Cube contains! " << pos.x << " " << pos.y << " " << pos.z << " " << std::endl;
		return true;
	}

	return false;
}

bool Deformation::MeshVolume::CubeContains(const glm::vec3& pos, const glm::vec3& cubeCenter) const
{
	auto min = cubeCenter - 0.5f * glm::vec3(1.0f);
	auto max = cubeCenter + 0.5f * glm::vec3(1.0f);

	for (int i = 0; i < 3; i++)
	{
		if (!(min[i] - cRadius() <= pos[i] && pos[i] <= max[i] + cRadius()))
			return false;
	}

	return true;
}

