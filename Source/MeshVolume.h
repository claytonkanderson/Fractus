#pragma once

#include <array>
#include <vector>
#include <glm/vec3.hpp>

////////////////////////////////////////////////////////////////////////////////

namespace Deformation
{
	class MeshVolume
	{
	public:
		static float cRadius() { return 0.01f; };

		class Tetrahedron
		{
		public:
			bool Contains(const glm::vec3& pos) const;
		public:
			std::array<glm::vec3, 4> mVertices;
		};

		// Contains means that a small sphere at the position is fully inside the volume.
		bool Contains(const glm::vec3& pos) const;

	private:
		bool CubeContains(const glm::vec3& pos, const glm::vec3& cubeCenter) const;

	public:
		std::vector<glm::vec3> mCubeCenters;
		std::vector<Tetrahedron> mTetrahedra;
	};
}

////////////////////////////////////////////////////////////////////////////////