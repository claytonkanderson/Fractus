#pragma once
#include "Deformation.hpp"

#include <glm/vec3.hpp>
#include <array>

////////////////////////////////////////////////////////////////////////////////

namespace Deformation
{
	class ConvexIntersection
	{
	public:

		void static ResolveCollisions(
			std::vector<Vertex>& vertices,
			std::unordered_map<size_t, Tetrahedra>& tetrahedra);
	};
}

