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

		struct CollisionResults
		{
			float mVolume = 0;
			glm::vec3 mForce;
		};

		bool static TetrahedraIntersection(
			const std::array<glm::vec3, 4>& tetVertices,
			const std::array<glm::vec3, 4>& otherVertices,
			CollisionResults& results);

		void static ResolveCollisions(
			std::vector<Vertex>& vertices,
			std::unordered_map<size_t, Tetrahedra>& tetrahedra);
	};
}

