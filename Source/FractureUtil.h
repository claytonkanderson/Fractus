#pragma once

#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <glm/mat3x3.hpp>
#include <array>

////////////////////////////////////////////////////////////////////////////////

namespace Deformation
{
	glm::ivec2 GetEdgeId(size_t idx1, size_t idx2);
	glm::ivec3 GetFaceId(size_t idx1, size_t idx2, size_t idx3);

	size_t GetCommonVertexFromEdges(const glm::ivec2& edge0, const glm::ivec2& edge1);
    void ComputeMMat(const glm::dvec3& eigenVector, glm::dmat3& outMat);
	double DistPointPlane(const glm::vec3& point, const glm::vec3& planeNormal, const glm::vec3& planePosition);
}

namespace Deformation
{
	static const std::array<glm::ivec3, 12> CubeTriangleIndices ({
	glm::ivec3(0, 4, 3),
	glm::ivec3(0, 2, 1),
	glm::ivec3(1, 2, 5),
	glm::ivec3(0, 3, 2),
	glm::ivec3(2, 3, 4),
	glm::ivec3(2, 4, 5),
	glm::ivec3(0, 7, 4),
	glm::ivec3(1, 5, 6),
	glm::ivec3(4, 7, 5),
	glm::ivec3(5, 7, 6),
	glm::ivec3(0, 6, 7),
	glm::ivec3(0, 1, 6)
	});

	static const std::array<glm::vec3, 8> CubeVertexPositions(
	{
		glm::vec3({0, 0, 0}),
		glm::vec3({1, 0, 0}),
		glm::vec3({1, 1, 0}),
		glm::vec3({0, 1, 0}),
		glm::vec3({0, 1, 1}),
		glm::vec3({1, 1, 1}),
		glm::vec3({1, 0, 1}),
		glm::vec3({0, 0, 1})
	});

	static const std::array<glm::ivec3, 4> TetrahedraIndices (
		{ glm::ivec3(0, 2, 1),
		glm::ivec3(0, 1, 3),
		glm::ivec3(0, 3, 2),
		glm::ivec3(1, 2, 3) }
	);
}

////////////////////////////////////////////////////////////////////////////////