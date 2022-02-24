#pragma once
#include "Deformation.hpp"
#include <unordered_map>
#include <glm/vec3.hpp>

////////////////////////////////////////////////////////////////////////////////

namespace Deformation
{
	class FractureContext
	{
	public:
		FractureContext(
			const glm::vec3& fracturePlaneNormal,
			size_t fractureVertexIdx,
			std::unordered_map<size_t, Tetrahedra>& idToTetrahedra,
			std::vector<Vertex>& vertices,
			size_t& tetIdCounter
		)
			:
			mFracturePlaneNormal(fracturePlaneNormal),
			mFractureNodePosition(vertices[fractureVertexIdx].mPosition),
			mFractureNodeIdx(fractureVertexIdx),
			mVertices(vertices),
			mTetIdCounter(tetIdCounter),
			mIdToTetrahedra(idToTetrahedra)
		{

		}

		// Returns true if fracturing occurred
		bool Fracture();

	private:
		size_t CloneVertex(size_t vertexId);
		size_t CreateEdgeVertex(const glm::ivec2& edgeId, double parametricDistance);

		void SeparateNodes(size_t& outPositiveNode, size_t& outNegativeNode, const std::array<size_t, 2>& nodes, const glm::dvec3& planePos, const glm::dvec3& planeNormal) const;
		void SeparateNodes(std::vector<size_t>& outPositiveNodes, std::vector<size_t>& outNegativeNodes, const std::vector<size_t>& nodes, const glm::dvec3& planePos, const glm::dvec3& planeNormal) const;
		bool PlaneIntersectEdge(const glm::dvec3& planePos, const glm::dvec3& planeNormal, const glm::dvec3& edgePos0, const glm::dvec3& edgePos1, double& d, glm::dvec3* intersectionPos) const;
		bool IsIsolatedEdge(const glm::ivec2& edgeId) const;
		std::vector<size_t> GetTetrahedraNeighbors(size_t nodeIdx) const;
		size_t GetNonFractureNode(const glm::ivec2& edge) const;

	private:
		std::vector<Tetrahedra> mNewTetrahedra;

		glm::dvec3 mFracturePlaneNormal;
		const glm::dvec3 mFractureNodePosition;
		size_t mFractureNodeIdx = -1;
		size_t& mTetIdCounter;
		std::unordered_map<size_t, Tetrahedra>& mIdToTetrahedra;
		std::unordered_set<size_t> mTetrahedraIdsToDelete;
		std::vector<Vertex>& mVertices;
	};
}

