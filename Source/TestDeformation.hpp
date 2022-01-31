#define _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS

#include <glm/common.hpp>
#include <glm/mat4x4.hpp>

#include <vector>
#include <array>

////////////////////////////////////////////////////////////////////////////////

namespace TestDeformation
{
	class Vertex
	{
	public:
		
		glm::vec3 mPosition;
		glm::vec3 mMaterialCoordinates;
		glm::vec3 mVelocity;
		float mMass = -1;
		float mInvMass = -1;

		glm::vec3 mForce;
		float mLargestEigenvalue = -1;
		glm::vec3 mPrincipalEigenVector;
		std::vector<glm::vec3> mCompressiveForces;
		std::vector<glm::vec3> mTensileForces;
	};

	class Tetrahedra
	{
	public:
		Tetrahedra() {}
		Tetrahedra(size_t id0, size_t id1, size_t id2, size_t id3)
		{
			mIndices[0] = id0;
			mIndices[1] = id1;
			mIndices[2] = id2;
			mIndices[3] = id3;
		}

		bool ContainsVertexIndex(size_t idx) const;
		bool ContainsEdgeIndex(const glm::ivec2& edgeId) const;
		std::array<glm::ivec2, 6> GetEdges() const;
		void ReplaceVertex(size_t oldVertexId, size_t newVertexId);

	public:
		float mMass = -1;
		float mVolume = -1;
		glm::mat4 mBeta = glm::mat4(0);
		std::array<size_t, 4> mIndices;
	};

	class TetraGroup
	{
	public:

		void Update(float timestep);
		void ComputeDerivedQuantities();
		void FractureNode(size_t nodeIdx, const glm::dvec3& fracturePlaneNormal);
		size_t SplitEdge(const glm::ivec2 & edgeIdx, size_t fractureNodeIdx, const glm::vec3& planeNormal);

		std::vector<size_t> GetTetrahedrasFromNode(size_t nodeIdx) const;
		bool EdgeIntersectsPlane(const glm::ivec2& edgeIdx, int fracutreNodeIdx, const glm::vec3& planeNormal) const;
		bool PlaneIntersectEdge(const glm::vec3& planePos, const glm::vec3& planeNormal,
			const glm::vec3& edgePos0, const glm::vec3& edgePos1, glm::vec3* intersectionPos = nullptr) const;
		float GetSignedDistanceToPlane(int nodeIdx, int fractureNodeIdx, const glm::vec3& planeNormal) const;
		std::vector<size_t> GetTetrahedraNeighbors(size_t tetrahedraIdx) const;
		size_t GetCommonVertexFromEdges(const glm::ivec2& edge0, const glm::ivec2& edge1) const;

	public:
		float mLambda = 1e7;
		float mPsi = 100;
		float mMu = 1;
		float mPhi = 1;
		float mDensity = 1;
		float mToughness = 10000.0f;
		int mMaxNumVertices = 0;
		int mMaxNumTetrahedra = 0;

		std::vector<Tetrahedra> mTetrahedra;
		std::vector<Vertex> mVertices;
	};
}