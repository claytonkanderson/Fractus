#define _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS

#include "core.pb.h"

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
		
		glm::dvec3 mPosition;
		glm::dvec3 mMaterialCoordinates;
		glm::dvec3 mVelocity;
		double mMass = -1;
		double mInvMass = -1;

		glm::dvec3 mForce;
		double mLargestEigenvalue = -1;
		glm::dvec3 mPrincipalEigenVector;
		std::vector<glm::dvec3> mCompressiveForces;
		std::vector<glm::dvec3> mTensileForces;
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

		std::array<size_t, 3> GetOtherVertices(size_t vertexId) const;
		std::array<size_t, 2> GetOtherVertices(size_t vertexId1, size_t vertexId2) const;
		bool ContainsVertexIndex(size_t idx) const;
		bool ContainsEdgeIndex(const glm::ivec2& edgeId) const;
		std::array<glm::ivec2, 6> GetEdges() const;
		void ReplaceVertex(size_t oldVertexId, size_t newVertexId);

	public:
		double mMass = -1;
		double mVolume = -1;
		glm::dmat4 mBeta = glm::dmat4(0);
		std::array<size_t, 4> mIndices;
	};

	class FractureContext
	{
	public:
		FractureContext(
			const glm::vec3 & fracturePlaneNormal,
			const glm::vec3 & fracturePlanePosition,
			size_t fractureVertexIdx,
			std::vector<Tetrahedra> & tetrahedra,
			std::vector<Vertex> & vertices)
			: 
			mFracturePlaneNormal(fracturePlaneNormal),
			mFractureNodePosition(fracturePlanePosition),
			mFractureNodeIdx(fractureVertexIdx),
			mTetrahedra(tetrahedra),
			mVertices(vertices)
		{

		}

		void Fracture();

	private:
		void DetermineSnapping();

		void FaceSnappedFracture();
		void EdgeSnappedFracture();

		std::vector<glm::ivec2> PlaneIntersectTetrahedraEdges() const;

		void FractureTetrahedra(size_t tetIdx);
		void AssignTetToSide();
		void RegularFracture();
		void NeighborFaceFracture();
		void NeighborEdgeFracture();

	private:
		size_t mTetrahedraIdx = -1;

		bool mSnapToFace = false;
		glm::ivec3 mSnappingFaceId;

		bool mSnapToEdge = false;
		glm::ivec2 mSnappingEdgeId;

		const glm::dvec3& mFracturePlaneNormal;
		const glm::dvec3& mFractureNodePosition;
		size_t mFractureNodeIdx = -1;
		std::vector<Tetrahedra>& mTetrahedra;
		std::vector<Vertex>& mVertices;
	};

	class TetraGroup
	{
	public:
		void Update(double timestep);
		void ComputeDerivedQuantities();
		void FractureNode(size_t nodeIdx, const glm::dvec3& fracturePlaneNormal);
		size_t SplitEdge(const glm::ivec2 & edgeIdx, size_t fractureNodeIdx, const glm::dvec3& planeNormal);

		std::vector<size_t> GetTetrahedrasFromNode(size_t nodeIdx) const;
		bool EdgeIntersectsPlane(const glm::ivec2& edgeIdx, int fractureNodeIdx, const glm::dvec3& planeNormal) const;
		bool PlaneIntersectEdge(const glm::dvec3& planePos, const glm::dvec3& planeNormal,
			const glm::dvec3& edgePos0, const glm::dvec3& edgePos1, double & d, glm::dvec3* intersectionPos = nullptr) const;
		double GetSignedDistanceToPlane(int nodeIdx, int fractureNodeIdx, const glm::dvec3& planeNormal) const;
		std::vector<size_t> GetTetrahedraNeighbors(size_t tetrahedraIdx) const;
		size_t GetCommonVertexFromEdges(const glm::ivec2& edge0, const glm::ivec2& edge1) const;

	public:
		double mLambda = 1e7;
		double mPsi = 100;
		double mMu = 1;
		double mPhi = 1;
		double mDensity = 1;
		double mToughness = 10000.0f;
		int mMaxNumVertices = 0;
		int mMaxNumTetrahedra = 0;
		int mStepNum = 0;
		int mStepsSinceLastSave = 100000;
		int mSaveEveryXSteps = 0;
		double mSimulationTime = 0;

		std::vector<Tetrahedra> mTetrahedra;
		std::vector<Vertex> mVertices;

		IronGames::SimulationSummary* mSummary = nullptr;
	};
}