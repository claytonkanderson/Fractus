#pragma once

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
		size_t GetOtherVertex(const std::array<size_t, 3>& vertices) const;
		bool ContainsVertexIndex(size_t idx) const;
		bool ContainsEdgeIndex(const glm::ivec2& edgeId) const;
		bool ContainsFaceIndex(const glm::ivec3& faceId) const;
		std::array<glm::ivec2, 6> GetEdges() const;
		void ReplaceVertex(size_t oldVertexId, size_t newVertexId);

	public:
		bool mFracturedThisFrame = false;

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
			size_t fractureVertexIdx,
			std::unordered_map<size_t, Tetrahedra>& idToTetrahedra,
			std::vector<Vertex> & vertices,
			size_t & tetIdCounter
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
		size_t CreateEdgeVertex(const glm::ivec2 & edgeId, double parametricDistance);
		
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
		size_t & mTetIdCounter;
		std::unordered_map<size_t, Tetrahedra>& mIdToTetrahedra;
		std::unordered_set<size_t> mTetrahedraIdsToDelete;
		std::vector<Vertex>& mVertices;
	};

	class TetraGroup
	{
	public:
		
		void Initialize(double lambda, double psi, double phi, double mu, double density, double toughness);
		void Update(double timestep);
		void ComputeDerivedQuantities();
		void FractureNode(size_t nodeIdx, const glm::dvec3& fracturePlaneNormal);
		size_t SplitEdge(const glm::ivec2 & edgeIdx, size_t fractureNodeIdx, const glm::dvec3& planeNormal);
		void OutputSaveFile();

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
		int mSaveEveryXSteps = 1;
		double mSimulationTime = 0;

		size_t mTetIdCounter = 0;
		std::unordered_map<size_t, Tetrahedra> mIdToTetrahedra;
		std::vector<Vertex> mVertices;

		IronGames::SimulationSummary mSummary;
	};
}