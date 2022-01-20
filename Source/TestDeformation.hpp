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
		std::vector<glm::vec3> mCompressiveForces;
		std::vector<glm::vec3> mTensileForces;
	};

	class Tetrahedra
	{
	public:
		
		float mMass = -1;
		float mVolume = -1;
		glm::mat4 mBeta;
		std::array<size_t, 4> mIndices;
	};

	class TetraGroup
	{
	public:

		void Update(float timestep);

		float mLambda = 1e7;
		float mPsi = 100;
		float mMu = 1;
		float mPhi = 1;

		std::vector<Tetrahedra> mTetrahedra;
		std::vector<Vertex> mVertices;
	};
}