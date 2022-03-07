#include "TestFractureManager.h"
#include "Deformation.hpp"
#include "DeformationAPI.h"
#include "FractureContext.h"
#include "ProtoConverter.hpp"
#include <fstream>

////////////////////////////////////////////////////////////////////////////////

namespace Deformation
{
	void SaveFrame(IronGames::SimulationSummary* summary, const TetraGroup& group)
	{
		auto frame = summary->add_frames();

		for (size_t i = 0; i < group.mVertices.size(); i++)
		{
			const auto& vertex = group.mVertices[i];
			auto vert = frame->add_vertices();
			*vert->mutable_position() = ProtoConverter::Convert(vertex.mPosition);
			*vert->mutable_material_coordinates() = ProtoConverter::Convert(vertex.mMaterialCoordinates);
			*vert->mutable_velocity() = ProtoConverter::Convert(vertex.mVelocity);
			vert->set_mass(vertex.mMass);
		}

		for (const auto& pair : group.mIdToTetrahedra)
		{
			const auto& tetrahedra = pair.second;
			auto& tet = *frame->add_tetrahedra();
			tet.set_mass(tetrahedra.mMass);
			tet.set_volume(tetrahedra.mVolume);
			for (auto& idx : tetrahedra.mIndices)
				tet.add_indices(idx);
		}
	}

	void InsertFracture(IronGames::SimulationFrame* frame, size_t vertId, const glm::vec3& fracturePlane)
	{
		double largeEigenValue = 1000;

		for (size_t i = 0; i < frame->mutable_vertices()->size(); i++)
		{
			if (i == vertId)
			{
				auto &vertex = (*frame->mutable_vertices())[i];
				vertex.set_largest_eigenvalue(largeEigenValue);
				*vertex.mutable_principal_eigenvector() = ProtoConverter::Convert(fracturePlane);
			}
		}
	}

	void GetTwoTetrahedraConfiguration(std::vector<float>& outPositions, std::vector<int>& outIndices, float yOffset = 0.0f)
	{
		outPositions = {
			0, yOffset + 0, 0,
			1, yOffset + 0.5f, 0,
			0, yOffset + 1, 0,
			0.4f, yOffset + 0.5f, 1,
			0.45f, yOffset + 0.5f, -1.2f
		};

		outIndices = {
			0, 1, 2, 3,
			0, 1, 4, 2
		};
	}
	
	void GetSingleTetrahedraConfiguration(std::vector<float>& outPositions, std::vector<int>& outIndices, float yOffset = 0.0f)
	{
		outPositions = {
			0, yOffset + 0, 0,
			1, yOffset + 0.5f, 0,
			0, yOffset + 1, 0,
			0.4f, yOffset + 0.5f, 1,
		};

		outIndices = {
			0, 1, 2, 3
		};
	}

	// Two tetrahedra with a shared face.
	// Fracture occurs on all nodes in the shared face with planes such that they will be snapped to that face.
	class SimpleCaseThreeFractures : public TestCase
	{
		using TestCase::TestCase;

	public:
		void Run(IronGames::SimulationSummary* summary) override
		{
			std::vector<float> positions;
			std::vector<int> indices;
			GetTwoTetrahedraConfiguration(positions, indices);
			positions.resize(3 * (mMaxNumVertices + 10));
			indices.resize(4 * (mMaxNumTetrahedra + 10));

			Initialize(positions.data(), 5, mMaxNumVertices, indices.data(), 2, mMaxNumTetrahedra, mLambda, mPsi, mMu, mPhi, mToughness, mDensity);
			Deformation::TetraGroup* group = (Deformation::TetraGroup*)mData;

			const glm::vec3 fracturePlane(0.0f, 0.0f, 1.0f);

			{
				SaveFrame(summary, *group);
				InsertFracture(&summary->mutable_frames()->at(0), 0, fracturePlane);
				Deformation::FractureContext context(fracturePlane, 0, group->mIdToTetrahedra, group->mVertices, group->mTetIdCounter);
				context.Fracture();
			}
			{
				SaveFrame(summary, *group);
				InsertFracture(&summary->mutable_frames()->at(1), 1, fracturePlane);
				Deformation::FractureContext context(fracturePlane, 1, group->mIdToTetrahedra, group->mVertices, group->mTetIdCounter);
				context.Fracture();
			}
			{
				SaveFrame(summary, *group);
				InsertFracture(&summary->mutable_frames()->at(2), 2, fracturePlane);
				Deformation::FractureContext context(fracturePlane, 2, group->mIdToTetrahedra, group->mVertices, group->mTetIdCounter);
				context.Fracture();
			}

			SaveFrame(summary, *group);
		}
	};

	class UpdateCaseTwoTetrahedra : public TestCase
	{
		using TestCase::TestCase;

	public:

		void Run(IronGames::SimulationSummary* summary) override
		{
			mMaxNumVertices = 20;
			mMaxNumTetrahedra = 10;
			mToughness = 100;

			std::vector<float> positions;
			std::vector<int> indices;
			GetTwoTetrahedraConfiguration(positions, indices, 1.0f);
			positions.resize(3 * (mMaxNumVertices + 10));
			indices.resize(4 * (mMaxNumTetrahedra + 10));

			Initialize(positions.data(), 5, mMaxNumVertices, indices.data(), 2, mMaxNumTetrahedra, mLambda, mPsi, mMu, mPhi, mToughness, mDensity);
			Deformation::TetraGroup* group = (Deformation::TetraGroup*)mData;

			for (auto& vert : group->mVertices)
				vert.mVelocity += glm::vec3(0, -1, 0);

			group->mSaveEveryXSteps = 1;

			float maxEigenvalue = -1;
			float maxEigenvalueTime = 0.0f;

			for (int i = 0; i < 6000; i++)
			{
				if (i == 3684)
					std::cout << "Here" << std::endl;

				if (!Update(*group))
					break;

				for (const auto& vertex : group->mVertices)
				{
					if (vertex.mLargestEigenvalue > maxEigenvalue)
					{
						maxEigenvalue = vertex.mLargestEigenvalue;
						maxEigenvalueTime = i * mTimestep;
					}
				}
			}

			*summary = group->mSummary;
			std::cout << "Max Eigenvalue : " << maxEigenvalue << " time " << maxEigenvalueTime << "s." << std::endl;
		}
	};

	class UpdateCaseBowl : public TestCase
	{
		using TestCase::TestCase;

	public:
		void Run(IronGames::SimulationSummary* summary) override
		{
			std::ifstream file;
			//file.open("D:/UnityProjects/3D_Template/Assets/Resources/Bowl.obj");
			file.open("D:/UnityProjects/3D_Template/Assets/Resources/BasicMace.obj");

			std::vector<float> positions;
			std::vector<int> indices;

const float positionScale = 1.0f;
const float heightOffset = 8.0f;

			std::string line;
			float x, y, z;
			int i0, i1, i2, i3;
			char c;
			while (std::getline(file, line))
			{
				std::istringstream iss(line);

				iss >> c;

				if (c == 'v')
				{
					iss >> x >> y >> z;

					positions.push_back(positionScale * x);
					positions.push_back(positionScale * y + heightOffset);
					positions.push_back(positionScale * z);
				}
				else if (c == 't')
				{
					iss >> i0 >> i1 >> i2 >> i3;
					indices.push_back(i0 - 1);
					indices.push_back(i1 - 1);
					indices.push_back(i2 - 1);
					indices.push_back(i3 - 1);
				}
			}

			int numVertices = positions.size() / 3;
			int numTetrahedra = indices.size() / 4;

			int maxNumVertices = 2 * numVertices;
			int maxNumTetrahedra = 2 * numTetrahedra;

			Initialize(positions.data(), numVertices, maxNumVertices, indices.data(), numTetrahedra, maxNumTetrahedra, mLambda, mPsi, mMu, mPhi, mToughness, mDensity);
			Deformation::TetraGroup* group = (Deformation::TetraGroup*)mData;

			float maxEigenvalue = -1;
			float maxEigenvalueTime = 0.0f;

			SaveFrame(summary, *group);

			int numSteps = 10000;
			for (int i = 0; i < numSteps; i++)
			{
				if (!Update(*group))
					break;

				std::cout << "Completed step " << i << " out of " << numSteps << " total steps." << std::endl;
			}

			//*summary = group->mSummary;

			SaveFrame(summary, *group);
		}
	};

	// Two tetrahedra with a shared face.
	// A single fracture event occurs on a non-shared node to split the shared face in two with edge snapping.
	// so that 4 tetrahedra are expected as the result.
	class SimpleCaseEdgeSnappedFracture : public TestCase
	{
		using TestCase::TestCase;
	public:
		void Run(IronGames::SimulationSummary* summary) override
		{
			std::vector<float> positions;
			std::vector<int> indices;
			GetTwoTetrahedraConfiguration(positions, indices);

			Initialize(positions.data(), 5, mMaxNumVertices, indices.data(), 2, mMaxNumTetrahedra, mLambda, mPsi, mMu, mPhi, mToughness, mDensity);
			Deformation::TetraGroup* group = (Deformation::TetraGroup*)mData;

			const glm::vec3 fracturePlane(0.0f, 1.0f, 0.0f);

			auto fracturingVertId = 3;
			{
				SaveFrame(summary, *group);
				InsertFracture(&summary->mutable_frames()->at(0), fracturingVertId, fracturePlane);
				Deformation::FractureContext context(fracturePlane, fracturingVertId, group->mIdToTetrahedra, group->mVertices, group->mTetIdCounter);
				context.Fracture();
			}

			SaveFrame(summary, *group);
		}
	};
	
	// Two tetrahedra with a shared face.
	// A single fracture event occurs on a non-shared node to split the shared face in two.
	// so that 6 tetrahedra are expected as the result.
	class SimpleCaseRegularFracture : public TestCase
	{
		using TestCase::TestCase;
	public:
		void Run(IronGames::SimulationSummary* summary) override
		{
			std::vector<float> positions;
			std::vector<int> indices;
			GetTwoTetrahedraConfiguration(positions, indices);

			Initialize(positions.data(), 5, mMaxNumVertices, indices.data(), 2, mMaxNumTetrahedra, mLambda, mPsi, mMu, mPhi, mToughness, mDensity);
			Deformation::TetraGroup* group = (Deformation::TetraGroup*)mData;

			const glm::vec3 fracturePlane = glm::normalize(glm::vec3(0.0f, 0.8f, 0.2f));

			auto fracturingVertId = 3;
			{
				SaveFrame(summary, *group);
				InsertFracture(&summary->mutable_frames()->at(0), fracturingVertId, fracturePlane);
				Deformation::FractureContext context(fracturePlane, fracturingVertId, group->mIdToTetrahedra, group->mVertices, group->mTetIdCounter);
				context.Fracture();
			}

			SaveFrame(summary, *group);
		}
	};

	// Two tetrahedra with a shared face.
	// A single fracture event occurs on a shared node that does not split the shared face.
	// The plane does not intersect one of the tetrahedra.
	class SimpleCaseSharedNodeSingleFracture : public TestCase
	{
		using TestCase::TestCase;
	public:
		void Run(IronGames::SimulationSummary* summary) override
		{
			std::vector<float> positions;
			std::vector<int> indices;
			GetTwoTetrahedraConfiguration(positions, indices);

			Initialize(positions.data(), 5, mMaxNumVertices, indices.data(), 2, mMaxNumTetrahedra, mLambda, mPsi, mMu, mPhi, mToughness, mDensity);
			Deformation::TetraGroup* group = (Deformation::TetraGroup*)mData;

			const glm::vec3 fracturePlane = glm::normalize(glm::vec3(0.0f, 0.8f, 0.5f));

			auto fracturingVertId = 0;
			{
				SaveFrame(summary, *group);
				InsertFracture(&summary->mutable_frames()->at(0), fracturingVertId, fracturePlane);
				Deformation::FractureContext context(fracturePlane, fracturingVertId, group->mIdToTetrahedra, group->mVertices, group->mTetIdCounter);
				context.Fracture();
			}

			SaveFrame(summary, *group);
		}
	};

	// Two tetrahedra with a shared face.
	// A single fracture event occurs on a shared node that splits the shared face.
	// Six tetrahedra are expected.
	class SimpleCaseSharedNodeRegularFracture : public TestCase
	{
		using TestCase::TestCase;
	public:
		void Run(IronGames::SimulationSummary* summary) override
		{
			std::vector<float> positions;
			std::vector<int> indices;
			GetTwoTetrahedraConfiguration(positions, indices);

			Initialize(positions.data(), 5, mMaxNumVertices, indices.data(), 2, mMaxNumTetrahedra, mLambda, mPsi, mMu, mPhi, mToughness, mDensity);
			Deformation::TetraGroup* group = (Deformation::TetraGroup*)mData;

			const glm::vec3 fracturePlane = glm::normalize(glm::vec3(-0.2f, 0.8f, 0.0f));

			auto fracturingVertId = 1;
			{
				SaveFrame(summary, *group);
				InsertFracture(&summary->mutable_frames()->at(0), fracturingVertId, fracturePlane);
				Deformation::FractureContext context(fracturePlane, fracturingVertId, group->mIdToTetrahedra, group->mVertices, group->mTetIdCounter);
				context.Fracture();
			}

			SaveFrame(summary, *group);
		}
	};

	class UpdateCaseTwoTetVelocityTensile : public TestCase
	{
		using TestCase::TestCase;
	public:
		void Run(IronGames::SimulationSummary* summary) override
		{
			std::vector<float> positions;
			std::vector<int> indices;
			GetTwoTetrahedraConfiguration(positions, indices);

			Initialize(positions.data(), 5, mMaxNumVertices, indices.data(), 2, mMaxNumTetrahedra, mLambda, mPsi, mMu, mPhi, mToughness, mDensity);
			Deformation::TetraGroup* group = (Deformation::TetraGroup*)mData;

			mTimestep = 0.001f;

			double speed = 1;
			group->mVertices[3].mVelocity += glm::vec3(0, 0, speed);
			group->mVertices[4].mVelocity += glm::vec3(0, 0, -speed);

			group->mSaveEveryXSteps = 1;

			float maxEigenvalue = -1;
			float maxEigenvalueTime = 0.0f;

			for (int i = 0; i < 10000; i++)
			{
				if (!Update(*group))
					break;

				for (const auto& vertex : group->mVertices)
				{
					if (vertex.mLargestEigenvalue > maxEigenvalue)
					{
						maxEigenvalue = vertex.mLargestEigenvalue;
						maxEigenvalueTime = i * mTimestep;
					}
				}
			}

			*summary = group->mSummary;
			std::cout << "Max Eigenvalue : " << maxEigenvalue << " time " << maxEigenvalueTime << "s." << std::endl;
		}
	};

	class SimpleCaseSingleTetTensile : public TestCase
	{
		using TestCase::TestCase;
	public:
		void Run(IronGames::SimulationSummary* summary) override
		{
			std::vector<float> positions;
			std::vector<int> indices;
			GetSingleTetrahedraConfiguration(positions, indices);

			Initialize(positions.data(), 4, mMaxNumVertices, indices.data(), 1, mMaxNumTetrahedra, mLambda, mPsi, mMu, mPhi, mToughness, mDensity);
			Deformation::TetraGroup* group = (Deformation::TetraGroup*)mData;

			double speed = 1;
			group->mVertices[3].mVelocity += glm::vec3(0, 0, speed);

			group->mSaveEveryXSteps = 1;

			float maxEigenvalue = -1;
			float maxEigenvalueTime = 0.0f;

			for (int i = 0; i < 100; i++)
			{
				if (!Update(*group))
					break;

				for (const auto& vertex : group->mVertices)
				{
					if (vertex.mLargestEigenvalue > maxEigenvalue)
					{
						maxEigenvalue = vertex.mLargestEigenvalue;
						maxEigenvalueTime = i * mTimestep;
					}
				}
			}

			*summary = group->mSummary;
			std::cout << "Max Eigenvalue : " << maxEigenvalue << " time " << maxEigenvalueTime << "s." << std::endl;
		}
	};

	class SimpleCaseSingleTetCompressive : public TestCase
	{
		using TestCase::TestCase;
	public:
		void Run(IronGames::SimulationSummary* summary) override
		{
			std::vector<float> positions;
			std::vector<int> indices;
			GetSingleTetrahedraConfiguration(positions, indices);

			Initialize(positions.data(), 4, mMaxNumVertices, indices.data(), 1, mMaxNumTetrahedra, mLambda, mPsi, mMu, mPhi, mToughness, mDensity);
			Deformation::TetraGroup* group = (Deformation::TetraGroup*)mData;

			double speed = 1;
			group->mVertices[3].mVelocity += glm::vec3(0, 0, -speed);

			group->mSaveEveryXSteps = 1;

			float maxEigenvalue = -1;
			float maxEigenvalueTime = 0.0f;

			for (int i = 0; i < 100; i++)
			{
				if (!Update(*group))
					break;

				for (const auto& vertex : group->mVertices)
				{
					if (vertex.mLargestEigenvalue > maxEigenvalue)
					{
						maxEigenvalue = vertex.mLargestEigenvalue;
						maxEigenvalueTime = i * mTimestep;
					}
				}
			}

			*summary = group->mSummary;
			std::cout << "Max Eigenvalue : " << maxEigenvalue << " time " << maxEigenvalueTime << "s." << std::endl;
		}
	};

	class UpdateCaseTwoTetVelocityBowTie : public TestCase
	{
		using TestCase::TestCase;
	public:
		void Run(IronGames::SimulationSummary* summary) override
		{
			std::vector<float> positions = {
				-0.5f, 0.0f, 1.0f,
				0.5f, 0.0f, 1.0f,
				0.0f, 1.0f, 1.0f,
				0.0f, 0.5f, 0.0f,
				-0.5f, 0.0f, -1.0f,
				0.5f, 0.0f, -1.0f,
				0.0f, 1.0f, -1.0f,
			};

			std::vector<int> indices = {
				0, 1, 2, 3,
				3, 4, 5, 6
			};

			Initialize(positions.data(), 7, mMaxNumVertices, indices.data(), 2, mMaxNumTetrahedra, mLambda, mPsi, mMu, mPhi, mToughness, mDensity);
			Deformation::TetraGroup* group = (Deformation::TetraGroup*)mData;

			group->mVertices[0].mVelocity += glm::vec3(0, 0, 1);
			group->mSaveEveryXSteps = 1;

			float maxEigenvalue = -1;
			float maxEigenvalueTime = 0.0f;

			for (int i = 0; i < 100; i++)
			{
				if (!Update(*group))
					break;

				for (const auto& vertex : group->mVertices)
				{
					if (vertex.mLargestEigenvalue > maxEigenvalue)
					{
						maxEigenvalue = vertex.mLargestEigenvalue;
						maxEigenvalueTime = i * mTimestep;
					}
				}
			}

			*summary = group->mSummary;
			std::cout << "Max Eigenvalue : " << maxEigenvalue << " time " << maxEigenvalueTime << "s." << std::endl;
		}
	};
	
	////////////////////////////////////////////////////////////////////////////////

	class UpdateCaseTwoTetCollision : public TestCase
	{
		using TestCase::TestCase;
	public:
		void Run(IronGames::SimulationSummary* summary) override
		{
			std::vector<float> positions = {
				-0.5f, 0.0f, 1.0f,
				0.5f, 0.0f, 1.0f,
				0.0f, 1.0f, 1.0f,
				0.0f, 0.5f, 0.0f,
				-0.5f, 0.0f, -1.0f,
				0.5f, 0.0f, -1.0f,
				0.0f, 1.0f, -1.0f,
				0.0f, 0.5f, 0.5f
			};

			std::vector<int> indices = {
				0, 1, 2, 3,
				4, 5, 6, 7
			};

			Initialize(positions.data(), 8, mMaxNumVertices, indices.data(), 2, mMaxNumTetrahedra, mLambda, mPsi, mMu, mPhi, mToughness, mDensity);
			Deformation::TetraGroup* group = (Deformation::TetraGroup*)mData;

			group->mVertices[0].mVelocity += glm::vec3(0, 0, 1);
			group->mSaveEveryXSteps = 1;

			float maxEigenvalue = -1;
			float maxEigenvalueTime = 0.0f;

			for (int i = 0; i < 10000; i++)
			{
				if (!Update(*group))
					break;

				for (const auto& vertex : group->mVertices)
				{
					if (vertex.mLargestEigenvalue > maxEigenvalue)
					{
						maxEigenvalue = vertex.mLargestEigenvalue;
						maxEigenvalueTime = i * mTimestep;
					}
				}
			}

			*summary = group->mSummary;
			std::cout << "Max Eigenvalue : " << maxEigenvalue << " time " << maxEigenvalueTime << "s." << std::endl;
		}
	};

	////////////////////////////////////////////////////////////////////////////////

	bool TestCase::Update(Deformation::TetraGroup& group)
	{
		if (mBreakOnFailure)
		{
			group.Update(mTimestep);
			return true;
		}

		try {
			group.Update(mTimestep);
		}
		catch (const std::exception& e)
		{
			std::cout << "Error while updating : " << e.what() << std::endl;
			return false;
		}

		return true;
	}
};

////////////////////////////////////////////////////////////////////////////////

Deformation::TestFractureManager::TestFractureManager(IronGames::SimulationSummaries* summaries)
	: mSummaries(summaries)
{
	bool breakOnFailure = true;

	mTestCases.push_back(new SimpleCaseThreeFractures(breakOnFailure));
	mTestCases.push_back(new UpdateCaseTwoTetrahedra(breakOnFailure));
	
	mTestCases.push_back(new SimpleCaseEdgeSnappedFracture(breakOnFailure));
	mTestCases.push_back(new SimpleCaseRegularFracture(breakOnFailure));
	mTestCases.push_back(new SimpleCaseSharedNodeSingleFracture(breakOnFailure));
	mTestCases.push_back(new SimpleCaseSharedNodeRegularFracture(breakOnFailure));
	mTestCases.push_back(new UpdateCaseTwoTetVelocityTensile(breakOnFailure));
	mTestCases.push_back(new SimpleCaseSingleTetTensile(breakOnFailure));
	mTestCases.push_back(new SimpleCaseSingleTetCompressive(breakOnFailure));
	mTestCases.push_back(new UpdateCaseTwoTetVelocityBowTie(breakOnFailure));
	mTestCases.push_back(new UpdateCaseTwoTetCollision(breakOnFailure));

	mTestCases.push_back(new UpdateCaseBowl(breakOnFailure));
}

////////////////////////////////////////////////////////////////////////////////

Deformation::TestFractureManager::~TestFractureManager()
{
	for (auto testCase : mTestCases)
		delete testCase;
}

////////////////////////////////////////////////////////////////////////////////

#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

void Deformation::TestFractureManager::RunAllTestCases()
{
	for (int i = 0; i < mTestCases.size(); i++)
	{
		std::cout << std::endl;
		std::cout << "Beginning test case " << i << "." << std::endl;
		
		auto t1 = Clock::now();
		RunTestCase(i);
		auto t2 = Clock::now();
		auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
		std::cout << "Test Case " << i << " took " << dt << " ms." << std::endl;
	}
}

////////////////////////////////////////////////////////////////////////////////

void Deformation::TestFractureManager::RunTestCase(int testNum)
{
	mTestCases[testNum]->Run(mSummaries->add_summaries());
}

////////////////////////////////////////////////////////////////////////////////
