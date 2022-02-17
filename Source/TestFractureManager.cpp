#include "TestFractureManager.h"
#include "TestDeformation.hpp"
#include "DeformationAPI.h"
#include "ProtoConverter.hpp"
#include <fstream>

////////////////////////////////////////////////////////////////////////////////

namespace TestDeformation
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
	
	// Two tetrahedra with a shared face.
	// Fracture occurs on all nodes in the shared face with planes such that they will be snapped to that face.
	class SimpleCaseThreeFractures : public TestCase
	{
		using TestCase::TestCase;

	public:
		void Run(IronGames::SimulationSummary* summary) override
		{
			std::vector<float> positions = {
			0, 0, 0,
			1, 0.5f, 0,
			0, 1, 0,
			0, 0.5f, 1,
			-0.5f, 0.5f, -1.2f
			};
			positions.resize(3 * (mMaxNumVertices + 10));

			std::vector<int> indices = {
				0, 1, 2, 3,
				0, 1, 4, 2
			};
			indices.resize(4 * (mMaxNumTetrahedra + 10));

			Initialize(positions.data(), 5, mMaxNumVertices, indices.data(), 2, mMaxNumTetrahedra, mLambda, mPsi, mMu, mPhi, mToughness, mDensity);
			TestDeformation::TetraGroup* group = (TestDeformation::TetraGroup*)mData;

			const glm::vec3 fracturePlane(0.0f, 0.0f, 1.0f);

			{
				SaveFrame(summary, *group);
				InsertFracture(&summary->mutable_frames()->at(0), 0, fracturePlane);
				TestDeformation::FractureContext context(fracturePlane, 0, group->mIdToTetrahedra, group->mVertices, group->mTetIdCounter);
				context.Fracture();
			}
			{
				SaveFrame(summary, *group);
				InsertFracture(&summary->mutable_frames()->at(1), 1, fracturePlane);
				TestDeformation::FractureContext context(fracturePlane, 1, group->mIdToTetrahedra, group->mVertices, group->mTetIdCounter);
				context.Fracture();
			}
			{
				SaveFrame(summary, *group);
				InsertFracture(&summary->mutable_frames()->at(2), 2, fracturePlane);
				TestDeformation::FractureContext context(fracturePlane, 2, group->mIdToTetrahedra, group->mVertices, group->mTetIdCounter);
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
			int maxNumVertices = 20;
			int maxNumTetrahedra = 10;

			std::vector<float> positions = {
			0,1+ 0, 0,
			1,1+ 0.5f, 0,
			0,1+ 1, 0,
			0.4f,1+ 0.5f, 1,
			-0.5f,1+ 0.5f, -1.2f
			};
			positions.resize(3 * (maxNumVertices + 10));

			std::vector<int> indices = {
				0, 1, 2, 3,
				0, 1, 4, 2
			};
			indices.resize(4 * (maxNumTetrahedra + 10));

			Initialize(positions.data(), 5, maxNumVertices, indices.data(), 2, maxNumTetrahedra, mLambda, mPsi, mMu, mPhi, mToughness, mDensity);
			TestDeformation::TetraGroup* group = (TestDeformation::TetraGroup*)mData;

			for (auto& vert : group->mVertices)
				vert.mVelocity += glm::vec3(0, -1, 0);

			group->mSaveEveryXSteps = 1;

			float maxEigenvalue = -1;
			float maxEigenvalueTime = 0.0f;

			try {
				for (int i = 0; i < 6000; i++)
				{
					group->Update(mTimestep);

					for (const auto& vertex : group->mVertices)
					{
						if (vertex.mLargestEigenvalue > maxEigenvalue)
						{
							maxEigenvalue = vertex.mLargestEigenvalue;
							maxEigenvalueTime = i * mTimestep;
						}
					}
				}
			}
			catch (const std::exception& e)
			{
				std::cout << "Error occurred during Test Case Two : " << e.what() << std::endl;
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
			file.open("D:/UnityProjects/3D_Template/Assets/Resources/Bowl.obj");

			std::vector<float> positions;
			std::vector<int> indices;

const float positionScale = 2.0f;
const float heightOffset = 2.0f;

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
			TestDeformation::TetraGroup* group = (TestDeformation::TetraGroup*)mData;

			float maxEigenvalue = -1;
			float maxEigenvalueTime = 0.0f;

			//try {
				for (int i = 0; i < 1; i++)
				{
					group->Update(mTimestep);
				}
			//}
			//catch (const std::exception& e)
			//{
			//	std::cout << "Error occurred during Test Case Three : " << e.what() << std::endl;
			//}
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
			int maxNumVertices = 10;
			int maxNumTetrahedra = 10;

			std::vector<float> positions = {
			0, 0, 0,
			1, 0.5f, 0,
			0, 1, 0,
			0, 0.5f, 1,
			-0.5f, 0.5f, -1.2f
			};
			positions.resize(3 * (maxNumVertices + 10));

			std::vector<int> indices = {
				0, 1, 2, 3,
				0, 1, 4, 2
			};
			indices.resize(4 * (maxNumTetrahedra + 10));

			Initialize(positions.data(), 5, maxNumVertices, indices.data(), 2, maxNumTetrahedra, mLambda, mPsi, mMu, mPhi, mToughness, mDensity);
			TestDeformation::TetraGroup* group = (TestDeformation::TetraGroup*)mData;

			const glm::vec3 fracturePlane(0.0f, 1.0f, 0.0f);

			auto fracturingVertId = 3;
			{
				SaveFrame(summary, *group);
				InsertFracture(&summary->mutable_frames()->at(0), fracturingVertId, fracturePlane);
				TestDeformation::FractureContext context(fracturePlane, fracturingVertId, group->mIdToTetrahedra, group->mVertices, group->mTetIdCounter);
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
			int maxNumVertices = 10;
			int maxNumTetrahedra = 10;

			std::vector<float> positions = {
			0, 0, 0,
			1, 0.5f, 0,
			0, 1, 0,
			0, 0.5f, 1,
			-0.5f, 0.5f, -1.2f
			};
			positions.resize(3 * (maxNumVertices + 10));

			std::vector<int> indices = {
				0, 1, 2, 3,
				0, 1, 4, 2
			};
			indices.resize(4 * (maxNumTetrahedra + 10));

			Initialize(positions.data(), 5, maxNumVertices, indices.data(), 2, maxNumTetrahedra, mLambda, mPsi, mMu, mPhi, mToughness, mDensity);
			TestDeformation::TetraGroup* group = (TestDeformation::TetraGroup*)mData;

			const glm::vec3 fracturePlane = glm::normalize(glm::vec3(0.0f, 0.8f, 0.2f));

			auto fracturingVertId = 3;
			{
				SaveFrame(summary, *group);
				InsertFracture(&summary->mutable_frames()->at(0), fracturingVertId, fracturePlane);
				TestDeformation::FractureContext context(fracturePlane, fracturingVertId, group->mIdToTetrahedra, group->mVertices, group->mTetIdCounter);
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
			int maxNumVertices = 10;
			int maxNumTetrahedra = 10;

			std::vector<float> positions = {
			0, 0, 0,
			1, 0.5f, 0,
			0, 1, 0,
			0, 0.5f, 1,
			-0.5f, 0.5f, -1.2f
			};
			positions.resize(3 * (maxNumVertices + 10));

			std::vector<int> indices = {
				0, 1, 2, 3,
				0, 1, 4, 2
			};
			indices.resize(4 * (maxNumTetrahedra + 10));

			Initialize(positions.data(), 5, maxNumVertices, indices.data(), 2, maxNumTetrahedra, mLambda, mPsi, mMu, mPhi, mToughness, mDensity);
			TestDeformation::TetraGroup* group = (TestDeformation::TetraGroup*)mData;

			const glm::vec3 fracturePlane = glm::normalize(glm::vec3(0.0f, 0.8f, 0.5f));

			auto fracturingVertId = 0;
			{
				SaveFrame(summary, *group);
				InsertFracture(&summary->mutable_frames()->at(0), fracturingVertId, fracturePlane);
				TestDeformation::FractureContext context(fracturePlane, fracturingVertId, group->mIdToTetrahedra, group->mVertices, group->mTetIdCounter);
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
			int maxNumVertices = 10;
			int maxNumTetrahedra = 10;

			std::vector<float> positions = {
			0, 0, 0,
			1, 0.5f, 0,
			0, 1, 0,
			0, 0.5f, 1,
			-0.5f, 0.5f, -1.2f
			};
			positions.resize(3 * (maxNumVertices + 10));

			std::vector<int> indices = {
				0, 1, 2, 3,
				0, 1, 4, 2
			};
			indices.resize(4 * (maxNumTetrahedra + 10));

			Initialize(positions.data(), 5, maxNumVertices, indices.data(), 2, maxNumTetrahedra, mLambda, mPsi, mMu, mPhi, mToughness, mDensity);
			TestDeformation::TetraGroup* group = (TestDeformation::TetraGroup*)mData;

			const glm::vec3 fracturePlane = glm::normalize(glm::vec3(-0.2f, 0.8f, 0.0f));

			auto fracturingVertId = 1;
			{
				SaveFrame(summary, *group);
				InsertFracture(&summary->mutable_frames()->at(0), fracturingVertId, fracturePlane);
				TestDeformation::FractureContext context(fracturePlane, fracturingVertId, group->mIdToTetrahedra, group->mVertices, group->mTetIdCounter);
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
			int maxNumVertices = 10;
			int maxNumTetrahedra = 10;

			std::vector<float> positions = {
			0, 0, 0,
			1, 0.5f, 0,
			0, 1, 0,
			0, 0.5f, 1,
			-0.5f, 0.5f, -1.2f
			};
			positions.resize(3 * (maxNumVertices + 10));

			std::vector<int> indices = {
				0, 1, 2, 3,
				0, 1, 4, 2
			};
			indices.resize(4 * (maxNumTetrahedra + 10));

			Initialize(positions.data(), 5, maxNumVertices, indices.data(), 2, maxNumTetrahedra, mLambda, mPsi, mMu, mPhi, mToughness, mDensity);
			TestDeformation::TetraGroup* group = (TestDeformation::TetraGroup*)mData;

			double speed = 1;
			group->mVertices[3].mVelocity += glm::vec3(0, 0, speed);
			group->mVertices[4].mVelocity += glm::vec3(0, 0, -speed);

			group->mSaveEveryXSteps = 1;

			float maxEigenvalue = -1;
			float maxEigenvalueTime = 0.0f;

			try {
				for (int i = 0; i < 100; i++)
				{
					group->Update(mTimestep);

					for (const auto& vertex : group->mVertices)
					{
						if (vertex.mLargestEigenvalue > maxEigenvalue)
						{
							maxEigenvalue = vertex.mLargestEigenvalue;
							maxEigenvalueTime = i * mTimestep;
						}
					}
				}
			}
			catch (const std::exception& e)
			{
				std::cout << "Error occurred during Test Case Two : " << e.what() << std::endl;
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
			std::vector<float> positions = {
			0, 0, 0,
			1, 0.5f, 0,
			0, 1, 0,
			0.4f, 0.5f, 1
			};
			positions.resize(3 * (mMaxNumVertices + 10));

			std::vector<int> indices = {
				0, 1, 2, 3
			};
			indices.resize(4 * (mMaxNumTetrahedra + 10));

			Initialize(positions.data(), 4, mMaxNumVertices, indices.data(), 1, mMaxNumTetrahedra, mLambda, mPsi, mMu, mPhi, mToughness, mDensity);
			TestDeformation::TetraGroup* group = (TestDeformation::TetraGroup*)mData;

			double speed = 1;
			group->mVertices[3].mVelocity += glm::vec3(0, 0, speed);

			group->mSaveEveryXSteps = 1;

			float maxEigenvalue = -1;
			float maxEigenvalueTime = 0.0f;

			try {
				for (int i = 0; i < 100; i++)
				{
					group->Update(mTimestep);

					for (const auto& vertex : group->mVertices)
					{
						if (vertex.mLargestEigenvalue > maxEigenvalue)
						{
							maxEigenvalue = vertex.mLargestEigenvalue;
							maxEigenvalueTime = i * mTimestep;
						}
					}
				}
			}
			catch (const std::exception& e)
			{
				std::cout << "Error occurred during Test Case Two : " << e.what() << std::endl;
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
			int maxNumVertices = 10;
			int maxNumTetrahedra = 10;

			std::vector<float> positions = {
			0, 0, 0,
			1, 0.5f, 0,
			0, 1, 0,
			0.4f, 0.5f, 1
			};
			positions.resize(3 * (maxNumVertices + 10));

			std::vector<int> indices = {
				0, 1, 2, 3
			};
			indices.resize(4 * (maxNumTetrahedra + 10));

			Initialize(positions.data(), 4, maxNumVertices, indices.data(), 1, maxNumTetrahedra, mLambda, mPsi, mMu, mPhi, mToughness, mDensity);
			TestDeformation::TetraGroup* group = (TestDeformation::TetraGroup*)mData;

			double speed = 1;
			group->mVertices[3].mVelocity += glm::vec3(0, 0, -speed);

			group->mSaveEveryXSteps = 1;

			float maxEigenvalue = -1;
			float maxEigenvalueTime = 0.0f;

			try {
				for (int i = 0; i < 100; i++)
				{
					group->Update(mTimestep);

					for (const auto& vertex : group->mVertices)
					{
						if (vertex.mLargestEigenvalue > maxEigenvalue)
						{
							maxEigenvalue = vertex.mLargestEigenvalue;
							maxEigenvalueTime = i * mTimestep;
						}
					}
				}
			}
			catch (const std::exception& e)
			{
				std::cout << "Error occurred during Test Case Two : " << e.what() << std::endl;
			}

			*summary = group->mSummary;
			std::cout << "Max Eigenvalue : " << maxEigenvalue << " time " << maxEigenvalueTime << "s." << std::endl;
		}
	};
};

////////////////////////////////////////////////////////////////////////////////

TestDeformation::TestFractureManager::TestFractureManager(IronGames::SimulationSummaries* summaries)
	: mSummaries(summaries)
{
	bool breakOnFailure = true;

	mTestCases.push_back(new SimpleCaseThreeFractures(breakOnFailure));
	mTestCases.push_back(new UpdateCaseTwoTetrahedra(breakOnFailure));

	mTestCases.push_back(new UpdateCaseBowl(breakOnFailure));
	
	mTestCases.push_back(new SimpleCaseEdgeSnappedFracture(breakOnFailure));
	mTestCases.push_back(new SimpleCaseRegularFracture(breakOnFailure));
	mTestCases.push_back(new SimpleCaseSharedNodeSingleFracture(breakOnFailure));
	mTestCases.push_back(new SimpleCaseSharedNodeRegularFracture(breakOnFailure));
	mTestCases.push_back(new UpdateCaseTwoTetVelocityTensile(breakOnFailure));
	mTestCases.push_back(new SimpleCaseSingleTetTensile(breakOnFailure));
	mTestCases.push_back(new SimpleCaseSingleTetCompressive(breakOnFailure));
}

////////////////////////////////////////////////////////////////////////////////

TestDeformation::TestFractureManager::~TestFractureManager()
{
	for (auto testCase : mTestCases)
		delete testCase;
}

////////////////////////////////////////////////////////////////////////////////

void TestDeformation::TestFractureManager::RunAllTestCases()
{
	for (int i = 0; i < mTestCases.size(); i++)
		RunTestCase(i);
}

////////////////////////////////////////////////////////////////////////////////

void TestDeformation::TestFractureManager::RunTestCase(int testNum)
{
	mTestCases[testNum]->Run(mSummaries->add_summaries());
}

////////////////////////////////////////////////////////////////////////////////
