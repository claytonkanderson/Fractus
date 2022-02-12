#include "TestFractureManager.h"
#include "TestDeformation.hpp"
#include "DeformationAPI.h"
#include "ProtoConverter.hpp"

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

		for (const auto& tetrahedra : group.mTetrahedra)
		{
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

	class TestCaseOne : public TestCase
	{
	public:
		void Run(IronGames::SimulationSummary* summary) override
		{
			int maxNumVertices = 10;
			int maxNumTetrahedra = 10;

			double lambda = 2.65e6f;
			double psi = 397.f;
			double phi = 264.f;
			double mu = 3.97e6f;
			double density = 5013.f;
			double timestep = 0.0001f;
			double toughness = 10.f;

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

			Initialize(positions.data(), 5, maxNumVertices, indices.data(), 2, maxNumTetrahedra, lambda, psi, mu, phi, toughness, density);
			TestDeformation::TetraGroup* group = (TestDeformation::TetraGroup*)mData;

			const glm::vec3 fracturePlane(0.0f, 0.0f, 1.0f);

			{
				SaveFrame(summary, *group);
				InsertFracture(&summary->mutable_frames()->at(0), 0, fracturePlane);
				TestDeformation::FractureContext context(fracturePlane, 0, group->mTetrahedra, group->mVertices);
				context.Fracture();
			}
			{
				SaveFrame(summary, *group);
				InsertFracture(&summary->mutable_frames()->at(1), 1, fracturePlane);
				TestDeformation::FractureContext context(fracturePlane, 1, group->mTetrahedra, group->mVertices);
				context.Fracture();
			}
			{
				SaveFrame(summary, *group);
				InsertFracture(&summary->mutable_frames()->at(2), 2, fracturePlane);
				TestDeformation::FractureContext context(fracturePlane, 2, group->mTetrahedra, group->mVertices);
				context.Fracture();
			}

			SaveFrame(summary, *group);
		}
	};

	class TestCaseTwo : public TestCase
	{
	public:
		void Run(IronGames::SimulationSummary* summary) override
		{
			int maxNumVertices = 10;
			int maxNumTetrahedra = 10;

			double lambda = 2.65e6f;
			double psi = 397.f;
			double phi = 264.f;
			double mu = 3.97e6f;
			double density = 5013.f;
			double timestep = 0.0001f;
			double toughness = 10.f;

			std::vector<float> positions = {
			0,1+ 0, 0,
			1,1+ 0.5f, 0,
			0,1+ 1, 0,
			0,1+ 0.5f, 1,
			-0.5f,1+ 0.5f, -1.2f
			};
			positions.resize(3 * (maxNumVertices + 10));

			std::vector<int> indices = {
				0, 1, 2, 3,
				0, 1, 4, 2
			};
			indices.resize(4 * (maxNumTetrahedra + 10));

			Initialize(positions.data(), 5, maxNumVertices, indices.data(), 2, maxNumTetrahedra, lambda, psi, mu, phi, toughness, density);
			TestDeformation::TetraGroup* group = (TestDeformation::TetraGroup*)mData;

			for (auto& vert : group->mVertices)
				vert.mVelocity += glm::vec3(0, -10, 0);

			group->mSaveEveryXSteps = 1;

			float maxEigenvalue = -1;
			float maxEigenvalueTime = 0.0f;

			for (int i = 0; i < 3000; i++)
			{
				group->Update(timestep);

				for (const auto& vertex : group->mVertices)
				{
					if (vertex.mLargestEigenvalue > maxEigenvalue)
					{
						maxEigenvalue = vertex.mLargestEigenvalue;
						maxEigenvalueTime = i * timestep;
					}
				}
			}

			std::cout << "Max Eigenvalue : " << maxEigenvalue << " time " << maxEigenvalueTime << "s." << std::endl;
		}
	};
};

////////////////////////////////////////////////////////////////////////////////

TestDeformation::TestFractureManager::TestFractureManager(IronGames::SimulationSummaries* summaries)
	: mSummaries(summaries)
{
	mTestCases.push_back(new TestCaseOne());
	mTestCases.push_back(new TestCaseTwo());
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
