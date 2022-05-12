#include "TestFractureManager.h"
#include "Deformation.hpp"
#include "DeformationAPI.h"
#include "FractureContext.h"
#include "ProtoConverter.hpp"
#include "ImplicitIntegrator.h"
#include <fstream>
#include <sstream>
#include <chrono>

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

			*vert->mutable_force() = ProtoConverter::Convert(vertex.mForce);
			vert->set_largest_eigenvalue(vertex.mLargestEigenvalue);
			*vert->mutable_principal_eigenvector() = ProtoConverter::Convert(vertex.mPrincipalEigenVector);

			for (const auto& compressiveForce : vertex.mCompressiveForces)
				*vert->add_compressive_forces() = ProtoConverter::Convert(compressiveForce);
			for (const auto& tensileForce : vertex.mTensileForces)
				*vert->add_tensile_forces() = ProtoConverter::Convert(tensileForce);
			for (const auto& collisionForce : vertex.mCollisionForces)
				*vert->add_collision_forces() = ProtoConverter::Convert(collisionForce);
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
		void Run(IronGames::SimulationSummaries* summaries) override
		{
			auto summary = summaries->add_summaries();
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

		void Run(IronGames::SimulationSummaries* summaries) override
		{
			auto summary = summaries->add_summaries();

			mMaxNumVertices = 20;
			mMaxNumTetrahedra = 10;
			mToughness = 100;
			mTimestep = 0.03f;
			std::vector<float> positions;
			std::vector<int> indices;
			GetTwoTetrahedraConfiguration(positions, indices, 1.0f);
			positions.resize(3 * (mMaxNumVertices + 10));
			indices.resize(4 * (mMaxNumTetrahedra + 10));

			Initialize(positions.data(), 5, mMaxNumVertices, indices.data(), 2, mMaxNumTetrahedra, mLambda, mPsi, mMu, mPhi, mToughness, mDensity);
			Deformation::TetraGroup* group = (Deformation::TetraGroup*)mData;

			for (auto& vert : group->mVertices)
				vert.mVelocity += glm::vec3(0, -16, 0);

			group->mSaveEveryXSteps = 1;

			float maxEigenvalue = -1;
			float maxEigenvalueTime = 0.0f;

			SaveFrame(summary, *group);

			for (int i = 0; i < 200; i++)
			{
				if (!Update(*group))
					break;

				SaveFrame(summary, *group);

				for (const auto& vertex : group->mVertices)
				{
					if (vertex.mLargestEigenvalue > maxEigenvalue)
					{
						maxEigenvalue = vertex.mLargestEigenvalue;
						maxEigenvalueTime = i * mTimestep;
					}
				}
			}

			SaveFrame(summary, *group);

			std::cout << "Max Eigenvalue : " << maxEigenvalue << " time " << maxEigenvalueTime << "s." << std::endl;
		}
	};

	class UpdateCaseBowl : public TestCase
	{
		using TestCase::TestCase;

	public:
		void Run(IronGames::SimulationSummaries* summaries) override
		{
			auto summary = summaries->add_summaries();

			std::ifstream file;
			//file.open("D:/UnityProjects/3D_Template/Assets/Resources/Bowl.obj");
			file.open("D:/UnityProjects/3D_Template/Assets/Resources/BasicMace.obj");

			std::vector<float> positions;
			std::vector<int> indices;

const float positionScale = 10.0f;
			const float heightAboveZero = 0.0001f;

			glm::vec3 boundingBoxMin(0);
			glm::vec3 boundingBoxMax(0);

			std::string line;
			float x, y, z;
			int i0, i1, i2, i3;
			char c;
			bool firstVert = false;
			while (std::getline(file, line))
			{
				std::istringstream iss(line);

				iss >> c;

				if (c == 'v')
				{
					iss >> x >> y >> z;

					glm::vec3 pos(x * positionScale, y * positionScale, z * positionScale);

					positions.push_back(pos[0]);
					positions.push_back(pos[1]);
					positions.push_back(pos[2]);

					if (firstVert)
					{
						firstVert = false;
						boundingBoxMin = pos;
						boundingBoxMax = pos;
					}
					else
					{
						for (int i = 0; i < 3; i++)
						{
							boundingBoxMin[i] = std::min(boundingBoxMin[i], pos[i]);
							boundingBoxMax[i] = std::max(boundingBoxMax[i], pos[i]);
						}
					}
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

			std::cout << "Box has dimensions : (" << boundingBoxMax.x - boundingBoxMin.x << ", " << boundingBoxMax.y - boundingBoxMin.y << ", " << boundingBoxMax.z - boundingBoxMin.z << ")" << std::endl;

			for (int i = 0; i < numVertices; i++)
				positions[3 * i + 1] = positions[3 * i + 1] - boundingBoxMin[1] + heightAboveZero;

			mTimestep = 1e-4;
			//mTimestep = 1e-3;
			mLambda = 0.0f;
			mMu = 5.29e7f;
			mPhi = 0.0f;
			mPsi = 198.0f;
			mDensity = 1013.0f;
			mToughness = 73.6;
			Initialize(positions.data(), numVertices, maxNumVertices, indices.data(), numTetrahedra, maxNumTetrahedra, mLambda, mPsi, mMu, mPhi, mToughness, mDensity);
			Deformation::TetraGroup* group = (Deformation::TetraGroup*)mData;

			//for (const auto& pair : group->mIdToTetrahedra)
			//{
			//	auto angle = pair.second.GetMinDihedralAngle(group->mVertices);
			//	std::cout << "Min angle : " << angle << std::endl;
			//}

			std::cout << "Min Vertex Mass : " << group->GetMinVertexMass() << std::endl;
			std::cout << "Max Vertex Mass : " << group->GetMaxVertexMass() << std::endl;

			for (auto& vert : group->mVertices)
				vert.mVelocity = glm::vec3(0, -5.4f, 0);

			float maxEigenvalue = -1;
			float maxEigenvalueTime = 0.0f;

			SaveFrame(summary, *group);

			int numSteps = 3000;
			for (int i = 0; i < numSteps; i++)
			{
				if (!Update(*group))
					break;

				if (i % (numSteps/10) == 0)
				{
					double minHeight = DBL_MAX;
					for (const auto& vert : group->mVertices)
						minHeight = std::min(vert.mPosition.y, minHeight);

					std::cout << "Progress Report at step : " << i << std::endl;
					std::cout << "Num Vertices : " << group->mVertices.size() << " out of max " << group->mMaxNumVertices << std::endl;
					std::cout << "Num Tetrahedra : " << group->mIdToTetrahedra.size() << " out of max " << group->mMaxNumTetrahedra << std::endl;
					std::cout << "Average max eigenvalue across all vertices : " << group->GetAverageMaxEigenvalue() << std::endl;
					std::cout << "Max eigenvalue across all vertices : " << group->GetMaxEigenvalue() << std::endl;
					std::cout << "Total Mass : " << group->GetTotalMass() << std::endl;
					std::cout << "Total Volume : " << group->GetTotalVolume() << std::endl;
					std::cout << "Min Vertex Height : " << minHeight << std::endl;
					SaveFrame(summary, *group);
				}

				if (i % 1000 == 0)
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
		void Run(IronGames::SimulationSummaries* summaries) override
		{
			auto summary = summaries->add_summaries();

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
		void Run(IronGames::SimulationSummaries* summaries) override
		{
			auto summary = summaries->add_summaries();

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
		void Run(IronGames::SimulationSummaries* summaries) override
		{
			auto summary = summaries->add_summaries();

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
		void Run(IronGames::SimulationSummaries* summaries) override
		{
			auto summary = summaries->add_summaries();

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
		void Run(IronGames::SimulationSummaries* summaries) override
		{
			auto summary = summaries->add_summaries();

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
		void Run(IronGames::SimulationSummaries* summaries) override
		{
			auto summary = summaries->add_summaries();

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
		void Run(IronGames::SimulationSummaries* summaries) override
		{
			auto summary = summaries->add_summaries();

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

	class SimpleCaseSingleTetDeformed : public TestCase
	{
		using TestCase::TestCase;
	public:
		void Run(IronGames::SimulationSummaries* summaries) override
		{
			float a = 2.0f;
			int numSims = 20;
			float delta = a / 2 / numSims;

			for (int sim = 0; sim <= numSims; sim++)
			{
				auto summary = summaries->add_summaries();

				std::vector<float> positions = { -a / 2, 0, 0,
												a / 2, 0, 0,
												0, 0, sqrt(3.0f) / 2.0f * a,
												0, a / 2, sqrt(3.0f) / 4.0f * a };
				std::vector<int> indices = { 0,1,2,3 };

				Initialize(positions.data(), 4, mMaxNumVertices, indices.data(), 1, mMaxNumTetrahedra, mLambda, mPsi, mMu, mPhi, mToughness, mDensity);
				Deformation::TetraGroup* group = (Deformation::TetraGroup*)mData;

				mTimestep = 0.01f;
				//group->mVertices[3].mPosition -= glm::dvec3(0, 0.1f, 0);
				group->mVertices[3].mPosition -= glm::dvec3(0, sim *delta, 0);

				std::cout << "initial deformation amount : " << sim * delta << std::endl;

				group->mSaveEveryXSteps = 1;

				SaveFrame(summary, *group);

				auto numSteps = 250;
				int i = 0;
				for (; i < numSteps; i++)
				{
					if (!Update(*group))
						break;

					if (i % (numSteps / numSteps) == 0)
						SaveFrame(summary, *group);
				}

				std::cout << "Sim " << sim << " Finished " << i << " steps out of " << numSteps << "." << std::endl << std::endl;

				SaveFrame(summary, *group);
			}
		}
	};

	class UpdateCaseTwoTetVelocityBowTie : public TestCase
	{
		using TestCase::TestCase;
	public:
		void Run(IronGames::SimulationSummaries* summaries) override
		{
			auto summary = summaries->add_summaries();

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
		void Run(IronGames::SimulationSummaries* summaries) override
		{
			auto summary = summaries->add_summaries();

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
			mTimestep = 0.000001f;
			float maxEigenvalue = -1;
			float maxEigenvalueTime = 0.0f;
			SaveFrame(summary, *group);

			for (int i = 0; i < 1000; i++)
			{
				if (!Update(*group))
					break;

				//for (const auto& vertex : group->mVertices)
				//{
				//	if (vertex.mLargestEigenvalue > maxEigenvalue)
				//	{
				//		maxEigenvalue = vertex.mLargestEigenvalue;
				//		maxEigenvalueTime = i * mTimestep;
				//	}
				//}

				SaveFrame(summary, *group);
			}

			SaveFrame(summary, *group);

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
	bool breakOnFailure = false;

	// Search
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
	mTestCases.push_back(new SimpleCaseSingleTetDeformed(breakOnFailure));

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
	mTestCases[testNum]->Run(mSummaries);
}

////////////////////////////////////////////////////////////////////////////////
