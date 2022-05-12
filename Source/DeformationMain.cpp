#include "Deformation.hpp"
#include "DeformationAPI.h"
#include "TestFractureManager.h"
#include "MeshVolume.h"
#include "FractureUtil.h"
#include "ConvexIntersection.h"

#include <glm/vec4.hpp>
#include <glm/common.hpp>
#include <iostream>
#include <fstream>

namespace
{
	void RunTestCases()
	{
		IronGames::SimulationSummaries summaries;

		Deformation::TestFractureManager testManager(&summaries);
		//testManager.RunAllTestCases();
		//testManager.RunTestCase(testManager.GetIndexOfLastCase());
		testManager.RunTestCase(1);

		// write to file...
		std::ofstream ofs("D:/UnityProjects/3D_Template/Assets/Resources/simulation.summary", std::ios_base::out | std::ios_base::binary);
		summaries.SerializeToOstream(&ofs);
	}

	void TestSurfaceMeshCreation()
	{		
		std::vector<float> cubeCenters = { 0,0,0, 1,0,0 };
		int numCubes = 2;
		int numTet = 0;
		int maxNumVertices = 1000;
		int maxNumTriangles = 1000;
		std::vector<float> outVertexPositions(maxNumVertices);
		std::vector<int> outTriangleIndices(3 * maxNumTriangles);
		int outNumVertices;
		int outNumTriangles;
		auto errVal = CreateSurfaceMesh(cubeCenters.data(), numCubes, nullptr, numTet, outVertexPositions.data(), &outNumVertices, maxNumVertices, outTriangleIndices.data(), &outNumTriangles, maxNumTriangles);

		Deformation::MeshVolume vol;
		vol.mCubeCenters.push_back(glm::vec3(-0.5f, 0, 0));
		Deformation::MeshVolume::Tetrahedron tet;
		tet.mVertices[0] = glm::vec3(0.0f, -0.5f, -0.5f);
		tet.mVertices[1] = glm::vec3(0.0f, 0.5f, -0.5f);
		tet.mVertices[2] = glm::vec3(0.0f, -0.5f, 0.5f);
		tet.mVertices[3] = glm::vec3(1.0f, -0.5f, -0.5f);
		vol.mTetrahedra.push_back(tet);

		vol.Contains(glm::vec3(0.0f, -0.5f / 3.0f, -0.5f / 3.0f));
	}

	void TestVolumeMesh()
	{
		std::vector<float> vertexPositions;

		for (const auto& pos : Deformation::CubeVertexPositions)
		{
			vertexPositions.push_back(pos.x - 0.5f);
			vertexPositions.push_back(pos.y - 0.5f);
			vertexPositions.push_back(pos.z - 0.5f);
		}

		int numVertices = vertexPositions.size() / 3;

		std::vector<int> triangleIndices;
		for (const auto& tri : Deformation::CubeTriangleIndices)
		{
			triangleIndices.push_back(tri[0]);
			triangleIndices.push_back(tri[1]);
			triangleIndices.push_back(tri[2]);
		}

		int numTriangles = triangleIndices.size() / 3;
		int maxNumVertices = 1000;
		int maxNumTetrahedra = 1000;
		std::vector<float> outVertexPositions(maxNumVertices);
		std::vector<int> outTriangleIndices(4 * maxNumTetrahedra);
		int outNumVertices;
		int outNumTetrahedra;
		CreateTetrahedralMesh(vertexPositions.data(), numVertices, triangleIndices.data(), numTriangles, outVertexPositions.data(),
			&outNumVertices, maxNumVertices, outTriangleIndices.data(), &outNumTetrahedra, maxNumTetrahedra);

		for (int i = 0; i < outNumVertices; i++)
			std::cout << outVertexPositions[3 * i] << ", " << outVertexPositions[3 * i + 1] << ", " << outVertexPositions[3 * i + 2] << std::endl;
	}

	void TestCollisionDetection()
	{
		std::array<glm::vec3, 4> verts
		{
			glm::vec3(-0.5f, 0.0f, 1.0f),
			glm::vec3(0.5f, 0.0f, 1.0f),
			glm::vec3(0.0f, 1.0f, 1.0f),
			glm::vec3(0.0f, 0.5f, 0.0f)
		};

		std::array<glm::vec3, 4> otherVerts
		{
			glm::vec3(-0.5f, 0.0f, -1.0f),
			glm::vec3(0.5f, 0.0f, -1.0f),
			glm::vec3(0.0f, 1.0f, -1.0f),
			glm::vec3(0.0f, 0.4f, 0.5f)
		};
		//Deformation::ConvexIntersection::CollisionResults results;
		// needs to be updated to use ConvexIntersection::ResolveCollisions
		//if (Deformation::ConvexIntersection::TetrahedraIntersection(verts, otherVerts, results))
		//	std::cout << "Collision detected!" << std::endl;
		//else
		//	std::cout << "No Collision detected!" << std::endl;

		const float cCollisionStrength = 1e7f;
		//std::cout << "Collision mag : " << cCollisionStrength * results.mVolume << std::endl;
	}
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char* argv[]) 
{
	auto start = std::chrono::high_resolution_clock::now();
	RunTestCases();
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "Finished running after " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s." << std::endl;
	return 0;
}