#include "Deformation.hpp"
#include "DeformationAPI.h"
#include "TestFractureManager.h"
#include "MeshVolume.h"
#include <glm/vec4.hpp>
#include <glm/common.hpp>
#include <iostream>
#include <fstream>

////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char* argv[]) 
{
	//{
	//	IronGames::SimulationSummaries summaries;

	//	Deformation::TestFractureManager testManager(&summaries);
	//	testManager.RunAllTestCases();

	//	// write to file...
	//	std::ofstream ofs("simulation.summary", std::ios_base::out | std::ios_base::binary);
	//	summaries.SerializeToOstream(&ofs);
	//}

	//extern "C" int __declspec(dllexport) __stdcall CreateSurfaceMesh(
	//float* cubeCenters, int numCubes, float* tetrahedraPositions, int numTetrahedra,
	//float* outVertexPositions, int* outNumVertices, int numMaxVertices,
	//int* outTriangleIndices, int* outNumTriangles, int numMaxTriangles)
	std::vector<float> cubeCenters = { 0,0,0, 1,0,0 };
	int numCubes = 2;
	int numTet = 0;
	int maxNumVertices = 1000;
	int maxNumTriangles = 1000;
	std::vector<float> outVertexPositions(maxNumVertices);
	std::vector<int> outTriangleIndices(3* maxNumTriangles);
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

	vol.Contains(glm::vec3(0.0f, -0.5f/3.0f, -0.5f/3.0f));

	return 0;
}