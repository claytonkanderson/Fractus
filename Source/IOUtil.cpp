//
//  Vertex.cpp
//  Fracturing
//
//  Created by Clayton Anderson on 4/17/17.
//  Copyright © 2017 Clayton Anderson. All rights reserved.
//

#include "IOUtil.hpp"

#include <fstream>
#include <string>

////////////////////////////////////////////////////////////////////////////////
// Duplicated with StressGeometry.
namespace
{
	void WriteVertex(std::ofstream & outFile, double x, double y, double z) {
		outFile << "v " << x << " " << y << " " << z << std::endl;
	}

	void WriteFace(std::ofstream & outFile, unsigned int i, unsigned int j, unsigned int k) {
		outFile << "f " << i + 1 << " " << j + 1 << " " << k + 1 << std::endl;
	}

	void WriteTetrahedron(std::ofstream & outFile, unsigned int i, unsigned int j, unsigned int k, unsigned int l) {
		outFile << "t " << i + 1 << " " << j + 1 << " " << k + 1 << " " << l + 1 << std::endl;
	}
}

void IOUtil::LoadTetrahedronObj(TetraGroup & outGroup, const char * filename)
{
	std::ifstream file(filename);
	std::string str;

	char * pch;
	glm::vec3 tempPos;

	std::vector<size_t> indices;

	bool oneBasedIndices = true;

	while (std::getline(file, str))
	{
		if (str.empty())
			continue;
		else if (str.front() == 'v')
		{
			std::vector<char> cstr(str.c_str() + 1, str.c_str() + str.size() + 1);
			pch = strtok(&cstr[0], " \t");

			tempPos.x = (float)atof(pch);
			pch = strtok(NULL, " ");
			tempPos.y = (float)atof(pch);
			pch = strtok(NULL, " ");
			tempPos.z = (float)atof(pch);

			outGroup.Vertices.emplace_back();
			outGroup.Vertices.back().setPos(tempPos);

		}
		else if (str.front() == 't')
		{
			std::vector<char> cstr(str.c_str() + 1, str.c_str() + str.size() + 1);
			pch = strtok(&cstr[0], " \t");

			indices.emplace_back(atoi(pch));
			if (oneBasedIndices)
				indices.back()--;
			
			pch = strtok(NULL, " ");
			indices.emplace_back(atoi(pch));
			if (oneBasedIndices)
				indices.back()--;
			
			pch = strtok(NULL, " ");
			indices.emplace_back(atoi(pch));
			if (oneBasedIndices)
				indices.back()--;

			pch = strtok(NULL, " ");
			indices.emplace_back(atoi(pch));
			if (oneBasedIndices)
				indices.back()--;
		}
	}

	for (unsigned int i = 0; i < indices.size(); i+=4)
	{
		outGroup.AddTetrahedron(
			new Tetrahedron(
				outGroup.GetVertex(indices[i]),
				outGroup.GetVertex(indices[i+1]),
				outGroup.GetVertex(indices[i+2]),
				outGroup.GetVertex(indices[i+3])
			)
		);
	}
	
	outGroup.Initialize();

	//for (const auto & tet : outGroup.GetTetrahedra())
	//	std::cout << "Volume : " << tet->GetVolume() << std::endl;

	std::cout << "Finished reading " << filename << std::endl;
	std::cout << "Number of Vertices: " << outGroup.Vertices.size() << std::endl;
	std::cout << "Number of Tetrahedra: " << outGroup.GetTetrahedra().size() << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

void IOUtil::WriteTetrahedronObj(const TetraGroup & group, const char * filename)
{
	std::ofstream outfile;
	outfile.open(filename);

	// Write Vertices
	for (const auto & vert : group.Vertices)
		WriteVertex(outfile, vert.getPos().x, vert.getPos().y, vert.getPos().z);

	// Write indices
	for (const auto & tet : group.GetTetrahedra())
	{
		// Slow .. should be using indices
		WriteTetrahedron(outfile,
			group.GetVertexIndex(tet->GetVertex(0)),
			group.GetVertexIndex(tet->GetVertex(1)),
			group.GetVertexIndex(tet->GetVertex(2)),
			group.GetVertexIndex(tet->GetVertex(3))
		);
	}

	outfile.close();
}

////////////////////////////////////////////////////////////////////////////////