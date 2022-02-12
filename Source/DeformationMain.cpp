#include "TestDeformation.hpp"
#include "DeformationAPI.h"
#include "TestFractureManager.h"
#include <glm/vec4.hpp>
#include <glm/common.hpp>
#include <iostream>
#include <fstream>

#include <Mathematics/Delaunay3.h>
#include <Mathematics/IntrSegment3Plane3.h>
#include <Mathematics/AlignedBox.h>
#include <Mathematics/DistPointHyperplane.h>

////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char* argv[]) 
{
	IronGames::SimulationSummaries summaries;

	try
	{
		TestDeformation::TestFractureManager testManager(&summaries);
		testManager.RunAllTestCases();
	}
	catch (const std::exception& e)
	{
		// nah
	};


	// write to file...
	std::ofstream ofs("simulation.summary", std::ios_base::out | std::ios_base::binary);
	summaries.SerializeToOstream(&ofs);

}