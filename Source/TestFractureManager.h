#pragma once
#include "core.pb.h"
#include <vector>

////////////////////////////////////////////////////////////////////////////////

namespace TestDeformation
{
	class TestCase
	{
	public:

		virtual void Run(IronGames::SimulationSummary* /*summary*/) {}
	};

	class TestFractureManager
	{
	public:
		TestFractureManager(IronGames::SimulationSummaries *summaries = nullptr);

		~TestFractureManager();

		void RunAllTestCases();

		void RunTestCase(int testNum);

	private:
		IronGames::SimulationSummaries* mSummaries = nullptr;
		std::vector<TestCase*> mTestCases;
	};

}