#pragma once
#include "core.pb.h"
#include <vector>

////////////////////////////////////////////////////////////////////////////////

namespace TestDeformation
{
	class TestCase
	{
	public:
		TestCase(bool breakOnFailure)
			: mBreakOnFailure(breakOnFailure) {}

		virtual void Run(IronGames::SimulationSummary* /*summary*/) {}

	protected:
		double mLambda = 2.65e6f;
		double mPsi = 397.f;
		double mPhi = 264.f;
		double mMu = 3.97e6f;
		double mDensity = 5013.f;
		double mTimestep = 0.0001f;
		double mToughness = 10.f;

		int mMaxNumVertices = 10;
		int mMaxNumTetrahedra = 10;

	private:
		bool mBreakOnFailure = true;
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