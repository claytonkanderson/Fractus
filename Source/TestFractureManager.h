#pragma once
#include "core.pb.h"
#include "Deformation.hpp"
#include <vector>

////////////////////////////////////////////////////////////////////////////////

namespace Deformation
{
	class TestCase
	{
	public:
		TestCase(bool breakOnFailure)
			: mBreakOnFailure(breakOnFailure) {}

		virtual void Run(IronGames::SimulationSummary* /*summary*/) {}
		
		bool Update(Deformation::TetraGroup& group);

	protected:
		float mLambda = 2.65e6f;
		float mPsi = 397.f;
		float mPhi = 264.f;
		float mMu = 3.97e6f;
		float mDensity = 1013.f;
		float mTimestep = 0.0001f;
		float mToughness = 10.f;

		int mMaxNumVertices = 10;
		int mMaxNumTetrahedra = 10;

	private:
		bool mBreakOnFailure = false;
	};

	class TestFractureManager
	{
	public:
		TestFractureManager(IronGames::SimulationSummaries *summaries = nullptr);

		~TestFractureManager();

		void RunAllTestCases();

		void RunTestCase(int testNum);

		int GetIndexOfLastCase() const { return mTestCases.size() - 1; }

	private:
		IronGames::SimulationSummaries* mSummaries = nullptr;
		std::vector<TestCase*> mTestCases;
	};

}