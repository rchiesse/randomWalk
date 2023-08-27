#pragma once
#include "defs.h"
#include <iomanip>		// ----> std::setprecision()
#include <time.h>
#include <chrono>
#include <ctime>

namespace sim {		// ----> Simulator's namespace.
class Reporter {
private:
	static rtl timeElapsed;
	static std::time_t timestamp;
	static std::vector<clock_t> timeScope;			// ----> Allows chronometers to be nested: a chronometer does not need to stop for another one to start.

public:
	static void highlight		 (const std::string& message);
	static void tell			 (const std::string& message);
	static void logTimestamp	 (const std::string& message);
	static void startChronometer (const std::string& message = "");
	static void stopChronometer	 (const std::string& message);
	static void durationInfo	 (const rtl& duration);
	static void progress		 (const uint& round, const uint& granularity = 10);
	static void avSimTimeInfo	 (const rtl& avDuration);
	static void networkInfo		 (const uint& n, const uint& m, const rtl& averageDegree, const uint& largestDegree, const uint& smallestDegree, const uint& lccSize);
	//static void simulationInfo	 (const uint& ROUNDS, const uint& T, const uint& numAg, const uint& itotal, const uint& n, const rtl& TAU_aa, const rtl& TAU_al, const rtl& TAU_la, const rtl& GAMMA_a, const rtl& GAMMA_l, const rtl& LAMBDA);
	static void simulationInfo(const uint& itotal, const rtl& ROUNDS, const rtl& T, const rtl& NUM_AGENTS, const rtl& TAU_aa, const rtl& GAMMA_a, const rtl& LAMBDA, const rtl& Wi, const rtl& Ws, const rtl& avDegree, const rtl& _2ndMoment, const rtl& Eag, const rtl& bk, const rtl& maxKbnb, const std::vector<rtl>& block_prob);
	static void errorOpening	 (const std::string& fileName);
	static void openedWithSuccess(const std::string& fileName);
#ifdef OCCUPANCY
	static void iniOccRatioInfo	 (const uint& n, const uint& numAg);
	static void minMaxOccInfo	 (const uint& min, const uint& max, const graph::node& whereMin, const graph::node& whereMax, const uint& ntwSize, const uint& numAgents);
#endif
};
}