#pragma once
#include "defs.h"
#include <iomanip>		// ----> std::setprecision()

namespace sim {		// ----> Simulator's namespace.
class Reporter {
private:
	static real timeElapsed;
	static std::time_t timestamp;
	static std::vector<clock_t> timeScope;			// ----> Allows chronometers to be nested: a chronometer does not need to stop for another one to start.

public:
	static void highlight		 (const std::string& message);
	static void tell			 (const std::string& message);
	static void logTimestamp	 (const std::string& message);
	static void startChronometer (const std::string& message = "");
	static void stopChronometer	 (const std::string& message);
	static void durationInfo	 (const real& duration);
	static void progress		 (const uint& round, const uint& granularity = 10);
	static void avSimTimeInfo	 (const real& avDuration);
	static void networkInfo		 (const uint& n, const uint& m, const real& averageDegree, const uint& largestDegree, const uint& smallestDegree, const uint& lccSize);
	static void simulationInfo	 (const uint& ROUNDS, const uint& T, const uint& numAg, const uint& itotal, const uint& n, const real& TAU, const real& GAMMA, const real& LAMBDA, const real& Ws, const real& Wi);
	static void errorOpening	 (const std::string& fileName);
	static void openedWithSuccess(const std::string& fileName);
#ifdef OCCUPANCY
	static void iniOccRatioInfo	 (const uint& n, const uint& numAg);
	static void minMaxOccInfo	 (const uint& min, const uint& max, const graph::node& whereMin, const graph::node& whereMax, const uint& ntwSize, const uint& numAgents);
#endif
};
}