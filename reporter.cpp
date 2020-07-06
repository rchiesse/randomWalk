#include "reporter.h"

namespace sim{

real				 Reporter::timeElapsed;
std::time_t			 Reporter::timestamp;
std::vector<clock_t> Reporter::timeScope;

void Reporter::highlight(const std::string& message) {
	using std::cout; using std::endl;
	cout << "\n\n\n" << "**************************************************************";
	cout << message << "\n";
}
void Reporter::tell(const std::string& message) {
	std::cout << message;
}
void Reporter::logTimestamp(const std::string& message) {
	const char TIMESTAMP_COLUMNS = 26;
	char timestamp_string[TIMESTAMP_COLUMNS];
	timestamp = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	ctime_s(timestamp_string, sizeof timestamp_string, &timestamp);
	std::cout << std::endl << timestamp_string << message;
}
void Reporter::startChronometer(const std::string& message) {
	std::cout << message;
	timeScope.emplace_back(clock());
}
void Reporter::stopChronometer(const std::string& message) {
	if (timeScope.size() == 0) {
		std::cout << "Invalid 'stopChronometer()' call. No chronometer exist to be stopped (i.e. no 'startChronometer()' call remains left for being stopped) in this scope.";
		return;
	}
	timeElapsed = (real)clock() - timeScope[timeScope.size() - 1];
	timeScope.resize(timeScope.size() - 1);

	//User-friendly "elapsed-time" message (i.e. in the format "hh:mm:ss"):
	real s = timeElapsed / CLOCKS_PER_SEC;	//seconds
	uint m = (uint)s / 60;					//minutes
	uint h = 0;								//hours
	if (m) {
		s = (uint)s % (m * 60);
		h = m / 60;
	}
	if (h) m /= (h * 60);
	std::stringstream friendly;
	friendly.precision(2);
	if (h)		friendly << h << "h " << m << "m " << s << "s";
	else if (m)	friendly << m << "m " << s << "s";
	else		friendly << s << "s";
	std::cout << message << " (" << timeElapsed / CLOCKS_PER_SEC << "s == " << friendly.str() << ").";
}

void Reporter::progress(const uint& round, const uint& granularity) {
	std::cout << ((round % granularity == 0) ? std::to_string(round) : ".");
}

void Reporter::durationInfo(const real& duration) {
	std::cout << std::fixed << std::setprecision(2);
	std::cout << " Sim. time: " << duration;
}
void Reporter::avSimTimeInfo(const real& avDuration) {
	std::cout << "\nAverage simulation time = " << avDuration;
}
void Reporter::networkInfo(const uint& n, const uint& m, const real& averageDegree, const uint& largestDegree, const uint& smallestDegree, const uint& lccSize) {
	using std::endl;
	std::cout << endl
		<< "\t ---> "					<< n << " nodes."	<< endl
		<< "\t ---> "					<< m << " edges."	<< endl
		<< "\t ---> Average degree: "	<< averageDegree	<< endl
		<< "\t ---> Largest degree: "	<< largestDegree	<< endl
		<< "\t ---> Smallest degree: "	<< smallestDegree	<< endl
		<< "\t ---> LCC size: "			<< lccSize << " (" << ((real)lccSize / n) * 100 << "%)" << std::endl;
}
void Reporter::simulationInfo(const uint& ROUNDS, const uint& T, const uint& numAg, const uint& itotal, const uint& n, const real& TAU, const real& GAMMA, const real& LAMBDA, const real& Ws, const real& Wi) {
	using std::endl;
	std::cout << endl
		<< "\tROUNDS: "								<< ROUNDS		<< '\n'
		<< "\tT: "									<< T			<< '\n'
		<< "\tNUM_AGENTS: "							<< numAg		<< '\n'
		<< "\tInitially inf: "						<< itotal		<< '\n'
		<< "\tN: "									<< n			<< '\n'
		<< "\tTAU (Infect): "						<< TAU			<< '\n'
		<< "\tGAMMA (Recover): "					<< GAMMA		<< '\n'
		<< "\tLAMBDA (Walk): "						<< LAMBDA		<< '\n'
		<< "\tWs (dangerous-destiny tolerance): "	<< Ws			<< '\n'
		<< "\tWi (safe-destiny tolerance): "		<< Wi			<< '\n'
#ifdef i_t_FROM_MODEL
		<< "\tR0 (Reprod. number = BETA/GAMMA): "			<< BETA/GAMMA			<< '\n'
		<< "\tR0_pfx (Protection FX considered): "			<< Ro_pfx << '\n'
		<< "\ti_inf (Estimated % infected by the end): "	<< 1.0 - GAMMA/BETA		<< '\n'
		<< "\ti_inf_pfx (from \\beta(i_0)): "		<< i_inf_pfx << '\n'
#endif
		;
}
void Reporter::errorOpening(const std::string& fileName) {
	std::cout << "\nError while trying to open the file " << fileName << std::endl;
}
void Reporter::openedWithSuccess(const std::string& fileName) {
	std::cout << std::endl << "File " << fileName << " opened with success." << std::endl;
}
#ifdef OCCUPANCY
void Reporter::iniOccRatioInfo(const uint& n, const uint& numAg) {
	std::cout << "Initial occupancy ratio (NUM_AGENTS/n) = " << numAg << "/" << n << " = " << (real)numAg / n << std::endl;
}
void Reporter::minMaxOccInfo(const uint& min, const uint& max, const graph::node& whereMin, const graph::node& whereMax, const uint& ntwSize, const uint& numAgents) {
	std::cout << "\nMax occupancy over all rounds: node " << whereMax << " with " << max << " agents (" << ((real)max / numAgents) * 100 << "% of the agents)" << std::endl;
	std::cout   << "Min occupancy over all rounds: node " << whereMin << " with " << min << " agents (" << ((real)min / numAgents) * 100 << "% of the agents)" << std::endl;
}

#endif
} // Namespace sim