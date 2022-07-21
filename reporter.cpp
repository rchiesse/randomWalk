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
		std::cout << "Invalid 'stopChronometer()' call. No chronometer exists to be stopped (i.e. no 'startChronometer()' call remains left for being stopped) in this scope.";
		return;
	}
	timeElapsed = (real)clock() - timeScope[timeScope.size() - 1];
	timeScope.resize(timeScope.size() - 1);

	//User-friendly "elapsed-time" message ("hh:mm:ss" format):
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
	if		(h)	friendly << h << "h " << m << "m " << s << "s";
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
		<< "\t ---> Network name: "		<< NWTK_LABEL		<< '\n'
		<< "\t ---> "					<< n << " nodes."	<< '\n'
		<< "\t ---> "					<< m << " edges."	<< '\n'
		<< "\t ---> Average degree: "	<< averageDegree	<< '\n'
		<< "\t ---> Largest degree: "	<< largestDegree	<< '\n'
		<< "\t ---> Smallest degree: "	<< smallestDegree	<< '\n'
		<< "\t ---> LCC size: "			<< lccSize << " (" << ((real)lccSize / n) * 100 << "%)";
}
void Reporter::simulationInfo(const uint& itotal, const real& ROUNDS, const real& T, const real& NUM_AGENTS, const real& TAU_aa, const real& GAMMA_a, const real& LAMBDA) {
	std::cout << '\n'
		<< "\tROUNDS: "								<< ROUNDS		<< '\n'
		<< "\tT: "									<< T			<< '\n'
		<< "\tNUM_AGENTS: "							<< NUM_AGENTS	<< '\n'
		<< "\tInitially infected: "					<< itotal		<< '\n'
		<< "\tN: "									<< N			<< '\n'
		<< "\tTAU (Infect): "						<< TAU_aa		<< '\n'
		<< "\tGAMMA (Recover): "					<< GAMMA_a		<< '\n'
		<< "\tLAMBDA (Walk): "						<< LAMBDA		<< '\n'
		//<< "\tTAU_al (Infect): "					<< TAU_al		<< '\n'
		//<< "\tTAU_la (Infect): "					<< TAU_la		<< '\n'
		//<< "\tGAMMA_l (Recover): "				<< GAMMA_l		<< '\n'
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