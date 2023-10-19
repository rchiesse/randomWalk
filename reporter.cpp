#include "reporter.h"

namespace sim{

rtl				 Reporter::timeElapsed;
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
	timeElapsed = (rtl)clock() - timeScope[timeScope.size() - 1];
	timeScope.resize(timeScope.size() - 1);

	//User-friendly "elapsed-time" message ("hh:mm:ss" format):
	rtl s = timeElapsed / CLOCKS_PER_SEC;	//seconds
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

void Reporter::durationInfo(const rtl& duration) {
	std::cout << std::fixed << std::setprecision(2);
	std::cout << " Sim. time: " << duration;
}
void Reporter::avSimTimeInfo(const rtl& avDuration) {
	std::cout << "\nAverage simulation time = " << avDuration;
}
void Reporter::networkInfo(const uint& n, const uint& m, const rtl& averageDegree, const uint& largestDegree, const uint& smallestDegree, const uint& lccSize) {
	using std::endl;
	std::cout << endl
		<< "\t ---> Network name: "		<< NWTK_LABEL		<< '\n'
		<< "\t ---> "					<< n << " nodes."	<< '\n'
		<< "\t ---> "					<< m << " edges."	<< '\n'
		<< "\t ---> Average degree: "	<< averageDegree	<< '\n'
		<< "\t ---> Largest degree: "	<< largestDegree	<< '\n'
		<< "\t ---> Smallest degree: "	<< smallestDegree	<< '\n'
		<< "\t ---> LCC size: "			<< lccSize << " (" << ((rtl)lccSize / n) * 100 << "%)";
}
void Reporter::simulationInfo(const uint& itotal, const rtl& ROUNDS, const rtl& T, const rtl& NUM_AGENTS, const rtl& TAU_aa, const rtl& GAMMA_a, const rtl& LAMBDA, const rtl& Wi, const rtl& Ws, const rtl& avDegree, const rtl& _2ndMoment, const rtl& Eag, const rtl& bk, const rtl& maxKbnb, const std::vector<rtl>& block_prob) {
	rtl w = (Wi + Ws) / 2.0;
	rtl sigma = TAU_aa / (2.0 * LAMBDA + TAU_aa);

	rtl elambda = LAMBDA * (1.0 - (1.0 / avDegree));
	rtl esigma = TAU_aa / (2.0 * elambda + TAU_aa);

	//rtl sigma = 1;

	//TESTE!!!
	rtl sum_b2pb2 = 0, sum_bpb2 = 0, factor, val;
	for (uint b = (uint)block_prob.size() - 1; b > 0; --b) {
		sum_b2pb2 += (size_t)b * b * block_prob[b] * block_prob[b];
		sum_bpb2  += (size_t)b * block_prob[b] * block_prob[b];
	}
	factor = sum_b2pb2 / sum_bpb2;
	val = (NUM_AGENTS / (N * avDegree)) * factor;

	rtl tamExpBlock = (N / avDegree) * sum_bpb2;
	rtl expPopulation = ((NUM_AGENTS * _2ndMoment) / (pow(avDegree, 3))) * sum_bpb2;
	rtl _2sl = 2.0 * sigma * LAMBDA;
	rtl beta_ronald = (2.0 * esigma * NUM_AGENTS * elambda * _2ndMoment * w) / ( N * pow(avDegree, 2.0));
	rtl exp_agglom = (NUM_AGENTS * _2ndMoment) / (N * pow(avDegree, 2.0));
	rtl beta_asym = TAU_aa * w * exp_agglom;
	rtl w_bound		= std::min((rtl)1.0, (GAMMA_a * N * pow(avDegree, 2.0)) / (2.0 * esigma * elambda * NUM_AGENTS * _2ndMoment));
	rtl w_asym_lambda	= std::min((rtl)1.0, (GAMMA_a * N * pow(avDegree, 2.0)) / (TAU_aa * NUM_AGENTS * _2ndMoment));
	std::cout << '\n'
		<< "\tROUNDS: "								<< ROUNDS		<< '\n'
		<< "\tT: "									<< T			<< '\n'
		<< "\tNUM_AGENTS: "							<< NUM_AGENTS	<< '\n'
		<< "\tInitially infected: "					<< itotal		<< '\n'
		<< "\tWi: "									<< Wi			<< '\n'
		<< "\tWs: "									<< Ws			<< '\n'
		<< "\tTAU (Infect): "						<< TAU_aa		<< '\n'
		<< "\tGAMMA (Recover): "					<< GAMMA_a		<< '\n'
		<< "\tLAMBDA (Walk): "						<< LAMBDA		<< '\n'
		<< "\tBETA (Infection force): "				<< beta_ronald << '\n'
		<< "\tBETA (Asymptotic on lambda): "		<< beta_asym << '\n'
		<< "\t2nd moment (<b^2>): "					<< _2ndMoment << '\n'
		<< "\t1st moment squared (<b>^2): "			<< avDegree * avDegree << '\n'
		<< "\t<b^2> / <b>^2: "						<< _2ndMoment / (avDegree * avDegree) << '\n'
		<< "\t<b>_K: "								<< bk << '\n'
		<< "\tE[#Ag] per Node (uniform): "			<< NUM_AGENTS / N << '\n'
		<< "\tK<b^2>/(n<b>^2): "					<< (NUM_AGENTS * _2ndMoment) / (N * avDegree * avDegree) << '\n'
		<< "\tval: "								<< val << '\n'
		<< "\ttamExpBlock: "						<< tamExpBlock << '\n'
		<< "\texpPopulation: "						<< expPopulation << '\n'
		<< "\tE[#Ag] (q_b weighted): "				<< Eag << "; E[#Ag]/<b>_K = " << Eag/bk << '\n'
		<< "\tmaxKbnb: "							<< maxKbnb << '\n'
		
		<< "\tE[#hops as I]: "						<< LAMBDA / GAMMA_a << '\n'
		<< "\tEarly mobility (50 * (1.0/LAMBDA)): "	<< (50.0 * (1.0 / LAMBDA)) << '\n'
		<< "\tw ((Wi + Ws) / 2): "					<< w			<< '\n'
		<< "\tw_bound: "							<< w_bound		<< '\n'
		<< "\tw_asym_lambda: "						<< w_asym_lambda << '\n'
		<< "\tR0 (BETA/GAMMA): "					<< beta_ronald/GAMMA_a	<< '\n'
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
	std::cout << "Initial occupancy ratio (NUM_AGENTS/n) = " << numAg << "/" << n << " = " << (rtl)numAg / n << std::endl;
}
void Reporter::minMaxOccInfo(const uint& min, const uint& max, const graph::node& whereMin, const graph::node& whereMax, const uint& ntwSize, const uint& numAgents) {
	std::cout << "\nMax occupancy over all rounds: node " << whereMax << " with " << max << " agents (" << ((rtl)max / numAgents) * 100 << "% of the agents)" << std::endl;
	std::cout   << "Min occupancy over all rounds: node " << whereMin << " with " << min << " agents (" << ((rtl)min / numAgents) * 100 << "% of the agents)" << std::endl;
}

#endif
} // Namespace sim