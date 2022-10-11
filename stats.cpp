#include "stats.h"
using namespace sim;
std::ofstream Stats::infFracData;
real Stats::t			= 0;
uint Stats::numAgents	= 0;
uint Stats::rounds		= 0;
real Stats::tau			= 0;
real Stats::gamma		= 0;
real Stats::lambda		= 0;
real Stats::Wi			= 0;
real Stats::Ws			= 0;

#ifdef INFECTED_FRACTION
uint Stats::_bufferPos = 0;
std::vector<real> Stats::infectedAgFractionBuffer(BUFFER_SIZE);
std::vector<real> Stats::infectedSiteFractionBuffer(BUFFER_SIZE);
std::vector<real> Stats::timestampBuffer(BUFFER_SIZE);
std::vector<uint> Stats::agentBuffer(BUFFER_SIZE);
//std::vector<std::string> Stats::actionBuffer(BUFFER_SIZE);
#endif

#ifdef ESTIMATE_PROBS
real Stats::meetings = 0;
real Stats::siMeetings = 0;
real Stats::totalExposition = 0;
real Stats::totalFate = 0;
real Stats::totalInfections = 0;

real Stats::Psi;
real Stats::Pinf;
real Stats::meetingRate_;
real Stats::siMeetingRate;
std::ofstream Stats::probsData;
#endif


#ifdef OCCUPANCY
std::vector<real> Stats::avOcc;
std::vector<real> Stats::avSOcc;
std::vector<real> Stats::avIOcc;
std::vector<real> Stats::lastUpdate;
std::vector<real> Stats::sLastUpdate;
std::vector<real> Stats::iLastUpdate;
std::vector<uint> Stats::maxOcc;
std::vector<uint> Stats::maxSOcc;
std::vector<uint> Stats::maxIOcc;
std::vector<uint> Stats::minOcc;
std::vector<uint> Stats::minSOcc;
std::vector<uint> Stats::minIOcc;
std::ofstream Stats::occData;
std::ofstream Stats::netOccData;
#endif //OCCUPANCY
real Stats::avDur = 0;
uint Stats::partials = 0;
bool Stats::avDurComputed = false;
std::ofstream Stats::avDurData;

#ifdef ESTIMATE_PROBS
void Stats::probsToFile() {
	Psi = siMeetings / meetings;
	Pinf = totalInfections / siMeetings;
	meetingRate_ = meetings / T;
	siMeetingRate = siMeetings / T;
	probsData << Psi << '\t'
		<< Pinf << '\t'
		<< sim::SIGMA_aa << '\t' //modelPinf
		//<< meetingRate_		 								<< '\t' //meeting rate from simulation
		//<< meetingRate_ * Pinf * K							<< '\t' //meeting rate from PFX model
		//<< (sim::meetingRate * (pow(10,(uint)log10(N))))	<< '\t' //meeting rate from standard model
		//<< (meetingRate_pfx  * (pow(10,(uint)log10(N))))	<< '\t' //meeting rate from PFX model
		//<< SIGMA_aa * meetingRate_pfx * K				<< '\t' //meeting rate from simulation
		<< siMeetingRate << std::endl;
}
#endif
void Stats::setParams(const real& time, const uint& _numAgents, const uint& _rounds, const real& _tau, const real& _gamma, const real& _lambda, const real& _Wi, const real& _Ws) {
	t			= time;
	numAgents	= _numAgents;
	rounds		= _rounds;
	tau			= _tau;
	gamma		= _gamma;
	lambda		= _lambda;
	Wi			= _Wi;
	Ws			= _Ws;	
}

//void Stats::initStream(const real& Ws, const real& Wi, const streamType& s, const std::string& ntwLabel, const uint& ntwSize, const uint& numAgents, const double& tau, const double& gamma, const double& lambda, const uint& T, const uint& rounds) {
void Stats::initStream(const streamType& s) {
	setBasename();
	std::stringstream ss;
	ss.precision(2);
	std::string fileName;
	switch (s) {
#ifdef INFECTED_FRACTION
	case streamType::infFrac:
		ss.str("");			// ----> Clear content.
		ss << "fractionInfected_" << baseName;
		{
			std::string ref = ss.str();
#ifdef BYPASS_SIMULATION
			genPlotScript(ref, true);
#else
			genPlotScript(ref);
#endif
		}
		ss.str("");			// ----> Clear content.
		ss << "./stats/fractionInfected_" << baseName << ".csv";
		fileName = ss.str();
		infFracData.open(fileName);

		//Header
		infFracData << "Agent\tTime\ti_ag\ti_site\n";
		break;
#endif //INFECTED_FRACTION
#ifdef ESTIMATE_PROBS
	case streamType::probs:
		ss.str("");			// ----> Clear content.
		ss << EXE_DIR << "/stats/probs_" << ntwLabel << "_Ws" << Ws << "_Wi" << Wi << "_N" << ntwSize << "_AG" << numAgents << "_T" << tau << "_G" << gamma << "_L" << lambda << "_STime" << T << "_R" << rounds << ".csv";
		fileName = ss.str();
		probsData.open(fileName);
		//Header
		probsData << "Psi\tPinf\tmodelPinf\tmeetingRate\tstdModelMR\tpfxModelMR\tsiMeetingRate\n";
		break;
#endif 
#ifdef OCCUPANCY
	case streamType::occupancy:
		ss.str("");			// ----> Clear content.
		ss << EXE_DIR << "/stats/nodeOccupancy_" << ntwLabel << "_N" << ntwSize << "_T" << tau << "_G" << gamma << "_L" << lambda << "_STime" << T << "_R" << rounds << ".csv";
		fileName = ss.str();
		occData.open(fileName, std::ios::app);
		//Header
		if (isEmpty(occData))
			occData << "numAg\tWs\tWi\tnode\tavOcc\tavIOcc\tavSOcc\tmaxOcc\tmaxIOcc\tmaxSOcc\tminOcc\tminIOcc\tminSOcc" << std::endl;

		ss.str("");			// ----> Clear content.
		ss << EXE_DIR << "/stats/netInfAgMean_" << ntwLabel << "_N" << ntwSize << "_L" << tau << "_G" << gamma << "_T" << lambda << "_STime" << T << "_R" << rounds << ".csv";
		fileName = ss.str();
		netOccData.open(fileName, std::ios::app);
		//Header
		if (isEmpty(netOccData))
			netOccData << "numAg\tWs\tWi\tinfAgMean\tsusAgMean" << std::endl;

		break;
#endif //OCCUPANCY
	case streamType::avDuration:
		ss.str("");			// ----> Clear content.
		ss << "./stats/averageDuration_" << baseName << ".csv";
		fileName = ss.str();
		avDurData.open(fileName, std::ios::app);
		//Header
		if (isEmpty(avDurData))
			avDurData << "numAgents\tWs\tWi\tavDuration\n";
	default:;
	}
}
void Stats::endStream(const streamType& s) {
	switch (s) {
	case streamType::avDuration:
		avDurData.close();
		break;
#ifdef ESTIMATE_PROBS
	case streamType::probs:
		probsData.close();
		break;
#endif
#ifdef OCCUPANCY
	case streamType::occupancy:
		occData.close();
		netOccData.close();
		break;
#endif //OCCUPANCY
	default:
		infFracData.close();
	}
}
#ifdef INFECTED_FRACTION
void Stats::bufferizeIFrac(const int& ag, const real& now, const std::string& action, const uint& itotal, const uint& iltotal, const uint& NUM_AGENTS, const uint& EVT_GRANULARITY) {
	infectedAgFractionBuffer[_bufferPos] = (real)itotal / NUM_AGENTS;

#ifdef NORM_SITE_PER_AG
	infectedSiteFractionBuffer[_bufferPos] = (real)iltotal / NUM_AGENTS;
#else
	infectedSiteFractionBuffer[_bufferPos] = (real)iltotal / N;
#endif
	timestampBuffer[_bufferPos] = now;
	agentBuffer[_bufferPos] = ag;
	++_bufferPos;
	if (_bufferPos == BUFFER_SIZE) {	// ----> Flush
		for (uint i = 0; i < BUFFER_SIZE; ++i) {
			//infFracData << agentBuffer[i] << "\t" << actionBuffer[i] << "\t" << timestampBuffer[i] << "\t" << infectedAgFractionBuffer[i] << "\t" << infectedSiteFractionBuffer[i] /*<< "\t" << sim::i_t(timestampBuffer[i])*/ << '\n';
			infFracData << agentBuffer[i] << "\t" << timestampBuffer[i] << "\t" << infectedAgFractionBuffer[i] << "\t" << infectedSiteFractionBuffer[i] /*<< "\t" << sim::i_t(timestampBuffer[i])*/ << '\n';
			i += EVT_GRANULARITY;
		}
		_bufferPos = 0;
	}
}
void Stats::iFracToFile(const uint& overlook) {
	uint pos = 0;
	uint remainder = _bufferPos % overlook;
	if (remainder == _bufferPos) {
		//The entire buffer must be printed:
		while (pos < _bufferPos) {
			infFracData << agentBuffer[pos] << "\t" << timestampBuffer[pos] << "\t" << infectedAgFractionBuffer[pos] << "\t" << infectedSiteFractionBuffer[pos] << '\n';
			++pos;
		}
		return;
	}

	//For each 'overlook' number of events, only one gets printed:
	//pos = 0;
	for (uint i = 0; i < _bufferPos; ++i) {
		infFracData << agentBuffer[i] << "\t" << timestampBuffer[i] << "\t" << infectedAgFractionBuffer[i] << "\t" << infectedSiteFractionBuffer[i] << '\n';
		i += overlook;
		//pos += overlook;
		//++pos;
	}

	//Remainder:
	//for (uint i = _bufferPos - remainder - 1; i < _bufferPos; ++i) {
	//for (uint i = _bufferPos; i < agentBuffer.size(); ++i) {
	//	//infFracData << agentBuffer[i] << "\t" << actionBuffer[i] << "\t" << timestampBuffer[i] << "\t" << infectedFractionBuffer[i] << "\t" << sim::i_t(timestampBuffer[i]) << "\t" << sim::i_t_2ndMmt_naive(timestampBuffer[i]) << "\t" << sim::i_t_2ndMmt(timestampBuffer[i]) << '\n';
	//	
	//	//infFracData << agentBuffer[i] << "\t" << timestampBuffer[i] << "\t" << infectedAgFractionBuffer[i] << "\t" << infectedSiteFractionBuffer[i] << '\n';
	//	infFracData << agentBuffer[pos] << "\t" << timestampBuffer[pos] << "\t" << infectedAgFractionBuffer[pos] << "\t" << infectedSiteFractionBuffer[pos] << '\n';
	//	//++pos;
	//	//i += overlook;
	//}
}
void Stats::genPlotScript(const std::string& referenceFile, const bool&& numericOnly) {
	std::string fileName(std::string("./genplot.py"));
	std::ofstream of;
	of.open(fileName);
	of <<
		"import matplotlib\n" <<
		"#Prevents the plot from being shown in the screen when saving it to file\n" <<
		"matplotlib.use('Agg')\n\n" <<

		"import matplotlib.pyplot as plt\n" <<
		"import numpy as np\n" <<
		"import csv\n\n";

	if (!numericOnly) {
		of <<
			"#Import CSV data\n" <<
			"with open(\"./stats/" << referenceFile << ".csv\", \"r\") as i :\n" <<
			"\trawdata = list(csv.reader(i, delimiter = \"\\t\"))\n\n" <<
			"myData = np.array(rawdata[1:], dtype = np.float64)\n" <<
			"timeData = myData[:, 1]\n" <<
			"infAgSimul = myData[:, 2]\n" <<
			"#infSiteSimul = myData[:, 3]\n" <<
			"cumSumAg = np.cumsum(infAgSimul)\n" <<
			"cumMeanAg = cumSumAg / np.arange(1, len(timeData)+1)\n\n" <<
			"#cumSumSites = np.cumsum(infSiteSimul)\n" <<
			"#cumMeanSites = cumSumSites / np.arange(1, len(timeData)+1)\n\n";
	}
	of <<
		"with open(\"./stats/Runge-Kutta_" << baseName << ".csv\", \"r\") as j :\n" <<
		"\trawRK = list(csv.reader(j, delimiter = \"\\t\"))\n\n" <<

		"rkData = np.array(rawRK[1:], dtype = np.float64)\n" <<
		"timeRK = rkData[:, 0]\n" <<
		"infAgRK = rkData[:, 1]\n" <<
		"#infSiteRK = rkData[:, 2]\n" <<

		"#Plot\n" <<
		"plt.figure(1, dpi = 120)\n" <<
		"plt.title(\"Fraction of Infected Agents over Time\")\n" <<
		"plt.xlabel(\"Time\")\n" <<
		"plt.ylabel(\"Infected Fraction\")\n" <<
		"plt.xlim(0, " << t << ")\n" <<
		"plt.ylim(0, 1)\n";

	if (!numericOnly) {
		of <<
			"plt.plot(timeData, infAgSimul, label = \"Simulation\")\n" <<
			"#plt.plot(timeData, infSiteSimul, label = \"InfSites\")\n" <<
			"plt.plot(timeData, cumMeanAg, label = \"Cumul. Average\")\n" <<
			"#plt.plot(timeData, cumMeanSites, label = \"Cum.Av.#infSites\")\n" <<
			"#plt.plot(timeData, infSiteSimul, label = \"Model\")\n" <<
			"#plt.plot(timeData, [np.mean(infAgSimul) for i in range(len(timeData))], label = \"Av.#infAg\")\n";
	}

	of <<
		"plt.plot(timeRK, infAgRK, label = \"Model\")\n" <<
		"#plt.plot(timeRK, infSiteRK, label = \"Model-Site\")\n" <<
		"plt.legend()\n" <<
		"plt.grid()\n" <<
		"#plt.xscale(\"log\")\n" <<
		"#plt.yscale(\"log\")\n\n" <<

		"plt.savefig(" << "\"./plots/" << referenceFile << ".pdf\"" << ")\n" <<
		"#plt.show()\n";
	of.close();
}
#endif //INFECTED_FRACTION
void Stats::writeToFile(const streamType& s, const real& Ws, const real& Wi, const uint& numAg) {
	switch (s) {
	case streamType::avDuration:
		avDurData << numAg << "\t" << Ws << "\t" << Wi << "\t" << avDuration() << '\n';
		break;
#ifdef OCCUPANCY
	case streamType::occupancy:
		const uint n = (uint)avOcc.size();
		for (graph::node v = 0; v < n; ++v) {
			occData
				<< numAg << "\t"
				<< Ws << "\t"
				<< Wi << "\t"
				<< v << "\t\t"
				<< avOcc[v] << "\t"
				<< avIOcc[v] << "\t"
				<< avSOcc[v] << "\t"
				<< maxOcc[v] << "\t\t"
				<< maxIOcc[v] << "\t\t"
				<< maxSOcc[v] << "\t\t"
				<< minOcc[v] << "\t\t"
				<< minIOcc[v] << "\t\t"
				<< minSOcc[v]
				<< std::endl;
		}

		real avInfAgents = 0;
		real avSusAgents = 0;
		for (graph::node v = 0; v < n; ++v) avInfAgents += avIOcc[v];
		for (graph::node v = 0; v < n; ++v) avSusAgents += avSOcc[v];
		netOccData
			<< numAg << "\t"
			<< Ws << "\t"
			<< Wi << "\t"
			<< avInfAgents << "\t"
			<< avSusAgents
			<< std::endl;
		break;
#endif //OCCUPANCY
	}
}
void Stats::partialsAvDur(const real& duration) {
	avDur += duration;
	++partials;
}
void Stats::resetAvDur() {
	avDur = 0;
	partials = 0;
	avDurComputed = false;
}
const real& Stats::avDuration() {
	if (avDurComputed) 
		return avDur;
	avDur = avDur / partials;
	avDurComputed = true;
	return avDur;
}
bool Stats::isEmpty(std::ofstream& s) {
	s.seekp(0, std::ios::end);
	return s.tellp() == 0;
}
void Stats::setBasename() {
	std::stringstream name;
	name << SHORT_LABEL 
		<< "_N"		<< N 
		<< "_AG"	<< numAgents 
		<< "_T"		<< tau 
		<< "_G"		<< gamma 
		<< "_L"		<< lambda 
		<< "_Wi"	<< Wi 
		<< "_Ws"	<< Ws 
		<< "_STime" << t 
		<< "_R"		<< rounds;
		//<< "_Tal"	<< TAU_al 
		//<< "_Tla"	<< TAU_la 
		//<< "_Gl"	<< GAMMA_l 
	baseName = name.str();
}
void Stats::resetStats() {
#ifdef OCCUPANCY
	for (uint i = 0; i < lastUpdate.size(); ++i) lastUpdate[i] = 0;
	for (uint i = 0; i < sLastUpdate.size(); ++i) sLastUpdate[i] = 0;
	for (uint i = 0; i < iLastUpdate.size(); ++i) iLastUpdate[i] = 0;
#endif
}

#ifdef OCCUPANCY
void Stats::initOcc(const uint& n) {
	if (avOcc.empty()) {
		avOcc.resize(n, 0);
		avSOcc.resize(n, 0);
		avIOcc.resize(n, 0);
		lastUpdate.resize(n, 0);
		sLastUpdate.resize(n, 0);
		iLastUpdate.resize(n, 0);
		maxOcc.resize(n, 0);
		maxSOcc.resize(n, 0);
		maxIOcc.resize(n, 0);
		minOcc.resize(n, n);
		minSOcc.resize(n, n);
		minIOcc.resize(n, n);
	}
}
void Stats::updateInfOccAv(const graph::node& v, const uint& prevNumA, const uint& prevNumI, const real& now) {
	const double timeFrame = now - lastUpdate[v];
	const double iTimeFrame = now - iLastUpdate[v];
	avOcc[v] += (timeFrame * prevNumA);
	avIOcc[v] += (iTimeFrame * prevNumI);
	lastUpdate[v] = now;
	iLastUpdate[v] = now;
}
void Stats::updateSusOccAv(const graph::node& v, const uint& prevNumA, const uint& prevNumS, const real& now) {
	const double timeFrame = now - lastUpdate[v];
	const double sTimeFrame = now - sLastUpdate[v];
	avOcc[v] += (timeFrame * prevNumA);
	avSOcc[v] += (sTimeFrame * prevNumS);
	lastUpdate[v] = now;
	sLastUpdate[v] = now;
}
void Stats::checkInfMaxOcc(const graph::node& v, const uint& numA, const uint& numI) {
	if (numA > maxOcc[v])  maxOcc[v] = numA;
	if (numI > maxIOcc[v]) maxIOcc[v] = numI;
}
void Stats::checkInfMinOcc(const graph::node& v, const uint& numA, const uint& numI) {
	if (numA < minOcc[v])  minOcc[v] = numA;
	if (numI < minIOcc[v]) minIOcc[v] = numI;
}
void Stats::checkSusMaxOcc(const graph::node& v, const uint& numA, const uint& numS) {
	if (numA > maxOcc[v])  maxOcc[v] = numA;
	if (numS > maxSOcc[v]) maxSOcc[v] = numS;
}
void Stats::checkSusMinOcc(const graph::node& v, const uint& numA, const uint& numS) {
	if (numA < minOcc[v])  minOcc[v] = numA;
	if (numS < minSOcc[v]) minSOcc[v] = numS;
}
void Stats::updateInfOccRaised(const graph::node& v, const uint& numA, const uint& numI, const real& now) {
	updateInfOccAv(v, numA - 1, numI - 1, now);
	checkInfMaxOcc(v, numA, numI);
}
void Stats::updateInfOccLowered(const graph::node& v, const uint& numA, const uint& numI, const real& now) {
	updateInfOccAv(v, numA + 1, numI + 1, now);
	checkInfMinOcc(v, numA, numI);
}
void Stats::updateSusOccRaised(const graph::node& v, const uint& numA, const uint& numS, const real& now) {
	updateSusOccAv(v, numA - 1, numS - 1, now);
	checkSusMaxOcc(v, numA, numS);
}
void Stats::updateSusOccLowered(const graph::node& v, const uint& numA, const uint& numS, const real& now) {
	updateSusOccAv(v, numA + 1, numS + 1, now);
	checkSusMinOcc(v, numA, numS);
}
void Stats::getGlobalMaxOcc(uint& _max, graph::node& _where) {
	for (uint i = 0; i < maxOcc.size(); ++i) {
		if (maxOcc[i] > _max) { _max = maxOcc[i]; _where = i; }
	}
}
void Stats::getGlobalMinOcc(uint& _min, graph::node& _where) {
	for (uint i = 0; i < minOcc.size(); ++i) {
		if (minOcc[i] < _min) { _min = minOcc[i]; _where = i; }
	}
}
void Stats::computeOcc(const real& totalTime) {
	const uint n = (uint)avOcc.size();
	for (graph::node v = 0; v < n; ++v) avOcc[v] /= totalTime;
	for (graph::node v = 0; v < n; ++v) avIOcc[v] /= totalTime;
	for (graph::node v = 0; v < n; ++v) avSOcc[v] /= totalTime;
}
//void Stats::netInfAgMeanToFile(const real& Ws, const uint& numAg) {
//	const uint n = (uint)avOcc.size();
//	real avInfAgents = 0;
//	for (graph::node v = 0; v < n; ++v) avInfAgents += avIOcc[v];
//	netOccData
//		<< numAg << "\t"
//		<< Ws << "\t"
//		<< avInfAgents
//		<< std::endl;
//}
//void Stats::occToFile(const real& Ws, const uint& numAg) {
//	const uint n = (uint)avOcc.size();
//	for (graph::node v = 0; v < n; ++v){
//		occData
//			<< numAg << "\t"
//			<< Ws << "\t"
//			<< v << "\t\t"
//			<< avOcc  [v] << "\t"
//			<< avIOcc [v] << "\t"
//			<< avSOcc [v] << "\t"
//			<< maxOcc [v] << "\t\t"
//			<< maxIOcc[v] << "\t\t"
//			<< maxSOcc[v] << "\t\t"
//			<< minOcc [v] << "\t\t"
//			<< minIOcc[v] << "\t\t"
//			<< minSOcc[v]
//			<< std::endl;
//	}
//}

#endif //OCCUPANCY