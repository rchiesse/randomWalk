#include "stats.h"
#include "utils.h"
using namespace sim;
std::ofstream Stats::infFracData;
rtl Stats::t			= 0;
uint Stats::numAgents	= 0;
uint Stats::rounds		= 0;
rtl Stats::tau			= 0;
rtl Stats::gamma		= 0;
rtl Stats::lambda		= 0;
rtl Stats::Wi			= 0;
rtl Stats::Ws			= 0;
rtl Stats::W			= 0;

#ifdef INFECTED_FRACTION
uint Stats::_bufferPos = 0;
std::vector<rtl> Stats::infectedAgFractionBuffer(BUFFER_SIZE);
std::vector<rtl> Stats::infectedSiteFractionBuffer(BUFFER_SIZE);
std::vector<rtl> Stats::timestampBuffer(BUFFER_SIZE);
std::vector<uint> Stats::agentBuffer(BUFFER_SIZE);
//std::vector<std::string> Stats::actionBuffer(BUFFER_SIZE);
#endif

#ifdef ESTIMATE_PROBS
rtl Stats::meetings = 0;
rtl Stats::siMeetings = 0;
rtl Stats::totalExposition = 0;
rtl Stats::totalFate = 0;
rtl Stats::totalInfections = 0;

rtl Stats::Psi;
rtl Stats::Pinf;
rtl Stats::meetingRate_;
rtl Stats::siMeetingRate;
std::ofstream Stats::probsData;
#endif


#ifdef OCCUPANCY
std::vector<rtl> Stats::avOcc;
std::vector<rtl> Stats::avSOcc;
std::vector<rtl> Stats::avIOcc;
std::vector<rtl> Stats::lastUpdate;
std::vector<rtl> Stats::sLastUpdate;
std::vector<rtl> Stats::iLastUpdate;
std::vector<uint> Stats::maxOcc;
std::vector<uint> Stats::maxSOcc;
std::vector<uint> Stats::maxIOcc;
std::vector<uint> Stats::minOcc;
std::vector<uint> Stats::minSOcc;
std::vector<uint> Stats::minIOcc;
std::ofstream Stats::occData;
std::ofstream Stats::netOccData;
#endif //OCCUPANCY
rtl Stats::avDur = 0;
rtl Stats::lastRoundDuration = 0;
uint Stats::partials = 0;
bool Stats::avDurComputed = false;
std::ofstream Stats::avDurDataGroupL;
std::ofstream Stats::avDurDataGroupK;
std::ofstream Stats::avDurDataK;
std::ofstream Stats::avDurDataL;
std::string Stats::avDurSetNameL;
std::string Stats::avDurSetNameK;
std::string Stats::avDurKBaseName;
std::string Stats::avDurLBaseName;

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
void Stats::setParams(const rtl& time, const uint& _numAgents, const uint& _rounds, const rtl& _tau, const rtl& _gamma, const rtl& _lambda, const rtl& _Wi, const rtl& _Ws) {
	t			= time;
	numAgents	= _numAgents;
	rounds		= _rounds;
	tau			= _tau;
	gamma		= _gamma;
	lambda		= _lambda;
	Wi			= _Wi;
	Ws			= _Ws;	
	W           = (_Wi + _Ws) / 2.0;
}

//void Stats::initStream(const rtl& Ws, const rtl& Wi, const streamType& s, const std::string& ntwLabel, const uint& ntwSize, const uint& numAgents, const double& tau, const double& gamma, const double& lambda, const uint& T, const uint& rounds) {
void Stats::initStream(const streamType& s) {
	//setBasename();
	std::stringstream ss;
	ss.precision(2);
	std::string fileName;
	std::string ref;
	switch (s) {
#ifdef INFECTED_FRACTION
	case streamType::infFrac:
		ss.str("");			// ----> Clear content.
		ss << "fractionInfected_" << baseName;
		ref = ss.str();
#ifdef BYPASS_SIMULATION
		genPlotScript(ref, true);
#else
		genPlotScript(ref);
		ss.str("");			// ----> Clear content.
		ss << "./stats/fractionInfected_" << baseName << ".csv";
		fileName = ss.str();
		infFracData.open(fileName);

		//Header
		infFracData << "Agent\tTime\ti_ag\ti_site\n";
#endif
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
		ss << EXE_DIR << "/stats/nodeOccupancy_" << NWTK_LABEL << "_N" << N << "_T" << tau << "_G" << gamma << "_L" << lambda << "_STime" << t << "_R" << rounds << ".csv";
		fileName = ss.str();
		occData.open(fileName, std::ios::app);
		//Header
		if (sim::utils::isEmpty(occData))
			occData << "numAg\tWs\tWi\tnode\tavOcc\tavIOcc\tavSOcc\tmaxOcc\tmaxIOcc\tmaxSOcc\tminOcc\tminIOcc\tminSOcc" << std::endl;

		ss.str("");			// ----> Clear content.
		ss << EXE_DIR << "/stats/netInfAgMean_" << NWTK_LABEL << "_N" << N << "_L" << tau << "_G" << gamma << "_T" << lambda << "_STime" << t << "_R" << rounds << ".csv";
		fileName = ss.str();
		netOccData.open(fileName, std::ios::app);
		//Header
		if (sim::utils::isEmpty(netOccData))
			netOccData << "numAg\tWs\tWi\tinfAgMean\tsusAgMean" << std::endl;

		break;
#endif //OCCUPANCY
	case streamType::avDuration:
		ss.str("");			// ----> Clear content.
		ss.clear();
		ss << "./stats/averages/" << avDurKBaseName << ".csv";
		fileName = ss.str();
		avDurDataK.open(fileName, std::ios::app);
		//if (!avDurDataK.is_open()) {
		//	sim::Reporter::errorOpening(fileName);
		//	return;
		//}
		if (sim::utils::isEmpty(avDurDataK))
			avDurDataK << "k,duration\n";	//Header

		ss.str("");
		ss.clear();
		ss << "./stats/averages/" << avDurLBaseName << ".csv";
		fileName = ss.str();
		avDurDataL.open(fileName, std::ios::app);
		//if (!avDurDataL.is_open()) {
		//	sim::Reporter::errorOpening(fileName);
		//	return;
		//}
		if (sim::utils::isEmpty(avDurDataL))
			avDurDataL << "lambda,duration\n";
	default:;
	}
}
void Stats::endStream(const streamType& s) {
	switch (s) {
	case streamType::avDuration:
		avDurDataK.close();
		avDurDataL.close();
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
void Stats::bufferizeIFrac(const int& ag, const rtl& now, const std::string& action, const uint& itotal, const uint& iltotal, const uint& NUM_AGENTS, const uint& EVT_GRANULARITY) {
	infectedAgFractionBuffer[_bufferPos] = (rtl)itotal / NUM_AGENTS;

#ifdef NORM_SITE_PER_AG
	infectedSiteFractionBuffer[_bufferPos] = (rtl)iltotal / NUM_AGENTS;
#else
	infectedSiteFractionBuffer[_bufferPos] = (rtl)iltotal / N;
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

	for (uint i = 0; i < _bufferPos; ++i) {
		infFracData << agentBuffer[i] << "\t" << timestampBuffer[i] << "\t" << infectedAgFractionBuffer[i] << "\t" << infectedSiteFractionBuffer[i] << '\n';
		i += overlook;
	}
}
void Stats::genPlotScript(const std::string& referenceFile, const bool&& numericOnly) {
	std::string fileName(std::string("./genplot.py"));
	std::ofstream of;
	of.open(fileName);
	of <<
		"import matplotlib\n" <<
		"#Prevents the plot from being shown in the screen when saving it to file\n" <<
		"#matplotlib.use('Agg')\n\n" <<

		"import matplotlib.pyplot as plt\n" <<
		"import matplotlib as mpl\n" <<
		"import numpy as np\n" <<
		//"from IPython.display import clear_output\n" <<
		"import csv\n\n";

	of <<
		"plt.figure(1, dpi = 120)\n" <<
		"plt.title(\"Fraction of Infected Agents over Time\")\n" <<
		"plt.xlabel(\"Time\")\n" <<
		"plt.ylabel(\"Infected Fraction\")\n" <<
		"plt.xlim(0, " << t << ")\n" <<
		"plt.ylim(0, 1)\n\n";

	of << "def incluirPlot(nomeArq, ini, lbl) :\n" <<
		"\twith open(nomeArq, \"r\") as i :\n " <<
		"\t\trawdata = list(csv.reader(i, delimiter = \"\t\"))\n\n" <<

		"\tmyData = np.array(rawdata[1:], dtype = np.float64)\n" <<
		"\txData = myData[:, ini]\n" <<
		"\tyData = myData[:, ini + 1]\n\n" <<

		"\tplt.plot(xData, yData, label = lbl)\n\n\n" <<

		"#Import CSV data\n";


	//"INFECTED FRACTION X SIMULATION TIME":
	if (!numericOnly) 
		of << "incluirPlot(\"./stats/" << referenceFile << ".csv\", 1, \"w" << W << "\")\n";
	else
		of << "#incluirPlot(\"./stats/" << referenceFile << ".csv\", 1, \"w" << W << "\")\n";

	std::string method, lbl;
	method = "MASTER";
	lbl = "master";
#ifdef MASTER
	of << "incluirPlot(\"./stats/Runge-Kutta_" << method << "_" << baseName << ".csv\", 0, \"w" << W << " " << lbl << "\")\n";
#else
	of << "#incluirPlot(\"./stats/Runge-Kutta_" << method << "_" << baseName << ".csv\", 0, \"w" << W << " " << lbl << "\")\n";
#endif // MASTER
	method = "BLOCK";
	lbl = "block";
#ifdef PER_BLOCK
	of << "incluirPlot(\"./stats/Runge-Kutta_" << method << "_" << baseName << ".csv\", 0, \"w" << W << " " << lbl << "\")\n";
#else
	of << "#incluirPlot(\"./stats/Runge-Kutta_" << method << "_" << baseName << ".csv\", 0, \"w" << W << " " << lbl << "\")\n";
#endif // PER_BLOCK
	//method = "NODE";
	//lbl = "node";

	of << "\n\n" <<
		"plt.legend()\n" <<
		"plt.grid()\n" <<
		"#plt.xscale(\"log\")\n" <<
		"#plt.yscale(\"log\")\n\n" <<

		"plt.savefig(" << "\"./plots/" << referenceFile << ".pdf\"" << ")\n" <<
		"#plt.show()\n";




	//"AVERAGE DURATION X LAMBDA" & "AVERAGE DURATION X K":
	std::stringstream avLTxt;
	avLTxt << "./stats/averages/" << avDurSetNameL << ".csv";
	std::string avlName = avLTxt.str();
	std::ifstream arq__avL(avlName);
	//if (!arq__avL.is_open()) {
	//	sim::Reporter::errorOpening(avDurSetNameL);
	//	return;
	//}
	
	of << "\n\n#AVERAGE DURATION X LAMBDA & AVERAGE DURATION X K\n" <<
		"plt.pause(0.001)\n" <<
		"plt.clf()\n" <<
		//"clear_output(wait=True)\n" <<
		"plt.figure(1, dpi = 120)\n" <<
		"plt.title(\"Average Duration over Walk Rate\")\n" <<
		"plt.xlabel(\"Walk Rate\")\n" <<
		"plt.ylabel(\"Average duration\")\n" <<
		"plt.xlim(1, 10)\n" <<
		"plt.ylim(0, " << t << ")\n" << 
		"plt.xscale(\"log\")\n";

	std::string instanceName;
	do {
		getline(arq__avL, instanceName);
		if (instanceName.empty())
			break;

		vector<std::string> labels;
		utils::split(instanceName, '-', labels);
		std::string label = labels[1];
		
		//TESTE!!!
		of << "incluirPlot(\"./stats/averages/" << instanceName << ".csv\", 0, \"" << label << "\")\n";

		//of << "with open(\"./stats/averages/" << instanceName << ".csv\", \"r\") as j :\n";
		//
		//std::replace(instanceName.begin(), instanceName.end(), '.', '_');
		//std::replace(instanceName.begin(), instanceName.end(), '-', '_');
		//of <<
		//	"\traw_"<< instanceName <<" = list(csv.reader(j, delimiter = \",\"))\n\n" <<
		//
		//	instanceName << " = np.array(raw_"  << instanceName << "[1:], dtype = np.float64)\n" <<
		//	"lambda_"  << instanceName << " = " << instanceName << "[:, 0]\n" <<
		//	"duration_"<< instanceName << " = " << instanceName << "[:, 1]\n" <<
		//	"plt.plot(lambda_"  << instanceName << ", duration_" << instanceName << ", label = \"" << label << "\")\n\n";
	} while (arq__avL.good());
	arq__avL.close();

	of << "plt.legend()\n" <<
		"plt.grid()\n" <<
		"plt.savefig(" << "\"./plots/averages/" << referenceFile << ".pdf\"" << ")\n";
	
	of.close();
}

#endif //INFECTED_FRACTION
void Stats::writeToFile(const streamType& s, const rtl& Ws, const rtl& Wi, const uint& numAg) {
	switch (s) {
	case streamType::avDuration:
		//std::stringstream sstype;
		//sstype.precision(1);
		//sstype << W;
		//std::string type("N" + std::to_string(N) + "_w" + sstype.str());

		//avDurData << numAg << "\t" << Ws << "\t" << Wi << "\t" << avDuration() << '\n';
		avDurDataK << numAg  << ',' << avDuration() << '\n';
		avDurDataL << lambda << ',' << avDuration() << '\n';
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

		rtl avInfAgents = 0;
		rtl avSusAgents = 0;
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
void Stats::partialsAvDur(const rtl& duration) {
	lastRoundDuration = duration;
	avDur += duration;
	++partials;
}
void Stats::initAvDur() {
	{
		std::stringstream ss;
		ss << "./stats/averages/" << avDurSetNameK << ".csv";
		const std::string streamName = ss.str();
		if (!utils::containsLine(streamName, avDurKBaseName)) {
			avDurDataGroupK.open(streamName, std::ios::app);
			avDurDataGroupK << avDurKBaseName << '\n';
			avDurDataGroupK.close();
		}
	}
	{
		std::stringstream ss;
		ss << "./stats/averages/" << avDurSetNameL << ".csv";
		const std::string streamName = ss.str();
		if (!utils::containsLine(streamName, avDurLBaseName)) {
			avDurDataGroupL.open(streamName, std::ios::app);
			avDurDataGroupL << avDurLBaseName << '\n';
			avDurDataGroupL.close();
		}
	}
}
void Stats::resetAvDur() {
	avDur = 0;
	partials = 0;
	avDurComputed = false;
}
const rtl& Stats::avDuration() {
	if (avDurComputed) 
		return avDur;
	avDur = avDur / partials;
	avDurComputed = true;
	return avDur;
}

void Stats::setBasename() {
	{
		std::stringstream name;
		name << SHORT_LABEL
			<< "_N" << N
			<< "_AG" << numAgents
			<< "_T" << tau
			<< "_G" << gamma
			<< "_L" << lambda
			<< "_Wi" << Wi
			<< "_Ws" << Ws
			<< "_STime" << t
			<< "_R" << rounds;
		//<< "_Tal"	<< TAU_al 
		//<< "_Tla"	<< TAU_la 
		//<< "_Gl"	<< GAMMA_l 
		baseName = name.str();
	}
	{
		std::stringstream name;
		name << "group_"
			<< SHORT_LABEL
			<< "_AG" << numAgents
			<< "_T" << tau
			<< "_G" << gamma
			<< "_STime" << t
			<< "_R" << rounds;
		avDurSetNameL = name.str();
	}
	{
		std::stringstream name;
		name << "group_"
			<< SHORT_LABEL
			<< "_T" << tau
			<< "_G" << gamma
			<< "_L" << lambda
			<< "_STime" << t
			<< "_R" << rounds;
		avDurSetNameK = name.str();
	}
	{
		std::stringstream name;
		name << SHORT_LABEL
			<< "_T" << tau
			<< "_G" << gamma
			<< "_L" << lambda
			<< "_STime" << t
			<< "_R" << rounds
			<< "-N" << N
			<< "_w" << W;
		avDurKBaseName = name.str();
	}
	{
		std::stringstream name;
		name.clear();
		name.str() = "";
		name << SHORT_LABEL
			<< "_AG" << numAgents
			<< "_T" << tau
			<< "_G" << gamma
			<< "_STime" << t
			<< "_R" << rounds
			<< "-N" << N
			<< "_w" << W;
		avDurLBaseName = name.str();
	}
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
void Stats::updateInfOccAv(const graph::node& v, const uint& prevNumA, const uint& prevNumI, const rtl& now) {
	const double timeFrame = now - lastUpdate[v];
	const double iTimeFrame = now - iLastUpdate[v];
	avOcc[v] += (timeFrame * prevNumA);
	avIOcc[v] += (iTimeFrame * prevNumI);
	lastUpdate[v] = now;
	iLastUpdate[v] = now;
}
void Stats::updateSusOccAv(const graph::node& v, const uint& prevNumA, const uint& prevNumS, const rtl& now) {
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
void Stats::updateInfOccRaised(const graph::node& v, const uint& numA, const uint& numI, const rtl& now) {
	updateInfOccAv(v, numA - 1, numI - 1, now);
	checkInfMaxOcc(v, numA, numI);
}
void Stats::updateInfOccLowered(const graph::node& v, const uint& numA, const uint& numI, const rtl& now) {
	updateInfOccAv(v, numA + 1, numI + 1, now);
	checkInfMinOcc(v, numA, numI);
}
void Stats::updateSusOccRaised(const graph::node& v, const uint& numA, const uint& numS, const rtl& now) {
	updateSusOccAv(v, numA - 1, numS - 1, now);
	checkSusMaxOcc(v, numA, numS);
}
void Stats::updateSusOccLowered(const graph::node& v, const uint& numA, const uint& numS, const rtl& now) {
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
void Stats::computeOcc(const rtl& totalTime) {
	const uint n = (uint)avOcc.size();
	for (graph::node v = 0; v < n; ++v) avOcc[v] /= totalTime;
	for (graph::node v = 0; v < n; ++v) avIOcc[v] /= totalTime;
	for (graph::node v = 0; v < n; ++v) avSOcc[v] /= totalTime;
}
//void Stats::netInfAgMeanToFile(const rtl& Ws, const uint& numAg) {
//	const uint n = (uint)avOcc.size();
//	rtl avInfAgents = 0;
//	for (graph::node v = 0; v < n; ++v) avInfAgents += avIOcc[v];
//	netOccData
//		<< numAg << "\t"
//		<< Ws << "\t"
//		<< avInfAgents
//		<< std::endl;
//}
//void Stats::occToFile(const rtl& Ws, const uint& numAg) {
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