#include "graph.h"
#include "reporter.h"
#include "stats.h"


namespace sim { // Simulator's namespace.

//Main structures
enum class action{walk, recover, infect};
struct job {
	agent ag		= INT32_MAX;														// ----> Agent ID.
	real time		= 0;
	action a		= action::walk;

	//Infection only:
	uint validity_S			= UINT_MAX;		// ----> Snapshot. If the S-agent's current snapshot does not match this one, then the event is ignored and the infection does not occur.
	uint validity_I			= UINT_MAX;		
	agent infective			= UINT_MAX;
	//real delta			= UINT_MAX;		// ----> Time gap between the event being generated and then processed. If the S-agent's exposition time is smaller than this interval, then it does *not* get infected.
	//graph::node location	= UINT_MAX;		// ----> Node at which the infection event was generated.

	job(const agent&& _ag, const real&& _time, const action&& _action) : ag(_ag), time(_time), a(_action) {}
	job(const agent&& _ag, const real&& _time, const action&& _action, const agent&& _infective, const uint&& _snapshot_S, const uint&& _snapshot_I) : ag(_ag), time(_time), a(_action), infective(_infective), validity_S(_snapshot_S), validity_I(_snapshot_I) {}
	job() {}
	job(const agent& _ag, const real&& _time, const action&& _action) : ag(_ag), time(_time), a(_action) {}
	job(const agent& _ag, const real& _time, const action&& _action, const agent& _infective, const uint& _snapshot_S, const uint& _snapshot_I) : ag(_ag), time(_time), a(_action), infective(_infective), validity_S(_snapshot_S), validity_I(_snapshot_I) {}

	job(const job& other) { *this = other; }						// ----> Copy Contructor. Jobs will be handled by an STL vector (via priority-queue). It thus requires a copy contructor to be explicitly defined. 

};
struct earlier {
	bool operator()(const job& j1, const job& j2) {
		return j1.time > j2.time;
	}
};

std::priority_queue<job, std::vector<job>, earlier> schedule;


//Simulation variables:
uint stotal;														// ----> Up-to-date number of SUSCEPTIBLE agents during the simulation.
uint itotal;														// ----> Up-to-date number of INFECTED agents during the simulation.
real now;
std::vector<real> totalSimTime(ROUNDS, 0);							// ----> Total simulation time at each round, to average upon.
real beta2ndMmt;
real beta2ndMmt_logistic;
real beta2ndMmt_naive;
long real C_2ndMmt;
long real C_2ndMmt_logistic;
long real C_2ndMmt_naive;

//Agent control variables
using std::vector; using graph::node;
vector<real> exposure		(NUM_AGENTS);							// ----> The total time interval a susceptible agent remained exposed to infected individuals at a given node. Susceptible agents have this value "zeroed" every time they enter a new node.
vector<uint> snapshot		(NUM_AGENTS,0);							// ----> A counter assigned to each agent and incremented every time the agent walks. Its main purpose is to control whether or not an infection event should be applied: when the infection event is processed, even if the agent comes to be currently located in the node at which the event relates, we must ensure that it is not the case that between the event creation and its processing the agent went somewhere else and then came back to said node. We do this by checking the event's snapshot against the agent's. These must be the same for the infection to take place.
vector<real> iniExposureTime(NUM_AGENTS);							// ----> The moment a new exposition cicle starts for a susceptible agent at its current node.
vector<node> currentNode	(NUM_AGENTS);							// ----> Keeps track of the node an agent is currently located.
vector<uint> indexWithinNode(NUM_AGENTS);							// ----> By storing the agent's index in the list of its current node, we are able to always find any agent in O(1). **This is critical for the overall performance**
vector<bool> isInfected		(NUM_AGENTS);							// ----> Keeps track of each agent's current state.


//Node control variables
vector<agent> sInNode;												// ----> Keeps track, for each node v, of how many susceptible agents are in v at time t.
vector<agent> iInNode;												// ----> Keeps track, for each node v, of how many infected agents are in v at time t.
vector<vector<agent>> sAgents;										// ----> Up-to-date list of susceptible agents within each node v.
vector<vector<agent>> iAgents;										// ----> Up-to-date list of infected agents within each node v.
#ifdef CLIQUE
uint _neighbor;
#endif

// * GLOBAL UNIFORM(0,1) RANDOM-NUMBER GENERATOR * 
std::random_device rd;
std::mt19937_64 generator(rd());										// ----> "mt" = "Mersenne Twister".
std::uniform_real_distribution<real> distribution(0, 1);
inline static real U() { return distribution(generator); }
//auto U = bind(distribution, generator);


/* PROTOTYPES */
//Returns a random value uniformly drawn from the interval [0, openRange), i.e. a number between 0 and 'openRange - 1' (for ex., "randomInt(6)" will return an integer between 0 and 5). This is particularly useful when determining an agent's "next stop" during its random walk. In this case, "randomInt(<v's degree>)" provides a neighbor's index each time an agent walks out some node v.
static const uint randomInt(const uint& openRange);
//Exponential random-number generator.
static const real EXP(const real& param);
//Exponential random-number generator that determines whether or not an agent will become infected.
static const real EXPTau();
//Exponential random-number generator that sets the moment an infected agent will recover.
static const real EXPGamma();
//Exponential random-number generator that sets the next moment an agent walks.
static const real EXPLambda();
void resetVariables();
void readParams();
void runSimulation	(const uint& startingNumAg = 0, const uint& granularity = 5);
void walk		    (const agent& ag, const real& now);
void recover	    (const agent& ag, const real& now);

using graph::node;
void fate			(const agent& ag, const real& now, const agent& infective, const uint& validity_S, const uint& validity_I);		// ----> Determines whether or not an exposed, susceptible agent 'ag' will become infected.
void enterNodeAsSus (const agent& ag, const node& v, const real& now);
void checkinAsSus   (const agent& ag, const node& v);
void enterNodeAsInf (const agent& ag, const node& v, const real& now);
void checkinAsInf   (const agent& ag, const node& v);
void leaveNodeAsSus (const agent& ag, const node& v, const real& now);
void checkoutAsSus  (const agent& ag, const node& v);
void leaveNodeAsInf (const agent& ag, const node& v, const real& now);
void checkoutAsInf  (const agent& ag, const node& v);
//Defines the next node an agent is going to visit.
const node& nextNodeForSus(const node& _currNode);
const node& nextNodeForInf(const node& _currNode);
void nextJob(job& _j, real& _time);
void setBeta2ndMmt();
#ifndef CLIQUE
const node& randomLCCNode();
#endif
} // Namespace sim


/* MAIN */
int main() {
	sim::Reporter::logTimestamp("Simulation requested.");
	sim::Reporter::startChronometer();
#ifndef CLIQUE
	graph::Graph::readGraph(SOURCE_FILE);
	graph::Graph::set2ndMoment();
	sim::setBeta2ndMmt();
#endif
	sim::runSimulation(sim::STARTING_NUM_AG, sim::GRAN_NUM_AG);
	sim::Reporter::stopChronometer("\n\nSimulation completed");
	sim::Reporter::logTimestamp("End of simulation.");
	return 0;
}


/* IMPLEMENTATION */
static const uint sim::randomInt(const uint& openRange) { return (uint)(floor(openRange * U())); }
static const real sim::EXP(const real& param) { return -(1.0 / param) * log(U()); }
static const real sim::EXPTau()		{ return NEG_RECIPR_TAU		* log(U()); }
static const real sim::EXPGamma()	{ return NEG_RECIPR_GAMMA	* log(U()); }
static const real sim::EXPLambda()	{ return NEG_RECIPR_LAMBDA	* log(U()); }
#ifdef i_t_FROM_MODEL
real sim::i_t(const real& t) {
	double val = (t > 1000) ? 1000 : t;

	long real Ce = sim::C * exp(B_MINUS_G * val);
	return std::max((_1_MINUS_G_OVER_B) * (Ce / (1.0 + Ce)), (long real)0.0);
}
real sim::i_t_2ndMmt(const real& t) {
	double val = (beta2ndMmt - GAMMA) * t;
	if(val > 500) val = 500;

	//long real Ce = C_2ndMmt * exp((beta2ndMmt - GAMMA) * t);
	long real Ce = C_2ndMmt * exp(val);

	return std::max((1.0 - (GAMMA/beta2ndMmt)) * (Ce / (1.0 + Ce)), (long real)0.0);
}
real sim::i_t_2ndMmt_logistic(const real& t) {
	double val = (t > 1000) ? 1000 : t;

	long real Ce = C_2ndMmt_logistic * exp((beta2ndMmt_logistic - GAMMA) * val);
	return std::max((1.0 - (GAMMA / beta2ndMmt_logistic)) * (Ce / (1.0 + Ce)), (long real)0.0);
}
real sim::i_t_2ndMmt_naive(const real& t) {
	double val = (t > 1000) ? 1000 : t;

	long real Ce = C_2ndMmt_naive * exp((beta2ndMmt_naive - GAMMA) * val);
	return std::max((1.0 - (GAMMA / beta2ndMmt_naive)) * (Ce / (1.0 + Ce)), (long real)0.0);
}
//real sim::i_t_2ndMmt_sys(const real& t) {
//	for (size_t i = 0; i < sim::frequency.size(); i++) {
//
//		long real Ce = C_2ndMmt * exp((beta2ndMmt - GAMMA) * t);
//		return std::max((1.0 - (GAMMA / beta2ndMmt)) * (Ce / (1.0 + Ce)), (long real)0.0);
//	}
//}
real sim::i_t_pfx(const real& t) {
	long real Ce = sim::C_pfx * exp(B_MINUS_G_pfx * t);
	return std::max((_1_MINUS_G_OVER_B_pfx) * (Ce / (1.0 + Ce)), (long real)0.0);
}
void sim::setBeta2ndMmt() {
	////Excess:
	//vector<double> expAgBlock;				// ----> Expected number of agents within each degree block.
	//expAgBlock.resize(graph::Graph::largestDegree + 1, 0);
	//for (size_t _b = 1; _b < graph::Graph::frequency.size(); ++_b) {
	//	expAgBlock[_b] = (NUM_AGENTS * _b * graph::Graph::frequency[_b])/graph::Graph::averageDegree;
	//}
	//vector<double> expOcc;					// ----> Expected node occupancy within block _b, i.e. the expected number of agents within each node of some specific degree block.
	//expOcc.resize(graph::Graph::largestDegree + 1, 0);
	//for (size_t _b = 1; _b < graph::Graph::frequency.size(); ++_b) {
	//	if(graph::Graph::frequency[_b] > 0)
	//		expOcc[_b] = expAgBlock[_b] / (N * graph::Graph::frequency[_b]);
	//}
	//double generalExpOcc = 0;
	//double _generalExpOcc = 0;
	////uint activeBlocks = 0;
	//for (size_t b = 1; b < graph::Graph::frequency.size(); ++b) {
	//	generalExpOcc += expOcc[b] * ((double)b / N);
	//	_generalExpOcc += expOcc[b] * graph::Graph::frequency[b];
	//	//if (graph::Graph::frequency[b] > 0)
	//	//	++activeBlocks;
	//}
	//double excess = 0;
	//double countExcess = 0;
	//for (size_t b = 1; b < expOcc.size(); ++b) {
	//	if (expOcc[b] - 1.0 > 0) {
	//		//excess += expOcc[b] - (expOcc[b] - 1.0);
	//		excess += expOcc[b] - 1.0;
	//		++countExcess;
	//	}
	//}
	////excess /= NUM_AGENTS;

	//1st approach: 
	vector<double> beta_b(graph::Graph::frequency.size(), 0);
	for (size_t _b = 1; _b < graph::Graph::frequency.size(); ++_b){
		double _kb_ = (NUM_AGENTS * _b * N * graph::Graph::frequency[_b]) / (2 * graph::Graph::m);
		//if (_kb_ < N * graph::Graph::frequency[_b]) {
		if (_kb_ < N * graph::Graph::frequency[_b]) {
			beta_b[_b] = (LAMBDA * _b * (TAU / (2 * LAMBDA + TAU)) * _kb_) / graph::Graph::m;
		}
		else {
			beta_b[_b] = (LAMBDA * _b * (TAU / (LAMBDA + TAU)) * N * graph::Graph::frequency[_b]) / graph::Graph::m;
		}
	}
	beta2ndMmt = 0;
	for (size_t _b = 1; _b < graph::Graph::frequency.size(); ++_b) {
		beta2ndMmt += beta_b[_b];
	}

	//2nd approach:
	beta_b.resize(0);
	beta_b.resize(graph::Graph::frequency.size(), 0);
	double sigma = 0;
	double sigma_2 = 0;
	double psi = 0;
	uint validBlock = 0;
	double expBlock = graph::Graph::_2ndMmt / graph::Graph::averageDegree;
	double max_kb = (((double)NUM_AGENTS) * (double)(graph::Graph::frequency.size() - 1) * (double)N * graph::Graph::frequency[graph::Graph::frequency.size() - 1]) / (2.0 * graph::Graph::m);
	double variance = abs(graph::Graph::_2ndMmt - pow(graph::Graph::averageDegree, 2));
	double _diff = abs(graph::Graph::averageDegree - graph::Graph::_2ndMmt);
	double std_dev = abs(sqrt(variance));
	const double& avDg = graph::Graph::averageDegree;
	for (size_t _b = 1; _b < graph::Graph::frequency.size(); ++_b) {
		double _kb_ = (((double)NUM_AGENTS) * (double)_b * (double)N * graph::Graph::frequency[_b]) / (2.0 * graph::Graph::m);
		//double kbnb = _kb_ / N * graph::Graph::frequency[_b];
		double norm_kb = (_kb_/N)/ ((_kb_/N)+ std::exp(-(_kb_/N)));
		//beta_b[_b] = (LAMBDA * _b * (TAU / (LAMBDA + TAU + (LAMBDA * (1- norm_kbnb)) )) * (_kb_ * (1 - norm_kbnb))) / graph::Graph::m;
		beta_b[_b] = (LAMBDA * _b * (TAU / (2*LAMBDA + TAU) * norm_kb * N )) / graph::Graph::m;
		
		if (graph::Graph::frequency[_b] > 0) {
			++validBlock;
			//double dimmer = 1.0 - ((_kb_/(4.0 + _kb_)) / ((_kb_/ (4.0 + _kb_)) + exp(-(_kb_/ (4.0 + _kb_)))));
			//double dimmer = 1.0 - ((_kb_/(expBlock)) / ((_kb_/ (expBlock)) + exp(-(_kb_/ (expBlock)))));
			//double val = (_kb_) / (abs(std_dev - (log2(avDg) * _kb_)));
			//double val = (_kb_) / (1 + abs(variance - avDg));
			//double val = (_kb_) / (expBlock + (2 * _kb_));
			//double val = graph::Graph::averageDegree / graph::Graph::_2ndMmt;
			//double val = _kb_/(abs(expBlock - avDg) + _kb_);
			//double val = _kb_/(avDg + _kb_);
			//double val = _kb_/(log2(avDg)*max_kb);
			//double val = _kb_ / ((log2(expBlock) * max_kb) + _kb_);
			double val = _kb_/((log2(expBlock*max_kb))+_kb_);
			//double dimmer = 1.0 - (val / (val + exp(-val)));

			//Bom resultado:
			double nTau = (double)(TAU) / std::max(1.0, _kb_ * std::min(1.0, graph::Graph::_2ndMmt/avDg));
			
			//double dimmer = 1.0 - ((_kb_/(2*expBlock)) / ((_kb_/ (2*expBlock)) + exp(-(_kb_/ (2*expBlock)))));
			//double dimmer = 1.0 - ((1.0/(2*_kb_)) / (((1.0 / (2 * _kb_))) + exp(-(1.0 / (2 * _kb_)))));
			
			//double dimmer = (_kb_ < 2.0) ? 1.0 : 1.0 - (_kb_) / (_kb_ + exp(-(1.0 / (_kb_))));
			//sigma += ((double)TAU * dimmer) / (2 * LAMBDA + (TAU * dimmer));
			
			//sigma += (double)TAU / (LAMBDA + TAU + LAMBDA * (1.0 - (std::min(1.0,(_kb_/2.0)))));

			psi += 1.0 - pow(1.0 - (1.0 / ((double)N * graph::Graph::frequency[_b])), _kb_);

			//MUITO BOM NO BA (O NORMAL É MELHOR NO G(N,P)):
			sigma +=  (_kb_ < 1.5) ? TAU / (2*LAMBDA + TAU) : TAU / (LAMBDA + TAU);
			
			//sigma +=  (_kb_ < 1.5) ? TAU / (2*LAMBDA + TAU) : (TAU*dimmer) / (LAMBDA + (TAU*dimmer));
			sigma_2 += nTau / (2 * LAMBDA + nTau);
		}
	}
	//psi /= validBlock;
	sigma /= validBlock;
	sigma_2 /= validBlock;
	beta2ndMmt_logistic = 0;
	for (size_t _b = 1; _b < graph::Graph::frequency.size(); ++_b) {
		beta2ndMmt_logistic += beta_b[_b];
	}

	//double kOverN = (double)NUM_AGENTS / N;
	//double sigma = (double)TAU / (LAMBDA + TAU + LAMBDA * (1 - _k_n));
	//double _k_n = (kOverN) / (kOverN + std::exp(-(kOverN * (std::_Pi/2))));
	//beta2ndMmt_logistic = (double)((LAMBDA * _b * sigma * _k_n)) / graph::Graph::m;


	//double probC = 0;	// ----> Probability "correction".
	//probC =  TAU/(LAMBDA + TAU);

	//3rd approach:
	//double _aux = (excess * (1.0 / std::exp(log(N) - log(NUM_AGENTS))));
	//beta2ndMmt = (meetingRate * infectionProb * (std::floor((double)NUM_AGENTS - 2*generalExpOcc - _aux)) * graph::Graph::_2ndMmt) / pow(graph::Graph::averageDegree, 2);
	
	//beta2ndMmt_naive = (meetingRate * infectionProb * NUM_AGENTS * graph::Graph::_2ndMmt) / pow(graph::Graph::averageDegree, 2);
	
	//MUITO BOM NO BA (O NORMAL É MELHOR NO G(N,P)):
	//beta2ndMmt = (meetingRate * sigma * NUM_AGENTS * graph::Graph::_2ndMmt) / pow(graph::Graph::averageDegree, 2);
	beta2ndMmt = (meetingRate * infectionProb * NUM_AGENTS * graph::Graph::_2ndMmt) / pow(graph::Graph::averageDegree, 2);
	//beta2ndMmt_naive = (2*LAMBDA * NUM_AGENTS * psi * infectionProb) / graph::Graph::averageDegree;
	
	
	
	//Bom resultado na BA:
	//beta2ndMmt_naive = (meetingRate * sigma_2 * NUM_AGENTS * graph::Graph::_2ndMmt) / pow(graph::Graph::averageDegree, 2);
	
	//Matheus (prosseguir implementação):
	//beta2ndMmt_naive = (2 * LAMBDA * NUM_AGENTS * psi * sigma) / pow(graph::Graph::averageDegree,2);
	
	//beta2ndMmt_naive = (LAMBDA * infectionProb * N * graph::Graph::_2ndMmt) / (graph::Graph::m * avDg);

	C_2ndMmt = i_0 / (1.0 - i_0 - (GAMMA / beta2ndMmt));
	C_2ndMmt_logistic = i_0 / (1.0 - i_0 - (GAMMA / beta2ndMmt_logistic));
	C_2ndMmt_naive = i_0 / (1.0 - i_0 - (GAMMA / beta2ndMmt_naive));
}
#endif //i_t_FROM_MODEL

#ifdef SOLVE_NUMERICALLY
//real sim::didt(const real& i) { return  _g * (i * (C_1 * i + C_2)) / (i + _h); }
//real sim::didt(const real& i) {
//	//const real _crowd_increment = crowdFactor - ((1.0 - i) * i * (crowdFactor - 1.0));
//	//return  BETA_pfx * (1 - i) * i * _crowd_increment - GAMMA * i; 
//	return  BETA_pfx * (1 - i) * i - GAMMA * i;
//}
//real sim::didt(const real& i) { return  BETA_pfx * i * (1 - i) - (GAMMA * i); }
real sim::didt(const real& i) { return  beta2ndMmt * i * (1 - i) - (GAMMA * i); }



void sim::rungeKutta4thOrder(const real& t0, const real& i0, const real& t, const real& h, const real& epsilon, vector<real>& saveToFile, uint& outputSize, const uint& outputGranularity, const real& largerDetailUntil) {
	uint totalSteps = (uint)((t - t0) / h) + 1;
	saveToFile.resize((uint64_t)largerDetailUntil + (totalSteps - ((uint)largerDetailUntil) / outputGranularity) + 1);

	constexpr real one_sixth = 1.0 / 6.0;
	real k1, k2, k3, k4;
	real i = i0;
	bool end = false;
	saveToFile[0] = i0;
	++outputSize;

	//For the first 'largerDetailUntil' iterations every step is stored in a vector ('saveToFile'), for later being written to file.
	for (uint s = 1; s < largerDetailUntil; ++s) {
		k1 = h * didt(i);
		k2 = h * didt(i + 0.5 * k1);
		k3 = h * didt(i + 0.5 * k2);
		k4 = h * didt(i + k3);

		i = i + one_sixth * (k1 + 2 * k2 + 2 * k3 + k4);
		if (i < epsilon) {
			saveToFile[s] = 0;
			end = true;
			++outputSize;
			break;
		}
		saveToFile[s] = i;
		++outputSize;
	}
	if (end) return;

	//From the 'largerDetailUntil' iteration on, we afford to ignore 'outputGranularity'-size windows of values, so that the saved file does not grow explosively.
	for (uint s = (uint)largerDetailUntil; s < totalSteps; ++s) {
		k1 = h * didt(i);
		k2 = h * didt(i + 0.5 * k1);
		k3 = h * didt(i + 0.5 * k2);
		k4 = h * didt(i + k3);
		// ***  IF _T IS REQUIRED, USE THE VERSION BELOW  ***
		//k1 = h * didt(_t, i);
		//k2 = h * didt(_t + 0.5 * h, i + 0.5 * k1);
		//k3 = h * didt(_t + 0.5 * h, i + 0.5 * k2);
		//k4 = h * didt(_t + h, i + k3);
		//_t = t0 + h;
		i = i + one_sixth * (k1 + 2 * k2 + 2 * k3 + k4);
		if (i < epsilon) { 
			saveToFile[outputSize] = 0;
			++outputSize;
			break;
		}
		if (s % outputGranularity == 0) {
			saveToFile[outputSize] = i;
			++outputSize;
		}
	}
}
#endif

void sim::checkinAsInf   (const agent& ag, const node& v) {
	currentNode[ag] = v;
	vector<agent>& list = iAgents[v];
	if (list[ELEMS] == (list.size() - 1))
		list.resize(2 * (list.size()));
	uint& lastPos = list[ELEMS];
	++lastPos;
	indexWithinNode[ag] = lastPos;
	list[lastPos] = ag;
}
void sim::checkinAsSus   (const agent& ag, const node& v) {
	currentNode[ag] = v;
	vector<agent>& list = sAgents[v];
	if (list[ELEMS]  == (list.size() - 1))
		list.resize(2 * (list.size()));
	uint& lastPos = list[ELEMS];
	++lastPos;
	indexWithinNode[ag] = lastPos;
	list[lastPos] = ag;
}
void sim::checkoutAsSus  (const agent& ag, const node& v) {
	vector<agent>& list = sAgents[v];
	uint& lastPos = list[ELEMS];
	if (lastPos > 1) {	// ----> "Is it still going to remain any susceptible agent at 'v' once 'ag' has been removed?". If true, then the agent from the last position is copied to ag's. There's no problem if the outcoming agent happens to be the one at the 'top' position. In this case, both inner instructions become redundant, not wrong.
		list[indexWithinNode[ag]] = list[lastPos];
		indexWithinNode[list[lastPos]] = indexWithinNode[ag];
	}
	--lastPos;
}
void sim::checkoutAsInf  (const agent& ag, const node& v) {
	vector<agent>& list = iAgents[v];
	uint& lastPos = list[ELEMS];
	if (lastPos > 1) {	// ----> "Is it still going to remain any infected agent at 'v' once 'ag' has been removed?".  If true, then the agent from the last position is copied to ag's. There's no problem if the outcoming agent happens to be the one at the 'top' position. In this case, both inner instructions become redundant, not wrong.
		list[indexWithinNode[ag]] = list[lastPos];
		indexWithinNode[list[lastPos]] = indexWithinNode[ag];
	}
	--lastPos;
}
void sim::enterNodeAsSus (const agent& ag, const node& v, const real& now) {
	++snapshot[ag];
	checkinAsSus(ag, v);
	const uint& numI = iInNode[v];					// ----> Number of infected agents currently hosted in v.
	uint&		numS = sInNode[v];					// ----> Number of susceptible agents currently hosted in v.
	++numS;
#ifdef OCCUPANCY
	stat::Stats::updateSusOccRaised(v, numI + numS, numS, now);
#endif
#ifdef ESTIMATE_PROBS
	if (numS > 1) { 
		stat::Stats::meetings += (numS - 1); 
	}
	if (numI > 0) {
		stat::Stats::meetings += (numI); 
		stat::Stats::siMeetings += (numI);
	}
#endif
	//exposure[ag] = 0;
	if (numI > 0) {
		//iniExposureTime[ag] = now;
		const vector<uint>& list = iAgents[v];
		for (uint i = 1; i <= numI; ++i) {
			const real delta = EXPTau();
			//schedule.emplace(ag, now + delta, action::infect, snapshot[ag], delta, v);
			schedule.emplace(ag, now + delta, action::infect, list[i], snapshot[ag], snapshot[list[i]]);
		}
	}
#ifdef PROTECTION_FX
	if (numS == 1)
		graph::Graph::updateHasS(v);
#endif
}
void sim::enterNodeAsInf (const agent& ag, const node& v, const real& now) {
	++snapshot[ag];
	checkinAsInf(ag, v);
	const uint& numS = sInNode[v];				// ----> Number of susceptible agents currently hosted in v.
	uint&		numI = iInNode[v];				// ----> Number of infected agents currently hosted in v.
	++numI;
#ifdef OCCUPANCY
	stat::Stats::updateInfOccRaised(v, numI + numS, numI, now);
#endif
#ifdef ESTIMATE_PROBS
	if (numS > 0) {
		stat::Stats::siMeetings += (numS);
		stat::Stats::meetings += (numS);
	}
	if (numI > 1) {
		stat::Stats::meetings += (numI - 1);
	}
#endif
		//for (uint i = 1; i <= numS; ++i) {		// ----> We start by idx 1 since sAgents' position 0 is NOT an agent, but the list's actual size.
		//	iniExposureTime[list[i]] = now;
		//}
#ifdef PROTECTION_FX
	if (numI == 1) {
		graph::Graph::updateHasI(v);
	}
#endif
	const vector<uint>& list = sAgents[v];
	for (uint i = 1; i <= numS; ++i) {
		const real delta = EXPTau();
		schedule.emplace(list[i], now + delta, action::infect, ag, snapshot[list[i]], snapshot[ag]);
		//job(_ag, _time, _action, _infective, _snapshot_S, _snapshot_I)
	}
}
void sim::leaveNodeAsInf (const agent& ag, const node& v, const real& now) {
	++snapshot[ag];
	checkoutAsInf(ag, v);
	//const uint& numS = sInNode[v];				// ----> Number of susceptible agents currently hosted in v.
	uint&		numI = iInNode[v];				// ----> Number of infected agents currently hosted in v.
	--numI;
#ifdef OCCUPANCY
	stat::Stats::updateInfOccLowered(v, numI + numS, numI, now);
#endif
		//const vector<uint>& list = sAgents[v];
		//for (uint i = 1; i <= numS; ++i) {		// ----> We start by idx 1 since sAgents' position 0 is NOT an agent, but the list's actual size.
		//	exposure[list[i]] += now - iniExposureTime[list[i]];
		//}
#ifdef PROTECTION_FX
	if (numI == 0) {
		graph::Graph::updateNoI(v);
	}
#endif
}
void sim::leaveNodeAsSus (const agent& ag, const node& v, const real& now) {
	++snapshot[ag];
	checkoutAsSus(ag, v);
	//const uint& numI = iInNode[v];				// ----> Number of infected agents currently hosted in v.
	uint&		numS = sInNode[v];				// ----> Number of susceptible agents currently hosted in v.
	--numS;
#ifdef OCCUPANCY
	stat::Stats::updateSusOccLowered(v, numI + numS, numS, now);
#endif
	//if (numI > 0) 
	//	exposure[ag] += now - iniExposureTime[ag]; 
#ifdef PROTECTION_FX
	if (numS == 0)
		graph::Graph::updateNoS(v);
#endif
}
void sim::walk(const agent& ag, const real& now) {
	node v = (isInfected[ag]) ? nextNodeForInf(currentNode[ag]) : nextNodeForSus(currentNode[ag]);
	if (v == currentNode[ag])
		return;

	if (isInfected[ag]) {
		leaveNodeAsInf(ag, currentNode[ag], now);
		enterNodeAsInf(ag, v, now);
	} else {
		leaveNodeAsSus(ag, currentNode[ag], now);
		enterNodeAsSus(ag, v, now);
	}
}
void sim::recover		 (const agent& ag, const real& now) {
	isInfected[ag] = false;
	++stotal;
	--itotal;
#ifdef INFECTED_FRACTION
	stat::Stats::bufferizeIFrac(ag, now, 'R', itotal, NUM_AGENTS, OVERLOOK);
#endif
	const graph::node& v = currentNode[ag];
	leaveNodeAsInf(ag, v, now);
	enterNodeAsSus(ag, v, now);
}
void sim::fate(const agent& ag, const real& now, const agent& infective, const uint& validity_S, const uint& validity_I) {
#ifdef ESTIMATE_PROBS
	++stat::Stats::totalFate;
	if (exposure[ag] > 0) 
		stat::Stats::totalExposition += exposure[ag]; 
#endif
	if (validity_S != snapshot[ag] || validity_I != snapshot[infective]) {
		return;	// ----> Event became obsolete.
	}
	//if(iInNode[v] > 0)
	//	exposure[ag] += now - iniExposureTime[ag];
	//if (delta <= exposure[ag]) {
		isInfected[ag] = true;
		--stotal;
		++itotal;
#ifdef ESTIMATE_PROBS
		++stat::Stats::totalInfections;
#endif
#ifdef INFECTED_FRACTION
		stat::Stats::bufferizeIFrac(ag, now, 'I', itotal, NUM_AGENTS, OVERLOOK);
#endif
		leaveNodeAsSus(ag, currentNode[ag], now);
		enterNodeAsInf(ag, currentNode[ag], now);
		schedule.emplace(ag, now + EXPGamma(), action::recover); // ----> 'Recover' event is scheduled.
	//} 
}
const graph::node& sim::nextNodeForSus(const node& _currNode) {
	using graph::Graph;
#ifdef CLIQUE
#ifdef PROTECTION_FX
#ifdef PROPORTIONAL
	return Graph::nextNodeForS(U(), itotal, stotal);
#else
	return Graph::nextNodeForS(U());
#endif //PROPORTIONAL
#else
	return Graph::g[randomInt(graph::Graph::n)];
#endif //PROTECTION_FX
#else
#ifdef PROTECTION_FX
	return Graph::nextNodeForS(_currNode, U());
#else
	return Graph::g[_currNode][randomInt((uint)Graph::g[_currNode].size())];
#endif //PROTECTION_FX
#endif //CLIQUE
}
const graph::node& sim::nextNodeForInf(const node& _currNode) {
	using graph::Graph;
#ifdef CLIQUE
#ifdef PROTECTION_FX
#ifdef PROPORTIONAL
	return Graph::nextNodeForI(U(), itotal, stotal);
#else
	return Graph::nextNodeForI(U());
#endif //PROPORTIONAL
#else
	return Graph::g[randomInt(Graph::n)];
#endif //PROTECTION_FX
#else
#ifdef PROTECTION_FX
	return Graph::nextNodeForI(_currNode, U());
#else
	return Graph::g[_currNode][randomInt((uint)Graph::g[_currNode].size())];
#endif //PROTECTION_FX
#endif //CLIQUE
}
#ifndef CLIQUE
const graph::node& sim::randomLCCNode() { 
	return graph::Graph::lcc[(size_t)floor(graph::Graph::lcc.size() * sim::U())]; 
}
#endif
void sim::nextJob(job& _j, real& _time) {
	_j = schedule.top();
	schedule.pop();
	_time = _j.time;
}
void sim::readParams() {
	//TODO
}
void sim::resetVariables() {
	itotal = (I_0 > NUM_AGENTS) ? NUM_AGENTS : I_0;
	stotal = NUM_AGENTS - itotal;
	now = 0;
	while (!schedule.empty()) schedule.pop();	// ----> Needed from the 2nd round on.
	for (agent a = 0; a < isInfected.size(); ++a) isInfected[a] = false;
	snapshot.resize(NUM_AGENTS, 0);
	using graph::Graph; using graph::node;
	for (node v = 0; v < Graph::n; ++v) iAgents[v][ELEMS] = 0;
	for (node v = 0; v < Graph::n; ++v) sAgents[v][ELEMS] = 0;
	for (node v = 0; v < Graph::n; ++v) sInNode[v] = 0;
	for (node v = 0; v < Graph::n; ++v) iInNode[v] = 0;
}
void sim::runSimulation(const uint& startingNumAg, const uint& granularity) {
	using graph::Graph; using graph::node; using stat::Stats; 
#ifdef CLIQUE
	Graph::n = N;
#ifdef PROTECTION_FX
	Graph::gs.resize(Graph::n);
	Graph::gi.resize(Graph::n);
#else
	Graph::g.resize(Graph::n);
#endif //PROTECTION_FX
#endif
	sInNode.resize(Graph::n);									// ----> Keeps track, for each node v, of how many susceptible agents are in v at time t.
	iInNode.resize(Graph::n);									// ----> Keeps track, for each node v, of how many infected agents are in v at time t.
	sAgents.resize(Graph::n);									// ----> Up-to-date list of susceptible agents within each node v.
	iAgents.resize(Graph::n);									// ----> Up-to-date list of infected agents within each node v.
	for (node v = 0; v < Graph::n; ++v) sAgents[v].resize((size_t)LIST_INI_SZ + 1);	// ----> The extra spot is to store the actual number of elements in the list, which may differ from the container size.
	for (node v = 0; v < Graph::n; ++v) iAgents[v].resize((size_t)LIST_INI_SZ + 1);	// ----> The extra spot is to store the actual number of elements in the list, which may differ from the container size.

#ifdef PROTECTION_FX
	Graph::setProbs();
#endif

	uint _startingNumAg = ((startingNumAg != 0) && (startingNumAg < NUM_AGENTS)) ? startingNumAg : NUM_AGENTS;
	uint step = (_startingNumAg == NUM_AGENTS) ? 1 : (granularity == 0) ? 1 : granularity;
	const uint span = NUM_AGENTS - _startingNumAg + 1;
	uint numScenarios = (span / step) + 1 * (span % step);
	uint scenario = 0;
	uint _numAgents;
	while (scenario < numScenarios){
		_numAgents = _startingNumAg + (scenario * step);
		
		//Agent-wise reset:
#ifdef CLIQUE
#ifdef PROTECTION_FX
		for (node v = 0; v < Graph::n; ++v) Graph::gs[v] = v;		// ----> This will be important during the "protection effect" process, as the nodes will be re-arranged. 
		for (node v = 0; v < Graph::n; ++v) Graph::gi[v] = v;		// ----> This will be important during the "protection effect" process, as the nodes will be re-arranged. 
#else
		for (node v = 0; v < Graph::n; ++v) Graph::g[v] = v;		// ----> This will be important during the "protection effect" process, as the nodes will be re-arranged. 
#endif //PROTECTION_FX
#ifdef PROTECTION_FX
		Graph::resetAgentIdx();
#endif //PROTECTION_FX
#endif //CLIQUE
		Stats::resetAvDur();
		for (uint r = 0; r < ROUNDS; ++r) totalSimTime[r] = 0;
		itotal = (I_0 > NUM_AGENTS) ? NUM_AGENTS : I_0;
		
		//Initializing the applicable stats:
#ifdef INFECTED_FRACTION
		Stats::initStream(Ws, Wi, stat::streamType::infFrac, SHORT_LABEL, Graph::n, _numAgents, TAU, GAMMA, LAMBDA, T, ROUNDS);
#endif
#ifdef ESTIMATE_PROBS
		Stats::initStream(Ws, Wi, stat::streamType::probs, SHORT_LABEL, Graph::n, _numAgents, TAU, GAMMA, LAMBDA, T, ROUNDS);
#endif
#ifdef OCCUPANCY
		Stats::initOcc(Graph::n);
		Stats::initStream(Ws, Wi, stat::streamType::occupancy,  SHORT_LABEL, Graph::n, _numAgents, TAU, GAMMA, LAMBDA, T, ROUNDS);
#endif
		Stats::initStream(Ws, Wi, stat::streamType::avDuration, SHORT_LABEL, Graph::n, _numAgents, TAU, GAMMA, LAMBDA, T, ROUNDS);
		Reporter::startChronometer("\n\n\nRunning scenario " + std::to_string(scenario + 1) + "/" + std::to_string(numScenarios) + "...");
		Reporter::simulationInfo(ROUNDS, T, _numAgents, itotal, Graph::n, TAU, GAMMA, LAMBDA, Ws, Wi);
#ifdef PROTECTION_FX
		//long real s1 = 0, s2 = 0;
		//stat::Stats::roots(C13, C14, C15, s1, s2);
		//std::cout << "\ti_inf_pfx from di/dt = 0: " << ((s1 > 0 && s1 < 1)? s1 : s2) << '\n';
#endif
		for (uint round = 0; round < ROUNDS; ++round) {
#ifdef MEASURE_ROUND_EXE_TIME
			Reporter::startChronometer("\n  Round " + std::to_string(round + 1) + "...");
#else
			Reporter::progress(round);
#endif
			//Round-wise reset:
			resetVariables();
#ifdef PROTECTION_FX
			Graph::resetSchema();
#endif
			//Initially, we schedule a single 'walk' event for each agent. A new 'walk' job will then be created and scheduled for an agent the moment its current 'walk' job is processed.
			for (agent i = 0; i < _numAgents; ++i)
				schedule.emplace(i, EXPLambda(), action::walk);

			//Now we must define which agents are infected from the start, and then schedule their 'recover' event.
			//It is not a problem to always infect the first 'itotal' agents, since the starting node for each of them will be randomly set.
			for (agent i = 0; i < itotal; ++i) {
				isInfected[i] = true;
				schedule.emplace(i, EXPGamma(), action::recover);
			}

			// * DISTRIBUTING THE AGENTS ACROSS THE NETWORK *
			//Random distribution of the INFECTED agents:
#ifdef CLIQUE
			for (agent i = 0; i < itotal; ++i)
				enterNodeAsInf(i, randomInt(Graph::n), TIME_ZERO);
#else
			for (agent i = 0; i < itotal; ++i)
				enterNodeAsInf(i, randomLCCNode(), TIME_ZERO);
#endif

			//Random distribution of the SUSCEPTIBLE agents:
#ifdef CLIQUE
			for (agent i = itotal; i < _numAgents; ++i)
				enterNodeAsSus(i, randomInt(Graph::n), TIME_ZERO);
#else
			for (agent i = itotal; i < _numAgents; ++i)
				enterNodeAsSus(i, randomLCCNode(), TIME_ZERO);
#endif

			// * MAIN LOOP *
			bool earlyStop = true;
			real& roundDuration = totalSimTime[round]; 
			double timeLimit = T;
			job j;
			nextJob(j, now);
			while (now < timeLimit) {
				roundDuration += (now - roundDuration);
				switch (j.a) {
				case action::walk:
					walk(j.ag, now);
					schedule.emplace(j.ag, now + EXPLambda(), action::walk);
					break;
				case action::infect:
					fate(j.ag, now, j.infective, j.validity_S, j.validity_I);
					break;
				default:
					recover(j.ag, now);
					if (itotal == 0) { 
						earlyStop = false; 
						timeLimit = now; // ----> Exits the 'while'-loop by the next conditional test.
					}
					break;
				}
				nextJob(j, now);
			}
			Stats::partialsAvDur(roundDuration);
#ifdef MEASURE_ROUND_EXE_TIME
			Reporter::stopChronometer((earlyStop) ? "done-TL" : "done-IV");	// ----> TL == Time Limit; IV == Infection Vanished.
			Reporter::durationInfo(roundDuration);
#endif
#ifdef INFECTED_FRACTION
			Reporter::startChronometer(" STF...");	// ----> STF == "Saving To File"
			Stats::iFracToFile(OVERLOOK);
			Stats::endStream(stat::streamType::infFrac);
			Reporter::stopChronometer(" done");
#endif
#ifdef ESTIMATE_PROBS
			Stats::probsToFile();
			Stats::endStream(stat::streamType::probs);
#endif
		} // ** for (uint round = 0; round < ROUNDS; ++round)
#ifdef OCCUPANCY
		uint _max = 0, _min = _numAgents;
		graph::node whereMax, whereMin;
		Stats::getGlobalMaxOcc(_max, whereMax);
		Stats::getGlobalMinOcc(_min, whereMin);
		Reporter::minMaxOccInfo(_min, _max, whereMin, whereMax, Graph::n, _numAgents);
		real _totalTime = 0;
		for (uint i = 0; i < totalSimTime.size(); ++i)
			_totalTime += totalSimTime[i];	// ----> TODO: Verify whether it is the case of using an accumulator here instead of a for-loop. 
		Stats::computeOcc(_totalTime);
		Stats::writeToFile(stat::streamType::occupancy , Ws, Wi, _numAgents);
		Stats::endStream  (stat::streamType::occupancy);
#endif
		Stats::writeToFile(stat::streamType::avDuration, Ws, Wi, _numAgents);
		Stats::endStream  (stat::streamType::avDuration);
#ifndef MEASURE_ROUND_EXE_TIME
		Reporter::tell(" All rounds completed.\n");
#endif
		Reporter::stopChronometer("Scenario " + std::to_string(scenario + 1) + "/" + std::to_string(numScenarios) + " completed");
		Reporter::avSimTimeInfo(Stats::avDuration());
		++scenario;
	} // ** while (scenario < numScenarios)

	//TEMP:
	//std::ofstream R0data;
	//std::stringstream ss;
	//ss.precision(4);
	//std::string fileName;
	//ss << EXE_DIR << "/stats/R0_" << NWTK_LABEL << "_N" << N << "_AG" << K << "_T" << TAU << "_G" << GAMMA << "_L" << LAMBDA << "_STime" << T << "_R" << ROUNDS << ".csv";
	//fileName = ss.str();
	//R0data.open(fileName, std::ios::app);
	//if (Stats::isEmpty(R0data))
	//	R0data << "Ws\tWi\tbeta\tbetaApprox\tR0\tR0approx\ti_inf\ti_inf_approx\n";
	//R0data << Ws << '\t' << Wi << '\t' << BETA_pfx << '\t' << BETA_pfx_approx << '\t' << Ro_pfx << '\t' << Ro_pfx_approx << '\t' << std::max(i_inf_pfx, (long double)0.0) << '\t' << std::max(i_inf_pfx_approx, (long double)0.0) << '\n';
	//R0data.close();

#ifdef SOLVE_NUMERICALLY
	//Runge-Kutta:
	constexpr uint outputGranularity = 500;
	constexpr uint largerDetailUntil = 100;
	constexpr real stepSize = 0.1;
	constexpr real epsilon = 1.0 / N ;
	constexpr real timeIncrement = stepSize * outputGranularity;
	vector<real> saveToFile;
	uint outputSize = 0;
	rungeKutta4thOrder(0, FRAC_AG_INFECTED, T, stepSize, epsilon, saveToFile, outputSize, outputGranularity, largerDetailUntil);

	//Saving to file:
	std::ofstream RKdata;
	std::stringstream _ss;
	_ss.precision(4);
	_ss << EXE_DIR << "/stats/Runge-Kutta_" << NWTK_LABEL << "_Wi" << Wi << "_Ws" << Ws << "_N" << N << "_AG" << K << "_T" << TAU << "_G" << GAMMA << "_L" << LAMBDA << "_STime" << T << ".csv";
	std::string _fileName = _ss.str();
	RKdata.open(_fileName);
	RKdata << "Time i\n"; // ----> Header
	real _time = 0;
	for (size_t s = 0; s < largerDetailUntil; ++s) {
		RKdata << _time << ' ' << saveToFile[s] << '\n';
		_time += stepSize;
	}
	for (size_t s = largerDetailUntil; s < outputSize; ++s){
		RKdata << _time << ' ' << saveToFile[s] << '\n';
		_time += timeIncrement;
	}
	RKdata.close();
#endif //SOLVE_NUMERICALLY
}
