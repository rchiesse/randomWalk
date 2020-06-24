#include "graph.h"
#include "reporter.h"
#include "stats.h"

namespace sim { // Simulator's namespace.

//Main structures
enum class action{walk, recover};
struct job {
	agent ag;														// ----> Agent ID.
	real time;
	action a;

	job(const agent&& _ag, const real&& _time, const action&& _action) : ag(_ag), time(_time), a(_action) {}
	job() : job(INT32_MAX, INT32_MAX, action::walk) {}
	job(const agent& _ag, const real&& _time, const action&& _action) : ag(_ag), time(_time), a(_action) {}
	job(const job& other) { *this = other; }						// ----> Copy Contructor. Jobs will be handled by an STL vector (under a priority-queue approach --- see the "schedule" declaration). It thus requires a copy contructor to be explicitly defined. 

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


//Agent control variables
using std::vector; using graph::node;
vector<real> exposure		(NUM_AGENTS);							// ----> The total time interval a susceptible agent remained exposed to infected individuals at a given node. Susceptible agents have this value "zeroed" every time they enter a new node.
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
void enterNodeAsSus (const agent& ag, const node& v, const real& now);
void checkinAsSus   (const agent& ag, const node& v);
void enterNodeAsInf (const agent& ag, const node& v, const real& now);
void checkinAsInf   (const agent& ag, const node& v);
void leaveNodeAsSus (const agent& ag, const node& v, const real& now);
void checkoutAsSus  (const agent& ag, const node& v);
void leaveNodeAsInf (const agent& ag, const node& v, const real& now);
void checkoutAsInf  (const agent& ag, const node& v);
// Determines (i) whether or not an exposed, susceptible agent 'ag' will become infected and (ii) the next node 'ag' is going to visit.
void fateAndNextNode(const agent& ag, const real& now);
//Defines the next node an agent is going to visit.
const node& nextNodeForSus(const node& _currNode);
const node& nextNodeForInf(const node& _currNode);
void nextJob(job& _j, real& _time);
const node& randomLCCNode();
} // Namespace sim


/* MAIN */
int main() {
	sim::Reporter::logTimestamp("Simulation requested.");
	sim::Reporter::startChronometer();
#ifndef GENERATE_NETWORK
	graph::Graph::readGraph(SOURCE_FILE);
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
	long real Ce = sim::C * exp(B_MINUS_G * t);
	return std::max((_1_MINUS_G_OVER_B) * (Ce / (1.0 + Ce)), (long real)0.0);
}
real sim::i_t_pfx(const real& t) {
	long real Ce = sim::C_pfx * exp(B_MINUS_G_pfx * t);
	return std::max((_1_MINUS_G_OVER_B_pfx) * (Ce / (1.0 + Ce)), (long real)0.0);
}
#endif

#ifdef SOLVE_NUMERICALLY
//real sim::didt(const real& i) { return  _g * (i * (C_1 * i + C_2)) / (i + _h); }
real sim::didt(const real& i) {
	const real _crowd_increment = crowdFactor - ((1.0 - i) * i * (crowdFactor - 1.0));
	return  BETA_pfx * (1 - i) * i * _crowd_increment - GAMMA * i; 
}
//real sim::didt(const real& i) { return  BETA_pfx * i * (1 - i) - (GAMMA * i); }
void sim::rungeKutta4thOrder(const real& t0, const real& i0, const real& t, const real& h, const real& epsilon, vector<real>& saveToFile, uint& outputSize, const uint& outputGranularity, const real& largerDetailUntil) {
	uint totalSteps = (uint)((t - t0) / h) + 1;
	saveToFile.resize((uint64_t)largerDetailUntil + (totalSteps - (uint)largerDetailUntil)/outputGranularity);

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
	if (lastPos > 1) {	// ----> "will there be still any susceptible agent at 'v' once 'ag' has been removed?". If true, then the agent from the last position is copied to ag's. There's no problem if the outcoming agent comes to be the one at the 'top' position. In this case, both inner instructions become redundant, not wrong.
		list[indexWithinNode[ag]] = list[lastPos];
		indexWithinNode[list[lastPos]] = indexWithinNode[ag];
	}
	--lastPos;
}
void sim::checkoutAsInf  (const agent& ag, const node& v) {
	vector<agent>& list = iAgents[v];
	uint& lastPos = list[ELEMS];
	if (lastPos > 1) {	// ----> "will there be still any  infected  agent at 'v' once 'ag' has been removed?".  If true, then the agent from the last position is copied to ag's. There's no problem if the outcoming agent comes to be the one at the 'top' position. In this case, both inner instructions become redundant, not wrong.
		list[indexWithinNode[ag]] = list[lastPos];
		indexWithinNode[list[lastPos]] = indexWithinNode[ag];
	}
	--lastPos;
}
void sim::enterNodeAsSus (const agent& ag, const node& v, const real& now) {
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
	exposure[ag] = 0;
	if (numI > 0)
		iniExposureTime[ag] = now;
#ifdef PROTECTION_FX
	if (numS == 1)
		graph::Graph::updateHasS(v);
#endif
}
void sim::enterNodeAsInf (const agent& ag, const node& v, const real& now) {
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
	if (numI == 1) {
		const vector<uint>& list = sAgents[v];
		for (uint i = 1; i <= numS; ++i) {		// ----> We start by idx 1 since sAgents' position 0 is NOT an agent, but the list's actual size.
			iniExposureTime[list[i]] = now;
		}
#ifdef PROTECTION_FX
		graph::Graph::updateHasI(v);
#endif
	}
}
void sim::leaveNodeAsInf (const agent& ag, const node& v, const real& now) {
	checkoutAsInf(ag, v);
	const uint& numS = sInNode[v];				// ----> Number of susceptible agents currently hosted in v.
	uint&		numI = iInNode[v];				// ----> Number of infected agents currently hosted in v.
	--numI;
#ifdef OCCUPANCY
	stat::Stats::updateInfOccLowered(v, numI + numS, numI, now);
#endif
	if (numI == 0) {
		const vector<uint>& list = sAgents[v];
		for (uint i = 1; i <= numS; ++i) {		// ----> We start by idx 1 since sAgents' position 0 is NOT an agent, but the list's actual size.
			exposure[list[i]] += now - iniExposureTime[list[i]];
		}
#ifdef PROTECTION_FX
		graph::Graph::updateNoI(v);
#endif
	}
}
void sim::leaveNodeAsSus (const agent& ag, const node& v, const real& now) {
	checkoutAsSus(ag, v);
	const uint& numI = iInNode[v];				// ----> Number of infected agents currently hosted in v.
	uint&		numS = sInNode[v];				// ----> Number of susceptible agents currently hosted in v.
	--numS;
#ifdef OCCUPANCY
	stat::Stats::updateSusOccLowered(v, numI + numS, numS, now);
#endif
	if (iInNode[v] > 0) 
		exposure[ag] += now - iniExposureTime[ag]; 
#ifdef PROTECTION_FX
	if (numS == 0)
		graph::Graph::updateNoS(v);
#endif
}
void sim::walk			 (const agent& ag, const real& now) {
	if (isInfected[ag]) {
		leaveNodeAsInf(ag, currentNode[ag], now);
		enterNodeAsInf(ag, nextNodeForInf(currentNode[ag]), now);
	} else {
		leaveNodeAsSus(ag, currentNode[ag], now);
		fateAndNextNode(ag, now);
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
void sim::fateAndNextNode(const agent& ag, const real& now) {
#ifdef ESTIMATE_PROBS
	++stat::Stats::totalFate;
	if (exposure[ag] > 0) 
		stat::Stats::totalExposition += exposure[ag]; 
#endif
	if ((exposure[ag] > 0) && (EXPTau() < exposure[ag])) {
		isInfected[ag] = true;
		--stotal;
		++itotal;
#ifdef ESTIMATE_PROBS
		++stat::Stats::totalInfections;
#endif
#ifdef INFECTED_FRACTION
		stat::Stats::bufferizeIFrac(ag, now, 'I', itotal, NUM_AGENTS, OVERLOOK);
#endif
		enterNodeAsInf(ag, nextNodeForInf(currentNode[ag]), now);
		schedule.emplace(ag, now + EXPGamma(), action::recover); // ----> 'Recover' event is scheduled.
	} else { 
		enterNodeAsSus(ag, nextNodeForSus(currentNode[ag]), now);
	}
}
const graph::node& sim::nextNodeForSus(const node& _currNode) {
	using graph::Graph;
#ifdef CLIQUE
#ifdef PROTECTION_FX
	return Graph::nextNodeForS(U());
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
	return Graph::nextNodeForI(U());
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
const graph::node& sim::randomLCCNode() { 
	return graph::Graph::lcc[(size_t)floor(graph::Graph::lcc.size() * sim::U())]; 
}
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
		long real s1 = 0, s2 = 0;
		stat::Stats::roots(C13, C14, C15, s1, s2);
		std::cout << "\ti_inf_pfx from di/dt = 0: " << ((s1 > 0 && s1 < 1)? s1 : s2) << '\n';
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
			bool timeLimit = true;
			real& roundDuration = totalSimTime[round]; 
			job j;
			nextJob(j, now);
			while (now < T) {
				roundDuration += (now - roundDuration);
				if (j.a == action::recover) {
					recover(j.ag, now);
					if (itotal == 0) { timeLimit = false; break; }
				} else {
					walk(j.ag, now);
					schedule.emplace(j.ag, now + EXPLambda(), action::walk);
				}
				nextJob(j, now);
			}
			Stats::partialsAvDur(roundDuration);
#ifdef MEASURE_ROUND_EXE_TIME
			Reporter::stopChronometer((timeLimit) ? "done-TL" : "done-IV");	// ----> TL == Time Limit; IV == Infection Vanished.
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
			_totalTime += totalSimTime[i];
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
	RKdata << "Time\ti\n"; // ----> Header
	real _time = 0;
	for (size_t s = 0; s < largerDetailUntil; ++s) {
		RKdata << _time << '\t' << saveToFile[s] << '\n';
		_time += stepSize;
	}
	for (size_t s = largerDetailUntil; s < outputSize; ++s){
		RKdata << _time << '\t' << saveToFile[s] << '\n';
		_time += timeIncrement;
	}
	RKdata.close();
#endif //SOLVE_NUMERICALLY
}
