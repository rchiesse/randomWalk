#include "graph.h"
#include "reporter.h"
#include "stats.h"
#include "solver.h"
#include <numeric>		//std::accumulate()
#include <iomanip>		//std::fixed & std::setprecision


namespace sim { // Simulator's namespace.
// Input parameters
rtl T;																	// ----> Simulation time.
uint NUM_AGENTS;														// ----> Total number of agents in a simulation.
uint STARTING_NUM_AG;
uint GRAN_NUM_AG;
uint ROUNDS;															// ----> Number of simulation runs for a given setup. 
rtl TAU_aa;																// ----> Agent-to-agent transmissibility rate.
rtl TAU_al;																// ----> Agent-to-location transmissibility rate.
rtl TAU_la;																// ----> Location-to-agent transmissibility rate.
rtl GAMMA_a;															// ----> Recovery rate. 
rtl GAMMA_l;															// ----> Recovery rate. 
rtl LAMBDA;																// ----> Walk rate. 
rtl FRAC_AG_INFECTED;													// ----> Fraction of AGENTS initially infected (i.e. when the simulation starts).
rtl FRAC_ST_INFECTED;													// ----> Fraction of SITES initially infected (i.e. when the simulation starts).
uint ABS_INFECTED;														// ----> Absolute number of agents initially infected (i.e. when the simulation starts). This value is used whenever set to any value > 0, in which case it overrides 'FRAC_AG_INFECTED'. To use 'FRAC_AG_INFECTED' instead, set 'ABS_INFECTED = 0'.
rtl EARLY_MOBILITY;														// ----> Time dedicated for agents to freely walk through the network before the epidemic starts.


rtl Ws = 1.0;															// ----> Susceptible-agents' tolerance to enter nodes that contain infected agents, such that 0 <= Ws <= 1. This is the "s-protection-effect" single parameter.
rtl Wi = 1.0;															// ----> Infected-agents' tolerance to enter nodes that contain susceptible agents, such that 0 <= Wi <= 1. This is the "i-protection-effect" single parameter.

#ifdef PROTECTION_FX
//#define PROPORTIONAL													// ----> Promotes risk-tolerance proportional to the current number of infectives, according to the function w(i) = (1-i)^r, for some "rejection force" r.
#ifdef PROPORTIONAL
//If defined, we consider that the risk-tolerance is proportional to the number of infectives: w(i) = (1-i)^r, for some "rejection force" r.
static constexpr rtl _r = 1000.0;		// Rejection force.
#endif //PROPORTIONAL
#endif //PROTECTION_FX


// Auxiliary constants
uint I_0;
rtl i_0;
rtl meetingRate;
rtl SIGMA_aa;
rtl SIGMA_al;
rtl SIGMA_la;
rtl NEG_RECIPR_LAMBDA;
rtl NEG_RECIPR_GAMMA_a;
rtl NEG_RECIPR_GAMMA_l;
rtl NEG_RECIPR_TAU_aa;
rtl NEG_RECIPR_TAU_al;
rtl NEG_RECIPR_TAU_la;
uint ELEMS;																// ----> Zero here means "the first position of the container". Both 'sAgents' and 'iAgents' lists store their current number of elements in their respective first positions. 
rtl TIME_ZERO;
uint K;
uint LIST_INI_SZ;														// ----> Initial size of both 'sAgents' and 'iAgents' lists. Every time a node's list become full, its size gets doubled. Although arbitrary, the initial value provided here aims at reducing both the number of times a doubling operation is required and the vector's final size.

// Output control:
uint OVERLOOK;															// ----> Depending on the initial settings, the simulation may generate a firehose of data, which in turn becomes highly inconvenient (or even prohibitive) for being written to a file (in terms of either space and time). In such cases, it is strongly adviseable to purposely overlook a portion of the events generated per time unit.

//Main structures
enum class action{walk, recoverAg, recoverSite, agInfectAg, agInfectSite, siteInfectAg };
struct job {
	uint target		= INT32_MAX;										// ----> Target ID (either an agent or a site).
	rtl time		= 0;
	action a		= action::walk;

	//Infection only:
	uint validity_S			= UINT_MAX;									// ----> S-agent's snapshot from the moment the event was generated. If the S-agent's snapshot by the time the infection event is processed does not match this one, then the event is ignored and the infection does not occur.
	uint validity_I			= UINT_MAX;									// ----> I-agent's snapshot from the moment the event was generated. If the I-agent's snapshot by the time the infection event is processed does not match this one, then the event is ignored and the infection does not occur.		
	uint infective			= UINT_MAX;									// ----> ID of the infective actor (either an agent or a site). It will be necessary to later retrieve its snapshot by the time the infection event is processed. If this does not match 'validity_I', then the event is ignored and the infection does not occur.
	uint streamer			= UINT_MAX;									// ----> ID of the agent whose queue this event come from. This information makes it possible to get the next event from the correct queue (as it may come either from the S-agent or the I-agent) and insert it into the main queue.

	job() {}
	job(const uint&& _target, const rtl&& _time, const action&& _action) : target(_target), time(_time), a(_action) {}
	job(const uint& _target, const rtl&& _time, const action&& _action) : target(_target), time(_time), a(_action) {}
	job(const uint&& _target, const rtl&& _time, const action&& _action, const uint&& _infective, const uint&& _snapshot_S, const uint&& _snapshot_I, const uint&& _streamer) : target(_target), time(_time), a(_action), infective(_infective), validity_S(_snapshot_S), validity_I(_snapshot_I), streamer(_streamer){}
	job(const uint& _target, const rtl& _time, const action&& _action, const uint& _infective, const uint& _snapshot_S, const uint& _snapshot_I, const uint& _streamer) : target(_target), time(_time), a(_action), infective(_infective), validity_S(_snapshot_S), validity_I(_snapshot_I), streamer(_streamer) {}

	//Sites:
	job(const uint&& _target, const rtl&& _time, const action&& _action, const uint&& _infective, const uint&& _snapshot_S, const uint&& _snapshot_I) : target(_target), time(_time), a(_action), infective(_infective), validity_S(_snapshot_S), validity_I(_snapshot_I) {}
	job(const uint& _target, const rtl& _time, const action&& _action, const uint& _infective, const uint& _snapshot_S, const uint& _snapshot_I) : target(_target), time(_time), a(_action), infective(_infective), validity_S(_snapshot_S), validity_I(_snapshot_I) {}

	job(const job& other) { *this = other; }							// ----> Copy Contructor. Jobs will be handled by an STL vector (via priority-queue). It thus requires a copy contructor to be explicitly defined. 

};
struct earlier {
	bool operator()(const job& j1, const job& j2) {
		return j1.time > j2.time;
	}
};

std::vector<std::priority_queue<job, std::vector<job>, earlier>> myEvents;
std::priority_queue<job, std::vector<job>, earlier> schedule;


//Simulation variables:
uint saTotal;															// ----> Up-to-date number of SUSCEPTIBLE AGENTS during the simulation.
uint iaTotal;															// ----> Up-to-date number of INFECTED AGENTS during the simulation.
uint ilTotal = 0;														// ----> Up-to-date number of INFECTED SITES during the simulation.
rtl now;
std::vector<rtl> totalSimTime;											// ----> Total simulation time at each round, to average upon.


//Agent control variables
using std::vector; using graph::node;
vector<uint> snapshot_a		;											// ----> A counter assigned to each AGENT and incremented every time the agent walks. Its main purpose is to control whether or not an infection event should be applied: when the infection event is processed, even if the agent comes to be currently located in the node at which the event relates, we must ensure that it is not the case that between the event creation and its processing the agent went somewhere else and then came back to said node. We do this by checking the event's snapshot against the agent's. These must be the same for the infection to take place.
vector<uint> snapshot_l		;											// ----> A counter assigned to each LOCATION (NODE). Its purpose is the same of 'snapshot_a'.
vector<node> currentNode	;											// ----> Keeps track of the node an agent is currently located.
vector<uint> indexWithinNode;											// ----> By storing the agent's index in the list of its current node, we are able to always find any agent in O(1). **This is critical for the overall performance**
vector<bool> isInfectedAg	;											// ----> Keeps track of each agent's current state.
vector<bool> isInfectedSite ;											// ----> Keeps track of each site's current state.
#ifdef QUIET_INFECTION
vector<bool> quiet			(NUM_AGENTS, false);						// ----> Once an S-ag a is infected within the node v, this control variable governs whether a is allowed to further infect other S-agents inside v even before moving to somewhere else first. If set to 'false', then a will only become a propagator once it moves to another node first.
#endif
//Node control variables
vector<agent> sInNode;													// ----> Keeps track, for each node v, of how many susceptible agents are in v at time t.
vector<agent> iInNode;													// ----> Keeps track, for each node v, of how many infected agents are in v at time t.
vector<vector<agent>> sAgents;											// ----> Up-to-date list of susceptible agents within each node v.
vector<vector<agent>> iAgents;											// ----> Up-to-date list of infected agents within each node v.

//Average number of hops until changing status:
//vector<rtl> IavNumHops;
//vector<rtl> SavNumHops;
//vector<uint> ITransitions;
//vector<uint> STransitions;
//vector<uint> hopsUntilI;
//vector<uint> hopsUntilS;


// * GLOBAL UNIFORM(0,1) RANDOM-NUMBER GENERATOR * 
std::random_device rd;
std::mt19937_64 generator(rd());										// ----> "mt" = "Mersenne Twister".
std::uniform_real_distribution<rtl> distribution(0, 1);
inline static rtl U() { return distribution(generator); }


/* PROTOTYPES */
//Returns a random value uniformly drawn from the interval [0, openRange), i.e. a number between 0 and 'openRange - 1' (for ex., "randomInt(6)" returns an integer between 0 and 5). Likewise, for some node v, "randomInt(<v's degree>)" uniformly chooses an agent's next-hop up from v's neighbors.
static const uint randomInt(const uint& openRange);
//Exponential random-number generators.
static const rtl EXP(const rtl& param);
static const rtl EXPTau_aa();
static const rtl EXPTau_al();
static const rtl EXPTau_la();
static const rtl EXPGamma_a();
static const rtl EXPGamma_l();
static const rtl EXPLambda();
void resetVariables();
void readParams();
void runSimulation	(const uint& startingNumAg = 0, const uint& granularity = 5);
void walk		    (const agent& ag, const rtl& now);
void quietWalk		(const agent& ag, const rtl& now);
void infectAg	    (const agent& ag, const rtl& now);
void recoverAg	    (const agent& ag, const rtl& now);
void recoverSite    (const uint& site, const rtl& now);

using graph::node;
void agFate_fromAg	(const agent& ag, const rtl& now, const uint& infective, const uint& validity_S, const uint& validity_I, const uint& streamer);		// ----> Determines whether or not an exposed, susceptible agent 'ag' will become infected.
void agFate_fromSite(const agent& ag, const rtl& now, const uint& infective, const uint& validity_S, const uint& validity_I);		// ----> Determines whether or not an exposed, susceptible agent 'ag' will become infected.
void siteFate		(const uint& v, const rtl& now, const uint& infective, const uint& validity_L, const uint& validity_I);		// ----> Determines whether or not an exposed, susceptible agent 'ag' will become infected.
void check_in		(const agent& ag, const node& v, vector<vector<agent>>& _where);
void check_out		(const agent& ag, const node& v, vector<vector<agent>>& _where);
void enterNodeAsSus (const agent& ag, const node& v, const rtl& now);
void enterNodeAsInf (const agent& ag, const node& v, const rtl& now);
void leaveNodeAsSus (const agent& ag, const node& v, const rtl& now);
void leaveNodeAsInf (const agent& ag, const node& v, const rtl& now);
const node& nextNodeForSus(const node& _currNode);
const node& nextNodeForInf(const node& _currNode);
void nextJob(job& _j, rtl& _time);
void setEnvironment();
const node& randomLCCNode();
} // Namespace sim


/* MAIN */
int main() {
	sim::Reporter::logTimestamp("Simulation requested.");
	sim::Reporter::startChronometer();
	graph::Graph::readGraph(SOURCE_FILE);
	sim::setEnvironment();

	////Cria um anel:
	//std::ofstream rede;
	//rede.open("ring-10k.txt");
	//uint tam = 1e4, n = tam + 1;
	//while (--n)
	//	rede << n << ' ' << n - 1 << '\n';
	//rede << n << ' ' << tam << '\n';
	//rede.close();
	//return 0;
	sim::runSimulation(sim::STARTING_NUM_AG, sim::GRAN_NUM_AG);
	sim::Reporter::stopChronometer("\n\nSimulation completed");

	//=============================================================
	//Temp:
	//rtl IgeneralHops = 0;
	//rtl SgeneralHops = 0;
	//for (int i = sim::NUM_AGENTS - 1; i >= 0; --i){
	//	sim::IavNumHops[i] /= (sim::ITransitions[i] > 0) ? sim::ITransitions[i] : 1;
	//	sim::SavNumHops[i] /= (sim::STransitions[i] > 0) ? sim::STransitions[i] : 1;
	//}
	//for (int i = sim::NUM_AGENTS - 1; i >= 0; --i) {
	//	IgeneralHops += sim::IavNumHops[i];
	//	SgeneralHops += sim::SavNumHops[i];
	//}
	//IgeneralHops /= sim::NUM_AGENTS;
	//SgeneralHops /= sim::NUM_AGENTS;
	//
	//std::cout << "\n\nAv. num. hops as I-agent: " << IgeneralHops;
	//std::cout << "\nAv. num. hops as S-agent: " << SgeneralHops;
	//
	//std::sort(sim::hopsUntilI.begin(), sim::hopsUntilI.end());
	//std::cout << "\n\n10 minimum #hops until becoming infected: \n";
	//for (size_t i = 0; i < 10; ++i) {
	//	std::cout << sim::hopsUntilI[i] << " ";
	//}
	//std::cout << "\n\n10 maximum #hops until becoming infected: \n";
	//for (size_t i = 0; i < 10; ++i){
	//	std::cout << sim::hopsUntilI[sim::hopsUntilI.size() - 1 - i] << " ";
	//}
	//std::cout << '\n';
	//=============================================================

	sim::Reporter::logTimestamp("End of simulation.");
	return 0;
}


/* IMPLEMENTATION */
static const uint sim::randomInt(const uint& openRange) { return (uint)(floor((double)openRange * U())); }
static const rtl sim::EXP(const rtl& param) { return -(1.0 / param) * log(U()); }
static const rtl sim::EXPTau_aa()	{ return NEG_RECIPR_TAU_aa	* log(U()); }
static const rtl sim::EXPTau_al()	{ return NEG_RECIPR_TAU_al	* log(U()); }
static const rtl sim::EXPTau_la()	{ return NEG_RECIPR_TAU_la	* log(U()); }
static const rtl sim::EXPGamma_a()	{ return NEG_RECIPR_GAMMA_a	* log(U()); }
static const rtl sim::EXPGamma_l()	{ return NEG_RECIPR_GAMMA_l	* log(U()); }
static const rtl sim::EXPLambda()	{ return NEG_RECIPR_LAMBDA	* log(U()); }

void sim::setEnvironment() {
	T					= 0.02;											// ----> Simulation time.
	NUM_AGENTS			= 20000;										// ----> Total number of agents in a simulation.
	STARTING_NUM_AG		= 1000000;
	GRAN_NUM_AG			= 1;
	ROUNDS				= 1;											// ----> Number of simulation runs for a given setup. 
	TAU_aa				= 100.0;										// ----> Agent-to-agent transmissibility rate.
	GAMMA_a				= 15000.0;										// ----> Recovery rate. 
	LAMBDA				= 1.0;											// ----> Walk rate.
	FRAC_AG_INFECTED	= 0.5;											// ----> Fraction of AGENTS initially infected (i.e. when the simulation starts).
	FRAC_ST_INFECTED	= 0.0;											// ----> Fraction of SITES initially infected (i.e. when the simulation starts).
	ABS_INFECTED		= 0;											// ----> Absolute number of agents initially infected (i.e. when the simulation starts). This value is used whenever it is set to any value > 0, in which case it overrides 'FRAC_AG_INFECTED'. To use 'FRAC_AG_INFECTED' instead, set 'ABS_INFECTED = 0'.
	
	//TAU_al				= 0.000001;										// ----> Agent-to-location transmissibility rate.
	//TAU_la				= 0.000001;										// ----> Location-to-agent transmissibility rate.
	//GAMMA_l				= 20000.0;										// ----> Recovery rate. 
	
	//1K:
	//10K:
	//FIXAR R_0!
	// 
	//10k
	//T = 2000.0; NUM_AGENTS = 125; TAU_aa = 1.0; GAMMA_a = 0.007; LAMBDA = 1.0;
	
	//100k
	//T = 2000.0; NUM_AGENTS = 200; TAU_aa = 1.0; GAMMA_a = 0.0014; LAMBDA = 1.0;

	//1M
	T = 5000.0; NUM_AGENTS = 100; TAU_aa = 1.0; GAMMA_a = 0.001; LAMBDA = 1.0;

#ifdef PROTECTION_FX
	Wi = 1.0;
	Ws = 1.0;
	//Wi = 0.65;
	//Ws = 0.65;
#else
	Wi = Ws = 1.0;	// ----> Do not change this line.
#endif
	//Other parameters:
	EARLY_MOBILITY = std::max((rtl)1.0, 50.0 * (1.0/LAMBDA));
	OVERLOOK			= 1;
	LIST_INI_SZ			= (uint)(round(std::max((rtl)2.0, (rtl)NUM_AGENTS / (3.0 * N))));
	
	I_0					= (ABS_INFECTED > 0) ? ABS_INFECTED : (uint)((rtl)NUM_AGENTS * FRAC_AG_INFECTED);
	i_0					= (rtl)I_0 / NUM_AGENTS;
	meetingRate			= 2.0 * LAMBDA / N;
	SIGMA_aa			= (TAU_aa / (2.0 * LAMBDA + TAU_aa));
	SIGMA_al			= (TAU_al / (LAMBDA + TAU_al));
	SIGMA_la			= (TAU_la / (LAMBDA + TAU_la));
	NEG_RECIPR_LAMBDA	= -(1.0 / LAMBDA);	
	NEG_RECIPR_GAMMA_a	= -(1.0 / GAMMA_a);	
	NEG_RECIPR_GAMMA_l	= -(1.0 / GAMMA_l);	
	NEG_RECIPR_TAU_aa	= -(1.0 / TAU_aa);	
	NEG_RECIPR_TAU_al	= -(1.0 / TAU_al);	
	NEG_RECIPR_TAU_la	= -(1.0 / TAU_la);	
	ELEMS				= 0;													// ----> Zero here means "the first position of the container". Both 'sAgents' and 'iAgents' lists store their current number of elements in their respective first positions. 
	TIME_ZERO			= 0;
	K					= NUM_AGENTS;

	totalSimTime	.resize(ROUNDS, 0);
	snapshot_a		.resize(NUM_AGENTS, 0);
	snapshot_l		.resize(N, 0);
	currentNode		.resize(NUM_AGENTS);
	indexWithinNode	.resize(NUM_AGENTS);
	isInfectedAg	.resize(NUM_AGENTS);
	isInfectedSite	.resize(N, false);
	myEvents		.resize(NUM_AGENTS);

#ifdef SOLVE_NUMERICALLY
	nT = TAU_aa;
	nL = LAMBDA;
	nG = GAMMA_a;
#endif 
	graph::Graph::setParams(NUM_AGENTS, Ws, Wi);
	Solver::setParams(nT, nL, nG, NUM_AGENTS, Wi, Ws);
	Stats::setParams(T, NUM_AGENTS, ROUNDS, TAU_aa, GAMMA_a, LAMBDA, Wi, Ws);
	graph::Graph::setBlockData();
	Solver::setBlockData();
	Stats::setBasename();
	Stats::initAvDur();
}


void sim::check_in (const agent& ag, const node& v, vector<vector<agent>>& _where) {
	currentNode[ag] = v;
	vector<agent>& list = _where[v];
	if (list[ELEMS] == (list.size() - 1))
		list.resize(2 * (list.size()));
	uint& lastPos = list[ELEMS];
	++lastPos;
	indexWithinNode[ag] = lastPos;
	list[lastPos] = ag;
}
void sim::check_out  (const agent& ag, const node& v, vector<vector<agent>>& _where) {
	vector<agent>& list = _where[v];
	uint& lastPos = list[ELEMS];
	if (lastPos > 1) {	// ----> "Will any susceptible agent still remain at 'v' once 'ag' is gone?". If true, then the agent from the last position is copied to ag's. There's no problem if the outcoming agent happens to be the one at the 'top' position. In this case, both inner instructions simply become redundant, but NOT WRONG.
		list[indexWithinNode[ag]] = list[lastPos];
		indexWithinNode[list[lastPos]] = indexWithinNode[ag];
	}
	--lastPos;
}
void sim::enterNodeAsSus (const agent& ag, const node& v, const rtl& now) {
	++snapshot_a[ag];
	check_in(ag, v, sAgents);
	const uint& numI = iInNode[v];					// ----> Number of infected agents currently hosted in v.
	uint&		numS = sInNode[v];					// ----> Number of susceptible agents currently hosted in v.
	++numS;
	
	if (numI) {
		const vector<uint>& list = iAgents[v];
		rtl delta;
		uint i = numI + 1;
		while (--i) {
			delta = EXPTau_aa();
			myEvents[ag].emplace(ag, now + delta, action::agInfectAg, list[i], snapshot_a[ag], snapshot_a[list[i]], ag);
		}

		//The most recent event is copied to the main queue and deleted from the agent's:
		schedule.emplace(myEvents[ag].top());
		myEvents[ag].pop();
	}
	//if (isInfectedSite[v]) {
	//	const rtl delta = EXPTau_la();
	//	schedule.emplace(ag, now + delta, action::siteInfectAg, v, snapshot_a[ag], snapshot_l[v]);
	//}
#ifdef PROTECTION_FX
	if (numS == 1)
		graph::Graph::updateHasS(v);
#endif

}
void sim::enterNodeAsInf (const agent& ag, const node& v, const rtl& now) {
	++snapshot_a[ag];
	check_in(ag, v, iAgents);
	const uint& numS = sInNode[v];				// ----> Number of susceptible agents currently hosted in v.
	uint&		numI = iInNode[v];				// ----> Number of infected agents currently hosted in v.
	++numI;

	if (numS) {
		const vector<uint>& list = sAgents[v];
		rtl delta;
		uint i = numS + 1;
		while(--i) {
			delta = EXPTau_aa();
			myEvents[ag].emplace(list[i], now + delta, action::agInfectAg, ag, snapshot_a[list[i]], snapshot_a[ag], ag);
		}

		//The most recent event is copied to the main queue and deleted from the agent's:
		schedule.emplace(myEvents[ag].top());
		myEvents[ag].pop();
	}
#ifdef PROTECTION_FX
	if (numI == 1)
		graph::Graph::updateHasI(v);
#endif
	//if (!isInfectedSite[v]) {
	//	const rtl delta = EXPTau_al();
	//	schedule.emplace(v, now + delta, action::agInfectSite, ag, snapshot_l[v], snapshot_a[ag]);
	//}
}
void sim::leaveNodeAsInf (const agent& ag, const node& v, const rtl& now) {
	++snapshot_a[ag];
	myEvents[ag] = {};		// ----> Any events remaining into ag's queue becomes invalid, so that we may simply erase them.
	check_out(ag, v, iAgents);
	uint&		numI = iInNode[v];				// ----> Number of infected agents currently hosted in v.
	--numI;

#ifdef PROTECTION_FX
	if (numI == 0)
		graph::Graph::updateNoI(v);
#endif
}
void sim::leaveNodeAsSus (const agent& ag, const node& v, const rtl& now) {
	++snapshot_a[ag];
	myEvents[ag] = {};							// ----> Any events remaining into ag's queue becomes invalid, so that we may simply erase them.
	check_out(ag, v, sAgents);
	uint&		numS = sInNode[v];				// ----> Number of susceptible agents currently hosted in v.
	--numS;

#ifdef PROTECTION_FX
	if (numS == 0)
		graph::Graph::updateNoS(v);
#endif
}
void sim::quietWalk(const agent& ag, const rtl& now) {
	node v = nextNodeForSus(currentNode[ag]);
#ifdef SELF_LOOPS
	if (v == currentNode[ag])
		return;
#endif
	leaveNodeAsSus(ag, currentNode[ag], now);
	enterNodeAsSus(ag, v, now);
}
void sim::walk(const agent& ag, const rtl& now) {
	node v = (isInfectedAg[ag]) ? nextNodeForInf(currentNode[ag]) : nextNodeForSus(currentNode[ag]);
#ifdef SELF_LOOPS
	if (v == currentNode[ag])
		return;
#endif
	if (isInfectedAg[ag]) {
		leaveNodeAsInf(ag, currentNode[ag], now);
		enterNodeAsInf(ag, v, now);
	} else {
		leaveNodeAsSus(ag, currentNode[ag], now);
		enterNodeAsSus(ag, v, now);
		//if (ag == 0)
		//	++hopsUntilI[hopsUntilI.size() - 1];
	}
}

void sim::infectAg(const agent& ag, const rtl& now) {
	isInfectedAg[ag] = true;
	--saTotal;
	++iaTotal;
	leaveNodeAsSus(ag, currentNode[ag], now);
	enterNodeAsInf(ag, currentNode[ag], now);
	schedule.emplace(ag, now + EXPGamma_a(), action::recoverAg); // ----> 'Recover' event is scheduled.
}

void sim::recoverAg		 (const agent& ag, const rtl& now) {
	isInfectedAg[ag] = false;
	++saTotal;
	--iaTotal;
#ifdef INFECTED_FRACTION
	Stats::bufferizeIFrac(ag, now, "Ra", iaTotal, ilTotal, NUM_AGENTS, OVERLOOK);
#endif
	const graph::node& v = currentNode[ag];
	leaveNodeAsInf(ag, v, now);
	enterNodeAsSus(ag, v, now);
}

void sim::recoverSite	(const uint& v, const rtl& now) {
	isInfectedSite[v] = false;
	--ilTotal;
#ifdef INFECTED_FRACTION
	Stats::bufferizeIFrac(v, now, "Rl", iaTotal, ilTotal, NUM_AGENTS, OVERLOOK);
#endif
	++snapshot_l[v];
}

void sim::agFate_fromAg(const agent& ag, const rtl& now, const uint& infective, const uint& validity_S, const uint& validity_I, const uint& streamer) {
	if (!myEvents[streamer].empty()) {
		schedule.emplace(myEvents[streamer].top());
		myEvents[streamer].pop();
	}
	if (validity_S != snapshot_a[ag] || validity_I != snapshot_a[infective])
		return;	// ----> Event became obsolete.
	
	infectAg(ag, now);

#ifdef INFECTED_FRACTION
	Stats::bufferizeIFrac(ag, now, "Ia", iaTotal, ilTotal, NUM_AGENTS, OVERLOOK);
#endif
}
void sim::agFate_fromSite(const agent& ag, const rtl& now, const uint& infective, const uint& validity_S, const uint& validity_I) {
	if (validity_S != snapshot_a[ag] || validity_I != snapshot_l[infective]) 
		return;	// ----> Event became obsolete.
	
	infectAg(ag, now);

#ifdef INFECTED_FRACTION
	Stats::bufferizeIFrac(ag, now, "Ia", iaTotal, ilTotal, NUM_AGENTS, OVERLOOK);
#endif
}
void sim::siteFate(const uint& v, const rtl& now, const uint& infective, const uint& validity_S, const uint& validity_I) {
	if (validity_S != snapshot_l[v] || validity_I != snapshot_a[infective]) 
		return;	// ----> Event became obsolete.
	
	isInfectedSite[v] = true;
	++ilTotal;
#ifdef INFECTED_FRACTION
	Stats::bufferizeIFrac(v, now, "Il", iaTotal, ilTotal, NUM_AGENTS, OVERLOOK);
#endif
	schedule.emplace(v, now + EXPGamma_l(), action::recoverSite); // ----> 'Recover' event is scheduled.
	const uint& numS = sInNode[v];				// ----> Number of susceptible agents currently hosted in v.
	const vector<uint>& list = sAgents[v];
	for (uint i = numS; i > 0; --i) {
		const rtl delta = EXPTau_la();
		schedule.emplace(list[i], now + delta, action::siteInfectAg, v, snapshot_a[list[i]], snapshot_l[v]);
	}
}
const graph::node& sim::nextNodeForSus(const node& _currNode) {
	using graph::Graph;
#ifdef PROTECTION_FX
	return Graph::nextNodeForS(_currNode, U());
#else
	return Graph::g[_currNode][randomInt((uint)Graph::g[_currNode].size())];
#endif //PROTECTION_FX
}
const graph::node& sim::nextNodeForInf(const node& _currNode) {
	using graph::Graph;
#ifdef PROTECTION_FX
	return Graph::nextNodeForI(_currNode, U());
#else
	return Graph::g[_currNode][randomInt((uint)Graph::g[_currNode].size())];
#endif //PROTECTION_FX
}
const graph::node& sim::randomLCCNode() { 
	return graph::Graph::lcc[(size_t)floor(graph::Graph::lcc.size() * sim::U())]; 
}
void sim::nextJob(job& _j, rtl& _time) {
	_j = schedule.top();
	schedule.pop();
	_time = _j.time;
}
void sim::readParams() {
	//TODO
}
void sim::resetVariables() {
	iaTotal = 0;
	saTotal = NUM_AGENTS;
	//iaTotal = (I_0 > NUM_AGENTS) ? NUM_AGENTS : I_0;
	//saTotal = NUM_AGENTS - iaTotal;
	now = 0;
	schedule = {};	// ----> Needed from the 2nd round on.
	for (agent a = 0; a < isInfectedAg.size(); ++a) isInfectedAg[a] = false;
	snapshot_a.resize(NUM_AGENTS, 0);
	using graph::Graph; using graph::node;
	for (node v = 0; v < Graph::n; ++v) iAgents[v][ELEMS] = 0;
	for (node v = 0; v < Graph::n; ++v) sAgents[v][ELEMS] = 0;
	for (node v = 0; v < Graph::n; ++v) sInNode[v] = 0;
	for (node v = 0; v < Graph::n; ++v) iInNode[v] = 0;

#ifdef PROTECTION_FX
	graph::Graph::resetSchema();
#endif

}
void sim::runSimulation(const uint& startingNumAg, const uint& granularity) {
	using graph::Graph; using graph::node; 
	sInNode.resize(Graph::n);									// ----> Keeps track, for each node v, of how many susceptible agents are in v at time t.
	iInNode.resize(Graph::n);									// ----> Keeps track, for each node v, of how many infected agents are in v at time t.
	sAgents.resize(Graph::n);									// ----> Up-to-date list of susceptible agents within each node v.
	iAgents.resize(Graph::n);									// ----> Up-to-date list of infected agents within each node v.
	for (node v = 0; v < Graph::n; ++v) sAgents[v].resize((size_t)LIST_INI_SZ + 1);	// ----> The extra spot is to store the actual number of elements in the list, which may differ from the container size.
	for (node v = 0; v < Graph::n; ++v) iAgents[v].resize((size_t)LIST_INI_SZ + 1);	// ----> The extra spot is to store the actual number of elements in the list, which may differ from the container size.

	//Temp:
	//IavNumHops.resize(NUM_AGENTS, 0);
	//SavNumHops.resize(NUM_AGENTS, 0);
	//ITransitions.resize(NUM_AGENTS, 0);
	//STransitions.resize(NUM_AGENTS, 0);
	//
	//hopsUntilI.reserve((size_t)(T * nL));
	//hopsUntilS.reserve((size_t)(T * nL));
	//hopsUntilI.emplace_back(0);
	//hopsUntilS.emplace_back(0);

#ifdef PROTECTION_FX
	Graph::setProbs();
#endif

	uint _startingNumAg = ((startingNumAg != 0) && (startingNumAg < NUM_AGENTS)) ? startingNumAg : NUM_AGENTS;
	uint step = (_startingNumAg == NUM_AGENTS) ? 1 : (granularity == 0) ? 1 : granularity;
	const uint span = NUM_AGENTS - _startingNumAg + 1;
	uint numScenarios = (span / step) + 1 * (span % step);
	uint scenario = 0;
	uint _numAgents;

//#ifdef PER_BLOCK 
	vector<rtl> v_Iab(Graph::block_prob.size(), 0.0), v_ilb(Graph::block_prob.size(), 0.0), v_Sab(Graph::block_prob.size(), 0.0);
//#else
//	vector<rtl> v_Iv(Graph::n, 0.0), v_Sv(Graph::n, 0.0);
//#endif
	while (scenario < numScenarios){
		_numAgents = _startingNumAg + (scenario * step);
		
		//Agent-wise reset:
		Stats::resetAvDur();
		for (uint r = 0; r < ROUNDS; ++r) totalSimTime[r] = 0;
		I_0 = (I_0 > NUM_AGENTS) ? NUM_AGENTS : I_0;
		
		//Initializing applicable stats:
		Stats::initStream(streamType::avDuration);
#ifdef INFECTED_FRACTION
		Stats::initStream(streamType::infFrac);
#endif
		Reporter::startChronometer("\n\nRunning scenario " + std::to_string(scenario + 1) + "/" + std::to_string(numScenarios) + "...");

		//E[#ag] (q_b weighted):
		rtl Eag = 0;
		for (size_t b = 1; b < Graph::block_prob.size(); ++b) {
			Eag += NUM_AGENTS * Graph::q_b[b] * Graph::q_b[b];
		}		

		Reporter::simulationInfo(I_0, ROUNDS, T, NUM_AGENTS, TAU_aa, GAMMA_a, LAMBDA, Wi, Ws, Graph::averageDegree, Graph::_2ndMmt, Eag, Graph::bk, Graph::maxKbnb, Graph::block_prob);
		for (uint round = 0; round < ROUNDS; ++round) {
#ifdef MEASURE_ROUND_EXE_TIME
			//Reporter::startChronometer("\n  Round " + std::to_string(round + 1) + "...");
			Reporter::startChronometer("");
#else
			//Reporter::progress(round);
#endif
			//Round-wise reset:
			resetVariables();


			// * DISTRIBUTING THE AGENTS ACROSS THE NETWORK *
			//Random distribution of agents (initially susceptible):
			for (agent i = 0; i < _numAgents; ++i) {
				const node v = randomLCCNode();
				enterNodeAsSus(i, v, -EARLY_MOBILITY);
//#ifndef PER_BLOCK
//				++v_Sv[v];
//#endif
			}

			//EARLY MOBILITY:
			//Initially, we schedule a single 'walk' event for each agent. A new 'walk' job will then be created and scheduled for an agent the moment its current 'walk' job is processed.
			{
				std::cout << '\n' << "Applying early mobility... ";
				for (agent i = 0; i < _numAgents; ++i)
					schedule.emplace(i, EXPLambda() - EARLY_MOBILITY, action::walk);
				//Now, all agents are put to walk for a certain time window, so they can get somewhat mixed in the network:

#ifndef BYPASS_SIMULATION
				job j;
				rtl _now;
				do {
					nextJob(j, _now);
					quietWalk(j.target, _now);		// ----> A "quiet walk" does not produce any statistical data, as opposed to the 'walk()' method.
					schedule.emplace(j.target, _now + EXPLambda(), action::walk);
				} while (_now < TIME_ZERO);
				std::cout << "done.";
#endif //BYPASS_SIMULATION

				//Now we must define which agents are infected from the start, and then schedule their 'recover' event.
				//It is not a problem to always infect the first 'itotal' agents, since the starting node for each of them will be randomly set.
				for (agent ag = 0; ag < I_0; ++ag)
					infectAg(ag, TIME_ZERO);

			} //EARLY MOBILITY

			//Initializing the numerical solution:
			{
				const rtl ia = (rtl)I_0 / NUM_AGENTS, sa = (rtl)(NUM_AGENTS - I_0) / NUM_AGENTS;
				//Self-contained blocks:
				//for (uint b = (uint)graph::Graph::block_prob.size() - 1; b > 0; --b) {
				//	v_Iab[b] = ia * graph::Graph::Kb[b];
				//	v_Sab[b] = sa * graph::Graph::Kb[b];
				//}

				//Instance-wise:
				for (node v = 0; v < (uint)graph::Graph::n; ++v) {
					const size_t myBlock = graph::Graph::gs[v].size();
					v_Iab[myBlock] += iInNode[v];
					v_Sab[myBlock] += sInNode[v];
				}

#ifdef DEBUG
				rtl sum = 0.0;
				for (int b = (uint)graph::Graph::block_prob.size() - 1; b > 0; --b) {
					const rtl nb = (rtl)graph::Graph::n * graph::Graph::block_prob[b];
					sum += nb * (ivb[b] + svb[b]);
				}
				assertm(abs(sum - NUM_AGENTS) < epsilon, "The computed *expected number of S-/I-agents per node* is incorrect.");
#endif //DEBUG
			}
			
			// * MAIN LOOP *
			bool earlyStop = false;
			rtl& roundDuration = totalSimTime[round]; 
			rtl timeLimit = T;
			job j;
			nextJob(j, now);
			const rtl logInterval = 0.01 * T; // ----> Log progress to the console at one-percent increments.
			rtl prevLog = 0;
			std::cout << std::fixed;
			std::cout << std::setprecision(1);
			std::cout << '\n' << "Round " << (round + 1) << '/' << ROUNDS << ": 0% complete...     ";
#ifdef BYPASS_SIMULATION
			while (false) {
#else
			while (now < timeLimit) {
#endif
				roundDuration += (now - roundDuration);
				switch (j.a) {
				case action::walk:
					walk(j.target, now);
					schedule.emplace(j.target, now + EXPLambda(), action::walk);
					break;
				case action::agInfectAg:
					agFate_fromAg(j.target, now, j.infective, j.validity_S, j.validity_I, j.streamer);
					break;
				//case action::siteInfectAg:
				//	agFate_fromSite(j.target, now, j.infective, j.validity_S, j.validity_I);
				//	break;
				//case action::agInfectSite:
				//	siteFate(j.target, now, j.infective, j.validity_S, j.validity_I);
				//	break;
				//case action::recoverSite:
				//	recoverSite(j.target, now);
				//	break;
				default:
					recoverAg(j.target, now);
					if (iaTotal == 0) {
						earlyStop = true; 
						timeLimit = now; // ----> Exits the 'while'-loop by the next conditional test.
					}
					break;
				}
				nextJob(j, now);
				if (now - prevLog > logInterval) {
					prevLog = now;
					std::cout << '\r' << "Round " << (round + 1) << '/' << ROUNDS << ": " << std::round((now / T) * 100) << "% complete...";
				}
			}
			std::cout << '\r' << "Round " << (round + 1)  << '/' << ROUNDS << ": " << "100% complete!    ";
			Stats::partialsAvDur(roundDuration);
#ifdef MEASURE_ROUND_EXE_TIME
			Reporter::stopChronometer((earlyStop) ? "done-IV" : "done-TL");	// ----> TL == Time Limit; IV == Infection Vanished.
			Reporter::durationInfo(roundDuration);
#endif
#ifdef INFECTED_FRACTION
			Reporter::startChronometer(" Saving To File...");	
			Stats::iFracToFile(OVERLOOK);
			Stats::endStream(streamType::infFrac);
			Reporter::stopChronometer(" done");
#endif
		} // ** for (uint round = 0; round < ROUNDS; ++round)
		Stats::writeToFile(streamType::avDuration, Ws, Wi, _numAgents);
		Stats::endStream(streamType::avDuration);

#ifndef MEASURE_ROUND_EXE_TIME
		Reporter::tell("\nAll rounds completed.\n");
#endif
		Reporter::stopChronometer("Scenario " + std::to_string(scenario + 1) + "/" + std::to_string(numScenarios) + " completed");
		Reporter::avSimTimeInfo(Stats::avDuration());
		++scenario;
	} // ** while (scenario < numScenarios)

#ifdef SOLVE_NUMERICALLY
	//Runge-Kutta:
	constexpr uint outputGranularity = 50;
	constexpr rtl stepSize = 0.001;
	constexpr uint largerDetailUntil = 1000;
	constexpr rtl epsilon = 0.01 ;
	constexpr rtl timeIncrement = stepSize * outputGranularity;
	vector<rtl> saveToFile_diadt;
	vector<rtl> saveToFile_dildt;
	uint outputSize = 0;
	std::stringstream name;
	name << SHORT_LABEL 
		<< "_N"		<< N 
		<< "_AG"	<< NUM_AGENTS 
		<< "_T"		<< TAU_aa 
		<< "_G"		<< GAMMA_a 
		<< "_L"		<< LAMBDA 
		<< "_Wi"	<< Wi
		<< "_Ws"	<< Ws
		<< "_STime" << T 
		<< "_R"		<< ROUNDS;
		//<< "_Tal"	<< TAU_al 
		//<< "_Tla"	<< TAU_la 
		//<< "_Gl"	<< GAMMA_l 
	baseName = name.str();

	std::string method;
#ifdef MASTER	
	Reporter::startChronometer("\nSolving at network level...");
	Solver::rkMaster(0, v_Iab, v_Sab, T, stepSize, epsilon, saveToFile_diadt, saveToFile_dildt, outputSize, outputGranularity, largerDetailUntil);
	//Reporter::stopChronometer(" done", "execTime_master.txt");
	Reporter::stopChronometer(" done");
	method = "MASTER";
	
	//Saving to file:
	{
		std::ofstream RKdata;
		std::stringstream _ss;
		_ss.precision(4);
		_ss << "./stats/Runge-Kutta_" << method << "_" << baseName << ".csv";
		std::string _fileName = _ss.str();
		RKdata.open(_fileName);
		RKdata << "Time\tia\til\n"; // ----> Header
		rtl _time = 0;
		for (size_t s = 0; s < largerDetailUntil; ++s) {
			RKdata << _time << '\t' << saveToFile_diadt[s] << '\t' << saveToFile_diadt[s] << '\n';
			_time += stepSize;
		}
		for (size_t s = largerDetailUntil; s < outputSize; ++s) {
			RKdata << _time << '\t' << saveToFile_diadt[s] << '\t' << saveToFile_diadt[s] << '\n';
			_time += timeIncrement;
		}
		RKdata.close();
	}
#endif //MASTER	

#ifdef PER_BLOCK
	saveToFile_diadt.resize(0);
	saveToFile_dildt.resize(0);
	outputSize = 0;
	Reporter::startChronometer("\nSolving at block level...");
	Solver::rungeKutta4thOrder(0, v_Iab, v_Sab, v_ilb, T, stepSize, epsilon, saveToFile_diadt, saveToFile_dildt, outputSize, outputGranularity, largerDetailUntil);
	//Reporter::stopChronometer(" done", "execTime_block.txt");
	Reporter::stopChronometer(" done");
	method = "BLOCK";

	//Saving to file:
	{
		std::ofstream RKdata;
		std::stringstream _ss;
		_ss.precision(4);
		_ss << "./stats/Runge-Kutta_" << method << "_" << baseName << ".csv";
		std::string _fileName = _ss.str();
		RKdata.open(_fileName);
		RKdata << "Time\tia\til\n"; // ----> Header
		rtl _time = 0;
		for (size_t s = 0; s < largerDetailUntil; ++s) {
			RKdata << _time << '\t' << saveToFile_diadt[s] << '\t' << saveToFile_diadt[s] << '\n';
			_time += stepSize;
		}
		for (size_t s = largerDetailUntil; s < outputSize; ++s) {
			RKdata << _time << '\t' << saveToFile_diadt[s] << '\t' << saveToFile_diadt[s] << '\n';
			_time += timeIncrement;
		}
		RKdata.close();
	}
//#else
//	Solver::rungeKutta4thOrder(0, v_Iv, v_Sv, v_ilb, T, stepSize, epsilon, saveToFile_diadt, saveToFile_dildt, outputSize, outputGranularity, largerDetailUntil);
//	method = "NODE";
//#endif
#endif //PER_BLOCK
#endif //SOLVE_NUMERICALLY
}
