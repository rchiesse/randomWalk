#include "graph.h"
#include "reporter.h"
#include "stats.h"
#include <numeric>		//std::accumulate()
#include <iomanip>		//std::fixed & std::setprecision

namespace sim { // Simulator's namespace.

//Main structures
enum class action{walk, recoverAg, recoverSite, agInfectAg, agInfectSite, siteInfectAg };
struct job {
	uint target		= INT32_MAX;			// ----> Target ID (either an agent or a site).
	real time		= 0;
	action a		= action::walk;

	//Infection only:
	uint validity_S			= UINT_MAX;		// ----> S-agent's snapshot from the moment the event was generated. If the S-agent's snapshot by the time the infection event is processed does not match this one, then the event is ignored and the infection does not occur.
	uint validity_I			= UINT_MAX;		// ----> I-agent's snapshot from the moment the event was generated. If the I-agent's snapshot by the time the infection event is processed does not match this one, then the event is ignored and the infection does not occur.		
	uint infective			= UINT_MAX;		// ----> ID of the infective actor (either an agent or a site). It will be necessary to later retrieve its snapshot by the time the infection event is processed. If this does not match 'validity_I', then the event is ignored and the infection does not occur.
	uint streamer			= UINT_MAX;		// ----> ID of the agent whose queue this event come from. This information makes it possible to get the next event from the correct queue (as it may come either from the S-agent or the I-agent) and insert it into the main queue.

	job() {}
	job(const uint&& _target, const real&& _time, const action&& _action) : target(_target), time(_time), a(_action) {}
	job(const uint& _target, const real&& _time, const action&& _action) : target(_target), time(_time), a(_action) {}
	job(const uint&& _target, const real&& _time, const action&& _action, const uint&& _infective, const uint&& _snapshot_S, const uint&& _snapshot_I, const uint&& _streamer) : target(_target), time(_time), a(_action), infective(_infective), validity_S(_snapshot_S), validity_I(_snapshot_I), streamer(_streamer){}
	job(const uint& _target, const real& _time, const action&& _action, const uint& _infective, const uint& _snapshot_S, const uint& _snapshot_I, const uint& _streamer) : target(_target), time(_time), a(_action), infective(_infective), validity_S(_snapshot_S), validity_I(_snapshot_I), streamer(_streamer) {}

	//Sites:
	job(const uint&& _target, const real&& _time, const action&& _action, const uint&& _infective, const uint&& _snapshot_S, const uint&& _snapshot_I) : target(_target), time(_time), a(_action), infective(_infective), validity_S(_snapshot_S), validity_I(_snapshot_I) {}
	job(const uint& _target, const real& _time, const action&& _action, const uint& _infective, const uint& _snapshot_S, const uint& _snapshot_I) : target(_target), time(_time), a(_action), infective(_infective), validity_S(_snapshot_S), validity_I(_snapshot_I) {}

	job(const job& other) { *this = other; }						// ----> Copy Contructor. Jobs will be handled by an STL vector (via priority-queue). It thus requires a copy contructor to be explicitly defined. 

};
struct earlier {
	bool operator()(const job& j1, const job& j2) {
		return j1.time > j2.time;
	}
};

std::vector<std::priority_queue<job, std::vector<job>, earlier>> agJob(NUM_AGENTS);
std::priority_queue<job, std::vector<job>, earlier> schedule;


//Simulation variables:
uint saTotal;														// ----> Up-to-date number of SUSCEPTIBLE AGENTS during the simulation.
uint iaTotal;														// ----> Up-to-date number of INFECTED AGENTS during the simulation.
uint ilTotal = 0;													// ----> Up-to-date number of INFECTED SITES during the simulation.
real now;
std::vector<real> totalSimTime(ROUNDS, 0);							// ----> Total simulation time at each round, to average upon.

//Agent control variables
using std::vector; using graph::node;
vector<real> exposure		(NUM_AGENTS);							// ----> The total time interval a susceptible agent remained exposed to infected individuals at a given node. Susceptible agents have this value "zeroed" every time they enter a new node.
vector<uint> snapshot_a		(NUM_AGENTS,0);							// ----> A counter assigned to each AGENT and incremented every time the agent walks. Its main purpose is to control whether or not an infection event should be applied: when the infection event is processed, even if the agent comes to be currently located in the node at which the event relates, we must ensure that it is not the case that between the event creation and its processing the agent went somewhere else and then came back to said node. We do this by checking the event's snapshot against the agent's. These must be the same for the infection to take place.
vector<uint> snapshot_l		(N, 0);									// ----> A counter assigned to each LOCATION (NODE). Its purpose is the same of 'snapshot_a'.
vector<real> iniExposureTime(NUM_AGENTS);							// ----> The moment a new exposition cicle starts for a susceptible agent at its current node.
vector<node> currentNode	(NUM_AGENTS);							// ----> Keeps track of the node an agent is currently located.
vector<uint> indexWithinNode(NUM_AGENTS);							// ----> By storing the agent's index in the list of its current node, we are able to always find any agent in O(1). **This is critical for the overall performance**
vector<bool> isInfectedAg	(NUM_AGENTS);							// ----> Keeps track of each agent's current state.
vector<bool> isInfectedSite (N, false);								// ----> Keeps track of each site's current state.
#ifdef QUIET_INFECTION
vector<bool> quiet			(NUM_AGENTS, false);					// ----> Once an S-ag a is infected within the node v, this control variable governs whether a is allowed to further infect other S-agents inside v even before moving to somewhere else first. If set to 'false', then a will only become a propagator once it moves to another node first.
#endif
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

//std::string baseName;
//real Sim::beta_a;														// ----> Force of infection from an I-agent to an S-agent
//real Sim::beta_al;														// ----> Force of infection from an I-agent to a site
//real Sim::beta_la;														// ----> Force of infection from a site to an I-agent

/* PROTOTYPES */
//Returns a random value uniformly drawn from the interval [0, openRange), i.e. a number between 0 and 'openRange - 1' (for ex., "randomInt(6)" will return an integer between 0 and 5). This is particularly useful when determining an agent's "next stop" during its random walk. In this case, "randomInt(<v's degree>)" provides a neighbor's index each time an agent walks out some node v.
static const uint randomInt(const uint& openRange);
//Exponential random-number generator.
static const real EXP(const real& param);
//Exponential random-number generator based on the transmissibility parameter.
static const real EXPTau_aa();
static const real EXPTau_al();
static const real EXPTau_la();
//Exponential random-number generator that sets the moment an infected agent will recover.
static const real EXPGamma_a();
static const real EXPGamma_l();
//Exponential random-number generator that sets the next moment an agent walks.
static const real EXPLambda();
void resetVariables();
void readParams();
void runSimulation	(const uint& startingNumAg = 0, const uint& granularity = 5);
void walk		    (const agent& ag, const real& now);
void recoverAg	    (const agent& ag, const real& now);
void recoverSite    (const uint& site, const real& now);

using graph::node;
void agFate_fromAg	(const agent& ag, const real& now, const uint& infective, const uint& validity_S, const uint& validity_I, const uint& streamer);		// ----> Determines whether or not an exposed, susceptible agent 'ag' will become infected.
void agFate_fromSite(const agent& ag, const real& now, const uint& infective, const uint& validity_S, const uint& validity_I);		// ----> Determines whether or not an exposed, susceptible agent 'ag' will become infected.
void siteFate		(const uint& v, const real& now, const uint& infective, const uint& validity_L, const uint& validity_I);		// ----> Determines whether or not an exposed, susceptible agent 'ag' will become infected.
void enterNodeAsSus (const agent& ag, const node& v, const real& now);
void check_in		(const agent& ag, const node& v, vector<vector<agent>>& _where);
void check_out		(const agent& ag, const node& v, vector<vector<agent>>& _where);
void enterNodeAsInf (const agent& ag, const node& v, const real& now);
void leaveNodeAsSus (const agent& ag, const node& v, const real& now);
void leaveNodeAsInf (const agent& ag, const node& v, const real& now);
const node& nextNodeForSus(const node& _currNode);
const node& nextNodeForInf(const node& _currNode);
void nextJob(job& _j, real& _time);
void setBeta2ndMmt();
//void getSigma_a(const double& Ia, double& sigma_as, double& sigma_ai);
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
#endif
	sim::setBeta2ndMmt();
	sim::runSimulation(sim::STARTING_NUM_AG, sim::GRAN_NUM_AG);
	sim::Reporter::stopChronometer("\n\nSimulation completed");
	sim::Reporter::logTimestamp("End of simulation.");
	return 0;
}


/* IMPLEMENTATION */
static const uint sim::randomInt(const uint& openRange) { return (uint)(floor(openRange * U())); }
static const real sim::EXP(const real& param) { return -(1.0 / param) * log(U()); }
static const real sim::EXPTau_aa()	{ return NEG_RECIPR_TAU_aa	* log(U()); }
static const real sim::EXPTau_al()	{ return NEG_RECIPR_TAU_al	* log(U()); }
static const real sim::EXPTau_la()	{ return NEG_RECIPR_TAU_la	* log(U()); }
static const real sim::EXPGamma_a()	{ return NEG_RECIPR_GAMMA_a	* log(U()); }
static const real sim::EXPGamma_l()	{ return NEG_RECIPR_GAMMA_l	* log(U()); }
static const real sim::EXPLambda()	{ return NEG_RECIPR_LAMBDA	* log(U()); }

void sim::setBeta2ndMmt() {
#ifndef CLIQUE
	beta_a = (double)(LAMBDA * SIGMA_aa * graph::Graph::psi) / graph::Graph::averageDegree;
	beta_al = (LAMBDA * NUM_AGENTS *  SIGMA_al)/N;
	beta_la = (LAMBDA * graph::Graph::_2ndMmt * N * SIGMA_la) / (N * pow(graph::Graph::averageDegree, 2));
#endif
	//Normalization
	//nT = TAU_aa;
	//nL = LAMBDA;
	//nG = GAMMA_a;
	if (TAU_aa <= LAMBDA) {
		nT = 1.0;
		nL = LAMBDA / TAU_aa;
		nG = GAMMA_a / TAU_aa;
	}
	else {
		nL = 1.0;
		nT = TAU_aa / LAMBDA;
		nG = GAMMA_a / LAMBDA;
	}
}

//void sim::getSigma_a(const double& Ia, double& sigma_as, double& sigma_ai) {
//	sigma_as = 0;
//	sigma_ai = 0;
//	for (uint b = 1; b < graph::Graph::block_prob.size(); ++b) {
//		if (graph::Graph::kb[b] == 0)
//			continue;
//		const double expNumAg = graph::Graph::kb[b];
//		const double expNumInfAg = expNumAg * Ia;			// ----> Talvez seja errado fazer dessa forma...
//		const double minInfAg = std::min(1.0, expNumInfAg);	// ----> REVER! Talvesz o mínimo unitário faça mais sentido (fora que evita erros de número muito pequeno em cenários extremamente esparsos, onde esse número seria mt próx de zero).
//		const double expNumSusAg = expNumAg * (1.0 - Ia);
//		constexpr double euler = 0.5772156649;
//		const double digamma_i = log(expNumInfAg) - 1.0 / (2 * expNumInfAg);
//		const double digamma_s = log(expNumSusAg) - 1.0 / (2 * expNumSusAg);
//		const double H_i = std::max(1.0, euler + digamma_i);
//		const double H_s = std::max(1.0, euler + digamma_s);
//		const double prob_NoAcq = (H_i * LAMBDA + H_i * GAMMA_a + LAMBDA) / (H_i * LAMBDA + LAMBDA + H_i * GAMMA_a + expNumInfAg * TAU_aa);
//		const double prob_acq = (std::max(expNumInfAg, 1.0) * TAU_aa) / ((1.0 + H_i) * LAMBDA + std::max(expNumInfAg, 1.0) * TAU_aa);
//		const double prob_inf = (TAU_aa) / (2 * LAMBDA + std::max(expNumInfAg, 1.0) * TAU_aa);
//		const double prob_NoTransmission = (H_s * LAMBDA + LAMBDA + GAMMA_a) / (H_s * LAMBDA + LAMBDA + GAMMA_a + TAU_aa);
//		sigma_as += prob_acq * graph::Graph::block_prob[b];
//		sigma_ai += prob_inf * graph::Graph::block_prob[b];
//	}
//
//}

#ifdef SOLVE_NUMERICALLY



#ifdef CLIQUE
real sim::diabdt(const real& Ia, const real& Sa) {
	const double pb = 1.0;
	const double nb = graph::Graph::n;
	const double qb = 1.0;
	const double ibnb = Ia/nb;
	const double sbnb = Sa/nb;
	const double kbnb = ibnb + sbnb;
	//const double Sa = (double)NUM_AGENTS - Ia;

	//RONALD (BEST SO FAR):
	if (Sa == 0.0)
		return - (GAMMA_a * Ia);

	//long double H = EULER * log(sbnb + 1.0) / (2.0 * nL);
	long double H = EULER * log(std::min(sbnb, ibnb) * ibnb * (1.0 / nT)) / (2.0 * nL);
	long double ii = ibnb + 1.0;
	const double prob_inf = ii / (H + ii);
	//const double p = std::min(sbnb, 1.0) * std::min(ibnb, 1.0);
	return 
		Ia * sbnb * prob_inf * nT 
		- (nG * Ia);
}

real sim::dsabdt(const real& Ia, const real& Sa) {
	const double pb = 1.0;
	const double nb = graph::Graph::n;
	const double qb = 1.0;
	const double ibnb = Ia / nb;
	const double sbnb = Sa / nb;
	const double kbnb = ibnb + sbnb;
	//const double Sa = (double)NUM_AGENTS - Ia;

	//BEST SO FAR:
	if (Sa == 0.0) {
		return 
			GAMMA_a * Ia;
	}
	//long double H = EULER * log(sbnb + ibnb + 1.0) / (2.0 * nL);
	long double H = EULER * log(std::min(sbnb, ibnb) * ibnb * (1.0 / nT)) / (2.0 * nL);
	long double ii = ibnb + 1.0;
	const double prob_inf = ii / (H + ii);
	//const double p = std::min(sbnb, 1.0) * std::min(ibnb, 1.0);
	return 
		- Ia * sbnb * prob_inf * nT
		+ (nG * Ia);
}

void sim::step(const real& h, real& Ia, real& Sa) {
	constexpr real one_sixth = 1.0 / 6.0;
	vector<real> k1(2, 0), k2(2, 0), k3(2, 0), k4(2, 0);

	lookAhead(h, Ia, Sa, k1);
	lookAhead(h, Ia, Sa, k2, k1, 0.5);
	lookAhead(h, Ia, Sa, k3, k2, 0.5);
	lookAhead(h, Ia, Sa, k4, k3);

	//Take step:
	Ia += one_sixth * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
	Sa += one_sixth * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
}

void sim::lookAhead(const real& h, real& Ia, real& Sa, std::vector<real>& target) {
	target[0] = h * diabdt(Ia, Sa);
	target[1] = h * dsabdt(Ia, Sa);
}

void sim::lookAhead(const real& h, real& Ia, real& Sa, std::vector<real>& target, std::vector<real>& base, const double& fraction) {
	target[0] = h * diabdt(Ia + fraction * base[0], Sa + fraction * base[1]);
	target[1] = h * dsabdt(Ia + fraction * base[0], Sa + fraction * base[1]);
}

void sim::rungeKutta4thOrder(const real& t0, real& Ia, real& Sa, const real& t, const real& h, const real& epsilon, vector<real>& saveToFile_diadt, vector<real>& saveToFile_dildt, uint& outputSize, const uint& outputGranularity, const real& largerDetailUntil) {
	uint totalSteps = (uint)((t - t0) / h) + 1;
	saveToFile_diadt.resize((uint64_t)largerDetailUntil + (totalSteps - ((uint)largerDetailUntil) / outputGranularity) + 1);

	saveToFile_diadt[0] = Ia / NUM_AGENTS;
	bool end = false;
	++outputSize;

	//For the first 'largerDetailUntil' iterations every step is stored in a vector ('saveToFile'), for later being written to file.
	for (uint s = 1; s < largerDetailUntil; ++s) {
		real prevIA = Ia;
		step(h, Ia, Sa);
		if (Ia < epsilon) {
			saveToFile_diadt[outputSize] = 0;
			end = true;
			++outputSize;
			break;
		}
		if (Ia == prevIA) {
			while (s < largerDetailUntil) {
				saveToFile_diadt[outputSize] = Ia / NUM_AGENTS;
				++s;
				++outputSize;
			}
			end = true;
			break;
		}
		saveToFile_diadt[outputSize] = Ia / NUM_AGENTS;
		++outputSize;
	}
	if (end) return;

	//From the 'largerDetailUntil' iteration on, we afford to ignore 'outputGranularity'-size windows of values, so that the saved file does not grow explosively.
	for (uint s = (uint)largerDetailUntil; s < totalSteps; ++s) {
		real prevIA = Ia;
		step(h, Ia, Sa);
		if (Ia < epsilon) {
			saveToFile_diadt[outputSize] = 0;
			++outputSize;
			break;
		}
		if (Ia == prevIA) {
			while (s < totalSteps) {
				if (s % outputGranularity == 0) {
					saveToFile_diadt[outputSize] = Ia / NUM_AGENTS;
					++outputSize;
				}
				++s;
			}
			break;
		}
		if (s % outputGranularity == 0) {
			saveToFile_diadt[outputSize] = Ia / NUM_AGENTS;
			++outputSize;
		}
	}
}

#else //CLIQUE
#ifdef PER_BLOCK
real sim::diabdt(const real& Ia, const real& Iab, const real& Sab, const uint& b) {
	const double& pb = graph::Graph::block_prob[b];
	const double nb = graph::Graph::n * pb;
	const double& qb = graph::Graph::q_b[b];
	const double ibnb = Iab / nb;
	const double sbnb = Sab / nb;
	const double kbnb = ibnb + sbnb;
	const double Sa = (double)NUM_AGENTS - Ia;
	const double l = ((double)b - 1) / ((double)b);

	//RONALD (BEST SO FAR):
	if (Sab == 0.0) {
		return Ia * LAMBDA * qb - Iab * LAMBDA
			- (GAMMA_a * Iab);
	}

	long double H = EULER * log(sbnb + 1.0) / (2.0 * nL);
	long double ii = ibnb + 1.0;
	const double prob_inf = ii / (H + ii);
	return Ia * nL * qb - Iab * nL
		+ Iab * sbnb * prob_inf * nT
		- (nG * Iab);
}

real sim::dsabdt(const real& Ia, const real& Iab, const real& Sab, const uint& b) {
	const double& pb = graph::Graph::block_prob[b];
	const double nb = graph::Graph::n * pb;
	const double& qb = graph::Graph::q_b[b];
	const double ibnb = Iab / nb;
	const double sbnb = Sab / nb;
	const double kbnb = ibnb + sbnb;
	const double Sa = (double)NUM_AGENTS - Ia;
	const double l = ((double)b - 1) / ((double)b);

	//DON 2:
	//return nb * (
	//	(LAMBDA/nb) * (Sa * qb - Sab) 
	//	- LAMBDA * sbnb * ibnb * prob_inf
	//	- LAMBDA * sbnb * ibnb * prob_acq
	//	+ GAMMA_a * ibnb
	//	);

	//Ronald (sparse):
	//return (Sa - Sab) * LAMBDA * qb - Sab * LAMBDA * (1.0 - qb)
	//	- 2 * ((Sab * Iab) / nb) * LAMBDA * SIGMA_aa
	//	+ (GAMMA_a * Iab);

	////RONALD (excelente quando LAMBDA == TAU, mesmo em cenários densos):
	//const double C = TAU_aa - LAMBDA + 1.0;
	//C -= TAU_aa / ((2.0 * LAMBDA) / TAU_aa);
	//const double prob_inf = TAU_aa / (2 * LAMBDA + std::max(ibnb, 1.0) * TAU_aa);
	//const double prob_inf_2nd = TAU_aa / (2 * LAMBDA + ((Iab + Sab) / nb) * TAU_aa);
	//double C = (LAMBDA/TAU_aa);
	//const double prob_inf = 1.0 / ((1.0 + std::min(sbnb, 1.0)) + std::max(ibnb, 1.0));
	//double prob_inf_2nd = 1.0 / (2.0 + ((Iab + Sab) / nb));
	//prob_inf_2nd *= std::max(0.0, sbnb - ibnb) * prob_inf;
	//prob_inf_2nd *= (Ia > Sa) ? 1.0 : -1.0;
	//const double prob_acq = ibnb * prob_inf + (1.0 - (ibnb * prob_inf)) * prob_inf_2nd * (std::max(0.0, sbnb) * prob_inf);
	//if (ibnb < 1.0 || sbnb < 1.0) {
	//	return (Sa - Sab) * LAMBDA * qb - Sab * LAMBDA * (1.0 - qb)
	//		- 2 * ((Sab * Iab) / nb) * LAMBDA * SIGMA_aa
	//		+ (GAMMA_a * Iab);
	//}

	//RONALD CANDIDATO
	//const double lt = (LAMBDA > TAU_aa) ? 0.0 : 1.0;	// ----> lt is a flag for "LAMBDA larger than TAU".
	//const double prob_inf = 1.0 / (1.0 + lt + ibnb);
	////const double prob_inf = 1.0 / (2.0 + ibnb);
	//const double prob_acq = (ibnb * prob_inf);
	//const double s_ag = (LAMBDA > TAU_aa) ? sbnb : std::min(ibnb, sbnb);
	//return 
	//	Sa * LAMBDA * qb - Sab * LAMBDA 
	//	- Iab * s_ag * prob_inf * TAU_aa
	//	- Sab * ibnb * prob_acq * TAU_aa
	//	+ (GAMMA_a * Iab);


	//BEST SO FAR:
	if (Sab == 0.0) {
		return Sa * LAMBDA * qb - Sab * LAMBDA
			+ (GAMMA_a * Iab);
	}
	//const double s_ag = (LAMBDA > TAU_aa) ? sbnb : std::min(ibnb, sbnb);
	long double H = EULER * log(sbnb + 1.0) / (2.0 * nL);
	long double ii = ibnb + 1.0;
	const double prob_inf = ii / (H + ii);
	return Sa * nL * qb - Sab * nL
		- Iab * sbnb * prob_inf * nT
		+ (nG * Iab);

	//DON : Good for dense scenarios or large walk rate
	//if (kb == 0.0)
	//	return Sa * LAMBDA * qb - Sab * LAMBDA;
	//return 
	//	Sa * LAMBDA * qb - Sab * LAMBDA
	//	- Iab * sbnb * TAU_aa
	//	+ (GAMMA_a * Iab);
}

void sim::step(const real& h, real& Ia, std::vector<real>& v_Iab, std::vector<real>& v_Sab) {
	constexpr real one_sixth = 1.0 / 6.0;
	const uint blocks = static_cast<uint>(graph::Graph::block_prob.size());
	vector<real> k1(2 * blocks, 0), k2(2 * blocks, 0), k3(2 * blocks, 0), k4(2 * blocks, 0);

	lookAhead(h, Ia, v_Iab, v_Sab, k1);
	lookAhead(h, Ia, v_Iab, v_Sab, k2, k1, 0.5);
	lookAhead(h, Ia, v_Iab, v_Sab, k3, k2, 0.5);
	lookAhead(h, Ia, v_Iab, v_Sab, k4, k3);

	//Take step:
	for (uint b = (uint)graph::Graph::block_prob.size() - 1; b > 0; --b) {
		if (graph::Graph::block_prob[b] == 0)
			continue;

		v_Iab[b] += one_sixth * (k1[2 * b] + 2 * k2[2 * b] + 2 * k3[2 * b] + k4[2 * b]);
		v_Sab[b] += one_sixth * (k1[2 * b + 1] + 2 * k2[2 * b + 1] + 2 * k3[2 * b + 1] + k4[2 * b + 1]);
	}
}

void sim::lookAhead(const real& h, real& Ia, const std::vector<real>& v_Iab, const std::vector<real>& v_Sab, std::vector<real>& target) {
	update_Ia(Ia, v_Iab);
	for (uint b = (uint)graph::Graph::block_prob.size() - 1; b > 0; --b) {
		if (graph::Graph::block_prob[b] == 0)
			continue;

		target[2 * b] = h * diabdt(Ia, v_Iab[b], v_Sab[b], b);
		target[2 * b + 1] = h * dsabdt(Ia, v_Iab[b], v_Sab[b], b);
	}
}

void sim::lookAhead(const real& h, real& Ia, const std::vector<real>& v_Iab, const std::vector<real>& v_Sab, std::vector<real>& target, std::vector<real>& base, const double& fraction) {
	update_Ia(Ia, v_Iab, base, fraction);
	for (uint b = (uint)graph::Graph::block_prob.size() - 1; b > 0; --b) {
		if (graph::Graph::block_prob[b] == 0)
			continue;

		target[2 * b] = h * diabdt(Ia, v_Iab[b] + fraction * base[2 * b], v_Sab[b] + fraction * base[2 * b + 1], b);
		target[2 * b + 1] = h * dsabdt(Ia, v_Iab[b] + fraction * base[2 * b], v_Sab[b] + fraction * base[2 * b + 1], b);
	}
}

void sim::update_Ia(real& Ia, const std::vector<real>& v_Iab) {
	Ia = 0;
	for (uint b = (uint)graph::Graph::block_prob.size() - 1; b > 0; --b)
		Ia += v_Iab[b];
}

void sim::update_Ia(real& Ia, const std::vector<real>& v_Iab, const std::vector<real>& base, const double& fraction) {
	Ia = 0;
	for (uint b = (uint)graph::Graph::block_prob.size() - 1; b > 0; --b)
		Ia += v_Iab[b] + (fraction * base[2 * b]);
}
#else //PER_BLOCK
real sim::divbdt(const real& Ia, const real& Iv, const real& Sv, const uint& b) {
	const double& pb = graph::Graph::block_prob[b];
	const double nb = graph::Graph::n * pb;
	const double& qb = graph::Graph::q_b[b];
	const double qbnb = qb / nb;
	const double Sa = (double)NUM_AGENTS - Ia;
	const double prob_inf = (TAU_aa / (2.0 * LAMBDA + TAU_aa));
	return ((Ia - Iv) * LAMBDA * qbnb - Iv * LAMBDA * (1.0 - qbnb))
		+ (Sa - Sv) * LAMBDA * qbnb * Iv * prob_inf
		+ (Ia - Iv) * LAMBDA * qbnb * Sv * prob_inf
		- (GAMMA_a * Iv);
}
real sim::dsvbdt(const real& Ia, const real& Iv, const real& Sv, const uint& b) {
	const double& pb = graph::Graph::block_prob[b];
	const double nb = graph::Graph::n * pb;
	const double& qb = graph::Graph::q_b[b]; const double qbnb = qb / nb;
	const double Sa = (double)NUM_AGENTS - Ia;
	const double prob_inf = (TAU_aa / (2.0 * LAMBDA + TAU_aa));
	return ((Sa - Sv) * LAMBDA * qbnb - Sv * LAMBDA * (1.0 - qbnb))
		- (Sa - Sv) * LAMBDA * qbnb * Iv * prob_inf
		- (Ia - Iv) * LAMBDA * qbnb * Sv * prob_inf
		+ (GAMMA_a * Iv);
}
void sim::step(const real& h, real& Ia, std::vector<real>& v_Iv, std::vector<real>& v_Sv) {
	constexpr real one_sixth = 1.0 / 6.0;
	const uint blocks = static_cast<uint>(graph::Graph::block_prob.size());
	vector<real> k1(2 * graph::Graph::n, 0), k2(2 * graph::Graph::n, 0), k3(2 * graph::Graph::n, 0), k4(2 * graph::Graph::n, 0);

	lookAhead(h, Ia, v_Iv, v_Sv, k1);
	lookAhead(h, Ia, v_Iv, v_Sv, k2, k1, 0.5);
	lookAhead(h, Ia, v_Iv, v_Sv, k3, k2, 0.5);
	lookAhead(h, Ia, v_Iv, v_Sv, k4, k3);

	//Take step:
	for (int v = (uint)graph::Graph::n - 1; v >= 0; --v) {
		//if (graph::Graph::block_prob[b] == 0)
		//	continue;

		v_Iv[v] += one_sixth * (k1[2 * v] + 2 * k2[2 * v] + 2 * k3[2 * v] + k4[2 * v]);
		v_Sv[v] += one_sixth * (k1[2 * v + 1] + 2 * k2[2 * v + 1] + 2 * k3[2 * v + 1] + k4[2 * v + 1]);
	}
}
void sim::lookAhead(const real& h, real& Ia, const std::vector<real>& v_Iv, const std::vector<real>& v_Sv, std::vector<real>& target) {
	update_Ia(Ia, v_Iv);
	for (int v = (uint)graph::Graph::n - 1; v >= 0; --v) {
		target[2 * v] = h * divbdt(Ia, v_Iv[v], v_Sv[v], (uint)graph::Graph::g[v].size());
		target[2 * v + 1] = h * dsvbdt(Ia, v_Iv[v], v_Sv[v], (uint)graph::Graph::g[v].size());
	}
}
void sim::lookAhead(const real& h, real& Ia, const std::vector<real>& v_Iv, const std::vector<real>& v_Sv, std::vector<real>& target, std::vector<real>& base, const double& fraction) {
	update_Ia(Ia, v_Iv, base, fraction);
	for (int v = (uint)graph::Graph::n - 1; v >= 0; --v) {
		target[2 * v] = h * divbdt(Ia, v_Iv[v] + fraction * base[2 * v], v_Sv[v] + fraction * base[2 * v + 1], (uint)graph::Graph::g[v].size());
		target[2 * v + 1] = h * dsvbdt(Ia, v_Iv[v] + fraction * base[2 * v], v_Sv[v] + fraction * base[2 * v + 1], (uint)graph::Graph::g[v].size());
	}
}
void sim::update_Ia(real& Ia, const std::vector<real>& v_Iv) {
	Ia = 0;
	for (int v = (uint)graph::Graph::n - 1; v >= 0; --v)
		Ia += v_Iv[v];
}
void sim::update_Ia(real& Ia, const std::vector<real>& v_Iv, const std::vector<real>& base, const double& fraction) {
	Ia = 0;
	for (int v = (uint)graph::Graph::n - 1; v >= 0; --v)
		Ia += v_Iv[v] + (fraction * base[2 * v]);
}
#endif //PER_BLOCK

real sim::dilbdt(const real& Ia, const real& il, const real& Iab, const real& ilb, const uint& b) {
	return 0;
}

void sim::rungeKutta4thOrder(const real& t0, std::vector<real>& v_Iab, std::vector<real>& v_Sab, std::vector<real>& v_ilb, const real& t, const real& h, const real& epsilon, vector<real>& saveToFile_diadt, vector<real>& saveToFile_dildt, uint& outputSize, const uint& outputGranularity, const real& largerDetailUntil) {
	uint totalSteps = (uint)((t - t0) / h) + 1;
	saveToFile_diadt.resize((uint64_t)largerDetailUntil + (totalSteps - ((uint)largerDetailUntil) / outputGranularity) + 1);
	saveToFile_dildt.resize((uint64_t)largerDetailUntil + (totalSteps - ((uint)largerDetailUntil) / outputGranularity) + 1);

	double Ia = 0;
	for (uint b = (uint)v_Iab.size() - 1; b > 0; --b)
		Ia += v_Iab[b];
	real il = 0.0;

	saveToFile_diadt[0] = Ia / NUM_AGENTS;
	saveToFile_dildt[0] = il;
	bool end = false;
	++outputSize;

	//For the first 'largerDetailUntil' iterations every step is stored in a vector ('saveToFile'), for later being written to file.
	//bool stationary = false;
	for (uint s = 1; s < largerDetailUntil; ++s) {
		real prevIA = Ia;
		step(h, Ia, v_Iab, v_Sab);
		if (Ia < epsilon) {
			saveToFile_diadt[outputSize] = 0;
			end = true;
			++outputSize;
			break;
		}
		if (Ia == prevIA) {
			while (s < largerDetailUntil) {
				saveToFile_diadt[outputSize] = Ia / NUM_AGENTS;
				++s;
				++outputSize;
			}
			end = true;
			break;
		}
		saveToFile_diadt[outputSize] = Ia / NUM_AGENTS;

#ifdef NORM_SITE_PER_AG
		saveToFile_dildt[outputSize] = il * (N / NUM_AGENTS);
#else
		saveToFile_dildt[s] = il;
#endif
		++outputSize;
	}
	if (end) return;

	//From the 'largerDetailUntil' iteration on, we afford to ignore 'outputGranularity'-size windows of values, so that the saved file does not grow explosively.
	for (uint s = (uint)largerDetailUntil; s < totalSteps; ++s) {
		real prevIA = Ia;
		step(h, Ia, v_Iab, v_Sab);
		if (Ia < epsilon) {
			saveToFile_diadt[outputSize] = 0;
#ifdef NORM_SITE_PER_AG
			saveToFile_dildt[outputSize] = il * (N / NUM_AGENTS);
#else
			saveToFile_dildt[outputSize] = il;
#endif
			++outputSize;
			break;
		}
		if (Ia == prevIA) {
			while (s < totalSteps) {
				if (s % outputGranularity == 0) {
					saveToFile_diadt[outputSize] = Ia / NUM_AGENTS;
					++outputSize;
				}
				++s;
			}
			break;
		}
		if (s % outputGranularity == 0) {
			saveToFile_diadt[outputSize] = Ia / NUM_AGENTS;
#ifdef NORM_SITE_PER_AG
			saveToFile_dildt[outputSize] = il * (N / NUM_AGENTS);
#else
			saveToFile_dildt[outputSize] = il;
#endif
			++outputSize;
		}
	}
}
#endif //CLIQUE






//real sim::diadt(const real& Ia, const real& il) {
real sim::diadt(const real& Ia, const double& sumSbIb) {
	//return graph::Graph::psi * TAU_aa * sumSbIb - GAMMA_a * Ia;
	return - GAMMA_a * Ia;
}
real sim::dildt(const real& Ia, const real& il) {
	return  beta_al * (1.0 - il) * Ia - GAMMA_l * il;
}

#endif //SOLVE_NUMERICALLY

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
	if (lastPos > 1) {	// ----> "Is it still going to remain any susceptible agent at 'v' once 'ag' has been removed?". If true, then the agent from the last position is copied to ag's. There's no problem if the outcoming agent happens to be the one at the 'top' position. In this case, both inner instructions become redundant, not wrong.
		list[indexWithinNode[ag]] = list[lastPos];
		indexWithinNode[list[lastPos]] = indexWithinNode[ag];
	}
	--lastPos;
}
void sim::enterNodeAsSus (const agent& ag, const node& v, const real& now) {
	++snapshot_a[ag];
	check_in(ag, v, sAgents);
	const uint& numI = iInNode[v];					// ----> Number of infected agents currently hosted in v.
	uint&		numS = sInNode[v];					// ----> Number of susceptible agents currently hosted in v.
	++numS;
	
	if (numI > 0) {
		const vector<uint>& list = iAgents[v];
		real delta;
		for (uint i = numI; i > 0; --i) {
			delta = EXPTau_aa();
			agJob[ag].emplace(ag, now + delta, action::agInfectAg, list[i], snapshot_a[ag], snapshot_a[list[i]], ag);
		}

		//The most recent event is copied to the main queue and deleted from the agent's:
		schedule.emplace(agJob[ag].top());
		agJob[ag].pop();
	}
	//if (isInfectedSite[v]) {
	//	const real delta = EXPTau_la();
	//	schedule.emplace(ag, now + delta, action::siteInfectAg, v, snapshot_a[ag], snapshot_l[v]);
	//}
#ifdef PROTECTION_FX
	if (numS == 1)
		graph::Graph::updateHasS(v);
#endif
}
void sim::enterNodeAsInf (const agent& ag, const node& v, const real& now) {
	++snapshot_a[ag];
	check_in(ag, v, iAgents);
	const uint& numS = sInNode[v];				// ----> Number of susceptible agents currently hosted in v.
	uint&		numI = iInNode[v];				// ----> Number of infected agents currently hosted in v.
	++numI;

	if (numS > 0) {
		const vector<uint>& list = sAgents[v];
		real delta;
		for (uint i = numS; i > 0; --i) {
			delta = EXPTau_aa();
			agJob[ag].emplace(list[i], now + delta, action::agInfectAg, ag, snapshot_a[list[i]], snapshot_a[ag], ag);
		}

		//The most recent event is copied to the main queue and deleted from the agent's:
		schedule.emplace(agJob[ag].top());
		agJob[ag].pop();
	}
	//if (!isInfectedSite[v]) {
	//	const real delta = EXPTau_al();
	//	schedule.emplace(v, now + delta, action::agInfectSite, ag, snapshot_l[v], snapshot_a[ag]);
	//}
}
void sim::leaveNodeAsInf (const agent& ag, const node& v, const real& now) {
	++snapshot_a[ag];
	agJob[ag] = {};		// ----> Any events remaining into ag's queue becomes invalid, so that we may simply erase them.
	check_out(ag, v, iAgents);
	uint&		numI = iInNode[v];				// ----> Number of infected agents currently hosted in v.
	--numI;
}
void sim::leaveNodeAsSus (const agent& ag, const node& v, const real& now) {
	++snapshot_a[ag];
	agJob[ag] = {};		// ----> Any events remaining into ag's queue becomes invalid, so that we may simply erase them.
	check_out(ag, v, sAgents);
	uint&		numS = sInNode[v];				// ----> Number of susceptible agents currently hosted in v.
	--numS;
}
void sim::walk(const agent& ag, const real& now) {
	node v = (isInfectedAg[ag]) ? nextNodeForInf(currentNode[ag]) : nextNodeForSus(currentNode[ag]);
	if (v == currentNode[ag])
		return;
	if (isInfectedAg[ag]) {
		leaveNodeAsInf(ag, currentNode[ag], now);
		enterNodeAsInf(ag, v, now);
	} else {
		leaveNodeAsSus(ag, currentNode[ag], now);
		enterNodeAsSus(ag, v, now);
	}
}
void sim::recoverAg		 (const agent& ag, const real& now) {
	isInfectedAg[ag] = false;
	++saTotal;
	--iaTotal;
#ifdef INFECTED_FRACTION
	stat::Stats::bufferizeIFrac(ag, now, "Ra", iaTotal, ilTotal, NUM_AGENTS, OVERLOOK);
#endif
	const graph::node& v = currentNode[ag];
	leaveNodeAsInf(ag, v, now);
	enterNodeAsSus(ag, v, now);
}

void sim::recoverSite	(const uint& v, const real& now) {
	isInfectedSite[v] = false;
	--ilTotal;
#ifdef INFECTED_FRACTION
	stat::Stats::bufferizeIFrac(v, now, "Rl", iaTotal, ilTotal, NUM_AGENTS, OVERLOOK);
#endif
	++snapshot_l[v];
}

void sim::agFate_fromAg(const agent& ag, const real& now, const uint& infective, const uint& validity_S, const uint& validity_I, const uint& streamer) {
	if (!agJob[streamer].empty()) {
		schedule.emplace(agJob[streamer].top());
		agJob[streamer].pop();
	}
	if (validity_S != snapshot_a[ag] || validity_I != snapshot_a[infective])
		return;	// ----> Event became obsolete.
	
	isInfectedAg[ag] = true;
	--saTotal;
	++iaTotal;
#ifdef INFECTED_FRACTION
	stat::Stats::bufferizeIFrac(ag, now, "Ia", iaTotal, ilTotal, NUM_AGENTS, OVERLOOK);
#endif
	leaveNodeAsSus(ag, currentNode[ag], now);
	enterNodeAsInf(ag, currentNode[ag], now);
	schedule.emplace(ag, now + EXPGamma_a(), action::recoverAg); // ----> 'Recover' event is scheduled.
}
void sim::agFate_fromSite(const agent& ag, const real& now, const uint& infective, const uint& validity_S, const uint& validity_I) {
	if (validity_S != snapshot_a[ag] || validity_I != snapshot_l[infective]) 
		return;	// ----> Event became obsolete.
	
	isInfectedAg[ag] = true;
	--saTotal;
	++iaTotal;
#ifdef INFECTED_FRACTION
	stat::Stats::bufferizeIFrac(ag, now, "Ia", iaTotal, ilTotal, NUM_AGENTS, OVERLOOK);
#endif
	leaveNodeAsSus(ag, currentNode[ag], now);
	enterNodeAsInf(ag, currentNode[ag], now);
	schedule.emplace(ag, now + EXPGamma_a(), action::recoverAg); // ----> 'Recover' event is scheduled.
}
void sim::siteFate(const uint& v, const real& now, const uint& infective, const uint& validity_S, const uint& validity_I) {
	if (validity_S != snapshot_l[v] || validity_I != snapshot_a[infective]) 
		return;	// ----> Event became obsolete.
	
	isInfectedSite[v] = true;
	++ilTotal;
#ifdef INFECTED_FRACTION
	stat::Stats::bufferizeIFrac(v, now, "Il", iaTotal, ilTotal, NUM_AGENTS, OVERLOOK);
#endif
	schedule.emplace(v, now + EXPGamma_l(), action::recoverSite); // ----> 'Recover' event is scheduled.
	const uint& numS = sInNode[v];				// ----> Number of susceptible agents currently hosted in v.
	const vector<uint>& list = sAgents[v];
	for (uint i = numS; i > 0; --i) {
		const real delta = EXPTau_la();
		schedule.emplace(list[i], now + delta, action::siteInfectAg, v, snapshot_a[list[i]], snapshot_l[v]);
	}
}
const graph::node& sim::nextNodeForSus(const node& _currNode) {
	using graph::Graph;
#ifdef CLIQUE
	return Graph::g[randomInt(graph::Graph::n)];
#else
	return Graph::g[_currNode][randomInt((uint)Graph::g[_currNode].size())];
#endif //CLIQUE
}
const graph::node& sim::nextNodeForInf(const node& _currNode) {
	using graph::Graph;
#ifdef CLIQUE
	return Graph::g[randomInt(Graph::n)];
#else
	return Graph::g[_currNode][randomInt((uint)Graph::g[_currNode].size())];
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
	iaTotal = (I_0 > NUM_AGENTS) ? NUM_AGENTS : I_0;
	saTotal = NUM_AGENTS - iaTotal;
	now = 0;
	while (!schedule.empty()) schedule.pop();	// ----> Needed from the 2nd round on.
	for (agent a = 0; a < isInfectedAg.size(); ++a) isInfectedAg[a] = false;
	snapshot_a.resize(NUM_AGENTS, 0);
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
	Graph::g.resize(Graph::n);
	Graph::averageDegree = N;
	Graph::largestDegree = N;
	Graph::lccSize = N;
	Graph::m = N * N / 2;	// ----> Note that the total number of edges here is NOT (n*(n-1))/2 because each node contains an implicit self loop.
	Graph::selfLoops = N;
	Graph::smallestDegree = N;
	Graph::validBlocks = 1;
	sim::Reporter::networkInfo(Graph::n, Graph::m, Graph::averageDegree, Graph::largestDegree, Graph::smallestDegree, Graph::lccSize);
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

#ifdef CLIQUE
#ifdef PER_BLOCK
	real Ia = 0.0, Sa = 0.0;
#else
	vector<real> v_Iv(Graph::n, 0.0), v_Sv(Graph::n, 0.0);
#endif
#else //CLIQUE
#ifdef PER_BLOCK
	vector<real> v_Iab(Graph::block_prob.size(), 0.0), v_ilb(Graph::block_prob.size(), 0.0), v_Sab(Graph::block_prob.size(), 0.0);
#else
	vector<real> v_Iv(Graph::n, 0.0), v_Sv(Graph::n, 0.0);
#endif
#endif //CLIQUE
	while (scenario < numScenarios){
		_numAgents = _startingNumAg + (scenario * step);
		
		//Agent-wise reset:
#ifdef CLIQUE
		for (node v = 0; v < Graph::n; ++v) Graph::g[v] = v;		// ----> This will be important during the "protection effect" process, as the nodes will be re-arranged. 
#endif //CLIQUE
		Stats::resetAvDur();
		for (uint r = 0; r < ROUNDS; ++r) totalSimTime[r] = 0;
		iaTotal = (I_0 > NUM_AGENTS) ? NUM_AGENTS : I_0;
		
		//Initializing the applicable stats:
#ifdef INFECTED_FRACTION
		Stats::initStream(stat::streamType::infFrac);
#endif
		Stats::initStream(stat::streamType::avDuration);
		Reporter::startChronometer("\n\nRunning scenario " + std::to_string(scenario + 1) + "/" + std::to_string(numScenarios) + "...");
		Reporter::simulationInfo(iaTotal);
		for (uint round = 0; round < ROUNDS; ++round) {
#ifdef MEASURE_ROUND_EXE_TIME
			Reporter::startChronometer("\n  Round " + std::to_string(round + 1) + "...");
#else
			//Reporter::progress(round);
#endif
			//Round-wise reset:
			resetVariables();

			//Initially, we schedule a single 'walk' event for each agent. A new 'walk' job will then be created and scheduled for an agent the moment its current 'walk' job is processed.
			for (agent i = 0; i < _numAgents; ++i)
				schedule.emplace(i, EXPLambda(), action::walk);

			//Now we must define which agents are infected from the start, and then schedule their 'recover' event.
			//It is not a problem to always infect the first 'itotal' agents, since the starting node for each of them will be randomly set.
			for (agent i = 0; i < iaTotal; ++i) {
				isInfectedAg[i] = true;
				schedule.emplace(i, EXPGamma_a(), action::recoverAg);
			}

			// * DISTRIBUTING THE AGENTS ACROSS THE NETWORK *
			//Random distribution of the INFECTED agents:
#ifdef CLIQUE
			for (agent i = 0; i < iaTotal; ++i)
				enterNodeAsInf(i, randomInt(Graph::n), TIME_ZERO);
#else //CLIQUE
			for (agent i = 0; i < iaTotal; ++i) {
				const node v = randomLCCNode();
				enterNodeAsInf(i, v, TIME_ZERO);
#ifndef PER_BLOCK
				++v_Iv[v];
#endif
			}
#endif //CLIQUE

			//Random distribution of the SUSCEPTIBLE agents:
#ifdef CLIQUE
			for (agent i = iaTotal; i < _numAgents; ++i)
				enterNodeAsSus(i, randomInt(Graph::n), TIME_ZERO);
#else //CLIQUE
			for (agent i = iaTotal; i < _numAgents; ++i) {
				const node v = randomLCCNode();
				enterNodeAsSus(i, v, TIME_ZERO);
#ifndef PER_BLOCK
				++v_Sv[v];
#endif
			}
#endif //CLIQUE

#ifdef CLIQUE
#ifdef PER_BLOCK
			Ia = iaTotal;
			Sa = saTotal;
#endif //PER_BLOCK
#else
#ifdef PER_BLOCK
			//Expected number of S-/I-agents within each node from each block:
			{
				const double ia = (double)iaTotal / NUM_AGENTS, sa = (double)saTotal / NUM_AGENTS;
				for (uint b = (uint)graph::Graph::block_prob.size() - 1; b > 0; --b) {
					v_Iab[b] = ia * graph::Graph::kb[b];
					v_Sab[b] = sa * graph::Graph::kb[b];
				}
#ifdef DEBUG
				double sum = 0.0;
				for (uint b = (uint)graph::Graph::block_prob.size() - 1; b > 0; --b) {
					const double nb = (double)graph::Graph::n * graph::Graph::block_prob[b];
					sum += nb * (ivb[b] + svb[b]);
				}
				assertm(abs(sum - NUM_AGENTS) < epsilon, "The computed *expected number of S-/I-agents per node* is incorrect.");
#endif //DEBUG
			}
#endif //PER_BLOCK
#endif //CLIQUE
			
			// * MAIN LOOP *
			bool earlyStop = false;
			real& roundDuration = totalSimTime[round]; 
			double timeLimit = T;
			job j;
			nextJob(j, now);
			const real logInterval = 0.01 * T; // ----> Log progress to the console at one-percent increments.
			real prevLog = 0;
			std::cout << std::fixed;
			std::cout << std::setprecision(1);
			std::cout << "\r0% complete...";
#ifdef ONLY_NUMERIC
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
				case action::siteInfectAg:
					agFate_fromSite(j.target, now, j.infective, j.validity_S, j.validity_I);
					break;
				case action::agInfectSite:
					siteFate(j.target, now, j.infective, j.validity_S, j.validity_I);
					break;
				case action::recoverSite:
					recoverSite(j.target, now);
					break;
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
					std::cout << '\r' << std::round((now/T) * 100) << "% complete...";
				}
			}
			std::cout << '\r' << "100% complete!  ";
			Stats::partialsAvDur(roundDuration);
#ifdef MEASURE_ROUND_EXE_TIME
			Reporter::stopChronometer((earlyStop) ? "done-IV" : "done-TL");	// ----> TL == Time Limit; IV == Infection Vanished.
			Reporter::durationInfo(roundDuration);
#endif
#ifdef INFECTED_FRACTION
			Reporter::startChronometer(" Saving To File...");	
			Stats::iFracToFile(OVERLOOK);
			Stats::endStream(stat::streamType::infFrac);
			Reporter::stopChronometer(" done");
#endif
		} // ** for (uint round = 0; round < ROUNDS; ++round)

		Stats::writeToFile(stat::streamType::avDuration, Ws, Wi, _numAgents);
		Stats::endStream  (stat::streamType::avDuration);
#ifndef MEASURE_ROUND_EXE_TIME
		Reporter::tell(" All rounds completed.\n");
#endif
		Reporter::stopChronometer("Scenario " + std::to_string(scenario + 1) + "/" + std::to_string(numScenarios) + " completed");
		Reporter::avSimTimeInfo(Stats::avDuration());
		++scenario;
	} // ** while (scenario < numScenarios)

#ifdef SOLVE_NUMERICALLY
	//Runge-Kutta:
	constexpr uint outputGranularity = 50;
	constexpr uint largerDetailUntil = 100;
	constexpr real stepSize = 0.00001;
	constexpr real epsilon = 1.0 / N ;
	constexpr real timeIncrement = stepSize * outputGranularity;
	vector<real> saveToFile_diadt;
	vector<real> saveToFile_dildt;
	uint outputSize = 0;
	//Stats::setBasename();
	std::stringstream name;
	name << SHORT_LABEL 
		<< "_N"		<< N 
		<< "_AG"	<< NUM_AGENTS 
		<< "_T"	<< TAU_aa 
		<< "_G"	<< GAMMA_a 
		<< "_L"		<< LAMBDA 
		<< "_STime" << T 
		<< "_R"		<< ROUNDS;
		//<< "_Tal"	<< TAU_al 
		//<< "_Tla"	<< TAU_la 
		//<< "_Gl"	<< GAMMA_l 
	baseName = name.str();

#ifdef CLIQUE
#ifdef PER_BLOCK
	rungeKutta4thOrder(0, Ia, Sa, T, stepSize, epsilon, saveToFile_diadt, saveToFile_dildt, outputSize, outputGranularity, largerDetailUntil);
#else
	rungeKutta4thOrder(0, v_Iv, v_Sv, v_ilb, T, stepSize, epsilon, saveToFile_diadt, saveToFile_dildt, outputSize, outputGranularity, largerDetailUntil);
#endif
#else //CLIQUE
#ifdef PER_BLOCK
	rungeKutta4thOrder(0, v_Iab, v_Sab, v_ilb, T, stepSize, epsilon, saveToFile_diadt, saveToFile_dildt, outputSize, outputGranularity, largerDetailUntil);
#else
	rungeKutta4thOrder(0, v_Iv, v_Sv, v_ilb, T, stepSize, epsilon, saveToFile_diadt, saveToFile_dildt, outputSize, outputGranularity, largerDetailUntil);
#endif
#endif //CLIQUE
	//Saving to file:
	std::ofstream RKdata;
	std::stringstream _ss;
	_ss.precision(4);
	_ss << "./stats/Runge-Kutta_" << baseName << ".csv";
	std::string _fileName = _ss.str();
	RKdata.open(_fileName);
	RKdata << "Time\tia\til\n"; // ----> Header
	real _time = 0;
	for (size_t s = 0; s < largerDetailUntil; ++s) {
		RKdata << _time << '\t' << saveToFile_diadt[s] << '\t' << saveToFile_diadt[s] << '\n';
		_time += stepSize;
	}
	for (size_t s = largerDetailUntil; s < outputSize; ++s){
		RKdata << _time << '\t' << saveToFile_diadt[s] << '\t' << saveToFile_diadt[s] << '\n';
		_time += timeIncrement;
	}
	RKdata.close();
#endif //SOLVE_NUMERICALLY
}
