#include "graph.h"
#include "reporter.h"
#include "stats.h"
#include <numeric>		//std::accumulate()

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

	job(const uint&& _target, const real&& _time, const action&& _action) : target(_target), time(_time), a(_action) {}
	job(const uint&& _target, const real&& _time, const action&& _action, const uint&& _infective, const uint&& _snapshot_S, const uint&& _snapshot_I) : target(_target), time(_time), a(_action), infective(_infective), validity_S(_snapshot_S), validity_I(_snapshot_I) {}
	job() {}
	job(const uint& _target, const real&& _time, const action&& _action) : target(_target), time(_time), a(_action) {}
	job(const uint& _target, const real& _time, const action&& _action, const uint& _infective, const uint& _snapshot_S, const uint& _snapshot_I) : target(_target), time(_time), a(_action), infective(_infective), validity_S(_snapshot_S), validity_I(_snapshot_I) {}

	job(const job& other) { *this = other; }						// ----> Copy Contructor. Jobs will be handled by an STL vector (via priority-queue). It thus requires a copy contructor to be explicitly defined. 

};
struct earlier {
	bool operator()(const job& j1, const job& j2) {
		return j1.time > j2.time;
	}
};

std::priority_queue<job, std::vector<job>, earlier> schedule;


//Simulation variables:
uint saTotal;														// ----> Up-to-date number of SUSCEPTIBLE AGENTS during the simulation.
uint iaTotal;														// ----> Up-to-date number of INFECTED AGENTS during the simulation.
uint ilTotal = 0;													// ----> Up-to-date number of INFECTED SITES during the simulation.
real now;
std::vector<real> totalSimTime(ROUNDS, 0);							// ----> Total simulation time at each round, to average upon.
//real sim::beta_a;													// ----> Force of infection from an I-agent to an S-agent
//real sim::beta_al;													// ----> Force of infection from an I-agent to a site
//real sim::beta_la;													// ----> Force of infection from a site to an I-agent
//std::string sim::baseName;
//real beta2ndMmt_logistic;
//real beta2ndMmt_naive;
//long real C_2ndMmt;
//long real C_2ndMmt_logistic;
//long real C_2ndMmt_naive;

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
void agFate_fromAg	(const agent& ag, const real& now, const uint& infective, const uint& validity_S, const uint& validity_I);		// ----> Determines whether or not an exposed, susceptible agent 'ag' will become infected.
void agFate_fromSite(const agent& ag, const real& now, const uint& infective, const uint& validity_S, const uint& validity_I);		// ----> Determines whether or not an exposed, susceptible agent 'ag' will become infected.
void siteFate		(const uint& v, const real& now, const uint& infective, const uint& validity_L, const uint& validity_I);		// ----> Determines whether or not an exposed, susceptible agent 'ag' will become infected.
void enterNodeAsSus (const agent& ag, const node& v, const real& now);
//void checkinAsSus   (const agent& ag, const node& v);
//void checkoutAsSus  (const agent& ag, const node& v);
//void checkinAsInf   (const agent& ag, const node& v);
//void checkoutAsInf  (const agent& ag, const node& v);
void check_in		(const agent& ag, const node& v, vector<vector<agent>>& _where);
void check_out		(const agent& ag, const node& v, vector<vector<agent>>& _where);
void enterNodeAsInf (const agent& ag, const node& v, const real& now);
void leaveNodeAsSus (const agent& ag, const node& v, const real& now);
void leaveNodeAsInf (const agent& ag, const node& v, const real& now);
//Defines the next node an agent is going to visit.
const node& nextNodeForSus(const node& _currNode);
const node& nextNodeForInf(const node& _currNode);
void nextJob(job& _j, real& _time);
void setBeta2ndMmt();
void getSigma_a(const double& Ia, double& sigma_as, double& sigma_ai);
//double getSigma_a();
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
static const real sim::EXPTau_aa()	{ return NEG_RECIPR_TAU_aa	* log(U()); }
static const real sim::EXPTau_al()	{ return NEG_RECIPR_TAU_al	* log(U()); }
static const real sim::EXPTau_la()	{ return NEG_RECIPR_TAU_la	* log(U()); }
static const real sim::EXPGamma_a()	{ return NEG_RECIPR_GAMMA_a	* log(U()); }
static const real sim::EXPGamma_l()	{ return NEG_RECIPR_GAMMA_l	* log(U()); }
static const real sim::EXPLambda()	{ return NEG_RECIPR_LAMBDA	* log(U()); }
#ifdef i_t_FROM_MODEL

//real sim::i_t_2ndMmt(const real& t) {
//	double val = (beta2ndMmt - GAMMA_a) * t;
//	if(val > 500) val = 500;
//
//	//long real Ce = C_2ndMmt * exp((beta2ndMmt - GAMMA) * t);
//	long real Ce = C_2ndMmt * exp(val);
//
//	return std::max((1.0 - (GAMMA_a/beta2ndMmt)) * (Ce / (1.0 + Ce)), (long real)0.0);
//}
//
//real sim::i_t_2ndMmt_naive(const real& t) {
//	double val = (t > 1000) ? 1000 : t;
//
//	long real Ce = C_2ndMmt_naive * exp((beta2ndMmt_naive - GAMMA) * val);
//	return std::max((1.0 - (GAMMA / beta2ndMmt_naive)) * (Ce / (1.0 + Ce)), (long real)0.0);
//}
//real sim::i_t_2ndMmt_sys(const real& t) {
//	for (size_t i = 0; i < sim::block_prob.size(); i++) {
//
//		long real Ce = C_2ndMmt * exp((beta2ndMmt - GAMMA) * t);
//		return std::max((1.0 - (GAMMA / beta2ndMmt)) * (Ce / (1.0 + Ce)), (long real)0.0);
//	}
//}
//real sim::i_t_pfx(const real& t) {
//	long real Ce = sim::C_pfx * exp(B_MINUS_G_pfx * t);
//	return std::max((_1_MINUS_G_OVER_B_pfx) * (Ce / (1.0 + Ce)), (long real)0.0);
//}
//void sim::setBeta2ndMmt() {
//	
//	//1st approach: 
//	vector<double> beta_b(graph::Graph::block_prob.size(), 0);
//	for (size_t _b = 1; _b < graph::Graph::block_prob.size(); ++_b){
//		double _kb_ = (NUM_AGENTS * _b * N * graph::Graph::block_prob[_b]) / (2 * graph::Graph::m);
//		//if (_kb_ < N * graph::Graph::block_prob[_b]) {
//		if (_kb_ < N * graph::Graph::block_prob[_b]) {
//			beta_b[_b] = (LAMBDA * _b * (TAU / (2 * LAMBDA + TAU)) * _kb_) / graph::Graph::m;
//		}
//		else {
//			beta_b[_b] = (LAMBDA * _b * (TAU / (LAMBDA + TAU)) * N * graph::Graph::block_prob[_b]) / graph::Graph::m;
//		}
//	}
//	beta2ndMmt = 0;
//	for (size_t _b = 1; _b < graph::Graph::block_prob.size(); ++_b) {
//		beta2ndMmt += beta_b[_b];
//	}
//
//	//2nd approach:
//	beta_b.resize(0);
//	beta_b.resize(graph::Graph::block_prob.size(), 0);
//	double sigma = 0;
//	double sigma_2 = 0;
//	uint validBlock = 0;
//	double expBlock = graph::Graph::_2ndMmt / graph::Graph::averageDegree;
//	double max_kb = (((double)NUM_AGENTS) * (double)(graph::Graph::block_prob.size() - 1) * (double)N * graph::Graph::block_prob[graph::Graph::block_prob.size() - 1]) / (2.0 * graph::Graph::m);
//	double variance = abs(graph::Graph::_2ndMmt - pow(graph::Graph::averageDegree, 2));
//	double _diff = abs(graph::Graph::averageDegree - graph::Graph::_2ndMmt);
//	double std_dev = abs(sqrt(variance));
//	const double& avDg = graph::Graph::averageDegree;
//	for (size_t _b = 1; _b < graph::Graph::block_prob.size(); ++_b) {
//		double _kb_ = (((double)NUM_AGENTS) * (double)_b * (double)N * graph::Graph::block_prob[_b]) / (2.0 * graph::Graph::m);
//		//double kbnb = _kb_ / N * graph::Graph::block_prob[_b];
//		double norm_kb = (_kb_/N)/ ((_kb_/N)+ std::exp(-(_kb_/N)));
//		//beta_b[_b] = (LAMBDA * _b * (TAU / (LAMBDA + TAU + (LAMBDA * (1- norm_kbnb)) )) * (_kb_ * (1 - norm_kbnb))) / graph::Graph::m;
//		beta_b[_b] = (LAMBDA * _b * (TAU / (2*LAMBDA + TAU) * norm_kb * N )) / graph::Graph::m;
//		
//		if (graph::Graph::block_prob[_b] > 0) {
//			++validBlock;
//			//double dimmer = 1.0 - ((_kb_/(4.0 + _kb_)) / ((_kb_/ (4.0 + _kb_)) + exp(-(_kb_/ (4.0 + _kb_)))));
//			//double dimmer = 1.0 - ((_kb_/(expBlock)) / ((_kb_/ (expBlock)) + exp(-(_kb_/ (expBlock)))));
//			//double val = (_kb_) / (abs(std_dev - (log2(avDg) * _kb_)));
//			//double val = (_kb_) / (1 + abs(variance - avDg));
//			//double val = (_kb_) / (expBlock + (2 * _kb_));
//			//double val = graph::Graph::averageDegree / graph::Graph::_2ndMmt;
//			//double val = _kb_/(abs(expBlock - avDg) + _kb_);
//			//double val = _kb_/(avDg + _kb_);
//			//double val = _kb_/(log2(avDg)*max_kb);
//			//double val = _kb_ / ((log2(expBlock) * max_kb) + _kb_);
//			double val = _kb_/((log2(expBlock*max_kb))+_kb_);
//			//double dimmer = 1.0 - (val / (val + exp(-val)));
//
//			//Bom resultado:
//			double nTau = (double)(TAU) / std::max(1.0, _kb_ * std::min(1.0, graph::Graph::_2ndMmt/avDg));
//			
//			//double dimmer = 1.0 - ((_kb_/(2*expBlock)) / ((_kb_/ (2*expBlock)) + exp(-(_kb_/ (2*expBlock)))));
//			//double dimmer = 1.0 - ((1.0/(2*_kb_)) / (((1.0 / (2 * _kb_))) + exp(-(1.0 / (2 * _kb_)))));
//			
//			//double dimmer = (_kb_ < 2.0) ? 1.0 : 1.0 - (_kb_) / (_kb_ + exp(-(1.0 / (_kb_))));
//			//sigma += ((double)TAU * dimmer) / (2 * LAMBDA + (TAU * dimmer));
//			
//			//sigma += (double)TAU / (LAMBDA + TAU + LAMBDA * (1.0 - (std::min(1.0,(_kb_/2.0)))));
//
//
//			//MUITO BOM NO BA (O NORMAL É MELHOR NO G(N,P)):
//			sigma +=  (_kb_ < 1.5) ? TAU / (2*LAMBDA + TAU) : TAU / (LAMBDA + TAU);
//			
//			//sigma +=  (_kb_ < 1.5) ? TAU / (2*LAMBDA + TAU) : (TAU*dimmer) / (LAMBDA + (TAU*dimmer));
//			sigma_2 += nTau / (2 * LAMBDA + nTau);
//		}
//	}
//	sigma /= validBlock;
//	sigma_2 /= validBlock;
//	beta2ndMmt_logistic = 0;
//	for (size_t _b = 1; _b < graph::Graph::block_prob.size(); ++_b) {
//		beta2ndMmt_logistic += beta_b[_b];
//	}
//
//	//double kOverN = (double)NUM_AGENTS / N;
//	//double sigma = (double)TAU / (LAMBDA + TAU + LAMBDA * (1 - _k_n));
//	//double _k_n = (kOverN) / (kOverN + std::exp(-(kOverN * (std::_Pi/2))));
//	//beta2ndMmt_logistic = (double)((LAMBDA * _b * sigma * _k_n)) / graph::Graph::m;
//
//
//	//double probC = 0;	// ----> Probability "correction".
//	//probC =  TAU/(LAMBDA + TAU);
//
//	//3rd approach:
//	//double _aux = (excess * (1.0 / std::exp(log(N) - log(NUM_AGENTS))));
//	//beta2ndMmt = (meetingRate * SIGMA_aa * (std::floor((double)NUM_AGENTS - 2*generalExpOcc - _aux)) * graph::Graph::_2ndMmt) / pow(graph::Graph::averageDegree, 2);
//	
//	//beta2ndMmt_naive = (meetingRate * SIGMA_aa * NUM_AGENTS * graph::Graph::_2ndMmt) / pow(graph::Graph::averageDegree, 2);
//	
//	//MUITO BOM NO BA (O NORMAL É MELHOR NO G(N,P)):
//	//beta2ndMmt = (meetingRate * sigma * NUM_AGENTS * graph::Graph::_2ndMmt) / pow(graph::Graph::averageDegree, 2);
//	beta2ndMmt = (meetingRate * SIGMA_aa * NUM_AGENTS * graph::Graph::_2ndMmt) / pow(graph::Graph::averageDegree, 2);
//	
//	
//	//Bom resultado na BA:
//	//beta2ndMmt_naive = (meetingRate * sigma_2 * NUM_AGENTS * graph::Graph::_2ndMmt) / pow(graph::Graph::averageDegree, 2);
//	//beta2ndMmt_naive = (LAMBDA * SIGMA_aa * N * graph::Graph::_2ndMmt) / (graph::Graph::m * avDg);
//
//	C_2ndMmt = i_0 / (1.0 - i_0 - (GAMMA / beta2ndMmt));
//	C_2ndMmt_logistic = i_0 / (1.0 - i_0 - (GAMMA / beta2ndMmt_logistic));
//	C_2ndMmt_naive = i_0 / (1.0 - i_0 - (GAMMA / beta2ndMmt_naive));
//}
#endif //i_t_FROM_MODEL

void sim::setBeta2ndMmt() {
	//beta_a = (double)((2.0 * LAMBDA * SIGMA_aa * NUM_AGENTS * graph::Graph::_2ndMmt)) / (N * pow(graph::Graph::averageDegree, 2));
	//MELHOR COM SIGMA FIXO:
	//beta_a = (double)(2 * LAMBDA * SIGMA_aa * graph::Graph::psi) / graph::Graph::averageDegree;
	
	beta_a = (double)(LAMBDA * SIGMA_aa * graph::Graph::psi) / graph::Graph::averageDegree;
	
	beta_al = (LAMBDA * NUM_AGENTS *  SIGMA_al)/N;
	beta_la = (LAMBDA * graph::Graph::_2ndMmt * N * SIGMA_la) / (N * pow(graph::Graph::averageDegree, 2));
}

void sim::getSigma_a(const double& Ia, double& sigma_as, double& sigma_ai) {
	sigma_as = 0;
	sigma_ai = 0;
	for (uint b = 1; b < graph::Graph::block_prob.size(); ++b) {
		if (graph::Graph::kb[b] == 0)
			continue;
		const double expNumAg = graph::Graph::kb[b];
		const double expNumInfAg = expNumAg * Ia;			// ----> Talvez seja errado fazer dessa forma...
		const double minInfAg = std::min(1.0, expNumInfAg);	// ----> REVER! Talvesz o mínimo unitário faça mais sentido (fora que evita erros de número muito pequeno em cenários extremamente esparsos, onde esse número seria mt próx de zero).
		const double expNumSusAg = expNumAg * (1.0 - Ia);
		constexpr double euler = 0.5772156649;
		const double digamma_i = log(expNumInfAg) - 1.0 / (2 * expNumInfAg);
		const double digamma_s = log(expNumSusAg) - 1.0 / (2 * expNumSusAg);
		const double H_i = std::max(1.0, euler + digamma_i);
		const double H_s = std::max(1.0, euler + digamma_s);
		const double prob_NoAcq = (H_i * LAMBDA + H_i * GAMMA_a + LAMBDA) / (H_i * LAMBDA + LAMBDA + H_i * GAMMA_a + expNumInfAg * TAU_aa);
		const double prob_acq = (std::max(expNumInfAg, 1.0) * TAU_aa) / ((1.0 + H_i) * LAMBDA + std::max(expNumInfAg, 1.0) * TAU_aa);
		const double prob_inf = (TAU_aa) / (2 * LAMBDA + std::max(expNumInfAg, 1.0) * TAU_aa);
		const double prob_NoTransmission = (H_s * LAMBDA + LAMBDA + GAMMA_a) / (H_s * LAMBDA + LAMBDA + GAMMA_a + TAU_aa);
		sigma_as += prob_acq * graph::Graph::block_prob[b];
		sigma_ai += prob_inf * graph::Graph::block_prob[b];
	}

}

#ifdef SOLVE_NUMERICALLY
real sim::diabdt(const real& Ia, const real& Iab, const real& Sab, const uint& b) {
	

	//const double& rho = graph::Graph::rho_b[b];
	//const double notInNode = (double)b / (2 * graph::Graph::m);
	//double noOne;
	//double rho_s;
	//double rho_i;
	//if (Iab > Sab) {
	//	noOne = pow(1.0 - notInNode, round(Sab));
	//	rho_s = 1.0 - noOne;
	//	rho_i = 1.0 - pow(noOne, round(Iab) - round(Sab));
	//}
	//else if (Sab > Iab) {
	//	noOne = pow(1.0 - notInNode, round(Iab));
	//	rho_i = 1.0 - noOne;
	//	rho_s = 1.0 - pow(noOne, round(Sab) - round(Iab));
	//}
	//else {	// ----> Iab == Sab
	//	noOne = pow(1.0 - notInNode, round(Iab));
	//	rho_s = rho_i = 1.0 - noOne;
	//}

	//rho_s = 1.0 - pow(1.0 - ((double)b / (2 * graph::Graph::m)), Sab);
	//rho_i = 1.0 - pow(1.0 - ((double)b / (2 * graph::Graph::m)), Iab);

	const double& pb = graph::Graph::block_prob[b];
	const double nb = graph::Graph::n * pb;
	const double& qb = graph::Graph::q_b[b];
	//const double out = ((double)(b - 1) / b) * (1.0 - pb);
	//const double out = 1.0 - (((double)b * nb) / (2.0 * graph::Graph::m));	// ----> bem legal perto de 50%. Subestima saturação
	//const double out = 1.0 - (((double)(b-1) * nb) / (2.0 * (graph::Graph::m - graph::Graph::n)));	// --> bom perto da saturação. Superrestima regime intermediário.
	//const double out = 1.0 - pb; // ----> superestima um pouco...
	const double out = ((double)graph::Graph::m - (graph::Graph::n - nb) - (b * nb)) / ((double)graph::Graph::m - (graph::Graph::n - nb));

	//DON:
	if (Sab + Iab == 0)
		return LAMBDA * Ia * qb;
	//return (Ia - Iab) * LAMBDA * qb - Iab * LAMBDA * (1.0 - qb) + ((Sab * Iab * TAU_aa) / nb) - (GAMMA_a * Iab);
	//TESTE!!!
	return (Ia - Iab) * LAMBDA * qb - Iab * LAMBDA * out + ((Sab * Iab * TAU_aa) / nb) - (GAMMA_a * Iab);
	//return Ia * LAMBDA * qb - Iab * LAMBDA * (1.0 - qb) + ((Sab * Iab * TAU_aa) / nb) - (GAMMA_a * Iab);
	//return (LAMBDA * (Ia * qb - Iab)) + ((Sab * Iab * TAU_aa) / nb) - (GAMMA_a * Iab);
	//return (LAMBDA * (Ia * qb - Iab)) - (GAMMA_a * Iab);
	 
	////const double _kb = (Sab + Iab);
	////const double saturation = 1.0 - (0.25 - (Sab / _kb) * (Iab / _kb));
	////return LAMBDA * (Ia * qb - Iab) + (Sab * Iab * (saturation * TAU_aa)) / nb - GAMMA_a * Iab;
	
	//const double Ib = NUM_AGENTS * Iab;
	//const double Sb = NUM_AGENTS * Sab;
	////const double rateAllIagsRecover = Ib * GAMMA_a;
	////const double sigma_ai = TAU_aa / (2 * LAMBDA + Ib * TAU_aa + rateAllIagsRecover);
	
	//double _1to1 = (TAU_aa / (2 * LAMBDA + Iab * TAU_aa));
	//double noTransmission = pow(1.0 - _1to1, Sab);
	//const double sigma_ai = 1.0 - noTransmission;
	//if (Iab == 0) {
	//	if (Sab == 0)
	//		return;
	//	return (Ia * LAMBDA * qb * (Sab / nb) * sigma_ai);
	//}
	//
	//const double probNoInf = pow(1.0 - _1to1, Iab);
	//const double sigma_as = 1.0 - probNoInf;
	//
	////return LAMBDA * (Ia * qb - Iab) + (Ia * LAMBDA * qb * (Sab/nb) * TAU_aa) + ((NUM_AGENTS - Ia) * LAMBDA * qb * (Iab / nb) * TAU_aa) - GAMMA_a * Iab;
	//
	////return (Ia * LAMBDA * qb * (Sab/nb) * TAU_aa) + ((NUM_AGENTS - Ia) * LAMBDA * qb * (Iab / nb) * TAU_aa) - GAMMA_a * Iab;
	//return (Ia * LAMBDA * qb * (Sab / nb) * sigma_ai) + ((NUM_AGENTS - Ia) * LAMBDA * qb * (Iab / nb) * sigma_as) - GAMMA_a * Iab;

}

real sim::dsabdt(const real& Ia, const real& Iab, const real& Sab, const uint& b) {
	
	//const double& rho = graph::Graph::rho_b[b];

	//const double rho_s = 1.0 - pow(1.0 - ((double)b / (2 * graph::Graph::m)), Sab);
	//const double rho_i = 1.0 - pow(1.0 - ((double)b / (2 * graph::Graph::m)), Iab);

	const double& pb = graph::Graph::block_prob[b];
	const double nb = graph::Graph::n * pb;
	const double& qb = graph::Graph::q_b[b];
	const double out = ((double)graph::Graph::m - (graph::Graph::n - nb) - (b * nb)) / ((double)graph::Graph::m - (graph::Graph::n - nb));

	//DON:
	if (Sab + Iab == 0)
		return LAMBDA * ((double)NUM_AGENTS - Ia) * qb;
	//return (((double)NUM_AGENTS - Ia) - Sab) * LAMBDA * qb - Sab * LAMBDA * (1.0 - qb) - ((Sab * Iab * TAU_aa) / nb) + (GAMMA_a * Iab);
	//TESTE!!!
	return (((double)NUM_AGENTS - Ia) - Sab) * LAMBDA * qb - Sab * LAMBDA * out - ((Sab * Iab * TAU_aa) / nb) + (GAMMA_a * Iab);
	//return ((double)NUM_AGENTS - Ia)  * LAMBDA * qb - Sab * LAMBDA * (1.0 - qb) - ((Sab * Iab * TAU_aa) / nb) + (GAMMA_a * Iab);
	//return (LAMBDA * (((double)NUM_AGENTS - Ia) * qb - Sab)) - ((Sab * Iab * TAU_aa) / nb) + (GAMMA_a * Iab);
	//return (LAMBDA * (((double)NUM_AGENTS - Ia) * qb - Sab)) + (GAMMA_a * Iab);

	////const double _kb = (Sab + Iab);
	////const double saturation = 1.0 - (0.25 - (Sab / _kb) * (Iab / _kb));
	////return LAMBDA * ((NUM_AGENTS - Ia) * qb - Sab) - (Sab * Iab * (saturation * TAU_aa)) / nb + GAMMA_a * Iab;

	//const double Ib = NUM_AGENTS * Iab;
	//const double Sb = NUM_AGENTS * Sab;
	////const double rateAllIagsRecover = Ib * GAMMA_a;
	////const double sigma_ai = TAU_aa / (2 * LAMBDA + Ib * TAU_aa + rateAllIagsRecover);
	
	//double _1to1 = (TAU_aa / (2 * LAMBDA + Iab * TAU_aa));
	//double noTransmission = pow(1.0 - _1to1, Sab);
	//const double sigma_ai = 1.0 - noTransmission;
	//if (Iab == 0) {
	//	if (Sab == 0)
	//		return;
	//	return -(Ia * LAMBDA * qb * (Sab / nb) * sigma_ai);
	//}
	//
	//const double probNoInf = pow(1.0 - _1to1, Iab);
	//const double sigma_as = 1.0 - probNoInf;
	//
	////return LAMBDA * ((NUM_AGENTS - Ia) * qb - Sab) - (Ia * LAMBDA * qb * (Sab / nb) * TAU_aa) - ((NUM_AGENTS - Ia) * LAMBDA * qb * (Iab / nb) * TAU_aa) + GAMMA_a * Iab;
	//
	//return -(Ia * LAMBDA * qb * (Sab / nb) * sigma_ai) - ((NUM_AGENTS - Ia) * LAMBDA * qb * (Iab / nb) * sigma_as) + GAMMA_a * Iab;

}

void sim::takeStep(const real& h, real& Ia, std::vector<real>& v_Iab, std::vector<real>& v_Sab) {
	constexpr real one_sixth = 1.0 / 6.0;
	const uint blocks = static_cast<uint>(graph::Graph::block_prob.size());
	vector<real> k1(2 * blocks, 0), k2(2 * blocks, 0), k3(2 * blocks, 0), k4(2 * blocks, 0);
	
	Ia = 0;
	uint i;
	for (uint b = (uint)graph::Graph::block_prob.size() - 1; b > 0; --b) {
		Ia += v_Iab[b];
	}
	for (uint b = (uint)graph::Graph::block_prob.size() - 1; b > 0; --b) {
		if (graph::Graph::block_prob[b] == 0)
			continue;

		i = 2 * b;
		k1[i]		= h * diabdt(Ia, v_Iab[b], v_Sab[b], b);
		k1[i + 1]	= h * dsabdt(Ia, v_Iab[b], v_Sab[b], b);
	}

	Ia = 0;
	for (uint b = (uint)graph::Graph::block_prob.size() - 1; b > 0; --b) 
		Ia += v_Iab[b] + 0.5 * k1[2 * b];
	for (uint b = (uint)graph::Graph::block_prob.size() - 1; b > 0; --b) {
		if (graph::Graph::block_prob[b] == 0)
			continue;

		i = 2 * b;
		k2[i]		= h * diabdt(Ia, v_Iab[b] + 0.5 * k1[i], v_Sab[b] + 0.5 * k1[i + 1], b);
		k2[i + 1]	= h * dsabdt(Ia, v_Iab[b] + 0.5 * k1[i], v_Sab[b] + 0.5 * k1[i + 1], b);
	}

	Ia = 0;
	for (uint b = (uint)graph::Graph::block_prob.size() - 1; b > 0; --b)
		Ia += v_Iab[b] + 0.5 * k2[2 * b];
	for (uint b = (uint)graph::Graph::block_prob.size() - 1; b > 0; --b) {
		if (graph::Graph::block_prob[b] == 0)
			continue;

		i = 2 * b;
		k3[i]		= h * diabdt(Ia, v_Iab[b] + 0.5 * k2[i], v_Sab[b] + 0.5 * k2[i + 1], b);
		k3[i + 1]	= h * dsabdt(Ia, v_Iab[b] + 0.5 * k2[i], v_Sab[b] + 0.5 * k2[i + 1], b);
	}

	Ia = 0;
	for (uint b = (uint)graph::Graph::block_prob.size() - 1; b > 0; --b)
		Ia += v_Iab[b] + k3[2 * b];
	for (uint b = (uint)graph::Graph::block_prob.size() - 1; b > 0; --b) {
		if (graph::Graph::block_prob[b] == 0)
			continue;

		i = 2 * b;
		k4[i]		= h * diabdt(Ia, v_Iab[b] + k3[i], v_Sab[b] + k3[i + 1], b);
		k4[i + 1]	= h * dsabdt(Ia, v_Iab[b] + k3[i], v_Sab[b] + k3[i + 1], b);
	}

	for (uint b = (uint)graph::Graph::block_prob.size() - 1; b > 0; --b) {
		if (graph::Graph::block_prob[b] == 0)
			continue;

		i = 2 * b;
		v_Iab[b] += one_sixth * (k1[i]			+ 2 * k2[i]			+ 2 * k3[i]			+ k4[i]);
		v_Sab[b] += one_sixth * (k1[i + 1]		+ 2 * k2[i + 1]		+ 2 * k3[i + 1]		+ k4[i + 1]);
	}
}

real sim::dilbdt(const real& Ia, const real& il, const real& Iab, const real& ilb, const uint& b) {
	return 0;
}

//real sim::diadt(const real& Ia, const real& il) {
real sim::diadt(const real& Ia, const double& sumSbIb) {

	return graph::Graph::psi * TAU_aa * sumSbIb - GAMMA_a * Ia;
	//return - GAMMA_a * Ia;
}
real sim::dildt(const real& Ia, const real& il) {
	return  beta_al * (1.0 - il) * Ia - GAMMA_l * il;
}

void sim::rungeKutta4thOrder(const real& t0, std::vector<real>& v_Iab, std::vector<real>& v_Sab, std::vector<real>& v_ilb, const real& t, const real& h, const real& epsilon, vector<real>& saveToFile_diadt, vector<real>& saveToFile_dildt, uint& outputSize, const uint& outputGranularity, const real& largerDetailUntil) {
	uint totalSteps = (uint)((t - t0) / h) + 1;
	saveToFile_diadt.resize((uint64_t)largerDetailUntil + (totalSteps - ((uint)largerDetailUntil) / outputGranularity) + 1);
	saveToFile_dildt.resize((uint64_t)largerDetailUntil + (totalSteps - ((uint)largerDetailUntil) / outputGranularity) + 1);

	double Ia = 0;
	for (size_t i = 1; i < v_Iab.size(); ++i)
		Ia += v_Iab[i];
	//double Ia = std::accumulate(v_Iab.begin(), v_Iab.end(), 0.0);
	//real il = std::accumulate(v_ilb.begin(), v_ilb.end(), 0);
	real il = 0.0;

	saveToFile_diadt[0] = Ia / NUM_AGENTS;
	saveToFile_dildt[0]	= il;
	bool end = false;
	++outputSize;

	//For the first 'largerDetailUntil' iterations every step is stored in a vector ('saveToFile'), for later being written to file.
	for (uint s = 1; s < largerDetailUntil; ++s) {
		takeStep(h, Ia, v_Iab, v_Sab);
		if (Ia < epsilon) {
			saveToFile_diadt[s] = 0;
			end = true;
			++outputSize;
			break;
		}
		saveToFile_diadt[s] = Ia / NUM_AGENTS;

#ifdef NORM_SITE_PER_AG
		saveToFile_dildt[s] = il * (N / NUM_AGENTS);
#else
		saveToFile_dildt[s] = il;
#endif
		++outputSize;
	}
	if (end) return;

	//From the 'largerDetailUntil' iteration on, we afford to ignore 'outputGranularity'-size windows of values, so that the saved file does not grow explosively.
	for (uint s = (uint)largerDetailUntil; s < totalSteps; ++s) {
		takeStep(h, Ia, v_Iab, v_Sab);
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
#endif

//void sim::checkinAsInf   (const agent& ag, const node& v) {
void sim::check_in (const agent& ag, const node& v, vector<vector<agent>>& _where) {
	currentNode[ag] = v;
	//vector<agent>& list = iAgents[v];
	vector<agent>& list = _where[v];
	if (list[ELEMS] == (list.size() - 1))
		list.resize(2 * (list.size()));
	uint& lastPos = list[ELEMS];
	++lastPos;
	indexWithinNode[ag] = lastPos;
	list[lastPos] = ag;
}
//void sim::checkinAsSus   (const agent& ag, const node& v) {
//	currentNode[ag] = v;
//	vector<agent>& list = sAgents[v];
//	if (list[ELEMS]  == (list.size() - 1))
//		list.resize(2 * (list.size()));
//	uint& lastPos = list[ELEMS];
//	++lastPos;
//	indexWithinNode[ag] = lastPos;
//	list[lastPos] = ag;
//}
//void sim::checkoutAsSus  (const agent& ag, const node& v) {
void sim::check_out  (const agent& ag, const node& v, vector<vector<agent>>& _where) {
	//vector<agent>& list = sAgents[v];
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
	//checkinAsSus(ag, v);
	check_in(ag, v, sAgents);
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
	if (numI > 0) {
		const vector<uint>& list = iAgents[v];
		for (uint i = 1; i <= numI; ++i) {
			const real delta = EXPTau_aa();
			schedule.emplace(ag, now + delta, action::agInfectAg, list[i], snapshot_a[ag], snapshot_a[list[i]]);
		}
	}
	if (isInfectedSite[v]) {
		const real delta = EXPTau_la();
		schedule.emplace(ag, now + delta, action::siteInfectAg, v, snapshot_a[ag], snapshot_l[v]);
	}
#ifdef PROTECTION_FX
	if (numS == 1)
		graph::Graph::updateHasS(v);
#endif
}
void sim::enterNodeAsInf (const agent& ag, const node& v, const real& now) {
	++snapshot_a[ag];
	//checkinAsInf(ag, v);
	check_in(ag, v, iAgents);
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
#ifdef PROTECTION_FX
	if (numI == 1) {
		graph::Graph::updateHasI(v);
	}
#endif
	const vector<uint>& list = sAgents[v];
	for (uint i = 1; i <= numS; ++i) {
		const real delta = EXPTau_aa();
		schedule.emplace(list[i], now + delta, action::agInfectAg, ag, snapshot_a[list[i]], snapshot_a[ag]);
	}
	if (!isInfectedSite[v]) {
		const real delta = EXPTau_al();
		schedule.emplace(v, now + delta, action::agInfectSite, ag, snapshot_l[v], snapshot_a[ag]);
	}
}
void sim::leaveNodeAsInf (const agent& ag, const node& v, const real& now) {
	++snapshot_a[ag];
	//checkoutAsInf(ag, v);
	check_out(ag, v, iAgents);
	uint&		numI = iInNode[v];				// ----> Number of infected agents currently hosted in v.
	--numI;
#ifdef OCCUPANCY
	stat::Stats::updateInfOccLowered(v, numI + numS, numI, now);
#endif
#ifdef PROTECTION_FX
	if (numI == 0) {
		graph::Graph::updateNoI(v);
	}
#endif
}
void sim::leaveNodeAsSus (const agent& ag, const node& v, const real& now) {
	++snapshot_a[ag];
	//checkoutAsSus(ag, v);
	check_out(ag, v, sAgents);
	uint&		numS = sInNode[v];				// ----> Number of susceptible agents currently hosted in v.
	--numS;
#ifdef OCCUPANCY
	stat::Stats::updateSusOccLowered(v, numI + numS, numS, now);
#endif
#ifdef PROTECTION_FX
	if (numS == 0)
		graph::Graph::updateNoS(v);
#endif
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
	//const vector<uint>& slist = sAgents[v];
	//const vector<uint>& ilist = iAgents[v];
	//for (uint i = slist[0]; i > 0; --i) {
	//	++snapshot_a[slist[i]];
	//}
	//for (uint i = ilist[0]; i > 0; --i) {
	//	++snapshot_a[ilist[i]];
	//}
}

void sim::agFate_fromAg(const agent& ag, const real& now, const uint& infective, const uint& validity_S, const uint& validity_I) {
#ifdef ESTIMATE_PROBS
	++stat::Stats::totalFate;
	if (exposure[ag] > 0) 
		stat::Stats::totalExposition += exposure[ag]; 
#endif
	if (validity_S != snapshot_a[ag] || validity_I != snapshot_a[infective]) 
		return;	// ----> Event became obsolete.
	
	isInfectedAg[ag] = true;
	--saTotal;
	++iaTotal;
#ifdef ESTIMATE_PROBS
	++stat::Stats::totalInfections;
#endif
#ifdef INFECTED_FRACTION
	stat::Stats::bufferizeIFrac(ag, now, "Ia", iaTotal, ilTotal, NUM_AGENTS, OVERLOOK);
#endif
	leaveNodeAsSus(ag, currentNode[ag], now);
	enterNodeAsInf(ag, currentNode[ag], now);
	schedule.emplace(ag, now + EXPGamma_a(), action::recoverAg); // ----> 'Recover' event is scheduled.
}
void sim::agFate_fromSite(const agent& ag, const real& now, const uint& infective, const uint& validity_S, const uint& validity_I) {
#ifdef ESTIMATE_PROBS
	++stat::Stats::totalFate;
	if (exposure[ag] > 0)
		stat::Stats::totalExposition += exposure[ag];
#endif
	if (validity_S != snapshot_a[ag] || validity_I != snapshot_l[infective]) 
		return;	// ----> Event became obsolete.
	
	isInfectedAg[ag] = true;
	--saTotal;
	++iaTotal;
#ifdef ESTIMATE_PROBS
	++stat::Stats::totalInfections;
#endif
#ifdef INFECTED_FRACTION
	stat::Stats::bufferizeIFrac(ag, now, "Ia", iaTotal, ilTotal, NUM_AGENTS, OVERLOOK);
#endif
	leaveNodeAsSus(ag, currentNode[ag], now);
	enterNodeAsInf(ag, currentNode[ag], now);
	schedule.emplace(ag, now + EXPGamma_a(), action::recoverAg); // ----> 'Recover' event is scheduled.
}
void sim::siteFate(const uint& v, const real& now, const uint& infective, const uint& validity_S, const uint& validity_I) {
#ifdef ESTIMATE_PROBS
	++stat::Stats::totalFate;
	if (exposure[ag] > 0)
		stat::Stats::totalExposition += exposure[ag];
#endif
	if (validity_S != snapshot_l[v] || validity_I != snapshot_a[infective]) 
		return;	// ----> Event became obsolete.
	
	isInfectedSite[v] = true;
	++ilTotal;
#ifdef ESTIMATE_PROBS
	++stat::Stats::totalInfections;
#endif
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

	//v_Iab, v_ilb e v_Sab:
	vector<real> v_Iab(Graph::block_prob.size(), 0.0), v_ilb(Graph::block_prob.size(), 0.0), v_Sab(Graph::block_prob.size(), 0.0);

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
		iaTotal = (I_0 > NUM_AGENTS) ? NUM_AGENTS : I_0;
		
		//Initializing the applicable stats:
#ifdef INFECTED_FRACTION
		Stats::initStream(stat::streamType::infFrac);
#endif
#ifdef ESTIMATE_PROBS
		Stats::initStream(Ws, Wi, stat::streamType::probs, SHORT_LABEL, Graph::n, _numAgents, TAU, GAMMA, LAMBDA, T, ROUNDS);
#endif
#ifdef OCCUPANCY
		Stats::initOcc(Graph::n);
		Stats::initStream(Ws, Wi, stat::streamType::occupancy,  SHORT_LABEL, Graph::n, _numAgents, TAU, GAMMA, LAMBDA, T, ROUNDS);
#endif
		Stats::initStream(stat::streamType::avDuration);
		Reporter::startChronometer("\n\n\nRunning scenario " + std::to_string(scenario + 1) + "/" + std::to_string(numScenarios) + "...");
		Reporter::simulationInfo(iaTotal);
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
			for (agent i = 0; i < iaTotal; ++i) {
				isInfectedAg[i] = true;
				schedule.emplace(i, EXPGamma_a(), action::recoverAg);
			}

			// * DISTRIBUTING THE AGENTS ACROSS THE NETWORK *
			//Random distribution of the INFECTED agents:
#ifdef CLIQUE
			for (agent i = 0; i < itotal; ++i)
				enterNodeAsInf(i, randomInt(Graph::n), TIME_ZERO);
#else
			for (agent i = 0; i < iaTotal; ++i) {
				const node v = randomLCCNode();
				enterNodeAsInf(i, v, TIME_ZERO);
				++v_Iab[Graph::g[v].size()];
				++graph::Graph::kb[Graph::g[v].size()];
			}
#endif

			//Random distribution of the SUSCEPTIBLE agents:
#ifdef CLIQUE
			for (agent i = itotal; i < _numAgents; ++i)
				enterNodeAsSus(i, randomInt(Graph::n), TIME_ZERO);
#else
			for (agent i = iaTotal; i < _numAgents; ++i) {
				const node v = randomLCCNode();
				enterNodeAsSus(i, v, TIME_ZERO);
				++v_Sab[Graph::g[v].size()];
				++graph::Graph::kb[Graph::g[v].size()];
			}
#endif
			// * MAIN LOOP *
			bool earlyStop = false;
			real& roundDuration = totalSimTime[round]; 
			double timeLimit = T;
			job j;
			nextJob(j, now);
			while (now < timeLimit) {
				roundDuration += (now - roundDuration);
				switch (j.a) {
				case action::walk:
					walk(j.target, now);
					schedule.emplace(j.target, now + EXPLambda(), action::walk);
					break;
				case action::agInfectAg:
					agFate_fromAg(j.target, now, j.infective, j.validity_S, j.validity_I);
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
			}
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

#ifdef SOLVE_NUMERICALLY
	//Runge-Kutta:
	constexpr uint outputGranularity = 50;
	constexpr uint largerDetailUntil = 100;
	constexpr real stepSize = 0.001;
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
		<< "_Taa"	<< TAU_aa 
		<< "_Tal"	<< TAU_al 
		<< "_Tla"	<< TAU_la 
		<< "_Ga"	<< GAMMA_a 
		<< "_Gl"	<< GAMMA_l 
		<< "_L"		<< LAMBDA 
		<< "_STime" << T 
		<< "_R"		<< ROUNDS;
	baseName = name.str();

	rungeKutta4thOrder(0, v_Iab, v_Sab, v_ilb, T, stepSize, epsilon, saveToFile_diadt, saveToFile_dildt, outputSize, outputGranularity, largerDetailUntil);

	//Saving to file:
	std::ofstream RKdata;
	std::stringstream _ss;
	_ss.precision(4);
	_ss << EXE_DIR << "/stats/Runge-Kutta_" << baseName << ".csv";
	std::string _fileName = _ss.str();
	RKdata.open(_fileName);
	RKdata << "Time\tia\til\n"; // ----> Header
	real _time = 0;
	for (size_t s = 0; s < largerDetailUntil; ++s) {
		RKdata << _time << '\t' << saveToFile_diadt[s] << '\t' << saveToFile_dildt[s] << '\n';
		_time += stepSize;
	}
	for (size_t s = largerDetailUntil; s < outputSize; ++s){
		RKdata << _time << '\t' << saveToFile_diadt[s] << '\t' << saveToFile_dildt[s] << '\n';
		_time += timeIncrement;
	}
	RKdata.close();
#endif //SOLVE_NUMERICALLY
}
