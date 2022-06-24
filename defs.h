#pragma once

#include <stdio.h>
#include <tchar.h>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <queue>
#include <random>		// ----> uniform_real_distribution<>
#include <cassert>
#include <functional>	// ----> Function as an argument

// * DEBUG * 
//#define DEBUG
#ifdef DEBUG
namespace sim {
	constexpr double epsilon = 1.0 / 10e7;
}
#define assertm(exp, msg) assert(((void)msg, exp))
#endif 

// ---------------------------//----------------------------- //

// * DATA TYPES *
#define uint uint32_t
#define ulong uint64_t
#define ushort uint16_t
#define real double
//#define real float
typedef uint32_t agent;
namespace graph {
	typedef uint32_t node;
}

// ---------------------------//----------------------------- //

// * BASE DIRECTORY AND FILE NAME * 
#define EXE_DIR "C:/Users/Ronald/source/repos/randomWalkTimeFrame/x64/Release"
namespace sim {
	static std::string baseName;
	static real beta_a	= 0;														// ----> Force of infection from an I-agent to an S-agent
	static real beta_al	= 0;														// ----> Force of infection from an I-agent to a site
	static real beta_la	= 0;														// ----> Force of infection from a site to an I-agent
}

// ---------------------------//----------------------------- //


// * NETWORK *
//#define CLIQUE
//#define READ_NTWK_FROM_FILE
#define GNP
//#define STAR

#ifdef STAR
static constexpr uint N = 12008;
#define SOURCE_FILE "", 12008
#define NWTK_LABEL "Star"
#define SHORT_LABEL "star"
#endif

#ifdef CLIQUE
//static constexpr uint N = 12008;										// ----> Network size
//#define NWTK_LABEL "Clique"
//#define SHORT_LABEL "CL"
#endif

#ifdef GNP
static constexpr uint N = 12008;
#define SOURCE_FILE "", 12008
#define NWTK_LABEL "Gnp"
#define SHORT_LABEL "gnp"
#endif
#ifdef READ_NTWK_FROM_FILE
//static constexpr uint N = 55;										// ----> Network size
//#define SOURCE_FILE std::string(std::string(EXE_DIR) + "\\redes\\grafoDeTestes.txt"), 55
//#define NWTK_LABEL "Ronald"
//#define SHORT_LABEL "ron"
//
//static constexpr uint N = 4039;										// ----> Network size
//#define SOURCE_FILE std::string(std::string(EXE_DIR) + "\\redes\\facebook_combined.txt"), 4039
//#define NWTK_LABEL "Fb"
//#define SHORT_LABEL "fb"
//
//static constexpr uint N = 12008;										// ----> Network size
//#define SOURCE_FILE std::string(std::string(EXE_DIR) + "\\redes\\CA-HepPh.txt"), 12008
//#define NWTK_LABEL "HepPh"
//#define SHORT_LABEL "hep"
//
static constexpr uint N = 12008;										// ----> Network size
#define SOURCE_FILE std::string(std::string(EXE_DIR) + "\\synthetic\\BA.txt"), 12008
#define NWTK_LABEL "BA"
#define SHORT_LABEL "BA"
//
//static constexpr uint N = 15233;										// ----> Network size
//#define SOURCE_FILE std::string(std::string(EXE_DIR) + "\\redes\\netHEPT.txt"), 15233
//#define NWTK_LABEL "net"
//#define SHORT_LABEL "net"
//
//static constexpr uint N = 18772;										// ----> Network size
//#define SOURCE_FILE std::string(std::string(EXE_DIR) + "\\redes\\CA-AstroPh.txt"), 18772
//#define NWTK_LABEL "AstroPh"
//#define SHORT_LABEL "astro"
//
//static constexpr uint N = 23133;										// ----> Network size
//#define SOURCE_FILE std::string(std::string(EXE_DIR) + "\\redes\\CA-CondMat.txt"), 23133
//#define NWTK_LABEL "CondMat"
//#define SHORT_LABEL "cmat"
//
//static constexpr uint N = 36692;										// ----> Network size
//#define SOURCE_FILE std::string(std::string(EXE_DIR) + "\\redes\\Email-Enron.txt"), 36692
//#define NWTK_LABEL "Enron"
//#define SHORT_LABEL "enron"
//
//static constexpr uint N = 58228;										// ----> Network size
//#define SOURCE_FILE std::string(std::string(EXE_DIR) + "\\redes\\Brightkite_edges.txt"), 58228
//#define NWTK_LABEL "Brightkite"
//#define SHORT_LABEL "bk"
//
//static constexpr uint N = 196591;										// ----> Network size
//#define SOURCE_FILE std::string(std::string(EXE_DIR) + "\\redes\\Gowalla_edges.txt"), 196591
//#define NWTK_LABEL "Gowalla"
//#define SHORT_LABEL "gw"
//
//static constexpr uint N = 317080;										// ----> Network size
//#define SOURCE_FILE std::string(std::string(EXE_DIR) + "\\redes\\com-dblp.ungraph.txt"), 317080
//#define NWTK_LABEL "dblp"
//#define SHORT_LABEL "dblp"
//
//static constexpr uint N = 334863;										// ----> Network size
//#define SOURCE_FILE std::string(std::string(EXE_DIR) + "\\redes\\com-amazon.ungraph.txt"), 334863
//#define NWTK_LABEL "amazon"
//#define SHORT_LABEL "amz"
//
//static constexpr uint N = 1696415;										// ----> Network size
//#define SOURCE_FILE std::string(std::string(EXE_DIR) + "\\redes\\as-skitter.txt"), 1696415
//#define NWTK_LABEL "as-skitter"
//#define SHORT_LABEL "as"
#endif
// ---------------------------//----------------------------- //

// * AGENTS' BEHAVIOR *
#define AUTO_RELATION					// ----> Gives agents the option of staying at their current node upon their walk event. If not enabled, agents will necessarily change their current node when their walk event is processed.
//#define PROTECTION_FX
#ifndef PROTECTION_FX
static constexpr real Ws = 1.0;	// !DO NOT CHANGE THIS LINE! To set Ws to 1.0 here means "no protection effect", which is the desired behaviour when the pre-processor macro "PROTECTION_FX" is not defined.
static constexpr real Wi = 1.0;	// !DO NOT CHANGE THIS LINE! To set Wi to 1.0 here means "no protection effect", which is the desired behaviour when the pre-processor macro "PROTECTION_FX" is not defined.
#endif

// ---------------------------//----------------------------- //

// * STATS * 
//#define OCCUPANCY
#define INFECTED_FRACTION
//#define ESTIMATE_PROBS
#define i_t_FROM_MODEL
#define SOLVE_NUMERICALLY
#define NORM_SITE_PER_AG

// ---------------------------//----------------------------- //

// * REPORTER * 
//#define MEASURE_ROUND_EXE_TIME

// ---------------------------//----------------------------- //

namespace sim{		// ----> Simulator's namespace.

// * SIMULATION PARAMETERS *
#ifdef PROTECTION_FX
static constexpr real Ws  = 1.0;										// ----> Susceptible-agents' tolerance to enter nodes that contain infected agents, such that 0 <= Ws <= 1. This is the "s-protection-effect" single parameter.
static constexpr real Wi  = 1.0;										// ----> Infected-agents' tolerance to enter nodes that contain susceptible agents, such that 0 <= Wi <= 1. This is the "i-protection-effect" single parameter.

//#define PROPORTIONAL													// ----> Promotes risk-tolerance proportional to the current number of infectives, according to the function w(i) = (1-i)^r, for some "rejection force" r.
#ifdef PROPORTIONAL
//If defined, we consider that the risk-tolerance is proportional to the number of infectives: w(i) = (1-i)^r, for some "rejection force" r.
static constexpr real _r  = 1000.0;		// Rejection force.
#endif //PROPORTIONAL
#endif //PROTECTION_FX

static constexpr uint T					= 2000;						// ----> Simulation time.
static constexpr uint NUM_AGENTS		= 200;							// ----> Total number of agents in a simulation.
static constexpr uint STARTING_NUM_AG	= 1000000;							
static constexpr uint GRAN_NUM_AG		= 1;							
static constexpr uint ROUNDS			= 1;							// ----> Number of simulation runs for a given setup. 
static constexpr real TAU_aa			= 4.0;							// ----> Admits two different views: 1) "Resistance to exposure": the larger, the harder it gets to infect an exposed, susceptible agent; 2) "Propagator's 'Infectivity'": in this case, SMALLER values yield LARGER transmission probability. Parameter of an exponentially-distributed random-number generator.
static constexpr real TAU_al			= 0.000001;							// ----> Admits two different views: 1) "Resistance to exposure": the larger, the harder it gets to infect an exposed, susceptible agent; 2) "Propagator's 'Infectivity'": in this case, SMALLER values yield LARGER transmission probability. Parameter of an exponentially-distributed random-number generator.
static constexpr real TAU_la			= 0.000001;							// ----> Admits two different views: 1) "Resistance to exposure": the larger, the harder it gets to infect an exposed, susceptible agent; 2) "Propagator's 'Infectivity'": in this case, SMALLER values yield LARGER transmission probability. Parameter of an exponentially-distributed random-number generator.
static constexpr real GAMMA_a			= 0.01;							// ----> Recovery rate. The higher, the faster. Parameter of an exponentially-distributed random-number generator.
static constexpr real GAMMA_l			= 2000.0;							// ----> Recovery rate. The higher, the faster. Parameter of an exponentially-distributed random-number generator.
static constexpr real LAMBDA			= 4.0;							// ----> Walking speed. The higher, the faster. Parameter of an exponentially-distributed random-number generator.
static constexpr real FRAC_AG_INFECTED	= 0.5;							// ----> Fraction of AGENTS initially infected (i.e. when the simulation starts).
static constexpr real FRAC_ST_INFECTED	= 0.0;							// ----> Fraction of SITES initially infected (i.e. when the simulation starts).
static constexpr uint ABS_INFECTED		= 0;							// ----> Absolute number of agents initially infected (i.e. when the simulation starts). This value is used whenever set to any value > 0, in which case it overrides 'FRAC_AG_INFECTED'. To use 'FRAC_AG_INFECTED' instead, set 'ABS_INFECTED = 0'.

// * CONSTANTS * 
static constexpr uint I_0 = (ABS_INFECTED > 0) ? ABS_INFECTED : (uint)((real)NUM_AGENTS * FRAC_AG_INFECTED);
static constexpr long real i_0 = (real)I_0 / NUM_AGENTS;
static constexpr long real meetingRate = 2.0 * LAMBDA / N;
//static constexpr long real SIGMA_aa = (TAU_aa / (2.0 * LAMBDA + TAU_aa + GAMMA_a));
//static constexpr long real SIGMA_al = (TAU_al / (LAMBDA + TAU_al + GAMMA_a));
//static constexpr long real SIGMA_la = (TAU_la / (LAMBDA + TAU_la + GAMMA_l));

static constexpr long real SIGMA_aa = (TAU_aa / (2.0 * LAMBDA + TAU_aa));
static constexpr long real SIGMA_al = (TAU_al / (LAMBDA + TAU_al));
static constexpr long real SIGMA_la = (TAU_la / (LAMBDA + TAU_la));

//static constexpr long real BETA = meetingRate * SIGMA_aa * NUM_AGENTS;
//static constexpr long real B_MINUS_G = BETA - GAMMA;
//static constexpr long real G_OVER_B = GAMMA / BETA;
//static constexpr long real _1_MINUS_G_OVER_B = 1.0 - G_OVER_B;
//static constexpr long real C = i_0 / (1 - i_0 - G_OVER_B);
static constexpr long real NEG_RECIPR_LAMBDA	= -(1.0 / LAMBDA);				// ----> Preprocessing. The negative reciprocal of LAMBDA, to be used by the exponential random-number generator.
static constexpr long real NEG_RECIPR_GAMMA_a	= -(1.0 / GAMMA_a);				// ----> Preprocessing. The negative reciprocal of GAMMA, to be used by the exponential random-number generator.
static constexpr long real NEG_RECIPR_GAMMA_l	= -(1.0 / GAMMA_l);				// ----> Preprocessing. The negative reciprocal of GAMMA, to be used by the exponential random-number generator.
static constexpr long real NEG_RECIPR_TAU_aa	= -(1.0 / TAU_aa);				// ----> Preprocessing. The negative reciprocal of TAU, to be used by the exponential random-number generator.
static constexpr long real NEG_RECIPR_TAU_al	= -(1.0 / TAU_al);				// ----> Preprocessing. The negative reciprocal of TAU, to be used by the exponential random-number generator.
static constexpr long real NEG_RECIPR_TAU_la	= -(1.0 / TAU_la);				// ----> Preprocessing. The negative reciprocal of TAU, to be used by the exponential random-number generator.
static constexpr uint ELEMS = 0;												// ----> Zero here means "the first position of the container". Both 'sAgents' and 'iAgents' lists store their current number of elements in their respective first positions. 
static constexpr long real TIME_ZERO = 0;

static constexpr uint K = NUM_AGENTS;
//static constexpr long real a = Ws;
//static constexpr long real b = Wi;
//static constexpr long real meetingRate_pfx = LAMBDA * (a / (i_0 * (a - 1) + N) + b / ((1 - i_0) * K * (b - 1) + N));
//static constexpr long real meetingRate_pfx_approx = LAMBDA * ((a+b)/N);
//static constexpr long real BETA_pfx = SIGMA_aa * K * meetingRate_pfx;
//static constexpr long real BETA_pfx_approx = SIGMA_aa * K * meetingRate_pfx_approx;
//static constexpr real B_MINUS_G_pfx = BETA_pfx_approx - GAMMA;
//static constexpr real G_OVER_B_pfx = GAMMA / BETA_pfx_approx;
//static constexpr real _1_MINUS_G_OVER_B_pfx = 1.0 - G_OVER_B_pfx;
//static constexpr real C_pfx = i_0 / (1 - i_0 - G_OVER_B_pfx);
//static constexpr long real Ro_pfx = BETA_pfx / GAMMA;
//static constexpr long real Ro_pfx_approx = BETA_pfx_approx / GAMMA;
//static constexpr long real i_inf_pfx = 1 - (GAMMA / BETA_pfx);
//static constexpr long real i_inf_pfx_approx = 1 - (GAMMA / BETA_pfx_approx);

#ifdef INFECTED_FRACTION
static constexpr real OVERLOOK_RATE = 1.0;									// ----> Depending on the initial settings, the simulation may generate a firehose of data, which in turn becomes highly inconvenient (or even prohibitive) for being written to a file (in terms of either space and time). For such cases, it is strongly adviseable to purposely overlook a portion of the events generated per time unit.
static constexpr uint OVERLOOK = (uint)((NUM_AGENTS * OVERLOOK_RATE)/1);	// ----> Depending on the initial settings, the simulation may generate a firehose of data, which in turn becomes highly inconvenient (or even prohibitive) for being written to a file (in terms of either space and time). For such cases, it is strongly adviseable to purposely overlook a portion of the events generated per time unit.
#endif

#ifdef CLIQUE
static const uint LIST_INI_SZ = (uint)(round(std::max((real)2.0, (real)NUM_AGENTS / (3 * N)))); // ----> Initial size of both 'sAgents' and 'iAgents' lists. Every time a node's list become full, its size gets doubled. Although arbitrary, the initial value provided here aims at reducing both the number of times a doubling operation is required and the vector's final size.
#else
const uint LIST_INI_SZ = (uint)(round(std::max((real)2.0, (real)NUM_AGENTS / (3 * N)))); // ----> Initial size of both 'sAgents' and 'iAgents' lists. Every time a node's list become full, its size gets doubled. Although arbitrary, the initial value provided here aims at reducing both the number of times a doubling operation is required and the vector's final size.
#endif

// * SIMULATION UTILS *
//#ifdef i_t_FROM_MODEL
//real i_t_2ndMmt(const real& t);
//real i_t_2ndMmt_naive(const real& t);
//real i_t_pfx(const real& t);
//#endif

#ifdef SOLVE_NUMERICALLY
static constexpr long real crowdFactor = std::min(2.0, std::max((real)K / N, 1.0));
real diadt(const real& ia, const double& sumSbIb);
//real diadt(const real& ia, const real& il);
real dildt(const real& ia, const real& il);
real diabdt(const real& Ia, const real& Iab, const real& Sab, const uint& block);
real dsabdt(const real& Ia, const real& Iab, const real& Sab, const uint& block);
void step(const real& h, real& Ia, std::vector<real>& v_Iab, std::vector<real>& v_Sab);

//Returns 'true' if all the range of steps is performed; returns 'false' if the value for a new step happens to be 0, which means that the numerical solution has been completed prior to taking all the steps predicted at the input. Important: the function considers the interval [from, to[ i.e. 'from' is included but 'to' is excluded. The first index is "from" and the last is "to - 1".
//bool takeSteps(const uint& from, const uint& to);

real dilbdt(const real& ia, const real& il, const real& iab, const real& ilb, const uint& block);
void update_Ia(real& Ia, const std::vector<real>& v_Iab);
void update_Ia(real& Ia, const std::vector<real>& v_Iab, const std::vector<real>& base, const double& fraction = 1.0);
void lookAhead(const real& h, real& Ia, const std::vector<real>& v_Iab, const std::vector<real>& v_Sab, std::vector<real>& target);
void lookAhead(const real& h, real& Ia, const std::vector<real>& v_Iab, const std::vector<real>& v_Sab, std::vector<real>& target, std::vector<real>& base, const double& fraction = 1.0);
void rungeKutta4thOrder(const real& t0, std::vector<real>& v_Iab, std::vector<real>& v_Sab, std::vector<real>& v_ilb, const real& t, const real& h, const real& epsilon, std::vector<real>& saveToFile_diadt, std::vector<real>& saveToFile_dildt, uint& outputSize, const uint& outputGranularity = 50, const real& largerDetailUntil = 1000);
#endif

} //namespace sim

// ---------------------------//----------------------------- //