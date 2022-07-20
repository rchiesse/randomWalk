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
	static real nT;		// ----> Normalized TAU
	static real nG;		// ----> Normalized GAMMA
	static real nL;		// ----> Normalized LAMBDA
}

// ---------------------------//----------------------------- //


// * NETWORK *
#define CLIQUE
//#define READ_NTWK_FROM_FILE
//#define GNP
//#define STAR

#ifdef STAR
static constexpr uint N = 12008;
#define SOURCE_FILE "", 12008
#define NWTK_LABEL "Star"
#define SHORT_LABEL "star"
#endif

#ifdef CLIQUE
static constexpr uint N = 200;										// ----> Network size
#define NWTK_LABEL "Clique"
#define SHORT_LABEL "CL"
#endif

#ifdef GNP
//static constexpr uint N = 12008;
static constexpr uint N = 200;
#define SOURCE_FILE "", 12008
#define NWTK_LABEL "Gnp"
#define SHORT_LABEL "gnp"
#endif
#ifdef READ_NTWK_FROM_FILE
//static constexpr uint N = 55;										// ----> Network size
//#define SOURCE_FILE std::string(".\\redes\\grafoDeTestes.txt"), 55
//#define NWTK_LABEL "Ronald"
//#define SHORT_LABEL "ron"
//
//static constexpr uint N = 4039;										// ----> Network size
//#define SOURCE_FILE std::string(".\\redes\\facebook_combined.txt"), 4039
//#define NWTK_LABEL "Fb"
//#define SHORT_LABEL "fb"
//
//static constexpr uint N = 12008;										// ----> Network size
//#define SOURCE_FILE std::string(".\\redes\\CA-HepPh.txt"), 12008
//#define NWTK_LABEL "HepPh"
//#define SHORT_LABEL "hep"
//
static constexpr uint N = 12008;										// ----> Network size
#define SOURCE_FILE std::string(".\\synthetic\\BA.txt"), 12008
#define NWTK_LABEL "BA"
#define SHORT_LABEL "BA"
//
//static constexpr uint N = 15233;										// ----> Network size
//#define SOURCE_FILE std::string(".\\redes\\netHEPT.txt"), 15233
//#define NWTK_LABEL "net"
//#define SHORT_LABEL "net"
//
//static constexpr uint N = 18772;										// ----> Network size
//#define SOURCE_FILE std::string(".\\redes\\CA-AstroPh.txt"), 18772
//#define NWTK_LABEL "AstroPh"
//#define SHORT_LABEL "astro"
//
//static constexpr uint N = 23133;										// ----> Network size
//#define SOURCE_FILE std::string(".\\redes\\CA-CondMat.txt"), 23133
//#define NWTK_LABEL "CondMat"
//#define SHORT_LABEL "cmat"
//
//static constexpr uint N = 36692;										// ----> Network size
//#define SOURCE_FILE std::string(".\\redes\\Email-Enron.txt"), 36692
//#define NWTK_LABEL "Enron"
//#define SHORT_LABEL "enron"
//
//static constexpr uint N = 58228;										// ----> Network size
//#define SOURCE_FILE std::string(".\\redes\\Brightkite_edges.txt"), 58228
//#define NWTK_LABEL "Brightkite"
//#define SHORT_LABEL "bk"
//
//static constexpr uint N = 196591;										// ----> Network size
//#define SOURCE_FILE std::string(".\\redes\\Gowalla_edges.txt"), 196591
//#define NWTK_LABEL "Gowalla"
//#define SHORT_LABEL "gw"
//
//static constexpr uint N = 317080;										// ----> Network size
//#define SOURCE_FILE std::string(".\\redes\\com-dblp.ungraph.txt"), 317080
//#define NWTK_LABEL "dblp"
//#define SHORT_LABEL "dblp"
//
//static constexpr uint N = 334863;										// ----> Network size
//#define SOURCE_FILE std::string(".\\redes\\com-amazon.ungraph.txt"), 334863
//#define NWTK_LABEL "amazon"
//#define SHORT_LABEL "amz"
//
//static constexpr uint N = 1696415;										// ----> Network size
//#define SOURCE_FILE std::string(".\\redes\\as-skitter.txt"), 1696415
//#define NWTK_LABEL "as-skitter"
//#define SHORT_LABEL "as"
#endif
// ---------------------------//----------------------------- //


// * AGENTS' BEHAVIOR *
#define AUTO_RELATION					// ----> Gives an agent the option of staying at its current node upon every other walk event. If not enabled, agents will necessarily move to a different node when their walk event is processed.
//#define PROTECTION_FX
#ifndef PROTECTION_FX
static constexpr real Ws = 1.0;	// !DO NOT CHANGE THIS LINE! To set Ws to 1.0 here means "no protection effect", which is the desired behaviour when the pre-processor macro "PROTECTION_FX" is not defined.
static constexpr real Wi = 1.0;	// !DO NOT CHANGE THIS LINE! To set Wi to 1.0 here means "no protection effect", which is the desired behaviour when the pre-processor macro "PROTECTION_FX" is not defined.
#endif

// ---------------------------//----------------------------- //


// * STATS * 
#define INFECTED_FRACTION
//#define OCCUPANCY
//#define ESTIMATE_PROBS

// ---------------------------//----------------------------- //


// * REPORTER * 
//#define MEASURE_ROUND_EXE_TIME

// ---------------------------//----------------------------- //


// * SIMULATION PARAMETERS *
namespace sim{		// ----> Simulator's namespace.
#ifdef PROTECTION_FX
static constexpr real Ws  = 1.0;										// ----> Susceptible-agents' tolerance to enter nodes that contain infected agents, such that 0 <= Ws <= 1. This is the "s-protection-effect" single parameter.
static constexpr real Wi  = 1.0;										// ----> Infected-agents' tolerance to enter nodes that contain susceptible agents, such that 0 <= Wi <= 1. This is the "i-protection-effect" single parameter.

//#define PROPORTIONAL													// ----> Promotes risk-tolerance proportional to the current number of infectives, according to the function w(i) = (1-i)^r, for some "rejection force" r.
#ifdef PROPORTIONAL
//If defined, we consider that the risk-tolerance is proportional to the number of infectives: w(i) = (1-i)^r, for some "rejection force" r.
static constexpr real _r  = 1000.0;		// Rejection force.
#endif //PROPORTIONAL
#endif //PROTECTION_FX

// Input parameters
static real T;															// ----> Simulation time.
static uint NUM_AGENTS;													// ----> Total number of agents in a simulation.
static uint STARTING_NUM_AG;											
static uint GRAN_NUM_AG;	
static uint ROUNDS;														// ----> Number of simulation runs for a given setup. 
static real TAU_aa;														// ----> Agent-to-agent transmissibility rate.
static real TAU_al;														// ----> Agent-to-location transmissibility rate.
static real TAU_la;														// ----> Location-to-agent transmissibility rate.
static real GAMMA_a;													// ----> Recovery rate. 
static real GAMMA_l;													// ----> Recovery rate. 
static real LAMBDA;														// ----> Walk rate. 
static real FRAC_AG_INFECTED;											// ----> Fraction of AGENTS initially infected (i.e. when the simulation starts).
static real FRAC_ST_INFECTED;											// ----> Fraction of SITES initially infected (i.e. when the simulation starts).
static uint ABS_INFECTED;												// ----> Absolute number of agents initially infected (i.e. when the simulation starts). This value is used whenever set to any value > 0, in which case it overrides 'FRAC_AG_INFECTED'. To use 'FRAC_AG_INFECTED' instead, set 'ABS_INFECTED = 0'.

// Auxiliary constants
static const uint I_0 = (ABS_INFECTED > 0) ? ABS_INFECTED : (uint)((real)NUM_AGENTS * FRAC_AG_INFECTED);
static const long real i_0 = (real)I_0 / NUM_AGENTS;
static const long real meetingRate = 2.0 * LAMBDA / N;
static const long real SIGMA_aa = (TAU_aa / (2.0 * LAMBDA + TAU_aa));
static const long real SIGMA_al = (TAU_al / (LAMBDA + TAU_al));
static const long real SIGMA_la = (TAU_la / (LAMBDA + TAU_la));
static const long real EULER = 0.57721566490153286060651209008240243104215933593992;	// ----> The Euler–Mascheroni constant.
static const long real NEG_RECIPR_LAMBDA	= -(1.0 / LAMBDA);				// ----> Preprocessing. The negative reciprocal of LAMBDA, to be used by the exponential random-number generator.
static const long real NEG_RECIPR_GAMMA_a	= -(1.0 / GAMMA_a);				// ----> Preprocessing. The negative reciprocal of GAMMA, to be used by the exponential random-number generator.
static const long real NEG_RECIPR_GAMMA_l	= -(1.0 / GAMMA_l);				// ----> Preprocessing. The negative reciprocal of GAMMA, to be used by the exponential random-number generator.
static const long real NEG_RECIPR_TAU_aa	= -(1.0 / TAU_aa);				// ----> Preprocessing. The negative reciprocal of TAU, to be used by the exponential random-number generator.
static const long real NEG_RECIPR_TAU_al	= -(1.0 / TAU_al);				// ----> Preprocessing. The negative reciprocal of TAU, to be used by the exponential random-number generator.
static const long real NEG_RECIPR_TAU_la	= -(1.0 / TAU_la);				// ----> Preprocessing. The negative reciprocal of TAU, to be used by the exponential random-number generator.
static const uint ELEMS = 0;												// ----> Zero here means "the first position of the container". Both 'sAgents' and 'iAgents' lists store their current number of elements in their respective first positions. 
static const long real TIME_ZERO = 0;
static const uint K = NUM_AGENTS;
#ifdef CLIQUE
static const uint LIST_INI_SZ = (uint)(round(std::max((real)2.0, (real)NUM_AGENTS / (3 * N)))); // ----> Initial size of both 'sAgents' and 'iAgents' lists. Every time a node's list become full, its size gets doubled. Although arbitrary, the initial value provided here aims at reducing both the number of times a doubling operation is required and the vector's final size.
#else
const uint LIST_INI_SZ = (uint)(round(std::max((real)2.0, (real)NUM_AGENTS / (3 * N)))); // ----> Initial size of both 'sAgents' and 'iAgents' lists. Every time a node's list become full, its size gets doubled. Although arbitrary, the initial value provided here aims at reducing both the number of times a doubling operation is required and the vector's final size.
#endif

// Output control:
#ifdef INFECTED_FRACTION
static constexpr real OVERLOOK_RATE = 1.0;									// ----> Depending on the initial settings, the simulation may generate a firehose of data, which in turn becomes highly inconvenient (or even prohibitive) for being written to a file (in terms of either space and time). For such cases, it is strongly adviseable to purposely overlook a portion of the events generated per time unit.
static const uint OVERLOOK = (uint)((NUM_AGENTS * OVERLOOK_RATE)/1);
#endif

// ---------------------------//----------------------------- //


// * NUMERICAL SOLUTION *
#define SOLVE_NUMERICALLY
//#define ONLY_NUMERIC
#define NORM_SITE_PER_AG
#define PER_BLOCK		// ----> If defined, then the numerical solution is based on 2 equations per degree-block. If otherwise, then a fine-grained system of 2 equations per NODE is solved (computationally expensive).

#ifdef SOLVE_NUMERICALLY
//static constexpr long real crowdFactor = std::min(2.0, std::max((real)K / N, 1.0));
real diadt(const real& ia, const double& sumSbIb);
real dildt(const real& ia, const real& il);

#ifdef CLIQUE
#ifdef PER_BLOCK
real diabdt(const real& Ia, const real& Sa);
real dsabdt(const real& Ia, const real& Sa);
#else //PER_BLOCK
real divbdt(const real& Ia, const real& Iv, const real& Sv);
real dsvbdt(const real& Ia, const real& Iv, const real& Sv);
#endif //PER_BLOCK
void step(const real& h, real& Ia, real& Sa);
//real dilbdt(const real& ia, const real& il, const real& iab, const real& ilb, const uint& block);
void lookAhead(const real& h, real& Ia, real& Sa, std::vector<real>& target);
void lookAhead(const real& h, real& Ia, real& Sa, std::vector<real>& target, std::vector<real>& base, const double& fraction = 1.0);
void rungeKutta4thOrder(const real& t0, real& Ia, real& Sa, const real& t, const real& h, const real& epsilon, std::vector<real>& saveToFile_diadt, std::vector<real>& saveToFile_dildt, uint& outputSize, const uint& outputGranularity = 50, const real& largerDetailUntil = 1000);

#else //CLIQUE
#ifdef PER_BLOCK
real diabdt(const real& Ia, const real& Iab, const real& Sab, const uint& block);
real dsabdt(const real& Ia, const real& Iab, const real& Sab, const uint& block);
#else //PER_BLOCK
real divbdt(const real& Ia, const real& Iv, const real& Sv, const uint& block);
real dsvbdt(const real& Ia, const real& Iv, const real& Sv, const uint& block);
#endif //PER_BLOCK
void step(const real& h, real& Ia, std::vector<real>& v_Iab, std::vector<real>& v_Sab);
real dilbdt(const real& ia, const real& il, const real& iab, const real& ilb, const uint& block);
void update_Ia(real& Ia, const std::vector<real>& v_Iab);
void update_Ia(real& Ia, const std::vector<real>& v_Iab, const std::vector<real>& base, const double& fraction = 1.0);
void lookAhead(const real& h, real& Ia, const std::vector<real>& v_Iab, const std::vector<real>& v_Sab, std::vector<real>& target);
void lookAhead(const real& h, real& Ia, const std::vector<real>& v_Iab, const std::vector<real>& v_Sab, std::vector<real>& target, std::vector<real>& base, const double& fraction = 1.0);
void rungeKutta4thOrder(const real& t0, std::vector<real>& v_Iab, std::vector<real>& v_Sab, std::vector<real>& v_ilb, const real& t, const real& h, const real& epsilon, std::vector<real>& saveToFile_diadt, std::vector<real>& saveToFile_dildt, uint& outputSize, const uint& outputGranularity = 50, const real& largerDetailUntil = 1000);
#endif //CLIQUE
#endif //SOLVE_NUMERICALLY

} //namespace sim

// ---------------------------//----------------------------- //