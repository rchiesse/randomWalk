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
#include <random>		// uniform_real_distribution<>
#include <cassert>

// * DEBUG * 
//#define assertm(exp, msg) assert(((void)msg, exp))
//#define DEBUG

// ---------------------------//----------------------------- //

// * BASE DIRECTORY * 
#define EXE_DIR "C:/Users/Ronald/source/repos/randomWalkTimeFrame/x64/Release"

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

// * NETWORK *
//#define CLIQUE
//#define READ_NTWK_FROM_FILE
#define GNP

// ---------------------------//----------------------------- //

// * AGENTS' BEHAVIOR *
#define AUTO_RELATION					// ----> Gives agents the option of staying at their current node upon their walk event. If not enabled, agents will necessarily change their current node when their walk event is processed.
#define PROTECTION_FX
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
//#define SOLVE_NUMERICALLY

// ---------------------------//----------------------------- //

// * REPORTER * 
//#define MEASURE_ROUND_EXE_TIME

// ---------------------------//----------------------------- //

namespace sim{		// ----> Simulator's namespace.

// * SIMULATION PARAMETERS *
static constexpr uint N = 12008;										// ----> Network size
#ifdef PROTECTION_FX
static constexpr real Ws  = 1.0;										// ----> Susceptible-agents' tolerance to enter nodes that contain infected agents, such that 0 <= Ws <= 1. This is the "s-protection-effect" single parameter.
static constexpr real Wi  = 1.0;										// ----> Infected-agents' tolerance to enter nodes that contain susceptible agents, such that 0 <= Wi <= 1. This is the "i-protection-effect" single parameter.

//#define PROPORTIONAL													// ----> Promotes risk-tolerance proportional to the current number of infectives, according to the function w(i) = (1-i)^r, for some "rejection force" r.
#ifdef PROPORTIONAL
//If defined, we consider that the risk-tolerance is proportional to the number of infectives: w(i) = (1-i)^r, for some "rejection force" r.
static constexpr real _r  = 1000.0;		// Rejection force.
#endif //PROPORTIONAL
#endif //PROTECTION_FX

static constexpr uint T					= 20000;						// ----> Simulation time.
static constexpr uint NUM_AGENTS		= 800;							// ----> Total number of agents in a simulation.
static constexpr uint STARTING_NUM_AG	= 1000000;							
static constexpr uint GRAN_NUM_AG		= 1;							
static constexpr uint ROUNDS			= 1;							// ----> Number of simulation runs for a given setup. 
static constexpr real TAU				= 1.0;							// ----> Admits two different views: 1) "Resistance to exposure": the larger, the harder it gets to infect an exposed, susceptible agent; 2) "Propagator's 'Infectivity'": in this case, SMALLER values yield LARGER transmission probability. Parameter of an exponentially-distributed random-number generator.
static constexpr real GAMMA				= 0.0131;						// ----> Recovery rate. The higher, the faster. Parameter of an exponentially-distributed random-number generator.
static constexpr real LAMBDA			= 1.0;							// ----> Walking speed. The higher, the faster. Parameter of an exponentially-distributed random-number generator.
static constexpr real FRAC_AG_INFECTED	= 0.5;							// ----> Fraction of agents initially infected (i.e. when the simulation starts).
static constexpr uint ABS_INFECTED		= 0;							// ----> Absolute number of agents initially infected (i.e. when the simulation starts). This value is used whenever set to any value > 0, in which case it overrides 'FRAC_AG_INFECTED'. To use 'FRAC_AG_INFECTED' instead, set 'ABS_INFECTED = 0'.

// * CONSTANTS * 
static constexpr uint I_0 = (ABS_INFECTED > 0) ? ABS_INFECTED : (uint)((real)NUM_AGENTS * FRAC_AG_INFECTED);
static constexpr long real i_0 = (real)I_0 / NUM_AGENTS;
static constexpr long real meetingRate = 2.0 * LAMBDA / N;
static constexpr long real infectionProb = (TAU / (2.0 * LAMBDA + TAU));
static constexpr long real BETA = meetingRate * infectionProb * NUM_AGENTS;
static constexpr long real B_MINUS_G = BETA - GAMMA;
static constexpr long real G_OVER_B = GAMMA / BETA;
static constexpr long real _1_MINUS_G_OVER_B = 1.0 - G_OVER_B;
static constexpr long real C = i_0 / (1 - i_0 - G_OVER_B);
static constexpr long real NEG_RECIPR_LAMBDA = -(1.0 / LAMBDA);				// ----> Preprocessing. The negative reciprocal of LAMBDA, to be used by the exponential random-number generator.
static constexpr long real NEG_RECIPR_GAMMA	= -(1.0 / GAMMA);				// ----> Preprocessing. The negative reciprocal of GAMMA, to be used by the exponential random-number generator.
static constexpr long real NEG_RECIPR_TAU	= -(1.0 / TAU);					// ----> Preprocessing. The negative reciprocal of TAU, to be used by the exponential random-number generator.
static constexpr uint ELEMS = 0;											// ----> Zero here means "the first position of the container". Both 'sAgents' and 'iAgents' lists store their current number of elements in their respective first positions. 
static constexpr long real TIME_ZERO = 0;

static constexpr uint K = NUM_AGENTS;
static constexpr long real a = Ws;
static constexpr long real b = Wi;
#ifdef PROTECTION_FX


//static constexpr long real C10 = (GAMMA * K) * (a - a * b + b - 1);
//static constexpr long real C11 = (GAMMA * K) * (a * b - a - b + 1) - (GAMMA * N) * (b * K + K + a - 1);
//static constexpr long real C12 = (GAMMA * N) * (N - K + b * K);
//static constexpr long real C13 = -C5 - C10;	//a
//static constexpr long real C14 = C6 - C11;	//b
//static constexpr long real C15 = C4 - C12;	//c
//static constexpr long real p_ = C14/C13;	//p
//static constexpr long real q_ = C15/C13;	//q
#endif

static constexpr long real meetingRate_pfx = LAMBDA * (a / (i_0 * (a - 1) + N) + b / ((1 - i_0) * K * (b - 1) + N));
static constexpr long real meetingRate_pfx_approx = LAMBDA * ((a+b)/N);
static constexpr long real BETA_pfx = infectionProb * K * meetingRate_pfx;
static constexpr long real BETA_pfx_approx = infectionProb * K * meetingRate_pfx_approx;
static constexpr real B_MINUS_G_pfx = BETA_pfx_approx - GAMMA;
static constexpr real G_OVER_B_pfx = GAMMA / BETA_pfx_approx;
static constexpr real _1_MINUS_G_OVER_B_pfx = 1.0 - G_OVER_B_pfx;
static constexpr real C_pfx = i_0 / (1 - i_0 - G_OVER_B_pfx);
static constexpr long real Ro_pfx = BETA_pfx / GAMMA;
static constexpr long real Ro_pfx_approx = BETA_pfx_approx / GAMMA;
static constexpr long real i_inf_pfx = 1 - (GAMMA / BETA_pfx);
static constexpr long real i_inf_pfx_approx = 1 - (GAMMA / BETA_pfx_approx);

#ifdef INFECTED_FRACTION
static constexpr real OVERLOOK_RATE = 1;									// ----> Depending on the initial settings, the simulation may generate a firehose of data, which in turn becomes highly inconvenient (or even prohibitive) for being written to a file (in terms of either space and time). For such cases, it is strongly adviseable to purposely overlook a portion of the events generated per time unit.
static constexpr uint OVERLOOK = (uint)((NUM_AGENTS * OVERLOOK_RATE)/1);	// ----> Depending on the initial settings, the simulation may generate a firehose of data, which in turn becomes highly inconvenient (or even prohibitive) for being written to a file (in terms of either space and time). For such cases, it is strongly adviseable to purposely overlook a portion of the events generated per time unit.
#endif

#ifdef CLIQUE
static const uint LIST_INI_SZ = (uint)(round(std::max((real)2.0, (real)NUM_AGENTS / (3 * N)))); // ----> Initial size of both 'sAgents' and 'iAgents' lists. Every time a node's list become full, its size gets doubled. Although arbitrary, the initial value provided here aims at reducing both the number of times a doubling operation is required and the vector's final size.
#else
const uint LIST_INI_SZ = (uint)(round(std::max((real)2.0, (real)NUM_AGENTS / (3 * N)))); // ----> Initial size of both 'sAgents' and 'iAgents' lists. Every time a node's list become full, its size gets doubled. Although arbitrary, the initial value provided here aims at reducing both the number of times a doubling operation is required and the vector's final size.
#endif

// * SIMULATION UTILS *
#ifdef i_t_FROM_MODEL
real i_t(const real& t);
real i_t_2ndMmt(const real& t);
real i_t_2ndMmt_naive(const real& t);
real i_t_2ndMmt_logistic(const real& t);
//real i_t_2ndMmt_sys(const real& t);
real i_t_pfx(const real& t);
#endif

#ifdef SOLVE_NUMERICALLY
static constexpr long real crowdFactor = std::min(2.0, std::max((real)K / N, 1.0));
real didt(const real& i);
void rungeKutta4thOrder(const real& t0, const real& i0, const real& t, const real& h, const real& epsilon, std::vector<real>& saveToFile, uint& outputSize, const uint& outputGranularity = 50, const real& largerDetailUntil = 1000);
#endif

} //namespace sim

// ---------------------------//----------------------------- //