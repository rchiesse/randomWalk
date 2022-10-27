#pragma once

#include <stdio.h>
#include <tchar.h>
#include <iostream>
#include <iomanip>
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
namespace sim {
	constexpr double epsilon = 1.0 / 10e7;
}
#define assertm(exp, msg) assert(((void)msg, exp))
//#define DEBUG

// ---------------------------//----------------------------- //

// * DATA TYPES *
#define uint uint32_t
#define ulong uint64_t
#define ushort uint16_t
#define real long double
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
//#define CLIQUE
//#define READ_NTWK_FROM_FILE
#define GNP
//#define STAR

#ifdef STAR
static constexpr uint N = 200;
#define SOURCE_FILE "", 200
#define NWTK_LABEL "Star"
#define SHORT_LABEL "star"
#endif

#ifdef CLIQUE
static constexpr uint N = 10;										// ----> Network size
#define NWTK_LABEL "Clique"
#define SHORT_LABEL "CL"
#endif

#ifdef GNP
//static constexpr uint N = 12008;
static constexpr uint N = 10000;
#define SOURCE_FILE "", 10000
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
//static constexpr uint N = 120;										// ----> Network size
//#define SOURCE_FILE std::string(std::string(EXE_DIR) + std::string("/synthetic/BA-120-3.txt")), 120
//#define NWTK_LABEL "BA-120-3"
//#define SHORT_LABEL "BA-120-3" 

static constexpr uint N = 1200;										// ----> Network size
#define SOURCE_FILE std::string(std::string(EXE_DIR) + std::string("/synthetic/BA-1200-3.txt")), 1200
#define NWTK_LABEL "BA-1200-3"
#define SHORT_LABEL "BA-1200-3" 

//static constexpr uint N = 12008;										// ----> Network size
//#define SOURCE_FILE std::string(std::string(EXE_DIR) + std::string("/synthetic/BA.txt")), 12008
//#define NWTK_LABEL "BA-12k-10"
//#define SHORT_LABEL "BA-12k-10" 

//static constexpr uint N = 500;										// ----> Network size
//#define SOURCE_FILE std::string(std::string(EXE_DIR) + std::string("/synthetic/BA-500-3.txt")), 500
//#define NWTK_LABEL "BA-500-3"
//#define SHORT_LABEL "BA-500-3" 

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
#define AUTO_RELATION					// ----> Gives agents the option of staying at their current node upon their walk event. If not enabled, agents will necessarily change their current node when their walk event is processed.
#define PROTECTION_FX
#ifndef PROTECTION_FX
static constexpr real Ws = 1.0;	// !DO NOT CHANGE THIS LINE! To set Ws to 1.0 here means "no protection effect", which is the desired behaviour when the pre-processor macro "PROTECTION_FX" is not defined.
static constexpr real Wi = 1.0;	// !DO NOT CHANGE THIS LINE! To set Wi to 1.0 here means "no protection effect", which is the desired behaviour when the pre-processor macro "PROTECTION_FX" is not defined.
#endif

// ---------------------------//----------------------------- //


// * STATS * 
#define INFECTED_FRACTION
#define SI_PROPORTION
//#define OCCUPANCY
//#define ESTIMATE_PROBS

// ---------------------------//----------------------------- //


// * REPORTER * 
#define MEASURE_ROUND_EXE_TIME

// ---------------------------//----------------------------- //


// * NUMERICAL SOLUTION *
namespace sim{		// ----> Simulator's namespace.
#define SOLVE_NUMERICALLY
//#define BYPASS_SIMULATION
#define NORM_SITE_PER_AG
#define PER_BLOCK													// ----> If defined, then the numerical solution is based on 2 equations per degree-block. If otherwise, then a fine-grained system of 2 equations per NODE is solved (computationally expensive).
} //namespace sim

// ---------------------------//----------------------------- //