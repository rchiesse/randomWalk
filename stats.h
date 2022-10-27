#pragma once
#include "defs.h"
#include "reporter.h"

namespace sim {
	enum class streamType { infFrac, occupancy, avDuration, probs };
	using std::vector;
	class Stats {
	private:
		static real t;
		static uint numAgents;
		static uint rounds;
		static real tau;
		static real gamma;
		static real lambda;
		static real Wi;
		static real Ws;
		static real W;
		static std::ofstream infFracData;
#ifdef INFECTED_FRACTION
	private:
		static const uint BUFFER_SIZE = (uint)10e5;
		static uint _bufferPos;

		//Structures to measure the fraction of infected agents per time:
		static std::vector<real> infectedAgFractionBuffer;
		static std::vector<real> infectedSiteFractionBuffer;
		static std::vector<real> timestampBuffer;
		static std::vector<uint> agentBuffer;
		//static std::vector<std::string> actionBuffer;
	public:
		//File manipulation:
		static void iFracToFile(const uint& overlook);
		static void bufferizeIFrac(const int& ag, const real& now, const std::string& action, const uint& itotal, const uint& iltotal, const uint& NUM_AGENTS, const uint& EVT_GRANULARITY);
#endif

#ifdef OCCUPANCY
		//Structures to measure the average number of agents per node over all rounds:
	private:
		static std::ofstream occData;
		static std::ofstream netOccData;
		static void updateSusOccAv(const graph::node& v, const uint& prevNumA, const uint& prevNumS, const real& now);
		static void checkSusMaxOcc(const graph::node& v, const uint& numA, const uint& numS);
		static void checkSusMinOcc(const graph::node& v, const uint& numA, const uint& numS);
		static void updateInfOccAv(const graph::node& v, const uint& prevNumA, const uint& prevNumI, const real& now);
		static void checkInfMaxOcc(const graph::node& v, const uint& numA, const uint& numI);
		static void checkInfMinOcc(const graph::node& v, const uint& numA, const uint& numI);
	public:
		static std::vector<real> avOcc;											// ----> Average occupancy of each node by agents.                   There are 'n' accumulators in total (one for each node).
		static std::vector<real> avSOcc;										// ----> Average occupancy of each node by susceptible agents only.  There are 'n' accumulators in total (one for each node).
		static std::vector<real> avIOcc;										// ----> Average occupancy of each node by infected agents only.     There are 'n' accumulators in total (one for each node).
		static std::vector<real> lastUpdate;
		static std::vector<real> sLastUpdate;
		static std::vector<real> iLastUpdate;
		static std::vector<uint> maxOcc;										// ----> A "max-occupancy" counter for every node, to measure up to how many agents happen to be at the same node the same time.
		static std::vector<uint> maxSOcc;										// ----> A "max-occupancy" counter for every susceptible node, to measure up to how many agents happen to be at the same node the same time.
		static std::vector<uint> maxIOcc;										// ----> A "max-occupancy" counter for every infected node, to measure up to how many agents happen to be at the same node the same time.
		static std::vector<uint> minOcc;										// ----> A "min-occupancy" counter for every node, to measure down to how many agents a node happen to be left with, considering the entire simulation time T.
		static std::vector<uint> minSOcc;										// ----> A "min-occupancy" counter for every susceptible node, to measure down to how many agents a node happen to be left with, considering the entire simulation time T.
		static std::vector<uint> minIOcc;										// ----> A "min-occupancy" counter for every infected node, to measure down to how many agents a node happen to be left with, considering the entire simulation time T.

		static void computeOcc(const real& totalTime);
		//static void occToFile(const real& Ws, const uint& numAg);
		//static void netInfAgMeanToFile(const real& Ws, const uint& numAg);
		static void initOcc(const uint& ntwSize);
		static void getGlobalMaxOcc(uint& _max, graph::node& _where);			// ----> Determines the maximum number of agents (S, I and both) observed in the same node during a simulation run.
		static void getGlobalMinOcc(uint& _min, graph::node& _where);			// ----> Determines the minimum number of agents (S, I and both) observed in the same node during a simulation run.

		//Two-purpose functions. They depart from the premise that v's new number of infected (resp. susceptible) agents has increased (resp. decreased) by 1 unit. Then, the functions will (i) update the average of both "general" and "infected (resp. susceptible)" occupancy of v accordingly, and (ii) update (if necessary) v's maximum (resp. minimum) occupancy for both the general and the infected (resp. susceptible) types of agents.
		static void updateInfOccRaised(const graph::node& v, const uint& numA, const uint& numI, const real& now);
		static void updateInfOccLowered(const graph::node& v, const uint& numA, const uint& numI, const real& now);
		static void updateSusOccRaised(const graph::node& v, const uint& numA, const uint& numS, const real& now);
		static void updateSusOccLowered(const graph::node& v, const uint& numA, const uint& numS, const real& now);


		/*					***************  IMPORTANT  ****************
			The templates below play the important role of not letting the PUBLIC functions above
			be compiled whenever the type for ANY of their actual arguments do not explicitly match
			their formal counterparts. This reduces the risk of user errors.

			More specifically, no implicit cast will be allowed when calling such functions.
			For example, a 'uint_32' cannot be passed as an argument for "graph::node", nor can a 'double'
			be passed as 'uint', and so on (all this would be possible otherwise).

			Indeed, consider the following declarations:

				graph::node _v = 4;
				uint _A = 30
				uint _I = 5;
				real _time = 1.034;

			Now let's say we want to call 'updateInfOccRaised()' with these parameters.
			THE CORRECT CALL WOULD BE

				updateInfOccRaised(_v, _A, _I, _time);

			but no compile error, and possibly---and most dangerously---NO RUNTIME ERROR would
			be thrown if we were to inadvertently call

				updateInfOccRaised(_time, _A, _I, _v);

			and this could in turn lead to REALLY-HARD-TO-DEBUG ERRORS.

			**OBS: "= delete" since C++11.  */
		template <class arg_UINT, class arg_NODE>
		static void getGlobalMaxOcc(arg_UINT, arg_NODE) = delete;
		template <class arg_UINT, class arg_NODE>
		static void getGlobalMinOcc(arg_UINT, arg_NODE) = delete;
		template <class arg_NODE, class arg_UINT, class arg_DOUBLE>
		void updateInfOccRaised(arg_NODE, arg_UINT, arg_UINT, arg_DOUBLE) = delete;
		template <class arg_NODE, class arg_UINT, class arg_DOUBLE>
		void updateInfOccLowered(arg_NODE, arg_UINT, arg_UINT, arg_DOUBLE) = delete;
		template <class arg_NODE, class arg_UINT, class arg_DOUBLE>
		void updateSusOccRaised(arg_NODE, arg_UINT, arg_UINT, arg_DOUBLE) = delete;
		template <class arg_NODE, class arg_UINT, class arg_DOUBLE>
		void updateSusOccLowered(arg_NODE, arg_UINT, arg_UINT, arg_DOUBLE) = delete;

#endif //OCCUPANCY

		//Structures to measure the epidemic's average duration, over the number of ROUNDS provided:
	private:
		static real avDur;										// ----> Epidemics average duration;
		static real lastRoundDuration;							// ----> Stores the simulation time for only the most recent round;
		static uint partials;									// ----> Total number of values to average upon.
		static bool avDurComputed;
	public:
		static std::ofstream avDurDataGroupL;	    // Lists all the time series to be later included at the same plot.
		static std::ofstream avDurDataGroupK;	    // Lists all the time series to be later included at the same plot.
		static std::ofstream avDurDataK;	// For some fixed PrE, average duration of epidemics for different number of agents (K)
		static std::ofstream avDurDataL;	// For some fixed PrE, average duration of epidemics for different values of the walk rate (\lambda)

		//Average Duration (AD) partials. It simply assigns the informed 'duration' to the 'averageDuration' accumulator. The actual average is later computed by calling 'avDuration()'.
		static void initAvDur();
		static void partialsAvDur(const real& duration);
		static const real& avDuration();
		static void resetAvDur();

#ifdef ESTIMATE_PROBS
		static std::ofstream probsData;
		static real meetings;
		static real siMeetings;
		static real totalExposition;
		static real totalFate;
		static real totalInfections;

		static real Psi;										// ----> Probability of SI meetings (= siMeetings/meetings)
		static real Pinf;										// ----> Probability of an S agent become I (= Psi * (totalInfections/siMeetings))
		static real meetingRate_;								// ----> The rate at which agents meet, irrespective of their state (= meetings/T)
		static real siMeetingRate;								// ----> The rate of SI meetings (= siMeetings/T)
		static void probsToFile();								// ----> Computes 'Psi', 'Pinf', 'meetingRate' and 'siMeetingRate' based on the simulation-generated data.

#endif

//Utils
		static void setParams(const real& _t, const uint& _numAgents, const uint& _rounds, const real& _tau, const real& _gamma, const real& _lambda, const real& _Wi = 1.0, const real& _Ws = 1.0);
		static void initStream(const streamType& s);
		static void endStream(const streamType& s);
		static void writeToFile(const streamType& s, const real& Ws, const real& Wi, const uint& numAgents);
		static bool isEmpty(std::ofstream& s);
		static void resetStats();
		static void setBasename();
#ifdef INFECTED_FRACTION
		static void genPlotScript(const std::string& referenceFile, const bool&& numericOnly = false);
#endif
	};
} //namespace sim
