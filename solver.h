#pragma once
#include "defs.h"
namespace sim {
	class Solver
	{
	private:
		static real nT, nL, nG;
		static real sigma;
		static real EULER;													// ----> The Euler–Mascheroni constant.
		static uint numAgents;
	public:
		static void setParams(const real& tau, const real& lambda, const real& gamma, const uint& NUM_AGENTS);

		static real diadt(const real& ia, const real& sumSbIb);
		static real dildt(const real& ia, const real& il);
#ifdef CLIQUE
#ifdef PER_BLOCK
		static real diabdt(const real& Ia, const real& Sa);
		static real dsabdt(const real& Ia, const real& Sa);
#else //PER_BLOCK
		static real divbdt(const real& Ia, const real& Iv, const real& Sv);
		static real dsvbdt(const real& Ia, const real& Iv, const real& Sv);
#endif //PER_BLOCK
		static void step(const real& h, real& Ia, real& Sa);
		//static real dilbdt(const real& ia, const real& il, const real& iab, const real& ilb, const uint& block);
		static void lookAhead(const real& h, real& Ia, real& Sa, std::vector<real>& target);
		static void lookAhead(const real& h, real& Ia, real& Sa, std::vector<real>& target, std::vector<real>& base, const real& fraction = 1.0);
		static void rungeKutta4thOrder(const real& t0, real& Ia, real& Sa, const real& t, const real& h, const real& epsilon, std::vector<real>& saveToFile_diadt, std::vector<real>& saveToFile_dildt, uint& outputSize, const uint& outputGranularity = 50, const real& largerDetailUntil = 1000);

#else //CLIQUE
#ifdef PER_BLOCK
		static real diabdt(const real& Ia, const real& Iab, const real& Sab, const uint& block);
		static real dsabdt(const real& Ia, const real& Iab, const real& Sab, const uint& block);
#else //PER_BLOCK
		static real divbdt(const real& Ia, const real& Iv, const real& Sv, const uint& block);
		static real dsvbdt(const real& Ia, const real& Iv, const real& Sv, const uint& block);
#endif //PER_BLOCK
		static void step(const real& h, real& Ia, std::vector<real>& v_Iab, std::vector<real>& v_Sab);
		static real dilbdt(const real& ia, const real& il, const real& iab, const real& ilb, const uint& block);
		static void update_Ia(real& Ia, const std::vector<real>& v_Iab);
		static void update_Ia(real& Ia, const std::vector<real>& v_Iab, const std::vector<real>& base, const double& fraction = 1.0);
		static void lookAhead(const real& h, real& Ia, const std::vector<real>& v_Iab, const std::vector<real>& v_Sab, std::vector<real>& target);
		static void lookAhead(const real& h, real& Ia, const std::vector<real>& v_Iab, const std::vector<real>& v_Sab, std::vector<real>& target, std::vector<real>& base, const double& fraction = 1.0);
		static void rungeKutta4thOrder(const real& t0, std::vector<real>& v_Iab, std::vector<real>& v_Sab, std::vector<real>& v_ilb, const real& t, const real& h, const real& epsilon, std::vector<real>& saveToFile_diadt, std::vector<real>& saveToFile_dildt, uint& outputSize, const uint& outputGranularity = 50, const real& largerDetailUntil = 1000);
#endif //CLIQUE

		static real dIdt(const real& Ia);
		static void rkMaster(const real& t0, std::vector<real>& v_Iab, std::vector<real>& v_Sab, const real& t, const real& h, const real& epsilon, std::vector<real>& saveToFile_diadt, std::vector<real>& saveToFile_dildt, uint& outputSize, const uint& outputGranularity = 50, const real& largerDetailUntil = 1000);
		static void stepMaster(const real& h, real& Ia);
		static void lookAheadMaster(const real& h, real& Ia, real& target);
		static void lookAheadMaster(const real& h, real& Ia, real& target, real& base, const double& fraction = 1.0);
	};
} //namespace sim

