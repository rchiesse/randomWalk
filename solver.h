#pragma once
#include "defs.h"
namespace sim {
	class Solver
	{
	private:
		static rtl nT, nL, nG;
		static rtl Wi, Ws, w;
		static rtl sigma;
		static uint numAgents;
		static std::vector<rtl> block_walkRate;			
		static rtl master_esigma;
		static rtl master_walkRate;

	public:
		static void setParams(const rtl& tau, const rtl& lambda, const rtl& gamma, const uint& NUM_AGENTS, const rtl& _Wi = 1.0, const rtl& _Ws = 1.0);
		static void setBlockData();
		//static rtl diadt(const rtl& ia, const rtl& sumSbIb);
		//static rtl dildt(const rtl& ia, const rtl& il);

#ifdef PER_BLOCK
		static rtl diabdt(const rtl& Ia, const rtl& Iab, const rtl& Sab, const uint& block);
		static rtl dsabdt(const rtl& Ia, const rtl& Iab, const rtl& Sab, const uint& block);
#else //PER_BLOCK
		static rtl divbdt(const rtl& Ia, const rtl& Iv, const rtl& Sv, const uint& block);
		static rtl dsvbdt(const rtl& Ia, const rtl& Iv, const rtl& Sv, const uint& block);
#endif //PER_BLOCK
		static void step(const rtl& h, rtl& Ia, std::vector<rtl>& v_Iab, std::vector<rtl>& v_Sab);
		//static rtl dilbdt(const rtl& ia, const rtl& il, const rtl& iab, const rtl& ilb, const uint& block);
		static void update_Ia(rtl& Ia, const std::vector<rtl>& v_Iab);
		static void update_Ia(rtl& Ia, const std::vector<rtl>& v_Iab, const std::vector<rtl>& base, const double& fraction = 1.0);
		static void lookAhead(const rtl& h, rtl& Ia, const std::vector<rtl>& v_Iab, const std::vector<rtl>& v_Sab, std::vector<rtl>& target);
		static void lookAhead(const rtl& h, rtl& Ia, const std::vector<rtl>& v_Iab, const std::vector<rtl>& v_Sab, std::vector<rtl>& target, std::vector<rtl>& base, const double& fraction = 1.0);
		static void rungeKutta4thOrder(const rtl& t0, std::vector<rtl>& v_Iab, std::vector<rtl>& v_Sab, std::vector<rtl>& v_ilb, const rtl& t, const rtl& h, const rtl& epsilon, std::vector<rtl>& saveToFile_diadt, std::vector<rtl>& saveToFile_dildt, uint& outputSize, const uint& outputGranularity = 50, const rtl& largerDetailUntil = 1000);

		static rtl dIdt(const rtl& Ia);
		static void rkMaster(const rtl& t0, std::vector<rtl>& v_Iab, std::vector<rtl>& v_Sab, const rtl& t, const rtl& h, const rtl& epsilon, std::vector<rtl>& saveToFile_diadt, std::vector<rtl>& saveToFile_dildt, uint& outputSize, const uint& outputGranularity = 50, const rtl& largerDetailUntil = 1000);
		static void stepMaster(const rtl& h, rtl& Ia);
		//static void lookAheadMaster(const rtl& h, rtl& Ia, rtl& target);
		//static void lookAheadMaster(const rtl& h, rtl& Ia, rtl& target, rtl& base, const double& fraction = 1.0);
	};
} //namespace sim

