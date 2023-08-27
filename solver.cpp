#include "solver.h"
#include "graph.h"
using namespace sim;

rtl Solver::nT = 0;
rtl Solver::nL = 0;
rtl Solver::nG = 0;
rtl Solver::Wi = 0;
rtl Solver::Ws = 0;
rtl Solver::w = 0;
rtl Solver::sigma = 0;
rtl Solver::EULER = 0.57721566490153286060651209008240243104215933593992;		// ----> The Euler–Mascheroni constant.
uint Solver::numAgents = 0;

void Solver::setParams(const rtl& tau, const rtl& lambda, const rtl& gamma, const uint& NUM_AGENTS, const rtl& _Wi, const rtl& _Ws) {
	nT = tau;
	nL = lambda;
	nG = gamma;
#ifdef PROTECTION_FX
	Wi = _Wi;
	Ws = _Ws;
#else
	Wi = 1.0;
	Ws = 1.0;
#endif
	w = (Wi + Ws) / 2.0;
	numAgents = NUM_AGENTS;
	sigma = nT / (2.0 * nL + nT); 
}

#ifdef CLIQUE
rtl Solver::diabdt(const rtl& Ia, const rtl& Sa) {
	const rtl pb = 1.0;
	const rtl nb = graph::Graph::n;
	const rtl qb = 1.0;
	const rtl ibnb = Ia / nb;
	const rtl sbnb = Sa / nb;
	const rtl kbnb = ibnb + sbnb;
	const rtl b = N;
	const rtl K = Ia + Sa;
	const rtl p = b / (2.0 * graph::Graph::m - N);
	const rtl E_X = K * p;
	const rtl E_X2 = K*K * p*p;
	const rtl Var_X = E_X * (1.0 - p);
	const rtl stdDev = sqrt(Var_X);
	//const rtl cs = (1.0 - (sbnb / kbnb) * (ibnb / kbnb)) * stdDev;
	//const rtl Sa = (double)NUM_AGENTS - Ia;

	//RONALD (BEST SO FAR):
	if (Sa == 0.0)
		return -(nG * Ia);

	//RONALD v2 (ótimo em regime esparso):
	//return (Ia - Iab) * LAMBDA * qb - Iab * LAMBDA * (1.0 - qb)
	//	+ 2 * ((Sab * Iab)/nb) * LAMBDA * SIGMA_aa
	//	- (GAMMA_a * Iab);
	rtl sigma = nT / (2.0 * nL + nT);
	return
		2.0 * ((Sa * Ia) / nb) * nL * sigma
		- (nG * Ia);

	//long double H = EULER * log(sbnb + 1.0) / (2.0 * nL);
	//long double H = EULER * log(sbnb * ibnb * (1.0 / nT)) / (2.0 * nL);
	//long double dFactor = (sbnb * ibnb / nT) * std::min(1.0, abs(nT - nL));
	//long double H = EULER * log(ibnb + 1.0 + dFactor) / (2.0 * nL);
	//long double H = EULER * log((ibnb * (std::min(sbnb, ibnb) / (2.0))) + 1.0) / (2.0 * nL);
	//long double H = EULER * log((sbnb)+1.0) / (2.0 * nL) + EULER * log((ibnb)+1.0) / (2.0 * nL);
	//long double ii = ibnb + 1.0;
	//const double prob_inf = ii / (H + ii);
	//const double p = std::min(sbnb, 1.0) * std::min(ibnb, 1.0);
	//rtl Hs = std::max(1.0l, EULER * log((sbnb)+1.0) / (2.0 * nL));
	//rtl Hi = std::max(1.0l, EULER * log((ibnb)+1.0) / (2.0 * nL));
	//const rtl infRate = (sbnb > ibnb) ? (nT / Hs) : nT;
	//const rtl sigma = (ibnb * nT) / (nL + (ibnb * nT));
	//const rtl ve = Var_X / E_X;
	//const rtl ev = E_X / Var_X;
	//const rtl sd = sqrt(Var_X); // ----> Standard devation.
	//const rtl sdi = sd * (ibnb / kbnb);
	//const rtl sds = sd * (sbnb / kbnb);
	//const rtl dimm = (ve * ibnb * nG) / ((ev * sbnb) * nT + ve * ibnb * nG);
	//const rtl dimm = (ibnb * nG) / ((sbnb) * nT + ibnb * nG);
	//rtl ii = sbnb + 1.0;
	//rtl H = std::max(1.0l, EULER * log((sbnb)+1.0) / (2.0 * nL));
	//const rtl prob_inf = ii / (H + ii);
	
	//const rtl prob_inf = ii / (H + ii);
	//rtl H = (2.0 * nL) / EULER * log(sbnb);
	//rtl H = EULER * log(sbnb) / (2.0 * nL);
	//rtl H = std::max(1.0l, EULER * log(sbnb+1.0) / (nL));
	//rtl H = EULER * log(sbnb);
	//rtl ii = std::max(ibnb, 1.0l);
	//rtl H = EULER * log(sbnb+1.0) / (nL);
	//rtl H = std::max(1.0l, EULER * log(kbnb) / (2.0 * nL));	//LEGAL TB!
	//const rtl prob_inf = (ii * nT) / (H + (ii * nT));	//RON (figura plot) (aqui a presença de nT e nL n faz diferença. Avaliar outros cenários!)
	//rtl H = EULER * log(sbnb * ibnb) / (2.0 * nL); //INTERESSANTE...
	//rtl H = EULER * log(std::min(sbnb, ibnb)) / (2.0 * nL); // MT INTERESSANTE TB...
	
	//REGIME DENSO:
	////rtl H = EULER * log(kbnb) / (2.0 * nL);	//RON
	////rtl H = EULER * ((log(sbnb) / (2.0 * nL)) + (log(ibnb) / (2.0 * nL)));
	//rtl H = EULER * (log(sbnb * ibnb) / (2.0 * nL));
	//rtl ii = ibnb * nT;
	////rtl ii = ibnb;
	//const rtl prob_inf = (ii ) / (H + (ii ));
	////const rtl prob_inf = 1.0;
	//return
	//	//E-M:
	//	//+ Ia * sbnb * prob_inf * nT
	//	//- (nG * Ia);
	//
	//	//POR NÓ:
	//	//nb * (
	//	//	(Ia - ibnb) * nL / nb
	//	//	- ibnb * nL * (b - 1) / b
	//	//	+ ibnb * sbnb * nT
	//	//	- (nG * ibnb)
	//	//);
	//
	//	//EQUIVALENTE MASTER:
	//	((nT * Sa * Ia ) / (N)) - nG * Ia;
}

rtl Solver::dsabdt(const rtl& Ia, const rtl& Sa) {
	const rtl pb = 1.0;
	const rtl nb = graph::Graph::n;
	const rtl qb = 1.0;
	const rtl ibnb = Ia / nb;
	const rtl sbnb = Sa / nb;
	const rtl kbnb = ibnb + sbnb;
	const rtl b = N;
	const rtl pIn  = (b - 1) / (2.0 * graph::Graph::m - 1.0);		// ----> Probability that a randomly chosen link leads an agent to an specific node v_b, provided that the agent is outside.
	const rtl pOut = (b - 1) / b;									// ----> Probability of leaving a node v_b, provided that the agent is inside.
	const rtl K = Ia + Sa;
	const rtl p = b / (2.0 * graph::Graph::m - N);
	const rtl E_X = K * p;
	const rtl E_X2 = K * K * p * p;
	const rtl Var_X = E_X * (1.0 - p);
	const rtl stdDev = sqrt(Var_X);
	//const rtl cs = (1.0 - (sbnb / kbnb) * (ibnb / kbnb)) * stdDev;
	//const rtl Sa = (double)NUM_AGENTS - Ia;
	//BEST SO FAR:
	if (Sa == 0.0) 
		return nG * Ia;

	//rtl H = EULER * log(sbnb + ibnb + 1.0) / (2.0 * nL);
	//rtl H = EULER * log(sbnb * ibnb * (1.0 / nT)) / (2.0 * nL);
	//rtl dFactor = (sbnb * ibnb / nT) * std::min(1.0, abs(nT - nL));
	//rtl H = EULER * log(ibnb + 1.0 + dFactor) / (2.0 * nL);
	//rtl H = EULER * log((ibnb * (std::min(sbnb, ibnb) / (2.0))) + 1.0) / (2.0 * nL);
	//rtl H = EULER * log((sbnb)+1.0) / (2.0 * nL) + EULER * log((ibnb)+1.0) / (2.0 * nL);
	//rtl Hi = std::max(1.0l, EULER * log((ibnb)+1.0) / (2.0*nL));
	//const rtl infRate = (sbnb > ibnb) ? (nT - Hs) : nT;
	//const rtl sigma = (ibnb * nT) / (nL + (ibnb * nT));
	//const rtl ve = Var_X / E_X;
	//const rtl ev = E_X / Var_X;
	//const rtl dimm = (ibnb * nG) / ((sbnb) * nT + ibnb * nG);
	//rtl ii = ibnb + 1.0;
	//rtl H = std::max(1.0l, EULER * log((ibnb)+1.0) / (2.0 * nL));
	//const rtl prob_inf = ii / (H + ii);

	//rtl H = (2.0 * nL) / EULER * log(sbnb);
	//rtl H = std::max(1.0l, EULER * log(sbnb+1.0) / (nL));
	//rtl H = EULER * log(sbnb);
	//rtl H = EULER * log(sbnb) / (2.0 * nL);
	//rtl ii = std::max(ibnb, 1.0l);
	//rtl H = EULER * log(sbnb+1.0) / (nL);	//Ótimo qnd LAMBDA >> TAU; SUPERESTIMA C.C (e H = 1.0 acaba sendo melhor)
	//rtl H = EULER * log(kbnb) / (2.0 * nL);	//RON
	//rtl H = std::max(1.0l, EULER * log(kbnb) / (2.0 * nL));	//LEGAL TB!
	//rtl H = EULER * log(ibnb) / (2.0 * nL);
	//const rtl prob_inf = (ii * nT) / (H + (ii * nT));
	//rtl H = EULER * log(sbnb * ibnb) / (2.0 * nL);	// INTERESSANTE...
	//rtl H = EULER * log(std::min(sbnb, ibnb)) / (2.0 * nL); // MT INTERESSANTE TB...
	
	//RON
	//rtl H = EULER * log(kbnb) / (2.0 * nL);
	//rtl ii = ibnb;
	
	//RONALD v2 (ótimo em regime esparso):
	//return (Ia - Iab) * LAMBDA * qb - Iab * LAMBDA * (1.0 - qb)
	//	+ 2 * ((Sab * Iab)/nb) * LAMBDA * SIGMA_aa
	//	- (GAMMA_a * Iab);
	rtl sigma = nT / (2.0 * nL + nT);
	return 
		- 2.0 * ((Sa * Ia)/nb) * nL * sigma
		+ (nG * Ia);

	////ALTA DENSIDADE:
	////Muito preciso em alta densidade:
	//rtl H = EULER * (log(sbnb * ibnb) / (2.0 * nL));
	//rtl ii = ibnb * nT;
	//
	//const rtl prob_inf = (ii ) / (H + (ii ));
	////const rtl prob_inf = (ii) / (1.0 + (ii));	//Ótimo qnd TAU == LAMBDA; SUPERESTIMA qnd LAMBDA >> TAU (e EULER acaba sendo melhor)
	////const rtl prob_inf = 1.0;
	//return
	//	//POR NÓ:
	//	//-Ia * sbnb * prob_inf * nT
	//
	//	//E-M:
	//	//-Ia * sbnb * prob_inf * nT
	//	//+ (nG * Ia);
	//
	//	//OFFICIAL:
	//	//nb * (
	//	//	(Sa - sbnb) * nL / nb
	//	//	- sbnb * nL * (b - 1) / b
	//	//	- ibnb * sbnb * nT
	//	//	+ (nG * ibnb)
	//	//);
	//
	//	//EQUIVALENTE MASTER:
	//	-((nT * Sa * Ia ) / (N)) + nG * Ia;
}

void Solver::step(const rtl& h, rtl& Ia, rtl& Sa) {
	constexpr rtl one_sixth = 1.0 / 6.0;
	std::vector<rtl> k1(2, 0), k2(2, 0), k3(2, 0), k4(2, 0);

	lookAhead(h, Ia, Sa, k1);
	lookAhead(h, Ia, Sa, k2, k1, 0.5);
	lookAhead(h, Ia, Sa, k3, k2, 0.5);
	lookAhead(h, Ia, Sa, k4, k3);

	//Take step:
	Ia += one_sixth * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
	Sa += one_sixth * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
}

void Solver::lookAhead(const rtl& h, rtl& Ia, rtl& Sa, std::vector<rtl>& target) {
	target[0] = h * diabdt(Ia, Sa);
	target[1] = h * dsabdt(Ia, Sa);
}

void Solver::lookAhead(const rtl& h, rtl& Ia, rtl& Sa, std::vector<rtl>& target, std::vector<rtl>& base, const rtl& fraction) {
	target[0] = h * diabdt(Ia + fraction * base[0], Sa + fraction * base[1]);
	target[1] = h * dsabdt(Ia + fraction * base[0], Sa + fraction * base[1]);
}

void Solver::rungeKutta4thOrder(const rtl& t0, rtl& Ia, rtl& Sa, const rtl& t, const rtl& h, const rtl& epsilon, std::vector<rtl>& saveToFile_diadt, std::vector<rtl>& saveToFile_dildt, uint& outputSize, const uint& outputGranularity, const rtl& largerDetailUntil) {
	uint totalSteps = (uint)((t - t0) / h) + 1;
	saveToFile_diadt.resize((uint64_t)largerDetailUntil + (totalSteps - ((uint)largerDetailUntil) / outputGranularity) + 1);

	saveToFile_diadt[0] = Ia / numAgents;
	bool end = false;
	++outputSize;

	rtl sumH = 0;
	//For the first 'largerDetailUntil' iterations every step is stored in a vector ('saveToFile'), for later being written to file.
	for (uint s = 1; s < largerDetailUntil; ++s) {
		rtl prevIA = Ia;
		step(h, Ia, Sa);
		sumH += h;
		if (Ia < epsilon) {
			saveToFile_diadt[outputSize] = 0;
			end = true;
			++outputSize;
			break;
		}
		if (Ia == prevIA) {
			while (s < largerDetailUntil) {
				saveToFile_diadt[outputSize] = Ia / numAgents;
				++s;
				++outputSize;
			}
			end = true;
			break;
		}
		saveToFile_diadt[outputSize] = Ia / numAgents;
		++outputSize;
	}
	if (end) return;

	//From the 'largerDetailUntil' iteration on, we afford to ignore 'outputGranularity'-size windows of values, so that the saved file does not grow explosively.
	for (uint s = (uint)largerDetailUntil; s < totalSteps; ++s) {
		rtl prevIA = Ia;
		step(h, Ia, Sa);
		sumH += h;
		if (Ia < epsilon) {
			saveToFile_diadt[outputSize] = 0;
			++outputSize;
			break;
		}
		if (Ia == prevIA) {
			while (s < totalSteps) {
				if (s % outputGranularity == 0) {
					saveToFile_diadt[outputSize] = Ia / numAgents;
					++outputSize;
				}
				++s;
			}
			break;
		}
		if (s % outputGranularity == 0) {
			saveToFile_diadt[outputSize] = Ia / numAgents;
			++outputSize;
		}
	}
}

#else //CLIQUE
#ifdef PER_BLOCK
rtl Solver::diabdt(const rtl& Ia, const rtl& Iab, const rtl& Sab, const uint& b) {
	const rtl& pb = graph::Graph::block_prob[b];
	const rtl nb = graph::Graph::n * pb;
	const rtl& qb = graph::Graph::q_b[b];
	const rtl ibnb = Iab / nb;
	const rtl sbnb = Sab / nb;
	const rtl kbnb = ibnb + sbnb;
	const rtl Sa = (rtl)numAgents - Ia;
	const rtl& _k_b = graph::Graph::kb[b];
	const rtl _kvb_ = _k_b / nb;
	const rtl ii = (ibnb < 1.0) ? 1.0 : ibnb;
	const rtl ss = (sbnb < 1.0) ? 1.0 : sbnb;
	

	
	//RONALD PER BLOCK (BEST SO FAR):
	return nL * (Ia * qb - Iab)
		+ 2.0 * ((Sab * Iab) / nb) * nL * sigma * w - (nG * Iab);






	
	////rtl H = EULER * (log(sbnb * ibnb) / (2.0 * nL));		// ----> Ótimo em regime denso, mas falha se (sbnb * ibnb) < 1.0.
	//rtl H = EULER * (log((sbnb * ibnb) + 1.0) / (2.0 * nL));
	//rtl ii = ibnb;
	//const rtl prob_inf = ii / (H + ii);
	//return //Sa * nL * qb - Sab * nL
	//	Iab * sbnb * prob_inf * nT
	//	- (nG * Iab);

	////DENSO:
	////RONALD (BEST SO FAR):
	////if (Sab == 0.0) {
	////	return Ia * nL * qb - Iab * nL
	////		- (nG * Iab);
	////} 
	////rtl H = EULER * log(sbnb) / (2.0 * nL);
	////rtl H = EULER * log(kbnb) / (2.0 * nL);	//RON
	//rtl H = EULER * (log(sbnb * ibnb) / (2.0 * nL));
	////rtl H = EULER * log(sbnb * ibnb) / nL;	
	//rtl ii = ibnb;
	//const rtl prob_inf = (ii) / (H + ii);
	//return //Ia * nL * qb - Iab * nL
	//	+ Iab * sbnb * prob_inf * nT
	//	- (nG * Iab);
	//
	////NOVO TESTE (denso):
	////VERSÃO DE BLOCO EQUIVALENTE AO DE NÓ:
	////return nb * (
	////	(Ia - ibnb) * nL * (((rtl)b - 1.0) / (2.0 * graph::Graph::m - N - (b - 1))) - ibnb * nL * l
	////	+ ibnb * sbnb * nT
	////	- (nG * ibnb)
	////);
	////VERSÃO DE NÓ EQUIVALENTE AO DE BLOCO:
	////return Ia * nL * qb - Iab * nL
	////	+ Iab * sbnb * nT
	////	- (nG * Iab);
	//
	////ISOLANDO CADA BLOCO:
	////return nb * (
	////	+ ibnb * sbnb * nT
	////	- (nG * ibnb)
	////);
}

rtl Solver::dsabdt(const rtl& Ia, const rtl& Iab, const rtl& Sab, const uint& b) {
	const rtl& pb = graph::Graph::block_prob[b];
	const rtl nb = graph::Graph::n * pb;
	const rtl& qb = graph::Graph::q_b[b];
	const rtl ibnb = Iab / nb;
	const rtl sbnb = Sab / nb;
	const rtl kbnb = ibnb + sbnb;
	const rtl Sa = (rtl)numAgents - Ia;
	const rtl& _k_b = graph::Graph::kb[b];
	const rtl _kvb_ = _k_b / nb;
	const rtl ii = (ibnb < 1.0) ? 1.0 : ibnb;
	const rtl ss = (sbnb < 1.0) ? 1.0 : sbnb;
	
	//RONALD PER BLOCK (BEST SO FAR):
	return nL * (Sa * qb - Sab)
		- 2.0 * ((Sab * Iab) / nb) * nL * sigma * w + (nG * Iab);



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
	//if (Sab == 0.0) {
	//	return Sa * nL * qb - Sab * nL
	//		+ (nG * Iab);
	//}
	//rtl H = EULER * log(sbnb) / (2.0 * nL);
	//rtl H = EULER * log(kbnb) / (2.0 * nL);	//RON
	
	//if (Sab == 0.0 || Iab == 0.0)
	//	return (nG * Iab);
	//
	//rtl H = EULER * (log((sbnb * ibnb) + 1.0) / (2.0 * nL));
	//rtl ii = ibnb;
	//const rtl prob_inf = ii / (H + ii);
	//return //Sa * nL * qb - Sab * nL
	//	- Iab * sbnb * prob_inf * nT
	//	+ (nG * Iab);

	//DON : Good for dense scenarios or large walk rate
	//if (kb == 0.0)
	//	return Sa * LAMBDA * qb - Sab * LAMBDA;
	//return 
	//	Sa * LAMBDA * qb - Sab * LAMBDA
	//	- Iab * sbnb * TAU_aa
	//	+ (GAMMA_a * Iab);

	//NOVO TESTE (denso):
	//VERSÃO DE BLOCO EQUIVALENTE AO DE NÓ:
	//return Sa * nL * qb - Sab * nL
	//	- Iab * sbnb * nT
	//	+ (nG * Iab);
	//VERSÃO DE NÓ EQUIVALENTE AO DE BLOCO:
	//return nb * (
	//	(Sa - sbnb) * nL * (((rtl)b - 1.0) / (2.0 * graph::Graph::m - N - (b-1))) - sbnb * nL * l
	//	- ibnb * sbnb * nT
	//	+ (nG * ibnb)
	//);

	//ISOLANDO CADA BLOCO:
	//return nb * (
	//	- ibnb * sbnb * nT
	//	+ (nG * ibnb)
	//);
}

void Solver::step(const rtl& h, rtl& Ia, std::vector<rtl>& v_Iab, std::vector<rtl>& v_Sab) {
	constexpr rtl one_sixth = 1.0 / 6.0;
	const uint blocks = static_cast<uint>(graph::Graph::block_prob.size());
	std::vector<rtl> k1(2 * blocks, 0), k2(2 * blocks, 0), k3(2 * blocks, 0), k4(2 * blocks, 0);

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

void Solver::lookAhead(const rtl& h, rtl& Ia, const std::vector<rtl>& v_Iab, const std::vector<rtl>& v_Sab, std::vector<rtl>& target) {
	update_Ia(Ia, v_Iab);
	for (uint b = (uint)graph::Graph::block_prob.size() - 1; b > 0; --b) {
		if (graph::Graph::block_prob[b] == 0)
			continue;

		target[2 * b]		= h * diabdt(Ia, v_Iab[b], v_Sab[b], b);
		target[2 * b + 1]	= h * dsabdt(Ia, v_Iab[b], v_Sab[b], b);
	}
}

void Solver::lookAhead(const rtl& h, rtl& Ia, const std::vector<rtl>& v_Iab, const std::vector<rtl>& v_Sab, std::vector<rtl>& target, std::vector<rtl>& base, const double& fraction) {
	update_Ia(Ia, v_Iab, base, fraction);
	for (uint b = (uint)graph::Graph::block_prob.size() - 1; b > 0; --b) {
		if (graph::Graph::block_prob[b] == 0)
			continue;

		target[2 * b]		= h * diabdt(Ia, v_Iab[b] + fraction * base[2 * b], v_Sab[b] + fraction * base[2 * b + 1], b);
		target[2 * b + 1]	= h * dsabdt(Ia, v_Iab[b] + fraction * base[2 * b], v_Sab[b] + fraction * base[2 * b + 1], b);
	}
}

void Solver::update_Ia(rtl& Ia, const std::vector<rtl>& v_Iab) {
	Ia = 0;
	for (uint b = (uint)graph::Graph::block_prob.size() - 1; b > 0; --b)
		Ia += v_Iab[b];
}

void Solver::update_Ia(rtl& Ia, const std::vector<rtl>& v_Iab, const std::vector<rtl>& base, const double& fraction) {
	Ia = 0;
	for (uint b = (uint)graph::Graph::block_prob.size() - 1; b > 0; --b)
		Ia += v_Iab[b] + (fraction * base[2 * b]);
}
#else //PER_BLOCK
rtl Solver::divbdt(const rtl& Ia, const rtl& Iv, const rtl& Sv, const uint& b) {
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
rtl Solver::dsvbdt(const rtl& Ia, const rtl& Iv, const rtl& Sv, const uint& b) {
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
void Solver::step(const rtl& h, rtl& Ia, std::vector<rtl>& v_Iv, std::vector<rtl>& v_Sv) {
	constexpr rtl one_sixth = 1.0 / 6.0;
	const uint blocks = static_cast<uint>(graph::Graph::block_prob.size());
	vector<rtl> k1(2 * graph::Graph::n, 0), k2(2 * graph::Graph::n, 0), k3(2 * graph::Graph::n, 0), k4(2 * graph::Graph::n, 0);

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
void Solver::lookAhead(const rtl& h, rtl& Ia, const std::vector<rtl>& v_Iv, const std::vector<rtl>& v_Sv, std::vector<rtl>& target) {
	update_Ia(Ia, v_Iv);
	for (int v = (uint)graph::Graph::n - 1; v >= 0; --v) {
		target[2 * v] = h * divbdt(Ia, v_Iv[v], v_Sv[v], (uint)graph::Graph::g[v].size());
		target[2 * v + 1] = h * dsvbdt(Ia, v_Iv[v], v_Sv[v], (uint)graph::Graph::g[v].size());
	}
}
void Solver::lookAhead(const rtl& h, rtl& Ia, const std::vector<rtl>& v_Iv, const std::vector<rtl>& v_Sv, std::vector<rtl>& target, std::vector<rtl>& base, const double& fraction) {
	update_Ia(Ia, v_Iv, base, fraction);
	for (int v = (uint)graph::Graph::n - 1; v >= 0; --v) {
		target[2 * v] = h * divbdt(Ia, v_Iv[v] + fraction * base[2 * v], v_Sv[v] + fraction * base[2 * v + 1], (uint)graph::Graph::g[v].size());
		target[2 * v + 1] = h * dsvbdt(Ia, v_Iv[v] + fraction * base[2 * v], v_Sv[v] + fraction * base[2 * v + 1], (uint)graph::Graph::g[v].size());
	}
}
void Solver::update_Ia(rtl& Ia, const std::vector<rtl>& v_Iv) {
	Ia = 0;
	for (int v = (uint)graph::Graph::n - 1; v >= 0; --v)
		Ia += v_Iv[v];
}
void Solver::update_Ia(rtl& Ia, const std::vector<rtl>& v_Iv, const std::vector<rtl>& base, const double& fraction) {
	Ia = 0;
	for (int v = (uint)graph::Graph::n - 1; v >= 0; --v)
		Ia += v_Iv[v] + (fraction * base[2 * v]);
}
#endif //PER_BLOCK

rtl Solver::dilbdt(const rtl& Ia, const rtl& il, const rtl& Iab, const rtl& ilb, const uint& b) {
	return 0;
}

void Solver::rungeKutta4thOrder(const rtl& t0, std::vector<rtl>& v_Iab, std::vector<rtl>& v_Sab, std::vector<rtl>& v_ilb, const rtl& t, const rtl& h, const rtl& epsilon, std::vector<rtl>& saveToFile_diadt, std::vector<rtl>& saveToFile_dildt, uint& outputSize, const uint& outputGranularity, const rtl& largerDetailUntil) {
	uint totalSteps = (uint)((t - t0) / h) + 1;
	saveToFile_diadt.resize((uint64_t)largerDetailUntil + (totalSteps - ((uint)largerDetailUntil) / outputGranularity) + 1);
	saveToFile_dildt.resize((uint64_t)largerDetailUntil + (totalSteps - ((uint)largerDetailUntil) / outputGranularity) + 1);

	rtl Ia = 0;
	for (uint b = (uint)v_Iab.size() - 1; b > 0; --b)
		Ia += v_Iab[b];
	rtl il = 0.0;

	saveToFile_diadt[0] = Ia / numAgents;
	saveToFile_dildt[0] = il;
	bool end = false;
	++outputSize;

	//For the first 'largerDetailUntil' iterations every step is stored in a vector ('saveToFile'), for later being written to file.
	//bool stationary = false;
	for (uint s = 1; s < largerDetailUntil; ++s) {
		rtl prevIA = Ia;
		step(h, Ia, v_Iab, v_Sab);
		if (Ia < epsilon) {
			saveToFile_diadt[outputSize] = 0;
			end = true;
			++outputSize;
			break;
		}
		if (Ia == prevIA) {
			while (s < largerDetailUntil) {
				saveToFile_diadt[outputSize] = Ia / numAgents;
				++s;
				++outputSize;
			}
			end = true;
			break;
		}
		saveToFile_diadt[outputSize] = Ia / numAgents;

#ifdef NORM_SITE_PER_AG
		saveToFile_dildt[outputSize] = il * (N / numAgents);
#else
		saveToFile_dildt[s] = il;
#endif
		++outputSize;
	}
	if (end) return;

	//From the 'largerDetailUntil' iteration on, we afford to ignore 'outputGranularity'-size windows of values, so that the saved file does not grow explosively.
	for (uint s = (uint)largerDetailUntil; s < totalSteps; ++s) {
		rtl prevIA = Ia;
		step(h, Ia, v_Iab, v_Sab);
		if (Ia < epsilon) {
			saveToFile_diadt[outputSize] = 0;
#ifdef NORM_SITE_PER_AG
			saveToFile_dildt[outputSize] = il * (N / numAgents);
#else
			saveToFile_dildt[outputSize] = il;
#endif
			++outputSize;
			break;
		}
		if (Ia == prevIA) {
			while (s < totalSteps) {
				if (s % outputGranularity == 0) {
					saveToFile_diadt[outputSize] = Ia / numAgents;
					++outputSize;
				}
				++s;
			}
			break;
		}
		if (s % outputGranularity == 0) {
			saveToFile_diadt[outputSize] = Ia / numAgents;
#ifdef NORM_SITE_PER_AG
			saveToFile_dildt[outputSize] = il * (N / numAgents);
#else
			saveToFile_dildt[outputSize] = il;
#endif
			++outputSize;
		}
	}
}
#endif //CLIQUE

//rtl sim::diadt(const rtl& Ia, const rtl& il) {
rtl Solver::diadt(const rtl& Ia, const rtl& sumSbIb) {
	//return graph::Graph::psi * TAU_aa * sumSbIb - GAMMA_a * Ia;
	return -nG * Ia;
}
rtl Solver::dildt(const rtl& Ia, const rtl& il) {
	//return  beta_al * (1.0 - il) * Ia - GAMMA_l * il;
	return 0.0;
}

rtl Solver::dIdt(const rtl& Ia) {
	const rtl Sa = (rtl)numAgents - Ia;
	const rtl _b_ = (rtl)graph::Graph::averageDegree;
	const rtl& _b2_ = graph::Graph::_2ndMmt;
	//const uint& B = graph::Graph::B;
	return (2.0 * nL * sigma * w * _b2_ * Sa * Ia ) / (N * pow(_b_, 2.0)) - nG * Ia;	//----> Ótimo no esparso pra qq rede!
	//return (2.0 * nL * sigma * w * Sa * Ia ) / (N) - nG * Ia;	
	//return (2.0 * nL * sigma * w * (44.) * Sa * Ia ) / (N * _b_) - nG * Ia;	

	//rtl avKB = 0;
	//for (uint b = 0; b < graph::Graph::kb.size(); ++b){
	//	const rtl& pb = graph::Graph::block_prob[b];
	//	const rtl& qb = graph::Graph::q_b[b];
	//	const rtl& kb = graph::Graph::kb[b];
	//	if (pb == 0)
	//		continue;
	//	avKB += (kb*qb)/(N*pb) ;
	//}
	//avKB /= graph::Graph::largestDegree;
	//if (avKB < 2.0)
	//	avKB = 2.0;
	//
	//return (avKB * nL * sigma * _b2_ * Sa * Ia) / (N * pow(_b_, 2.0)) - nG * Ia;
	//
}

void Solver::rkMaster(const rtl& t0, std::vector<rtl>& v_Iab, std::vector<rtl>& v_Sab, const rtl& t, const rtl& h, const rtl& epsilon, std::vector<rtl>& saveToFile_diadt, std::vector<rtl>& saveToFile_dildt, uint& outputSize, const uint& outputGranularity, const rtl& largerDetailUntil) {
	uint totalSteps = (uint)((t - t0) / h) + 1;
	saveToFile_diadt.resize((uint64_t)largerDetailUntil + (totalSteps - ((uint)largerDetailUntil) / outputGranularity) + 1);
	saveToFile_dildt.resize((uint64_t)largerDetailUntil + (totalSteps - ((uint)largerDetailUntil) / outputGranularity) + 1);

	rtl Ia = 0;
	for (uint b = (uint)v_Iab.size() - 1; b > 0; --b)
		Ia += v_Iab[b];
	rtl il = 0.0;

	saveToFile_diadt[0] = Ia / numAgents;
	saveToFile_dildt[0] = il;
	bool end = false;
	++outputSize;

	//For the first 'largerDetailUntil' iterations every step is stored in a vector ('saveToFile'), for later being written to file.
	//bool stationary = false;
	for (uint s = 1; s < largerDetailUntil; ++s) {
		rtl prevIA = Ia;
		stepMaster(h, Ia);
		if (Ia < epsilon) {
			saveToFile_diadt[outputSize] = 0;
			end = true;
			++outputSize;
			break;
		}
		if (Ia == prevIA) {
			while (s < largerDetailUntil) {
				saveToFile_diadt[outputSize] = Ia / numAgents;
				++s;
				++outputSize;
			}
			end = true;
			break;
		}
		saveToFile_diadt[outputSize] = Ia / numAgents;

#ifdef NORM_SITE_PER_AG
		saveToFile_dildt[outputSize] = il * (N / numAgents);
#else
		saveToFile_dildt[s] = il;
#endif
		++outputSize;
	}
	if (end) return;

	//From the 'largerDetailUntil' iteration on, we afford to ignore 'outputGranularity'-size windows of values, so that the saved file does not grow explosively.
	for (uint s = (uint)largerDetailUntil; s < totalSteps; ++s) {
		rtl prevIA = Ia;
		stepMaster(h, Ia);
		if (Ia < epsilon) {
			saveToFile_diadt[outputSize] = 0;
#ifdef NORM_SITE_PER_AG
			saveToFile_dildt[outputSize] = il * (N / numAgents);
#else
			saveToFile_dildt[outputSize] = il;
#endif
			++outputSize;
			break;
		}
		if (Ia == prevIA) {
			while (s < totalSteps) {
				if (s % outputGranularity == 0) {
					saveToFile_diadt[outputSize] = Ia / numAgents;
					++outputSize;
				}
				++s;
			}
			break;
		}
		if (s % outputGranularity == 0) {
			saveToFile_diadt[outputSize] = Ia / numAgents;
#ifdef NORM_SITE_PER_AG
			saveToFile_dildt[outputSize] = il * (N / numAgents);
#else
			saveToFile_dildt[outputSize] = il;
#endif
			++outputSize;
		}
	}
}

void Solver::stepMaster(const rtl& h, rtl& Ia) {
	constexpr rtl one_sixth = 1.0 / 6.0;
	const uint blocks = static_cast<uint>(graph::Graph::block_prob.size());
	rtl k1(0), k2(0), k3(0), k4(0);

	k1 = h * dIdt(Ia);
	k2 = h * dIdt(Ia + 0.5 * k1);
	k3 = h * dIdt(Ia + 0.5 * k2);
	k4 = h * dIdt(Ia + k3);

	//Take step:
	Ia += one_sixth * (k1 + 2.0l * k2 + 2.0l * k3 + k4);
}
