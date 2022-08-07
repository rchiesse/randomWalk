#include "solver.h"
#include "graph.h"
using namespace sim;

real Solver::nT = 0;
real Solver::nL = 0;
real Solver::nG = 0;
real Solver::EULER = 0.57721566490153286060651209008240243104215933593992;		// ----> The Euler�Mascheroni constant.
uint Solver::numAgents = 0;

void Solver::setParams(const real& tau, const real& lambda, const real& gamma, const uint& NUM_AGENTS) {
	nT = tau;
	nL = lambda;
	nG = gamma;
	numAgents = NUM_AGENTS;
}

#ifdef CLIQUE
real Solver::diabdt(const real& Ia, const real& Sa) {
	const real pb = 1.0;
	const real nb = graph::Graph::n;
	const real qb = 1.0;
	const real ibnb = Ia / nb;
	const real sbnb = Sa / nb;
	const real kbnb = ibnb + sbnb;
	const real b = N;
	const real K = Ia + Sa;
	const real p = b / (2.0 * graph::Graph::m - N);
	const real E_X = K * p;
	const real E_X2 = K*K * p*p;
	const real Var_X = E_X * (1.0 - p);
	const real stdDev = sqrt(Var_X);
	//const real cs = (1.0 - (sbnb / kbnb) * (ibnb / kbnb)) * stdDev;
	//const real Sa = (double)NUM_AGENTS - Ia;

	//RONALD (BEST SO FAR):
	if (Sa == 0.0)
		return -(nG * Ia);

	//long double H = EULER * log(sbnb + 1.0) / (2.0 * nL);
	//long double H = EULER * log(sbnb * ibnb * (1.0 / nT)) / (2.0 * nL);
	//long double dFactor = (sbnb * ibnb / nT) * std::min(1.0, abs(nT - nL));
	//long double H = EULER * log(ibnb + 1.0 + dFactor) / (2.0 * nL);
	//long double H = EULER * log((ibnb * (std::min(sbnb, ibnb) / (2.0))) + 1.0) / (2.0 * nL);
	//long double H = EULER * log((sbnb)+1.0) / (2.0 * nL) + EULER * log((ibnb)+1.0) / (2.0 * nL);
	//long double ii = ibnb + 1.0;
	//const double prob_inf = ii / (H + ii);
	//const double p = std::min(sbnb, 1.0) * std::min(ibnb, 1.0);
	//real Hs = std::max(1.0l, EULER * log((sbnb)+1.0) / (2.0 * nL));
	//real Hi = std::max(1.0l, EULER * log((ibnb)+1.0) / (2.0 * nL));
	//const real infRate = (sbnb > ibnb) ? (nT / Hs) : nT;
	//const real sigma = (ibnb * nT) / (nL + (ibnb * nT));
	//const real ve = Var_X / E_X;
	//const real ev = E_X / Var_X;
	//const real sd = sqrt(Var_X); // ----> Standard devation.
	//const real sdi = sd * (ibnb / kbnb);
	//const real sds = sd * (sbnb / kbnb);
	//const real dimm = (ve * ibnb * nG) / ((ev * sbnb) * nT + ve * ibnb * nG);
	//const real dimm = (ibnb * nG) / ((sbnb) * nT + ibnb * nG);
	//real ii = sbnb + 1.0;
	//real H = std::max(1.0l, EULER * log((sbnb)+1.0) / (2.0 * nL));
	//const real prob_inf = ii / (H + ii);
	
	//const real prob_inf = ii / (H + ii);
	//real H = (2.0 * nL) / EULER * log(sbnb);
	//real H = EULER * log(sbnb) / (2.0 * nL);
	//real H = std::max(1.0l, EULER * log(sbnb+1.0) / (nL));
	//real H = EULER * log(sbnb);
	//real ii = std::max(ibnb, 1.0l);
	//real H = EULER * log(sbnb+1.0) / (nL);
	//real H = std::max(1.0l, EULER * log(kbnb) / (2.0 * nL));	//LEGAL TB!
	//const real prob_inf = (ii * nT) / (H + (ii * nT));	//RON (figura plot) (aqui a presen�a de nT e nL n faz diferen�a. Avaliar outros cen�rios!)
	//real H = EULER * log(sbnb * ibnb) / (2.0 * nL); //INTERESSANTE...
	//real H = EULER * log(std::min(sbnb, ibnb)) / (2.0 * nL); // MT INTERESSANTE TB...
	real H = EULER * log(kbnb) / (2.0 * nL);	//RON

	real ii = ibnb;
	const real prob_inf = (ii ) / (H + (ii ));
	//const real prob_inf = 1.0;
	return
		//E-M:
		//+ Ia * sbnb * prob_inf * nT
		//- (nG * Ia);

		//POR N�:
		//nb * (
		//	(Ia - ibnb) * nL / nb
		//	- ibnb * nL * (b - 1) / b
		//	+ ibnb * sbnb * nT
		//	- (nG * ibnb)
		//);

		//EQUIVALENTE MASTER:
		((nT * Sa * Ia * prob_inf) / (N)) - nG * Ia;
}

real Solver::dsabdt(const real& Ia, const real& Sa) {
	const real pb = 1.0;
	const real nb = graph::Graph::n;
	const real qb = 1.0;
	const real ibnb = Ia / nb;
	const real sbnb = Sa / nb;
	const real kbnb = ibnb + sbnb;
	const real b = N;
	const real pIn  = (b - 1) / (2.0 * graph::Graph::m - 1.0);		// ----> Probability that a randomly chosen link leads an agent to an specific node v_b, provided that the agent is outside.
	const real pOut = (b - 1) / b;									// ----> Probability of leaving a node v_b, provided that the agent is inside.
	const real K = Ia + Sa;
	const real p = b / (2.0 * graph::Graph::m - N);
	const real E_X = K * p;
	const real E_X2 = K * K * p * p;
	const real Var_X = E_X * (1.0 - p);
	const real stdDev = sqrt(Var_X);
	//const real cs = (1.0 - (sbnb / kbnb) * (ibnb / kbnb)) * stdDev;
	//const real Sa = (double)NUM_AGENTS - Ia;
	//BEST SO FAR:
	if (Sa == 0.0) 
		return nG * Ia;

	//real H = EULER * log(sbnb + ibnb + 1.0) / (2.0 * nL);
	//real H = EULER * log(sbnb * ibnb * (1.0 / nT)) / (2.0 * nL);
	//real dFactor = (sbnb * ibnb / nT) * std::min(1.0, abs(nT - nL));
	//real H = EULER * log(ibnb + 1.0 + dFactor) / (2.0 * nL);
	//real H = EULER * log((ibnb * (std::min(sbnb, ibnb) / (2.0))) + 1.0) / (2.0 * nL);
	//real H = EULER * log((sbnb)+1.0) / (2.0 * nL) + EULER * log((ibnb)+1.0) / (2.0 * nL);
	//real Hi = std::max(1.0l, EULER * log((ibnb)+1.0) / (2.0*nL));
	//const real infRate = (sbnb > ibnb) ? (nT - Hs) : nT;
	//const real sigma = (ibnb * nT) / (nL + (ibnb * nT));
	//const real ve = Var_X / E_X;
	//const real ev = E_X / Var_X;
	//const real dimm = (ibnb * nG) / ((sbnb) * nT + ibnb * nG);
	//real ii = ibnb + 1.0;
	//real H = std::max(1.0l, EULER * log((ibnb)+1.0) / (2.0 * nL));
	//const real prob_inf = ii / (H + ii);

	//real H = (2.0 * nL) / EULER * log(sbnb);
	//real H = std::max(1.0l, EULER * log(sbnb+1.0) / (nL));
	//real H = EULER * log(sbnb);
	//real H = EULER * log(sbnb) / (2.0 * nL);
	//real ii = std::max(ibnb, 1.0l);
	//real H = EULER * log(sbnb+1.0) / (nL);	//�timo qnd LAMBDA >> TAU; SUPERESTIMA C.C (e H = 1.0 acaba sendo melhor)
	//real H = EULER * log(kbnb) / (2.0 * nL);	//RON
	//real H = std::max(1.0l, EULER * log(kbnb) / (2.0 * nL));	//LEGAL TB!
	//real H = EULER * log(ibnb) / (2.0 * nL);
	//const real prob_inf = (ii * nT) / (H + (ii * nT));
	//real H = EULER * log(sbnb * ibnb) / (2.0 * nL);	// INTERESSANTE...
	//real H = EULER * log(std::min(sbnb, ibnb)) / (2.0 * nL); // MT INTERESSANTE TB...
	real H = EULER * log(kbnb) / (2.0 * nL);	//RON
	real ii = ibnb;
	const real prob_inf = (ii ) / (H + (ii ));
	//const real prob_inf = (ii) / (1.0 + (ii));	//�timo qnd TAU == LAMBDA; SUPERESTIMA qnd LAMBDA >> TAU (e EULER acaba sendo melhor)
	//const real prob_inf = 1.0;
	return
		//POR N�:
		//-Ia * sbnb * prob_inf * nT

		//E-M:
		//-Ia * sbnb * prob_inf * nT
		//+ (nG * Ia);

		//OFFICIAL:
		//nb * (
		//	(Sa - sbnb) * nL / nb
		//	- sbnb * nL * (b - 1) / b
		//	- ibnb * sbnb * nT
		//	+ (nG * ibnb)
		//);

		//EQUIVALENTE MASTER:
		-((nT * Sa * Ia * prob_inf) / (N)) + nG * Ia;
}

void Solver::step(const real& h, real& Ia, real& Sa) {
	constexpr real one_sixth = 1.0 / 6.0;
	std::vector<real> k1(2, 0), k2(2, 0), k3(2, 0), k4(2, 0);

	lookAhead(h, Ia, Sa, k1);
	lookAhead(h, Ia, Sa, k2, k1, 0.5);
	lookAhead(h, Ia, Sa, k3, k2, 0.5);
	lookAhead(h, Ia, Sa, k4, k3);

	//Take step:
	Ia += one_sixth * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
	Sa += one_sixth * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
}

void Solver::lookAhead(const real& h, real& Ia, real& Sa, std::vector<real>& target) {
	target[0] = h * diabdt(Ia, Sa);
	target[1] = h * dsabdt(Ia, Sa);
}

void Solver::lookAhead(const real& h, real& Ia, real& Sa, std::vector<real>& target, std::vector<real>& base, const real& fraction) {
	target[0] = h * diabdt(Ia + fraction * base[0], Sa + fraction * base[1]);
	target[1] = h * dsabdt(Ia + fraction * base[0], Sa + fraction * base[1]);
}

void Solver::rungeKutta4thOrder(const real& t0, real& Ia, real& Sa, const real& t, const real& h, const real& epsilon, std::vector<real>& saveToFile_diadt, std::vector<real>& saveToFile_dildt, uint& outputSize, const uint& outputGranularity, const real& largerDetailUntil) {
	uint totalSteps = (uint)((t - t0) / h) + 1;
	saveToFile_diadt.resize((uint64_t)largerDetailUntil + (totalSteps - ((uint)largerDetailUntil) / outputGranularity) + 1);

	saveToFile_diadt[0] = Ia / numAgents;
	bool end = false;
	++outputSize;

	real sumH = 0;
	//For the first 'largerDetailUntil' iterations every step is stored in a vector ('saveToFile'), for later being written to file.
	for (uint s = 1; s < largerDetailUntil; ++s) {
		real prevIA = Ia;
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
		real prevIA = Ia;
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
real Solver::diabdt(const real& Ia, const real& Iab, const real& Sab, const uint& b) {
	const real& pb = graph::Graph::block_prob[b];
	const real nb = graph::Graph::n * pb;
	const real& qb = graph::Graph::q_b[b];
	const real ibnb = Iab / nb;
	const real sbnb = Sab / nb;
	const real kbnb = ibnb + sbnb;
	const real Sa = (real)numAgents - Ia;
	const real l = ((real)b - 1) / ((real)b);

	//RONALD (BEST SO FAR):
	//if (Sab == 0.0) {
	//	return Ia * nL * qb - Iab * nL
	//		- (nG * Iab);
	//}
	//real H = EULER * log(sbnb) / (2.0 * nL);
	real H = EULER * log(kbnb) / (2.0 * nL);	//RON
	real ii = ibnb;
	const real prob_inf = ii / (H + ii);
	return //Ia * nL * qb - Iab * nL
		+ Iab * sbnb * prob_inf * nT
		- (nG * Iab);

	//NOVO TESTE (denso):
	//VERS�O DE BLOCO EQUIVALENTE AO DE N�:
	//return nb * (
	//	(Ia - ibnb) * nL * (((real)b - 1.0) / (2.0 * graph::Graph::m - N - (b - 1))) - ibnb * nL * l
	//	+ ibnb * sbnb * nT
	//	- (nG * ibnb)
	//);
	//VERS�O DE N� EQUIVALENTE AO DE BLOCO:
	//return Ia * nL * qb - Iab * nL
	//	+ Iab * sbnb * nT
	//	- (nG * Iab);

	//ISOLANDO CADA BLOCO:
	//return nb * (
	//	+ ibnb * sbnb * nT
	//	- (nG * ibnb)
	//);
}

real Solver::dsabdt(const real& Ia, const real& Iab, const real& Sab, const uint& b) {
	const real& pb = graph::Graph::block_prob[b];
	const real nb = graph::Graph::n * pb;
	const real& qb = graph::Graph::q_b[b];
	const real ibnb = Iab / nb;
	const real sbnb = Sab / nb;
	const real kbnb = ibnb + sbnb;
	const real Sa = (real)numAgents - Ia;
	const real l = ((real)b - 1) / ((real)b);

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

	////RONALD (excelente quando LAMBDA == TAU, mesmo em cen�rios densos):
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
	//real H = EULER * log(sbnb) / (2.0 * nL);
	real H = EULER * log(kbnb) / (2.0 * nL);	//RON
	real ii = ibnb;
	const real prob_inf = ii / (H + ii);
	return //Sa * nL * qb - Sab * nL
		- Iab * sbnb * prob_inf * nT
		+ (nG * Iab);

	//DON : Good for dense scenarios or large walk rate
	//if (kb == 0.0)
	//	return Sa * LAMBDA * qb - Sab * LAMBDA;
	//return 
	//	Sa * LAMBDA * qb - Sab * LAMBDA
	//	- Iab * sbnb * TAU_aa
	//	+ (GAMMA_a * Iab);

	//NOVO TESTE (denso):
	//VERS�O DE BLOCO EQUIVALENTE AO DE N�:
	//return Sa * nL * qb - Sab * nL
	//	- Iab * sbnb * nT
	//	+ (nG * Iab);
	//VERS�O DE N� EQUIVALENTE AO DE BLOCO:
	//return nb * (
	//	(Sa - sbnb) * nL * (((real)b - 1.0) / (2.0 * graph::Graph::m - N - (b-1))) - sbnb * nL * l
	//	- ibnb * sbnb * nT
	//	+ (nG * ibnb)
	//);

	//ISOLANDO CADA BLOCO:
	//return nb * (
	//	- ibnb * sbnb * nT
	//	+ (nG * ibnb)
	//);
}

void Solver::step(const real& h, real& Ia, std::vector<real>& v_Iab, std::vector<real>& v_Sab) {
	constexpr real one_sixth = 1.0 / 6.0;
	const uint blocks = static_cast<uint>(graph::Graph::block_prob.size());
	std::vector<real> k1(2 * blocks, 0), k2(2 * blocks, 0), k3(2 * blocks, 0), k4(2 * blocks, 0);

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

void Solver::lookAhead(const real& h, real& Ia, const std::vector<real>& v_Iab, const std::vector<real>& v_Sab, std::vector<real>& target) {
	update_Ia(Ia, v_Iab);
	for (uint b = (uint)graph::Graph::block_prob.size() - 1; b > 0; --b) {
		if (graph::Graph::block_prob[b] == 0)
			continue;

		target[2 * b]		= h * diabdt(Ia, v_Iab[b], v_Sab[b], b);
		target[2 * b + 1]	= h * dsabdt(Ia, v_Iab[b], v_Sab[b], b);
	}
}

void Solver::lookAhead(const real& h, real& Ia, const std::vector<real>& v_Iab, const std::vector<real>& v_Sab, std::vector<real>& target, std::vector<real>& base, const double& fraction) {
	update_Ia(Ia, v_Iab, base, fraction);
	for (uint b = (uint)graph::Graph::block_prob.size() - 1; b > 0; --b) {
		if (graph::Graph::block_prob[b] == 0)
			continue;

		target[2 * b]		= h * diabdt(Ia, v_Iab[b] + fraction * base[2 * b], v_Sab[b] + fraction * base[2 * b + 1], b);
		target[2 * b + 1]	= h * dsabdt(Ia, v_Iab[b] + fraction * base[2 * b], v_Sab[b] + fraction * base[2 * b + 1], b);
	}
}

void Solver::update_Ia(real& Ia, const std::vector<real>& v_Iab) {
	Ia = 0;
	for (uint b = (uint)graph::Graph::block_prob.size() - 1; b > 0; --b)
		Ia += v_Iab[b];
}

void Solver::update_Ia(real& Ia, const std::vector<real>& v_Iab, const std::vector<real>& base, const double& fraction) {
	Ia = 0;
	for (uint b = (uint)graph::Graph::block_prob.size() - 1; b > 0; --b)
		Ia += v_Iab[b] + (fraction * base[2 * b]);
}
#else //PER_BLOCK
real Solver::divbdt(const real& Ia, const real& Iv, const real& Sv, const uint& b) {
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
real Solver::dsvbdt(const real& Ia, const real& Iv, const real& Sv, const uint& b) {
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
void Solver::step(const real& h, real& Ia, std::vector<real>& v_Iv, std::vector<real>& v_Sv) {
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
void Solver::lookAhead(const real& h, real& Ia, const std::vector<real>& v_Iv, const std::vector<real>& v_Sv, std::vector<real>& target) {
	update_Ia(Ia, v_Iv);
	for (int v = (uint)graph::Graph::n - 1; v >= 0; --v) {
		target[2 * v] = h * divbdt(Ia, v_Iv[v], v_Sv[v], (uint)graph::Graph::g[v].size());
		target[2 * v + 1] = h * dsvbdt(Ia, v_Iv[v], v_Sv[v], (uint)graph::Graph::g[v].size());
	}
}
void Solver::lookAhead(const real& h, real& Ia, const std::vector<real>& v_Iv, const std::vector<real>& v_Sv, std::vector<real>& target, std::vector<real>& base, const double& fraction) {
	update_Ia(Ia, v_Iv, base, fraction);
	for (int v = (uint)graph::Graph::n - 1; v >= 0; --v) {
		target[2 * v] = h * divbdt(Ia, v_Iv[v] + fraction * base[2 * v], v_Sv[v] + fraction * base[2 * v + 1], (uint)graph::Graph::g[v].size());
		target[2 * v + 1] = h * dsvbdt(Ia, v_Iv[v] + fraction * base[2 * v], v_Sv[v] + fraction * base[2 * v + 1], (uint)graph::Graph::g[v].size());
	}
}
void Solver::update_Ia(real& Ia, const std::vector<real>& v_Iv) {
	Ia = 0;
	for (int v = (uint)graph::Graph::n - 1; v >= 0; --v)
		Ia += v_Iv[v];
}
void Solver::update_Ia(real& Ia, const std::vector<real>& v_Iv, const std::vector<real>& base, const double& fraction) {
	Ia = 0;
	for (int v = (uint)graph::Graph::n - 1; v >= 0; --v)
		Ia += v_Iv[v] + (fraction * base[2 * v]);
}
#endif //PER_BLOCK

real Solver::dilbdt(const real& Ia, const real& il, const real& Iab, const real& ilb, const uint& b) {
	return 0;
}

void Solver::rungeKutta4thOrder(const real& t0, std::vector<real>& v_Iab, std::vector<real>& v_Sab, std::vector<real>& v_ilb, const real& t, const real& h, const real& epsilon, std::vector<real>& saveToFile_diadt, std::vector<real>& saveToFile_dildt, uint& outputSize, const uint& outputGranularity, const real& largerDetailUntil) {
	uint totalSteps = (uint)((t - t0) / h) + 1;
	saveToFile_diadt.resize((uint64_t)largerDetailUntil + (totalSteps - ((uint)largerDetailUntil) / outputGranularity) + 1);
	saveToFile_dildt.resize((uint64_t)largerDetailUntil + (totalSteps - ((uint)largerDetailUntil) / outputGranularity) + 1);

	real Ia = 0;
	for (uint b = (uint)v_Iab.size() - 1; b > 0; --b)
		Ia += v_Iab[b];
	real il = 0.0;

	saveToFile_diadt[0] = Ia / numAgents;
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
		real prevIA = Ia;
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

//real sim::diadt(const real& Ia, const real& il) {
real Solver::diadt(const real& Ia, const real& sumSbIb) {
	//return graph::Graph::psi * TAU_aa * sumSbIb - GAMMA_a * Ia;
	return -nG * Ia;
}
real Solver::dildt(const real& Ia, const real& il) {
	//return  beta_al * (1.0 - il) * Ia - GAMMA_l * il;
	return 0.0;
}