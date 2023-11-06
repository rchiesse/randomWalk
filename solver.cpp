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
rtl Solver::master_esigma = 0;
uint Solver::numAgents = 0;
rtl Solver::master_walkRate = 0;
std::vector<rtl> Solver::block_walkRate = {0};

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

void Solver::setBlockData() {
	block_walkRate	.resize(graph::Graph::block_prob.size(), 0);
	
	for (const auto& b : graph::Graph::validBlocks) 
		block_walkRate[b] = (((rtl)(b - 1)) / b) * nL;
	
	for (const auto& b : graph::Graph::validBlocks)
		master_walkRate += block_walkRate[b] * graph::Graph::q_b[b];
		//master_walkRate += block_walkRate[b] * graph::Graph::block_prob[b];

	//rtl avEscapeRate = 0;
	for (const auto& b : graph::Graph::validBlocks) {
		const rtl out = (rtl)(b - 1);
		//const rtl escapeRate = (out / (out + w)) * nL;
		const rtl escapeRate = (out / b) * nL;
		master_esigma += (nT / (2 * escapeRate + nT)) * graph::Graph::q_b[b];
		//avEscapeRate += ((out / (out + w)) * nL) * graph::Graph::block_prob[b];
	}
	//master_esigma += nT / (2 * avEscapeRate + nT);
}

#ifdef PER_BLOCK
rtl Solver::diabdt(const rtl& Ia, const rtl& Iab, const rtl& Sab, const uint& b) {
	const rtl nb = graph::Graph::n * graph::Graph::block_prob[b];
	const rtl& qb = graph::Graph::q_b[b];
	const rtl Sa = (rtl)numAgents - Ia;
	const rtl out = (rtl)(b - 1);

	//const rtl safe_s = (graph::Graph::n - Ia) / graph::Graph::n;
	//const rtl hostile_s = Ia / graph::Graph::n;
	//const rtl safe_i = (graph::Graph::n - Sa) / graph::Graph::n;
	//const rtl hostile_i = Sa / graph::Graph::n;

	//const rtl escapeRate_s = (1.0 - (Ws / (out * (safe_s + hostile_s * Ws) + Ws))) * nL;
	//const rtl escapeRate_i = (1.0 - (Wi / (out * (safe_i + hostile_i * Wi) + Wi))) * nL;
	//const rtl escapeRate = ((out * safe + out * hostile * w) / (out * safe + out * hostile * w + w)) * nL;
	//const rtl esigma = nT / (escapeRate_i + escapeRate_s + nT);
	//const rtl escapeRate = (out / (out + w)) * nL;
	
	const rtl escapeRate = (out / b) * nL;
	const rtl etau = (nT / (escapeRate + nT)) * nT;
	//const rtl ii = Iab / nb;
	//const rtl esigma = etau / (escapeRate + ii * etau);
	const rtl esigma = nT / (2 * escapeRate + nT);
	
#ifdef SELF_LOOPS
	return nL * (Ia * qb - Iab)
		+ 2.0 * ((Sab * Iab) / nb) * escapeRate * esigma * w - nG * Iab;
		//+ ((Sab * Iab) / nb) * escapeRate * ii * esigma * w
		//+ ((Sab * Iab) / nb) * escapeRate * esigma * w 
		//- nG * Iab;
#else
	return nL * (Ia * qb - Iab)
		+ 2.0 * ((Sab * Iab) / nb) * nL * sigma * w - (nG * Iab);
#endif
	
	//TESTE DENSO - Bom apenas quando endemia é menor q 50% da população:
	//rtl sbnb = Sab / nb;
	//rtl ibnb = Iab / nb;
	////rtl maxL = (ibnb < 1.0) ? nL : ibnb*nL;
	//
	//rtl probinf = nT / ((ibnb * nL) + nT);
	//rtl probacq = ibnb * (nT / ((ibnb * nL) + nT));
	//rtl inf2 = sbnb * ibnb * probinf * probinf;
	//rtl acq2 = sbnb * ibnb * probacq * probacq
	//	+ 2.0 * sbnb * ibnb * probacq * probinf
	//	+ sbnb * ibnb * probinf * probinf;
	//
	//return nL * (Ia * qb - Iab)
	//	+ ((Sab * Iab) / nb) * nL * probinf * w
	//	+ ((Sab * Iab) / nb) * nL * probacq * w
	//	+ acq2
	//	+ acq2
	//	- (nG * Iab);
}

rtl Solver::dsabdt(const rtl& Ia, const rtl& Iab, const rtl& Sab, const uint& b) {
	const rtl nb = graph::Graph::n * graph::Graph::block_prob[b];
	const rtl& qb = graph::Graph::q_b[b];
	const rtl Sa = (rtl)numAgents - Ia;
	const rtl out = (rtl)(b - 1);
	
	//const rtl safe_s = (graph::Graph::n - Ia) / graph::Graph::n;
	//const rtl hostile_s = Ia / graph::Graph::n;
	//const rtl safe_i = (graph::Graph::n - Sa) / graph::Graph::n;
	//const rtl hostile_i = Sa / graph::Graph::n;

	//const rtl escapeRate_s = (1.0 - (Ws / (out * (safe_s + hostile_s * Ws) + Ws))) * nL;
	//const rtl escapeRate_i = (1.0 - (Wi / (out * (safe_i + hostile_i * Wi) + Wi))) * nL;
	//const rtl escapeRate = ((out * safe + out * hostile * w) / (out * safe + out * hostile * w + w)) * nL;
	//const rtl esigma = nT / (escapeRate_i + escapeRate_s + nT);
	//const rtl escapeRate = (out / (out + w)) * nL;
	
	const rtl escapeRate = (out / b) * nL;
	const rtl etau = (nT / (escapeRate + nT)) * nT;
	//const rtl ii = Iab / nb;
	//const rtl esigma = etau / (escapeRate + ii * etau);
	const rtl esigma = nT / (2 * escapeRate + nT);

#ifdef SELF_LOOPS
	return nL * (Sa * qb - Sab)
		- 2.0 * ((Sab * Iab) / nb) * escapeRate * esigma * w + nG * Iab;
		//- ((Sab * Iab) / nb) * escapeRate * ii * esigma * w
		//- ((Sab * Iab) / nb) * escapeRate * esigma * w 
		//+ nG * Iab;
#else
	return nL * (Sa * qb - Sab)
		- 2.0 * ((Sab * Iab) / nb) * nL * sigma * w + (nG * Iab);
#endif

	
	 
	//TESTE DENSO - Bom apenas quando endemia é menor q 50% da população:
	//rtl sbnb = Sab / nb;
	//rtl ibnb = Iab / nb;
	////rtl maxL = (ibnb < 1.0) ? nL : ibnb * nL;
	//
	//rtl probinf = nT / ((ibnb * nL)+nT);
	//rtl probacq = ibnb * (nT / ((ibnb * nL) + nT));
	//rtl inf2 = sbnb * ibnb * probinf * probinf;
	//rtl acq2 = sbnb * ibnb * probacq * probacq
	//	+ 2.0 * sbnb * ibnb * probacq * probinf
	//	+ sbnb * ibnb * probinf * probinf;
	//
	//return nL * (Sa * qb - Sab)
	//	- ((Sab * Iab) / nb) * nL * probinf * w
	//	- ((Sab * Iab) / nb) * nL * probacq * w
	//	- acq2
	//	- acq2
	//	+ (nG * Iab);
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
	uint b = (uint)graph::Graph::block_prob.size();
	while (--b) {
		if (graph::Graph::block_prob[b] == 0)
			continue;

		v_Iab[b] += one_sixth * (k1[2 * b] + 2 * k2[2 * b] + 2 * k3[2 * b] + k4[2 * b]);
		v_Sab[b] += one_sixth * (k1[2 * b + 1] + 2 * k2[2 * b + 1] + 2 * k3[2 * b + 1] + k4[2 * b + 1]);
	}
}

void Solver::lookAhead(const rtl& h, rtl& Ia, const std::vector<rtl>& v_Iab, const std::vector<rtl>& v_Sab, std::vector<rtl>& target) {
	update_Ia(Ia, v_Iab);
	uint b = (uint)graph::Graph::block_prob.size();
	while (--b) {
		if (graph::Graph::block_prob[b] == 0)
			continue;

		target[2 * b]		= h * diabdt(Ia, v_Iab[b], v_Sab[b], b);
		target[2 * b + 1]	= h * dsabdt(Ia, v_Iab[b], v_Sab[b], b);
	}
}

void Solver::lookAhead(const rtl& h, rtl& Ia, const std::vector<rtl>& v_Iab, const std::vector<rtl>& v_Sab, std::vector<rtl>& target, std::vector<rtl>& base, const double& fraction) {
	update_Ia(Ia, v_Iab, base, fraction);
	uint b = (uint)graph::Graph::block_prob.size();
	while (--b) {
		if (graph::Graph::block_prob[b] == 0)
			continue;

		target[2 * b]		= h * diabdt(Ia, v_Iab[b] + fraction * base[2 * b], v_Sab[b] + fraction * base[2 * b + 1], b);
		target[2 * b + 1]	= h * dsabdt(Ia, v_Iab[b] + fraction * base[2 * b], v_Sab[b] + fraction * base[2 * b + 1], b);
	}
}

void Solver::update_Ia(rtl& Ia, const std::vector<rtl>& v_Iab) {
	Ia = 0;
	uint b = (uint)graph::Graph::block_prob.size();
	while (--b)
		Ia += v_Iab[b];
}

void Solver::update_Ia(rtl& Ia, const std::vector<rtl>& v_Iab, const std::vector<rtl>& base, const double& fraction) {
	Ia = 0;
	uint b = (uint)graph::Graph::block_prob.size();
	while (--b)
		Ia += v_Iab[b] + (fraction * base[2 * b]);
}
#else //PER_BLOCK
rtl Solver::divbdt(const rtl& Ia, const rtl& Iv, const rtl& Sv, const uint& b) {
	//const double& pb = graph::Graph::block_prob[b];
	//const double nb = graph::Graph::n * pb;
	//const double& qb = graph::Graph::q_b[b];
	//const double qbnb = qb / nb;
	//const double Sa = (double)NUM_AGENTS - Ia;
	//const double prob_inf = (TAU_aa / (2.0 * LAMBDA + TAU_aa));
	//return ((Ia - Iv) * LAMBDA * qbnb - Iv * LAMBDA * (1.0 - qbnb))
	//	+ (Sa - Sv) * LAMBDA * qbnb * Iv * prob_inf
	//	+ (Ia - Iv) * LAMBDA * qbnb * Sv * prob_inf
	//	- (GAMMA_a * Iv);
}
rtl Solver::dsvbdt(const rtl& Ia, const rtl& Iv, const rtl& Sv, const uint& b) {
	//const double& pb = graph::Graph::block_prob[b];
	//const double nb = graph::Graph::n * pb;
	//const double& qb = graph::Graph::q_b[b]; const double qbnb = qb / nb;
	//const double Sa = (double)NUM_AGENTS - Ia;
	//const double prob_inf = (TAU_aa / (2.0 * LAMBDA + TAU_aa));
	//return ((Sa - Sv) * LAMBDA * qbnb - Sv * LAMBDA * (1.0 - qbnb))
	//	- (Sa - Sv) * LAMBDA * qbnb * Iv * prob_inf
	//	- (Ia - Iv) * LAMBDA * qbnb * Sv * prob_inf
	//	+ (GAMMA_a * Iv);
}
void Solver::step(const rtl& h, rtl& Ia, std::vector<rtl>& v_Iv, std::vector<rtl>& v_Sv) {
	constexpr rtl one_sixth = 1.0 / 6.0;
	const uint blocks = static_cast<uint>(graph::Graph::block_prob.size());
	std::vector<rtl> k1(2 * graph::Graph::n, 0), k2(2 * graph::Graph::n, 0), k3(2 * graph::Graph::n, 0), k4(2 * graph::Graph::n, 0);

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
		target[2 * v] = h * divbdt(Ia, v_Iv[v], v_Sv[v], (uint)graph::Graph::gs[v].size());
		target[2 * v + 1] = h * dsvbdt(Ia, v_Iv[v], v_Sv[v], (uint)graph::Graph::gs[v].size());
	}
}
void Solver::lookAhead(const rtl& h, rtl& Ia, const std::vector<rtl>& v_Iv, const std::vector<rtl>& v_Sv, std::vector<rtl>& target, std::vector<rtl>& base, const double& fraction) {
	update_Ia(Ia, v_Iv, base, fraction);
	for (int v = (uint)graph::Graph::n - 1; v >= 0; --v) {
		target[2 * v] = h * divbdt(Ia, v_Iv[v] + fraction * base[2 * v], v_Sv[v] + fraction * base[2 * v + 1], (uint)graph::Graph::gs[v].size());
		target[2 * v + 1] = h * dsvbdt(Ia, v_Iv[v] + fraction * base[2 * v], v_Sv[v] + fraction * base[2 * v + 1], (uint)graph::Graph::gs[v].size());
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


void Solver::rungeKutta4thOrder(const rtl& t0, std::vector<rtl>& v_Iab, std::vector<rtl>& v_Sab, std::vector<rtl>& v_ilb, const rtl& t, const rtl& h, const rtl& epsilon, std::vector<rtl>& saveToFile_diadt, std::vector<rtl>& saveToFile_dildt, uint& outputSize, const uint& outputGranularity, const rtl& largerDetailUntil) {
	uint totalSteps = (uint)((t - t0) / h) + 1;
	saveToFile_diadt.resize((uint64_t)largerDetailUntil + (totalSteps - ((uint)largerDetailUntil) / outputGranularity) + 1);
	//saveToFile_dildt.resize((uint64_t)largerDetailUntil + (totalSteps - ((uint)largerDetailUntil) / outputGranularity) + 1);

	rtl Ia = 0;
	for (uint b = (uint)v_Iab.size() - 1; b > 0; --b)
		Ia += v_Iab[b];
	rtl il = 0.0;

	saveToFile_diadt[0] = Ia / numAgents;
	//saveToFile_dildt[0] = il;
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

//#ifdef NORM_SITE_PER_AG
//		saveToFile_dildt[outputSize] = il * (N / numAgents);
//#else
//		saveToFile_dildt[s] = il;
//#endif
		++outputSize;
	}
	if (end) return;

	//From the 'largerDetailUntil' iteration on, we afford to ignore 'outputGranularity'-size windows of values, so that the saved file does not grow explosively.
	for (uint s = (uint)largerDetailUntil; s < totalSteps; ++s) {
		rtl prevIA = Ia;
		step(h, Ia, v_Iab, v_Sab);
		if (Ia < epsilon) {
			saveToFile_diadt[outputSize] = 0;
//#ifdef NORM_SITE_PER_AG
//			saveToFile_dildt[outputSize] = il * (N / numAgents);
//#else
//			saveToFile_dildt[outputSize] = il;
//#endif
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
//#ifdef NORM_SITE_PER_AG
//			saveToFile_dildt[outputSize] = il * (N / numAgents);
//#else
//			saveToFile_dildt[outputSize] = il;
//#endif
			++outputSize;
		}
	}
}

rtl Solver::dIdt(const rtl& Ia) {
	const rtl Sa = (rtl)numAgents - Ia;
	const rtl _b_ = (rtl)graph::Graph::averageDegree;
	const rtl& _b2_ = graph::Graph::_2ndMmt;
	
#ifdef SELF_LOOPS
	return (2.0 * master_walkRate * master_esigma * w * _b2_ * Sa * Ia) / (N * pow(_b_, 2.0)) - nG * Ia;
#else
	return (2.0 * nL * sigma * w * _b2_ * Sa * Ia ) / (N * pow(_b_, 2.0)) - nG * Ia;
#endif
}

void Solver::rkMaster(const rtl& t0, std::vector<rtl>& v_Iab, std::vector<rtl>& v_Sab, const rtl& t, const rtl& h, const rtl& epsilon, std::vector<rtl>& saveToFile_diadt, std::vector<rtl>& saveToFile_dildt, uint& outputSize, const uint& outputGranularity, const rtl& largerDetailUntil) {
	uint totalSteps = (uint)((t - t0) / h) + 1;
	saveToFile_diadt.resize((uint64_t)largerDetailUntil + (totalSteps - ((uint)largerDetailUntil) / outputGranularity) + 1);
	//saveToFile_dildt.resize((uint64_t)largerDetailUntil + (totalSteps - ((uint)largerDetailUntil) / outputGranularity) + 1);

	rtl Ia = 0;
	for (uint b = (uint)v_Iab.size() - 1; b > 0; --b)
		Ia += v_Iab[b];
	rtl il = 0.0;

	saveToFile_diadt[0] = Ia / numAgents;
	//saveToFile_dildt[0] = il;
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

//#ifdef NORM_SITE_PER_AG
//		saveToFile_dildt[outputSize] = il * (N / numAgents);
//#else
//		saveToFile_dildt[s] = il;
//#endif
		++outputSize;
	}
	if (end) return;

	//From the 'largerDetailUntil' iteration on, we afford to ignore 'outputGranularity'-size windows of values, so that the saved file does not grow explosively.
	for (uint s = (uint)largerDetailUntil; s < totalSteps; ++s) {
		rtl prevIA = Ia;
		stepMaster(h, Ia);
		if (Ia < epsilon) {
			saveToFile_diadt[outputSize] = 0;
//#ifdef NORM_SITE_PER_AG
//			saveToFile_dildt[outputSize] = il * (N / numAgents);
//#else
//			saveToFile_dildt[outputSize] = il;
//#endif
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
//#ifdef NORM_SITE_PER_AG
//			saveToFile_dildt[outputSize] = il * (N / numAgents);
//#else
//			saveToFile_dildt[outputSize] = il;
//#endif
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
