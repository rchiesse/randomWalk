#include "graph.h"
#include "reporter.h"
using namespace graph;

std::unordered_map<uint, uint> Graph::idMap;
uint Graph::n;													// ----> Network size, i.e. the number of nodes.
#ifdef CLIQUE
#ifdef PROTECTION_FX
vector<node> Graph::gs;
vector<node> Graph::gi;
#else
vector<node> Graph::g;
#endif //PROTECTION_FX
#else  //CLIQUE
#ifdef PROTECTION_FX
vector<vector<node>> Graph::gs;									// ----> The graph, as an edge list.
vector<vector<node>> Graph::gi;									// ----> The graph, as an edge list.
#else
vector<vector<node>> Graph::g;									// ----> The graph, as an edge list.
#endif //PROTECTION_FX
uint Graph::m;												
uint Graph::largestDegree;
uint Graph::smallestDegree;
uint Graph::selfLoops;
uint Graph::largestDgNode;
uint Graph::lccSize;											// ---> The size of the Largest Connected Component (LCC).
vector<node> Graph::lcc;										// ---> List of nodes that belong to the LCC.
real Graph::averageDegree;
real Graph::originalAvDeg;
real Graph::original2ndMmt;
real Graph::_2ndMmt;
vector<double> Graph::block_prob;
vector<double> Graph::q_b;
vector<double> Graph::originalFreq;
vector<double> Graph::kb;
vector<double> Graph::rho_b;
//vector<double> Graph::rho_bs;
//vector<double> Graph::rho_bi;
//real Graph::sumKB;
real Graph::psi;
real Graph::validBlocks;
#endif //CLIQUE

#ifdef PROTECTION_FX
#ifdef CLIQUE
uint			Graph::_schema_s;
uint			Graph::_schema_i;
vector<real>	Graph::_ps;
vector<real>	Graph::_pi;
vector<uint>	Graph::_myIdx_s;
vector<uint>	Graph::_myIdx_i;
#else
vector<uint>			Graph::schema_s;
vector<uint>			Graph::schema_i;
vector<vector<real>>	Graph::ps;
vector<vector<real>>	Graph::pi;
vector<vector<uint>>	Graph::myForeignIdx_s;
vector<vector<uint>>	Graph::myForeignIdx_i;
#endif //CLIQUE

void Graph::setProbs() {
#ifdef CLIQUE
	_ps.resize((size_t)n + 1);
	_pi.resize((size_t)n + 1);
	vector<real>& probsS = _ps;	// ----> It is important to create the reference only AFTER _ps BEING RESIZED. The need for contiguous space may lead '_ps' to be reallocated in memory upon resizing, what in turn would invalidate any previously created reference (as it would still be pointing to _ps's old address).
	vector<real>& probsI = _pi;	// ----> It is important to create the reference only AFTER _pi BEING RESIZED. The need for contiguous space may lead '_pi' to be reallocated in memory upon resizing, what in turn would invalidate any previously created reference (as it would still be pointing to _pi's old address).
	real sumW = 0;
	for (uint schema = 0; schema < probsS.size(); ++schema) {
		probsS[schema] = sumW / (n - schema + sumW);
		sumW += sim::Ws;
	}
	sumW = 0;
	for (uint schema = 0; schema < probsI.size(); ++schema) {
		probsI[schema] = sumW / (n - schema + sumW);
		sumW += sim::Wi;
	}
#else
	ps.resize(n);
	pi.resize(n);
	for (node v = 0; v < n; ++v) ps[v].resize((size_t)gs[v].size() + 1);
	for (node v = 0; v < n; ++v) pi[v].resize((size_t)gi[v].size() + 1);

	for (node v = 0; v < n; ++v) {
		vector<real>& vProbsS = ps[v];
		vector<real>& vProbsI = pi[v];
		const real vDegree = (real)(gs[v].size());
		real sumW = 0;
		for (uint schema = 0; schema < vProbsS.size(); ++schema) {
			vProbsS[schema] = sumW / (vDegree - schema + sumW);
			sumW += sim::Ws;
		}
		sumW = 0;
		for (uint schema = 0; schema < vProbsI.size(); ++schema) {
			vProbsI[schema] = sumW / (vDegree - schema + sumW);
			sumW += sim::Wi;
		}
	}
#endif //#ifdef CLIQUE
}

void Graph::resetSchema() {
#ifdef CLIQUE
	_schema_s = 0;
	_schema_i = 0;
#else
	for (size_t i = 0; i < schema_s.size(); ++i) schema_s[i] = 0;
	for (size_t i = 0; i < schema_i.size(); ++i) schema_i[i] = 0;
#endif
}
#ifdef CLIQUE
void Graph::resetAgentIdx() {
	_myIdx_s = Graph::gs;
	_myIdx_i = Graph::gi;
}
void Graph::raiseSchema	(uint& _schema) { ++_schema; }
void Graph::lowerSchema	(uint& _schema) { --_schema; }
void Graph::updateHasI(const node& v) {
	updateSBound(v);
	raiseSchema(_schema_s);
}
void Graph::updateNoI(const node& v) {
	//Opposite to 'updateHasI()', here we must first modify v's neighbors' schema and only then update their respective bounds:
	lowerSchema(_schema_s);
	updateSBound(v);
}
void Graph::updateHasS(const node& v) {
	updateIBound(v);
	raiseSchema(_schema_i);
}
void Graph::updateNoS(const node& v) {
	//Opposite to 'updateHasS()', here we must first modify v's neighbors' schema and only then update their respective bounds:
	lowerSchema(_schema_i);
	updateIBound(v);
}
void Graph::updateSBound(const node& v) {
	//Swap v and u. Note that u is the node at the clique's schema bound, to be swapped with v on 'Graph::g' vector. Also, v and u exchange their indexes and update their respective '_myIdx_s' accordingly.
	uint& vIdx = _myIdx_s[v];
	const node u = gs[_schema_s];
	gs[_schema_s] = v;
	gs[vIdx] = u;
	_myIdx_s[u] = vIdx;
	vIdx = _schema_s;
}
void Graph::updateIBound(const node& v) {
	//Swap v and u. Note that u is the node at the clique's schema bound, to be swapped with v on 'Graph::g' vector. Also, v and u exchange their indexes and update their respective '_myIdx_i' accordingly.
	uint& vIdx = _myIdx_i[v];
	const node u = gi[_schema_i];
	gi[_schema_i] = v;
	gi[vIdx] = u;
	_myIdx_i[u] = vIdx;
	vIdx = _schema_i;
}
#ifdef PROPORTIONAL
const node& Graph::nextNodeForS(const real& randBase, const real& itotal, const real& stotal) {
	real ps = pow(1.0 - (itotal / (itotal + stotal)), sim::_r) * _schema_s;
	ps /= (ps + (Graph::n - _schema_s));
#else
const node& Graph::nextNodeForS(const real& randBase) {
	const real& ps = _ps[_schema_s];
#endif
	return (randBase < ps) ? gs[(size_t)(floor((randBase * _schema_s) / ps))] : gs[(size_t)floor(  randBase * (gs.size() - _schema_s)  ) + _schema_s];
}
#ifdef PROPORTIONAL
const node& Graph::nextNodeForI(const real& randBase, const real& itotal, const real& stotal) {
	real pi = pow(1.0 - (itotal / (itotal + stotal)), sim::_r) * _schema_i;
	pi /= (pi + (Graph::n - _schema_s));
#else
const node& Graph::nextNodeForI(const real& randBase) {
	const real& pi = _pi[_schema_i];
#endif
	return (randBase < pi) ? gi[(uint)(floor((randBase * _schema_i) / pi))] : gi[(uint)floor(  randBase * (gi.size() - _schema_i) ) + _schema_i];
}

#else //CLIQUE

void Graph::updateHasI(const node& v) {
	updateNeighborsBound(v, gs, myForeignIdx_s, schema_s);
	raiseSchema(v, schema_s, gs[v]);
}
void Graph::updateNoI(const node& v) {
	//Opposite to 'updateHasI()', here we must first modify v's neighbors' schema and only then update their respective bounds:
	lowerSchema(v, schema_s, gs[v]);
	updateNeighborsBound(v, gs, myForeignIdx_s, schema_s);
}
void Graph::updateHasS(const node& v) {
	updateNeighborsBound(v, gi, myForeignIdx_i, schema_i);
	raiseSchema(v, schema_i, gi[v]);
}
void Graph::updateNoS(const node& v) {
	//Opposite to 'updateHasS()', here we must first modify v's neighbors' schema and only then update their respective bounds:
	lowerSchema(v, schema_i, gi[v]);
	updateNeighborsBound(v, gi, myForeignIdx_i, schema_i);
}
void Graph::updateNeighborsBound(const node& v, vector<vector<node>>& g, vector<vector<uint>>& foreignIdx, const vector<uint>& schema) {
	//For each node w, swaps v and u in w's list of neighbors. Note that u is the node at w's schema bound. Nodes v and u must exchange their indexes accordingly, thus their respective 'myForeignIdx' lists are also updated.
	for (size_t idxW = 0; idxW < g[v].size(); ++idxW) {
		const node w = g[v][idxW];
		vector<node>& wNeighbors = g[w];
		const uint wSchema = schema[w];
		const node u = wNeighbors[wSchema];	// ----> Opposite to 'w', node 'u' cannot be simply referenced (&) here since its position in 'wNeighbors' is exchanged with 'v' (and the reference would then be to 'v' instead of 'u').
#ifdef AUTO_RELATION
		if (w == v && u != v) {
			const uint w_index_in_u = foreignIdx[w][wSchema];
			std::swap(wNeighbors[idxW], wNeighbors[wSchema]);
			foreignIdx[u][w_index_in_u] = (uint)idxW;
			std::swap(foreignIdx[w][idxW], foreignIdx[w][wSchema]);
		}
		else if (u != v && u != w) {
#else
		if (u != v) {
#endif
			const uint w_index_in_u = foreignIdx[w][wSchema];
			uint& u_index_in_w = foreignIdx[u][w_index_in_u];
			uint& v_index_in_w = foreignIdx[v][idxW];
#ifdef DEBUG
			std::string message = "The index of 'u' in the neighborhood list of 'w' is " + std::to_string(u_index_in_w) + " but should be equal to w's schema (which is " + std::to_string(wSchema) + " at the moment).";
			assertm(u_index_in_w == wSchema, message);
#endif
			std::swap(wNeighbors[v_index_in_w], wNeighbors[u_index_in_w]);
			std::swap(foreignIdx[w][v_index_in_w], foreignIdx[w][u_index_in_w]);	// ----> Node w (which has received the alert from v for an s-bound update) has to update its 'myForeignIdx_s' to reflect that the indexes at which v and u poll their own position at w's neighborhood were exchanged.
			std::swap(v_index_in_w, u_index_in_w);
		}
	}
}
void Graph::raiseSchema(const node& v, vector<uint>& schema, const vector<node>& neighbors) {
	const uint vNeighbSz = (uint)neighbors.size();
	for (uint idxW = 0; idxW < vNeighbSz; ++idxW)
		++(schema[neighbors[idxW]]);
}
void Graph::lowerSchema(const node& v, vector<uint>& schema, const vector<node>& neighbors) {
	const uint vNeighbSz = (uint)neighbors.size();
	for (uint idxW = 0; idxW < vNeighbSz; ++idxW)
		--(schema[neighbors[idxW]]);
}
const node& Graph::nextNodeForS(const node& _currNode, const real& p) {
	const uint& schema = schema_s[_currNode];
	const real& _ps = ps[_currNode][schema];
	return (p < _ps) ? gs[_currNode][(uint)(floor(p * schema))] : gs[_currNode][(uint)floor(  p * (gs[_currNode].size() - schema) ) + schema];
}
const node& Graph::nextNodeForI(const node& _currNode, const real& p) {
	const uint& schema = schema_i[_currNode];
	const real& _pi = pi[_currNode][schema];
	return (p < _pi) ? gi[_currNode][(uint)(floor((p * schema) / _pi))] : gi[_currNode][(uint)floor(p * (gi[_currNode].size() - schema)) + schema];
}
#endif //CLIQUE
#endif //PROTECTION_FX

#ifndef CLIQUE

void Graph::set2ndMoment() {
	block_prob.resize(largestDegree + 1, 0);		// ----> Each position refers to a node degree, hence this "+ 1" happening. Vectors in C++ are indexed from 0 to n-1 (where n is the size of the vector). If the largest degree is, say, 5, then we need to acess the position 'block_prob[5]' instead of 'block_prob[4]'. Note that block_prob[0] will always be 0 (since no 0-degree nodes exist in the LCC).
	q_b.resize(largestDegree + 1, 0);
	originalFreq.resize(largestDegree + 1, 0);
	kb.resize(largestDegree + 1, 0);
	rho_b.resize(largestDegree + 1, 0);

	//Frequencies:
	for (uint i = 0; i < lccSize; ++i) {
#ifdef PROTECTION_FX
		++(block_prob[gs[lcc[i]].size()]);		// ----> '-1' applied since every node was given an extra edge, which should not be counted here. This extra edge is an auto-relation, which allows an agent to remain at its current node upon a walk event (See the AUTORELATION macro definition and its associated commentary for more info). 
		++(originalFreq[gs[lcc[i]].size()-1]);		// ----> '-1' applied since every node was given an extra edge, which should not be counted here. This extra edge is an auto-relation, which allows an agent to remain at its current node upon a walk event (See the AUTORELATION macro definition and its associated commentary for more info). 
#else
		++(block_prob[g[lcc[i]].size()]);		
		++(originalFreq[g[lcc[i]].size() - 1]);		// ----> '-1' applied since every node was given an extra edge, which should not be counted here. This extra edge is an auto-relation, which allows an agent to remain at its current node upon a walk event (See the AUTORELATION macro definition and its associated commentary for more info). 
#endif
	}
	//Probabilities:
	for (uint b = 0; b < block_prob.size(); ++b) {	// ----> Equivalent to "p_b" in [1].
		block_prob[b] /= n;
		originalFreq[b] /= n;
	}
	for (uint b = 1; b < block_prob.size(); ++b) {	
		//q_b[b] = (b * block_prob[b]) / averageDegree;

		//TESTE!!!
		//q_b[b] = ((double)(b-1) * block_prob[b]) / originalAvDeg;	// ----> interessante para grandes endemias.
		//q_b[b] = ((double)(b-1) * block_prob[b]) / averageDegree;
		q_b[b] = ((double)(b - 1) * n * block_prob[b]) / (2.0 * (m - n * block_prob[b]));
	}
	//for (uint b = 1; b < block_prob.size(); ++b) {	
	//	kb[b] = sim::NUM_AGENTS * q_b[b];
	//}
	//sumKB = 0;
	//for (uint b = 1; b < block_prob.size(); ++b) {
	//	if (block_prob[b] > 0)
	//		sumKB += pow(kb[b],2) / block_prob[b];
	//}

	for (uint b = 1; b < block_prob.size(); ++b) {
		rho_b[b]	= 1.0 - pow(1.0 - ((double)b / (2 * m)), sim::NUM_AGENTS);
	}
	psi = 0;
	for (uint b = 1; b < block_prob.size(); ++b) {
		if (block_prob[b] == 0)
			continue;
		psi += 1.0 / (n * block_prob[b]);
		//psi += b * block_prob[b] * rho_b[b];
		//psi += b * originalFreq[b] * rho_b[b];
	}

	validBlocks = 0;
	for (uint b = 1; b < block_prob.size(); ++b) {
		if (block_prob[b] > 0) ++validBlocks;
	}
	//2nd moment:
	_2ndMmt = 0;
	original2ndMmt = 0;
	for (uint b = 1; b < block_prob.size(); ++b) {
		_2ndMmt += pow(b, 2) * block_prob[b];
		original2ndMmt += pow(b, 2) * originalFreq[b];
	}
	
}

void Graph::readGraph(const string& fileName, const size_t& totalNodes) {
#ifdef PROTECTION_FX
	gs.resize(totalNodes);
	gi.resize(totalNodes);
	myForeignIdx_s.resize(totalNodes);
	myForeignIdx_i.resize(totalNodes);
#else
	g.resize(totalNodes);
#endif //PROTECTION_FX
#ifdef GNP
	/*The following procedure to build a G(n,p) graph is based on the pseudocode from [1].

		[1] Efficient generation of large random networks
			Vladimir Batagelj and Ulrik Brandes
			Phys. Rev. E 71, 036113
	*/
	
	{// ----> Extra scope.
		std::random_device _rd;
		std::mt19937_64 _gen(_rd());										// ----> Generator.
		std::uniform_real_distribution<real> _U(0, 1);

		constexpr real _avDegree = 19.74;									// ----> The average degree <d> of a G(n,p) graph is (n-1)p. Here we first fix <d> at some value (say, the same <d> of a real network, for the sake of comparisons) and then adjust 'p' accordingly.
		constexpr uint _n = N;
		constexpr real _p = _avDegree / (_n - 1);
		const real log_q = log(1 - _p);

		int v = 1, w = -1;
		real r;
		while (v < _n) {
			r = _U(_gen);													// ----> Random number uniformly drawn from the interval (0,1).
			w = w + 1 + (int)std::floor(log(1 - r) / log_q);
			while (w >= v && v < _n) {
				w -= v;
				++v;
			}
			if (v < _n) {
#ifdef PROTECTION_FX
				gs[v].emplace_back(w);
				gs[w].emplace_back(v);
				gi[v].emplace_back(w);
				gi[w].emplace_back(v);
				myForeignIdx_s[v].emplace_back((uint)(gs[w].size() - 1));
				myForeignIdx_s[w].emplace_back((uint)(gs[v].size() - 1));
				myForeignIdx_i[v].emplace_back((uint)(gs[w].size() - 1));
				myForeignIdx_i[w].emplace_back((uint)(gs[v].size() - 1));
#else
				g[v].emplace_back(w);
				g[w].emplace_back(v);
#endif
				++m;
			}
		}
		n = _n;
	} //Extra scope.
#endif

#ifdef STAR
	{
#ifdef PROTECTION_FX
		for (size_t i = 1; i < N; ++i) {
			gs[0].emplace_back(i);
			gs[i].emplace_back(0);
			myForeignIdx_s[0].emplace_back((uint)(gs[i].size() - 1));
			myForeignIdx_s[i].emplace_back((uint)(gs[0].size() - 1));
		}
		
		for (size_t i = 1; i < N; ++i) {
			gi[0].emplace_back(i);
			gi[i].emplace_back(0);
			myForeignIdx_i[0].emplace_back((uint)(gs[i].size() - 1));
			myForeignIdx_i[i].emplace_back((uint)(gs[0].size() - 1));
		}
		
		//gs[v].emplace_back(w);
		//gs[w].emplace_back(v);
		//gi[v].emplace_back(w);
		//gi[w].emplace_back(v);

#else
		for (node i = 1; i < N; ++i) 
			g[0].emplace_back(i);
		for (node i = 1; i < N; ++i) 
			g[i].emplace_back(0);
#endif
		m = 2*(N-1) + 1;
		n = N;
	}
#endif


#ifdef READ_NTWK_FROM_FILE
	string serialFileName = string(EXE_DIR) + string("\\serializadas\\") + string(NWTK_LABEL) + string(".txt");
	string msg = "\nReading the network " + string(NWTK_LABEL) + "...";
	sim::Reporter::highlight(msg);

#ifdef DE_SERIALIZE
	ifstream serializedFile(serialFileName);
	if (serializedFile.is_open()) {
		Graph::startChronometer("De-serializing the network... ");
		uint numViz;
		serializedFile >> n;
		g.resize(n);
		for (uint i = 0; i < n; ++i) {
			serializedFile >> numViz;
			g[i].resize(numViz);
			for (uint j = 0; j < numViz; ++j)
				serializedFile >> g[i][j];
		}
		serializedFile >> m;
		serializedFile >> averageDegree;

		stopChronometer("Done");
	}
	else {
#endif
		n = 0;
		m = 0;
		selfLoops = 0;
		std::ifstream arq__G(fileName);
		string line;
		if (!arq__G.is_open()) {
			sim::Reporter::errorOpening(fileName);
			return;
		}
		sim::Reporter::openedWithSuccess(fileName);

		do {
			getline(arq__G, line);
		} while (line.at(0) == COMMENTARY && arq__G.good());

		sim::Reporter::startChronometer("Reading the network...");
		//The method below guarantees the network will be internally dealt with through sequential id's (0, 1, 2, ...)
		uint sequentialIDs = 0;
		uint idV1, idV2;
		std::stringstream ss(line);
		ss >> idV1;
		ss >> idV2;

		if (idMap.count(idV1) == 0) {
			++n;
			idMap[idV1] = sequentialIDs++;
		}
		if (idMap.count(idV2) == 0) {
			++n;
			idMap[idV2] = sequentialIDs++;
		}
		idV1 = idMap[idV1];
		idV2 = idMap[idV2];

		// Creates an edge connecting idV1 and idV2, if such an edge is not a self-loop:
		if (idV1 != idV2) {
#ifdef PROTECTION_FX
			gs[idV1].emplace_back(idV2);
			gs[idV2].emplace_back(idV1);
			gi[idV1].emplace_back(idV2);
			gi[idV2].emplace_back(idV1);
#else
			g[idV1].emplace_back(idV2);
			g[idV2].emplace_back(idV1);
#endif
			++m;
#ifdef PROTECTION_FX
			myForeignIdx_s[idV1].emplace_back((uint)(gs[idV2].size() - 1));
			myForeignIdx_s[idV2].emplace_back((uint)(gs[idV1].size() - 1));
			myForeignIdx_i[idV1].emplace_back((uint)(gs[idV2].size() - 1));
			myForeignIdx_i[idV2].emplace_back((uint)(gs[idV1].size() - 1));
#endif
		}

		while (arq__G.good()) {
			arq__G >> idV1;
			arq__G >> idV2;

			//IDs' translation:
			if (idMap.count(idV1) == 0) {
				++n;
				idMap[idV1] = sequentialIDs++;
			}
			if (idMap.count(idV2) == 0) {
				++n;
				idMap[idV2] = sequentialIDs++;
			}
			idV1 = idMap[idV1];
			idV2 = idMap[idV2];

			// Identifying (and ignoring) self-loops (edges connecting a node to itself):
			if (idV1 == idV2) {
				++selfLoops;
				continue;
			}

			// Creates an edge between idV1 and idV2 if it still does not exist.
#ifdef PROTECTION_FX
			if (find(gs[idV1].begin(), gs[idV1].end(), idV2) == gs[idV1].end()) {
				gs[idV1].emplace_back(idV2);
				gs[idV2].emplace_back(idV1);
				gi[idV1].emplace_back(idV2);
				gi[idV2].emplace_back(idV1);
				++m;

				myForeignIdx_s[idV1].emplace_back((uint)(gs[idV2].size() - 1));
				myForeignIdx_s[idV2].emplace_back((uint)(gs[idV1].size() - 1));
				myForeignIdx_i[idV1].emplace_back((uint)(gs[idV2].size() - 1));
				myForeignIdx_i[idV2].emplace_back((uint)(gs[idV1].size() - 1));	
			}
#else
			if (find(g[idV1].begin(), g[idV1].end(), idV2) == g[idV1].end()) {
				g[idV1].emplace_back(idV2);
				g[idV2].emplace_back(idV1);
				++m;
		}
#endif
	}

#ifdef DEBUG
		// * ASSERTS BOTH myForeignIdx_s AND myForeignIdx_i ARE BUILT CORRECTLY *
		{ // ----> Extra scope, so that inner variables are destroyed upon assertion.
			uint errors = 0;
			for (node v = 0; v < n; ++v) {
				for (uint pos = 0; pos < gs[v].size(); ++pos) {
					const node& u = gs[v][pos];
					const uint v_index_in_u = myForeignIdx_s[v][pos];
					if (v != gs[u][v_index_in_u]) { ++errors; }
				}
			}
			assertm(errors == 0, "bad 'myForeignIdx_s': not all indexes map to node's actual position on 'gs'");
			for (node v = 0; v < n; ++v) {
				for (uint pos = 0; pos < gi[v].size(); ++pos) {
					const node& u = gi[v][pos];
					const uint v_index_in_u = myForeignIdx_i[v][pos];
					if (v != gi[u][v_index_in_u]) { ++errors; }
				}
			}
			assertm(errors == 0, "bad 'myForeignIdx_i': not all indexes map to node's actual position on 'gi'");
		} 
#endif //DEBUG
		sim::Reporter::stopChronometer("done");
#endif //READ_NTWK_FROM_FILE
#ifdef PROTECTION_FX
		schema_s.resize(n);
		schema_i.resize(n);
#endif
		originalAvDeg = averageDegree = 2 * (float)m / n;
#ifdef SERIALIZE
		//Serializa a rede:
		startChronometer("Serializing the network (for future use)... ");
		ofstream arq(serialFileName);
		arq << g.size();
		for (size_t i = 0; i < g.size(); i++) {
			arq << " " << g[i].size();
			for (size_t j = 0; j < g[i].size(); j++)
				arq << " " << g[i][j];
		}
		arq << " " << m << " " << averageDegree;
		arq.close();
		stopChronometer("done");
#endif
#ifdef DE_SERIALIZE
	}
#endif
#ifdef PROTECTION_FX
	const vector<vector<node>>& g = gs;
#endif
	smallestDegree = (uint)g[0].size();
	largestDegree = (uint)g[0].size();
	largestDgNode = 0;
	for (node i = 0; i < n; ++i) {
		if (g[i].size() > largestDegree) {
			largestDegree = (uint)g[i].size();
			largestDgNode = i;
		}
		else {
			if (g[i].size() < smallestDegree) {
				smallestDegree = (uint)g[i].size();
			}
		}
	}

	//Determining the LCC:
	vector<node> Q(n);
	vector<bool> notVisited(n, true);
	lcc.resize(n);
	lccSize = 0;
	for (node v = 0; v < n; ++v) {
		if (notVisited[v]) {
			node next;
			uint front = 0, end = 0, ccSize = 0;
			Q[end++] = v;
			while (front != end) {
				next = Q[front];
				notVisited[next] = false;
				++ccSize;
				const uint neighbSz = (uint)(g[next].size());
				for (uint j = 0; j < neighbSz; ++j) {
					const uint& neighbor = g[next][j];
					if (notVisited[neighbor]) {
						Q[end++] = neighbor;
						notVisited[neighbor] = false;
					}
				}
				++front; // ---> "pop".
			}
			if (ccSize > lccSize) {
				lccSize = ccSize;
				for (uint i = 0; i < end; ++i) {
					lcc[i] = Q[i];
				}
			}
		}
	}
	lcc.resize(lccSize);

	sim::Reporter::networkInfo(n, m, averageDegree, largestDegree, smallestDegree, lccSize);

	//TESTE!!!
	//std::ofstream graphData;
	//std::stringstream tt;
	//tt.str("");			// ----> Clear content.
	//tt << EXE_DIR << "/stats/graphData_" << NWTK_LABEL << ".csv";
	//string dataFile = tt.str();
	//graphData.open(dataFile);
	//graphData.seekp(0, std::ios::end);
	//if (graphData.tellp() == 0)
	//	cout << "new";
	//else
	//	cout << "appending...";
	////Header
	//graphData << "node\td(v)\td(v)/2m\tm" << std::endl;
	////Data:
	//for (size_t i = 0; i < g.size(); ++i){
	//	graphData 
	//		<< i << "\t"
	//		<< g[i].size() + 1 << "\t"
	//		<< ((real)(g[i].size() + 1))/(real)((2*m)+n) << "\t"
	//		<< m
	//		<< std::endl;
	//}
	//graphData.close();

	if (selfLoops > 0)
		cout << "\t ---> " << selfLoops << " native self-loops identified and removed." <<  "\n";

#ifndef CLIQUE
#ifdef AUTO_RELATION
#ifdef PROTECTION_FX
	for (node v = 0; v < n; ++v) gs[v].emplace_back(v);
	for (node v = 0; v < n; ++v) gi[v].emplace_back(v);
	m += n;
	++largestDegree;
	averageDegree = (2.0 * m) / n;

	for (node v = 0; v < n; ++v) myForeignIdx_s[v].emplace_back((uint)(gs[v].size() - 1));
	for (node v = 0; v < n; ++v) myForeignIdx_i[v].emplace_back((uint)(gi[v].size() - 1));
#else
	for (node v = 0; v < n; ++v) g[v].emplace_back(v);
	m += n;
	averageDegree = (2.0 * m) / n;
	++largestDegree;
#endif //PROTECTION_FX
	cout << "\t ---> AUTORELATION active. Self-loop added for each node. Updated stats: \n";
	cout << "\t\t" << m << " nodes.\n";
	cout << "\t\tAverage degree = " << averageDegree << "\n";
	cout << "\t\tLargest degree = " << largestDegree << "\n\n";
#endif //AUTO_RELATION
#endif //CLIQUE
}
#endif // CLIQUE
