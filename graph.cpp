#include "graph.h"
#include "reporter.h"
using namespace graph;

std::unordered_map<uint, uint> Graph::idMap;
uint Graph::n;													// ----> Network size, i.e. the number of nodes.
//std::random_device graph::_rd;
//std::mt19937_64 graph::_generator(graph::_rd());				// ----> "mt" = "Mersenne Twister".
//std::uniform_real_distribution<real> graph::_distribution(0, 1);
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
vector<vector<node>> Graph::gs_safe;
vector<vector<node>> Graph::gs_unsafe;
vector<vector<node>> Graph::gi_safe;
vector<vector<node>> Graph::gi_unsafe;
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
float Graph::averageDegree;
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
			//if ((vDegree - schema + sumW) == 0) std::cout << "v=" << v << ", d=" << vDegree << ", schema=" << schema << ", sumW=" << sumW << '\n';
			vProbsS[schema] = sumW / (vDegree - schema + sumW);
			sumW += sim::Ws;
		}
		sumW = 0;
		for (uint schema = 0; schema < vProbsI.size(); ++schema) {
			//if ((vDegree - schema + sumW) == 0) std::cout << "v=" << v << ", d=" << vDegree << ", schema=" << schema << ", sumW=" << sumW << '\n';
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
void Graph::lowSchema	(uint& _schema) { --_schema; }
void Graph::updateHasI(const node& v) {
	updateSBound(v);
	raiseSchema(_schema_s);
}
void Graph::updateNoI(const node& v) {
	//Opposite to 'updateHasI()', here we must first modify v's neighbors' schema and only then update their respective bounds:
	lowSchema(_schema_s);
	updateSBound(v);
}
void Graph::updateHasS(const node& v) {
	updateIBound(v);
	raiseSchema(_schema_i);
}
void Graph::updateNoS(const node& v) {
	//Opposite to 'updateHasS()', here we must first modify v's neighbors' schema and only then update their respective bounds:
	lowSchema(_schema_i);
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
const node& Graph::nextNodeForS(const real& randBase) {
	//return (randBase < _ps[_schema_s]) ? g[sim::randomInt(_schema_s)] : g[sim::randomInt(n - _schema_s) + _schema_s];
	const real& ps = _ps[_schema_s];
	//return (randBase < ps) ? g[(uint)(floor((randBase * _schema_s) / ps))] : g[(uint)(floor((randBase * g.size()) / (1.0 - ps)))];
	return (randBase < ps) ? gs[(uint)(floor((randBase * _schema_s) / ps))] : gs[(uint)floor(  randBase * (gs.size() - _schema_s)  ) + _schema_s];
}
const node& Graph::nextNodeForI(const real& randBase) {
	//return (randBase < _pi[_schema_i]) ? g[sim::randomInt(_schema_i)] : g[sim::randomInt(n - _schema_i) + _schema_i];
	const real& pi = _pi[_schema_i];
	//return (randBase < pi) ? g[(uint)(floor((randBase * _schema_i) / pi))] : g[(uint)(floor((randBase * g.size()) / (1.0 - pi)))];
	return (randBase < pi) ? gi[(uint)(floor((randBase * _schema_i) / pi))] : gi[(uint)floor(  randBase * (gi.size() - _schema_i) ) + _schema_i];
}

#else //CLIQUE
//void Graph::resetProtectionFX() {
//	for (size_t i = 0; i < schema_s.size(); ++i) schema_s[i] = 0;
//	for (size_t i = 0; i < schema_i.size(); ++i) schema_i[i] = 0;
//}
void Graph::updateHasI(const node& v) {
	updateSBound(v);
	raiseSchema(v, schema_s, gs[v]);
}
void Graph::updateNoI(const node& v) {
	//Opposite to 'updateHasI()', here we must first modify v's neighbors' schema and only then update their respective bounds:
	lowSchema(v, schema_s, gs[v]);
	updateSBound(v);
}
void Graph::updateHasS(const node& v) {
	updateIBound(v);
	raiseSchema(v, schema_i, gi[v]);
}
void Graph::updateNoS(const node& v) {
	//Opposite to 'updateHasS()', here we must first modify v's neighbors' schema and only then update their respective bounds:
	lowSchema(v, schema_i, gi[v]);
	updateIBound(v);
}
void Graph::updateSBound(const node& v) {
	//For each node w, swaps v and u in w's list of neighbors. Note that u is the node at w's schema bound. Nodes v and u must exchange their indexes accordingly, thus their respective 'myForeignIdx' lists are also updated.
	const uint vNeighbSz = (uint)gs[v].size();
	for (uint idxW = 0; idxW < vNeighbSz; ++idxW) {
		//Convenient references (for both efficiency and readability):
		const node& w = gs[v][idxW];
		vector<node>& wNeighbors = gs[w];
		const uint& wSchema = schema_s[w];

		const node u = wNeighbors[wSchema];	// ----> Opposite to 'w', node 'u' cannot be simply referenced (&) here since its position in 'wNeighbors' is exchanged with 'v' (and the reference would then be to 'v' instead of 'u').
		if (u != v) {
			const uint w_index_in_u = myForeignIdx_s[w][wSchema];
			uint& u_index_in_w = myForeignIdx_s[u][w_index_in_u];
			uint& v_index_in_w = myForeignIdx_s[v][idxW];
#ifdef DEBUG
			std::string message = "The index of 'u' in the neighborhood list of 'w' is " + std::to_string(u_index_in_w) + " but should be equal to w's schema (which is " + std::to_string(wSchema) + " at the moment).";
			assertm(u_index_in_w == wSchema, message);
#endif
			std::swap(wNeighbors[v_index_in_w], wNeighbors[u_index_in_w]);
			std::swap(myForeignIdx_s[w][v_index_in_w], myForeignIdx_s[w][u_index_in_w]);	// ----> Node w (which has received the alert from v for an s-bound update) has to update its 'myForeignIdx_s' to reflect that the indexes at which v and u poll their own position at w's neighborhood were exchanged.
			std::swap(v_index_in_w, u_index_in_w);
		}
	}
}
void Graph::updateIBound(const node& v) {
	//For each node w, swaps v and u in w's list of neighbors. Note that u is the node at w's schema bound. Nodes v and u must exchange their indexes accordingly, thus their respective 'myForeignIdx' lists are also updated.
	const uint vNeighbSz = (uint)gi[v].size();
	for (uint idxW = 0; idxW < vNeighbSz; ++idxW) {
		//Convenient references (for both efficiency and readability):
		const node& w = gi[v][idxW];
		vector<node>& wNeighbors = gi[w];
		const uint& wSchema = schema_i[w];

		const node u = wNeighbors[wSchema];	// ----> Opposite to 'w', node 'u' cannot be simply referenced (&) here since its position in 'wNeighbors' is exchanged with 'v' (and the reference would then be to 'v' instead of 'u').
		if (u != v) {
			const uint w_index_in_u = myForeignIdx_i[w][wSchema];
			uint& u_index_in_w = myForeignIdx_i[u][w_index_in_u];
			uint& v_index_in_w = myForeignIdx_i[v][idxW];
#ifdef DEBUG
			std::string message = "The index of 'u' in the neighborhood list of 'w' is " + std::to_string(u_index_in_w) + " but should be equal to w's schema (which is " + std::to_string(wSchema) + " at the moment).";
			assertm(u_index_in_w == wSchema, message);
#endif
			std::swap(wNeighbors[v_index_in_w], wNeighbors[u_index_in_w]);
			std::swap(myForeignIdx_i[w][v_index_in_w], myForeignIdx_i[w][u_index_in_w]);	// ----> Node w (which has received the alert from v for an s-bound update) has to update its 'myForeignIdx_s' to reflect that the indexes at which v and u poll their own position at w's neighborhood were exchanged.
			std::swap(v_index_in_w, u_index_in_w);
		}
	}
}
void Graph::raiseSchema(const node& v, vector<uint>& schema, const vector<node>& neighbors) {
	const uint vNeighbSz = (uint)neighbors.size();
	for (uint idxW = 0; idxW < vNeighbSz; ++idxW)
		++(schema[neighbors[idxW]]);
}
void Graph::lowSchema(const node& v, vector<uint>& schema, const vector<node>& neighbors) {
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
void Graph::readGraph(const string& fileName, const size_t& totalNodes) {
	string serialFileName = string(EXE_DIR) + string("\\serializadas\\") + string(NWTK_LABEL) + string(".txt");

	bool buildNetwork = false; // ----> TRUE if the network still needs to be built; FALSE otherwise (serialized file).
	string msg = "\nReading the network " + string(NWTK_LABEL) + "...";
	sim::Reporter::highlight(msg);

#ifdef DE_SERIALIZE
	ifstream serializedFile(serialFileName);
	if (serializedFile.is_open()) {
		Graph::startChronometer("De-serializing the network... ");
		uint numViz;
		serializedFile >> n;
		g.resize(n);
		for (uint i = 0; i < n; i++) {
			serializedFile >> numViz;
			g[i].resize(numViz);
			for (uint j = 0; j < numViz; j++)
				serializedFile >> g[i][j];
		}
		serializedFile >> m;
		serializedFile >> averageDegree;

		stopChronometer("Done");
	}
	else {
#endif
#ifdef READ_NTWK_FROM_FILE
		buildNetwork = true;
		n = 0;
		m = 0;
		selfLoops = 0;
#ifdef PROTECTION_FX
		gs.resize(totalNodes);
		gi.resize(totalNodes);

		gs_safe.resize	(totalNodes); // An extra spot for each node (the last one) will be later added, in order to store the actual size of the list (0 initially), which will often differ from the container size.
		gs_unsafe.resize(totalNodes); // An extra spot for each node (the last one) will be later added, in order to store the actual size of the list (0 initially), which will often differ from the container size.
		gi_safe.resize	(totalNodes); // An extra spot for each node (the last one) will be later added, in order to store the actual size of the list (0 initially), which will often differ from the container size.
		gi_unsafe.resize(totalNodes); // An extra spot for each node (the last one) will be later added, in order to store the actual size of the list (0 initially), which will often differ from the container size.
#ifndef CLIQUE
		myForeignIdx_s.resize(totalNodes);
		myForeignIdx_i.resize(totalNodes);
#endif //CLIQUE
#else
		g.resize(totalNodes);
#endif //PROTECTION_FX
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
			
			gs_safe[idV1].emplace_back(idV2);
			//++(gs_safe[idV1][LIST_SIZE_POS]);
			gs_safe[idV2].emplace_back(idV1);
			//++(gs_safe[idV2][LIST_SIZE_POS]);
			gi_safe[idV1].emplace_back(idV2);
			//++(gi_safe[idV1][LIST_SIZE_POS]);
			gi_safe[idV2].emplace_back(idV1);
			//++(gi_safe[idV2][LIST_SIZE_POS]);
#else
			g[idV1].emplace_back(idV2);
			g[idV2].emplace_back(idV1);
#endif
			m++;
#ifdef PROTECTION_FX
#ifndef CLIQUE
			myForeignIdx_s[idV1].emplace_back((uint)(gs_safe[idV2].size() - 1));
			myForeignIdx_s[idV2].emplace_back((uint)(gs_safe[idV1].size() - 1));
			myForeignIdx_i[idV1].emplace_back((uint)(gi_safe[idV2].size() - 1));
			myForeignIdx_i[idV2].emplace_back((uint)(gi_safe[idV1].size() - 1));
#endif
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
				selfLoops++;
				continue;
			}

			// Creates an edge between idV1 and idV2 if it still does not exist.
#ifdef PROTECTION_FX
			if (find(gs_safe[idV1].begin(), gs_safe[idV1].end(), idV2) == gs_safe[idV1].end()) {
				gs[idV1].emplace_back(idV2);
				gs[idV2].emplace_back(idV1);
				gi[idV1].emplace_back(idV2);
				gi[idV2].emplace_back(idV1);

				gs_safe[idV1].emplace_back(idV2);
				//++(gs_safe[idV1][LIST_SIZE_POS]);
				gs_safe[idV2].emplace_back(idV1);
				//++(gs_safe[idV2][LIST_SIZE_POS]);
				gi_safe[idV1].emplace_back(idV2);
				//++(gi_safe[idV1][LIST_SIZE_POS]);
				gi_safe[idV2].emplace_back(idV1);
				//++(gi_safe[idV2][LIST_SIZE_POS]);
				m++;
#ifndef CLIQUE
				myForeignIdx_s[idV1].emplace_back((uint)(gs_safe[idV2].size() - 1));
				myForeignIdx_s[idV2].emplace_back((uint)(gs_safe[idV1].size() - 1));
				myForeignIdx_i[idV1].emplace_back((uint)(gi_safe[idV2].size() - 1));
				myForeignIdx_i[idV2].emplace_back((uint)(gi_safe[idV1].size() - 1));
#endif			
			}
#else
			if (find(g[idV1].begin(), g[idV1].end(), idV2) == g[idV1].end()) {
				g[idV1].emplace_back(idV2);
				g[idV2].emplace_back(idV1);
				m++;
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
#ifdef PROTECTION_FX
		schema_s.resize(n);
		schema_i.resize(n);

		gs_unsafe = gs_safe;
		gi_unsafe = gi_safe;
		for (size_t i = 0; i < totalNodes; ++i) gs_safe[i].emplace_back((uint)(gs_safe[i].size()));
		for (size_t i = 0; i < totalNodes; ++i) gi_safe[i].emplace_back((uint)(gs_safe[i].size()));
		for (size_t i = 0; i < totalNodes; ++i) gs_unsafe[i].emplace_back(0);
		for (size_t i = 0; i < totalNodes; ++i) gi_unsafe[i].emplace_back(0);
#endif
		sim::Reporter::stopChronometer("done");

		averageDegree = 2 * (float)m / n;
#endif
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
	const vector<vector<node>>& g = gs_safe;
#endif
	smallestDegree = (uint)g[0].size();
	largestDegree = (uint)g[0].size();
	largestDgNode = 0;
	for (node i = 0; i < n; i++) {
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

	if (buildNetwork && selfLoops > 0)
		cout << "\t ---> Removed self-loops: " << selfLoops << "\n\n";
}
#endif // CLIQUE
