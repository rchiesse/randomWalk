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
	vector<real>& probsS = _ps;	// ----> It is important to create the reference only AFTER _ps BEING RESIZED. The need for contiguous space may lead '_ps' to be reallocated in memory upon resizing, what in turn would invalidate any previously-created reference (as it would still be pointing to _ps's old address).
	vector<real>& probsI = _pi;	// ----> It is important to create the reference only AFTER _pi BEING RESIZED. The need for contiguous space may lead '_pi' to be reallocated in memory upon resizing, what in turn would invalidate any previously-created reference (as it would still be pointing to _pi's old address).
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

#ifdef CLIQUE
void Graph::resetSchema() {
	_schema_s = 0;
	_schema_i = 0;
}
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
	const uint vNeighbSz = (uint)gs[v].size();
	for (uint idxW = 0; idxW < vNeighbSz; ++idxW) {
		//Not necessary but quite convenient REFERENCES (for both efficiency and readability):
		const node& w = gs[v][idxW];
		vector<node>& wNeighbors = gs[w];
		uint& wSchema = schema_s[w];
		uint& v_index_in_w = myForeignIdx_s[v][w];
		const uint& w_index_in_u = myForeignIdx_s[w][wSchema];

		//Swap v and u. Note that u is the node at w's schema bound, to be swapped with v on w's neighborhood list. Also, v and u exchange their indexes and update their respective 'myForeignIdx' lists accordingly.
		const node u = wNeighbors[wSchema];
		wNeighbors[wSchema] = v;
		wNeighbors[v_index_in_w] = u;
		myForeignIdx_s[u][w_index_in_u] = v_index_in_w;
		v_index_in_w = wSchema;
	}
}
void Graph::updateIBound(const node& v) {
	const uint vNeighbSz = (uint)gi[v].size();
	for (uint idxW = 0; idxW < vNeighbSz; ++idxW) {
		//Not necessary but quite convenient REFERENCES (for both efficiency and readability):
		const node& w = gi[v][idxW];
		vector<node>& wNeighbors = gi[w];
		uint& wSchema = schema_i[w];
		uint& v_index_in_w = myForeignIdx_i[v][w];
		const uint& w_index_in_u = myForeignIdx_i[w][wSchema];

		//Swap v and u. Note that u is the node at w's schema bound, to be swapped with v on w's neighborhood list. Also, v and u exchange their indexes and update their respective 'myForeignIdx' lists accordingly.
		const node u = wNeighbors[wSchema];
		wNeighbors[wSchema] = v;
		wNeighbors[v_index_in_w] = u;
		myForeignIdx_i[u][w_index_in_u] = v_index_in_w;
		v_index_in_w = wSchema;
	}
}
void Graph::raiseSchema(const node& v, vector<uint>& schema, const vector<node>& neighbors) {
	const uint vNeighbSz = (uint)neighbors.size();
	for (uint idxW = 0; idxW < vNeighbSz; ++idxW)
		++schema[neighbors[idxW]];
}
void Graph::lowSchema(const node& v, vector<uint>& schema, const vector<node>& neighbors) {
	const uint vNeighbSz = (uint)neighbors.size();
	for (uint idxW = 0; idxW < vNeighbSz; ++idxW)
		--schema[neighbors[idxW]];
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
void Graph::readGraph(const string& fileName, const uint& totalNodes) {
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
			n++;
			idMap[idV1] = sequentialIDs++;
		}
		if (idMap.count(idV2) == 0) {
			n++;
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
			m++;
#ifdef PROTECTION_FX
#ifndef CLIQUE
			myForeignIdx_s[idV1][gs[idV1].size() - 1] = gs[idV2].size() - 1;
			myForeignIdx_i[idV1][gi[idV1].size() - 1] = gi[idV2].size() - 1;
#endif
#endif
		}

		while (arq__G.good()) {
			arq__G >> idV1;
			arq__G >> idV2;

			//IDs' translation:
			if (idMap.count(idV1) == 0) {
				n++;
				idMap[idV1] = sequentialIDs++;
			}
			if (idMap.count(idV2) == 0) {
				n++;
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
			if (find(gs[idV1].begin(), gs[idV1].end(), idV2) == gs[idV1].end()) {
				gs[idV1].emplace_back(idV2);
				gs[idV2].emplace_back(idV1);
				gi[idV1].emplace_back(idV2);
				gi[idV2].emplace_back(idV1);
				m++;
#ifndef CLIQUE
				myForeignIdx_s[idV1][gs[idV1].size() - 1] = gs[idV2].size() - 1;
				myForeignIdx_i[idV1][gi[idV1].size() - 1] = gi[idV2].size() - 1;
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
	const vector<vector<node>>& g = gs;
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
		cout << "\t ---> Removed self-loops: " << selfLoops << endl << endl;
}
#endif // CLIQUE
