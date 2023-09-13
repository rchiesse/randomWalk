#pragma once
#include "defs.h"
#include <algorithm>	// swap(int, int)

//#define SERIALIZE
//#define DE_SERIALIZE
#define COMMENTARY '#'

namespace graph {
	using std::vector; using std::string; using std::cout; using std::endl;
	
	
class Graph {
private:
	static rtl numAgents;
	static rtl Ws;
	static rtl Wi;
public:
	static std::unordered_map<uint, uint> idMap;
	static uint n;													// ----> Network size, i.e. the number of nodes (n = |V|).
	static uint m;													// ----> Network's total number of edges.
	static uint largestDegree;
	static uint smallestDegree;
	static uint selfLoops;
	static uint largestDgNode;
	static uint lccSize;											// ---> The size (i.e. the number of nodes) of the Largest Connected Component (LCC).
	static rtl averageDegree;
	static rtl originalAvDeg;

#ifdef PROTECTION_FX
	static vector<vector<node>> gs;									// ----> The graph, as an edge list, with nodes varying their index according to the susceptible agents' random walk.
	static vector<vector<node>> gi;									// ----> The graph, as an edge list, with nodes varying their index according to the infected agents' random walk.
	struct largerDegreeGS {
		bool operator()(const node& v1, const node& v2) {
			return gs[v1].size() > gs[v2].size();
		}
	};
	struct largerDegreeGI {
		bool operator()(const node& v1, const node& v2) {
			return gi[v1].size() > gi[v2].size();
		}
	};
#else 
	static vector<vector<node>> g;									// ----> The graph, as an edge list.
	struct largerDegree {
		bool operator()(const node& v1, const node& v2) {
			return g[v1].size() > g[v2].size();
		}
	};
#endif // PROTECTION_FX
	static vector<node>lcc;											// ---> List of nodes that belong to the LCC.
	static rtl _2ndMmt;
	static rtl bk;													// ----> Expected block.
	static vector<rtl> q_b;											// ----> Probability q_b that a randomly chosen link points to a degree-b node.
	static vector<rtl> kb;											// ----> Expected number of agents in each block b.
	static vector<rtl> rho_b;										// ----> Probability that an specific node v_b from block b is NOT empty, i.e. the probability that v_b hosts at least one agent.
	static vector<rtl> Kbnb;
	static rtl maxKbnb;
	static rtl maxKbnbBlock;
	static vector<rtl> block_prob;
	static vector<uint> validBlocks;
	//static rtl avSelfLoop;
	static void setBlockData();
	static void readGraph(const string& fileName, const size_t& totalNodes);
	static void setParams(const uint& _numAgents, const rtl& _Ws, const rtl& _Wi);

#ifdef PROTECTION_FX
private:
	enum class direction{raise, lower};
	//static const uint LIST_SIZE_POS = 0;
	static vector<uint> schema_s;
	static vector<uint> schema_i;
	static vector<vector<rtl>> ps;
	static vector<vector<rtl>> pi;
	static vector<vector<uint>> myForeignIdx_s;
	static vector<vector<uint>> myForeignIdx_i;
public:
	static void resetSchema();
	static void setProbs();
	static const node& nextNodeForS(const node& _currNode, const rtl& p);
	static const node& nextNodeForI(const node& _currNode, const rtl& p);

	//Updates node v's neighbors' schema, to reflect that v wasn't hosting any infected agent, and now one of such agents has just arrived at v. It means v is no longer a safe spot and a new schema is necessary to reflect that the probability at which nearby susceptible agents choose v as their next hop becomes smaller.
	static void updateHasI	(const node& v);
	//Updates node v's neighbors' schema, to reflect that v was hosting an infected agent, and now none of such agents is located at v. It means v is now a safe spot and a new schema is necessary to reflect that the probability at which nearby susceptible agents choose v as their next hop becomes larger.
	static void updateNoI	(const node& v);
	static void updateHasS	(const node& v);
	static void updateNoS	(const node& v);

private:
	static void raiseSchema	(const node& v, vector<uint>& schema, const vector<node>& neighbors);
	static void lowerSchema	(const node& v, vector<uint>& schema, const vector<node>& neighbors);
	static void updateNeighborsBound(const node&, vector<vector<node>>& g, vector<vector<uint>>& foreignIdx, const vector<uint>& schema);
#endif //PROTECTION_FX
};
}