#pragma once
#include "defs.h"
#include <algorithm>	// swap(int, int)

//#define SERIALIZE
//#define DE_SERIALIZE
#define COMMENTARY '#'

namespace graph {
	using std::vector; using std::string; using std::cout; using std::endl;
	
	
class Graph {
public:
	static std::unordered_map<uint, uint> idMap;
	static uint n;													// ----> Network size, i.e. the number of nodes (n = |V|).
#ifdef CLIQUE
#ifdef PROTECTION_FX
	static vector<node> gs;											// ----> The graph, with nodes varying their index according to the susceptible agents' random walk.
	static vector<node> gi;											// ----> The graph, with nodes varying their index according to the infected agents' random walk.
#else 
	static vector<node> g;											// ----> The graph.
#endif // PROTECTION_FX
#else  // CLIQUE
#ifdef PROTECTION_FX
	static vector<vector<node>> gs;									// ----> The graph, as an edge list, with nodes varying their index according to the susceptible agents' random walk.
	static vector<vector<node>> gi;									// ----> The graph, as an edge list, with nodes varying their index according to the infected agents' random walk.
#else 
	static vector<vector<node>> g;									// ----> The graph, as an edge list.
#endif // PROTECTION_FX
	static uint m;													// ----> Network's total number of edges.
	static uint largestDegree;
	static uint smallestDegree;
	static uint selfLoops;
	static uint largestDgNode;
	static uint lccSize;											// ---> The size (i.e. the number of nodes) of the Largest Connected Component (LCC).
	static vector<node>lcc;											// ---> List of nodes that belong to the LCC.
	static real averageDegree;
	static real originalAvDeg;
	static real original2ndMmt;
	static real _2ndMmt;
	static vector<double> frequency;								// ----> Degree blocks.
	static vector<double> originalFreq;								// ----> Degree blocks.
	static vector<double> k_b;										// ----> Expected number of agents in each block b.
	static vector<double> rho_b;									// ----> Probability that an specific node v_b from block b is NOT empty, i.e. the probability that v_b hosts at least one agent.
	static real sumKB;
	static real psi;

#ifdef PROTECTION_FX
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
	struct largerDegree {
		bool operator()(const node& v1, const node& v2) {
			return g[v1].size() > g[v2].size();
		}
	};
#endif
	static void readGraph(const string& fileName, const size_t& totalNodes);
	static void set2ndMoment();
#endif // CLIQUE
#ifdef PROTECTION_FX
private:
	enum class direction{raise, lower};
	//static const uint LIST_SIZE_POS = 0;
#ifdef CLIQUE
	static uint _schema_s;
	static uint _schema_i;
	static vector<real> _ps;
	static vector<real> _pi;
	static vector<uint> _myIdx_s;
	static vector<uint> _myIdx_i;
#else
	static vector<uint> schema_s;
	static vector<uint> schema_i;
	static vector<vector<real>> ps;
	static vector<vector<real>> pi;
	static vector<vector<uint>> myForeignIdx_s;
	static vector<vector<uint>> myForeignIdx_i;
#endif
public:
	static void resetSchema();
	static void resetAgentIdx();
	static void setProbs();
#ifdef CLIQUE
#ifdef PROPORTIONAL
	static const node& nextNodeForS(const real& randBase, const real& itotal, const real& stotal);
	static const node& nextNodeForI(const real& randBase, const real& itotal, const real& stotal);
#else
	static const node& nextNodeForS(const real& randBase);
	static const node& nextNodeForI(const real& randBase);
#endif //PROPORTIONAL
#else
	static const node& nextNodeForS(const node& _currNode, const real& p);
	static const node& nextNodeForI(const node& _currNode, const real& p);
#endif //CLIQUE

	//Updates node v's neighbors' schema, to reflect that v wasn't hosting any infected agent, and now one of such agents has just arrived at v. It means v is no longer a safe spot and a new schema is necessary to reflect that the probability at which nearby susceptible agents choose v as their next hop becomes smaller.
	static void updateHasI	(const node& v);
	//Updates node v's neighbors' schema, to reflect that v was hosting an infected agent, and now none of such agents is located at v. It means v is now a safe spot and a new schema is necessary to reflect that the probability at which nearby susceptible agents choose v as their next hop becomes larger.
	static void updateNoI	(const node& v);
	static void updateHasS	(const node& v);
	static void updateNoS	(const node& v);

private:
#ifdef CLIQUE
	static void updateSBound(const node& v);
	static void updateIBound(const node& v);
	static void raiseSchema (uint& _schema);
	static void lowerSchema	(uint& _schema);
#else
	static void raiseSchema	(const node& v, vector<uint>& schema, const vector<node>& neighbors);
	static void lowerSchema	(const node& v, vector<uint>& schema, const vector<node>& neighbors);
	static void updateNeighborsBound(const node&, vector<vector<node>>& g, vector<vector<uint>>& foreignIdx, const vector<uint>& schema);
#endif //CLIQUE
#endif //PROTECTION_FX
};
}