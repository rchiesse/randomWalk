#pragma once
#include "defs.h"
#include <algorithm>	// swap(int, int)

//#define SERIALIZE
//#define DE_SERIALIZE
#define COMMENTARY '#'

#define NWTK_LABEL "Clique"
#define SHORT_LABEL "CL"
#ifdef GNP
#define NTWK_SIZE 23133
#define SOURCE_FILE "", NTWK_SIZE
#define NWTK_LABEL "Gnp"
#define SHORT_LABEL "gnp"
#endif
//#define NTWK_SIZE 55
//#define SOURCE_FILE std::string(std::string(EXE_DIR) + "\\redes\\grafoDeTestes.txt"), 55
//#define NWTK_LABEL "Ronald"
//#define SHORT_LABEL "ron"

//#define NTWK_SIZE 4039
//#define SOURCE_FILE std::string(std::string(EXE_DIR) + "\\redes\\facebook_combined.txt"), 4039
//#define NWTK_LABEL "Fb"
//#define SHORT_LABEL "fb"

//#define NTWK_SIZE 12008
//#define SOURCE_FILE std::string(std::string(EXE_DIR) + "\\redes\\CA-HepPh.txt"), 12008
//#define NWTK_LABEL "HepPh"
//#define SHORT_LABEL "hep"

//#define NTWK_SIZE 15233
//#define SOURCE_FILE std::string(std::string(EXE_DIR) + "\\redes\\netHEPT.txt"), 15233
//#define NWTK_LABEL "net"
//#define SHORT_LABEL "net"

//#define NTWK_SIZE 18772
//#define SOURCE_FILE std::string(std::string(EXE_DIR) + "\\redes\\CA-AstroPh.txt"), 18772
//#define NWTK_LABEL "AstroPh"
//#define SHORT_LABEL "astro"

//#define NTWK_SIZE 23133
//#define SOURCE_FILE std::string(std::string(EXE_DIR) + "\\redes\\CA-CondMat.txt"), 23133
//#define NWTK_LABEL "CondMat"
//#define SHORT_LABEL "cmat"

//#define NTWK_SIZE 36692
//#define SOURCE_FILE std::string(std::string(EXE_DIR) + "\\redes\\Email-Enron.txt"), 36692
//#define NWTK_LABEL "Enron"
//#define SHORT_LABEL "enron"

//#define NTWK_SIZE 58228
//#define SOURCE_FILE std::string(std::string(EXE_DIR) + "\\redes\\Brightkite_edges.txt"), 58228
//#define NWTK_LABEL "Brightkite"
//#define SHORT_LABEL "bk"

//#define NTWK_SIZE 196591
//#define SOURCE_FILE std::string(std::string(EXE_DIR) + "\\redes\\Gowalla_edges.txt"), 196591
//#define NWTK_LABEL "Gowalla"
//#define SHORT_LABEL "gw"

//#define NTWK_SIZE 317080
//#define SOURCE_FILE std::string(std::string(EXE_DIR) + "\\redes\\com-dblp.ungraph.txt"), 317080
//#define NWTK_LABEL "dblp"
//#define SHORT_LABEL "dblp"

//#define NTWK_SIZE 334863
//#define SOURCE_FILE std::string(std::string(EXE_DIR) + "\\redes\\com-amazon.ungraph.txt"), 334863
//#define NWTK_LABEL "amazon"
//#define SHORT_LABEL "amz"

//#define NTWK_SIZE 1696415
//#define SOURCE_FILE std::string(std::string(EXE_DIR) + "\\redes\\as-skitter.txt"), 1696415
//#define NWTK_LABEL "as-skitter"
//#define SHORT_LABEL "as"

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
	static float averageDegree;
	
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