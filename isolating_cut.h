#include <vector>
#include <list>
#include <limits>
#include <map>
#include <set>

#include <lemon/list_graph.h>
#include <lemon/concepts/graph.h>
#include <lemon/bits/map_extender.h>
#include <lemon/bits/graph_extender.h>
#include <lemon/preflow.h>
#include <lemon/tolerance.h>
#include <lemon/color.h>
#include <lemon/core.h>

// write python pseudocode and then c++

namespace lemon {

#ifdef DOXYGEN
  template <typename GR, typename CAP, typename TOL>
#else
  template <typename GR,
            typename CAP = typename GR::template ArcMap<int>,
            typename TOL = Tolerance<typename CAP::Value> >
#endif

    class IsolatingCut {
        
    public:
        /// The graph type of the algorithm
        typedef GR ListGraph;
        /// The capacity map type of the algorithm
        typedef CAP CapacityMap;
        // /// The tolerance type of the algorithm
        // typedef TOL Tolerance;

    private:

        typedef typename CapacityMap::Value Value;

        TEMPLATE_GRAPH_TYPEDEFS(ListGraph);

        const Graph& _graph;
        const CapacityMap& _capacity;
        const CutMap& _cutmap;

        std::vector<Node>& _subset;
        std::map<Node, int32_t> _labels;
        std::map<int, std::vector<Node>& > _colorings:
        std::set<Edge> _color_mincuts;
        
    public:
    IsolatingCut(const Graph& graph, const CapacityMap& capacity, const vector<Node>& subset) :
        _graph(graph), _capacity(&capacity), _subset(subset), _color_mincuts() {}

    ~IsolatingCut() {
        if (_graph) {
            delete _graph;
        }
        if (_capacity) {
            delete _capacity;
        }
        if (_subset) {
            delete subset;
        }
        if (_labels) {
            delete _labels;
        }
        if (_colorings) {
            delete _colorings;
        }
        if (_color_minuts) {
            delete _color_mincuts;
        }
    }
 
    private:
        void initLabels(std::vector<Node>& _subset) {
            int len = _subset.size()
            for (int i = 0; i < len; ++i) {
                _labels[_subset[i]] = i
            }
        }
        ListGraph& copyGraph() {
            ListGraph graphCopy;
            GraphCopy<Graph, ListGraph> cg(_graph, graphCopy);
            Graph::NodeMap<ListGraph::Node> nr(_graph);
            cg.nodeRef(nr);
            Graph::EdgeMap<int> oemap(_graph);
            ListGraph::EdgeMap<int> nemap(graphCopy);
            cg.edgeMap(oemap, nemap);
            cg.run();
            return &graphCopy;
        }
        void colorContraction(ListGraph& graph) {
            for (int color : [0, 1]) {
                std::vector<int, Node>& component = _colorings[color];
                for (int i = 1; i < component.size(); ++i) {
                    graph.contract(component[0], component[i]);
                }
            }
        }
        void colorGraph(ListGraph& graph, int i) {
            for (auto& v : _subset) {
                _colorings[findColor(_labels[v], i)].push_back(v);
            }
        }
        int findColor(int label, int i) {
            return (label >> i) & 1;
        }

        void findMincuts(ListGraph& graph, Node& s, Node& t) {
            Preflow<ListGraph, CapacityMap> pf(graph, _capacity, *s, *t);
            pf.run();
            pf.minCutMap(_cutmap);
        }

    public:
        void init() {
            initLabels(_subset);
        }

        void run() {

        }

    }; //class IsolatingCut
} //namespace lemon

/*
Public 

what we need mincut for i colorings - use s-t maxflow w/ color contraction
when finding these mincuts, denote which v in v is in the C_i and which side it is on
find connected component with v
union of C_i
graph complement
maxflow call for disjoint union of G_v
recover S_v and w(S_v)


class IsolatingCut

    def __init__(subset: R):
        self.ordering = R
        self.labels = {}
        self.colorings = {}
        for i in range(len(R)):
            label[R[i]] = binary(i)
        
        for i in range(log(len(R))):
            for v in R:
                color(v, labels[v], i)
            mincut(C_i between red and blue in G) - use maxflow
    
    def _color(v, label: binary, i):
        self.colorings[r] = blue if binary[i] == 1 else red

    once C_i are found, U_v is subset in G \ U C_i that contains v
    for each v
    contract V \ U_v ito a single vertex t
    compute min v-t cut on contracted graph (equivalent to S_v)
    do these in parallel on all contracted graphs w/ max flow on disjoint union of contracted
    (there are |v| of them)

    recover the S_v and minimum weight across the S_v
    return that then 

*/