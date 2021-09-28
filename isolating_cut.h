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
            typename CAP = typename GR::template EdgeMap<int>,
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

        const ListGraph& _graph;
        const CapacityMap& _capacity;
        typedef typename ListGraph::template NodeMap<bool> MinCutMap;
        const MinCutMap* _cutmap;

        const std::vector<Node>& _subset;
        std::map<Node, int32_t> _labels;
        std::map<int, std::vector<Node>& > _colorings;
        // we need some struct of saving the logV colorings
        // typedef typename ListGraph::template EdgeMap<Color> ColorMap;
        // ColorMap _color_mincuts;
        
    public:
    IsolatingCut(const ListGraph& graph, const CapacityMap& capacity, const vector<Node>& subset) :
        _graph(graph), _capacity(capacity), _subset(subset) {
            init();
        }

    ~IsolatingCut() {
        if (_graph) {
            delete _graph;
        }
        if (_capacity) {
            delete _capacity;
        }
        if (_subset) {
            delete _subset;
        }
        if (_labels) {
            delete _labels;
        }
        if (_colorings) {
            delete _colorings;
        }
    }
 
    private:
        void initLabels(std::vector<Node>& _subset) {
            int len = _subset.size();
            for (int i = 0; i < len; ++i) {
                _labels[_subset[i]] = i;
            }
        }
        void copyGraph(ListGraph& new_graph, ListGraph& orig_graph) {
            GraphCopy<ListGraph, ListGraph> cg(orig_graph, new_graph);
            cg.run();
        }
        std::list<Node> colorContraction(ListGraph& graph) {
            std::list<Node> st;
            std::vector<int> colors = {0, 1};
            for (const auto& color : colors) {
                std::vector<int, Node>& component = _colorings[color];
                st.push_back(component[0]);
                for (int i = 1; i < component.size(); ++i) {
                    graph.contract(component[0], component[i]);
                }
            }
            return st;
        }
        void colorGraph(ListGraph& graph, int i) {
            for (const auto& v : _subset) {
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

        void findIthMinCut(ListGraph& g, int i) {
            colorGraph(g, i);
            std::list<Node> st = colorContraction(g);
            findMincuts(g, st[0], st[1]);
        }

    public:
        void init() {
            initLabels(_subset);
        }

        void run() {
            int len = _subset.size();
            for (int i = 0; i < len; ++i) {
                ListGraph newGraph;
                copyGraph(newGraph, _graph);
                findIthMinCut(newGraph, i);
            }
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