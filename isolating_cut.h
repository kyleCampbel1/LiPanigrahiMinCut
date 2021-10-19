#include <vector>
#include <list>
#include <limits>
#include <map>

#include <math.h>

#include <lemon/list_graph.h>
#include <lemon/concepts/graph.h>
#include <lemon/concepts/graph_components.h>
#include <lemon/bits/graph_extender.h>
#include <lemon/concepts/bpgraph.h>

#include <lemon/preflow.h>
#include <lemon/dfs.h>
#include <lemon/core.h>

// TODO move func defs to .cc file

namespace lemon {
    template <typename GR, typename CAP>
    struct LiPanigrahiDefaultTraits {

    typedef GR ListGraph;

    typedef CAP CapacityMap;

#ifdef DOXYGEN
    typedef GR::NodeMap<Value> MinCutMap;
#else
    typedef typename ListGraph::template NodeMap<bool> MinCutMap;
#endif

    static MinCutMap* createCutMap(const ListGraph& graph) {
      return new MinCutMap(graph);
    }


#ifdef DOXYGEN
    typedef GR::EdgeMap<Value> MinCutEdgeMap;
#else
    typedef typename ListGraph::template EdgeMap<bool> MinCutEdgeMap;
#endif

    static MinCutEdgeMap* createCutEdgeMap(const ListGraph& graph) {
      return new MinCutEdgeMap(graph);
    }

typedef typename std::vector<int> NodeIds;

typedef struct SV {
    int weight;
    NodeIds sv;
} SV;

};



#ifdef DOXYGEN
  template <typename GR, typename CAP, typename TOL>
#else
  template <typename GR,
            typename CAP = typename GR::template EdgeMap<int>,
            typename TR = LiPanigrahiDefaultTraits<GR, CAP> >
#endif

    class IsolatingCut {
        
    public:
        typedef TR Traits;
        /// The graph type of the algorithm
        typedef GR ListGraph;
        /// The capacity map type of the algorithm
        typedef CAP CapacityMap;

        typedef typename Traits::SV SV;

        typedef typename Traits::MinCutMap MinCutMap;

        typedef typename Traits::MinCutEdgeMap MinCutEdgeMap;

        typedef typename Traits::NodeIds NodeIds;
       
    private:

        typedef typename CapacityMap::Value Value;

        TEMPLATE_GRAPH_TYPEDEFS(ListGraph);

    public:
    IsolatingCut(const ListGraph& graph, CapacityMap& capacity, vector<Node>& subset, ListGraph& c1, ListGraph& c2):
        _capacity(&capacity), _subset(&subset), _ogGraph(&copyGraph(c1, graph)), _graph(&copyGraph(c2, graph)),
        _s((*_graph).addNode()), _t((*_graph).addNode()), _maxWeight(graphSumWeight(capacity))
    {}

    ~IsolatingCut() {
        delete _allCutsUnion;
        delete _cutmap;
        delete isolatingCuts;
    }
    
    private: 
        // graph pointer to original graph
        ListGraph* _ogGraph;
        // graph we will find s-t cuts and create edges on
        ListGraph* _graph;
        CapacityMap* _capacity;
        Node _s, _t;
        int _maxWeight = 0;

        
        MinCutEdgeMap* _allCutsUnion;
        MinCutMap* _cutmap; 

        std::vector<SV>* isolatingCuts; // eventually make this thread safe

        std::vector<Node>* _subset;
 
    private:
        void createStructures() {

            // structures for storing results from log|V| maxflow comps
            _allCutsUnion = TR::createCutEdgeMap(*_ogGraph);
            _cutmap = TR::createCutMap(*_ogGraph);
            isolatingCuts = new std::vector<SV>;
        }

        void removeCuts(ListGraph& trimGraph) {
        
            for (EdgeIt e(*_ogGraph); e != INVALID; ++e) {
                if ((*_allCutsUnion)[e] == true) {
                    trimGraph.erase(e);
                }
            }
        }

        ListGraph& copyGraph(ListGraph& new_graph, const ListGraph& orig_graph) {
            ExtendedListGraphBase::NodeMap<ListGraphBase::Node> nr(orig_graph);
            ExtendedListGraphBase::EdgeMap<ListGraphBase::Edge> ecr(new_graph);
            GraphCopy<ListGraph, ListGraph> gc(orig_graph, new_graph);
            gc.nodeRef(nr).edgeCrossRef(ecr);
            gc.run();
            return new_graph;
        }

        bool findColor(int label, int i) {
            bool val = label & (1 << i);
            return val;
        }

        int graphSumWeight(CAP& capacityMap) {
            int weight = 0;
            for (EdgeIt e(*_graph); e != INVALID; ++e) {
                weight += capacityMap[e];
            }
            return weight;
        }

        int findMincuts(ListGraph& graph, Node s, Node t, MinCutMap& cutmap) {
            Preflow<ListGraph, CapacityMap> pf(graph, *_capacity, s, t);
            pf.run();
            pf.minCutMap(cutmap);
            return pf.flowValue();
        }

        // Uv is a mapping of the connected component of v
        // contractGraph is a copy of the G/ union(mincuts)
        Node contractUv(MinCutMap& Uv, ListGraph& contractGraph) {
            Node firstNode;
            bool fnSet = false;

            // we can contract contractGraph off of og _graph bc it is a copy
            for (NodeIt n(*_ogGraph); n != INVALID; ++n) {
                if (! Uv[n]) {
                    if (! fnSet) {
                        fnSet = true;
                        firstNode = n;
                    }
                    else {
                        contractGraph.contract(firstNode, n);
                    }
                }
            } 
            return firstNode;
        }

        void findUv(const Node& v, MinCutMap& compMap, const ListGraph& isolatedCuts) {
            Dfs<ListGraph> dfs(isolatedCuts);
            dfs.reachedMap(compMap);
            dfs.run(v);
        }

        void getCutNodeIds(std::vector<int>& nodeIds, MinCutMap& mincut, ListGraph& graph) {
            for (NodeIt n(*_ogGraph); n != INVALID; ++n) {
                if (mincut[n] == true) {
                    nodeIds.push_back(graph.id(n));
                }
            }
        }

        // once the isolated components are found, this finds isolating cuts
        void parallelContractMinCut(ListGraph& isolatedComponents) {
            // in parallel loop -- TODO
            MinCutMap* connectedComponent = TR::createCutMap(*_ogGraph);
            for (const auto& v : *_subset) {
                
                findUv(v, *connectedComponent, isolatedComponents);

                ListGraph cgraph;
                copyGraph(cgraph, *_ogGraph);

                Node s = contractUv(*connectedComponent, cgraph);
                int dSV = findMincuts(cgraph, s, v, *connectedComponent);
                std::vector<int> nIds;
                getCutNodeIds(nIds, *connectedComponent, *_ogGraph);
                (*isolatingCuts).push_back(SV{dSV, nIds});
                
            }
            delete connectedComponent;
        }

        void undoSTComponents() {
            for (OutArcIt a(*_graph, _s); a!=INVALID; ++a) {
                (*_graph).erase(a);
            };
            for (InArcIt a(*_graph, _t); a!=INVALID; ++a) {
                (*_graph).erase(a);
            }
        }

        void addToCutUnion() {
            for (typename ListGraph::ArcIt e(*_ogGraph); e != INVALID; ++e) {
                if ((*_cutmap)[(*_graph).source(e)] != (*_cutmap)[(*_graph).target(e)]) {
                    (*_allCutsUnion)[e] = true;
                } 
            }
        }

        void makeSTComponents(int i) {
            int len = (*_subset).size();
            for (int j = 0; j < len; ++j) {
                Node tmp = (*_subset)[j];
                if (findColor(j, i)) {
                    (*_capacity)[(*_graph).addEdge(_s, tmp)] = _maxWeight;
                }
                else {
                    (*_capacity)[(*_graph).addEdge(_t, tmp)] = _maxWeight;
                }
            }
        }
        
        // removes edges for each of the log|V| calls
        void findIthMinCut(int i) {
            makeSTComponents(i);
            int _ = findMincuts(*_graph, _s, _t, *_cutmap);
            addToCutUnion();
            undoSTComponents();
        }

    public:

        void init() {
            createStructures();
        }

        // must call IsolatingCut().init() before .run()
        void run() {
            int len = ceil(log2((*_subset).size()));
            for (int i = 0; i < len; ++i) {
                findIthMinCut(i);
            }
            ListGraph trimGraph;
            copyGraph(trimGraph, *_graph);
            removeCuts(trimGraph);
            parallelContractMinCut(trimGraph);
        }

        std::vector<SV> getIsolatingCuts() {
            return *isolatingCuts;
        }

        // test methods
        int testIthMinCut(int i) {
            makeSTComponents(i);
            int cut = findMincuts(*_graph, _s, _t, *_cutmap);
            return cut;
        }

        int getMaxWeight() {
            return _maxWeight;
        }

        void testPhase1(Node& v, MinCutMap& compMap, Node& t) {
            int len = ceil(log2((*_subset).size()));
            for (int i = 0; i < len; ++i) {
                findIthMinCut(i);
            }
            ListGraph trimGraph;
            copyGraph(trimGraph, *_ogGraph);
            removeCuts(trimGraph);
            ExtendedListGraphBase::NodeMap<ListGraphBase::Arc> pm(trimGraph);

            Dfs<ListGraph> dfs(trimGraph);
            dfs.reachedMap(compMap);
            dfs.run(v);
            PredMapPath<ListGraph, ExtendedListGraphBase::NodeMap<ExtendedListGraphBase::Arc>> a = dfs.path(t);
        }

    }; //class IsolatingCut
} //namespace lemon


