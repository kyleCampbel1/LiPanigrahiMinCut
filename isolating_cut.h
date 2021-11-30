#include <vector>
#include <list>
#include <limits>
#include <map>
#include <stdexcept>

#include <math.h>

#include <lemon/list_graph.h>
#include <lemon/concepts/graph.h>
#include <lemon/concepts/graph_components.h>
#include <lemon/bits/graph_extender.h>
#include <lemon/concepts/bpgraph.h>

#include <lemon/preflow.h>
#include <lemon/dfs.h>
#include <lemon/core.h>


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
};



#ifdef DOXYGEN
  template <typename GR, typename CAP, typename TOL>
#else
  template <typename GR,
            typename TR = LiPanigrahiDefaultTraits<GR, typename GR::template EdgeMap<int>> >
#endif

    class IsolatingCut {
        
    public:
        typedef TR Traits;
        /// The graph type of the algorithm
        // typedef GR ListGraph;
        /// The capacity map type of the algorithm
        typedef ListGraph::EdgeMap<int> CapacityMap;

        typedef typename Traits::MinCutMap MinCutMap;

        typedef typename Traits::MinCutEdgeMap MinCutEdgeMap;
       
    private:

        typedef typename CapacityMap::Value Value;

        static CapacityMap* createCapactiyMap(GR& graph) {
            return new CapacityMap(graph);
        }

        TEMPLATE_GRAPH_TYPEDEFS(ListGraph);

    public:
    // @graph the ListGraph to compute Mininimum cut of
    // @capacity the EdgeMap for weights of edges in graph
    // @subset the list of nodes to compute the minimum isolating cut between
    IsolatingCut(const ListGraph& graph, CapacityMap& capacity, vector<Node>& subset):
        _capacity(&capacity), _subset(&subset), _ogGraph(&graph), _maxWeight(graphSumWeight(capacity))
    {}

    ~IsolatingCut() {
        delete _allCutsUnion;
        delete _capMap;
        delete _graph;

    }
    
    private: 
        // graph pointer to original graph
        const ListGraph* _ogGraph;
        // graph we will find s-t cuts and create edges on
        ListGraph* _graph;
        CapacityMap* _capacity;
        Node _s, _t;
        int _maxWeight;
        int _numEdges;

        MinCutEdgeMap* _allCutsUnion;
        CapacityMap* _capMap;

       int minimalCutValue = INT_MAX;


        std::vector<Node>* _subset;

 
    private:
        void createStructures() {

            // structures for storing results from log|V| maxflow comps
            _numEdges = (*_ogGraph).maxEdgeId();
            _graph = new ListGraph;
            _capMap = createCapactiyMap(*_graph);
            copyGraph(*_graph, *_capMap, *_ogGraph, *_capacity);
            _allCutsUnion = TR::createCutEdgeMap(*_graph);
            _s = (*_graph).addNode();
            _t = (*_graph).addNode();
        }

        // removes the edges included in U C_i
        // induces |V| connected components U_v for each v in subset
        void removeCuts(ListGraph& trimGraph) {
            for (EdgeIt e(trimGraph); e != INVALID; ++e) {
                if ((*_allCutsUnion)[e] == true) {
                    trimGraph.erase(e);
                }
            }
        }

        void copyGraph(ListGraph& newGraph, CapacityMap& newMap, const ListGraph& oldGraph, CapacityMap& oldMap) {
            for (int i = 0; i <= oldGraph.maxNodeId(); ++i) {
                newGraph.addNode();
            }
            // ExtendedListGraphBase::NodeMap<ListGraphBase::Node> nr(oldGraph);
            // ExtendedListGraphBase::EdgeMap<ListGraphBase::Edge> ecr(newGraph);
            // GraphCopy<ListGraph, ListGraph> gc(oldGraph, newGraph);
            // gc.nodeRef(nr).edgeCrossRef(ecr);
            // gc.run();
            for (ListGraph::EdgeIt e(oldGraph); e != INVALID; ++e) {
                newMap[newGraph.addEdge(oldGraph.u(e), oldGraph.v(e))] = oldMap[e];
            }
        }

        // performs vertex labeling given vertex index and i of iteration
        bool findColor(int label, int i) {
            bool val = label & (1 << i);
            return val;
        }

        int graphSumWeight(CapacityMap& capacityMap) {
            int weight = 0;
            for (EdgeIt e(*_ogGraph); e != INVALID; ++e) {
                weight += capacityMap[e];
            }
            return weight;
        }

        // this finds a min s-t and returns value of the cut. sets cutmap to bools for nodes in/out of cut
        int findMincuts(ListGraph& graph, Node s, Node t, MinCutMap& cutmap, CapacityMap& cap) {
            Preflow<ListGraph, CapacityMap> pf(graph, cap, s, t);
            pf.run();
            pf.minCutMap(cutmap);
            return pf.flowValue();
        }

        // Uv is a mapping of the connected component of v
        // contractGraph is a copy of the G/ union(mincuts)
        // performs contraction on graph that is not reachable from v to speed up maxflow computation
        Node contractUv(MinCutMap& Uv, ListGraph& contractGraph) {
            Node firstNode;
            bool fnSet = false;

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
            if (fnSet) {
                return firstNode;
            }
            throw runtime_error("Contraction failed");
        }

        // performs dfs to set compMap to the nodes reachable from v in G\ U C_i
        void findUv(const Node& v, MinCutMap& compMap, const ListGraph& isolatedCuts) {
            Dfs<ListGraph> dfs(isolatedCuts); // could pass this in as arg because it gets called often
            dfs.reachedMap(compMap);
            dfs.run(v);
        }

        // retrieves the nodeIds from a NodeMap to store cut information
        void getCutNodeIds(vector<int>& nodeIds, MinCutMap& mincut, const ListGraph& graph) {
            for (NodeIt n(*_ogGraph); n != INVALID; ++n) {
                if (mincut[n] == true) {
                    nodeIds.push_back(graph.id(n));
                }
            }
        }

        // G\ U C_i is found, this finds isolating cuts
        void parallelContractMinCut(ListGraph& isolatedComponents) {
            // in parallel loop -- TODO
            // check runtime of this loop
            for (const auto& v : *_subset){
                MinCutMap* connectedComponent = TR::createCutMap(*_ogGraph);  // O(N)
                // finds connected component that contains v in graph
                // that is pruned to exclude the isolating cuts
                findUv(v, *connectedComponent, isolatedComponents);  // O(|Uv(n+e)|)

                ListGraph cgraph;
                ListGraph::EdgeMap<int> cmap(cgraph);
                copyGraph(cgraph, cmap, *_ogGraph, *_capacity);  // O(N+E) --> expensive, TODO: improve
                // copy a map over

                // if node is not in the connected component of v, it is contracted in the graph copy
                Node t = contractUv(*connectedComponent, cgraph);  // O(N)

                // returns value of min v-t cut
                // sets connected component as an node map for the mincut
                MinCutMap* vtCut = TR::createCutMap(*_ogGraph);

                int dSV = findMincuts(cgraph, v, t, *vtCut, cmap);  // O(Uv(n)^2*Uv(e)^0.5)

                if (dSV < minimalCutValue) {
                    minimalCutValue = dSV;
                    // markCutNodes(*vtCut);
                }
                // TODO set the mincut NodeMap
                
                delete connectedComponent;
                delete vtCut;
            }
        }

        // after finding s-t mincut for ith coloring, need to undo connections for the colorings
        void undoSTComponents() {
            for (IncEdgeIt e(*_graph, _s); e!=INVALID; ++e) {
                (*_graph).erase(e);
            };
            for (IncEdgeIt e(*_graph, _t); e!=INVALID; ++e) {
                (*_graph).erase(e);
            }
        }

        // after each coloring, need to add edges which have not been included in a cut before
        void addToCutUnion(MinCutMap& cutmap, ListGraph& graph) {
            for (EdgeIt e(graph); e != INVALID; ++e) {
                Node u = graph.u(e);
                Node v = graph.v(e);
                
                if (cutmap[u] != cutmap[v]) {
                    if (u == _s || v == _t) {
                        throw logic_error("DEBUG Edge should not be in cut");
                    }
                    (*_allCutsUnion)[e] = true; 
                }
            }
        }

        // adds connects one color of edges to source and other to sink
        void makeSTComponents(int i, ListGraph& graph, Node& s, Node& t, CapacityMap* cap) {
            int len = (*_subset).size();
            for (int j = 1; j <= len; ++j) {
                Node tmp = (*_subset)[j-1]; // make this and id and map back
                if (graph.id(tmp) > graph.maxNodeId() || graph.id(tmp) < 0) {
                    throw domain_error("Subset Node not part of graph.");
                }
                if (findColor(j, i)) {
                    (*cap)[graph.addEdge(s, tmp)] = _maxWeight;
                }
                else {
                    (*cap)[graph.addEdge(t, tmp)] = _maxWeight;
                }
            }
        }

        int _findNumEdges() {
            int numEdges = 0;
            for (EdgeIt e(*_graph); e != INVALID; ++e) {
                ++numEdges;
            }
            return numEdges;
        }
        
        // removes edges for each of the log|V| calls
        void findIthMinCut(int i) {
            makeSTComponents(i, *_graph, _s, _t, _capMap);  // O(|V|) (where V is the subset of nodes)
            MinCutMap* icutmap = TR::createCutMap(*_graph);  // O(N)
            int _ = findMincuts(*_graph, _s, _t, *icutmap, *_capMap);  // O(N^2*E^(0.5))
            addToCutUnion(*icutmap, *_graph); // O(E)
            undoSTComponents(); // O(|V|)
            delete icutmap;
        }

    public:

        void init() {
            // O(N+E)
            createStructures();
        }

        // must call IsolatingCut().init() before .run()
        void run() {
            int len = ceil(log2((*_subset).size()));

            // O(log(|V|)*N^2*E^(0.5))
            for (int i = 0; i < len; ++i) {
                findIthMinCut(i);  // O(N^2*E^(0.5))
            }

            // copyGraph(*_trimGraph, *_ogGraph);
            (*_graph).erase(_s);
            (*_graph).erase(_t);
            removeCuts(*_graph); // O(E)
            parallelContractMinCut(*_graph);  // O(|V|*|N^2e^0.5)
        }

        int getMinCut() {
            return minimalCutValue;
        }
        
        // test methods
        int testIthMinCut(int i) {
            MinCutMap* icutmap = TR::createCutMap(*_graph);
            makeSTComponents(i, *_graph, _s, _t, _capMap);
            int cut = findMincuts(*_graph, _s, _t, *icutmap, *_capMap);
            undoSTComponents();
            delete icutmap;
            return cut;
        }

        int testFindMincuts(ListGraph& graph, Node s, Node t, MinCutMap& cutmap) {
            return findMincuts(graph, s, t, cutmap, *_capacity);
        }

        int getMaxWeight() {
            return _maxWeight;
        }

        void runPhase1() {
            int len = ceil(log2((*_subset).size()));
            for (int i = 0; i < len; ++i) {
                findIthMinCut(i);
            }
            (*_graph).erase(_s);
            (*_graph).erase(_t);
            removeCuts(*_graph);
        }

        // must call in order init(), runPhase1(), then runPhase2()
        void runPhase2() {
            parallelContractMinCut(*_graph);
        }

        // must call in order init(), runPhase1()
        void runFindUv(Node& v, MinCutMap& connectedComponent) {
            findUv(v, connectedComponent, *_graph);
        }

        std::vector<int> runContractAndCut(Node& v, ListGraph& cgraph, MinCutMap& connectedComponent) {
            // this is to validate the S_v is properly computed for any v
            Node t = contractUv(*connectedComponent, cgraph);

            // sets connectedComponents to the edges on the v side of the cut in contracted graph
            int dSV = findMincuts(cgraph, v, t, *connectedComponent);
            std::vector<int> nIds;
            getCutNodeIds(nIds, *connectedComponent, *_ogGraph); 
            return nIds;
        }

    }; //class IsolatingCut
} //namespace lemon


