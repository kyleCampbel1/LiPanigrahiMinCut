#include <vector>
#include <unordered_map>
#include <list>
#include <limits>
#include <map>
#include <stdexcept>

#include <math.h>

#include <lemon/list_graph.h>
#include <lemon/concepts/graph.h>
#include <lemon/adaptors.h>
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
    typedef GR::NodeMap<Value> NodeCutMap;
#else
    typedef typename ListGraph::template NodeMap<bool> NodeCutMap;
#endif

    static NodeCutMap* createCutMap(const ListGraph& graph) {
      return new NodeCutMap(graph);
    }

#ifdef DOXYGEN
    typedef GR::EdgeMap<Value> EdgeCutMap;
#else
    typedef typename ListGraph::template EdgeMap<bool> EdgeCutMap;
#endif

    static EdgeCutMap* createCutEdgeMap(const ListGraph& graph) {
      return new EdgeCutMap(graph, true);
    }

typedef class FilterEdges<ListGraph, EdgeCutMap> FilterEdges;

    static FilterEdges* createFilterEdges(ListGraph& graph, EdgeCutMap& filter) {
       return new FilterEdges(graph, filter);
    }

#ifdef DOXYGEN
    typedef GR::NodeMap<Value> NodeCutMapF;
#else
    typedef typename FilterEdges::template NodeMap<bool> NodeCutMapF;
#endif

    static NodeCutMapF* createCutMapF(const FilterEdges& graph) {
      return new NodeCutMapF(graph);
    }

typedef FilterNodes<ListGraph, NodeCutMapF> FilterNodes;

    static FilterNodes* createFilterNodes(ListGraph& graph, NodeCutMapF& filter) {
        return new FilterNodes(graph, filter);
    }

typedef class Preflow<FilterNodes, CapacityMap> Maxflow;

    static Maxflow* createPreflow(FilterNodes& graph, CapacityMap& map, const typename ListGraph::Node& s, const typename ListGraph::Node& t) {
        return new Preflow<FilterNodes, CapacityMap>(graph, map, s, t);
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

        typedef typename Traits::NodeCutMap NodeCutMap;

        typedef typename Traits::NodeCutMapF NodeCutMapF;

        typedef typename Traits::EdgeCutMap EdgeCutMap;

        typedef typename Traits::FilterEdges FilterEdges;

        typedef typename Traits::FilterNodes FilterNodes;

        typedef typename Traits::Maxflow Maxflow;
       
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
        _capacity(&capacity), _subset(&subset), _graph(&graph), _max_weight(graphSumWeight(capacity))
    {}

    ~IsolatingCut() {
        delete _all_cuts_union;
        delete _cap_map;
        delete _graph_w;
        delete fe;
        delete _cut_map;

    }
    
    private: 
        // graph pointer to original graph
        const ListGraph* _graph;
        // graph we will find s-t cuts and create edges on
        ListGraph* _graph_w;
        CapacityMap* _capacity;
        Node _s, _t;
        int _max_weight;
        int last_edge;

        NodeCutMap* _cut_map;
        EdgeCutMap* _all_cuts_union;
        EdgeCutMap* st_connections;
        CapacityMap* _cap_map;
        FilterEdges* fe;
        FilterEdges* st_graph;


        int _min_cut = INT_MAX;

        std::vector<Node>* _subset;

 
    private:
        void createStructures() {

            // structures for storing results from log|V| maxflow comps
            _graph_w = new ListGraph;
            _cap_map = createCapactiyMap(*_graph_w);
            copyGraph(*_graph_w, *_cap_map, *_graph, *_capacity);
            _all_cuts_union = TR::createCutEdgeMap(*_graph_w);
            st_connections = TR::createCutEdgeMap(*_graph_w);
            _cut_map = TR::createCutMap(*_graph);

            _s = (*_graph_w).addNode();
            _t = (*_graph_w).addNode();
            fe = TR::createFilterEdges(*_graph_w, *_all_cuts_union);
            st_graph = TR::createFilterEdges(*_graph_w, *st_connections);
            last_edge = connectST(*_graph_w, _s, _t, _cap_map);
        }

        void copyGraph(ListGraph& new_graph, CapacityMap& new_map, const ListGraph& old_graph, CapacityMap& old_map) {
            for (int i = 0; i <= old_graph.maxNodeId(); ++i) {
                new_graph.addNode();
            }

            for (ListGraph::EdgeIt e(old_graph); e != INVALID; ++e) {
                new_map[new_graph.addEdge(old_graph.u(e), old_graph.v(e))] = old_map[e];
            }
        }

        // performs vertex labeling given vertex index and i of iteration
        bool findColor(int label, int i) {
            bool val = label & (1 << i);
            return val;
        }

        int graphSumWeight(CapacityMap& capacity_map) {
            int weight = 0;
            for (EdgeIt e(*_graph); e != INVALID; ++e) {
                weight += capacity_map[e];
            }
            return weight;
        }

        // this finds a min s-t and returns value of the cut. sets cutmap to bools for nodes in/out of cut
        int findSTCut(ListGraph& graph, Node s, Node t, NodeCutMap& cutmap, CapacityMap& cap) {
            Preflow<ListGraph, CapacityMap> pf(graph, cap, s, t);
            pf.run();
            pf.minCutMap(cutmap);
            return pf.flowValue();
        }

        int findSTCut(FilterEdges& graph, Node s, Node t, NodeCutMap& cutmap, CapacityMap& cap) {
            Preflow<FilterEdges, CapacityMap> pf(graph, cap, s, t);
            pf.run();
            pf.minCutMap(cutmap);
            return pf.flowValue();
        }

        // performs dfs to set compMap to the nodes reachable from v in G\ U C_i
        void dfs(const Node& v, NodeCutMapF& compMap, const FilterEdges& isolatedCuts) {
            Dfs<FilterEdges> dfs(isolatedCuts); // could pass this in as arg because it gets called often
            dfs.reachedMap(compMap);
            dfs.run(v);
        }

        void spanComponents(FilterNodes& graph, Node& t, NodeCutMapF& component, CapacityMap& capacity) {
            for (typename FilterNodes::NodeIt n(graph); n != INVALID; ++n) {
                // cout << graph.id(n) << " nidtrim " << (*_graph_w).id(n) << " nid" << endl;
                if (n != t) {
                    for (IncEdgeIt e(*_graph, n); e != INVALID; ++e) {
                        if (component[(*_graph).u(e)] != component[(*_graph).v(e)]) {
                            capacity[(*_graph_w).addEdge(n, t)] = (*_capacity)[e];  // be careful to avoid duplicates
                        }
                    }
                }
            }
        }

        // G\ U C_i is found, this finds isolating cuts
        void parallelContractMinCut(FilterEdges& isolated_components) {
            // in parallel loop -- TODO
            // check runtime of this loop
            NodeCutMapF* connected_component = TR::createCutMapF(isolated_components);
            FilterNodes* fn = TR::createFilterNodes(*_graph_w, *connected_component);
                
            for (const auto& v : *_subset) {
                // finds connected component that contains v in pruned graph
                dfs(v, *connected_component, isolated_components);  // O(|Uv(n+e)|)
                
                Node t = (*_graph_w).addNode();
                fn->status(t, true);
                // replicates graph contraction
                spanComponents(*fn, t, *connected_component, *_cap_map);

                Maxflow* pf = TR::createPreflow(*fn, *_cap_map, v, t);
                pf->run();
                int dSV = pf->flowValue();

                if (dSV < _min_cut) {
                    _min_cut = dSV;
                    (*pf).minCutMap(*_cut_map);
                }
                // TODO set the mincut NodeMap


                (*_graph_w).erase(t); 
                delete pf;
            }
            delete connected_component;
            delete fn;
            
        }

        // after each coloring, need to add edges which have not been included in a cut before
        void addToCutUnion(NodeCutMap& cutmap, FilterEdges& graph) {
            for (typename FilterEdges::EdgeIt e(graph); e != INVALID; ++e) {
                Node u = graph.u(e);
                Node v = graph.v(e);
                
                if (cutmap[u] != cutmap[v]) {
                    if (u == _s || v == _t) {
                        throw logic_error("DEBUG Edge should not be in cut");
                    }
                    (*_all_cuts_union)[e] = false;
                }
            }
        }

        int connectST(ListGraph& graph, Node& s, Node& t, CapacityMap* cap){
            int startId = (*_graph_w).maxEdgeId()+1;
            int len = (*_subset).size();
            for (int j = 1; j <= len; ++j) {
                Node tmp = (*_subset)[j-1]; // make this and id and map back
                if (graph.id(tmp) > graph.maxNodeId() || graph.id(tmp) < 0) {
                    throw domain_error("Subset Node not part of graph.");
                }
                (*cap)[graph.addEdge(s, tmp)] = _max_weight;
                (*cap)[graph.addEdge(t, tmp)] = _max_weight;
            }
            return startId;
        }

        // adds connects one color of edges to source and other to sink
        void makeSTComponents(int i, int startIndex) {
            int len = (*_subset).size();
            for (int j = 1; j <= len; ++j) {
                bool sflag = findColor(j, i);
                (*st_connections)[(*_graph_w).edgeFromId(startIndex+2*(j-1))] = sflag;
                (*st_connections)[(*_graph_w).edgeFromId(startIndex+2*(j-1)+1)] = !sflag;
            }
        }
        
        // removes edges for each of the log|V| calls
        void findIthMinCut(int i) {
            makeSTComponents(i, last_edge);  // O(|V|) (where V is the subset of nodes)
            NodeCutMap* icutmap = TR::createCutMap(*_graph_w);  // O(N)
            int _ = findSTCut(*st_graph, _s, _t, *icutmap, *_cap_map);  // O(N^2*E^(0.5))
            addToCutUnion(*icutmap, *st_graph); // O(E)
            delete icutmap;
        }

    public:

        void init() {
            // O(N+E)
            createStructures();
        }

        // must call IsolatingCut().init() before .run()
        // this call calculates the minimum cut
        void run() {
            int len = ceil(log2((*_subset).size()));

            // O(log(|V|)*N^2*E^(0.5))
            for (int i = 0; i < len; ++i) {
                findIthMinCut(i);  // O(N^2*E^(0.5))
            }
            (*_graph_w).erase(_s);
            (*_graph_w).erase(_t);
            
            parallelContractMinCut(*fe);  // O(|V|*|N^2e^0.5)
        }

        int minCutValue() const {
            return _min_cut;
        }

        int minCutMap(NodeCutMap& cutmap) const {
            for (ListGraph::NodeIt n(*_graph); n != INVALID; ++n) {
                cutmap.set(n, (*_cut_map)[n]);
            }
            return minCutValue();
        }
        
        // test methods
        int testIthMinCut(int i) {
            NodeCutMap* icutmap = TR::createCutMap(*_graph_w);
            makeSTComponents(i, last_edge);
            int cut = findSTCut(*_graph_w, _s, _t, *icutmap, *_cap_map);
            delete icutmap;
            return cut;
        }

        int testFindMincuts(ListGraph& graph, Node s, Node t, NodeCutMap& cutmap) {
            return findSTCut(graph, s, t, cutmap, *_capacity);
        }

        int getMaxWeight() {
            return _max_weight;
        }

        void runPhase1() {
            int len = ceil(log2((*_subset).size()));
            for (int i = 0; i < len; ++i) {
                findIthMinCut(i);
            }
            // remove auxilliary nodes before finding mincut in graph
            (*_graph_w).erase(_s);
            (*_graph_w).erase(_t);
        }

        // must call in order init(), runPhase1(), then runPhase2()
        void runPhase2() {
            parallelContractMinCut(*fe);
        }

        // must call in order init(), runPhase1()
        void runFindUv(Node& v, NodeCutMapF& connected_component) {
            dfs(v, connected_component, *fe);
        }

        FilterEdges& getFilterSubset() {
            return *fe;
        }

    }; //class IsolatingCut
} //namespace lemon


