#include "isolating_cut.h"

namespace lemon {

#ifdef DOXYGEN
  template <typename GR, typename CAP, typename TOL>
#else
  template <typename GR,
            typename CAP = typename GR::template ArcMap<int> >
#endif
    class LiPanigrahi { 
    public:

        typedef GR ListGraph;

        typedef CAP CapacityMap;

        // bring in types from isolating cut
    
    private:
        typedef typename CapacityMap::Value Value;

        TEMPLATE_GRAPH_TYPEDEFS(ListGraph);


    
    public:
        LiPanigrahi(const ListGraph& graph, CapacityMap& capacity, float epsilon = 0.1):
         _capacity(capacity), _ogGraph(graph), _epsilon(epsilon)
         { }

    private:
        float _epsilon;
        ListGraph* _ogGraph;
        CapacityMap* _capacity;
        LiPanigrahiDefaultTraits<ListGraph, ListGraph::EdgeMap<int>>::SV* _minCut;


    public:

        void run() {

        }
  
    }; // LiPanigrahi

} // namespace lemon
/*
figure out unbalanced case
balanced case just randomly choose - look at levine paper
write algo
write random graph
run it
*/