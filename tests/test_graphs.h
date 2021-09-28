#include <iostream>
#include <vector>

#include <lemon/maps.h>
#include <lemon/hao_orlin.h>
#include <lemon/nagamochi_ibaraki.h>
#include <lemon/list_graph.h>

using namespace lemon;
using namespace std;

class TestGraphs {

    private:
    ListGraph* _graph;
    ListGraph::EdgeMap<int>* _capacity;

    public:
    TestGraphs(ListGraph* graph, ListGraph::EdgeMap<int>* capacity) :
    _graph(graph), _capacity(capacity) { }

    ~TestGraphs() { }

    void createBicycleWheel(int n) {
        ListGraph::Node c1 = (*_graph).addNode();
        ListGraph::Node c2 = (*_graph).addNode(); // center nodes of wheel
        ListGraph::Node firstNodeOfWheel = (*_graph).addNode();
        (*_capacity)[(*_graph).addEdge(c1, c2)] = 1;
        (*_capacity)[(*_graph).addEdge(firstNodeOfWheel, c1)] = 1;
        ListGraph::Node prev = firstNodeOfWheel;
        ListGraph::Node cur;
        for (int i = 3; i < n; ++ i) {
          cur = (*_graph).addNode();
          (*_capacity)[(*_graph).addEdge(prev, cur)] = n/2;
          if (i % 2 == 0) {
            (*_capacity)[(*_graph).addEdge(c2, cur)] = 1;
          }
          else {
            (*_capacity)[(*_graph).addEdge(c1, cur)] = 1;
          }
          prev = cur;
        }
        (*_capacity)[(*_graph).addEdge(firstNodeOfWheel, prev)] = n/2;
    }

    ListGraph& getGraph() {
        return *_graph;
    }

    ListGraph::EdgeMap<int>& getCapacityMap() {
        return *_capacity;
    }
};