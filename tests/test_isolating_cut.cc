#include "test_graphs.h"
#include "../isolating_cut.h"


void testInitLabels() {

};

void testColorContradction() {

};

void testColorGraph() {

};

int main()
{
  ListGraph g;
  ListGraph::EdgeMap<int> map(g);
  TestGraphs test(&g, &map);
  test.createBicycleWheel(128);
  HaoOrlin<ListGraph, ListGraph::EdgeMap<int> >  haoOrlin(g, map);
  haoOrlin.run();
  cout << "Hello World! This is LEMON library here." << endl;
  cout << "We have a directed graph with HO mincut " << haoOrlin.minCutValue() << endl;
  NagamochiIbaraki<ListGraph, ListGraph::EdgeMap<int> > nagamochiIbaraki(test.getGraph(), test.getCapacityMap());
  nagamochiIbaraki.run();
  cout << "NI Mincut " << nagamochiIbaraki.minCutValue() << endl;
  for (ListGraph::NodeIt n(g); n != INVALID; ++n) {
    int cnt = 0;
    for (ListGraph::IncEdgeIt e(g, n); e != INVALID; ++e) {
        cnt++;
    }
    // std::cout << "deg(" << g.id(n) << ") = " << cnt << std::endl;
  };

  return 0;
}