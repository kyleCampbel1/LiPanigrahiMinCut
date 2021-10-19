#include "test_graphs.h"
#include "../isolating_cut.h"
#include <lemon/list_graph.h>
#include <lemon/concepts/graph.h>
#include <lemon/concepts/bpgraph.h>
#include <lemon/bits/graph_extender.h>
#include <vector>
#include <assert.h>
#include <iostream>

using namespace std;
using namespace lemon;

TEMPLATE_GRAPH_TYPEDEFS(ListGraph);
typedef ListGraph::EdgeMap<int> CapacityMap;
typedef IsolatingCut<ListGraph, CapacityMap> IsoCut;

void testInitLabels() {
  ListGraph g;
  ListGraph::EdgeMap<int> map(g);
  TestGraphs test(&g, &map);
  test.createBicycleWheel(8);
  vector<Node> subset = {g.nodeFromId(2), g.nodeFromId(4), g.nodeFromId(7)};
  ListGraph c1, c2;
  IsoCut ic(g, map, subset, c1, c2);
  ic.init();
  cout << ic.testIthMinCut(3) << endl;
  // vector<SV> isolatingCuts = ic.getIsolatingCuts();
  // cout << isolatingCuts.size() << endl;
  // for (auto& cut : isolatingCuts) {
  //   std::cout << cut.weight << endl;
  //   for (int nodeId : cut.sv) {
  //     std::cout << nodeId << endl;
  //   }
  // }
};

void testMaxWeight() {
  ListGraph g;
  ListGraph::EdgeMap<int> map(g);
  TestGraphs test(&g, &map);
  test.createBicycleWheel(8);
  vector<Node> subset = {g.nodeFromId(2), g.nodeFromId(4), g.nodeFromId(7)};
  ListGraph c1, c2;
  IsoCut ic(g, map, subset, c1, c2);
  ic.init();
  cout << ic.getMaxWeight() << endl;
};

void testRun() {
  ListGraph g;
  ListGraph::EdgeMap<int> map(g);
  TestGraphs test(&g, &map);
  test.createBicycleWheel(16);
  vector<Node> subset = {g.nodeFromId(7), g.nodeFromId(4), g.nodeFromId(12), g.nodeFromId(14)};
  ListGraph c1, c2;
  IsoCut ic(g, map, subset, c1, c2);
  ic.init();
  ic.run();
  vector<LiPanigrahiDefaultTraits<ListGraph, ListGraph::EdgeMap<int>>::SV> isolatingCuts = ic.getIsolatingCuts();
};

void testDisconnection() {
  ListGraph g;
  ListGraph::EdgeMap<int> map(g);
  TestGraphs test(&g, &map);
  test.createBicycleWheel(8);
  vector<Node> subset = { g.nodeFromId(2), g.nodeFromId(4), g.nodeFromId(1)};
  ListGraph c1, c2;
  IsoCut ic(g, map, subset, c1, c2);
  ic.init();
  ListGraph::NodeMap<bool> cutMap(g);
  Node v = g.nodeFromId(2);
  Node t = g.nodeFromId(4);
  ic.testPhase1(v, cutMap, t);
}

void testSetup() {
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
    std::cout << "deg(" << g.id(n) << ") = " << cnt << std::endl;
  }
};

void testConstructor() {
  ListGraph g;
  ListGraph::EdgeMap<int> map(g);
  TestGraphs test(&g, &map);
  test.createBicycleWheel(8);
  vector<Node> subset = {g.nodeFromId(2), g.nodeFromId(4), g.nodeFromId(7)};
  ListGraph c1, c2;
  IsoCut ic(g, map, subset, c1, c2);
  ic.init();
  // assert (ic.getMaxWeight > 0);
}

int main() {
  // testInitLabels();
  testRun();

  // ListGraph g;
  // ListGraph::EdgeMap<int> map(g);
  // TestGraphs test(&g, &map);
  // test.createBicycleWheel(8);
  // ListGraph* t;
  // ListGraph a;
  // ExtendedListGraphBase::NodeMap<ListGraphBase::Node> nr(g);
  // ExtendedListGraphBase::EdgeMap<ListGraphBase::Edge> ecr(a);
  // GraphCopy<ListGraph, ListGraph> gc(g, a);
  // gc.nodeRef(nr).edgeCrossRef(ecr);
  // gc.run();
  // t = &a;
  // int c = 0;
  return 0;
};