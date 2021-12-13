#include "test_graphs.h"
#include "../isolating_cut.h"

#include <lemon/list_graph.h>

#include <lemon/hao_orlin.h>
#include <lemon/nagamochi_ibaraki.h>
#include <vector>
#include <unordered_set>
#include <assert.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>


using namespace std;
using namespace lemon;

TEMPLATE_GRAPH_TYPEDEFS(ListGraph);
typedef ListGraph::EdgeMap<int> CapacityMap;
typedef IsolatingCut<ListGraph> IsoCut;

void testInitLabels() {
  ListGraph g;
  ListGraph::EdgeMap<int> map(g);
  TestGraphs test(&g, &map);
  test.createBicycleWheel(8);
  vector<Node> subset = {g.nodeFromId(2), g.nodeFromId(4), g.nodeFromId(7)};
  IsoCut ic(g, map, subset);
  ic.init();
  cout << ic.testIthMinCut(3) << endl;

};


void testMaxWeight() {
  ListGraph g;
  ListGraph::EdgeMap<int> map(g);
  TestGraphs test(&g, &map);
  test.createBicycleWheel(8);
  vector<Node> subset = {g.nodeFromId(2), g.nodeFromId(4), g.nodeFromId(7)};
  IsoCut ic(g, map, subset);
  ic.init();
  cout << ic.getMaxWeight() << endl;
};

void testDisconnection() {
  ListGraph g;
  ListGraph::EdgeMap<int> map(g);
  TestGraphs test(&g, &map);
  test.createBicycleWheel(8);
  vector<Node> subset = { g.nodeFromId(2), g.nodeFromId(4), g.nodeFromId(1)};
  IsoCut ic(g, map, subset);
  ic.init();
  ListGraph::NodeMap<bool> cutMap(g);
  Node v = g.nodeFromId(2);
  Node t = g.nodeFromId(4);
  // ic.testPhase1(v, cutMap, t);
}

void testSetupWheel(int n) {
  ListGraph g;
  ListGraph::EdgeMap<int> map(g);
  TestGraphs test(&g, &map);
  test.createBicycleWheel(n);
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
}

void testSetupRandomGraph() {
  ListGraph g;
  ListGraph::EdgeMap<int> map(g);
  TestGraphs test(&g, &map);
  int n = 128;
  test.createRandomGraph(n, 0.5);
  HaoOrlin<ListGraph, ListGraph::EdgeMap<int> >  haoOrlin(g, map);
  haoOrlin.run();
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
  cout << "Max Edge Id: " << g.maxEdgeId() << endl;
  cout << "max node id: " << g.maxNodeId() << endl;
  cout << "Expected nodes: " << n << endl;

};

void testPhase1(ListGraph& g, ListGraph::EdgeMap<int>& map, vector<Node>& subset) {
  IsoCut ic(g, map, subset);
  ic.init();
  ic.runPhase1();
  IsoCut::NodeCutMapF* connectedComponent = IsoCut::Traits::createCutMapF(ic.getFilterSubset());
  for (auto& v : subset) {
    ic.runFindUv(v, *connectedComponent);
    bool successPhase1 = true;
    for (auto& u : subset) {
      if (u != v && (*connectedComponent)[u] == true) {
        cout <<  "v = " << g.id(v) << ". node misclassified: " << g.id(u) << endl;
        successPhase1 = false;
      }
    }
    cout << "phase 1 success: " << successPhase1 << endl;
  }
  delete connectedComponent;
}

void testConstructor() {
  ListGraph g;
  ListGraph::EdgeMap<int> map(g);
  TestGraphs test(&g, &map);
  test.createBicycleWheel(8);
  vector<Node> subset = {g.nodeFromId(2), g.nodeFromId(4), g.nodeFromId(7)};
  IsoCut ic(g, map, subset);
  ic.init();
}

void testFindCut() {
  ListGraph g;
  ListGraph::EdgeMap<int> map(g);
  TestGraphs test(&g, &map);
  test.createBicycleWheel(16);
  vector<Node> subset = {g.nodeFromId(5), g.nodeFromId(8), g.nodeFromId(11)};
  testPhase1(g, map, subset);
  IsoCut ic(g, map, subset);
  ic.init();
  ic.run();
  int mincut = ic.minCutValue();
  cout << mincut << endl;
  IsoCut::NodeCutMap* nm = IsoCut::Traits::createCutMap(g);
  ic.minCutMap(*nm);
  for (NodeIt n(g); n != INVALID; ++n) {
    cout << (*nm)[n] << endl;
  }
  delete nm;

  ListGraph g1;
  ListGraph::EdgeMap<int> map1(g1);
  TestGraphs test1(&g1, &map1);
  test1.createRandomGraph(42, 0.5);
  vector<Node> subset1 = {g1.nodeFromId(16), g1.nodeFromId(40), g1.nodeFromId(30), g1.nodeFromId(28), g1.nodeFromId(12)};
  testPhase1(g1, map1, subset1);
  IsoCut ic1(g1, map1, subset1);
  ic1.init();
  ic1.run();
  int mincut1 = ic1.minCutValue();
  cout << mincut1 << endl;
  IsoCut::NodeCutMap* nm1 = IsoCut::Traits::createCutMap(g1);
  ic.minCutMap(*nm1);
  for (NodeIt n(g); n != INVALID; ++n) {
    cout << (*nm1)[n] << endl;
  }
  delete nm;
}

void testRun() {
  ListGraph g;
  ListGraph::EdgeMap<int> map(g);
  TestGraphs test(&g, &map);
  test.createBicycleWheel(16);
  vector<Node> subset = {g.nodeFromId(5), g.nodeFromId(8), g.nodeFromId(11)};
  testPhase1(g, map, subset);
  IsoCut ic(g, map, subset);
  ic.init();
  ic.run();
  int mincut = ic.minCutValue();
  cout << mincut << endl;

  ListGraph g1;
  ListGraph::EdgeMap<int> map1(g1);
  TestGraphs test1(&g1, &map1);
  test1.createRandomGraph(512, 0.5);
  vector<Node> subset1 = {g1.nodeFromId(16), g1.nodeFromId(264), g1.nodeFromId(409), g1.nodeFromId(166), g1.nodeFromId(306), g1.nodeFromId(90), g1.nodeFromId(76), g1.nodeFromId(28), g1.nodeFromId(458), g1.nodeFromId(322), g1.nodeFromId(226), g1.nodeFromId(285)};
  testPhase1(g1, map1, subset1);
  IsoCut ic1(g1, map1, subset1);
  ic1.init();
  ic1.run();
  int mincut1 = ic1.minCutValue();
  cout << mincut1 << endl;

  ListGraph g2;
  ListGraph::EdgeMap<int> map2(g2);
  TestGraphs test2(&g2, &map2);
  test2.createBicycleWheel(64);
  vector<Node> subset2 = {g2.nodeFromId(16), g2.nodeFromId(62), g2.nodeFromId(9), g2.nodeFromId(26), g2.nodeFromId(28)};
  testPhase1(g2, map2, subset2);
  IsoCut ic2(g2, map2, subset2);
  ic2.init();
  ic2.run();
  int mincut2 = ic2.minCutValue();
  cout << mincut2 << endl;

  ListGraph g3;
  ListGraph::EdgeMap<int> map3(g3);
  TestGraphs test3(&g3, &map3);
  test3.createRandomGraph(64, 0.5);
  vector<Node> subset3 = {g3.nodeFromId(7), g3.nodeFromId(22), g3.nodeFromId(19), g3.nodeFromId(45), g3.nodeFromId(60)};
  testPhase1(g3, map3, subset3);
  IsoCut ic3(g3, map3, subset3);
  ic3.init();
  ic3.run();
  int mincut3 = ic3.minCutValue();
  cout << mincut3 << endl;

};

int benchMarkLPTest(const int n, const int subsetSize, ListGraph& graph, CapacityMap& capacityMap) {

  vector<int> nodes(n);
  for (int i = 0; i < n; ++i) {
    nodes[i] = i;
  }

  // O(n) in number of nodes - looking for more efficient shuffling
  std::shuffle(nodes.begin(), nodes.end(), std::mt19937{std::random_device{}()});
  
  vector<Node> subset(subsetSize);
  for (int i = 0; i < subsetSize; ++i) {
    subset[i] = graph.nodeFromId(nodes[i]);
  }
  auto start = chrono::high_resolution_clock::now();
  IsoCut ic(graph, capacityMap, subset);
  ic.init();
  ic.run();
  auto stop = chrono::high_resolution_clock::now();
  chrono::duration<double, milli> runtime_ms = stop-start;
  cout << "Li Panigrahi exec. time (ms): " << runtime_ms.count() << " Nodes:" << n << " Subset Size" << subsetSize << endl;
  return ic.minCutValue();
}

int benchMarkHOTest(const int n, ListGraph& graph, CapacityMap& capacityMap) {
  auto start = chrono::high_resolution_clock::now();
  HaoOrlin<ListGraph, CapacityMap> ho(graph, capacityMap);
  ListGraph::NodeMap<bool> mincut(graph);
  ho.run();
  auto stop = chrono::high_resolution_clock::now();
  chrono::duration<double, milli> runtime_ms = stop-start;
  cout << "Hao Orlin exec. time (ms): " << runtime_ms.count() << " Nodes:" << n << endl;
  return ho.minCutValue();
}

int benchMarkNITest(const int n, ListGraph& graph, CapacityMap& capacityMap) {
  auto start = chrono::high_resolution_clock::now();
  NagamochiIbaraki<ListGraph, CapacityMap> ni(graph, capacityMap);
  ListGraph::NodeMap<bool> mincut(graph);
  ni.run();
  int cut = ni.minCutMap(mincut);
  auto stop = chrono::high_resolution_clock::now();
  chrono::duration<double, milli> runtime_ms = stop-start;
  cout << "Nagamochi Ibaraki exec. time (ms): " << runtime_ms.count() << " Nodes:" << n << endl;
  return ni.minCutValue();
}

void runBenchmarks() {
  ListGraph g1;
  ListGraph::EdgeMap<int> map1(g1);
  TestGraphs test(&g1, &map1);
  test.createBicycleWheel(1024);
  int nim = benchMarkNITest(1024, g1, map1);
  int hom = benchMarkHOTest(1024, g1, map1);
  int lpm = INT_MAX;
  for (int s : {12, 100, 200, 500, 700, 800}) {
    lpm = min(benchMarkLPTest(1024, s, g1, map1), lpm);
  }
  cout << "NI | HO | LP " << nim << " | " << hom << " | " << lpm << endl;

  ListGraph g2;
  ListGraph::EdgeMap<int> map2(g2);
  TestGraphs test2(&g2, &map2);
  test2.createBicycleWheel(10000);
  int ni2 = benchMarkNITest(10000, g2, map2);
  int ho2 = benchMarkHOTest(10000, g2, map2);
  int lpm2 = INT_MAX;
  for (int s : {50, 250, 1000, 4000, 5000, 7500}) {
    lpm2 = min(benchMarkLPTest(10000, s, g2, map2), lpm2);
  }
  cout << "NI | HO | LP " << ni2 << " | " << ho2 << " | " << lpm2 << endl;

  ListGraph g4;
  ListGraph::EdgeMap<int> map4(g4);
  TestGraphs test4(&g4, &map4);
  test4.createRandomGraph(1024, 0.5);
  int ni4 = benchMarkNITest(1024, g4, map4);
  int ho4 = benchMarkHOTest(1024, g4, map4);
  int lpm4 = INT_MAX;
  for (int s : {50, 250, 500, 750}) {
    lpm4 = min(benchMarkLPTest(1024, s, g4, map4), lpm4);
  }
  cout << "NI | HO | LP " << ni4 << " | " << ho4 << " | " << lpm4 << endl;

  ListGraph g5;
  ListGraph::EdgeMap<int> map5(g5);
  TestGraphs test5(&g5, &map5);
  test5.createRandomGraph(1024, 0.75);
  int ni5 = benchMarkNITest(1024, g5, map5);
  int ho5 = benchMarkHOTest(1024, g5, map5);
  int lpm5 = INT_MAX;
  for (int s : {12, 100, 500, 700}) {
    lpm5 = min(benchMarkLPTest(1024, s, g5, map5), lpm5);
  }
  cout << "NI | HO | LP " << ni5 << " | " << ho5 << " | " << lpm5 << endl;

  ListGraph g6;
  ListGraph::EdgeMap<int> map6(g6);
  TestGraphs test6(&g6, &map6);
  test6.createRandomGraph(1024, 0.25);
  int ni6 = benchMarkNITest(1024, g4, map4);
  int ho6 = benchMarkHOTest(1024, g4, map4);
  int lpm6 = INT_MAX;
  for (int s : {12, 100, 500, 700}) {
    lpm6 = min(benchMarkLPTest(1024, s, g4, map4), lpm6);
  }
  cout << "NI | HO | LP " << ni6 << " | " << ho6 << " | " << lpm6 << endl;
  
  ListGraph g3;
  ListGraph::EdgeMap<int> map3(g3);
  TestGraphs test3(&g3, &map3);
  test3.createBicycleWheel(100000);
  int ni3 = benchMarkNITest(100000, g3, map3);
  // int ho3 = benchMarkHOTest(100000, g3, map3);
  int lpm3 = INT_MAX;
  for (int s : {500, 5000, 10000, 40000}) {
    lpm3 = min(benchMarkLPTest(100000, s, g3, map3), lpm3);
  }
  cout << "NI | HO | LP " << ni3 << " | " << "NA" << " | " << lpm3 << endl;

  ListGraph g7;
  ListGraph::EdgeMap<int> map7(g7);
  TestGraphs test7(&g7, &map7);
  test7.createRandomGraph(10000, 0.75);
  int ni7 = benchMarkNITest(10000, g7, map7);
  // int ho7 = benchMarkHOTest(10000, g7, map7);
  int lpm7 = INT_MAX;
  for (int s : {100, 500, 2500, 5000}) {
    lpm7 = min(benchMarkLPTest(10000, s, g7, map7), lpm7);
  }
  cout << "NI | HO | LP " << ni7 << " | " << "NA" << " | " << lpm7 << endl;
}



void testSmallGraphs() {
  ListGraph g1;
  ListGraph::EdgeMap<int> map1(g1);
  TestGraphs test1(&g1, &map1);
  test1.createBicycleWheel(64);
  int mincut1 = INT_MAX;
  for (int s : {4, 10, 20, 40}) {
    mincut1 = benchMarkLPTest(64, s, g1, map1);
  }
  int ni1 = benchMarkNITest(64, g1, map1);

  cout << "Success: " << (ni1 == mincut1) << endl;
  cout << "NI: " << ni1 << " | LP: " << mincut1 << endl;

  ListGraph g5;
  ListGraph::EdgeMap<int> map5(g5);
  TestGraphs test5(&g5, &map5);
  test5.createRandomGraph(64, 0.5);
  int mincut2 = INT_MAX;
  for (int s : {4, 10, 20, 40}) {
    mincut2 = min(mincut2, benchMarkLPTest(64, s, g5, map5));
  }
  int ni2 = benchMarkNITest(64, g5, map5);
  cout << "Success: " << (ni2 == mincut2) << endl;
  cout << "NI: " << ni2 << " | LP: " << mincut2 << endl;
}

void copyGraphAndCapacity(ListGraph& newGraph, CapacityMap& newMap, ListGraph& oldGraph, CapacityMap& oldMap) {
  ExtendedListGraphBase::NodeMap<ListGraphBase::Node> nr(oldGraph);
  ExtendedListGraphBase::EdgeMap<ListGraphBase::Edge> ecr(newGraph);
  GraphCopy<ListGraph, ListGraph> gc(oldGraph, newGraph);
  gc.nodeRef(nr).edgeCrossRef(ecr);
  gc.run();
  // for (int i = 0; i <= oldGraph.maxNodeId(); ++i) {
  //   newGraph.addNode();
  // }

  for (ListGraph::EdgeIt e(oldGraph); e != INVALID; ++e) {
    newMap[newGraph.addEdge(oldGraph.u(e), oldGraph.v(e))] = oldMap[e];
  }
}

int graphSumWeight(CapacityMap& capacityMap, ListGraph& g) {
  int weight = 0;
  for (EdgeIt e(g); e != INVALID; ++e) {
      weight += capacityMap[e];
  }
  
  return weight;
}

void copyGraphFindSTCut() {
  ListGraph g;
  ListGraph::EdgeMap<int> map(g);
  TestGraphs test(&g, &map);
  test.createBicycleWheel(8);

  // copying graph
  ListGraph copy;
  ListGraph::EdgeMap<int> copymap(copy);
  copyGraphAndCapacity(copy, copymap, g, map);

  ListGraph::Node s = g.addNode();
  ListGraph::Node t = g.addNode();
  map[g.addEdge(t, g.nodeFromId(1))] = 31;
  map[g.addEdge(s, g.nodeFromId(5))] = 31;

  ListGraph::Node s1 = copy.addNode();
  ListGraph::Node t1 = copy.addNode();
  copymap[copy.addEdge(t1, copy.nodeFromId(1))] = 31;
  copymap[copy.addEdge(s1, copy.nodeFromId(5))] = 31;

  IsoCut::Traits::NodeCutMap* cutMap = new IsoCut::Traits::NodeCutMap(g);
  Preflow<ListGraph, CapacityMap> pf(g, map, s, t);
  pf.run();
  pf.minCutMap(*cutMap);
  int cutval = pf.flowValue();

  IsoCut::Traits::NodeCutMap* cutMap1 = new IsoCut::Traits::NodeCutMap(copy);
  Preflow<ListGraph, CapacityMap> pf1(copy, copymap, s1, t1);
  pf1.run();
  pf1.minCutMap(*cutMap1);
  int cutval1 = pf1.flowValue();

  bool success = (cutval == cutval1);
  cout << "Success: " << success << endl;
  cout << "first cut: " << cutval << endl;
  cout << "copy cut: " << cutval1 << endl;

  delete cutMap1;
  delete cutMap;
}



void testSTCut() {
  ListGraph g;
  ListGraph::EdgeMap<int> map(g);
  TestGraphs test(&g, &map);
  test.createBicycleWheel(8);
  ListGraph::Node s = g.addNode();
  ListGraph::Node t = g.addNode();
  map[g.addEdge(t, g.nodeFromId(1))] = 31;
  map[g.addEdge(s, g.nodeFromId(5))] = 31;

  IsoCut::Traits::NodeCutMap* cutMap = new IsoCut::Traits::NodeCutMap(g);
  Preflow<ListGraph, CapacityMap> pf(g, map, s, t);
  pf.run();
  pf.minCutMap(*cutMap);
  
  int mincutval = pf.flowValue();
  for (EdgeIt e(g); e != INVALID; ++e) {
      Node u = g.u(e);
      Node v = g.v(e);
      if ((*cutMap)[u] != (*cutMap)[v]) {
          g.erase(e);
      }
  }
  IsoCut::Traits::NodeCutMap* dcutmap = new IsoCut::Traits::NodeCutMap(g);
  Dfs<ListGraph> dfs(g);
  dfs.reachedMap(*dcutmap);
  dfs.run(g.nodeFromId(1));

  cout << mincutval << endl;

  delete cutMap;
  delete dcutmap;
}

void testSTComponents() {
  ListGraph g;
  ListGraph::EdgeMap<int> map(g);
  TestGraphs test(&g, &map);
  test.createBicycleWheel(64);

  vector<Node> subset = {g.nodeFromId(44), g.nodeFromId(11), g.nodeFromId(40), g.nodeFromId(56), g.nodeFromId(1), g.nodeFromId(53), g.nodeFromId(63), g.nodeFromId(34), g.nodeFromId(33)};
  IsoCut ic(g, map, subset);
  ic.init();

  int mincutval2 = ic.testIthMinCut(2);
  int mincutval = ic.testIthMinCut(0);
  int mincutval1 = ic.testIthMinCut(1);
  
  cout << mincutval1 << endl;
}

int main() {
  // testInitLabels();
  // testRun();
  // testFindCut();
  // testSTCut();
  // testSTComponents();
  // recreateSTCut();
  // testSetupRandomGraph();
  // testSmallGraphs();
  runBenchmarks();
  // // testSetupWheel(65);
  // copyGraphFindSTCut();
  // testSpanComponents();
  return 0;
};