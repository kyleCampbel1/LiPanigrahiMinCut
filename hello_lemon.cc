#include <iostream>
#include <lemon/list_graph.h>
#include <lemon/hao_orlin.h>
#include <lemon/bits/map_extender.h>
#include <lemon/bits/graph_extender.h>
#include <lemon/nagamochi_ibaraki.h>

using namespace lemon;
using namespace std;
int main()
{
  ListGraph g;
  ListGraph::Node u = g.addNode();
  ListGraph::Node v = g.addNode();
  
  ListGraph::EdgeMap<int> map(g);
  map[g.addEdge(u, v)] = 5;
  HaoOrlin<ListGraph, ListGraph::EdgeMap<int> >  haoOrlin(g, map);
  haoOrlin.run();
  cout << "Hello World! This is LEMON library here." << endl;
  cout << "We have a directed graph with HO mincut " << haoOrlin.minCutValue() << endl;
  NagamochiIbaraki<ListGraph, ListGraph::EdgeMap<int> > nagamochiIbaraki(g, map);
  nagamochiIbaraki.run();
  cout << "NI Mincut " << nagamochiIbaraki.minCutValue() << endl;
  return 0;
}