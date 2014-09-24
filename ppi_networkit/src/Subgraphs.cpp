#include "Subgraphs.h"

namespace tsppi {

NetworKit::Graph Subgraphs::getSubgraph(count_t id) const
{
  // Copy graph and remove all edges not expressed

  // get a copy of the graph (this performs a deep copy)
  NetworKit::Graph result(this->graph);

  result.forEdges([&](NetworKit::node u, NetworKit::node v){
    if (this->node_exists[u][id] && this->node_exists[v][id])
    {
      // do nothing, because edge is expressed
    }
    else
    {
      // remove edge
      result.removeEdge(u, v);
    }
  });
  return result;
}

} // namespace tsppi
