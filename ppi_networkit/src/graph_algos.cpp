#include "graph_algos.h"


namespace tsppi {
namespace algo {

/*
 * Create vector of degrees for the Graph and return this.
 */
std::vector<NetworKit::count> graph_degrees(const NetworKit::Graph& G)
{
    // init the output:
    std::vector<NetworKit::count> degrees(G.numberOfNodes());

    // simply get the degree from the graph for each node
    // each lookup is in constant time, especially not dependent on m=|E|
    // (because NetworKit holds the global degrees in an internal vector)
    G.parallelForNodes([&](NetworKit::node u) {
        degrees[u] = G.degree(u);
    });

    return degrees;
}

/*
 * Instead of doing a DFS till depth 2 and checking if those nodes link
 * back to the original node (i.e., the NetworKit algorithm), this algorithm
 * iterates over all unique pairs of neighbors of a node and checks if these
 * neighbors are connected. Naturally, this approach only counts every
 * triangle once and not twice, as is done by NetworKit.
 *
 * Additionally, this approach can be vectorized for application on
 * the Subgraph datastructure.
 */
std::vector<double> cc_neighbor_comb(const NetworKit::Graph& G)
{
  // use NetworKit defined types
  typedef NetworKit::node node;
  typedef NetworKit::count count;

  count n = G.numberOfNodes();
  std::vector<double> coefficient(n); // $c(u) := \frac{2 \cdot |E(N(u))| }{\deg(u) \cdot ( \deg(u) - 1)}$

  G.balancedParallelForNodes([&](node u) {
    // get neighbors
    std::vector<node> neighbors;
    G.forNeighborsOf(u, [&](node v){
      if (v != u) // ignore self loops
      {
        neighbors.push_back(v);
      }
    });
    count d = neighbors.size();

    if (d < 2) {
      coefficient[u] = 0.0;
    } else {
      count triangles = 0;
      // iterate over all unique pairs of neighbors
      for (count i = 0; i < d-1;++i) {
        for (count j = i+1; j < d; ++j) {
            // count triangles if the two neighbors are connected
            if (G.hasEdge(neighbors[i], neighbors[j])) {
              triangles += 1;
            }
        }
      }
      coefficient[u] = 2.0*(double)triangles / (double)(d * (d - 1)); // No division by 2 since triangles are counted twice as well!
    }
  });

  return coefficient;

}

} // namespace algo
} // namespace tsppi

