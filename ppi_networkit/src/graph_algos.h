
/*
 * Copyright (c) 2014 Patrick Flick <patrick.flick@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 *
*/

#ifndef GRAPHS_ALGOS_H
#define GRAPHS_ALGOS_H

#include <vector>

#include <graph/Graph.h>

namespace tsppi {
namespace algo {


/**
* @brief Returns a vector of degrees for each node.
*
* @param G The NetworKit::Graph instance.
*
* @return A vector of degrees.
*/
std::vector<NetworKit::count> graph_degrees(const NetworKit::Graph& G);

/**
 * @brief   Calculates the exact clustering coefficients for each node in G.
 *
 * This calculates the clustering coefficients for each node in the given
 * graph using a different approach from the implementation in NetworKit.
 *
 * Instead of doing a DFS till depth 2 and checking if those nodes link
 * back to the original node (i.e., the NetworKit algorithm), this algorithm
 * iterates over all unique pairs of neighbors of a node and checks if these
 * neighbors are connected. Naturally, this approach only counts every
 * triangle once and not twice, as is done by NetworKit.
 *
 * Additionally, this approach can be vectorized for application on
 * the Subgraph datastructure.
 *
 * @param G The input graph.
 *
 * @return  The clustering coefficients for each node in the graph.
 */
std::vector<double> cc_neighbor_comb(const NetworKit::Graph& G);

} // namespace algo
} // namespace tsppi

#endif // GRAPHS_ALGOS_H
