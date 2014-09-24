
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

#ifndef SUBGRAPHS_ALGOS_H
#define SUBGRAPHS_ALGOS_H

#include <vector>

#include "Subgraphs.h"

// NetworKit Paritition (for clustering)
#include <graph/Graph.h>
#include <structures/Partition.h>

namespace tsppi {
namespace algo {

/*********************************************************************
 *                              Degrees                              *
 *********************************************************************/

/**
 * @brief Returns the degrees of all nodes in all subgraphs
 *
 * @param subgraphs The subgraph datastructure
 *
 * @return The degrees for all nodes in all subgraphs.
 */
std::vector<std::vector<NetworKit::count> > subgraph_degrees(const Subgraphs& subgraphs);

/**
 * @brief Returns for each node the maximum degree of that node across all subgraphs.
 *
 * @param subgraphs The subgraph datastructure
 *
 * @return For each node the max degree of that node across all subgraphs.
 */
std::vector<NetworKit::count> subgraph_max_degrees(const Subgraphs& subgraphs);


/**
 * @brief Returns for each node the number of edges which exist in at least
 *        one subgraph. We call this the co-existence degree.
 *
 * @param subgraphs The subgraph datastructure
 *
 * @return For each node the number of edges which exist in at least
 *         one subgraph.
 */
std::vector<NetworKit::count> subgraph_edge_exist_degrees(const Subgraphs& subgraphs);


/*********************************************************************
 *     maximum and minumum exists overlap with neighboring nodes     *
 *********************************************************************/


/**
 * @brief Returns for each node the minimum co-existance count across all neighbors.
 *
 * @param subgraphs The subgraph datastructure.
 *
 * @return For each node the minimum co-existance count across all neighbors.
 */
std::vector<NetworKit::count> subgraph_neighbor_min_exists_count(const Subgraphs& subgraphs);

/**
 * @brief Returns for each node the maximum co-existance count across all neighbors.
 *
 * @param subgraphs The subgraph datastructure.
 *
 * @return For each node the maximum co-existance count across all neighbors.
 */
std::vector<NetworKit::count> subgraph_neighbor_max_exists_count(const Subgraphs& subgraphs);



/*********************************************************************
 *                         Aggregate graphs                          *
 *********************************************************************/

/**
 * @brief Creates a aggregate graph from all subgraphs. Edge weights are set
 *        to be equal to the number of subgraphs in which the edge exists.
 *
 * @param subgraphs The subgraph datastructure.
 *
 * @return The aggregate graph.
 */
NetworKit::Graph subgraph_edge_coexist_count_graph(const Subgraphs& subgraphs);

/**
 * @brief Creates a aggregate graph from all subgraphs. Edge weights are set
 *        to a correlation score of the two interacting node's existance
 *        vectors.
 *
 * @param subgraphs The subgraph datastructure.
 *
 * @return The aggregate graph.
 */
NetworKit::Graph subgraph_edge_correlation_graph(const Subgraphs& subgraphs);


/*********************************************************************
 *               Clustering Coefficients for Subgraphs               *
 *********************************************************************/

/*
 * @brief Calculates and returns the clustering coefficients for all nodes
 *        in all subgraphs for the given subgraph.
 *
 * Simplest implementation: create subgraphs and run NetworKit's CC algo
 *
 * @param subgraphs The subgraph datastructure
 *
 * @returns The clustering coeffiecients for all nodes in all subgraphs.
 */
std::vector<std::vector< double > > subgraph_cc(const Subgraphs& subgraphs);

/*
 * @brief Calculates and returns the clustering coefficients for all nodes
 *        in all subgraphs for the given subgraph.
 *
 * Create subgraphs and run our CC implementation (neighbor combinations)
 * for each one.
 *
 * @param subgraphs The subgraph datastructure
 *
 * @returns The clustering coeffiecients for all nodes in all subgraphs.
 */
std::vector<std::vector< double > > subgraph_cc_neighbor_comb(const Subgraphs& subgraphs);

/*
 * @brief Calculates and returns the clustering coefficients for all nodes
 *        in all subgraphs for the given subgraph.
 *
 * `Vectorized` version of the neighbor combination CC algorithm, that
 * calculates the clustering coefficients for all nodes in all subgraphs
 * at the same time, without having to create a graph instance for each
 * subgraph.
 *
 * @param subgraphs The subgraph datastructure
 *
 * @returns The clustering coeffiecients for all nodes in all subgraphs.
 */
std::vector<std::vector< double > > subgraph_cc_neighbor_comb_vec(const Subgraphs& subgraphs);




/*********************************************************************
 *                     Betweenness for Subgraphs                     *
 *********************************************************************/

/**
 * @brief Returns the betweenness for all nodes in all tissue specific
 *        subnetworks.
 *
 * Simplest implementation: create subgraphs and run NetworKit's betweenness
 * algorithm.
 *
 * @return The betweenness centrality in form for a vector of vectors
 */
std::vector<std::vector<double> > subgraph_betweenness(const Subgraphs& subgraphs);

/**
 * @brief Returns the betweenness for all nodes in all tissue specific
 *        subnetworks.
 *
 * This implements the betweenness calculation for the custom
 * TSPPI graph datastructure, so that it is not necessary to copy
 * and create subgraphs for every tissue.
 *
 * @return The betweenness centrality in form for a vector of vectors
 */
std::vector<std::vector<double> > subgraph_betweenness_fast(const Subgraphs& subgraphs);




/*********************************************************************
 *                   PLP clustering for Subgraphs                    *
 *********************************************************************/

/**
 * @brief Runs the PLP community detection algorithm and returns the
 *        detected partitions.
 *
 * Simplest implementation: create subgraphs and run NetworKit's betweenness
 * algorithm.
 *
 * @param subgraphs The Subgraphs datastructure instance.
 *
 * @return  The partitions for each subgraph.
 */
std::vector<NetworKit::Partition> subgraph_PLP(const Subgraphs& subgraphs);

/**
 * @brief Runs the PLP community detection algorithm and returns the
 *        detected partitions.
 *
 * Modified version of the PLP algorithm that runs directly on the
 * Subgraphs datastructure by checking if an edge exists for each traversal.
 *
 * @param subgraphs The Subgraphs datastructure instance.
 *
 * @return  The partitions for each subgraph.
 */
std::vector<NetworKit::Partition> subgraph_PLP_vec(const Subgraphs& subgraphs,
                            NetworKit::count updateThreshold = NetworKit::none);

} // namespace algo
} // namespace tsppi

#endif // SUBGRAPHS_ALGOS_H
