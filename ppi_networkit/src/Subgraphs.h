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

#ifndef SUBGRAPHS_H
#define SUBGRAPHS_H

// STL includes
#include <vector>

// include boost::dynamic_bitset, which is used for node labels
#include <boost/dynamic_bitset.hpp>

// NetworKit includes
#include <graph/Graph.h>

namespace tsppi
{

class Subgraphs
{
public:
    // type of a node label (bit vector), defining whether a node is expressed
    // in a given subgraph index
    typedef boost::dynamic_bitset<> node_label_t;
    typedef std::vector< node_label_t > node_labels_t;
    typedef NetworKit::count count_t;

    // the NetworKit graph instance
    NetworKit::Graph graph;

    // binary vectors per node, stating which nodes are existant in which
    // subgraph
    node_labels_t node_exists;

    /**
     * @brief Returns the subgraph with the given offset as a new graph instance.
     *
     * @param id The subgraph index
     *
     * @return  A new NetworKit::Graph instance for the given subgraph index.
     */
    NetworKit::Graph getSubgraph(count_t id) const;

    /**
     * @brief Returns the number of subgraphs.
     *
     * @return The number of subgraphs.
     */
    count_t numberOfSubgraphs() const {
        return nSubgraphs;
    }

    /**
     * @brief Returns the number of nodes in the full graph.
     */
    count_t numberOfNodes() const {
        return graph.numberOfNodes();
    }

    /**
     * @brief Returns the number of edges in the full graph.
     */
    count_t numberOfEdges() const {
        return graph.numberOfEdges();
    }

    /**
     * Creates new graph instances for each subgraph and executes the given
     * function on the subgraph.
     *
     * @tparam L  The function type with signature (NetworKit::Graph&, int)
     * @param handle The function to be executed for each subgraph.
     */
    template<typename L>
    void foreach_subgraph(L handle) const {
        for (count_t t = 0; t < nSubgraphs; ++t)
        {
          // get copy of subgraph with current offset
          NetworKit::Graph subgraph = getSubgraph(t);
          // call given handler
          handle(subgraph, t);
        }
    }

    /**
     * Creates new graph instances for each subgraph and executes the given
     * function on the subgraph.
     *
     * Uses OpenMP parallelism for creating the subgraphs and calling
     * the given function.
     *
     * @tparam L  The function type with signature (NetworKit::Graph&, int)
     * @param handle The function to be executed for each subgraph.
     */
    template<typename L>
    void parallel_foreach_subgraph(L handle) const {
#ifdef BENCHMARK_USE_OMP
#pragma omp parallel for schedule(runtime)
#endif
        for (count_t t = 0; t < nSubgraphs; ++t)
        {
          // get copy of subgraph with current offset
          NetworKit::Graph subgraph = getSubgraph(t);
          // call given handler
          handle(subgraph, t);
        }
    }

    /// Default constructor
    Subgraphs() = default;
    Subgraphs(const Subgraphs&) = default;

    /**
     * @brief Constuctor for the subgraphs class.
     *
     * @param graph         The NetworKit::Graph instance to wrap around.
     * @param node_exists   The node labels.
     */
    Subgraphs (NetworKit::Graph graph, node_labels_t node_exists)
        : graph(graph), node_exists(node_exists), nSubgraphs(node_exists[0].size()) {}

    /// Destructor
    virtual ~Subgraphs () {}

private:
    // The number of subgraphs
    count_t nSubgraphs;
};

} // namespace tsppi

#endif // SUBGRAPHS_H
