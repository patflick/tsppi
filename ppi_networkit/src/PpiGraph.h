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

#ifndef PPIGRAPH_H
#define PPIGRAPH_H

// C includes
#include <assert.h>

// STL includes
#include <map>
#include <vector>

// NetworKit includes
#include <graph/Graph.h>

// subgraphs
#include "Subgraphs.h"
#include "Utils.h"

namespace tsppi
{

typedef uint32_t nodeid_t;
typedef NetworKit::count count_t;


class name_id_mapping
{
    /// The forward name->id map
    std::map<std::string, int> id_map;
    /// The reverse id->name map
    std::map<int, std::string> rev_id_map;
public:
    // default constructor
    name_id_mapping() = default;
    name_id_mapping(const name_id_mapping&) = default;
    // constuctors, O(n) for reversing the map
    name_id_mapping(std::map<std::string, int> name2id)
        : id_map(name2id), rev_id_map(reverse_map(name2id)) {}
    name_id_mapping(std::map<int, std::string> id2name)
        : id_map(reverse_map(id2name)), rev_id_map(id2name) {}

    /**
     * @brief Returns the name matched to the given id.
     */
    std::string getName(const int id);

    /**
     * @brief Returns all names in the name-id-mapping.
     */
    std::vector<std::string> getAllNames();

    /**
     * @brief Returns the id matched to the given name
     */
    int getId(const std::string s);

    /**
     * @brief Returns the size of the map.
     */
    std::size_t size();
};

// class for PPI graphs (non-tissue specific subgraphs)
class PpiGraph
{
public:
    /// Default constructor
    PpiGraph() = default;

    PpiGraph(const PpiGraph& other) = default;

    PpiGraph(NetworKit::Graph graph, std::map<std::string, int> id_map, std::string name)
        : graph(graph), id_map(id_map), ppi_name(name) {}
    virtual ~PpiGraph();

    /// the networKit graph instance
    NetworKit::Graph graph;

    /// mapping from gene names to integer node ids
    name_id_mapping id_map;

    /// Holds the name of the Ppi network
    std::string ppi_name;
public:
    /// Returns the name of this PPI network
    std::string getPpiName() const;

    /// Returns the gene name for the given node id
    std::string getGeneName(const nodeid_t node_id);

    /**
     * @brief Returns all gene names in the graph.
     *
     * @returns A vector of all gene names.
     */
    std::vector<std::string> getAllGenes();
};


// class for the tissue specific ppi graphs
class TsPpiGraph
{
public:
    Subgraphs subgraphs;
    name_id_mapping id_map;
    name_id_mapping tissue_id_map;
    std::string ppi_name;
public:
    TsPpiGraph() = default;
    TsPpiGraph(const TsPpiGraph& other) = default;

    TsPpiGraph(Subgraphs subgraphs, std::map<std::string, int> id_map,
            std::map<std::string, int> tissue_id_map, std::string name)
        : subgraphs(subgraphs), id_map(id_map), tissue_id_map(tissue_id_map),
          ppi_name(name) {}
    virtual ~TsPpiGraph();
public:

    // returns number of tissues
    count_t numberOfTissues() const;

    /**
     * @brief Returns the gene name associated with the given integer id.
     */
    std::string getGeneName(const nodeid_t gene_id);

    /**
     * @brief Returns the name of the tissue with the given ID.
     *
     * @param tissue_id  The tissue id for which to return the name
     *
     * @return A string containing the tissue name.
     */
    std::string getTissueName(const nodeid_t tissue_id);

};

} // namespace tsppi

#endif // PPIGRAPH_H
