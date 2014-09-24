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

// STL includes
#include <algorithm>
#include <iostream>
#include <vector>
#include <deque>
#include <stack>

// own includes
#include "PpiGraph.h"
#include "Utils.h"

/* NetworKit includes */
// clustering
#include <community/PLP.h>
#include <structures/Partition.h>

using namespace tsppi;


PpiGraph::~PpiGraph()
{

}


TsPpiGraph::~TsPpiGraph()
{

}


std::string PpiGraph::getPpiName() const
{
  return this->ppi_name;
}

std::string name_id_mapping::getName(const int id)
{
    // check if the reverse map has been set up, if not: do so
    if (this->rev_id_map.size() != this->id_map.size())
    {
        this->rev_id_map = reverse_map(this->id_map);

        if (this->rev_id_map.size() != this->id_map.size())
        {
            // this can only happen if the mapping is not a one-to-one mapping
            // which would hint towards an error elsewhere
            throw std::runtime_error("the ID mapping is not one-to-one)");
        }
    }

    // return the answer
    return this->rev_id_map[id];
}

int name_id_mapping::getId(const std::string name)
{
    return this->id_map[name];
}

std::vector<std::string> name_id_mapping::getAllNames()
{
    // get all values from the std::map
    std::vector<std::string> result(id_map.size());

    auto map_begin = id_map.begin();
    auto vec_begin = result.begin();
    // copy all values from the map into the vector
    while (map_begin != id_map.end())
    {
      *(vec_begin++) = (map_begin++)->first;
    }

    return result;
}

std::size_t name_id_mapping::size()
{
  return id_map.size();
}

std::string PpiGraph::getGeneName(const nodeid_t node_id)
{
  return id_map.getName(node_id);
}




std::vector< std::string > PpiGraph::getAllGenes()
{
  return id_map.getAllNames();
}


std::string TsPpiGraph::getGeneName(const nodeid_t gene_id)
{
  return id_map.getName(gene_id);
}

std::string TsPpiGraph::getTissueName(const nodeid_t tissue_id)
{
  return tissue_id_map.getName(tissue_id);
}

count_t TsPpiGraph::numberOfTissues() const
{
  return subgraphs.numberOfSubgraphs();
}

