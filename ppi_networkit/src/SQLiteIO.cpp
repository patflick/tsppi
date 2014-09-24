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

#include "SQLiteIO.h"
#include "Subgraphs.h"


namespace tsppi
{


TsPpiGraph SQLiteIO::load_tsppi_graph(const std::string& ppi_name, const std::string& expr_name)
{
    // first load the shared ids
    std::map<std::string, int> id_map;
    get_tsppi_ids(id_map, ppi_name, expr_name);

    // load the graph, given the shared ids
    NetworKit::Graph graph;
    load_graph(graph, id_map, ppi_name);

    // load the tissue ids of the tissue expression data
    std::map<std::string, int> tissue_id_map;
    get_tissue_id_map(tissue_id_map, expr_name);

    // load the tissue expression data
    Subgraphs::node_labels_t tissue_expr;
    load_tissue_expression(tissue_expr, id_map, tissue_id_map, expr_name);

    // constuct tsppi and subgraph
    Subgraphs sg(graph, tissue_expr);
    TsPpiGraph result(sg, id_map, tissue_id_map, ppi_name + "_" + expr_name);

    return result;
}


PpiGraph SQLiteIO::load_ppi_graph(const std::string& ppi_name)
{
    // first load the shared ids
    std::map<std::string, int> id_map;
    get_ppi_ids(id_map, ppi_name);

    // load the graph, given the shared ids
    NetworKit::Graph g;
    load_graph(g, id_map, ppi_name);

    PpiGraph ppi(g, id_map, ppi_name);
    return ppi;
}

void SQLiteIO::get_ppi_ids(std::map< std::string, int >& ppi_ids, const std::string& ppi_name)
{
    ppi_ids.clear();
    try
    {
        std::stringstream ss;
        ss << "SELECT * FROM";
        ss << "(SELECT DISTINCT Gene1 AS Gene FROM " << ppi_name;
        ss << " UNION ";
        ss << "SELECT DISTINCT Gene2 AS Gene FROM " << ppi_name << ")";
        std::string ids_query = ss.str();

        // create map for ids
        SQLite::Statement query(db, ids_query);

        // Loop to execute the query step by step, to get rows of result
        int count = 0;
        while (query.executeStep())
        {
            std::string gene_id = query.getColumn(0);
            ppi_ids[gene_id] = count++; // start counting from 0
        }
    }
    catch (std::exception& ex)
    {
        std::cout << "SQLite exception: " << ex.what() << std::endl;
        throw std::runtime_error("unable to load PPI and expression data gene intersection");
    }
}


void SQLiteIO::get_tsppi_ids(std::map<std::string, int>& tsppi_ids, const std::string& ppi_name, const std::string& expr_name)
{
    tsppi_ids.clear();
    try
    {
        std::stringstream ss;
        ss << "SELECT DISTINCT Gene FROM " << expr_name << "_core";
        ss << " INTERSECT ";
        ss << "SELECT * FROM";
        ss << "(SELECT DISTINCT Gene1 AS Gene FROM " << ppi_name;
        ss << " UNION ";
        ss << "SELECT DISTINCT Gene2 AS Gene FROM " << ppi_name << ")";
        std::string ids_query = ss.str();

        // create map for ids
        SQLite::Statement query(db, ids_query);

        // Loop to execute the query step by step, to get rows of result
        int count = 0;
        while (query.executeStep())
        {
            std::string gene_id = query.getColumn(0);
            tsppi_ids[gene_id] = count++; // start counting from 0
        }
    }
    catch (std::exception& ex)
    {
        std::cout << "SQLite exception: " << ex.what() << std::endl;
        throw std::runtime_error("unable to load PPI and expression data gene intersection");
    }
}

void SQLiteIO::load_graph(NetworKit::Graph& G, const std::map<std::string, int>& id_map, const std::string& ppi_name)
{
    try
    {
        // get size of graph
        unsigned int graph_size = id_map.size();

        // create graph
        G = NetworKit::Graph(graph_size);

        // query to get edge list
        std::string edgelist_query_str = std::string("SELECT Gene1, Gene2 FROM `") + ppi_name + "`";
        SQLite::Statement edgelist_query(db, edgelist_query_str);

        // Loop to execute the query step by step, to get rows of result
        while (edgelist_query.executeStep())
        {
            // get the result columns
            std::string gene1_id = edgelist_query.getColumn(0);
            std::string gene2_id = edgelist_query.getColumn(1);

            auto gene1_it = id_map.find(gene1_id);
            auto gene2_it = id_map.find(gene2_id);
            // if both gene IDs are part of the ID map (i.e. part of the
            // intersection between the PPI and the expression data
            if (gene1_it != id_map.end() && gene2_it != id_map.end())
            {
                int u = gene1_it->second;
                int v = gene2_it->second;
                // check if this edge has already been added
                if (!G.hasEdge(u, v) && !G.hasEdge(v, u))
                {
                    G.addEdge(u, v);
                }
            }
        }
    }
    catch (std::exception& ex)
    {
        std::cout << "SQLite exception: " << ex.what() << std::endl;
        throw std::runtime_error("unable to load graph");
    }
}

void SQLiteIO::load_tissue_expression(Subgraphs::node_labels_t& tissue_expr, const std::map<std::string, int> id_map, std::map<std::string, int> tissue_id_map, const std::string& expr_name)
{
    try
    {
        // remember the number of tissues
        const int nTissues = tissue_id_map.size();

        // query to load the core of the expression dataset
        const std::string nodeexpr_query_str = std::string("SELECT Gene, Type, Expressed FROM ") + expr_name + "_core";
        SQLite::Statement nodeexpr_query(db, nodeexpr_query_str);

        // create tissue expression matrix of correct size
        int nGenes = id_map.size();
        //std::vector<std::vector<bool>> node_expression(nGenes, std::vector<bool>(nTissues));
        tissue_expr = std::vector< boost::dynamic_bitset<> >(nGenes, boost::dynamic_bitset<>(nTissues, false));

        // loop through all results and fill matrix
        while (nodeexpr_query.executeStep())
        {
            // get columns
            std::string gene_name = nodeexpr_query.getColumn(0);
            std::string tissue_name = nodeexpr_query.getColumn(1);
            int expressed_bit = nodeexpr_query.getColumn(2);

            // check if gene is used
            auto gene_it = id_map.find(gene_name);
            if (gene_it != id_map.end())
            {
                // get gene index
                int gene_id = gene_it->second;

                // get the tissue index
                int tissue_id = tissue_id_map[tissue_name];

                // set relevant entry of the matrix
                tissue_expr[gene_id][tissue_id] = expressed_bit == 0 ? false : true;
            }
            else
            {
                // ignore: genes that are not in the PPI network are not saved
                //         into the matrix data structure
            }
        }
    }
    catch (std::exception& e)
    {
        std::cout << "SQLite exception: " << e.what() << std::endl;
        throw std::runtime_error("unable to load tissue expression data");
    }
}


template <>
std::string SQLiteIO::getSqlType<std::string>()
{
    return "TEXT";
}

template <>
std::string SQLiteIO::getSqlType<unsigned long>()
{
    return "INTEGER";
}

template <>
std::string SQLiteIO::getSqlType<int>()
{
    return "INTEGER";
}

template <>
std::string SQLiteIO::getSqlType<float>()
{
    return "REAL";
}



void SQLiteIO::get_tissue_id_map(std::map<std::string, int>& tissue_id_map, const std::string& expr_name)
{
    try
    {
        // first create tissue table
        // we are using the `core` expression table, thus there is no need to filter the tissues at this point
        std::string tissue_query_str = std::string("SELECT DISTINCT Type FROM ") + expr_name + "_core ORDER BY Type";
        SQLite::Statement tissue_query(db, tissue_query_str);

        // clear tissue id map in case it is not properly initialized
        tissue_id_map.clear();
        // next available tissue ID
        int next_tissue_id = 0;

        // loop through results and fill the map with increasing, unique integer IDs
        while (tissue_query.executeStep())
        {
            std::string tissue_name = tissue_query.getColumn(0);
            tissue_id_map[tissue_name] = next_tissue_id++; // start counting from 0
        }
    }
    catch (std::exception& e)
    {
        std::cout << "SQLite exception: " << e.what() << std::endl;
        throw std::runtime_error("unable to load tissue id mapping");
    }
}

} // namespace tsppi
