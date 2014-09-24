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

#ifndef SQLITEIO_H
#define SQLITEIO_H

// STL includes
#include <map>

// SQLite includes
#include <SQLiteCpp/SQLiteCpp.h>

// NetworKit includes
#include <graph/Graph.h>

// own includes
#include <PpiGraph.h>

namespace tsppi
{

class SQLiteIO
{
private:
    /// The SQLite database filepath
    std::string db_filename;
    /// The SQLite database instance
    SQLite::Database db;

public:
    /**
     * @brief Constructor
     *
     * @param database_filename The filepath to the SQLite Database file.
     */
    SQLiteIO(const std::string & database_filename)
        : db_filename(database_filename), db(database_filename, SQLITE_OPEN_READWRITE) {}

    /// Destructor
    virtual ~SQLiteIO() {}

    // load the tissue specific graph
    TsPpiGraph load_tsppi_graph(const std::string& ppi_name, const std::string& expr_name);

    // load the global Ppi graph (no tissue specific)
    PpiGraph load_ppi_graph(const std::string& ppi_name);

    // Returns SQLite datatypes for the given C++ datatype
    template <typename T>
    static std::string getSqlType();

private:
    // returns all IDs used in the graph, mapped to an integer range
    void get_ppi_ids(std::map<std::string, int>& ppi_ids, const std::string& ppi_name);
    // returns all IDs used in the graph and the expression data set (intersection), mapped to an integer range
    void get_tsppi_ids(std::map<std::string, int>& tsppi_ids, const std::string& ppi_name, const std::string& expr_name);
    /// Loads the graph from the database.
    void load_graph(NetworKit::Graph& G, const std::map<std::string, int>& id_map, const std::string& ppi_name);

    /// Loads the tissue specific expression from the database
    void load_tissue_expression(Subgraphs::node_labels_t& tissue_expr, const std::map<std::string, int> id_map, const std::map<std::string, int> tissue_id_map, const std::string& expr_name);
    /// Loads the tissue-name -> integer mapping from the database
    void get_tissue_id_map(std::map<std::string, int>& tissue_id_map, const std::string& expr_name);
};

template <>
std::string SQLiteIO::getSqlType<std::string>();

template <>
std::string SQLiteIO::getSqlType<int>();
template <>
std::string SQLiteIO::getSqlType<unsigned long>();
template <>
std::string SQLiteIO::getSqlType<float>();


} // namespace tsppi

#endif // SQLITEIO_H
