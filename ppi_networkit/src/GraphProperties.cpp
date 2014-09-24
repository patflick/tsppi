
#include <GraphProperties.h>
#include <iostream>
#include <vector>
#include <string>

void calc_and_save_node_properties(std::string& ppi_name, std::string& expr_name, tsppi::SQLiteIO& sqlio)
{
  tsppi::TsPpiGraph graph;
  std::cout << "Load TSPPI graph from SQLite connection ..." << std::endl;
  sqlio.load_tsppi_graph(graph, ppi_name, expr_name);

  std::cout << "size of id_map: " << graph.id_map.size() << std::endl;
  std::cout << "size of tissue_id_map: " << graph.tissue_id_map.size() << std::endl;
  std::cout << "size of tissue_expr: " << graph.tissue_expr.size() << std::endl;
  std::cout << "size of graph: " << graph.graph.numberOfNodes() << std::endl;

  graph.printBasicGraphInfo();
  //get_degrees(graph);
  std::vector< std::vector< tsppi::count_t > > ts_degrees;

  // get max degrees of tissue specific network:
  std::cout << "getting global degrees ..." << std::endl;
  std::vector< tsppi::count_t > degrees;
  graph.getDegrees(degrees);
  std::cout << "saving global degrees ..." << std::endl;
  sqlio.save_node_property(graph, degrees, "degree", ppi_name, expr_name);

  // get max degrees of tissue specific network:
  std::cout << "getting max ts degrees ..." << std::endl;
  degrees.clear();
  graph.getMaxTsDegree(degrees);
  std::cout << "saving max ts degrees ..." << std::endl;
  sqlio.save_node_property(graph, degrees, "maxts_degree", ppi_name, expr_name);

  // get coexpr degrees of tissue specific network:
  std::cout << "getting coexpr degrees ..." << std::endl;
  degrees.clear();
  graph.getCoexprDegree(degrees);
  std::cout << "saving coexpr degrees ..." << std::endl;
  sqlio.save_node_property(graph, degrees, "coexpr_degree", ppi_name, expr_name);

  std::cout << "getting min/max neighbor expression count ..." << std::endl;
  std::vector<tsppi::count_t> mins;
  std::vector<tsppi::count_t> maxs;
  graph.getNeighborMinMaxExprCount(mins, maxs);
  // save to SQL
  std::cout << "saving min/max neighbor expression count ..." << std::endl;
  sqlio.save_node_property(graph, mins, "min_neighbor_expr_count", ppi_name, expr_name);
  sqlio.save_node_property(graph, maxs, "max_neighbor_expr_count", ppi_name, expr_name);
}

void getAllNodeProperties(tsppi::SQLiteIO& sqlio)
{
    std::vector<std::string> ppis = {"bossi", "string", "psicquic_all", "havu", "ccsb"};
    std::vector<std::string> exprs = {"emtab", "gene_atlas", "rnaseq_atlas", "hpa", "hpa_all"};

    for (auto ppi_name : ppis)
    {
        for (auto expr_name : exprs)
        {
            std::cout << "##################################################################" << std::endl;
            std::cout << "  PPI = `" << ppi_name << "`, EXPR = `" << expr_name << "`" << std::endl;
            std::cout << "##################################################################" << std::endl;
            calc_and_save_node_properties(ppi_name, expr_name, sqlio);
        }
    }
}
