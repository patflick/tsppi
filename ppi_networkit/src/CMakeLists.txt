cmake_minimum_required(VERSION 2.6)
project(ppi_networkit)

add_library(ppi_networkit PpiGraph.cpp SQLiteIO.cpp ModularityPerCluster.cpp
                          ClusterTools.cpp Subgraphs.cpp subgraph_algos.cpp
                          graph_algos.cpp)

add_dependencies(ppi_networkit networkit-lib)

target_link_libraries(ppi_networkit NetworKit-Core-Opt SQLiteCpp sqlite3)
