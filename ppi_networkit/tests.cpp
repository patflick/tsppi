// STL includes
#include <iostream>
#include <vector>

// own includes
#include <SQLiteIO.h>
#include <PpiGraph.h>
#include <cputimer.h>
#include <ClusterTools.h>
#include <subgraph_algos.h>
#include "ppi_utils.h"

// include auxiliary networKit
#include <auxiliary/Log.h>

#define LOG(msg) std::cerr << msg << std::endl;
//#define LOG(msg)


/*********************************************************************
 *                     Testing helper functions                      *
 *********************************************************************/

template<typename T>
bool vector_eq(const std::vector<T>& v1, const std::vector<T>& v2)
{
    if (v1.size() != v2.size())
    {
        std::cerr << "ERROR: differing in size" << std::endl;
        return false;
    }
    bool all_equal = true;
    for (unsigned int i = 0; i < v1.size(); ++i)
    {
        if (v1[i] != v2[i])
        {
            std::cerr << "ERROR: i=" << i << ": " << v1[i] << " != " << v2[i] << std::endl;
            all_equal = false;
        }
    }
    return all_equal;
}

template<typename T>
bool vector_vector_eq(const std::vector<std::vector<T> >& v1, const std::vector<std::vector<T> >& v2)
{
    if (v1.size() != v2.size())
    {
        std::cerr << "ERROR: differing in size" << std::endl;
        return false;
    }
    bool all_equal = true;
    for (unsigned int j = 0; j < v1.size(); ++j)
    {
        if (v1[j].size() != v2[j].size())
        {
            std::cerr << "ERROR: inner array j=" << j << " differing in size" << std::endl;
            all_equal = false;
        }
        else if (!vector_eq(v1[j], v2[j]))
        {
            std::cerr << "ERROR: j=" << j << ": v1[j] != v2[j]" << std::endl;
            all_equal = false;
        }
    }
    return all_equal;
}


/*********************************************************************
 *                 Testing betweenness for subgraphs                 *
 *********************************************************************/

void test_ts_betweenness(tsppi::TsPpiGraph& tsppi)
{
    LOG("Testing graph: " << tsppi.ppi_name);
    LOG("Testing Betweenness: Naive");
    std::vector<std::vector<double> > naive_ts_bw = tsppi::algo::subgraph_betweenness(tsppi.subgraphs);

    LOG("Testing Betweenness: Fast");
    std::vector<std::vector<double> > fast_ts_bw = tsppi::algo::subgraph_betweenness_fast(tsppi.subgraphs);

    LOG("Checking if results are equal");
    if(vector_eq(naive_ts_bw[0], fast_ts_bw[0]))
    {
        LOG("Testing SUCCESSFUL!");
    }
    else
    {
        LOG("Testing FAILED!"); 
    }
}

void test_all_bw(tsppi::SQLiteIO& sqlio)
{
    foreach_ppi_expr([&](std::string p, std::string e){
        tsppi::TsPpiGraph ts_graph = sqlio.load_tsppi_graph(p, e);
        test_ts_betweenness(ts_graph);
    });
}


/*********************************************************************
 *           Testing clustering coefficient for subgraphs            *
 *********************************************************************/

void test_ts_clustercoeff(tsppi::TsPpiGraph& tsppi)
{
    LOG("Testing graph: " << tsppi.ppi_name);
    LOG("Testing CC: Naive");
    std::vector<std::vector<double> > naive_ts_bw = tsppi::algo::subgraph_cc(tsppi.subgraphs);

    LOG("Testing CC: Fast NeighborComb");
    std::vector<std::vector<double> > fast_ts_bw = tsppi::algo::subgraph_cc_neighbor_comb(tsppi.subgraphs);

    LOG("Testing CC: Fast TsGraph");
    std::vector<std::vector<double> > faster_ts_bw = tsppi::algo::subgraph_cc_neighbor_comb_vec(tsppi.subgraphs);

    // check that the results of the new algorithms are identical
    LOG("Checking if results are equal");
    if(vector_vector_eq(naive_ts_bw, fast_ts_bw) && vector_vector_eq(naive_ts_bw, faster_ts_bw))
    {
        LOG("Testing SUCCESSFUL!");
    }
    else
    {
        LOG("Testing FAILED!");
    }
}

void test_all_cc(tsppi::SQLiteIO& sqlio)
{
    foreach_ppi_expr([&](std::string p, std::string e){
        tsppi::TsPpiGraph ts_graph = sqlio.load_tsppi_graph(p, e);
        test_ts_clustercoeff(ts_graph);
    });
}


/*********************************************************************
 *                        main() entry point                         *
 *********************************************************************/

void print_usage()
{
    std::cerr << "Usage: ./benchmark [OPTIONS] <sqlite_db_filename>" << std::endl;
    std::cerr << "  Options:" << std::endl;
    std::cerr << "      -a  The algorithm to test, can be: `cc`, or `bw`" << std::endl;
    std::cerr << "      -p  (optional) The PPI to test this on. `-e` has to be set as well" << std::endl;
    std::cerr << "      -e  (optional) The Expression data to test this on. `-p` has to be set as well" << std::endl;
    std::cerr << "           If neigher -p or -e is set, the given algorithm is tested" << std::endl;
    std::cerr << "           on all PPIs and expression datasets." << std::endl;
}

int main(int argc, char *argv[])
{
    --argc;
    ++argv;
    if (argc < 1){
        print_usage();
        exit(1);
    }

    std::string algo;
    std::string ppi = "";
    std::string expr = "";
    std::string sqlite_db_file = "";

    // parse optional arguments
    while(argc > 0){
        if (argv[0][0] == '-')
        {
            switch(argv[0][1])
            {
                case 'a':
                    // set algorithm to benchmark
                    algo = argv[1];
                    ++argv;
                    --argc;
                    break;
                case 'p':
                    // set PPI to load
                    ppi = argv[1];
                    ++argv;
                    --argc;
                    break;
                case 'e':
                    // set expr to load
                    expr = argv[1];
                    ++argv;
                    --argc;
                    break;
                default:
                    print_usage();
                    exit(1);
            }
            --argc;
            ++argv;
        }
        else
        {
            // this is the database file
            sqlite_db_file = argv[0];
            ++argv;
            --argc;
        }
    }

    if (sqlite_db_file == "" || algo == "")
    {
        print_usage();
        exit(1);
    }

    // open Sqlite database file
    tsppi::SQLiteIO sqlio(sqlite_db_file);

    if (ppi != "")
    {
        assert(expr != "");

        tsppi::TsPpiGraph ts_graph = sqlio.load_tsppi_graph(ppi, expr);
        if (algo == "bw")
        {
            test_ts_betweenness(ts_graph);
        } else if (algo == "cc")
        {
            test_ts_clustercoeff(ts_graph);
        } else {
            print_usage();
            exit(1);
        }
    }
    else
    {
        // benchmark all ppi and expr
        if (algo == "bw")
        {
            test_all_bw(sqlio);
        } else if (algo == "cc")
        {
            test_all_cc(sqlio);
        } else {
            print_usage();
            exit(1);
        }
    }

    return 0;
}
